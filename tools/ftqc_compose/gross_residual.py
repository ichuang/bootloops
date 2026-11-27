"""Pauli displacements at the open output of a BB144 logical Z measurement."""

from functools import lru_cache

from ftqc_faults.audit.record_model import RecordModel
from ftqc_codes.gf2 import Echelon, dot, rank
from ftqc_extractors.__main__ import verify
from ftqc_extractors.algebra import bits, kernel
from ftqc_extractors.gross import gross_z, logical_action, port_frame, shifted_z
from ftqc_extractors.gross_probes import experiment, quotient

from .supports import fields


def pack(value, rows):
    return sum(dot(value, row) << j for j, row in enumerate(rows))


def check_kernel(rows, trivial, n):
    if rank(rows) + rank(trivial) != n or any(dot(a, b) for a in rows for b in trivial):
        raise ValueError("terminal tests do not separate the projective Pauli quotient")


class Probe:
    def __init__(self, basis, rounds, mem, edges, rows, cycles, expected_hash):
        circ = experiment(rounds, basis)
        if circ.sha256() != expected_hash:
            raise ValueError("distance certificate covers a different probe")
        end = len(circ.meas) - mem.n
        split = end - edges
        expected = [(q, "X") for q in range(mem.n, mem.n + edges)]
        expected += [(q, basis) for q in range(mem.n)]
        if [(q, b) for _, _, q, b in circ.meas[split:]] != expected:
            raise ValueError("unexpected split or terminal record layout")
        self.model = RecordModel(circ, response_records=range(split, len(circ.meas)))
        self.offset = self.model.nd + len(circ.observables)
        self.checks = len(mem.hx if basis == "X" else mem.hz)
        self.operational = self.model.nd - self.checks
        self.history_mask = (1 << self.operational) - 1
        self.cycle_start = self.operational - len(cycles)
        self.cycles = cycles
        self.read_bit = int(basis == "Z")
        self.edges = edges
        self.projection = lru_cache(maxsize=32768)(lambda error: pack(error, rows))

    def raw(self, value):
        return divmod(value >> self.offset, 1 << self.edges)

    def check(self, value, corrected, split):
        terminal = (value & self.model.mask) >> self.operational
        observables = (value % (1 << self.offset)) >> (self.model.nd + self.read_bit)
        if terminal | (observables << self.checks) != self.projection(corrected):
            raise ValueError("frame-corrected residual disagrees with terminal probe")
        if self.cycles and (
            (value & self.history_mask) >> self.cycle_start != pack(split, self.cycles)
        ):
            raise ValueError("split-cycle displacement disagrees with the probe")


def residuals(rounds, distance):
    mem, _, ports, edges = gross_z()
    target = shifted_z()
    sites = list(bits(target))
    paths, cycles = port_frame(ports, mem.n, len(edges))
    z_rows = mem.hz + quotient([target, *kernel(mem.hx, mem.n)], mem.hz)
    x_rows = mem.hx + quotient(kernel([*mem.hz, target], mem.n), mem.hx)
    trivial_z = [*mem.hz, target]
    check_kernel(z_rows, mem.hx, mem.n)
    check_kernel(x_rows, trivial_z, mem.n)
    bx, bz = Echelon(mem.hx), Echelon(trivial_z)

    @lru_cache(maxsize=32768)
    def frame(split):
        return sum(dot(split, path) << q for q, path in zip(sites, paths, strict=True))

    @lru_cache(maxsize=32768)
    def representative(x, z):
        # Equivalent representatives suffice for a constructive upper bound.
        return min(
            (a | b).bit_count()
            for a in (x, bx.reduce(x))
            for b in (z, z ^ target, bz.reduce(z))
        )

    xp = Probe(
        "X",
        rounds,
        mem,
        len(edges),
        x_rows,
        cycles,
        distance["probes"]["X"]["circuit_sha256"],
    )
    zp = Probe(
        "Z",
        rounds,
        mem,
        len(edges),
        z_rows,
        [],
        distance["probes"]["Z"]["circuit_sha256"],
    )
    accepted = total = framed = nontrivial = largest = max_open = 0
    examples = {}
    first_x = xp.model.locations[0][0]
    first_z = zp.model.locations[0][0]
    for i, (a, b) in enumerate(
        zip(xp.model.locations, zp.model.locations, strict=True)
    ):
        if (a[0] - first_x, a[1], a[2].text()) != (b[0] - first_z, b[1], b[2].text()):
            raise ValueError("fault location mismatch between complementary probes")
        gx, gz = xp.model.groups[i], zp.model.groups[i]
        if gx.keys() != gz.keys():
            raise ValueError("fault mode mismatch")
        for mode, vx in gx.items():
            vz = gz[mode]
            raw_z, split = xp.raw(vx)
            raw_x, other_split = zp.raw(vz)
            if split != other_split:
                raise ValueError("split records disagree between probes")
            fz = frame(split)
            z = raw_z ^ fz
            xp.check(vx, z, split)
            zp.check(vz, raw_x, split)
            weight = representative(raw_x, z)
            total += 1
            largest = max(largest, weight)
            if (vx & xp.history_mask) or (vz & zp.history_mask):
                continue
            accepted += 1
            framed += bool(fz)
            nontrivial += bool(weight)
            labels = []
            if weight > max_open:
                max_open = weight
                labels.append("largest_accepted_representative")
            if fz and "accepted_frame" not in examples:
                labels.append("accepted_frame")
            if labels:
                witness = {
                    "fault": [i, mode],
                    "op": a[2].text(),
                    "split_flips": list(bits(split)),
                    "raw_x": list(bits(raw_x)),
                    "raw_z": list(bits(raw_z)),
                    "frame_z": list(bits(fz)),
                    "corrected_z": list(bits(z)),
                    "representative_weight": weight,
                    "readout_flip": (vz >> zp.model.nd) & 1,
                }
                for label in labels:
                    examples[label] = witness
    for example in examples.values():
        for probe in (xp, zp):
            probe.model.replay([example["fault"]])
    return {
        "noisy_locations": len(xp.model.locations),
        "single_location_modes": total,
        "operationally_accepted_modes": accepted,
        "accepted_with_frame_update": framed,
        "accepted_with_nontrivial_residual": nontrivial,
        "max_accepted_representative_weight": max_open,
        "all_modes_spread_upper": largest,
        "witnesses": examples,
        "zero_test_kernel_ranks": {"X": bx.rank, "Z": bz.rank},
    }


def gross_result(spec):
    fields(spec, "kind", "extractor", "distance_receipt")
    source, receipt = spec["extractor"], spec["distance_receipt"]
    if source["kind"] != "gross_circuit":
        raise ValueError("expected a gross logical Z circuit")
    verify(source, receipt)
    result = receipt["result"]
    if result["status"] != "proved":
        raise ValueError("gross circuit needs a proved distance lower bound")
    logical_action()
    return {
        "status": "proved",
        "scope": "clean encoded input, diagonal location-Pauli faults",
        "rounds": source["rounds"],
        "distance": [result["lower"], result["upper"]],
        "open_output": residuals(source["rounds"], result),
        "operational_acceptance": "syndrome history and split cycles",
        "residual_rule": "XOR raw data Paulis and split-induced Z frames; at t faults a representative has weight at most min(144, t * all_modes_spread_upper)",
        "closed_output_through": result["lower"] - 1,
        "closed_output_condition": "operational acceptance and zero ideal terminal syndrome",
        "closed_output_conclusion": "correct measured bit; residual is a memory stabilizer times a power of the measured Z on its projected branch",
    }
