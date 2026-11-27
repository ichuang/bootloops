"""CSS record probes of a one-round full block reading with explicit gates and idles."""

from ftqc_faults.audit.circuit import Circuit
from ftqc_faults.audit.record_model import RecordModel
from ftqc_codes.codes import CSSCode, load
from .algebra import bits, kernel, parity, transpose
from .block_reading import combine, quotient, reading
from .gross_probes import encoder
from . import css_cover


def linear_frame(section, checks):
    # Gaussian elimination with data corrections carried along. Extend the
    # right inverse linearly to all records; invalid split parities are detectors.
    basis = {}
    for value, word in section:
        for pivot in sorted(basis, reverse=True):
            if value >> pivot & 1:
                value ^= basis[pivot][0]
                word ^= basis[pivot][1]
        basis[value.bit_length() - 1] = (value, word)
    columns = []
    for i in range(checks):
        value, word = 1 << i, 0
        for pivot in sorted(basis, reverse=True):
            if value >> pivot & 1:
                value ^= basis[pivot][0]
                word ^= basis[pivot][1]
        columns.append(word)
    return columns


def meter_layers(checks, width, basis):
    """Greedy disjoint edge layers; every live idle during extraction is noisy."""
    meters = list(range(width, width + len(checks)))
    layers = [
        [("RX" if basis == "X" else "R", q) for q in meters]
        + [("I", q) for q in range(width)]
    ]
    slots = []
    used = []
    for meter, row in zip(meters, checks):
        for q in bits(row):
            k = next(
                (
                    j
                    for j, sites in enumerate(used)
                    if q not in sites and meter not in sites
                ),
                len(used),
            )
            if k == len(used):
                used.append(set())
                slots.append([])
            used[k].update((q, meter))
            slots[k].append(("CX", meter, q) if basis == "X" else ("CX", q, meter))
    for gates, sites in zip(slots, used):
        layers.append(
            gates + [("I", q) for q in range(width + len(meters)) if q not in sites]
        )
    layers.append(
        [("MX" if basis == "X" else "M", q) for q in meters]
        + [("I", q) for q in range(width)]
    )
    return layers


def probe(mem, axis, basis):
    # Work in Z-reading coordinates; the X experiment is its global Hadamard dual.
    n = mem.n
    read = reading((n, mem.hx, mem.hz), 2, [3], "Z")
    data = 2 * n
    width = read["width"]
    layers = []
    for tick in encoder(mem, basis):
        layers.append(tick + [(op, *(q + n for q in qs)) for op, *qs in tick])
    layers.append(
        [("RX", q) for q in range(data, width)] + [("I", q) for q in range(data)]
    )
    layers.extend(meter_layers(read["hx"], width, "X"))
    layers.extend(meter_layers(read["hz"], width, "Z"))
    layers.append(
        [("MX", q) for q in range(data, width)] + [("I", q) for q in range(data)]
    )
    layers.append([("!M" if basis == "Z" else "!MX", q) for q in range(data)])
    nx = len(read["hx"])
    oldz = 2 * len(mem.hz)
    split = nx + len(read["hz"])
    end = split + width - data
    det = [[i] for i in range(nx + oldz)]
    det.extend([nx + i for i in bits(m)] for m in read["metachecks"])
    det.extend(
        [split + i for i in bits(m)] for m in kernel(transpose(mem.hx, n), len(mem.hx))
    )
    frame = linear_frame(read["split_basis"], len(mem.hx))

    def terminal(row):
        records = [end + i for i in bits(row)]
        if basis == "X":
            records += [
                split + i for i, v in enumerate(frame) if parity(row, v | (v << n))
            ]
        return records

    old_x = [h << j for j in (0, n) for h in mem.hx]
    old_z = [h << j for j in (0, n) for h in mem.hz]
    det.extend(terminal(row) for row in (old_x if basis == "X" else old_z))
    if basis == "Z":
        obs = [[nx + oldz + i for i in bits(m)] for m in read["port_readouts"]]
        reps = quotient(old_x, old_z, data)
    else:
        obs = []
        reps = quotient(old_z + read["measured"], old_x, data)
    obs.extend(terminal(row) for row in reps)
    if axis == "X":
        dual = {"R": "RX", "RX": "R", "M": "MX", "MX": "M", "I": "I"}

        def swap(op, *qs):
            name = op.lstrip("!")
            if name == "CX":
                return op, *reversed(qs)
            return op.replace(name, dual[name]), *qs

        layers = [[swap(*op) for op in tick] for tick in layers]
    return (
        Circuit(
            layers,
            detectors=det,
            observables=obs,
            outputs=[],
            name=f"BB144 {axis} block reading {basis} canonical probe",
        ),
        split,
        end,
        frame,
    )


def audit_probe(mem, axis, basis):
    circ, split, end, frame = probe(mem, axis, basis)
    model = RecordModel(circ, response_records=range(split, len(circ.meas)))
    obs_mask = (1 << len(circ.observables)) - 1
    read_offset = model.nd + len(circ.observables)
    splitmask = (1 << (end - split)) - 1
    upper = None
    signatures = {}
    spread = 0
    for loc, group in enumerate(model.groups):
        current = []
        for fault, value in group.items():
            det = value & model.mask
            obs = (value >> model.nd) & obs_mask
            if not det and obs:
                upper = [(loc, fault)]
            elif upper is None:
                for prev, other, prev_obs in signatures.get(det, ()):
                    if obs != prev_obs:
                        upper = [(prev, other), (loc, fault)]
                        break
            current.append((det, obs, fault))
            split_data = value >> read_offset
            correction = combine(frame, split_data & splitmask) if basis == "X" else 0
            residual = (
                (split_data >> (end - split)) ^ correction ^ (correction << mem.n)
            )
            spread = max(spread, residual.bit_count())
        for det, obs, fault in current:
            seen = signatures.setdefault(det, [])
            if len(seen) < 2 and all(old != obs for _, _, old in seen):
                seen.append((loc, fault, obs))
    lower = len(upper) if upper is not None else 3
    if upper is None:
        # A full memory logical applied at the last idle evades terminal checks.
        reps = quotient(
            mem.hz if basis == "Z" else mem.hx,
            mem.hx if basis == "Z" else mem.hz,
            mem.n,
        )
        from .gross import shifted_z, LEFT_Z, RIGHT_Z, cell

        dual = sum(
            1 << (half * 72 + cell(-i, -j))
            for half, terms in enumerate((RIGHT_Z, LEFT_Z))
            for i, j in terms
        )
        physical_fault = "X" if (axis == "Z") == (basis == "Z") else "Z"
        candidate = dual if physical_fault == "X" else shifted_z()
        for rep in [candidate, *reps]:
            where = {
                qs[0]: i
                for i, (t, _, op) in enumerate(model.locations)
                if t == len(circ.ticks) - 2 and op.name == "I"
                for qs in [op.qs]
            }
            fault = "X" if (basis == "Z") == (axis == "Z") else "Z"
            witness = [(where[q], fault) for q in bits(rep)]
            det, obs = model.replay(witness)
            if det == 0 and obs & obs_mask:
                upper = witness
                break
    if upper is None:
        raise ValueError("missing circuit upper witness")
    det, obs = model.replay(upper)
    if det or not obs & obs_mask:
        raise ValueError("physical witness failed forward replay")
    return dict(
        circuit_sha256=circ.sha256(),
        locations=len(model.locations),
        lower=lower,
        upper=len(upper),
        witness=[list(p) for p in upper],
        residual_component_weight_per_fault=spread,
    )


def analyze(s, budget):
    if s["axis"] not in ("X", "Z") or s["target"] != 12 or s["padding"] != [12, 12]:
        raise ValueError(
            "BB144 benchmark requires distance 12 and two 12-round boundaries"
        )
    original = load("bb144")
    mem = (
        original if s["axis"] == "Z" else CSSCode("dual", 144, original.hz, original.hx)
    )
    for checks, stabs in ((mem.hx, mem.hz), (mem.hz, mem.hx)):
        short, nodes, _ = css_cover.cover(144, checks, stabs, 11, budget)
        budget -= nodes
        if short:
            raise ValueError("BB144 memory distance below twelve")
    read = reading((144, mem.hx, mem.hz), 2, [3], "Z")
    probes = {b: audit_probe(mem, s["axis"], b) for b in ("X", "Z")}
    return dict(
        status="proved",
        axis=s["axis"],
        measured=read["measured"],
        survivors=read["survivors"],
        added_data=read["added_data"],
        added_check_meters=read["added_check_meters"],
        phenomenological=dict(
            theorem="2510.14895 Proposition 3.9",
            lower=12,
            padding=[12, 12],
            merge_rounds=1,
        ),
        circuit=dict(
            model="one merge round, ideal encoded inputs and terminal checks",
            probes=probes,
            lower=min(p["lower"] for p in probes.values()),
            upper=min(p["upper"] for p in probes.values()),
        ),
    )
