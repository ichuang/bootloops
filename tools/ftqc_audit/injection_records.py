"""Physical location responses of Clifford logical-Y injection experiments."""

from collections import Counter, defaultdict
from fractions import Fraction


def parse(circuit):
    """Resolve absolute measurement records and retain each measured Pauli product."""
    records = []
    ops = []
    detectors = []
    obs = {}
    for inst in circuit.flattened():
        g = inst.name
        args = inst.gate_args_copy()
        targets = inst.targets_copy()
        if g in ("DETECTOR", "OBSERVABLE_INCLUDE"):
            rr = []
            for t in targets:
                if not t.is_measurement_record_target or not 0 <= len(
                    records
                ) + t.value < len(records):
                    raise ValueError("invalid measurement reference")
                rr.append(len(records) + t.value)
            if g == "DETECTOR":
                detectors.append(rr)
            else:
                obs.setdefault(int(args[0]), []).extend(rr)
            continue
        if g == "TICK":
            continue
        start = len(records)
        if g == "MPP":
            for group in inst.target_groups():
                records.append(
                    tuple(
                        (
                            "X" if t.is_x_target else "Y" if t.is_y_target else "Z",
                            t.value,
                        )
                        for t in group
                    )
                )
            ops.append((g, (), Fraction(0), start, len(records), len(ops)))
            continue
        qs = tuple(t.value for t in targets)
        if g in ("MX", "M", "MR"):
            records.extend((("X" if g == "MX" else "Z", q),) for q in qs)
        elif g not in (
            "H",
            "CX",
            "R",
            "RX",
            "RY",
            "DEPOLARIZE1",
            "DEPOLARIZE2",
            "X_ERROR",
            "Y_ERROR",
            "Z_ERROR",
        ):
            raise ValueError(f"unsupported physical operation {g}")
        ops.append(
            (
                g,
                qs,
                Fraction(str(args[0])) if args else Fraction(0),
                start,
                len(records),
                len(ops),
            )
        )
    if sorted(obs) != list(range(len(obs))):
        raise ValueError("expected dense observable IDs")
    return records, ops, detectors, [obs[i] for i in range(len(obs))]


def responses(circuit, parsed):
    records, ops, det, obs = parsed
    columns = [0] * len(records)
    for i, row in enumerate(det + obs):
        for j in row:
            columns[j] ^= 1 << i
    x, z = [0] * circuit.num_qubits, [0] * circuit.num_qubits
    for gate, qs, p, start, end, op in reversed(ops):
        if p and gate in ("M", "MX", "MR"):
            for i, q in enumerate(qs):
                yield op, (q,), "flip", p, columns[start + i]
        if gate == "MR":
            for q in qs:
                x[q] = z[q] = 0
        if end > start:
            for i in range(start, end):
                for axis, q in records[i]:
                    if axis in ("Z", "Y"):
                        x[q] ^= columns[i]
                    if axis in ("X", "Y"):
                        z[q] ^= columns[i]
        elif gate in ("R", "RX", "RY"):
            for q in qs:
                x[q] = z[q] = 0
        elif gate == "H":
            for q in qs:
                x[q], z[q] = z[q], x[q]
        elif gate == "CX":
            for c, t in reversed(list(zip(qs[::2], qs[1::2]))):
                x[c] ^= x[t]
                z[t] ^= z[c]
        elif p:
            import itertools

            width = 2 if gate == "DEPOLARIZE2" else 1
            modes = (
                [
                    "".join(w)
                    for w in itertools.product("IXYZ", repeat=width)
                    if any(a != "I" for a in w)
                ]
                if gate.startswith("DEPOLARIZE")
                else [gate[0]]
            )
            for i in range(0, len(qs), width):
                sites = qs[i : i + width]
                for mode in modes:
                    value = 0
                    for q, a in zip(sites, mode):
                        if a in "XY":
                            value ^= x[q]
                        if a in "ZY":
                            value ^= z[q]
                    yield op, sites, mode, p / len(modes), value


def boundary(circuit, parsed, source, selected, rounds):
    records, ops, det, obs = parsed
    n = source["memory"]["n"]
    hx, hz = (source["memory"][key] for key in ("hx", "hz"))
    uid_order = list(range(n, n + len(hx) + len(hz)))
    k = len(obs)
    if k != len(source["logicals"]) or rounds < 2:
        raise ValueError("logical count or round count mismatch")
    by_meter = defaultdict(list)
    for i, rr in enumerate(records):
        if len(rr) == 1 and rr[0][1] in uid_order:
            by_meter[rr[0][1]].append(i)
    if any(len(by_meter[q]) != rounds for q in uid_order):
        raise ValueError("missing syndrome measurement")
    first = {tuple([by_meter[n + i][0]]) for i in selected["X"]}
    first.update((by_meter[n + len(hx) + i][0],) for i in selected["Z"])
    width = len(first)
    if len(det) < width or {tuple(v) for v in det[:width]} != first:
        raise ValueError("incorrect first-round postselection mask")
    parity = {tuple(sorted(by_meter[q][:2])) for q in uid_order}
    if {tuple(sorted(v)) for v in det[width : width + len(uid_order)]} != parity:
        raise ValueError(
            "postselection must compare all checks across the first two rounds"
        )
    terminal = len(records) - k - len(uid_order)

    def masks(rr):
        return tuple(sum(1 << q for a, q in rr if a in axes) for axes in ("XY", "ZY"))

    expected = (
        [(w, 0) for w in hx]
        + [(0, w) for w in hz]
        + [tuple(v[:2]) for v in source["logicals"]]
    )
    if [masks(v) for v in records[terminal:]] != expected or obs != [
        [len(records) - k + j] for j in range(k)
    ]:
        raise ValueError("closed measurement differs from the source")
    closed = [
        tuple(sorted((by_meter[q][-1], terminal + i))) for i, q in enumerate(uid_order)
    ]
    remaining = [tuple(sorted(v)) for v in det[width + len(uid_order) :]]
    later = [
        tuple(sorted(by_meter[q][j : j + 2]))
        for j in range(1, rounds - 1)
        for q in uid_order
    ]
    if sorted(remaining) != sorted(later + closed):
        raise ValueError("invalid decoder syndrome boundary")
    ideal = circuit.without_noise()
    circuit.detector_error_model()  # Reject stochastic detectors and observables.
    reference = ideal.reference_sample()
    for i, (xx, zz, sx, sz) in enumerate(source["logicals"]):
        if reference[-k + i] != (sx ^ sz ^ (((xx & zz).bit_count() - 1) // 2 % 2)):
            raise ValueError("wrong ideal signed logical Y")
    return width + len(uid_order)


def audit(circuit, source, selected, rounds, p_ref, ledger=None):
    parsed = parse(circuit)
    mask_size = boundary(circuit, parsed, source, selected, rounds)
    records, ops, det, obs = parsed
    checks = (1 << len(det)) - 1
    rejected = (1 << mask_size) - 1
    counts = Counter()
    coefficient = Counter()
    raw = {}
    carriers = source["preparation"]["carriers"]
    # Mark each explicit X reset-flip channel following the source's RY.
    pending = set()
    carrier_errors = set()
    for gate, qs, p, _, _, op in ops:
        if gate == "RY":
            pending.update(qs)
        elif gate == "X_ERROR":
            carrier_errors.update((op, q) for q in pending.intersection(qs))
            pending.difference_update(qs)
        elif gate in ("H", "CX", "R", "RX", "M", "MX", "MR"):
            pending.difference_update(qs)
    for op, sites, mode, p, v in responses(circuit, parsed):
        ds = v & checks
        ls = v >> len(det)
        counts["physical_modes"] += 1
        if not ds & rejected:
            counts["accepted_modes"] += 1
            if not ds and ls:
                coefficient[ls] += p / p_ref
            if ledger is not None:
                ledger(op, sites, mode, ds, ls, p)
        if len(sites) == 1 and (op, sites[0]) in carrier_errors and mode == "X":
            j = carriers.index(sites[0])
            if ds or ls != 1 << j:
                raise ValueError("raw preparation fault has wrong logical response")
            raw[j] = dict(operation=op, site=sites[0], mode=mode, mask=ls)
    return dict(
        postselect_width=mask_size,
        detectors=len(det),
        outputs=len(obs),
        **counts,
        first_order_silent_by_mask={
            str(i): str(coefficient[i]) for i in range(1, 1 << len(obs))
        },
        first_order_total=str(sum(coefficient.values())),
        first_order_correlated=str(
            sum(v for i, v in coefficient.items() if i.bit_count() > 1)
        ),
        explicit_raw_flip_witnesses=[raw[i] for i in sorted(raw)],
    )


def experiment_cost(circuit, mask_size):
    """Full-attempt two-qubit layer count and the independent DEM discard proxy."""
    product = 1.0
    for instruction in circuit.detector_error_model().flattened():
        if instruction.type == "error" and any(
            t.is_relative_detector_id() and t.val < mask_size
            for t in instruction.targets_copy()
        ):
            product *= 1 - instruction.args_copy()[0]
    # Respect printed clock ticks. Account for any intratick wire reuse.
    depth = 0
    wire_depth = {}
    for inst in circuit.flattened():
        if inst.name == "TICK":
            depth += max(wire_depth.values(), default=0)
            wire_depth.clear()
        elif inst.name == "CX":
            targets = inst.targets_copy()
            for c, t in zip(targets[::2], targets[1::2]):
                layer = 1 + max(wire_depth.get(c.value, 0), wire_depth.get(t.value, 0))
                wire_depth[c.value] = wire_depth[t.value] = layer
    return dict(
        physical_qubits=circuit.num_qubits,
        two_qubit_layers=depth + max(wire_depth.values(), default=0),
        dem_discard_union_proxy=round(1 - product, 10),
    )
