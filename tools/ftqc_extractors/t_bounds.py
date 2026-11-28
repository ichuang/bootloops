"""Structural distance bounds for verified cat meters on equivalent HGP ports.

The proof and trust boundary are specified in T_PROOFS.md. This checker accepts
only the compiled meter topology and perfect terminal syndrome boundaries.
"""

from collections import Counter
from .algebra import bits, hgp_parameters, seed
from .t_faults import FaultCircuit
from .magic import check_readout_classifier
from .proof import verify_witness


def inspect_round(m, part, repetitions):
    start, end, x, z, record, flags = part
    sites = list(bits(x | z))
    cat = list(range(2 * m.n, 2 * m.n + len(sites)))
    meter = cat[-1] + 1
    # Independent parse: no live data interaction before the certified sweep.
    pos = start

    def take(op):
        nonlocal pos
        if pos >= end or m.ops[pos] != op:
            raise ValueError("unrecognized cat topology in structural certificate")
        pos += 1

    take(("RX", cat[0]))
    for u, v in zip(cat, cat[1:]):
        take(("RZ", v))
        take(("CX", u, v))
    checked = []
    for _ in range(repetitions):
        for u, v in zip(cat, cat[1:]):
            take(("RZ", meter))
            take(("CX", u, meter))
            take(("CX", v, meter))
            checked.append(pos)
            take(("MZ", meter))
    if checked != flags:
        raise ValueError("unbound cat flag")
    for u, q in zip(cat, sites):
        if (x & z) >> q & 1:
            take(("SDG", q))
            take(("CX", u, q))
            take(("S", q))
        elif x >> q & 1:
            take(("CX", u, q))
        else:
            take(("H", q))
            take(("CX", u, q))
            take(("H", q))
    reads = []
    for u in cat:
        reads.append(pos)
        take(("MX", u))
    if pos != end or sum(1 << t for t in reads) != record:
        raise ValueError("incomplete cat readout")
    return sites, reads[0]


def structural(s):
    if s["resource"] != "ideal-encoded-T":
        raise ValueError("structural bounds require a supplied ideal resource")
    check_readout_classifier()
    rounds = s["rounds"]
    verification = s.get("verification_rounds", 1)
    base = hgp_parameters(*seed(s["a"]), *seed(s["b"]))
    cases = {}
    for name, arm in [("gate", None), ("read_X", 0), ("read_Y", 1)]:
        circuit = FaultCircuit(s, arm)
        stage_bounds = []
        read_faults = []
        offset = 0
        expected_det = []
        for stage in range(2 if arm is None else 3):
            supports = Counter()
            first = None
            for part in circuit.parts[stage * rounds : (stage + 1) * rounds]:
                if part[0] != offset:
                    raise ValueError("unexpected interstage operation")
                offset = part[1]
                sites, read = inspect_round(circuit, part, verification)
                expected_det.extend((0, 0, 1 << f) for f in part[-1])
                if first is None:
                    first = part[4]
                else:
                    expected_det.append((0, 0, first ^ part[4]))
                supports.update(sites)
                if stage == 0:
                    read_faults.append([read, 1])
            width = max(supports.values())
            stage_bounds.append(
                dict(
                    rounds=rounds,
                    max_incidence=width,
                    lower=(rounds + width - 1) // width,
                )
            )
        if (
            offset != len(circuit.ops)
            or circuit.det[: len(expected_det)] != expected_det
        ):
            raise ValueError("unverified operations or checks")
        # With fewer than V+1 faults either one full verification sweep is clean,
        # or (at V faults) every sweep is faulty and cat preparation is perfect.
        # In the latter case all faults are after cat preparation and each spreads
        # onto at most one data site. The clean-sweep argument also covers V=1.
        lower = min(base["d"], verification + 1, *[r["lower"] for r in stage_bounds])
        verify_witness(circuit, rounds + 1, read_faults)
        cases[name] = dict(
            lower=lower,
            upper=rounds,
            upper_witness=read_faults,
            exact=(lower == rounds),
            stages=stage_bounds,
            locations=len(circuit.ops),
        )
    return dict(
        status="proved",
        method="verified-cat spread and representative coverage",
        resource=s["resource"],
        base=base,
        verification_sweeps=verification,
        cases=cases,
    )
