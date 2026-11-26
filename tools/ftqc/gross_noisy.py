"""Circuit distance of the finite noisy gross Z operation under CSS sector probes."""

from faultcert.record_model import RecordModel
from faultcert.projected import exclude

from .algebra import bits
from .gross import shifted_z
from .gross_cycle import PORT, cycle
from .gross_probes import experiment


def projection(model, basis, rounds):
    """Keep terminal syndrome, split consistency and protected terminal logicals.
    Also check the clean-round readout lemma on every elementary signature."""
    nz = 72
    nport = 12
    nterminal = nz if basis == "Z" else nz + 7
    terminal_shift = model.nd - nterminal
    phase_mask = (1 << nz) - 1
    first_round = model.locations[0][0] + 1
    length = len(cycle())
    columns = set()
    by_round = {}
    for (tick, _, _), group in zip(model.locations, model.groups, strict=True):
        epoch = (tick - first_round) // length
        values = set(group.values())
        for v in values:
            end = v >> terminal_shift
            if basis == "Z":
                # Observable zero is the port product, one is measured Z on
                # the final data, followed by the rest of the Z logical basis.
                end = (end & phase_mask) | ((end >> (nz + 1)) << nz)
            if end:
                columns.add(end)
        if basis == "Z":
            by_round.setdefault(epoch, set()).update(values)
    if basis == "Z":
        for epoch, values in by_round.items():
            for v in values:
                synd = v & phase_mask
                read = (v >> model.nd) & 1
                final_syn = (v >> terminal_shift) & phase_mask
                final_z = (v >> (model.nd + 1)) & 1
                for t in range(rounds):
                    if t != epoch and (synd, read) != (
                        (final_syn, final_z) if epoch < t else (0, 0)
                    ):
                        raise ValueError(
                            "fault violates clean-round syndrome/readout identity"
                        )
                    if t + 1 < rounds:
                        diff = v >> (nz + t * (nz + nport))
                        synd ^= diff & phase_mask
                        read ^= ((diff >> nz) & ((1 << nport) - 1)).bit_count() & 1
    return sorted(columns), nterminal


def time_witness(model):
    return [
        [i, "flip"]
        for i, (_, _, op) in enumerate(model.locations)
        if op.name == "M" and op.qs == (PORT,)
    ]


def report(spec, budget):
    rounds = spec["rounds"]
    models = {b: RecordModel(experiment(rounds, b)) for b in "XZ"}

    # Noisy sites correspond one-for-one; the ideal probes have different encoders.
    def physical(m):
        offset = m.locations[0][0]
        return [(t - offset, op.text()) for t, _, op in m.locations]

    if physical(models["X"]) != physical(models["Z"]):
        raise ValueError("basis probes use different physical fault locations")
    m = models["X"]
    split = m.locations[-1][0]
    sites = set(bits(shifted_z(0, 1)))
    spatial = [
        [i, "Z"]
        for i, (t, _, op) in enumerate(m.locations)
        if t == split and op.name == "I" and op.qs[0] in sites
    ]
    candidates = [time_witness(models["Z"]), spatial]
    bank = []
    for faults in candidates:
        signatures = {b: m.replay(faults) for b, m in models.items()}
        if any(d for d, _ in signatures.values()) or not any(
            o for _, o in signatures.values()
        ):
            raise ValueError("upper witness does not replay in both physical probes")
        bank.append(
            dict(
                faults=faults,
                flipped_tests={b: list(bits(o)) for b, (_, o) in signatures.items()},
            )
        )
    best = min(bank, key=lambda w: len(w["faults"]))
    upper = len(best["faults"])
    details = {}
    for basis, m in models.items():
        columns, nd = projection(m, basis, rounds)
        lower = 1
        evidence = []
        # Per-probe budget bounds recursion; a completed projection hit is
        # never reported as a physical counterexample.
        remaining = budget
        # Increasing order preserves proved intermediate bounds on timeout.
        for w in range(1, upper):
            ok, nodes, hit = exclude(columns, nd, w, remaining)
            evidence.append(
                dict(
                    weight=w,
                    nodes=nodes,
                    outcome="excluded"
                    if ok
                    else "projected witness"
                    if hit
                    else "budget",
                )
            )
            remaining -= nodes
            if not ok:
                break
            lower = w + 1
        details[basis] = dict(
            circuit_sha256=m.circuit.sha256(),
            locations=len(m.locations),
            detector_count=m.nd,
            logical_tests=len(m.circuit.observables),
            signatures=len(m.columns),
            projected_columns=len(columns),
            spatial_lower=lower,
            search=evidence,
        )
    lower = min(rounds, *(d["spatial_lower"] for d in details.values()))
    return dict(
        status="proved"
        if lower >= spec["target"]
        else "counterexample"
        if upper < spec["target"]
        else "inconclusive",
        model="gross-z-location-v1",
        rounds=rounds,
        lower=lower,
        upper=upper,
        exact=lower == upper,
        probes=details,
        witness=best,
    )
