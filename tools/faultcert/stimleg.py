"""Leg S: Stim as an independent oracle for the fault model.

The circuit is exported to Stim with the same noise (DEPOLARIZE1 after one-qubit
operations, DEPOLARIZE2 after CX, a flip probability on measurements) and the same
detectors and observables.  Stim then
  * refuses the circuit if a declared detector or observable is not deterministic,
    which checks this package's stabilizer simulation;
  * derives its own detector error model, whose set of distinct symptoms must equal
    the set of (detector, observable) signatures of this package's fault columns;
  * searches for an undetectable logical error; what it finds is an upper bound on the
    circuit distance and must not be below a lower bound proven here.
Stim (C. Gidney, Quantum 5, 497, 2021) is not trusted for any bound in a receipt; a
disagreement is reported as a refusal.
"""
P = 0.001


def available():
    try:
        import stim  # noqa: F401
        return None
    except ImportError:
        return "stim is not importable in this interpreter (pip install stim)"


def to_stim(circuit, model, noisy=True):
    lines = []
    count = 0
    for tick in circuit.ticks:
        after = []
        for op in tick:
            q = " ".join(map(str, op.qs))
            nz = noisy and op.noisy
            if op.name in ("M", "MX"):
                count += 1
                lines.append(f"{op.name}({P}) {q}" if nz else f"{op.name} {q}")
            else:
                if op.name == "Z_IF":
                    lines.append(f"CZ rec[{op.record - count}] {q}")
                else:
                    gate = "S_DAG" if op.name == "SDG" else op.name
                    lines.append(f"{gate} {q}")
                if nz:
                    after.append(("DEPOLARIZE2" if op.name == "CX" else "DEPOLARIZE1") + f"({P}) {q}")
        lines += after
        lines.append("TICK")
    nm = len(circuit.meas)
    for d in model.detectors:
        lines.append("DETECTOR " + " ".join(f"rec[{i - nm}]" for i in d))
    for k, o in enumerate(circuit.observables):
        lines.append(f"OBSERVABLE_INCLUDE({k}) " + " ".join(f"rec[{i - nm}]" for i in o))
    return "\n".join(lines) + "\n"


def compare(circuit, model, distance_lower=None, search_logical=True):
    """Returns an evidence dict; status 'agrees', 'DISAGREES' or 'skipped'."""
    why = available()
    ev = dict(leg="S", method="Stim detector error model and logical-error search")
    if why:
        ev.update(status="skipped", reason=why)
        return ev
    import stim
    ev["stim"] = stim.__version__
    try:
        sc = stim.Circuit(to_stim(circuit, model))
        dem = sc.detector_error_model(decompose_errors=False, allow_gauge_detectors=False)
    except ValueError as e:
        ev.update(status="DISAGREES", reason=f"Stim rejects the circuit: {str(e)[:300]}")
        return ev
    theirs = set()
    for inst in dem.flattened():
        if inst.type != "error":
            continue
        det = obs = 0
        for t in inst.targets_copy():
            if t.is_relative_detector_id():
                det ^= 1 << t.val
            elif t.is_logical_observable_id():
                obs ^= 1 << t.val
        if det or obs:
            theirs.add((det, obs))
    ours = {(d, o) for d, o, s in model.columns if d or o}
    ev.update(mechanisms_ours=len(ours), mechanisms_stim=len(theirs))
    if ours != theirs:
        ev.update(status="DISAGREES",
                  reason=f"{len(ours - theirs)} signatures only here, {len(theirs - ours)} only in Stim")
        return ev
    if circuit.observables and search_logical:
        try:
            found = sc.search_for_undetectable_logical_errors(
                dont_explore_detection_event_sets_with_size_above=6,
                dont_explore_edges_with_degree_above=6,
                dont_explore_edges_increasing_symptom_degree=False)
            ev["stim_logical_error_weight"] = len(found)
            if distance_lower is not None and len(found) < distance_lower:
                ev.update(status="DISAGREES",
                          reason=f"Stim found a logical error of weight {len(found)}, below the "
                                 f"bound {distance_lower} proven here")
                return ev
        except ValueError as e:
            ev["stim_search"] = f"no logical error found by Stim's search: {str(e)[:120]}"
    ev["status"] = "agrees"
    return ev
