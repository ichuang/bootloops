"""Questions about "at most j faults", answered on the columns of the model.

circuit_distance   least number of faults that fire no detector and flip a logical
                   observable (memory and syndrome-extraction experiments)
strict_ft          for state preparations with postselection: every accepted set of
                   j <= s faults leaves the output within Pauli weight j of the
                   intended state
undetected_harm    least number of faults that fire no detector and leave a harmful
                   residual (the "fault distance" of a preparation)

Leg A is exact enumeration (meet in the middle on the detector signature).
Leg B hands the same statement to SAT through qeccert.cnf and has the LRAT proof
re-checked by qeccert.lrat.
"""
from itertools import combinations


class Budget(Exception):
    pass


def _sums(vecs, j, cap):
    """Yield (xor of the tuples, index tuple) over all j-subsets of vecs."""
    from math import comb
    if comb(len(vecs), j) > cap:
        raise Budget(f"{comb(len(vecs), j)} subsets of size {j} exceed the cap {cap}")
    k = len(vecs[0]) if vecs else 0
    for idx in combinations(range(len(vecs)), j):
        acc = [0] * k
        for i in idx:
            v = vecs[i]
            for t in range(k):
                acc[t] ^= v[t]
        yield tuple(acc), idx


def min_undetected(cols, harmful, max_weight, cap=3_000_000):
    """Least j <= max_weight such that some j columns sum to (det = 0, payload with
    harmful(payload) true).  cols: list of (det, payload...) tuples.
    Returns dict(lower, upper, witness (column indices) or None, exact, note)."""
    lower, note = 1, None
    if not cols:
        return dict(lower=None, upper=None, witness=None, exact=True, note="no fault has any effect")
    tables = {0: {tuple([0] * len(cols[0])): ()}}

    def table(j):
        if j not in tables:
            t = {}
            for v, idx in _sums(cols, j, cap):
                t.setdefault(v, idx)
            tables[j] = t
        return tables[j]

    for w in range(1, max_weight + 1):
        a = w // 2
        b = w - a
        try:
            ta, tb = table(a), table(b)
        except Budget as e:
            note = f"stopped before weight {w}: {e}"
            break
        # index the smaller half by detector signature
        by_det = {}
        for v, idx in ta.items():
            by_det.setdefault(v[0], []).append((v, idx))
        for v, idx in tb.items():
            for u, jdx in by_det.get(v[0], ()):
                tot = tuple(p ^ q for p, q in zip(u, v))
                if harmful(tot[1:]):
                    wit = sorted(set(idx) ^ set(jdx))
                    return dict(lower=len(wit), upper=len(wit), witness=wit, exact=True, note=None)
        lower = w + 1
    return dict(lower=lower, upper=None, witness=None, exact=False,
                note=note or f"no harmful undetected fault set of weight <= {max_weight}")


def weight_table(model, j):
    """Syndromes of all Paulis of weight <= j on the output qubits."""
    outs = model.circuit.outputs
    singles = []
    for q in outs:
        b = 1 << q
        singles.append([model.syndrome(b, 0), model.syndrome(0, b), model.syndrome(b, b)])
    seen = {0}
    for r in range(1, j + 1):
        for qs in combinations(range(len(outs)), r):
            def rec(i, acc):
                if i == r:
                    seen.add(acc)
                    return
                for s in singles[qs[i]]:
                    rec(i + 1, acc ^ s)
            rec(0, 0)
    return seen


def _accepted_sets(cols, j, cap):
    """Yield index tuples of the j-subsets of cols (pairs (det, syn)) whose detector
    signatures cancel, without enumerating the other subsets: choose j - 1 columns,
    then look the last one up by its required detector signature."""
    from math import comb
    by_det = {}
    for i, (d, s) in enumerate(cols):
        by_det.setdefault(d, []).append(i)
    if j == 1:
        for i in by_det.get(0, ()):
            yield (i,)
        return
    if comb(len(cols), j - 1) > cap:
        raise Budget(f"{comb(len(cols), j - 1)} subsets of size {j - 1} exceed the cap {cap}")
    for idx in combinations(range(len(cols)), j - 1):
        d = 0
        for i in idx:
            d ^= cols[i][0]
        for last in by_det.get(d, ()):
            if last > idx[-1]:
                yield idx + (last,)


def strict_ft(model, order, cap=3_000_000):
    """Check, for j = 1..order: every set of j faults firing no detector leaves the
    output within Pauli weight j of the believed state.  Returns per-order counts and
    the first violation (lowest j) with its columns."""
    index = {}
    for i, (d, o, s) in enumerate(model.columns):
        if d or s:
            index.setdefault((d, s), i)
    cols = sorted(index)
    per, violation, note = [], None, None
    for j in range(1, order + 1):
        ok_syn = weight_table(model, j)
        accepted = bad = 0
        try:
            for idx in _accepted_sets(cols, j, cap):
                s = 0
                for i in idx:
                    s ^= cols[i][1]
                accepted += 1
                if s not in ok_syn:
                    bad += 1
                    if violation is None:
                        violation = dict(order=j, columns=[index[cols[i]] for i in idx])
        except Budget as e:
            note = f"stopped at order {j}: {e}"
            break
        per.append(dict(order=j, accepted_sets=accepted, violations=bad))
        if violation:
            break
    holds = max([p["order"] for p in per if p["violations"] == 0], default=0)
    if violation:
        holds = violation["order"] - 1
    return dict(holds_through=holds, per_order=per, violation=violation, note=note)


def min_heavier(model, weight, max_faults, cap=3_000_000):
    """Least number of faults that fire no detector and leave the output further than
    Pauli weight `weight` from the believed state ("an accepted error heavier than
    `weight`").  Same return shape as min_undetected, witness in model columns."""
    ok = weight_table(model, weight)
    index = {}
    for i, (d, o, s) in enumerate(model.columns):
        if d or s:
            index.setdefault((d, s), i)
    cols = sorted(index)
    r = min_undetected(cols, lambda p: p[0] not in ok, max_faults, cap)
    if r["witness"] is not None:
        r["witness"] = [index[cols[i]] for i in r["witness"]]
    return r


def census(model):
    """Per location, the worst thing any of its Paulis leaves on the output, ignoring
    detection: 0 harmless, 1 equivalent to a single-qubit error, 2 heavier ("hook")."""
    one = {}
    for q in model.circuit.outputs:
        b = 1 << q
        for x, z in ((b, 0), (0, b), (b, b)):
            one.setdefault(model.syndrome(x, z), set()).add(q)
    hook, local, harmless = [], {q: 0 for q in model.circuit.outputs}, 0
    for loc in model.locations:
        worst, where = 0, set()
        for eff in loc.faults.values():
            s = model.syndrome(*eff)
            if s == 0:
                continue
            if s in one:
                worst = max(worst, 1)
                where |= one[s]
            else:
                worst = 2
        if worst == 0:
            harmless += 1
        elif worst == 2:
            hook.append(loc.index)
        else:
            for q in where:
                local[q] += 1
    return dict(locations=len(model.locations), harmless=harmless, hook_capable=len(hook),
                hook_locations=hook, single_qubit_only=len(model.locations) - harmless - len(hook),
                max_local_per_output=max(local.values(), default=0),
                local_per_output=[local[q] for q in model.circuit.outputs])
