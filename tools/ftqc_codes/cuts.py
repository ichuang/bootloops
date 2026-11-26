"""Leg E: distance lower bounds from disjoint logical representatives (k = 1).

Let the code encode one qubit.  If s_1, ..., s_m are pairwise disjoint supports
of nontrivial Z-type logical operators (each s_i in ker H_X, none in
rowspace H_Z), then every nontrivial X-type logical operator has weight at
least m, so d_x >= m.  The mirror statement bounds d_z.

Proof.  For k = 1, ker H_X = rowspace H_Z + span(z0) for one vector z0.  The
overlap x . s of an X-type logical x (in ker H_Z) with s in ker H_X does not
change when s is shifted by rowspace H_Z, so x . s_i = x . z0 for every i.
If x . z0 were 0, x would be orthogonal to all of ker H_X and hence lie in
rowspace H_X, i.e. be a stabilizer.  So x has odd overlap with each of the m
disjoint sets and at least one coordinate in each.  (The counting step is
`weight_ge_of_cuts` in Lean-QEC, LeanQEC/LinearAlgebra/DistanceLB.lean.)

The sets are the certificate: they are stored in the evidence entry and
`certify.verify` re-checks them with dot products and one rank computation,
so this is the only leg whose lower bounds are re-verified in cheap mode.
The search for the sets is heuristic and untrusted.
"""
from math import isqrt

from . import certify, gf2


def _sides(code, q):
    """(checks the sets must satisfy, stabilizers they must avoid)."""
    if q == "d_x":
        return code.hx, code.hz     # Z-type representatives bound X-type weight
    if q == "d_z":
        return code.hz, code.hx
    raise ValueError(q)


def verify_cuts(code, q, sets):
    """Number of sets if they form a valid certificate, else raises Refusal."""
    n = code.n
    ex, ez = gf2.Echelon(code.hx), gf2.Echelon(code.hz)
    if n - ex.rank - ez.rank != 1:
        raise certify.Refusal("cut certificate needs k = 1")
    checks, stabs = _sides(code, q)
    stab_space = gf2.Echelon(stabs)
    used = 0
    for i, s in enumerate(sets):
        if any(j < 0 or j >= n for j in s) or len(set(s)) != len(s):
            raise certify.Refusal(f"cut {i}: bad support list")
        v = gf2.from_supp(s)
        if v == 0 or any(gf2.dot(v, h) for h in checks):
            raise certify.Refusal(f"cut {i}: not a logical operator of the required type")
        if v in stab_space:
            raise certify.Refusal(f"cut {i}: is a stabilizer")
        if v & used:
            raise certify.Refusal(f"cut {i}: overlaps an earlier cut")
        used |= v
    return len(sets)


def _grid_hints(n):
    d = isqrt(n)
    if d * d != n:
        return []
    rows = [[r * d + c for c in range(d)] for r in range(d)]
    cols = [[r * d + c for r in range(d)] for c in range(d)]
    return rows + cols


def find_cuts(code, q, hints=()):
    """Greedy search for disjoint representatives.  Heuristic and untrusted."""
    n = code.n
    checks, stabs = _sides(code, q)
    stab_space = gf2.Echelon(stabs)

    def ok(v):
        return v and not any(gf2.dot(v, h) for h in checks) and v not in stab_space
    chosen, used = [], 0
    cands = sorted({gf2.from_supp(s) for s in list(hints) + _grid_hints(n)}, key=gf2.wt)
    for v in cands:
        if ok(v) and not v & used:
            chosen.append(v)
            used |= v
    while True:
        # logical operators supported on the unused coordinates
        blocked = [1 << j for j in range(n) if (used >> j) & 1]
        basis = gf2.nullspace(list(checks) + blocked, n)
        pool = [b for b in basis if ok(b)]
        if not pool:
            break
        best = min(pool, key=gf2.wt)
        # cheap lightening: add stabilizers or trivial kernel vectors that reduce weight
        improved = True
        reducers = [b for b in basis if b != best]
        while improved:
            improved = False
            for r in reducers:
                t = best ^ r
                if gf2.wt(t) < gf2.wt(best) and ok(t):
                    best, improved = t, True
        chosen.append(best)
        used |= best
    return [gf2.supp(v) for v in chosen]


def run_leg(receipt, code, hints=None):
    hints = hints or {}
    for q in ("d_x", "d_z"):
        ev = dict(leg="E", quantity=q, method="disjoint logical representatives (k = 1), "
                                              "re-checked by certify.verify")
        if receipt["facts"]["k"] != 1:
            ev.update(status="skipped", reason="the cut argument is implemented for k = 1 only")
            receipt["evidence"].append(ev)
            continue
        sets = find_cuts(code, q, hints.get(q, ()))
        if not sets:
            ev.update(status="skipped", reason="no representative found")
            receipt["evidence"].append(ev)
            continue
        m = verify_cuts(code, q, sets)
        ev.update(status="proved", claim=f"{q} >= {m}", cuts=sets)
        certify.add_lower_bound(receipt, q, m, ev)
    return receipt
