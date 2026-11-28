"""Leg A: exact minimum-weight bounds by information-set enumeration.

Problem.  C is the row space of `gens` (dimension k) in GF(2)^n and S is a
subspace of C (possibly {0}).  Bound  m(C, S) = min{ wt(c) : c in C, c not in S }.

Method (Brouwer-Zimmermann).  Pick pairwise disjoint column sets P_1, P_2, ...
For each P_j build a basis B_j of C consisting of r_j "pivot rows" that
restrict to the identity on P_j, and k - r_j "extra rows" that vanish on P_j
(r_j is the rank of C restricted to the columns the set was drawn from).  A
codeword that is the sum of a pivot rows and b extra rows of B_j has exactly a
ones inside P_j.  Enumerate every sum of at most w rows of every B_j.  A
codeword that was not enumerated uses at least w + 1 rows of each B_j, of
which at most k - r_j are extra rows, so it has at least w + 1 - (k - r_j)
ones inside P_j, and since the P_j are disjoint its weight is at least

    L(w) = sum_j max(0, w + 1 - (k - r_j)).

Enumerated codewords are examined one by one and those in S are discarded.
Hence m(C, S) >= min(L(w), least weight seen outside S), and the bound is
exact as soon as L(w) reaches the least weight seen.

The result is a pair (lower, upper) with a witness for the upper bound.  The
bound is exact when lower == upper.  Everything is integer arithmetic.
"""
from . import gf2


def _info_sets(gens, n, max_sets=None):
    """Disjoint column sets with full k-row bases adapted to each."""
    k = len(gens)
    free = list(range(n))
    sets = []
    while free:
        rows, piv = gf2.systematic(gens, n, free)
        r = len(piv)
        if r == 0:
            break
        # complete to a basis of C with k - r rows that vanish on the pivots
        if r < k:
            ech = gf2.Echelon(rows)
            extra = []
            for g in gens:
                v = g
                for row, p in zip(rows, piv):
                    if (v >> p) & 1:
                        v ^= row
                if ech.add(v):
                    extra.append(v)
            rows = rows + extra
            assert len(rows) == k
        sets.append((rows, piv, r))
        pset = set(piv)
        free = [c for c in free if c not in pset]
        if max_sets and len(sets) >= max_sets:
            break
    return sets


def _combos_at_weight(rows, w, visit):
    """Call visit(v) for the sum of every w-subset of rows."""
    k = len(rows)

    def rec(start, left, acc):
        if left == 0:
            visit(acc)
            return
        for i in range(start, k - left + 1):
            rec(i + 1, left - 1, acc ^ rows[i])

    rec(0, w, 0)


def _ncomb(k, w):
    from math import comb
    return comb(k, w)


def min_weight(gens, n, exclude=(), target=None, max_combos=5_000_000, max_sets=None):
    """Bound min{wt(c): c in span(gens), c not in span(exclude)}.

    target: stop as soon as the lower bound reaches it (a claimed distance);
            None runs until the bound is exact.
    Returns a dict: lower, upper (None if no such word exists or none found),
    witness (int or None), exact (bool), depth, combos, info_sets [(size, rank)].
    If more than max_combos combinations would be needed the search stops and
    returns the bound proven so far, with the reason in "budget".
    """
    ech = gf2.Echelon(gens)
    basis = list(ech.basis)
    k = len(basis)
    S = gf2.Echelon(exclude)
    for b in S.basis:
        if b not in ech:
            raise ValueError("exclude is not a subspace of span(gens)")
    if k == S.rank:
        return dict(lower=None, upper=None, witness=None, exact=True, depth=0,
                    combos=0, info_sets=[], empty=True)
    sets = _info_sets(basis, n, max_sets)
    best = [None, None]  # weight, word
    combos = 0

    def visit(v):
        w = v.bit_count()
        if best[0] is None or w < best[0]:
            if v and v not in S:
                best[0], best[1] = w, v

    def L(w):
        return sum(max(0, w + 1 - (k - r)) for _, _, r in sets)

    lower, depth, budget = 1, 0, None
    # every word outside S is nonzero, so 1 is a valid starting bound
    for w in range(1, k + 1):
        need = _ncomb(k, w) * len(sets)
        if combos + need > max_combos:
            # stop with what is proven: unseen words have weight >= L(depth)
            if best[0] is not None:
                lower = min(lower, best[0])
            budget = (f"stopped before depth {w}: it needs {need} more combinations "
                      f"(cap {max_combos})")
            break
        for rows, _, _ in sets:
            _combos_at_weight(rows, w, visit)
        combos += need
        depth = w
        lower = max(lower, L(w))
        if best[0] is not None and lower >= best[0]:
            lower = best[0]
            break
        if target is not None and lower >= target:
            break
    else:
        # all 2^k - 1 nonzero words were enumerated in the first set
        lower = best[0]
    return dict(lower=lower, upper=best[0], witness=best[1],
                exact=(best[0] is not None and lower == best[0]),
                depth=depth, combos=combos,
                info_sets=[[len(p), r] for _, p, r in sets], empty=False, budget=budget)


def brute_min_weight(gens, n, exclude=()):
    """Reference implementation by full enumeration (k <= 22).  For tests."""
    basis = gf2.Echelon(gens).basis
    k = len(basis)
    if k > 22:
        raise ValueError("brute force limited to dimension 22")
    S = gf2.Echelon(exclude)
    best, word, v = None, None, 0
    for i in range(1, 1 << k):
        # Gray code step
        v ^= basis[(i & -i).bit_length() - 1]
        w = v.bit_count()
        if (best is None or w < best) and v not in S:
            best, word = w, v
    return best, word
