"""Two-sided logical-error enclosures from a linked-cluster expansion, for localized decoders.

A decoder is *localized at detector range r* when it splits the fired detectors into groups,
two detectors being grouped when their distance in the detector graph (detectors adjacent
when one mechanism flips both) is at most r, decodes every group on its own with a base
decoder, and adds the predictions.  `decode_localized` wraps any base decoder this way.

Additivity.  Call two mechanisms adjacent when they share a detector, and let G^R be the
R-th power of that graph with R = r + 1.  If a fault set F has components C_1, ..., C_m in
G^R, then the residual of F under a localized decoder is the sum of the residuals of the
C_j.  Proof: a detector fired by C_i and one fired by C_j at detector distance k give
mechanisms of C_i and C_j at mechanism distance at most k + 1, so for i != j the fired
detectors are farther apart than r and are never grouped together; the groups of H(F) are
the groups of the H(C_j), and the base decoder is a function of the group alone.

Enclosure.  Let lambda = E[(-1)^(kappa . residual)] under independent faults of probability
p, activity x = p / (1 - p), M mechanisms, and D the maximal degree of G^R.  The numerator
and denominator of lambda are hard-core polymer gases of G^R-connected fault sets, the first
with the sign of each set as activity sign.  By the Kotecky-Preiss criterion (Commun. Math.
Phys. 103, 491, 1986) with the count (e D)^(n-1) of connected sets of size n through a
vertex, both logarithms converge for x <= x_R = 1 / (e (D + 1) + e^2 D), and with
theta = x / x_R <= 1 the clusters of total size at least s contribute at most M theta^s to
each.  Clusters without a failing set cancel.  Hence, if every connected set of size at
most t is corrected and S_m is the sum of x^|C| over failing connected sets of size at
most m,

    | log lambda + 2 S_m |  <=  2 M theta^min(m + 1, t + 2).

The parity-flip probability is (1 - lambda) / 2.  With y_lo = 2 S_m - remainder and
y_hi = 2 S_m + remainder, the bounds 1 - exp(-y) <= y and 1 - exp(-y) >= y - y^2 / 2 give
the enclosure.  The connected sets through size m are enumerated and decoded, so t and S_m
are computed, not assumed.

Arithmetic is exact in `fractions.Fraction`; printed bounds are rounded outward.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Iterator

from .bp_osd import DecodeResult
from .decimal_bounds import RefusalError, fraction_to_outward_strs, parse_exact_probability
from .dem import DecodingProblem, bit_positions

EULER_E_UPPER = Fraction(27183, 10000)


def _detector_neighbors(problem: DecodingProblem) -> list[set[int]]:
    nbrs: list[set[int]] = [set() for _ in range(problem.num_checks)]
    for col in problem.mech_detectors:
        dets = bit_positions(col)
        for d in dets:
            nbrs[d].update(dets)
    for d, s in enumerate(nbrs):
        s.discard(d)
    return nbrs


def _mechanism_power_graph(problem: DecodingProblem, radius: int) -> list[frozenset[int]]:
    """Adjacency of G^radius: mechanisms at distance 1..radius, adjacent when sharing a detector."""
    members = [bit_positions(r) for r in problem.check_rows]
    base: list[set[int]] = []
    for j in range(problem.num_mechanisms):
        s: set[int] = set()
        for c in bit_positions(problem.mech_detectors[j]):
            s.update(members[c])
        s.discard(j)
        base.append(s)
    out = []
    for j in range(problem.num_mechanisms):
        seen = {j}
        frontier = [j]
        for _ in range(radius):
            nxt = []
            for v in frontier:
                for w in base[v]:
                    if w not in seen:
                        seen.add(w)
                        nxt.append(w)
            frontier = nxt
        seen.discard(j)
        out.append(frozenset(seen))
    return out


def syndrome_groups(problem: DecodingProblem, syndrome: int, detector_range: int,
                    _nbrs: list[set[int]] | None = None) -> list[int]:
    """Split a syndrome into groups of fired detectors linked at detector distance <= range."""
    if detector_range < 0:
        raise RefusalError("INVALID-DETECTOR-RANGE", "detector_range must be >= 0")
    nbrs = _nbrs if _nbrs is not None else _detector_neighbors(problem)
    fired = set(bit_positions(syndrome))
    groups = []
    while fired:
        root = min(fired)
        fired.discard(root)
        group = 1 << root
        todo = [root]
        while todo:
            d = todo.pop()
            ball, frontier = {d}, [d]
            for _ in range(detector_range):
                nxt = []
                for v in frontier:
                    for w in nbrs[v]:
                        if w not in ball:
                            ball.add(w)
                            nxt.append(w)
                frontier = nxt
            for w in ball & fired:
                fired.discard(w)
                group |= 1 << w
                todo.append(w)
        groups.append(group)
    return groups


def decode_localized(problem: DecodingProblem, syndrome: int, base: str | dict[str, Any] = "union_find",
                     detector_range: int = 1) -> DecodeResult:
    """Decode every group of fired detectors on its own with `base` and add the results."""
    from .effective_distance import run_decoder

    groups = syndrome_groups(problem, syndrome, detector_range)
    error = obs = 0
    converged = True
    for g in groups:
        r = run_decoder(problem, g, base)
        if r.converged and problem.syndrome_of(r.error_mask) == g:
            error ^= r.error_mask
            obs ^= r.observable_mask
        else:
            converged = False
    return DecodeResult(
        decoder="localized",
        error_mask=error,
        converged=converged,
        observable_mask=obs,
        hamming_weight=error.bit_count(),
        soft_weight=problem.soft_weight_of(error),
        stats={"groups": len(groups), "detector_range": detector_range,
               "base": base if isinstance(base, str) else base.get("algorithm")},
    )


def connected_sets(adj: list[frozenset[int]], max_size: int) -> Iterator[tuple[int, ...]]:
    """Every connected vertex set of size <= max_size, each exactly once."""

    def extend(sub: tuple[int, ...], ext: set[int], closed: set[int], root: int):
        yield sub
        if len(sub) == max_size:
            return
        ext = set(ext)
        while ext:
            w = ext.pop()
            new = {u for u in adj[w] if u > root and u not in closed}
            yield from extend(sub + (w,), ext | new, closed | new, root)

    for v in range(len(adj)):
        first = {u for u in adj[v] if u > v}
        yield from extend((v,), first, first | {v}, v)


class _Signs:
    def __init__(self, problem, base, detector_range, kappa):
        self.problem, self.base, self.r, self.kappa = problem, base, detector_range, kappa
        self.cache: dict[int, int] = {0: 0}

    def sign(self, faults) -> int:
        s = o = 0
        for i in faults:
            s ^= self.problem.mech_detectors[i]
            o ^= self.problem.mech_observables[i]
        if s not in self.cache:
            self.cache[s] = decode_localized(self.problem, s, self.base, self.r).observable_mask
        return -1 if ((o ^ self.cache[s]) & self.kappa).bit_count() & 1 else 1


def additivity_counterexample(problem: DecodingProblem, sign, radius: int, max_weight: int,
                              budget_cap: int = 200_000):
    """Search fault sets of weight <= max_weight for one whose sign is not the product of the
    signs of its components in G^radius.  `sign` maps a tuple of mechanisms to +1 or -1.
    Returns None or a dict with the set and its components."""
    from itertools import combinations

    adj = _mechanism_power_graph(problem, radius)
    seen = 0
    for w in range(2, max_weight + 1):
        for F in combinations(range(problem.num_mechanisms), w):
            seen += 1
            if seen > budget_cap:
                raise RefusalError("ADDITIVITY-BUDGET-EXCEEDED", f"more than {budget_cap} fault sets")
            left, comps = set(F), []
            while left:
                comp, todo = [], [left.pop()]
                while todo:
                    v = todo.pop()
                    comp.append(v)
                    for u in list(left):
                        if u in adj[v]:
                            left.discard(u)
                            todo.append(u)
                comps.append(tuple(sorted(comp)))
            if len(comps) < 2:
                continue
            prod = 1
            for c in comps:
                prod *= sign(c)
            if prod != sign(F):
                return {"faults": list(F), "components": [list(c) for c in sorted(comps)]}
    return None


def enclose_linked_cluster(
    problem: DecodingProblem,
    base_decoder: str | dict[str, Any],
    physical_p: str | int | Fraction,
    *,
    detector_range: int = 1,
    max_cluster_size: int = 2,
    kappa: int = 1,
    budget_cap: int = 200_000,
    digits: int = 30,
) -> dict[str, Any]:
    """Two-sided enclosure of the probability that the observable parity `kappa` is flipped
    after decoding with `base_decoder` localized at `detector_range`, under independent faults
    of probability `physical_p` on every mechanism."""
    p = parse_exact_probability(physical_p, allow_zero=True, allow_one=False, max_half=True)
    if max_cluster_size < 1:
        raise RefusalError("INVALID-CLUSTER-SIZE", "max_cluster_size must be >= 1")
    if kappa <= 0 or kappa >> problem.num_observables:
        raise RefusalError("INVALID-OBSERVABLE-MASK", "kappa must select at least one observable")
    radius = detector_range + 1
    adj = _mechanism_power_graph(problem, radius)
    M = problem.num_mechanisms
    degree = max(1, max((len(a) for a in adj), default=1))
    x_radius = 1 / (EULER_E_UPPER * (degree + 1) + EULER_E_UPPER ** 2 * degree)
    x = p / (1 - p)
    theta = x / x_radius

    signs = _Signs(problem, base_decoder, detector_range, kappa)
    failing = [0] * (max_cluster_size + 1)
    total = 0
    for C in connected_sets(adj, max_cluster_size):
        total += 1
        if total > budget_cap:
            raise RefusalError("CLUSTER-BUDGET-EXCEEDED",
                               f"more than {budget_cap} connected fault sets through size {max_cluster_size}")
        if signs.sign(C) < 0:
            failing[len(C)] += 1
    first = next((s for s in range(1, max_cluster_size + 1) if failing[s]), None)
    t = (first - 1) if first is not None else max_cluster_size

    out: dict[str, Any] = {
        "mode": "linked_cluster",
        "evidence_kind": "analytical_linked_cluster_enclosure",
        "decoder": {"algorithm": "localized", "base": base_decoder, "detector_range": detector_range},
        "observable_mask": kappa,
        "physical_p": f"{p.numerator}/{p.denominator}",
        "additivity": f"proved for a localized decoder at mechanism range {radius} (module docstring)",
        "mechanism_range": radius,
        "interaction_degree": degree,
        "connected_sets_decoded": total,
        "max_cluster_size": max_cluster_size,
        "failing_connected_sets_by_size": {str(s): failing[s] for s in range(1, max_cluster_size + 1)},
        "corrected_through_size": t,
        "activity_radius": f"{x_radius.numerator}/{x_radius.denominator}",
        "theta": fraction_to_outward_strs(theta, digits=digits)[1],
        "noise": "independent faults, the same probability on every mechanism",
    }
    if theta > 1:
        out.update(subthreshold_condition_met=False, ler_lo="0", ler_hi="1")
        return out
    leading = sum((failing[s] * x ** s for s in range(1, max_cluster_size + 1)), Fraction(0))
    remainder = 2 * M * theta ** min(max_cluster_size + 1, t + 2)
    y_hi = 2 * leading + remainder
    y_lo = 2 * leading - remainder
    hi = min(Fraction(1, 2), y_hi / 2)
    lo = max(Fraction(0), (y_lo - y_lo * y_lo / 2) / 2) if y_lo > 0 else Fraction(0)
    if not 0 <= lo <= hi <= Fraction(1, 2):
        raise RefusalError("ENCLOSURE-ORDER", "computed bounds are not ordered")
    lo_str, _ = fraction_to_outward_strs(lo, digits=digits)
    _, hi_str = fraction_to_outward_strs(hi, digits=digits)
    out.update(
        subthreshold_condition_met=True,
        leading_sum=fraction_to_outward_strs(leading, digits=digits)[1],
        remainder=fraction_to_outward_strs(remainder, digits=digits)[1],
        ler_lo=lo_str,
        ler_hi=hi_str,
        lower_bound_nontrivial=lo > 0,
    )
    return out
