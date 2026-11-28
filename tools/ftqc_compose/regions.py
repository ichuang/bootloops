"""Dominant failure mechanisms of a bad-support family under several noise rates.

A study of kind `regions` takes a checked source (a family, a level reduction or a contract,
as for `local-stochastic`) and a map from its physical sites to named kinds of location
(gate, idle, measurement, ...).  If a location of kind k is faulty with probability
c_k eps^(w_k), the union bound over the minimal bad supports is a sum of monomials
count * prod_k (c_k eps^(w_k))^(n_k), one for every vector n of fault counts by kind.  As
eps -> 0 the bound is governed by

    tau(w) = min over count vectors n of  w . n,

a tropical polynomial in the exponents w.  It is linear on the cones of w on which one vertex
of the Newton polytope of the count vectors is minimal; those cones are the regions, and the
vertex names the kinds of fault that limit the bound there.  With independent faults the
probability of the bad-support event itself has the same exponent and leading coefficient,
because it is at least the probability that exactly one minimal support is faulty.

The result lists the count vectors with multiplicities, the vectors that are minimal in the
componentwise order (every vertex that is minimal for some positive w is among them), for two
kinds the vertices of the lower hull and the exponent ratios at which the dominant vertex
changes, and for every requested exponent vector the value of tau, the dominant vectors and
their total multiplicity.  Arithmetic is exact.
"""
from collections import Counter
from fractions import Fraction

from .stochastic import source_cover
from .supports import fields


def _exponent(x):
    if type(x) not in (str, int):
        raise ValueError("an exponent is an exact number given as a string or an integer")
    v = Fraction(x)
    if v <= 0:
        raise ValueError("exponents must be positive")
    return v


def count_vectors(minimal, kind_of, kinds):
    out = Counter()
    for term in minimal:
        c = Counter(kind_of[s] for s in term)
        out[tuple(c.get(k, 0) for k in kinds)] += 1
    return out


def pareto_minimal(vectors):
    vs = sorted(vectors)
    return [v for v in vs if not any(u != v and all(a <= b for a, b in zip(u, v)) for u in vs)]


def lower_hull(points):
    """Vertices of a planar point set that minimise w . n for some w > 0, by increasing first
    coordinate, with the ratio w_2 / w_1 at which consecutive vertices tie."""
    pts = pareto_minimal(points)                       # first coordinate up, second down
    hull = []
    for q in pts:
        while len(hull) >= 2:
            (a1, b1), (a2, b2) = hull[-2], hull[-1]
            # drop hull[-1] when it lies on or above the segment from hull[-2] to q
            if (a2 - a1) * (q[1] - b1) - (b2 - b1) * (q[0] - a1) <= 0:
                hull.pop()
            else:
                break
        hull.append(q)
    ties = [Fraction(b[0] - a[0], a[1] - b[1]) for a, b in zip(hull, hull[1:])]
    return hull, ties


def regions_result(spec):
    keys = ["kind", "source", "kinds"]
    if "exponents" in spec:
        keys.append("exponents")
    fields(spec, *keys)
    cover, premises, event, _, ideal = source_cover(spec["source"])
    sites = [s for s in cover["sites"] if s not in set(ideal)]
    kind_of = spec["kinds"]
    if not isinstance(kind_of, dict) or set(kind_of) != set(sites) or any(
            not isinstance(v, str) or not v for v in kind_of.values()):
        raise ValueError("kinds must name a kind for every non-ideal physical site and no other")
    kinds = sorted(set(kind_of.values()))
    minimal = [t for t in cover["minimal"] if not set(t) & set(ideal)]
    vectors = count_vectors(minimal, kind_of, kinds)
    if not vectors:
        raise ValueError("the source has no bad support")
    out = {
        "status": "conditional" if premises else "proved",
        "rule": "tropical exponent of the union bound over minimal bad supports",
        "event": event,
        "kinds": kinds,
        "count_vectors": [{"counts": list(v), "supports": c} for v, c in sorted(vectors.items())],
        "minimal_vectors": [list(v) for v in pareto_minimal(vectors)],
        "premises": premises,
    }
    if len(kinds) == 2:
        hull, ties = lower_hull(list(vectors))
        out["lower_hull"] = [list(v) for v in hull]
        out["ties_at_exponent_ratio"] = {"ratio": f"w[{kinds[1]}] / w[{kinds[0]}]", "values": [str(t) for t in ties]}
    rows = []
    for w in spec.get("exponents", []):
        if not isinstance(w, dict) or set(w) != set(kinds):
            raise ValueError("an exponent vector gives one exponent for every kind")
        ws = [_exponent(w[k]) for k in kinds]
        value = {v: sum(a * b for a, b in zip(ws, v)) for v in vectors}
        tau = min(value.values())
        dom = sorted(v for v in vectors if value[v] == tau)
        rows.append({"exponents": {k: str(x) for k, x in zip(kinds, ws)}, "tau": str(tau),
                     "dominant": [list(v) for v in dom], "supports": sum(vectors[v] for v in dom)})
    if rows:
        out["exponent_vectors"] = rows
    return out
