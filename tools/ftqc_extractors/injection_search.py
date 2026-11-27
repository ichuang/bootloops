"""Finite search over raw sites, oriented Bell pairs and free preparation bases."""

from itertools import combinations, product
from functools import reduce
from operator import or_
from .algebra import bits, Limit
from .injection import analyze, logicals, mask


def pairings(group):
    if not group:
        yield []
        return
    a, *rest = group
    for j, b in enumerate(rest):
        for tail in pairings(rest[:j] + rest[j + 1 :]):
            yield [[a, b], *tail]
            yield [[b, a], *tail]


def search(s, budget):
    n, hx, hz, all_pairs = logicals(s)
    objective = s["objective"]
    count = s["count"]
    if (
        objective not in ("protected", "fixed")
        or type(count) is not int
        or not 1 <= count <= len(all_pairs)
    ):
        raise ValueError("invalid objective or subset size")
    best = None
    evaluated = 0

    def spend():
        nonlocal budget
        budget -= 1
        if budget < 0:
            raise Limit("injection search budget")

    for subset in combinations(range(len(all_pairs)), count):
        spend()
        pairs = [all_pairs[i] for i in subset]
        support = reduce(or_, (x | z for x, z, _, _ in pairs), 0)
        choices = [
            list(
                bits(
                    x
                    & z
                    & ~reduce(
                        or_,
                        (a | b for j, (a, b, _, _) in enumerate(pairs) if j != i),
                        0,
                    )
                )
            )
            for i, (x, z, _, _) in enumerate(pairs)
        ]
        for sites in product(*choices):
            spend()
            plus, zero, free, groups = [], [], [], {}
            for q in range(n):
                if q in sites:
                    continue
                sig = tuple((x >> q & 1, z >> q & 1) for x, z, _, _ in pairs)
                has_x = any(a for a, b in sig)
                has_z = any(b for a, b in sig)
                if has_x and has_z:
                    groups.setdefault(sig, []).append(q)
                else:
                    (plus if has_x else zero if has_z else free).append(q)
            if any(len(g) % 2 for g in groups.values()):
                continue

            # Traverse lazily: a large class must respect the finite search budget.
            def arrangements(todo, bells):
                spend()
                if not todo:
                    yield bells
                else:
                    for new in pairings(todo[0]):
                        yield from arrangements(todo[1:], bells + new)

            for bells in arrangements(list(groups.values()), []):
                for bases in product((0, 1), repeat=len(free)):
                    spend()
                    spec = dict(
                        schema=1,
                        kind="injection",
                        method="selected",
                        memory=s["memory"],
                        logicals=pairs,
                        preparation=dict(
                            carriers=list(sites),
                            bells=bells,
                            plus=plus + [q for q, b in zip(free, bases) if b],
                            zero=zero + [q for q, b in zip(free, bases) if not b],
                        ),
                    )
                    r = analyze(spec)
                    if r["status"] != "proved":
                        raise ValueError("candidate failed signed source test")
                    evaluated += 1
                    fixed = r["fixed_checks"]
                    xcover = reduce(or_, (hx[i] for i in fixed["X"]), 0)
                    zcover = reduce(or_, (hz[i] for i in fixed["Z"]), 0)
                    p = spec["preparation"]
                    protected = support & (
                        (xcover & mask(p["plus"] + [a for a, b in bells]))
                        | (zcover & mask(p["zero"] + [b for a, b in bells]))
                    )
                    score = (
                        protected.bit_count()
                        if objective == "protected"
                        else sum(map(len, fixed.values()))
                    )
                    if best is None or score > best["score"]:
                        best = dict(
                            score=score,
                            selected=list(subset),
                            source=spec,
                            protected=list(bits(protected)),
                        )
    if best is None:
        return dict(
            status="counterexample",
            obligation="no feasible selected subset at requested size",
        )
    return dict(status="proved", objective=objective, evaluated=evaluated, **best)
