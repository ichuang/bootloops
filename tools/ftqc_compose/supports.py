"""Antichains represent upward-closed bad events on named physical locations."""
from collections import Counter
from itertools import combinations, product



def fields(obj, *keys):
    if not isinstance(obj, dict) or set(obj) != set(keys):
        raise ValueError(f"expected fields {keys}")


def names(items):
    if not isinstance(items, list) or any(not isinstance(i, str) or not i for i in items):
        raise ValueError("location IDs must be nonempty strings in a list")
    if len(items) != len(set(items)):
        raise ValueError("duplicate location ID")
    return frozenset(items)


class Family:
    def __init__(self, sites, terms):
        self.sites = frozenset(sites)
        ordered = sorted(set(map(frozenset, terms)), key=lambda t: (len(t), sorted(t)))
        self.terms = []
        smaller = []
        last_size = -1
        for term in ordered:
            if not term <= self.sites:
                raise ValueError("bad support outside its location domain")
            if len(term) != last_size:
                smaller = list(self.terms)
                last_size = len(term)
            if not any(prior <= term for prior in smaller):
                self.terms.append(term)

    def bad(self, support):
        if not support <= self.sites:
            raise ValueError("support outside domain")
        return any(t <= support for t in self.terms)

    def covers(self, other):
        """Every event bad for other is bad for self (paper Def. 2.9)."""
        if self.sites != other.sites:
            raise ValueError("containment requires the same location domain")
        return all(self.bad(t) for t in other.terms)

    def describe(self):
        counts = Counter(map(len, self.terms))
        return dict(sites=sorted(self.sites), minimal=[sorted(t) for t in self.terms],
                    enumerator={str(k): v for k, v in sorted(counts.items())})

    def union(self, other):
        return Family(self.sites | other.sites, [*self.terms, *other.terms])

    def product(self, other):
        disjoint([self, other])
        return Family(self.sites | other.sites, (x | y for x in self.terms for y in other.terms))


def disjoint(families):
    used = set()
    for f in families:
        overlap = used.intersection(f.sites)
        if overlap:
            raise ValueError(f"shared physical locations: {sorted(overlap)}")
        used.update(f.sites)
    return used


def leaf(node):
    fields(node, "op", "sites", "minimal")
    if not isinstance(node["minimal"], list):
        raise ValueError("minimal must be a list of supports")
    return Family(names(node["sites"]), [names(t) for t in node["minimal"]])


def evaluate(node):
    op = node["op"]
    if op == "family":
        return leaf(node)
    if op in ("union", "product"):
        fields(node, "op", "left", "right")
        a, b = evaluate(node["left"]), evaluate(node["right"])
        return a.union(b) if op == "union" else a.product(b)
    if op == "substitute":
        fields(node, "op", "outer", "at")
        outer = evaluate(node["outer"])
        inner = {k: evaluate(v) for k, v in node["at"].items()}
        if inner.keys() != outer.sites:
            raise ValueError("substitution must cover every outer site exactly once")
        sites = disjoint(inner.values())
        terms = (frozenset().union(*choice) for t in outer.terms
                 for choice in product(*(inner[k].terms for k in sorted(t))))
        return Family(sites, terms)
    raise ValueError(f"unknown family operation {op!r}")


def predicate(node, support):
    """Independent Boolean semantics. Does not use antichain operations."""
    op = node["op"]
    if op == "family":
        return any(all(site in support for site in term) for term in node["minimal"])
    if op == "union":
        return predicate(node["left"], support) or predicate(node["right"], support)
    if op == "product":
        return predicate(node["left"], support) and predicate(node["right"], support)
    if op == "substitute":
        flags = {k for k, child in node["at"].items() if predicate(child, support)}
        return predicate(node["outer"], flags)
    raise ValueError("unknown predicate")


def family_result(spec):
    fields(spec, "kind", "expression")
    node = spec["expression"]
    family = evaluate(node)
    if len(family.sites) > 16:
        raise ValueError("truth-table checker supports at most 16 physical sites")
    for size in range(len(family.sites) + 1):
        for sites in combinations(sorted(family.sites), size):
            if family.bad(frozenset(sites)) != predicate(node, set(sites)):
                raise ValueError("family disagrees with Boolean semantics")
    return dict(status="verified", interpretation="exact minimal supports of expression",
                checked_subsets=2**len(family.sites), **family.describe())
