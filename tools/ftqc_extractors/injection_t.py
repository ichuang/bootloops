"""Transport an explicit joint raw-carrier Z model through a signed CSS isometry."""

from fractions import Fraction
from .injection import analyze
from .algebra import Limit


def accepted_t(m, budget):
    from .__main__ import digest, validate
    from .injection_export import named

    if set(m) != {"type", "injection", "events"}:
        raise ValueError("unknown raw T model fields")
    spec = m["injection"]
    validate(spec)
    if spec["kind"] != "injection":
        raise ValueError("expected an injection isometry")
    proof = analyze(spec)
    if proof["status"] != "proved":
        raise ValueError("open injection isometry")
    carriers = proof["carriers"]
    ids = set()
    mass = {0: Fraction(1)}
    for e in m["events"]:
        if (
            set(e) != {"id", "sites", "p", "reject"}
            or named(e["id"]) in ids
            or type(e["reject"]) is not bool
        ):
            raise ValueError("invalid or duplicate source event")
        ids.add(e["id"])
        sites = e["sites"]
        if (
            not isinstance(sites, list)
            or len(sites) != len(set(sites))
            or any(type(q) is not int or q not in carriers for q in sites)
        ):
            raise ValueError("Z fault outside raw carrier set")
        mask = sum(1 << carriers.index(q) for q in sites)
        if type(e["p"]) not in (str, int):
            raise ValueError("expected rational probability")
        p = Fraction(e["p"])
        if not 0 <= p <= 1:
            raise ValueError("invalid event probability")
        nxt = {}
        for k, v in mass.items():
            if p < 1:
                nxt[k] = nxt.get(k, 0) + v * (1 - p)
            if not e["reject"] and p:
                nxt[k ^ mask] = nxt.get(k ^ mask, 0) + v * p
        if len(nxt) > budget:
            raise Limit("raw-carrier distribution exceeds state budget")
        mass = nxt
    return spec, dict(
        evidence="exact-T-Z",
        target="T",
        phases=[0, 1],
        modulus=8,
        carriers=carriers,
        acceptance=str(sum(mass.values())),
        accepted_mass={str(i): str(v) for i, v in sorted(mass.items())},
        preparation_sha256=digest(spec),
        isometry=proof,
        events=m["events"],
        premises="independent Bernoulli source events; joint Z per event before ideal signed CSS encoding; flagged events reject the attempt",
    )
