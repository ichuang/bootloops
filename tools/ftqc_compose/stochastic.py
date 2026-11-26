"""Exact arithmetic for the uniform local-stochastic subset condition."""

from fractions import Fraction

from .contracts import contract_result
from .level import level_result
from .supports import Family, family_result, fields, finite_domain, names


def probability(value):
    if type(value) not in (str, int):
        raise ValueError(
            "probability must be an exact integer, decimal string or fraction string"
        )
    try:
        rate = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError("invalid exact probability") from error
    if not 0 <= rate <= 1:
        raise ValueError("probability outside [0,1]")
    return rate


def check_law(law, family, rate, ideal_sites):
    """Superset zeta transform of an explicit finite law, including normalization."""
    finite_domain(family.sites)
    if not isinstance(law, list):
        raise TypeError("expected an explicit list of fault-support masses")
    sites = sorted(family.sites)
    index = {s: 1 << i for i, s in enumerate(sites)}
    total = [Fraction(0)] * (1 << len(sites))
    event = Fraction(0)
    for atom in law:
        fields(atom, "support", "mass")
        support = names(atom["support"])
        mass = probability(atom["mass"])
        if mass and support.intersection(ideal_sites):
            raise ValueError("support law faults an ideal physical site")
        if family.bad(support):
            event += mass
        total[sum(index[s] for s in support)] += mass
    if sum(total) != 1:
        raise ValueError("support law must have total mass one")
    for bit in index.values():
        for mask in range(len(total)):
            if not mask & bit:
                total[mask] += total[mask | bit]
    powers = [rate**k for k in range(len(sites) + 1)]
    for mask, q in enumerate(total):
        if q > powers[mask.bit_count()]:
            witness = [s for s, bit in index.items() if mask & bit]
            raise ValueError(f"local-stochastic subset condition fails at {witness}")
    return {"subsets_checked": len(total), "bad_event_mass": str(event)}


def source_cover(source):
    kind = source["kind"]
    if kind == "family":
        return family_result(source), [], "bad-support event", None, []
    if kind == "level":
        r = level_result(source)
        consequence = None
        if r["output_correctness"] is not None:
            consequence = "unconditional output distance for every input in the checked logical study"
        return (
            r["bad_faults"],
            r["premises"],
            "lifted bad-support event",
            consequence,
            r["ideal_physical_sites"],
        )
    if kind == "contract":
        r = contract_result(source)
        return (
            r["bad_faults"],
            r["premises"],
            "bad-support event",
            "joint probability of acceptance and incorrect output on a promised input",
            [],
        )
    raise ValueError("bound needs a checked family, level reduction, or contract")


def stochastic_result(spec):
    keys = ["kind", "p", "source"]
    if "law" in spec:
        keys.append("law")
    fields(spec, *keys)
    rate = probability(spec["p"])
    cover, premises, event, consequence, ideal_sites = source_cover(spec["source"])
    family = Family(cover["sites"], cover["minimal"])
    coeff = cover["enumerator"]
    raw = sum((count * rate ** int(k) for k, count in coeff.items()), Fraction(0))
    bound = min(Fraction(1), raw)
    check = check_law(spec["law"], family, rate, ideal_sites) if "law" in spec else None
    if check and Fraction(check["bad_event_mass"]) > bound:
        raise ValueError("support event exceeds its computed union bound")
    return {
        "status": "conditional" if premises else "proved",
        "rule": "local-stochastic union bound",
        "noise_condition": "Pr[S contained in fault support] <= p^|S| for every physical subset S",
        "noise_evidence": check
        if check is not None
        else "assumed local-stochastic support law",
        "p": str(rate),
        "physical_sites": sorted(family.sites),
        "ideal_physical_sites": ideal_sites,
        "enumerator": coeff,
        "raw_bound": str(raw),
        "upper": str(bound),
        "event": event,
        "also_bounds": consequence,
        "premises": premises,
    }
