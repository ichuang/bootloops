"""Probabilistic closure of a two-level program: one bound from the noise condition alone.

A study of kind `closure` takes a two-level transversal program (format of
`ftqc_gadgets cones`), the inner and outer codes' n and t, a purification tree for every
inner ancilla kind (as a physical `filter-contract` study) and for every outer ancilla kind
(as an encoded one), and the weight p of a physical location.  It chains three studies
under local-stochastic noise:

    1. the contract of every inner ancilla, from its physical tree;
    2. the contract of every outer ancilla, from its tree of encoded blocks, with the
       inner contracts of step 1;
    3. the tail bound over the cones of the program, with the contracts of steps 1 and 2.

Every number in the chain is the total weight of a family of sets of physical locations
that covers an event: weights of covers on disjoint sources multiply, weights of
alternatives add, and a first accepted attempt costs 1 / (1 - weight of rejection).  No
contract is declared and nothing is assumed independent.  The result bounds the
probability that some outer decoding event of the program sees more than its radius of dirty
columns or that an outer ancilla delivered to the program is junk.  A column is dirty for an
outer event when an inner decoding event of its cone fails there, an inner ancilla that event
consumed is junk, or an outer ancilla in the cone is dirty there; failures of single inner
events are not excluded, they are what the outer code corrects.  The chain exists only when
the weight of rejection is below one at every check, which the study verifies and otherwise
refuses.  What remains is listed in two groups:
the noise condition, and the deterministic statements quoted from the other studies, whose
hypotheses those studies check on the instance.
"""
import hashlib
import json

from .filters import filter_contract_result
from .stochastic import probability
from .supports import fields
from .tails import cone_bound_result

NOISE = "Pr[all locations of S faulty] <= p^|S| for every finite set S of physical locations"
DETERMINISTIC = [
    "the cone theorem for the program, for every raw preparation and for every check (ftqc_gadgets.cones: "
    "hypotheses checked, induction quoted)",
    "the filter rule for every purification tree, with the classification of its locations (junk-capable, local "
    "to one coordinate, able to cause rejection) and of the steps of its checks (ftqc_gadgets.factories: quoted "
    "and tested against propagation)",
    "a check accepts when both incoming patterns of its sector are trivial and no location or decoding event of "
    "the step that can cause rejection is faulty, and a rejected check is repeated with fresh inputs",
    "preparations, steps of checks, attempts, ancillas, program gates and columns occupy pairwise distinct "
    "physical locations, each location belongs to one coordinate of one column, and the noise condition holds "
    "jointly on the locations of all attempts, executed or not",
    "in a tree of encoded blocks, a survivor shares sources with the partner's side of a check only through its "
    "latest correction and only at the same column (checked on one column of the tree)",
]
DISCHARGED = [
    "the ancilla contracts (derived in steps 1 and 2 instead of declared)",
    "independence of faults, of ancillas and of columns (replaced by disjointness of sources)",
    "conditioning on acceptance (replaced by the first accepted attempt)",
]
CONCLUSION = (
    "With probability at least 1 - upper, no outer ancilla delivered to the program is junk and every outer "
    "decoding event of the program has at most its radius of dirty columns; by the cone theorem every outer "
    "decoder then returns the pattern it faces.  Single inner decoding events may fail: their failures are the "
    "dirty columns that the outer code corrects.")


def _sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _tree(kind, frag, common, extra):
    if not isinstance(frag, dict):
        raise ValueError(f"ancilla {kind!r}: a tree description is an object")
    allowed = {"tree", "leaf", "factory", "raw", "modes", "order", "max_modes"}
    if set(frag) - allowed:
        raise ValueError(f"ancilla {kind!r}: unknown fields {sorted(set(frag) - allowed)}")
    r = filter_contract_result(dict(common, **extra, **frag))
    if r["status"] != "conditional":
        raise ValueError(f"ancilla {kind!r}: {r.get('reason', 'no bound')}")
    return r


def closure_result(spec):
    fields(spec, "kind", "p", "inner", "outer", "program", "inner_ancillas", "outer_ancillas")
    probability(spec["p"])
    for group in ("inner_ancillas", "outer_ancillas"):
        if not isinstance(spec[group], dict) or not spec[group]:
            raise ValueError(f"{group}: one tree for every ancilla kind")
    if set(spec["inner_ancillas"]) & set(spec["outer_ancillas"]):
        raise ValueError("an ancilla kind is inner or outer, not both")
    common = dict(kind="filter-contract", noise="local-stochastic", p=spec["p"])
    stages, contracts = [], {}
    for kind in sorted(spec["inner_ancillas"]):
        r = _tree(kind, spec["inner_ancillas"][kind], common, dict(level="physical"))
        contracts[kind] = r["contract"]
        stages.append(dict(step=1, ancilla=kind, level="physical", junk=r["junk"], main_mode=r["main_mode"],
                           accept_lower=r["accept_lower"], expected_cost=r["expected_cost"], sha256=_sha(r)))
    inner = dict(contracts)
    for kind in sorted(spec["outer_ancillas"]):
        r = _tree(kind, spec["outer_ancillas"][kind], common,
                  dict(level="encoded", inner=spec["inner"], contracts=inner))
        contracts[kind] = r["contract"]
        stages.append(dict(step=2, ancilla=kind, level="encoded", junk=r["junk"], main_mode=r["main_mode"],
                           accept_lower=r["accept_lower"], carry_last=r["carry_last"], sha256=_sha(r)))
    bound = cone_bound_result(dict(kind="cone-bound", noise="local-stochastic", p=spec["p"], inner=spec["inner"],
                                   outer=spec["outer"], contracts=contracts, program=spec["program"]))
    if bound["outer"] is None:
        raise ValueError("the program is not a two-level program")
    return dict(
        status="conditional",
        rule="covers by sets of physical locations, chained through both levels",
        p=spec["p"],
        inner=spec["inner"],
        outer=spec["outer"],
        program_sha256=bound["program_sha256"],
        cones=bound["cones"],
        stages=stages,
        by_tag=bound["by_tag"],
        upper=bound["upper"],
        upper_sci=bound["upper_sci"],
        conclusion=CONCLUSION,
        noise_premise=NOISE,
        deterministic_premises=DETERMINISTIC,
        no_longer_premises=DISCHARGED,
    )
