"""Finite level reduction, He–Nguyen–Pattison Theorem 4.11.

The circuit is defined by expanding the typed vertices in their supplied
topological order. Non-native strong local simulation lemmas are premises.
"""

from itertools import product
from math import prod

from .contracts import BIT, derive, join, theorem, validate_port
from .supports import Family, disjoint, fields, leaf, names

QBIT = dict(BIT, kind="Q", code="qubit", basis="computational")
SCOPE = "unconditional-level-reduction"


def local(g):
    op = g["op"]
    if op == "assume-4.8":
        fields(g, "op", "contract")
        h = g["contract"]
        if h.get("op") != "hypothesis":
            raise ValueError("expected an explicit Definition 4.8 hypothesis")
        result = derive(h)
        if result["filters"] is None:
            raise ValueError("level reduction requires bundle-wise friendly filters")
        result["evidence"] = (
            "premise: unconditional Def. 4.8 for all avoiding Pauli pairs with environment"
        )
        return result
    if op in ("magic-source", "rm15-ideal"):
        fields(g, "op", "site")
        sites = names([g["site"]])
        if op == "magic-source":
            return theorem(
                [],
                [QBIT],
                "prepare-T",
                Family(sites, [sites]),
                [],
                "unencoded fresh T|+>; any faulty preparation is bad",
            )
        return theorem(
            [QBIT] * 15,
            [QBIT],
            "RM15-inject-correct-decode",
            Family(sites, []),
            ["identity"] * 15,
            "fault-free trace-preserving distillation channel on trivial bundles",
        )
    if op == "prepare0":
        # Fresh preparations, then a classical readout at each of the three sites.
        # Good sites are exactly zero; the at-most-one touched bit is an allowed
        # residual on an encoded zero, also for off-diagonal Pauli terms.
        r = derive(g)
        r["ideal"] = "zero"
        r["evidence"] = (
            "two untouched zero bits suffice for every avoiding support; no input filter"
        )
        return r
    if op in ("majority", "bit-wire", "vote3"):
        fields(g, "op", "site")
        site = names([g["site"]])
        if op == "vote3":
            return theorem(
                [BIT] * 3,
                [BIT],
                "majority-of-three",
                Family(site, []),
                ["identity"] * 3,
                "ideal classical truth table on trivial bundles",
            )
        r = derive({"op": op})
        r["family"] = Family(site, [])
        # The physical decoder depends only on the decoded logical bit, even
        # for arbitrary input words. Its filter decodes and re-encodes.
        if op == "majority":
            r["filters"] = ["re-encode majority"]
            r["evidence"] = (
                "majority(repeat(majority(w))) == majority(w) on all 8 classical words"
            )
            r["ideal"] = "identity-bit"
        return r
    raise ValueError(
        "unsupported unconditional gadget; accepted-run receipts cannot discharge Theorem 4.11"
    )


def lift(outer, blocks):
    sites = disjoint(blocks.values())
    if set(blocks) != outer.sites:
        raise ValueError("outer fault domain must match all contracted vertices")
    # Each outer generator independently chooses one local bad generator at
    # each participating vertex (Def. 2.7). The empty Cartesian product is one.
    sizes = [prod(len(blocks[k].terms) for k in t) for t in outer.terms]
    if sum(sizes) > 200000:
        raise ValueError("explicit level-reduction expansion exceeds 200000 terms")
    terms = (
        frozenset().union(*choices)
        for t in outer.terms
        for choices in product(*(blocks[k].terms for k in sorted(t)))
    )
    return Family(sites, terms)


def bit_gate(op, inputs):
    if op == "zero":
        return 0
    if op == "identity-bit":
        return inputs[0]
    if op == "majority-of-three":
        return int(sum(inputs) >= 2)
    raise ValueError("classical output checker needs supported logical truth tables")


def check_classical(vertices, external, outputs, outer):
    if any(
        len(v["outputs"]) != 1
        or v["ideal"] not in ("zero", "identity-bit", "majority-of-three")
        for v in vertices
    ):
        raise ValueError("unsupported logical truth table for output check")
    if len(vertices) + len(external) > 16:
        raise ValueError("classical output check exceeds 16 assigned bits")
    tested = 0
    for initial in product((0, 1), repeat=len(external)):
        ideal = dict(zip(external, initial, strict=True))
        for v in vertices:
            ideal[v["outputs"][0]] = bit_gate(
                v["ideal"], [ideal[i] for i in v["inputs"]]
            )
        for assignment in product((0, 1), repeat=len(vertices)):
            actual = dict(zip(external, initial, strict=True))
            faulty = set()
            for v, bit in zip(vertices, assignment, strict=True):
                if bit != bit_gate(v["ideal"], [actual[i] for i in v["inputs"]]):
                    faulty.add(v["id"])
                actual[v["outputs"][0]] = bit
            tested += 1
            if any(actual[i] != ideal[i] for i in outputs) and not outer.bad(faulty):
                raise ValueError(
                    f"outer family misses a harmful logical execution: {sorted(faulty)}"
                )
    return {"method": "exhaustive logical wire assignments", "assignments": tested}


def trivial(port):
    # One unencoded bit or qubit per bundle. Multiple such bundles are allowed.
    return port in (BIT, QBIT)


def level_result(spec):
    fields(
        spec,
        "kind",
        "inputs",
        "outputs",
        "physical_sites",
        "vertices",
        "outer",
        "output_check",
        "require_closed",
    )
    if type(spec["require_closed"]) is not bool:
        raise ValueError("require_closed must be Boolean")
    if not isinstance(spec["inputs"], dict) or not isinstance(spec["vertices"], list):
        raise TypeError("expected input bundle map and ordered vertices")
    names(list(spec["inputs"]))
    live = dict(spec["inputs"])
    for port in live.values():
        validate_port(port)
        if not trivial(port):
            raise ValueError("level-reduction external input must be trivial")
    wires = set(live)
    blocks = {}
    vertices = []
    owners = {}
    ideal_sites = set()
    premises = []
    for v in spec["vertices"]:
        fields(v, "id", "inputs", "outputs", "gadget")
        key = v["id"]
        names([key])
        if key in blocks:
            raise ValueError("duplicate contracted vertex")
        names(v["inputs"])
        outs = names(v["outputs"])
        if not set(v["inputs"]) <= live.keys() or outs & wires:
            raise ValueError("missing, reused, or out-of-order bundle")
        r = local(v["gadget"])
        if v["gadget"]["op"] in ("bit-wire", "majority", "vote3", "rm15-ideal"):
            ideal_sites.update(r["family"].sites)
        join([live.pop(i) for i in v["inputs"]], r["inputs"])
        if len(outs) != len(r["outputs"]):
            raise ValueError("output arity mismatch")
        for wire, port in zip(v["outputs"], r["outputs"], strict=True):
            live[wire] = port
            owners[wire] = key
        wires.update(outs)
        blocks[key] = r["family"]
        premises.extend(r["premises"])
        vertices.append(
            {
                "id": key,
                "inputs": v["inputs"],
                "outputs": v["outputs"],
                "ideal": r["ideal"],
                "physical_sites": sorted(r["family"].sites),
                "filters": r["filters"],
                "evidence": r["evidence"],
            }
        )
    outputs = names(spec["outputs"])
    if outputs != live.keys() or any(not trivial(p) for p in live.values()):
        raise ValueError(
            "external output bundles must be exactly the unused, trivial bundles"
        )
    if any(blocks[owners[w]].terms for w in outputs if w in owners):
        raise ValueError("an output-adjacent gadget has bad fault paths")
    if names(spec["physical_sites"]) != disjoint(blocks.values()):
        raise ValueError("gadget domains must partition the physical sites")
    if premises and spec["require_closed"]:
        raise ValueError("undischarged Definition 4.8 premise")
    outer = leaf(spec["outer"])
    cover = lift(outer, blocks)
    method = spec["output_check"]
    if method == "classical":
        if premises or any(p != BIT for p in spec["inputs"].values()):
            raise ValueError("logical output check requires native classical gates")
        check = check_classical(vertices, spec["inputs"], spec["outputs"], outer)
    elif method == "distill15":
        from .distill import output_check

        if premises:
            raise ValueError("distillation requires native unconditional gadgets")
        check = output_check(vertices, spec["inputs"], spec["outputs"], outer)
    elif method == "none":
        check = None
    else:
        raise ValueError("unknown logical output check")
    return {
        "status": "conditional" if premises else "proved",
        "scope": SCOPE,
        "rule": "Theorem 4.11",
        "logical_vertices": vertices,
        "ideal_physical_sites": sorted(ideal_sites),
        "logical_faults": "general, with environment",
        "conclusion": "avoiding physical faults simulate some avoiding logical fault",
        "premises": premises,
        "output_correctness": check,
        "bad_faults": cover.describe(),
    }
