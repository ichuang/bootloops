"""Typed, ordered-bundle composition for accepted diagonal-Pauli branches.

A leaf is either a checked construction or an explicit mathematical premise.
All verification is relative to the recorded ideal boundaries and fault model.
"""

from itertools import combinations, product

from .supports import Family, fields, names, disjoint, leaf


SCOPE = "accepted-diagonal-pauli"
BIT = dict(
    kind="C", wires=["0"], code="bit", basis="Z", decoder="identity", bad=[["0"]]
)
REP = dict(
    kind="C",
    wires=["0", "1", "2"],
    code="repetition-3",
    basis="Z",
    decoder="majority",
    bad=[list(t) for t in combinations(["0", "1", "2"], 2)],
)


def validate_port(port):
    fields(port, "kind", "wires", "code", "basis", "decoder", "bad")
    if port["kind"] not in ("Q", "C") or not port["wires"]:
        raise ValueError("expected a nonempty classical or quantum bundle")
    for key in ("code", "basis", "decoder"):
        if not isinstance(port[key], str) or not port[key]:
            raise ValueError(f"invalid port {key}")
    return leaf(dict(op="family", sites=port["wires"], minimal=port["bad"]))


def join(outputs, inputs):
    if len(outputs) != len(inputs):
        raise ValueError("join: ordered port count differs")
    for i, (a, b) in enumerate(zip(outputs, inputs)):
        for key in ("kind", "wires", "code", "basis", "decoder"):
            if a[key] != b[key]:
                raise ValueError(f"join port {i}: {key} differs")
        if not validate_port(a).covers(validate_port(b)):
            raise ValueError(f"join port {i}: output residual guarantee too weak")


def theorem(inputs, outputs, ideal, family, filters, evidence, premises=()):
    if filters is not None and len(filters) != len(inputs):
        raise ValueError("friendly proof requires one filter per input bundle")
    return dict(
        inputs=inputs,
        outputs=outputs,
        ideal=ideal,
        family=family,
        filters=filters,
        evidence=evidence,
        premises=list(premises),
    )


def classical(node):
    op = node["op"]
    if op == "prepare0":
        fields(node, "op", "sites")
        sites = sorted(names(node["sites"]))
        if len(sites) != 3:
            raise ValueError("prepare0 requires three distinct bit locations")
        pairs = Family(sites, combinations(sites, 2))
        for word in product((0, 1), repeat=3):
            support = {site for site, value in zip(sites, word) if value}
            if pairs.bad(support) != (sum(word) >= 2):
                raise ValueError("repetition preparation check failed")
        return theorem(
            [],
            [REP],
            "prepare-encoded-zero",
            pairs,
            [],
            "checked all eight bit-flip supports",
        )
    fields(node, "op")
    if op == "majority":
        for logical, flips in product((0, 1), range(8)):
            word = (7 if logical else 0) ^ flips
            if flips.bit_count() <= 1 and (word.bit_count() >= 2) != logical:
                raise ValueError("classical decoder failed")
        return theorem(
            [REP],
            [BIT],
            "decode-bit",
            Family([], []),
            ["identity"],
            "ideal majority; every input word decodes into one good ball",
        )
    return theorem(
        [BIT],
        [BIT],
        "identity-bit",
        Family([], []),
        ["identity"],
        "ideal classical wire",
    )


def hypothesis(node):
    fields(node, "op", "id", "inputs", "outputs", "ideal", "faults", "filters")
    for key in ("id", "ideal"):
        if not isinstance(node[key], str) or not node[key]:
            raise ValueError(f"missing hypothesis {key}")
    for ports in (node["inputs"], node["outputs"]):
        if not isinstance(ports, list):
            raise ValueError("expected a list of ordered bundles")
        for p in ports:
            validate_port(p)
    filters = node["filters"]
    if filters is not None and (
        not isinstance(filters, list)
        or not all(isinstance(v, str) and v for v in filters)
    ):
        raise ValueError("expected named bundle filters or null")
    f = leaf(node["faults"])
    return theorem(
        node["inputs"],
        node["outputs"],
        node["ideal"],
        f,
        filters,
        "assumed simulation and listed filter obligations",
        [node["id"]],
    )


def derive(node):
    op = node["op"]
    if op in ("prepare0", "majority", "bit-wire"):
        return classical(node)
    if op == "hypothesis":
        return hypothesis(node)
    if op not in ("sequence", "parallel"):
        raise ValueError(f"unknown contract rule {op!r}")
    fields(node, "op", "left", "right")
    a, b = derive(node["left"]), derive(node["right"])
    disjoint([a["family"], b["family"]])
    if op == "sequence":
        join(a["outputs"], b["inputs"])
        inputs, outputs = a["inputs"], b["outputs"]
        filters = a["filters"]
    else:
        inputs, outputs = a["inputs"] + b["inputs"], a["outputs"] + b["outputs"]
        filters = (
            None
            if a["filters"] is None or b["filters"] is None
            else a["filters"] + b["filters"]
        )
    evidence = dict(rule=op, left=a["evidence"], right=b["evidence"])
    return theorem(
        inputs,
        outputs,
        [op, a["ideal"], b["ideal"]],
        a["family"].union(b["family"]),
        filters,
        evidence,
        a["premises"] + b["premises"],
    )


def contract_result(spec):
    fields(spec, "kind", "scope", "program", "require_friendly", "require_closed")
    if spec["scope"] != SCOPE:
        raise ValueError("unsupported fault scope")
    if any(type(spec[k]) is not bool for k in ("require_friendly", "require_closed")):
        raise ValueError("proof requirements must be booleans")
    result = derive(spec["program"])
    if spec["require_friendly"] and result["filters"] is None:
        raise ValueError("missing friendly input-filter proof")
    if spec["require_closed"] and result["premises"]:
        raise ValueError(f"unproved premises: {result['premises']}")
    family = result.pop("family")
    return dict(
        scope=SCOPE,
        status="conditional" if result["premises"] else "proved",
        bad_faults=family.describe(),
        enumerator_scope="minimal generators of sufficient cover",
        **result,
    )
