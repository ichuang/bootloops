"""A checked accepted-Pauli leaf for a standalone or placed projective measurement."""

from itertools import combinations
from math import comb

from ftqc_gadgets.surgery import verify
from .contracts import theorem, Violation
from .supports import Family, fields


def leaf(node):
    fields(
        node,
        "op",
        "instance",
        "receipt",
        "wires",
        *(["incoming"] if "incoming" in node else []),
    )
    r = verify(node["receipt"])
    s = r["spec"]
    wires = node["wires"]
    if (
        not isinstance(wires, list)
        or len(wires) != s["n"]
        or len(set(wires)) != len(wires)
        or not all(isinstance(w, str) and w for w in wires)
    ):
        raise ValueError("expected named physical data wires in order")
    instance = node["instance"]
    if not isinstance(instance, str) or not instance:
        raise ValueError("expected instance name")
    if r["faults"]["status"] == "counterexample":
        raise Violation(r["faults"])
    if r["faults"]["status"] != "proved":
        raise ValueError("surgery leaf requires a closed fault exclusion")
    import json

    def port(checks):
        return dict(
            kind="Q",
            wires=wires,
            code=json.dumps(checks),
            basis="physical signed Pauli",
            decoder="clean codespace, recorded projective branch",
            bad=[[q] for q in wires],
        )

    outputs = [port(s["checks"])]
    if "basis" in s:
        named = {p["name"]: p for p in s["basis"]["input"]}
        outputs[0]["basis"] = json.dumps(
            [named[name] for name in s["basis"]["surviving"]], sort_keys=True
        )
    for i in range(len(s["targets"])):
        w = f"{instance}/bit/{i}"
        outputs.append(
            dict(
                kind="C",
                wires=[w],
                code="bit",
                basis="Z",
                decoder="identity",
                bad=[[w]],
            )
        )
    sites = [
        f"{instance}/loc/{i}"
        for i in range(
            sum(not op[0].startswith("!") for row in s["ticks"] for op in row)
        )
    ]
    d = r["faults"]["lower"]
    if comb(len(sites), d) > 1000000:
        raise ValueError("surgery support family exceeds enumeration budget")
    input_port = port(s["checks"])
    if "basis" in s:
        input_port["basis"] = json.dumps(s["basis"]["input"], sort_keys=True)
    evidence = {}
    if "incoming" in node:
        from ftqc_gadgets.surgery import Instrument
        from ftqc_gadgets.surgery_faults import residual_profile

        weight = node["incoming"]
        profile = residual_profile(Instrument(s), weight, d - 1, r["cap"])
        if profile["status"] != "proved":
            raise ValueError("incoming and internal search incomplete")
        if any(
            row["wrong_public"] or row["residual_weight_upper"]
            for row in profile["completed"]
        ):
            raise Violation(profile)
        input_port["bad"] = [list(v) for v in combinations(wires, weight + 1)]
        evidence["joint_input_internal"] = profile
    return theorem(
        [input_port],
        outputs,
        dict(project=s["targets"]),
        Family(sites, combinations(sites, d)),
        None,
        dict(
            **evidence,
            circuit=r["circuit_sha256"],
            readouts=s["readouts"],
            terminal_syndrome=s["terminal_syndrome"],
            branches="residual identity modulo checks and measured targets on accepted branch",
        ),
    )


def channel(spec):
    """Exact accepted law and monotone upper envelope by admitted input weight."""
    from fractions import Fraction
    from ftqc_gadgets.surgery import Instrument
    from ftqc_gadgets.surgery_faults import columns, effects
    from ftqc_gadgets.surgery_operational import conditional_law

    fields(spec, "kind", "receipt", "incoming_weight", "p", "cap")
    r = verify(spec["receipt"])
    inst = Instrument(r["spec"])
    internal, inputs = columns(inst)
    radius, cap = spec["incoming_weight"], spec["cap"]
    if type(radius) != int or not 0 <= radius <= inst.n or type(cap) != int or cap < 1:
        raise ValueError("invalid input or state budget")
    # Budget also bounds distinct incoming signatures; repeated signatures reuse the same law.
    if sum(comb(inst.n, i) * 3**i for i in range(radius + 1)) > cap:
        return dict(status="inconclusive", reason="incoming enumeration")
    cache = {}
    envelope = []
    worst = Fraction(0)
    witness = []
    for weight in range(radius + 1):
        for incoming, sig in effects(inputs, weight):
            if sig not in cache:
                law = conditional_law(internal, spec["p"], cap, initial=sig)
                if law["status"] != "exact":
                    return law
                cache[sig] = law
            law = cache[sig]
            if Fraction(law["acceptance"]):
                bad = sum(
                    (
                        Fraction(v["p"])
                        for v in law["conditional"]
                        if v["readout"] or v["residual"]
                    ),
                    Fraction(0),
                )
                if bad > worst:
                    worst, witness = bad, incoming
        envelope.append(
            dict(
                maximum_input_weight=weight,
                accepted_damage_upper=str(worst),
                witness=witness,
            )
        )
    return dict(
        status="proved",
        circuit=r["circuit_sha256"],
        independent_physical_locations=True,
        categorical_depolarizing_probability=str(Fraction(spec["p"])),
        clean_input=cache[0, 0, 0],
        envelope=envelope,
        scope="fixed-rate joint channel; envelope over nonzero-acceptance inputs",
    )
