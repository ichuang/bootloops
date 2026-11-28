"""Replayable logical operations on ordered, encoded data blocks."""

import hashlib
import json
from ftqc_codes import codes


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def endpoint(code_name, names):
    code = codes.load(code_name)
    full = (1 << code.n) - 1
    return [
        dict(
            block=name,
            code=code_name,
            checks_sha256=code.sha256(),
            width=code.n,
            logical_x=full,
            logical_z=full,
            signs="++",
            frame="identity",
        )
        for name in names
    ]


class Counterexample(Exception):
    def __init__(self, obligation, **witness):
        self.evidence = dict(obligation=obligation, **witness)


def analyze(spec):
    if not isinstance(spec, dict) or spec.get("kind") != "staged-gadget":
        raise ValueError("expected a staged-gadget specification")
    names = spec["blocks"]
    if (
        not isinstance(names, list)
        or not names
        or any(type(n) is not str or not n for n in names)
        or len(set(names)) != len(names)
    ):
        raise ValueError("expected distinct ordered block names")
    rule = spec.get("terminal_recovery")
    if "terminal_recovery" in spec and not isinstance(rule, list):
        raise ValueError("terminal_recovery must be a list of syndrome table rows")
    protocol = {key: value for key, value in spec.items() if key != "terminal_recovery"}
    from . import staged_switch

    try:
        if spec["operation"] == "pieceable-ccz":
            from .staged_pieceable import check

            result = check(protocol, rule)
        elif spec["operation"] == "target-h":
            from .staged_h import check

            result = check(protocol, rule)
        else:
            result = staged_switch.check(protocol, rule)
        from .staged_work import workload

        result["workload"] = workload(spec, result)
    except Counterexample as failure:
        result = dict(status="counterexample", witness=failure.evidence)
    return dict(
        schema="ftqc_gadgets-staged-1",
        spec=spec,
        sha256=hashlib.sha256(canonical(spec).encode()).hexdigest(),
        result=result,
    )


def verify(record):
    if canonical(analyze(record["spec"])) != canonical(record):
        raise ValueError("staged gadget replay mismatch")
