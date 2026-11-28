"""Prepare an encoded T state by executing a replayed gate on a fresh logical plus."""

from fractions import Fraction
from ftqc_compose.staged import program


def prepare(study):
    result = program(study)
    if (
        result["status"] != "proved"
        or result["action"] != "T"
        or len(result["inputs"]) != 1
        or result["inputs"] != result["outputs"]
        or result["inputs"][0]["code"] != "steane"
    ):
        raise ValueError("expected a closed Steane logical T")
    # The fixed-length switch has no reject rule; define a recovery (identity)
    # on unseen terminal syndromes to keep it trace preserving at higher order.
    if any(
        r["workload"]["clean_ticks"] != r["workload"]["max_ticks"]
        for r in result["stages"]
    ):
        raise ValueError("state source requires total deterministic stages")
    return dict(
        state="T",
        port=result["outputs"][0],
        input="fresh ideal logical +",
        acceptance="1",
        conditional_infidelity_upper=str(Fraction(result["failure_upper"])),
        noise=result["noise"],
        premises=[
            dict(terminal=r["terminal"], resources=r["ideal_resources"])
            for r in result["stages"]
        ],
        basis="encoded resource; arbitrary logical noise on a bad attempt",
    )
