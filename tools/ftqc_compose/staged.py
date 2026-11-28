"""Serial composition of closed, order-one staged channels.

Location sets are disjoint by operation instance. Four syndrome rounds per
pieceable cut reserve a fixed domain, including slots reached after faults.
A failure to settle by the cap is a failed attempt.
"""

from fractions import Fraction
from math import comb
from ftqc_gadgets import staged
from .supports import fields, names


def program(spec):
    fields(spec, "kind", "target", "stages", "instances", "p")
    if spec["kind"] != "staged-program":
        raise ValueError("expected staged-program")
    stages = spec["stages"]
    if (
        not isinstance(stages, list)
        or not stages
        or len(names(spec["instances"])) != len(stages)
    ):
        raise ValueError("each stage needs a distinct location namespace")
    if type(spec["p"]) not in (str, int):
        raise ValueError("p must be an exact rational")
    p = Fraction(spec["p"])
    if not 0 <= p <= 1:
        raise ValueError("p outside [0,1]")
    rows = []
    previous = None
    actions = []
    for instance, item in zip(spec["instances"], stages):
        if item.get("schema") == "ftqc_gadgets-staged-1":
            staged.verify(item)
            rec = item
        else:
            rec = staged.analyze(item)
        result = rec["result"]
        if result["status"] != "proved":
            return dict(
                status="counterexample", instance=instance, witness=result["witness"]
            )
        if previous is not None and previous != result["inputs"]:
            raise ValueError(
                "stage interface mismatch: block order, code, basis or frame"
            )
        previous = result["outputs"]
        if result["action"] != "I":
            actions.append(result["action"])
        faults = result["faults"]
        rows.append(
            dict(
                instance=instance,
                stage_sha256=rec["sha256"],
                inputs=result["inputs"],
                outputs=previous,
                terminal=faults["terminal"],
                ideal_resources=faults["ideal_resources"],
                workload=result["workload"],
            )
        )
    if spec["target"] not in ("I", "T", "CCZ") or actions != (
        [] if spec["target"] == "I" else [spec["target"]]
    ):
        raise ValueError("composed logical action differs from requested target")
    count = sum(r["workload"]["max_locations"] for r in rows)
    # Clean corrected output at each join: under the stated model a bad stage
    # needs two local faults. No independence between stages is required.
    bound = min(
        Fraction(1), sum(comb(r["workload"]["max_locations"], 2) * p * p for r in rows)
    )
    return dict(
        status="proved",
        action=spec["target"],
        stages=rows,
        inputs=rows[0]["inputs"],
        outputs=previous,
        max_locations=count,
        failure_upper=str(bound),
        noise="joint pair probabilities <= p^2 on reserved location slots",
        boundary="clean encoded input; ideal declared resources and terminal recoveries",
    )
