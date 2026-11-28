"""Price reserved data/cat capacity, feedback and the ideal services of a staged study."""

from ftqc_compose.staged import program


def price(study, reaction_ticks, terminal_wire_ticks, seed_wire_ticks):
    if any(
        type(v) is not int or v < 0
        for v in (reaction_ticks, terminal_wire_ticks, seed_wire_ticks)
    ):
        raise ValueError("expected nonnegative integer service budgets")
    result = program(study)
    if result["status"] != "proved":
        raise ValueError("cannot price a failed gadget")
    total = 0
    for row in result["stages"]:
        w = row["workload"]
        active = w["max_ticks"] * (w["data"] + w["peak_cat_legs"])
        wait = w["data"] * w["max_reactions"] * reaction_ticks
        total += active + wait + terminal_wire_ticks + seed_wire_ticks
    return dict(
        action=result["action"],
        reserved_wire_ticks_upper=total,
        model="serial stages; upper service budgets per stage; supplied feedback noise ideal",
        terminal_wire_ticks_per_stage=terminal_wire_ticks,
        seed_wire_ticks_per_stage=seed_wire_ticks,
        feedback_ticks_per_round=reaction_ticks,
        conditional_model_failure_upper=result["failure_upper"],
    )
