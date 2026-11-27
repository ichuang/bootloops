"""Volume of a schedule: expected physical fault locations, from counts and ancilla costs.

`ftqc_gadgets.cones.ledger` counts what a schedule spends.  This module prices the counts:

    online        n1 * n2 * block-ticks: one location per physical qubit of a live block per tick
    corrections   every inner correction consumes one |+> and one |0> ancilla in each of the n2
                  columns
    resources     every delivered resource of a costed kind (outer ancillas, the two of every
                  outer correction included; inner ancillas used as data)

A cost is the expected number of physical locations spent to deliver one accepted ancilla of a
kind, repetitions after rejection included (`filter-contract` reports it as expected_cost).
Costs are inputs of the study; kinds named in `uncosted` are counted and not priced, and the
result says so.  Blocks of kind std are existing blocks and cost nothing.

Not counted, by this convention: the gate and measurement locations of an ancilla while it is
used by a correction (its cost covers its preparation only).
"""
from fractions import Fraction

from .filters import decimal
from .supports import fields
from .tails import sci, up

FREE = ("std",)


def cost_of(value, where):
    if type(value) not in (str, int):
        raise ValueError(f"{where}: a cost is an exact integer or decimal string")
    try:
        c = Fraction(value)
    except (ValueError, ZeroDivisionError):
        raise ValueError(f"{where}: invalid cost") from None
    if c < 0:
        raise ValueError(f"{where}: a cost is not negative")
    return c


def price(counts, n1, n2, costs, uncosted=()):
    """Parts and total of a ledger.  n2 is 1 for a single-level program."""
    inner = ("zero", "plus", "cat")
    parts = dict(online=Fraction(n1 * n2 * counts["block_ticks"]))
    if counts["inner_corrections"]:
        for k in ("plus", "zero"):
            if k not in costs:
                raise ValueError(f"no cost for the {k} ancillas of inner corrections")
        parts["inner_corrections"] = counts["inner_corrections"] * n2 * (costs["plus"] + costs["zero"])
    else:
        parts["inner_corrections"] = Fraction(0)
    res, skipped = {}, {}
    for kind, count in counts["resources"].items():
        if kind in FREE:
            continue
        if kind in uncosted:
            skipped[kind] = count
            continue
        if kind not in costs:
            raise ValueError(f"no cost for resource kind {kind!r}; name it in uncosted to leave it out")
        res[kind] = count * costs[kind] * (n2 if kind in inner else 1)
    parts["resources"] = res
    total = parts["online"] + parts["inner_corrections"] + sum(res.values(), Fraction(0))
    return parts, total, skipped


def program_volume(program, n1, n2, costs, uncosted=()):
    """(ledger, parts, total, uncosted counts) of a program in the format of ftqc_gadgets cones."""
    from ftqc_gadgets import cones

    flow = cones.Flow(cones.Program.from_json(program), two=True if n2 is not None else None)
    if flow.two and n2 is None:
        raise ValueError("a two-level program needs the outer code's n")
    counts = cones.ledger(flow)
    parts, total, skipped = price(counts, n1, n2 or 1, costs, uncosted)
    return counts, parts, total, skipped


def volume_result(spec):
    """Study kind `volume`."""
    keys = ["kind", "inner", "program", "costs"]
    for opt in ("outer", "uncosted"):
        if opt in spec:
            keys.append(opt)
    fields(spec, *keys)
    n = {}
    for name in ("inner", "outer"):
        if name in spec:
            fields(spec[name], "n")
            if type(spec[name]["n"]) is not int or spec[name]["n"] < 1:
                raise ValueError(f"{name}: n is a positive integer")
            n[name] = spec[name]["n"]
    if not isinstance(spec["costs"], dict):
        raise ValueError("costs: {kind: expected locations per delivered ancilla}")
    costs = {k: cost_of(v, f"cost of {k}") for k, v in spec["costs"].items()}
    uncosted = spec.get("uncosted", [])
    if not isinstance(uncosted, list) or set(uncosted) & set(costs):
        raise ValueError("uncosted is a list of kinds without a cost")
    from ftqc_gadgets import cones

    counts, parts, total, skipped = program_volume(spec["program"], n["inner"], n.get("outer"), costs, uncosted)
    return dict(
        status="computed", rule="online block-ticks plus the cost of every ancilla consumed",
        inner=spec["inner"], outer=spec.get("outer"),
        program_sha256=cones.Program.from_json(spec["program"]).sha256(), counts=counts,
        peak_physical_qubits=counts["peak_blocks"] * n["inner"] * n.get("outer", 1),
        parts=dict(online=decimal(parts["online"]), inner_corrections=decimal(up(parts["inner_corrections"], 15)),
                   resources={k: decimal(up(v, 15)) for k, v in sorted(parts["resources"].items())}),
        total=decimal(up(total, 15)), total_sci=sci(total, 4), uncosted=dict(sorted(skipped.items())),
        premises=["every cost is the expected number of physical locations per delivered ancilla of its kind",
                  "resources named in uncosted are not part of the total",
                  "the locations an ancilla occupies while a correction uses it are not counted"])
