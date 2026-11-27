"""Replay native inputs before contracting their joint law through RM15."""

from ftqc_extractors.__main__ import verify
from ftqc_magic.insitu import evaluate
from .supports import fields


def factory(s):
    fields(s, "kind", "inputs", "independent_attempts")
    batches = []
    for row in s["inputs"]:
        fields(row, "spec", "receipt", "ports")
        if row["spec"]["kind"] != "injection_export":
            raise ValueError("expected injection export")
        verify(row["spec"], row["receipt"])
        batches.append((row["receipt"]["result"], row["ports"]))
    return dict(status="proved", factory=evaluate(batches, s["independent_attempts"]))


def retries(s):
    from fractions import Fraction

    fields(
        s,
        "kind",
        "factory",
        "lane_ticks",
        "factory_ticks",
        "factory_width",
        "transcript",
    )
    source = s["factory"]
    f = factory(source)["factory"]
    lanes = s["lane_ticks"]
    if len(lanes) != len(source["inputs"]) or any(
        type(t) is not int or t < 1 for t in lanes
    ):
        raise ValueError("one positive duration per simultaneous preparation")
    specs = [row["spec"] for row in source["inputs"]]
    if len({a["block"] for a in specs}) != len(specs):
        raise ValueError("simultaneous preparations require distinct blocks")
    width, ticks = s["factory_width"], s["factory_ticks"]
    if any(type(t) is not int or t <= 0 for t in (width, ticks)):
        raise ValueError("invalid factory timing")
    end = max(lanes)
    # Lanes retain the full encoded block until the synchronized accept decision.
    prep = sum(row["model"]["injection"]["memory"]["n"] * end for row in specs)
    a, final = Fraction(f["source_accept"]), Fraction(f["final_accept"])
    chances = {"source-reject": 1 - a, "factory-reject": a - final, "delivered": final}
    elapsed = volume = 0
    records = []
    for epoch, outcome in enumerate(s["transcript"]):
        if (
            outcome not in chances
            or not chances[outcome]
            or (records and records[-1]["outcome"] == "delivered")
        ):
            raise ValueError("impossible retry transcript")
        advance = outcome != "source-reject"
        start = elapsed
        elapsed += end + ticks * advance
        volume += prep + width * ticks * advance
        records.append(
            dict(
                epoch=epoch,
                start=start,
                release=elapsed,
                outcome=outcome,
                attempts=[[epoch, a["attempt"]] for a in specs],
            )
        )
    return dict(
        status="proved",
        factory=f,
        records=records,
        wire_ticks=volume,
        elapsed_ticks=elapsed,
        buffered_block_ticks=[end - t for t in lanes],
        yield_per_group=str(final),
        expected_wire_ticks_per_delivery=str((prep + a * width * ticks) / final)
        if final
        else None,
        premise="independent repeated groups; discard all on rejection; declared timing; ideal input decoding and transfer",
    )
