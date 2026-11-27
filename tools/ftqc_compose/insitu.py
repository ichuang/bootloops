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
