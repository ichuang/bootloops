"""Replay the released, closed logical-Y experiments on verified CSS sources."""

import hashlib
import json
from pathlib import Path
from fractions import Fraction
from ftqc_audit.injection_records import audit, experiment_cost
from .injection import analyze

EXAMPLES = Path(__file__).parent / "examples"
CASES = {
    "bb-uniform": ("bb144-injection-4", "uniform", 2, "0.001"),
    "bb-uniform-three": ("bb144-injection-4", "uniform-three", 3, "0.001"),
    "bb-asymmetric": ("bb144-injection-4", "asymmetric", 3, "0.001"),
    "hgp-release": ("hgp225-injection-6", "uniform", 3, "0.0001"),
}


def paths(case):
    if case not in CASES:
        raise ValueError("unknown injection circuit fixture")
    name, suffix, rounds, p = CASES[case]
    return (
        EXAMPLES / (name + ".json"),
        EXAMPLES / "injection-circuits" / (name + "-" + suffix + ".stim"),
        rounds,
        Fraction(p),
    )


def run(s, ledger=None):
    try:
        import stim
    except ImportError as exc:
        from .algebra import Limit

        raise Limit("Stim is required for the closed logical-Y circuit") from exc
    source, path, rounds, p = paths(s["case"])
    spec = json.loads(source.read_text())
    proof = analyze(spec)
    if proof["status"] != "proved":
        raise ValueError("invalid source")
    data = path.read_bytes()
    report = audit(
        stim.Circuit(data.decode()), spec, proof["postselect"], rounds, p, ledger
    )
    cost = experiment_cost(stim.Circuit(data.decode()), report["postselect_width"])
    return dict(
        cost=cost,
        status="proved",
        scope="finite single-location responses of closed logical Y probes",
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        circuit_sha256=hashlib.sha256(data).hexdigest(),
        rounds=rounds,
        p_reference=str(p),
        **report,
    )
