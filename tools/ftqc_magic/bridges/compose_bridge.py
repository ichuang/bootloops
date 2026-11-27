"""Read-only bridge exporting `ftqc_compose` (`tails.py` / `cone-bound` / `logical_t`) contracts.

Strictly enforces correlation and catalyst-poisoning discipline at the boundary without
modifying a single file in `tools/ftqc_compose/`.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any

from ..ir import CorrelationMode, MagicRefusalError, MagicSourceContract
from ..rings import DEFAULT_SIG, round_down_str, round_up_str


def to_compose_ancilla_contract(
    contract: MagicSourceContract,
    *,
    atom_tag: str = "and1",
    interleaving: str | None = None,
    terminal_single_shot: bool = False,
    sig: int = DEFAULT_SIG,
) -> dict[str, Any]:
    """Export an `ftqc_compose` `cone-bound` / `filter-contract`-compatible ancilla mixture dict!

    Enforces:
    - Refuses `CATALYST_POISONED` streams unless `terminal_single_shot=True`!
    - Refuses `BATCH_CORRELATED` (`k_outputs > 1`) unless `interleaving == "distinct-instances"`!
    """
    if contract.correlation_mode == CorrelationMode.CATALYST_POISONED and not terminal_single_shot:
        raise MagicRefusalError(
            "CATALYST-POISONING-HAZARD",
            f"Cannot export CATALYST_POISONED source {contract.name!r} as an independent multi-shot "
            "ancilla stream for ftqc_compose cone-bound without terminal_single_shot=True!",
        )
    if contract.correlation_mode == CorrelationMode.BATCH_CORRELATED and interleaving != "distinct-instances":
        raise MagicRefusalError(
            "MULTI-OUTPUT-CORRELATION-HAZARD",
            f"Cannot export BATCH_CORRELATED multi-output source {contract.name!r} (k={contract.k_outputs}) "
            "as independent ancillas without interleaving='distinct-instances'!",
        )

    fail_q: Fraction = contract.conditional_infidelity_upper
    clean_q: Fraction = max(Fraction(0), Fraction(1) - fail_q)
    fail_hi = round_up_str(fail_q, sig)
    pi_lo = round_down_str(clean_q, sig)

    return {
        "atom_tag": atom_tag,
        "source_name": contract.name,
        "target": contract.target.label,
        "correlation_mode": contract.correlation_mode.value,
        "interleaving": interleaving,
        "independent_contract": {
            "fail": fail_hi,
            "pi": pi_lo,
            "modes": [[1, fail_hi]],
        },
        "local_stochastic_contract": {
            "equivalents": {atom_tag: fail_hi},
        },
    }


def to_compose_logical_t_summary(
    contract: MagicSourceContract,
    *,
    sig: int = DEFAULT_SIG,
) -> dict[str, Any]:
    """Export summary bounds compatible with `ftqc_compose/logical_t.py` resource ledgers."""
    return {
        "factory_name": contract.name,
        "target": contract.target.label,
        "distance": contract.distance,
        "acceptance_prob_lo": round_down_str(contract.acceptance_lower, sig),
        "conditional_failure_hi": round_up_str(contract.conditional_infidelity_upper, sig),
        "unnormalized_accepted_failure_hi": round_up_str(contract.unnormalized_failure_upper, sig),
        "expected_wire_ticks_per_output_hi": round_up_str(contract.expected_wire_ticks_per_output_upper(), sig),
    }
