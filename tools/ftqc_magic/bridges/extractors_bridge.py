"""Read-only bridge exporting `ftqc_extractors` PBC/magic records and `ftqc_uncompute` pricing kwargs.

Leaves `tools/ftqc_uncompute/` and `tools/ftqc_extractors/` 100% untouched.
"""

from __future__ import annotations

from typing import Any

from ..ir import MagicSourceContract
from ..rings import DEFAULT_SIG, round_up_str


def to_extractors_magic_record(
    contract: MagicSourceContract,
    *,
    sig: int = DEFAULT_SIG,
) -> dict[str, Any]:
    """Emit a magic-source descriptor consumable alongside `ftqc_extractors` PBC / injection schedules."""
    return {
        "kind": "ftqc_magic_source",
        "name": contract.name,
        "target_state": contract.target.label,
        "k_outputs": contract.k_outputs,
        "fault_distance": contract.distance,
        "peak_active_qubits": contract.active_qubits,
        "conditional_infidelity_hi": round_up_str(contract.conditional_infidelity_upper, sig),
        "expected_wire_ticks_per_output_hi": round_up_str(contract.expected_wire_ticks_per_output_upper(), sig),
        "correlation_mode": contract.correlation_mode.value,
    }


def to_uncompute_pricing_kwargs(
    t_contract: MagicSourceContract,
    *,
    ccz_contract: MagicSourceContract | None = None,
    sig: int = DEFAULT_SIG,
) -> dict[str, Any]:
    """Emit certified keyword arguments directly compatible with
    `ftqc_uncompute.adapters.ftqc_bridge.price_circuit_with_zx_contract(...)`!
    """
    wt_frac = t_contract.expected_wire_ticks_per_output_upper()
    # Round UP to integer wire-ticks (`ceil(num / den)`) so `t_factory_wire_ticks: int` is an exact upper bound!
    wt_int_ceil = (wt_frac.numerator + wt_frac.denominator - 1) // wt_frac.denominator
    out: dict[str, Any] = {
        "t_factory_wire_ticks": wt_int_ceil,
        "t_factory_escape_hi": round_up_str(t_contract.conditional_infidelity_upper, sig),
        "magic_metadata": {
            "t_source": t_contract.name,
            "t_peak_active_qubits": t_contract.active_qubits,
            "t_correlation_mode": t_contract.correlation_mode.value,
        },
    }
    if ccz_contract is not None:
        ccz_wt = ccz_contract.expected_wire_ticks_per_output_upper()
        out["magic_metadata"]["ccz_source"] = ccz_contract.name
        out["magic_metadata"]["ccz_factory_wire_ticks_ceil"] = (
            ccz_wt.numerator + ccz_wt.denominator - 1
        ) // ccz_wt.denominator
        out["magic_metadata"]["ccz_factory_escape_hi"] = round_up_str(
            ccz_contract.conditional_infidelity_upper, sig
        )
    return out
