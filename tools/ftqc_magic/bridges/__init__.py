"""Read-only bridges from `ftqc_magic` to sibling BootLoops FTQC tools.

Note: `tools/ftqc_compose/` and `tools/ftqc_uncompute/` are kept 100% untouched.
"""

from .compose_bridge import to_compose_ancilla_contract, to_compose_logical_t_summary
from .extractors_bridge import to_extractors_magic_record, to_uncompute_pricing_kwargs
from .faults_bridge import audit_ppr_check_frame, verify_zx_clifford_frame

__all__ = [
    "audit_ppr_check_frame",
    "to_compose_ancilla_contract",
    "to_compose_logical_t_summary",
    "to_extractors_magic_record",
    "to_uncompute_pricing_kwargs",
    "verify_zx_clifford_frame",
]
