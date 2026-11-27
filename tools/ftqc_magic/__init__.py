"""ftqc_magic: certified magic-state distillation factories, cultivation, and pipelines.

100% Python standard library, MIT-licensed, zero floating-point literals or float()/complex() calls.
"""

from .catalog import CATALOG_NAMES, get_catalog_factory, list_catalog
from .galois import (
    GF2Ext,
    GaloisFactorySpec,
    build_gf4_hermitian_ag_cs_factory,
    build_gf8_reed_solomon_ccz_factory,
    build_gprw_64t_2ccz_trace_factory,
    compare_global_vs_local_concatenation,
    three_twirl_ccz_weight_transfer,
)
from .ir import (
    CorrelationMode,
    MagicRefusalError,
    MagicSourceContract,
    PhasePolynomialTarget,
)
from .limits import (
    SignedStabilizerFactory,
    analytical_lower_bound_n,
    build_bgm_3t_to_cs_factory,
    build_bgm_6t_to_cccz_factory,
    build_bgm_9t_to_1t_factory,
    stabilizer_nullity,
)
from .matrix import (
    BinaryFactory,
    from_borrowed_identity,
    free_row_append,
    parent_capacity,
    parity_lift,
    puncture_base_code,
    verify_puncture_downward_closure,
)
from .rings import PolyQ, QSqrt2, ZOmega8, ZSqrt2, round_down_str, round_up_str

VERSION = "ftqc_magic-receipt-1"

__all__ = [
    "CATALOG_NAMES",
    "GF2Ext",
    "GaloisFactorySpec",
    "VERSION",
    "BinaryFactory",
    "CorrelationMode",
    "MagicRefusalError",
    "MagicSourceContract",
    "PhasePolynomialTarget",
    "PolyQ",
    "QSqrt2",
    "SignedStabilizerFactory",
    "ZOmega8",
    "ZSqrt2",
    "analytical_lower_bound_n",
    "build_bgm_3t_to_cs_factory",
    "build_bgm_6t_to_cccz_factory",
    "build_bgm_9t_to_1t_factory",
    "build_gf4_hermitian_ag_cs_factory",
    "build_gf8_reed_solomon_ccz_factory",
    "build_gprw_64t_2ccz_trace_factory",
    "compare_global_vs_local_concatenation",
    "free_row_append",
    "from_borrowed_identity",
    "get_catalog_factory",
    "list_catalog",
    "parent_capacity",
    "parity_lift",
    "puncture_base_code",
    "round_down_str",
    "round_up_str",
    "stabilizer_nullity",
    "three_twirl_ccz_weight_transfer",
    "verify_puncture_downward_closure",
]
