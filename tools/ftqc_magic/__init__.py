"""ftqc_magic: certified magic-state distillation factories, cultivation, and pipelines.

100% Python standard library, MIT-licensed, zero floating-point literals or float()/complex() calls.
"""

from .catalog import CATALOG_NAMES, get_catalog_factory, list_catalog
from .cultivate import (
    CliffordTableau,
    analyze_cultivation_double_check,
    build_css_encoder_pair,
    build_itogawa_short_check_schedule,
    build_msc3_double_check_schedule,
    clifford_error_acceptance_prob,
    colour_code_7_1_3_data,
    colour_code_t_conditional_fidelity,
    propagate_pauli_across_t_layer,
)
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
from .recycle import (
    RecycledSchedule,
    certify_active_qubits_sat,
    compress_active_qubits,
    evaluate_active_trace,
    minimize_generator_weights,
)
from .rings import PolyQ, QSqrt2, ZOmega8, ZSqrt2, round_down_str, round_up_str
from .search import (
    SpatialSATInstance,
    build_canonical_t_factory,
    canonical_sqrt_t_theoretical_params,
    canonical_t_theoretical_params,
    encode_spatial_nogo_cnf,
    enumerate_jws_symmetry_orbits,
    search_symmetric_sn_family,
    solve_or_verify_spatial_sat,
    verify_banked_spatial_nogo,
)

VERSION = "ftqc_magic-receipt-1"

__all__ = [
    "CATALOG_NAMES",
    "CliffordTableau",
    "GF2Ext",
    "GaloisFactorySpec",
    "RecycledSchedule",
    "SpatialSATInstance",
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
    "analyze_cultivation_double_check",
    "build_bgm_3t_to_cs_factory",
    "build_bgm_6t_to_cccz_factory",
    "build_bgm_9t_to_1t_factory",
    "build_canonical_t_factory",
    "build_css_encoder_pair",
    "build_gf4_hermitian_ag_cs_factory",
    "build_gf8_reed_solomon_ccz_factory",
    "build_gprw_64t_2ccz_trace_factory",
    "build_itogawa_short_check_schedule",
    "build_msc3_double_check_schedule",
    "canonical_sqrt_t_theoretical_params",
    "canonical_t_theoretical_params",
    "certify_active_qubits_sat",
    "clifford_error_acceptance_prob",
    "colour_code_7_1_3_data",
    "colour_code_t_conditional_fidelity",
    "compare_global_vs_local_concatenation",
    "compress_active_qubits",
    "encode_spatial_nogo_cnf",
    "enumerate_jws_symmetry_orbits",
    "evaluate_active_trace",
    "free_row_append",
    "from_borrowed_identity",
    "get_catalog_factory",
    "list_catalog",
    "minimize_generator_weights",
    "parent_capacity",
    "parity_lift",
    "propagate_pauli_across_t_layer",
    "puncture_base_code",
    "round_down_str",
    "round_up_str",
    "search_symmetric_sn_family",
    "solve_or_verify_spatial_sat",
    "stabilizer_nullity",
    "three_twirl_ccz_weight_transfer",
    "verify_banked_spatial_nogo",
    "verify_puncture_downward_closure",
]
