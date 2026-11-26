"""ftqc_uncompute — Measurement-based uncomputation, spooky pebbling, and superposition masking."""

from .gf2 import GF2Field, anf_and, anf_const, anf_degree, anf_eval, anf_not, anf_var, anf_xor
from .hmr import HMRFinding, HMRVerificationResult, verify_hmr_circuit
from .ir_dag import ArithDAG, DAGNode, RegDecl, lift_scoped_kmx_to_dag
from .ir_kickmix import (
    Bit,
    CircuitBuilder,
    KickMixCircuit,
    KickMixOp,
    MaskDirective,
    OperationStatistics,
    Qubit,
    RegisterSpec,
)
from .kickmix_sim import (
    FiatShamirPRNG,
    TrajectoryResult,
    clopper_pearson_upper_bound,
    coherent_trace_distance_bound,
    dqi_repetition_factor,
    fiat_shamir_seed,
    simulate_trajectories,
    z_pow_synthesis_cost,
)
from .kraus_oracle import (
    KrausVerificationResult,
    evolve_kraus_ensemble,
    exact_ensemble_trace_distance,
    verify_channel_kraus,
)

__all__ = [
    "ArithDAG",
    "Bit",
    "CircuitBuilder",
    "DAGNode",
    "FiatShamirPRNG",
    "GF2Field",
    "HMRFinding",
    "HMRVerificationResult",
    "KickMixCircuit",
    "KickMixOp",
    "KrausVerificationResult",
    "MaskDirective",
    "OperationStatistics",
    "Qubit",
    "RegDecl",
    "RegisterSpec",
    "TrajectoryResult",
    "anf_and",
    "anf_const",
    "anf_degree",
    "anf_eval",
    "anf_not",
    "anf_var",
    "anf_xor",
    "clopper_pearson_upper_bound",
    "coherent_trace_distance_bound",
    "dqi_repetition_factor",
    "evolve_kraus_ensemble",
    "exact_ensemble_trace_distance",
    "fiat_shamir_seed",
    "lift_scoped_kmx_to_dag",
    "simulate_trajectories",
    "verify_channel_kraus",
    "verify_hmr_circuit",
    "z_pow_synthesis_cost",
]
