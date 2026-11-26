"""ftqc_uncompute — Measurement-based uncomputation, spooky pebbling, and superposition masking."""

from .gf2 import GF2Field, anf_and, anf_const, anf_degree, anf_eval, anf_not, anf_var, anf_xor
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

__all__ = [
    "ArithDAG",
    "Bit",
    "CircuitBuilder",
    "DAGNode",
    "GF2Field",
    "KickMixCircuit",
    "KickMixOp",
    "MaskDirective",
    "OperationStatistics",
    "Qubit",
    "RegDecl",
    "RegisterSpec",
    "anf_and",
    "anf_const",
    "anf_degree",
    "anf_eval",
    "anf_not",
    "anf_var",
    "anf_xor",
    "lift_scoped_kmx_to_dag",
]
