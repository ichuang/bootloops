"""faultweb — ZX-calculus fault equivalence, Pauli-web analysis, and FTQC gadget certifier."""

from .diagram import (
    AtomicFaultGroup,
    CircuitOp,
    CliffordCircuit,
    RefusalError,
    ZXDiagram,
    ZXEdge,
    ZXNode,
    circuit_to_zx_fe,
    extract_circuit_fe,
    load_spec_file,
)
from .frame_oracle import SpackleOracleResult, analyze_circuit_spackle
from .pauliweb import (
    FaultSignatureTables,
    OutcomePushOutMap,
    PauliWebAnalysis,
    analyze_pauli_webs,
    build_fault_signature_tables,
    compile_outcome_map,
    match_measurement_branches,
)

__all__ = [
    "AtomicFaultGroup",
    "CircuitOp",
    "CliffordCircuit",
    "FaultSignatureTables",
    "OutcomePushOutMap",
    "PauliWebAnalysis",
    "RefusalError",
    "SpackleOracleResult",
    "ZXDiagram",
    "ZXEdge",
    "ZXNode",
    "analyze_circuit_spackle",
    "analyze_pauli_webs",
    "build_fault_signature_tables",
    "circuit_to_zx_fe",
    "compile_outcome_map",
    "extract_circuit_fe",
    "load_spec_file",
    "match_measurement_branches",
]
