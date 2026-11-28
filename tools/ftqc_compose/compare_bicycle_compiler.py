"""Run the optional Qiskit Gross compiler on the fifteen BootLoops RM15 PPRs."""
import argparse
import json
import math
from pathlib import Path
import sys


def main():
    from bicycle_architecture_compiler import estimate
    from qiskit.circuit import QuantumCircuit
    from qiskit.circuit.library import PauliProductRotationGate
    from qiskit.quantum_info import Pauli
    sys.path.append(str(Path(__file__).resolve().parents[1]))
    from ftqc_extractors.bicycle import mapping
    from ftqc_gadgets.bicycle_ports import route
    from ftqc_magic.insitu import rm15_program

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", required=True, help="cache path for the external measurement table")
    args = parser.parse_args()
    program = rm15_program()
    n = len(program["prepare_plus"])
    circuit = QuantumCircuit(n)
    for row in program["rotations"]:
        label = "".join("Z" if row["parity"] >> q & 1 else "I" for q in reversed(range(n)))
        angle = -math.pi * row["exponent_i_pi_P_over_8"] / 4
        circuit.append(PauliProductRotationGate(Pauli(label), angle), range(n))
    external = estimate(circuit, "gross", "1e-3", "1e-10", measurement_table=args.table)
    moves = 0
    for row in mapping()["rotations"]:
        dx, dy = row["shift"]
        moves += len(route(-dx, -dy)) + len(route(dx, dy))
    print(json.dumps({"ppr_columns": len(program["rotations"]), "compiler_counts": external["counts"],
                      "native_gross_port_shift_and_return_generators": moves,
                      "native_shift_and_return_cx_ticks": 6 * moves}, indent=2))


if __name__ == "__main__":
    main()

