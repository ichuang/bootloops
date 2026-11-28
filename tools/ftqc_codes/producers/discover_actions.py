"""Optional qLDPC/GAP or AutQEC/Bliss discovery. Run with PYTHONPATH unset.

Write one independently replayable native receipt per generator. The external
programs discover candidates. Only successful signed BootLoops replay is saved.
"""
import argparse
import json
import sys
from pathlib import Path


def qldpc_circuits(matrix, gates):
    from qldpc.codes import QuditCode
    from qldpc.circuits import get_transversal_ops
    for _, circuit in get_transversal_ops(QuditCode(matrix), local_gates=gates, remove_redundancies=False):
        yield circuit


def autqec_circuits(matrix):
    import numpy as np
    import stim
    from autqec.automorphisms import circ_from_aut, logical_circ_and_pauli_correct
    from autqec.graph_auts import valid_clifford_auts_B_rows
    from autqec.utils.qiskit import construct_circuit
    from qiskit import transpile

    n = matrix.shape[1] // 2
    # For small codes give Bliss the full unsigned stabilizer set so its answer
    # is independent of the initial check generating set. Bound exponential growth.
    rows = matrix
    if len(matrix) <= 6:
        selectors = np.arange(1, 1 << len(matrix))[:, None] >> np.arange(len(matrix)) & 1
        rows = (selectors @ matrix) % 2
    x, z = rows[:, :n], rows[:, n:]
    # circ_from_aut numbers each physical site's triplet (Z, X, X+Z).
    triples = np.column_stack((z, x, x ^ z))
    for generator in valid_clifford_auts_B_rows(triples, return_order=False):
        physical, _ = circ_from_aut(matrix, generator).circ()
        _, corrected = logical_circ_and_pauli_correct(matrix, physical).run()
        gates = transpile(construct_circuit(corrected, n),
                          basis_gates=["h", "s", "sdg", "x", "y", "z", "cx", "cz", "swap"],
                          optimization_level=0)
        circuit = stim.Circuit()
        for inst in gates.data:
            name = inst.operation.name
            circuit.append("S_DAG" if name == "sdg" else name.upper(),
                           [gates.find_bit(q).index for q in inst.qubits])
        yield circuit


def steps_of(circuit):
    import stim
    steps = []
    for inst in circuit.flattened():
        if inst.name in ("I", "TICK"):
            continue
        targets = [t.value for t in inst.targets_copy()]
        gate = "SDG" if inst.name == "S_DAG" else inst.name
        if gate not in {"H", "S", "SDG", "X", "Y", "Z", "CX", "CZ", "SWAP"}:
            # Stim synthesizes this Clifford with elementary H, S and CX gates.
            part = stim.Circuit()
            part.append(inst)
            steps.extend(steps_of(part.to_tableau().to_circuit()))
            continue
        width = 2 if gate in ("CX", "CZ", "SWAP") else 1
        for i in range(0, len(targets), width):
            steps.append({"layer": [[gate, *targets[i:i + width]]]})
    return steps


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("backend", choices=["qldpc", "autqec"])
    ap.add_argument("code")
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--gates", nargs="+", choices=["SWAP", "H", "S", "SQRT_X"],
                    help="restrict qLDPC discovery (default H S SWAP)")
    args = ap.parse_args()
    if args.gates and args.backend != "qldpc":
        ap.error("--gates applies to qLDPC")
    import numpy as np
    import stim
    if args.backend == "qldpc":
        import qldpc  # load PyPI galois before adding the BootLoops tools path
    sys.path.append(str(Path(__file__).resolve().parents[2]))
    from ftqc_codes import action_domain, automorphisms as aut, basis
    code = action_domain.load(args.code)
    paired = action_domain.default_basis(code)
    checks = basis._code_signed_checks(code)
    independent = basis._SignedStabilizerGroup(code.n)
    rows = []
    for x, z, sign in checks:
        if independent.add((x, z, sign)):
            rows.append([(v >> j) & 1 for v in (x, z) for j in range(code.n)])
    matrix = np.array(rows, dtype=int).reshape(-1, 2 * code.n)
    circuits = (qldpc_circuits(matrix, args.gates or ["H", "S", "SWAP"])
                if args.backend == "qldpc" else autqec_circuits(matrix))
    receipts = []
    for circuit in circuits:
        steps = steps_of(circuit)
        rec = aut.build(code, paired, steps)
        for q in range(code.n):
            for x, z in ((1 << q, 0), (0, 1 << q)):
                inp = basis.pauli_from_xz(code.n, x, z)
                out = basis.pauli_from_xz(code.n, *aut.conjugate((x, z, 0), steps))
                if stim.PauliString(inp).after(circuit) != stim.PauliString(out):
                    raise ValueError("Stim and native Clifford transport disagree")
        receipts.append(rec)
    args.outdir.mkdir(parents=True, exist_ok=True)
    for i, rec in enumerate(receipts):
        path = args.outdir / f"{args.backend}-{i}.action.json"
        path.write_text(json.dumps(rec, indent=2) + "\n")
    print(f"{len(receipts)} signed generators verified")


if __name__ == "__main__":
    main()

