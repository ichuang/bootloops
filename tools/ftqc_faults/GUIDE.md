# Circuit fault analysis and ZX equivalence

From the BootLoops root, export `PYTHONPATH="$PWD/tools"`. The audit engine accepts a scheduled circuit with preparations, gates, idles, measurements and protected outputs. Its unit of fault weight is a faulty physical location. The ZX engine accepts Clifford circuits and ZX diagrams, with explicit boundary, branch and noise conventions. Its equivalence checker compares two circuits. Its stochastic contract checker uses a scheduled physical circuit. [Audit engine](../ftqc_audit/GUIDE.md) · [ZX engine](../ftqc_zx/GUIDE.md)

## Five-ancilla Steane extractor

```sh
export PYTHONPATH="$PWD/tools"
out=$(mktemp -d)
python3 -m ftqc_faults zx equiv steane_5anc steane_spec --w 3 --noise-model circuit_css --code steane_7_1_3
python3 -m ftqc_faults steane -o "$out/steane.json"
python3 -m ftqc_faults steane --verify "$out/steane.json"
```

For each of the X and Z syntheses, the bridge checks the ZX rewrite proof and bidirectional equivalence on the Steane codespace, with the independent ZX frame oracle enabled. The `w=3` comparison uses the ZX CSS-sector convention and tests classes of weight strictly below three.

The separate physical comparison takes the synthesized 15-CNOT circuit and its ten-tick ASAP schedule. Every live idle is explicit. The audit probes begin with an ideal encoded zero and end with ideal measurements of all seven data qubits, once in Z and once in X. They declare the flag and three syndrome parities as detectors. The comparison retains all five raw record displacements before recovery, relative to the same fault-free branch. Audit numbers measurements in tick order, while ZX uses the original circuit order. The bridge matches the labeled records before comparing them.

Each sector has 95 noisy locations: five preparations, 15 CNOTs, five measurements and 70 idles. The physical model tests all 15 joint nonidentity Paulis after each CNOT, all three Paulis after each preparation or idle, and one classical outcome flip per measurement. For each of the 455 modes, the ZX scheduled propagator must agree with the audit adjoint signatures and forward replay on the four detectors, all raw record flips, and the X and Z masks of the seven output qubits. The receipt binds both probe circuits and the ordered response ledger by digest. A changed result fails replay.

The physical claim covers single faults with the stated boundaries and frame. A circuit-distance or multi-fault claim requires its own test. The eight established suites retain their registry names, including `ftqc_audit` and `ftqc_zx`; the additional `ftqc_faults` battery exercises the bridge and compatibility paths.

## Existing clients and receipts

The canonical imports are `ftqc_faults.audit` and `ftqc_faults.zx`. Existing `ftqc_audit` and `ftqc_zx` imports and command lines continue to work. Circuit receipts and stochastic contracts retain their original schemas and independent verifiers. Published injection cases and decoder sampling belong to `ftqc_extractors`. Run `python3 -m ftqc_extractors injection bb-uniform` for the closed Y experiment. The old audit injection and sampling commands forward to the same implementation. [Injection guide](../ftqc_extractors/INJECTION.md)
