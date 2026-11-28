# ftqc_faults

Physical faults in scheduled Clifford circuits. The audit engine
(`ftqc_faults audit`) traces every location fault through one circuit and
reports circuit and fault distance with replayable witnesses. The ZX engine
(`ftqc_faults zx`) checks fault equivalence between two circuits, verifies ZX
rewrite proofs, synthesizes gadgets and writes stochastic gadget contracts.
`ftqc_faults steane` runs both engines on the five-ancilla Steane extractor
and checks that they agree on every single fault.

    export PYTHONPATH="$PWD/tools"                       # from the BootLoops root
    python3 -m ftqc_faults audit analyze surface:3,2,z,bad
    python3 -m ftqc_faults zx equiv steane_5anc steane_spec --w 3 --noise-model circuit_css --code steane_7_1_3
    python3 -m ftqc_faults steane -o steane.json && python3 -m ftqc_faults steane --verify steane.json

[Guide](GUIDE.md) · [Audit engine](../ftqc_audit/GUIDE.md) ·
[ZX engine](../ftqc_zx/GUIDE.md) · [FTQC overview](../FTQC.md)
