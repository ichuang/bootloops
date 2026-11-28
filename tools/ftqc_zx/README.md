# ftqc_zx

The ZX fault-equivalence engine of [ftqc_faults](../ftqc_faults/GUIDE.md),
after Rodatz, Poór and Kissinger. Its code is in `tools/ftqc_faults/zx/`; this
directory provides the import name `ftqc_zx` and its command line.

It decides w-fault equivalence between Clifford circuits (two independent
computations must agree), checks constructive ZX rewrite proofs, synthesizes
cat-state, Shor and Steane gadgets with their proofs and flag decoder tables,
and writes exact-rational stochastic contracts for one scheduled gadget.

    export PYTHONPATH="$PWD/tools"                       # from the BootLoops root
    python3 -m ftqc_faults zx equiv steane_5anc steane_spec --w 3 --noise-model circuit_css --code steane_7_1_3
    python3 -m ftqc_faults zx verify-proof cat --n 64 --w 3
    python3 -m ftqc_faults zx contract steane_5anc --p 1/10000 --idle-ratio 1/10 --noise-model circuit_css --code steane_7_1_3 --json
    python3 tools/ftqc_zx/selftest.py                    # battery

## Package layout (`tools/ftqc_faults/zx/`)

- `gf2.py`: bit-packed GF(2) row reduction, null spaces, linear solves, signed
  Pauli row reduction.
- `diagram.py`: `CliffordCircuit` and `ZXDiagram`, idealisations
  (`noisy`/`ft`/`ff`), fault gadgets, guarded circuit extraction.
- `pauliweb.py`: first computation, Pauli-web detectors and boundary
  stabilizers, and the outcome map of Appendix A.
- `frame_oracle.py`: second computation, stabilizer dilation and spacetime
  propagation.
- `equivalence.py`: w-fault equivalence, circuit distance, codespace
  preconditions (`steane_7_1_3`, `shor_9_1_3`, `golay_23_1_7`).
- `rewrites.py`: fault-equivalent rewrites (Props. 6.12–6.16, B.5–B.8) and
  `verify_derivation`.
- `synthesize.py`: cat (Props. 7.1–7.3), Shor (7.5–7.7) and Steane (7.8–7.12)
  synthesizers and planted negative controls.
- `decoder_synth.py`: flag-conditioned decoder tables (Appendix C.2).
- `contracts.py`: tick packing, single-fault audit and `ftqc_zx-contract-v1`.
- `adapters/`: OpenQASM 2 import and export, Stim emission, optional `pyzx`
  and `stim` bridges.
- `battery.py`, `selftest.py`: the 10-leg battery.

See [`GUIDE.md`](GUIDE.md) for scope, gates and hazards, and
[`../FTQC.md`](../FTQC.md) for how this engine relates to the other tools.
