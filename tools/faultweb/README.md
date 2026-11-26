# faultweb

Exact $\mathbb{F}_2$ Pauli-web verifier, constructive fault-equivalent ZX-calculus
compiler, flag-conditioned syndrome decoder synthesizer, and zero-float outward-rounded
stochastic gadget-contract emitter (`faultweb-contract-v1`) for fault-tolerant
quantum computing (FTQC).

See [`GUIDE.md`](GUIDE.md) for scope, acceptance gates, footguns, and full credit
keys.

## Quick start

```bash
# Run self-contained 10-leg battery (including 4 scratch mutation controls)
python3 tools/faultweb/selftest.py

# Verify 3-fault equivalence of Kissinger's 5-ancilla, 15-CNOT Steane extractor
python3 tools/faultweb/cli.py equiv steane_5anc steane_spec \
  --w 3 --noise-model circuit_css --code steane_7_1_3

# Verify constructive ZX proof for a 64-qubit w=3 cat state in polynomial time
python3 tools/faultweb/cli.py verify-proof cat --n 64 --w 3

# Emit a zero-float outward-rounded stochastic ledger contract at p = 1e-4
python3 tools/faultweb/cli.py contract steane_5anc \
  --p 1e-4 --idle-ratio 1/10 --noise-model circuit_css --code steane_7_1_3 --json
```

## Package layout

- `gf2.py` — bit-packed `int` $\mathbb{F}_2$ RREF, nullspace, linear solver, and signed Pauli-phase RREF.
- `diagram.py` — `CliffordCircuit` and `ZXDiagram` IRs, idealisations (`noisy`/`ft`/`ff`), fault gadgets, and fail-closed `extract_circuit_fe`.
- `pauliweb.py` — Lineage A: $\mathbb{F}_2$ Pauli-web detector / boundary-stabiliser solver and `Appendix A` outcome push-out compiler.
- `frame_oracle.py` — Lineage B: disjoint Stinespring/Choi unitary stabiliser dilation and spacetime spackle oracle.
- `equivalence.py` — semantic $w$-fault equivalence, circuit distance, and CSS codespace preconditioning (`steane_7_1_3`, `shor_9_1_3`, `golay_23_1_7`).
- `rewrites.py` — fault-equivalent ZX rewrites (`Props. 6.12-6.16, B.5-B.8`) and polynomial-time `verify_derivation`.
- `synthesize.py` — cat (`Props. 7.1-7.3`), Shor (`Props. 7.5-7.7`), and Steane (`Props. 7.8-7.12`) synthesizers + planted negative controls.
- `decoder_synth.py` — automated `Appendix C.2` punctured check matrix $H'$ and flag-conditioned decoder lookup compiler.
- `contracts.py` — ASAP tick packer, single-fault coordinate/filter audit, and zero-float outward-rounded `faultweb-contract-v1` emitter.
- `adapters/` — OpenQASM 2.0 importer/exporter, `.stim` + `DETECTOR` emitter, and optional runtime `pyzx` / `stim` bridges.
- `battery.py` / `selftest.py` — 10-leg synthetic + 4-mutant acceptance battery (`selftest` verification class).
