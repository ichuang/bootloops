# ftqc_zx

Exact $\mathbb{F}_2$ Pauli-web verifier, constructive fault-equivalent ZX-calculus
compiler, flag-conditioned syndrome decoder synthesizer, and zero-float outward-rounded
stochastic gadget-contract emitter (`ftqc_zx-contract-v1`) for fault-tolerant
quantum computing (FTQC).

See [`GUIDE.md`](GUIDE.md) for scope, acceptance gates, footguns, sibling-tool
comparisons, and full credit keys.

## Choosing between `ftqc_zx`, `ftqc_codes`, `ftqc_audit`, `ftqc_gadgets`, and `ftqc_extractors`

| Package | What it certifies | Use when | Do **not** use for |
|---|---|---|---|
| **`ftqc_codes`** | Static CSS code parameters $[n,k,d_X,d_Z,t,w_{\min}]$ via enumeration, CaDiCaL/LRAT, Lean-QEC, BCH, cuts, and concatenation | You have check matrices $(H_X, H_Z)$ and need a machine-checkable proof of **code** distance | Circuit hook errors, schedules, flags, or gates |
| **`ftqc_zx`** (this package) | Relational $w$-fault equivalence ($C_1 \stackrel{w}{\widehat{=}} C_2$), polynomial-time ZX rewrite proofs, ancilla-reduced cat/Shor/Steane synthesis, flag-decoder tables, and zero-float outward-rounded stochastic contracts (`ftqc_zx-contract-v1`) | Designing/optimising Clifford gadgets *by construction* via ZX rewrites, quotienting classical Pauli-frame byproducts (`Appendix A`), certifying `circuit_css` + flag decoders (`Appendix C.2`), or bounding all-orders stochastic escape/dirt rates | Static code-distance proofs (`ftqc_codes`), unstructured high-weight SAT/LRAT searches (`ftqc_audit`), or non-Clifford magic-state logical $T$ / HGP surgery (`ftqc_extractors`) |
| **`ftqc_audit`** | Single-circuit adversarial fault model, circuit/fault distance, strict FT order, and location census via meet-in-the-middle, CaDiCaL/LRAT, and Stim DEM | Auditing one hand-written Clifford tick circuit (e.g. rotated surface-code CX schedules) with LRAT proofs or Stim DEM diffs | Relational equivalence between open circuits, ZX synthesis, `circuit_css` sector mode, or all-orders probability enclosures |
| **`ftqc_gadgets`** | Aliferis-Gottesman-Preskill `prep`/`gate`/`ec`/`close` properties on unknown inputs, truncated refinement, and composition receipts | Proving AGP rectangle properties or composing sub-gadget refinement receipts with lookup decoders at small fault order | Constructive polynomial-time ZX rewrite derivations or all-orders stochastic `Decimal` contracts |
| **`ftqc_extractors`** | Hypergraph-product (HGP) static/scheduled logical extractors and non-Clifford magic-state logical $T$ protocols ($\mathbb{Z}_8$ branch checks + Clifford arm bounds) | Proving HGP graph-surgery logical measurements, spectator preservation, or conditional $S^p Z^x$ logical $T$ channels | General ZX circuit synthesis/leg-bending, flag-decoder table compilation, or outward-rounded stochastic gadget contracts |

## Quick start

```bash
# Run self-contained 10-leg battery (including 4 scratch mutation controls)
python3 tools/ftqc_zx/selftest.py

# Verify 3-fault equivalence of Kissinger's 5-ancilla, 15-CNOT Steane extractor
python3 tools/ftqc_zx/cli.py equiv steane_5anc steane_spec \
  --w 3 --noise-model circuit_css --code steane_7_1_3

# Verify constructive ZX proof for a 64-qubit w=3 cat state in polynomial time
python3 tools/ftqc_zx/cli.py verify-proof cat --n 64 --w 3

# Emit a zero-float outward-rounded stochastic ledger contract at p = 1e-4
python3 tools/ftqc_zx/cli.py contract steane_5anc \
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
- `contracts.py` — ASAP tick packer, single-fault coordinate/filter audit, and zero-float outward-rounded `ftqc_zx-contract-v1` emitter.
- `adapters/` — OpenQASM 2.0 importer/exporter, `.stim` + `DETECTOR` emitter, and optional runtime `pyzx` / `stim` bridges.
- `battery.py` / `selftest.py` — 10-leg synthetic + 4-mutant acceptance battery (`selftest` verification class).
