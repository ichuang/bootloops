# ftqc_zx — GUIDE

This engine is part of [ftqc_faults](../ftqc_faults/GUIDE.md);
overview of the FTQC tools: [`../FTQC.md`](../FTQC.md).

KIND: engine (Python standard library; no floating-point literals). The code
lives in `tools/ftqc_faults/zx/`: `gf2.py`, `diagram.py`, `pauliweb.py`,
`frame_oracle.py`, `equivalence.py`, `rewrites.py`, `synthesize.py`,
`decoder_synth.py`, `contracts.py`, `cli.py`, `adapters/` (optional `pyzx` and
`stim` bridges), `battery.py`, `selftest.py`, `examples/`. `tools/ftqc_zx/`
forwards to it. The canonical command is `python3 -m ftqc_faults zx`;
`python3 -m ftqc_zx` and `python3 tools/ftqc_zx/cli.py` still work. Battery
class: selftest.

## PURPOSE

Relations between Clifford circuits under faults, after the *Fault Tolerance
by Construction* framework of Rodatz, Poór and Kissinger [RPK25].

- **Fault equivalence.** C₁ and C₂ are w-fault-equivalent when every fault set
  of weight below w in one has the same effect as a fault set of no greater
  weight in the other. `equiv` decides this in both directions by enumeration,
  optionally on a CSS codespace.
- **Rewrite proofs.** A derivation by fault-equivalent ZX rewrites proves
  equivalence for families too large to enumerate. `verify-proof` checks a
  derivation in time linear in its size.
- **Synthesis.** Cat states (Props. 7.1–7.3), Shor syndrome measurements
  (7.5–7.7), Steane syndrome extractors including the 15-CNOT five-ancilla
  circuit (7.8–7.12), with their proofs, and flag-conditioned decoder tables
  (Appendix C.2).
- **Stochastic contracts** (`ftqc_zx-contract-v1`). For one scheduled gadget
  and a rational fault rate: the schedule (ticks, locations, wire-tick
  volume), an exact audit of every single fault (rejected, harmless, or bad),
  and outward-rounded bounds on acceptance and on escape. Every event with two
  or more faults is charged as a failure.

Every equivalence and distance verdict on a `CliffordCircuit` is computed
twice, by graphical Pauli webs (`pauliweb.py`) and by an independent forward
frame propagation (`frame_oracle.py`), and the two must agree.

## SCOPE

The Clifford fragment: spiders with phases kπ/2, Hadamard boxes,
preparations and Pauli measurements. Noise models `circuit_depolarizing`
(joint weight), `circuit_css` (X and Z sectors counted separately) and
`fault_free`, with per-wire and per-operation idealisations (`noisy`, `ft`,
`ff`). Non-Clifford phases and unsupported QASM are refused.

## USE-WHEN

- A compact or rescheduled Clifford gadget should be shown equivalent under
  faults to an idealised specification, with or without a codespace
  precondition.
- A gadget family is too large for enumeration but has a constructive ZX
  derivation.
- A flag gadget needs a certified decoder table, or a ledger needs a
  per-gadget acceptance and escape bound.

## NOT-FOR

- Static code distance: use `ftqc_codes`. A `CSSCodeSpec` enters here only as
  a precondition.
- One fixed circuit audited at higher weight with SAT/LRAT or against Stim:
  use the audit engine.
- Gadgets on unknown inputs with lookup decoders and composition receipts:
  use `ftqc_gadgets`.
- Non-Clifford protocols, HGP surgery and logical T: use `ftqc_extractors`.
- Unchecked `pyzx.full_reduce()` or `clifford_simp()`: these preserve the
  noiseless map but not fault equivalence. Check external rewrites with
  `audit_pyzx_rewrite`.

## INVOKE

From the BootLoops root with `export PYTHONPATH="$PWD/tools"`:

    python3 -m ftqc_faults zx selftest [--mutation-controls] [--json]
    python3 -m ftqc_faults zx equiv steane_5anc steane_spec --w 3 --noise-model circuit_css --code steane_7_1_3
    python3 -m ftqc_faults zx verify-proof cat --n 64 --w 3
    python3 -m ftqc_faults zx synth steane_5anc --out c.json --proof-out cert.json --stim-out c.stim --qasm-out c.qasm
    python3 -m ftqc_faults zx contract steane_5anc --p 1/10000 --idle-ratio 1/10 --noise-model circuit_css --code steane_7_1_3 --json
    python3 -m ftqc_faults zx audit tools/ftqc_faults/zx/examples/shor_opt_zzzz.json --w-max 3

Python API: `from ftqc_faults.zx import CliffordCircuit, ZXDiagram,
check_fault_equiv, circuit_distance, verify_derivation, synth_cat_state,
synth_shor_syndrome, synth_steane_syndrome, synthesize_flag_decoder,
build_gadget_contract, make_steane_7_1_3`.

Exit codes: 0 certified or pass; 1 a finding (failed equivalence, invalid
proof); 2 input or side-condition refusal (`RefusalError`).

## OUTPUTS

Verdict lines such as `[FAULT-EQUIV-PASS] w=3 ...` or `[PROOF-VALID] ...`, JSON
with `--json`, circuit files from `synth` (JSON, Stim, OpenQASM 2), and
`ftqc_zx-contract-v1` documents.

## GATES

`python3 tools/ftqc_zx/selftest.py`: 10 legs in 5 s on a four-core Linux
virtual machine, in a scratch directory.

- **L1** no `float` literals or calls in the package; JSON contracts reject
  float tokens.
- **L2** cat states (4, 5 and 8 qubits, full and w-truncated) and the 64-qubit
  proof check.
- **L3** Shor extraction (`standard`, three-ancilla `optimised`,
  rail-shifting `teleport`) on `ZZZZ`, `XXXX`, `XZYZ`, and the outcome map of
  Appendix A.
- **L4** Steane extraction (eight- and five-ancilla, X and Z sectors), closed
  distance 3, and the flag decoder.
- **L5** outward-rounding containment, contracts, the packaged examples,
  OpenQASM round trip, Stim `DETECTOR` emission, optional `pyzx`/`stim`
  probes, CLI exit codes.
- **N1–N5** an unflagged 4-cat fails; a bare CNOT is refused; Steane
  extraction without the encoder fails; invalid rewrites and non-Clifford QASM
  are refused.
- **M1–M4** four single-line sabotages in scratch copies (grouped Y = XZ
  weight, ±π/2 spider parity, the CNOT extraction guard, ceiling versus floor
  rounding) must each be caught.

## FOOTGUNS

- **`circuit_css` versus `circuit_depolarizing`.** The bent-leg Steane
  constructions (Props. 7.10–7.12) are 3-fault-equivalent under `circuit_css`
  and not under joint counting, where a CNOT fault with X on one qubit and Z
  on the other has support two. Receipts record `noise_model` and
  `bad_css_sector_labels`.
- **Open boundary wires have circuit distance 1** (Definition 2.8): a fault on
  an open input or output is undetectable alone. Compare against an idealised
  specification, or close the circuit with an idealised verifier.
- **Steane extraction needs a codespace.** It is not fault-equivalent on
  unencoded inputs; pass `--code steane_7_1_3`.
- **Switching CSS sectors.** Do not wrap noisy Hadamards around data wires;
  `synth_steane_syndrome(..., sector="Z")` emits the colour-dual circuit.
- **Contracts are pessimistic above one fault.** The escape bound adds the
  probability of a bad single fault (zero when the audit passes) to the
  probability of two or more faults. It is a sound upper bound and loose
  whenever multi-fault events are mostly harmless.
- `equiv` enumerates fault classes up to a state cap and refuses with
  `FAULT-ENUM-BUDGET-EXCEEDED` beyond it.

## KNOWN GAPS

- The rewrite rules and propositions are implemented from the paper; the
  rules themselves are not machine-checked here.
- Context-dependent single-ancilla surface-code hook schedules (Remark 3.6)
  need their spacetime context; audit the full multi-round patch with the
  audit engine instead.

CREDIT: implements the fault-equivalence, fault-gadget, fault-equivalent ZX
rewrite, cat-state, Shor and five-ancilla Steane constructions of Rodatz, Poór
and Kissinger [RPK25] (building on distance-preserving rewrites [RPK24], the ZX
calculus of Coecke and Duncan [CD08], Pauli webs of Bombin et al. [BLNP24] and
Borghans [Bor19], spacetime codes and spackle maps of Bacon et al. [BFHS17] and
Delfosse and Paetznick [DP23], Goto's Steane state preparation [Goto16], and
Gidney's Stim detector model [Gid21]). The two-lineage verifier, the
flag-decoder compiler and the contract layer are ours; the `pyzx` and `stim`
bridges call upstream packages only when installed. Bracketed keys resolve in
`REFERENCES.md` at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz;
code written by Claude (Anthropic) under his supervision. Documentation is
licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository
root.
