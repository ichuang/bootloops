# faultweb — GUIDE

Tool page: https://bootloops.ai/tools/index.html (registry index)

KIND: package (`gf2.py`, `diagram.py`, `pauliweb.py`, `frame_oracle.py`,
`equivalence.py`, `rewrites.py`, `synthesize.py`, `decoder_synth.py`,
`contracts.py`, `cli.py`, `adapters/{pyzx_bridge,stim_bridge}.py`, plus the
in-tree synthetic + mutation battery `battery.py` / `selftest.py` and
`examples/`; 100% Python stdlib core, zero `float` literals)

## PURPOSE

Exact $\mathbb{F}_2$ Pauli-web verifier, constructive fault-equivalent ZX-calculus
compiler, flag-conditioned syndrome decoder synthesizer, and zero-float outward-rounded
stochastic gadget-contract emitter (`faultweb-contract-v1`) for Clifford
fault-tolerant quantum computing (FTQC) gadgets, implementing and extending the
*Fault Tolerance by Construction* framework of Rodatz, Poór, and Kissinger [RPK25].

## SCOPE

`faultweb` operates strictly on the Clifford fragment of quantum circuits and ZX
diagrams (spiders with phases $k \frac{\pi}{2}, k \in \{0, 1, 2, 3\}$, Hadamard
boxes, preparations, and Pauli measurements) under explicit adversarial and
stochastic circuit-level noise models (`circuit_depolarizing`, `circuit_css`,
`fault_free`, plus per-wire / per-op idealisations `noisy`, `ft`, `ff`). Every
direct `CliffordCircuit` equivalence or distance verdict is dual-lineage
cross-checked between Lineage A (`pauliweb.py`: graphical $\mathbb{F}_2$ Pauli
webs and fault gadgets) and Lineage B (`frame_oracle.py`: Stinespring/Choi
unitary stabiliser dilation and forward spacetime spackle propagation). Every
stochastic probability bound in `contracts.py` is computed from exact
`fractions.Fraction` polynomials and outward-rounded with `decimal.Decimal`
(`ROUND_FLOOR` for lower bounds, `ROUND_CEILING` for upper bounds) with a
runtime containment check (`Fraction(lo) <= exact <= Fraction(hi)`).

## USE-WHEN

- Certifying or refuting semantic $w$-fault equivalence ($C_1 \stackrel{w}{\widehat{=}} C_2$)
  between an implementable Clifford gadget (state preparation, flag gadget,
  syndrome measurement, teleportation routing, or full syndrome refresh) and an
  idealised specification, with or without CSS codespace preconditioning.
- Compiling or verifying constructive polynomial-time ZX fault-equivalence
  derivations (`verify_derivation`) at scales where direct adversarial fault
  enumeration is NP-hard (e.g. 64-qubit $w$-fault-tolerant cat states in < 0.5 s).
- Synthesizing fault-tolerant cat states (`Prop. 7.1-7.3`), Shor syndrome
  measurements (`standard`, 3-ancilla `optimised`, rail-shifting `teleport`;
  `Prop. 7.5-7.7`), Steane syndrome extractors (`standard_8anc`, bent-leg
  15-CNOT `optimised_5anc`; `Prop. 7.8-7.12`), and certified flag-conditioned
  syndrome lookup tables (`synthesize_flag_decoder`; `Appendix C.2`).
- Emitting certified zero-float outward-rounded stochastic contracts
  (`faultweb-contract-v1`) with ASAP tick schedules, wire-tick spacetime volume,
  per-coordinate active fault-location counts, and conditional escape upper
  bounds for whole-circuit FTQC overhead and fault-masking ledgers.

## NOT-FOR

- Non-Clifford gates ($T$, $\text{CCZ}$, arbitrary-angle $R_Z(\theta)$):
  non-Clifford phases and QASM instructions are refused fail-closed with
  `NON-CLIFFORD-PHASE-REFUSED` / `NON-CLIFFORD-OR-UNSUPPORTED-QASM` (magic-state
  distillation and Clifford+$T$ gadgets enter via their Clifford verification
  wrappers or external algebraic verifiers).
- Brute-force semantic enumeration (`check_fault_equiv` / `circuit_distance`)
  at large weight on large codes: semantic distance/equivalence checking is
  NP-hard (`Theorem 3.11`); the reachable syndrome BFS enforces a deterministic
  cap (`FAULT-ENUM-BUDGET-EXCEEDED`) and directs large-scale callers to
  polynomial-time proof verification (`verify_derivation`).
- Unchecked `pyzx.full_reduce()` / `pyzx.clifford_simp()` optimisation: standard
  ZX simplification preserves only the noise-free linear map and destroys fault
  equivalence (`§5`); always verify external graph rewrites with
  `audit_pyzx_rewrite`.
- Context-dependent surface-code single-ancilla hook schedules (`Remark 3.6`)
  without supplying the surrounding spacetime context.

## INVOKE

- CLI (from repo root with `PYTHONPATH=tools` via `python3 -m faultweb`, or direct script `python3 tools/faultweb/cli.py`):
  - `python3 tools/faultweb/cli.py selftest [--mutation-controls] [--json]`
  - `python3 tools/faultweb/cli.py equiv steane_5anc steane_spec --w 3 --noise-model circuit_css --code steane_7_1_3 [--json]`
  - `python3 tools/faultweb/cli.py verify-proof cat --n 64 --w 3 [--json]`
  - `python3 tools/faultweb/cli.py synth steane_5anc --out c.json --proof-out cert.json --stim-out c.stim --qasm-out c.qasm`
  - `python3 tools/faultweb/cli.py contract steane_5anc --p 1e-4 --idle-ratio 1/10 --noise-model circuit_css --code steane_7_1_3 [--json]`
  - `python3 tools/faultweb/cli.py audit tools/faultweb/examples/shor_opt_zzzz.json --w-max 3 [--json]`
- Python API:
  `from faultweb import CliffordCircuit, ZXDiagram, check_fault_equiv,
  circuit_distance, verify_derivation, synth_cat_state, synth_shor_syndrome,
  synth_steane_syndrome, synthesize_flag_decoder, build_gadget_contract,
  make_steane_7_1_3`.
- Exit codes: `0` = certified / pass, `1` = mathematical finding / falsification,
  `2` = input or side-condition refusal (`RefusalError`).

## GATES

1. **Zero-float AST & JSON gate (`L1`, `L5`)**: AST scan forbids every `float`
   literal and `float` builtin call across `tools/faultweb/**/*.py`; JSON
   contract parser rejects any float token via `parse_float` trap.
2. **Disjoint dual-lineage agreement (`L2-L4`)**: every `CliffordCircuit`
   equivalence and distance check cross-verifies `pauliweb` against
   `frame_oracle` (`spackle_crosschecked=True`) and raises on any disagreement.
3. **Fail-closed composite CNOT extraction (`N2`, `M3`)**: extracting a noisy
   CNOT under `circuit_depolarizing` without either split/merge termination
   (`Props. 6.9-6.10`) or all four cross-sector fault gadgets (`Prop. 6.8`)
   raises `UNGUARDED-CNOT-EXTRACTION`.
4. **Rewrite side-condition gates (`N4`)**: `fuse_1_fe` refuses 0-leg scalar
   collapses (`FUSE1-ISOLATED-SCALAR-REFUSED`) and multi-leg hub merges with
   $\min(\text{ext}_u, \text{ext}_v) > 1$ (`FUSE1-DEGREE-REFUSED`);
   `ideal_fuse_fe` refuses non-purple noisy edges (`NON-IDEAL-EDGE-IN-IDEAL-RULE`).

## FOOTGUNS

- **`circuit_css` vs `circuit_depolarizing` noise semantics**: Kissinger's
  Goto/bent-leg Steane constructions (`Props. 7.10-7.12`) satisfy strict
  3-fault equivalence under `circuit_css` (separate $X$ and $Z$ fault sectors),
  not joint-weight `circuit_depolarizing` (where cross-sector $X_i Z_j$ CNOT
  faults have joint support 2 even though each CSS sector has weight $\le 1$).
  Every receipt records BOTH `noise_model` and `bad_css_sector_labels` so the
  two regimes are never conflated.
- **Open boundary wires have circuit distance 1**: per `Definition 2.8` (`p. 9`),
  faults on open input/output wires are undetectable in isolation; evaluate
  `check_fault_equiv` against an idealised specification, or compose with an
  idealised verifier (`Theorem 3.11`) when testing `circuit_distance`.
- **Single-shot Steane requires `codespace`**: `Steane` extraction is NOT
  fault-equivalent on arbitrary unencoded inputs (`Appendix C.1`); pass
  `codespace=make_steane_7_1_3()` (`--code steane_7_1_3`).
- **Colour-dual Z-sector Steane**: never wrap noisy Hadamards around data wires
  to switch CSS sectors under `circuit_css` (doing so mixes $X$ and $Z$ fault
  sectors); `synth_steane_syndrome(..., sector="Z")` emits the exact colour-dual
  circuit (`Prop. 7.8`, zero data Hadamards).

## Battery (in-tree, synthetic, self-contained)

One command: `python3 tools/faultweb/selftest.py` (~3 s from a cold clone,
stdlib-only, executes inside a scratch temporary directory). Runs 10 legs:

- **L1 (`L1_purity_zero_float`)** — AST purity scan across all package modules.
- **L2 (`L2_cat_and_64q_proof`)** — `Props. 7.1-7.3` (4-, 5-, 8-qubit full &
  $w$-truncated cat states) + 64-qubit polynomial-time proof check.
- **L3 (`L3_shor_and_outcome_map`)** — `Props. 7.5-7.7` (`standard`, 3-ancilla
  `optimised`, rail-shifting `teleport`) on `ZZZZ`, `XXXX`, `XZYZ` + `Appendix A`
  outcome push-out map.
- **L4 (`L4_steane_and_flag_decoder`)** — `Props. 7.8-7.12` (`8-anc` & `5-anc`
  in both $X$ and $Z$ sectors), closed distance $== 3$, + `Appendix C.2`
  flag-conditioned decoder certification.
- **L5 (`L5_contracts_examples_cli`)** — outward-rounded `Fraction` containment,
  zero-float `faultweb-contract-v1` contracts, `examples/*.json`, OpenQASM 2.0
  round-trip, `.stim` `DETECTOR` emission, optional `pyzx`/`stim` probes, and
  CLI `0`/`1`/`2` exit codes.
- **N1-N5 (`N1_N5_negatives_and_refusals`)** — unflagged 4-cat failure, bare
  CNOT refusal, missing-encoder Steane failure, invalid ZX rewrites, and
  non-Clifford QASM refusal.
- **M1-M4 (`M1..M4_mutant_*`)** — scratch-copy single-line sabotages (grouped
  $Y=XZ$ fault weight, $\pm\frac{\pi}{2}$ spider parity coupling, CNOT extraction
  guard, and `ROUND_CEILING` $\to$ `ROUND_FLOOR` outward rounding), requiring all
  four mutants to be caught.

CREDIT: implements the fault-equivalence, fault-gadget, fault-equivalent ZX
rewrite, cat-state, Shor, and 5-ancilla Steane constructions of Rodatz, Poór,
and Kissinger [RPK25] (building on distance-preserving rewrites [RPK24], the ZX
calculus of Coecke and Duncan [CD08], Pauli webs of Bombin et al. [BLNP24] and
Borghans [Bor19], spacetime codes and spackle maps of Bacon et al. [BFHS17] and
Delfosse and Paetznick [DP23], Goto's Steane state preparation [Goto16], and
Gidney's Stim detector model [Gid21]); the disjoint-lineage verifier, automated
flag-decoder compiler, and zero-float outward-rounded stochastic contract layer
are ours (`pyzx` and `stim` optional bridges invoke upstream packages only when
installed). Bracketed keys resolve in `REFERENCES.md` at the repository root.

## License

MIT License (code) and CC BY 4.0 (documentation); see repository root `LICENSE`
and `LICENSE-CONTENT`.
