# ftqc_uncompute — measurement-based uncomputation, spooky pebbling, and superposition masking

`ftqc_uncompute` verifies logical measurement-based uncomputation, validates spooky pebbling and superposition-mask specifications, and estimates arithmetic rewrite costs. Its CLI and examples are in this package. [FTQC integration guide](../FTQC.md)

KIND: library + CLI (`python3 -m tools.ftqc_uncompute` / `python3 tools/ftqc_uncompute/selftest.py`), 100% Python standard library (`fractions.Fraction`, `decimal.Decimal`, `int` bit-packing; zero floating-point literals or `float()` calls across all certified modules and JSON receipts). Verification class: **selftest** (runs green in ~1 s from a cold clone with zero external packages). Optional bridges under `adapters/` connect to sibling `ftqc_audit` and `ftqc_zx` packages in-tree, pure-stdlib OpenQASM 3 dynamic circuits, and optional external `qualtran`, `cirq`, `qiskit`, and `kickmix` installations (which skip cleanly by name when absent).

## PURPOSE

Reversible quantum arithmetic routinely spends half its non-Clifford gates and workspace qubits uncomputing intermediate carries, products, table lookups, and modular powers. Three measurement- and state-assisted techniques circumvent naive Bennett compute-copy-uncompute overhead, each carrying subtle soundness traps that classical computational-basis testing misses:

1. **KickMix circuits & Hadamard-Measure-Reset (`HMR`)** (`hmr.py`, `kickmix_sim.py`, `kraus_oracle.py`):
   - Symbolic Algebraic Normal Form (ANF) + $(\mathbb{Q}/2\mathbb{Z})$ multilinear phase-polynomial solver (`verify_hmr_circuit`) proving exact cancellation of `HMR` phase kickbacks $(-1)^{r_t \cdot f(x)}$ against classically conditioned diagonal fixups (`Z`, `CZ`, `CCZ`, `NEG`) for `del_and`, 0-Toffoli `del_gf2_mul`, and LSB-folded `del_lookup`, while enforcing dirty-reset absence (`DIRTY-RESET-HAZARD`), outcome-independent basis outputs (`OUTCOME-DEPENDENT-BASIS-OUTPUT`), and branch injectivity (`NON-INJECTIVE-KICKMIX-BRANCH`) required for Khattar et al.'s coherent trace-distance theorem.
   - 64-lane bit-packed Fiat-Shamir trajectory engine (`simulate_trajectories`) seeded by `SHA256(canonical_bytes || salt)` with outward-rounded `Decimal` Clopper-Pearson bounds $\varepsilon^+$ ($\delta = 2^{-128}$), coherent trace-distance enclosures $D_{\mathrm{tr}} \le 2\sqrt{\varepsilon^+} + 2\pi\tau$, `Z_POW` $T$-synthesis budgets, and DQI repetition factors.
   - Disjoint-lineage small-width ($n \le 12$) dense Kraus channel oracle (`verify_channel_kraus`, `exact_ensemble_trace_distance`) evolving full quantum superpositions in 44-digit `Decimal` arithmetic with subspace Gram-Schmidt + algebraic Jacobi trace-norm diagonalization.
2. **Spooky & parallel spooky pebbling** (`pebbling.py`):
   - Exact validator (`validate_pebbling_schedule`) over `{EMPTY, PEBBLE, GHOST}` states on line chains and general DAGs (predecessor presence, parallel active-set disjointness, ancilla-weighted peak space, and `UNCLEANED-GHOST-AT-EXIT` refusal).
   - KMRV26 $A$-sequence `Blast`/`Unblast` compiler (`compile_kmrv_blast_schedule`) achieving exact optimal parallel depth $2\ell$ in $\le \lceil 2.47\log_2\ell - 2.77\rceil$ pebbles, Gidney's sequential divide-and-conquer spooky pebbler, and an admissible $A^*$ solver (`solve_optimal_pebbling_astar`) supporting heterogeneous Regev $w$-window ancilla weights.
3. **Superposition masking & streaming traveling-eigenstate ladders** (`masking.py`):
   - Validates exact ($\lvert M_N\rangle$, $\lvert +\rangle^{\otimes w}$, $\lvert\mathrm{grad}_w\rangle$) and approximate ($\lvert R_N\rangle$, carry-runway) masks, refusing non-coprime branch multipliers (`MASK-NON-BIJECTIVE-GARBAGE`) and average-case-truncated subroutines (`MASK-SUPPORT-OVERFLOW-HAZARD`, catching $1.413n$ Euclidean Dialog truncation vs worst-case $2n - 1$ steps on full-support masked eigenstates).
   - Symbolic $\mathrm{GL}_2(\mathbb{Z})$ prover (`verify_streaming_fibonacci_ladder`) certifying that `Fib`/`Fib^-1` conjugation decouples streamed IQPE controls from $\lvert M_N\rangle\otimes\lvert u_\theta\rangle$ with zero reverse uncomputation pass, paired with factor-free outward-rounded `Decimal` bounds $\sqrt{1 - \varphi(N)/N}$.
4. **Simplification advisor & verified rewriter** (`suggest.py`):
   - Scans `.kmx` circuits and Qualtran-style `ArithDAG`s for rules `R1_TEMP_AND_HMR`, `R2_BILINEAR_HMR`, `R3_LOOKUP_HMR`, `R4_SPOOKY_PEBBLE`, `R5_STREAM_MBUC`, `R6_SUPERPOS_MASK`, and `R7_ADDER_PHASE_GRAD`, prices physical FTQC wire-tick savings against `ftqc_zx-contract-v1` contracts, and reports which rules were applied. KickMix verification compares the public basis map and phase of both circuits. Verified DAG rewriting requires a concrete lowering and checks its HMR obligations.

A global proof that exceeds its symbolic budget produces `SYMBOLIC-BUDGET-EXCEEDED`. Affine maps use exact binary rank. Nonlinear maps with at most twelve input bits can use exhaustive enumeration. A reversible circuit with restored scratch can establish injectivity directly. An unresolved wider nonlinear map produces `INJECTIVITY-NOT-PROVED`. A refused result cannot establish coherent cleanup.

`verify` resolves the source or accepts an explicit `--target`, rechecks the digest and compares every regenerated field, including the trajectory and Kraus options. Trajectory bounds refer to the selected oracle and sampling assumptions. A trace-preservation probe has a separate mode in the receipt.

## USE-WHEN

- Validating any `.kmx`, OpenQASM 3 dynamic circuit, or `ArithDAG` that uses mid-circuit `HMR` / `R` resets, conditioned phase fixups, spooky pebbling, or masked accumulators before trusting computational-basis tests.
- Simplifying reversible arithmetic circuits (adders, $\mathrm{GF}(2^m)$ multipliers, QROM tables, modular exponentiation chains, DQI decoders) to reduce Toffoli/$T$ count, workspace qubits, or physical QEC spacetime volume.
- Computing certified outward-rounded `Decimal` trace-distance and total-variation budgets (`ftqc_uncompute-receipt-1`) for approximate masks ($\lvert R_N\rangle$, carry runways) or randomized KickMix trajectory audits.

## NOT-FOR

- Arbitrary non-monomial quantum circuits inside KickMix branches (general `H` or $\sqrt{X}$ gates that create superpositions mid-branch rather than `HMR` measurement-resets): `KickMixIR` purposefully restricts gates to computational-basis permutations (`X`, `CX`, `CCX`, `SWAP`, `CSWAP`), diagonal phases (`Z`, `CZ`, `CCZ`, `Z_POW`, `NEG`), classical logic, `R`, and `HMR` so every branch $V_r$ is monomial and efficiently verifiable.
- Physical fault-distance or syndrome-decoder proofs below the logical gate level: lower Clifford+`HMR` gadgets via `adapters.ftqc_bridge` into [`ftqc_audit`](../ftqc_audit/GUIDE.md) or [`ftqc_zx`](../ftqc_zx/GUIDE.md) for location-level fault-tolerance proofs.

## INVOKE

Run from the BootLoops repository root:

```sh
# Validate a .kmx circuit or JSON spec and emit a zero-float receipt
python3 -m tools.ftqc_uncompute check tools/ftqc_uncompute/examples/dqi_rs_tiny_mbuc.kmx -o /tmp/dqi.receipt.json
python3 -m tools.ftqc_uncompute verify /tmp/dqi.receipt.json --target tools/ftqc_uncompute/examples/dqi_rs_tiny_mbuc.kmx

# Suggest and apply MBUC simplifications (with automatic post-rewrite verification)
python3 -m tools.ftqc_uncompute suggest my_circuit.kmx --apply -o my_circuit_opt.kmx

# Synthesize KMRV26 depth-2*ell or space-capped A* parallel spooky pebbling schedules
python3 -m tools.ftqc_uncompute pebble --length 20 --mode kmrv
python3 -m tools.ftqc_uncompute pebble --length 12 --mode astar --pebbles 5 --window 2

# Audit superposition mask bijectivity, Dialog worst-case safety, and Decimal d_TV bounds
python3 -m tools.ftqc_uncompute mask --kind uniform_mod_n --modulus 15 --multipliers 2,4 --runs 4

# Convert between .kmx, OpenQASM 3 (.qasm), and ArithDAG (.json)
python3 -m tools.ftqc_uncompute convert tools/ftqc_uncompute/examples/gidney_adder_4bit.kmx -o /tmp/adder.qasm

# Run the cold-clone acceptance battery
python3 tools/ftqc_uncompute/selftest.py
```

## WIDE ARITHMETIC: CUTS AND BLOCK PROGRAMS

`verify_hmr_circuit` writes every wire as one polynomial in the circuit inputs. A ripple carry doubles at every bit, so a 16-bit adder already exceeds the budget and is refused with `SYMBOLIC-BUDGET-EXCEEDED`. `cuts.py` proves the same obligations locally.

A **cut** replaces the value of a wire by a fresh variable and keeps its polynomial as the definition of that variable. An identity that holds for every value of the cut variables holds for their actual values, so a proof is sound wherever the cuts are placed; cuts only decide whether a proof is found. `verify_with_cuts` places them itself (a wire of degree two or more that is about to control a non-linear gate, or a wire whose polynomial has grown long) and accepts `INFO_TAG "@cut:q7"`. A wire that depends on a measurement outcome is never cut. When an obligation is not an identity in the free variables, the latest cut variable in it is replaced by its definition, and so on until none is left; a polynomial that is still not zero then comes with an input at which it is not zero. The checks are the ones of `hmr.py`: dirty resets, scratch at exit, outputs independent of outcomes, cancellation of every outcome-dependent phase. Branch injectivity follows from the last: a kickback (-1)^(r q) can only be cancelled by diagonal fix-ups on the remaining wires, so q is a function of them.

`prove_relation` shows that an integer combination of products of wire values, taken at entry, at exit or at a probe (`INFO_TAG "@probe:name"`), is zero for every input. Values are lifted to exact integer multilinear polynomials and the cut variables are eliminated latest first. On a carry chain each elimination is a + b + k = (a xor b xor k) + 2k', so the polynomial stays small. A multilinear polynomial in the inputs that vanishes on all 0/1 points has all coefficients zero, so the final test is exact.

`blocks.py` applies both to **block programs**: lists of logical operations on named wires, `ZERO`, `CX`, `XF` (logical X), `AND` (onto a fresh wire, by consuming an AND resource), `ERASE t u v z` (measurement deletion, claimed valid because t = (u AND v) XOR z), `DROP` (retire a wire claimed to hold 0), with optional `CUT` hints and `PROBE` marks.

```sh
python3 -m tools.ftqc_uncompute blocks tools/ftqc_uncompute/examples/block_adder_4bit.json -o /tmp/add4.receipt.json
python3 -m tools.ftqc_uncompute blocks tools/ftqc_uncompute/examples/block_adder_4bit.json --verify /tmp/add4.receipt.json
```

A specification is `{"kind": "block-program", "program": ..., "spec": {"relation": [[coef, [factor, ...]], ...], "preserved": [wires], "reading": text}}`, with a factor `{"at": "in" | "out" | probe, "xor": [wires]}`. The result is `proved` when every deletion is exact, no phase is left, the relation is an identity and the preserved wires leave as they entered. The relation is a statement; its `reading` (for instance "B_out + 2^n C = B + A") is carried into the receipt unchecked. `blocks.adder`, `subtractor`, `comparator` and `multiplier` build standard programs and `blocks.arithmetic_spec` their relations. A comparator is stated as A + (2^n − 1 − B) = S + 2^n F with S a sum of n bits read at a probe, from which F = [A > B].

With `"layers"` (a scheduled program in the format of `ftqc_gadgets cones`) the result also **binds** the schedule to the program: without its error-correction steps it must have the same operations on every wire and every classical bit in the same order as the tool's own lowering, every controlled correction in a later layer than its measurement, and no wire used twice in a layer. Two words with equal projections on every wire and bit are the same circuit. The lowering replaces `AND` and `ERASE` by measurement gadgets, and `gadget_identities` runs exactly those lowered operations on every basis input and outcome. The identities hold only with their Pauli-frame terms (X on t when both outcomes of an AND are 1; Z on z when the outcome of an ERASE is 1). Those terms and `XF` have no location in the layers. They belong to the classical control, which is assumed perfect, and the AND resource is assumed ideal.

## KNOWN HAZARDS

1. **Non-injective KickMix branches invalidate $D_{\mathrm{tr}} \le 2\sqrt{\varepsilon} + 2\pi\tau$**: Khattar et al.'s coherent trace-distance bound requires each fixed-outcome branch operator $V_r$ to preserve state norm on superpositions (injective computational-basis action). Measuring an uncopied input qubit (`HMR q0 b0`) or conditioning `X` on `b0` without restoring injectivity passes naive single-basis checks on $|0\rangle$ while collapsing $|+\rangle \mapsto |0\rangle$ ($D_{\mathrm{tr}} > 0.707$). `hmr.py` catches this via `NON-INJECTIVE-KICKMIX-BRANCH`, and `kickmix_sim.py` sets `coherent_bound_valid = False` whenever branch injectivity or clean-reset invariants fail.
2. **Dirty `R` resets silently dephase superpositions**: Applying `R q` when $q$ can be $|1\rangle$ traces out $q$ into a mixed state; `hmr.py` refuses with `DIRTY-RESET-HAZARD` and minimal witness $(x^*, r^*)$, while `kickmix_sim.py` injects KickMix's random $\pm 1$ trajectory sign flip so dirty resets fail phase checks with probability $1/2$.
3. **Average-case-truncated arithmetic fails under superposition masking**: Masked registers ($\lvert M_N\rangle$, $\mathrm{Fib}^{-k}(\lvert M_N\rangle\lvert u_\theta\rangle)$) populate correlated full-support superpositions over $\mathbb{Z}_N^\times \times \mathbb{Z}_N$, hitting worst-case $2n - 1$ Euclidean chains in Kaliski/Dialog multipliers where $1.413n$ average-case truncations overflow (`MASK-SUPPORT-OVERFLOW-HAZARD`).
4. **Stranded spooky ghosts at exit**: Leaving a node in `GHOST` state frees its qubits (`0` ancilla leak!) yet leaves a random branch-dependent $(-1)^{\langle b_v, x_v\rangle}$ phase error on the output superposition (`UNCLEANED-GHOST-AT-EXIT`).

## ACCEPTANCE BATTERY

- `python3 tools/ftqc_uncompute/selftest.py` (registered in `tools/BATTERIES.json` under `ftqc_uncompute`, `cwd: root`) runs the following checks together with regressions for proof budgets, injectivity, evidence replay, feedback, rewrites, mask inputs and numerical bounds.
  - `L0_purity_and_lineage`: AST zero-float/zero-complex scan across all package modules and disjoint-lineage import check on `kraus_oracle.py`.
  - `L1_ir_gf2_qasm3`: $\mathrm{GF}(2^m)$ field/LUP/bilinear identities, `.kmx` canonical SHA-256 invariance, and OpenQASM 3 dynamic-circuit round-trip.
  - `L2_hmr_sim_kraus_dual`: Dual-lineage agreement (`hmr.py` + `kickmix_sim.py` + 44-digit `kraus_oracle.py` with $D_{\mathrm{tr}} = 0$) on `del_and`, 4-bit Gidney `iadd`, 0-Toffoli `del_gf2_mul`, and QROM `del_lookup`.
  - `L3_spooky_pebbling`: KMRV26 `Blast`/`Unblast` depth-$2\ell$ and space bounds ($\ell = 1..30$), Gidney sequential pebbler, and $A^*$ windowed solver.
  - `L4_superposition_masking`: $\lvert R_{15}\rangle$ vs $\lvert M_{15}\rangle$ and carry-runway bounds cross-checked to $20+$ digits against `kraus_oracle.py`, RSA-2048 factor-free bound, and streaming Fibonacci ladder proof.
  - `L5_advisor_r1_r7`: All seven simplification rules (`R1`–`R7`) and `--apply` equivalence verification.
  - `L6_receipts_and_bridges`: `ftqc_uncompute-receipt-1` round-trips on all 5 `examples/*` fixtures plus live `ftqc_audit` and `ftqc_zx` bridges.
  - `L7_cuts_and_block_programs`: the cut verifier against `hmr.py` wire by wire and against a run on every input; relations against integer arithmetic; every single-operation corruption of small adders, comparators and multipliers refused or shown to be the same function; hazard witnesses replayed; bindings with corrupted layers refused; 64-bit adder proved where the global method exceeds its budget. Thirty-nine planted defects in `cuts.py` and `blocks.py` each fail this leg; four needed a new check first.
  - `N1_N6_sabotage_controls`: Six planted-failure controls (`UNCOMPENSATED-HMR-PHASE`, `DIRTY-RESET-HAZARD`, `NON-INJECTIVE-KICKMIX-BRANCH`, `UNCLEANED-GHOST-AT-EXIT`, `MASK-NON-BIJECTIVE-GARBAGE` / `MASK-SUPPORT-OVERFLOW-HAZARD`, `FLOAT-LITERAL-IN-RECEIPT`).
  - `M1_M4_mutant_kills`: Four subprocess code-mutation self-audits on isolated temp copies confirming zero surviving mutants.
- Expected summary line: `SELF-TEST: ALL PASS (rc 0)`.

CREDIT: Clean-room MIT implementation authored from published specifications (zero third-party source lines or competition files copied). KickMix gate model, `.kmx` grammar, Fiat-Shamir trajectory verification, coherent trace-distance theorem, and 0-Toffoli $\mathrm{GF}(2^m)$ / DQI MBUC constructions follow T. Khattar, C. Gidney, M. P. Harrigan, R. Babbush and co-authors [Khattar26] (whose upstream `kickmix` package is Apache-2.0 and is never bundled here). Temporary logical-AND (`init_and`/`del_and`), ripple-carry HMR addition, and phase-gradient catalysis follow C. Gidney [Gid18]; QROM `del_lookup` follows D. W. Berry, C. Gidney, M. Motta, J. R. McClean and R. Babbush [BCB19] and M. B. Mansky et al.; sequential spooky pebbling and carry runways follow C. Gidney [Gid19a] and N. Kornerup, J. Sadun and D. Soloveichik [KSS25]; parallel spooky pebbling (`Blast`/`Unblast`, $A$-sequence, and canonical backward-$A^*$) follows G. D. Kahanamoku-Meyer, S. Ragavan and K. Van Kirk [KMRV26]; superposition masking ($\lvert M_N\rangle$, $\lvert R_N\rangle$), streaming Fibonacci traveling-eigenstate ladders, and the worst-case Euclidean Dialog support bound follow K. Van Kirk and C. Gidney [VKG26]. Bracketed keys resolve in `REFERENCES.md` at the repository root.
