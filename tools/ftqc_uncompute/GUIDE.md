# ftqc_uncompute — GUIDE

Overview of the FTQC tools: [`../FTQC.md`](../FTQC.md).

KIND: library and CLI (`python3 -m tools.ftqc_uncompute`, or
`python3 -m ftqc_uncompute` with `tools/` on `PYTHONPATH`), Python standard
library, no floating-point literals. Circuit forms: `ir.py`, `ir_dag.py`,
`ir_kickmix.py` (`.kmx`), `gf2.py`. Checks: `hmr.py` (symbolic HMR proofs),
`cuts.py` and `blocks.py` (wide arithmetic), `kickmix_sim.py` (sampled
trajectories), `kraus_oracle.py` (dense numerical check for small widths),
`pebbling.py`, `masking.py`, `decimal_bounds.py`. `suggest.py` (rewrite
advisor), `certify.py`, `cli.py`, `adapters/` (bridges to `ftqc_faults`,
OpenQASM 3, and optional `qualtran`, `cirq`, `qiskit`, `kickmix`, which skip by
name when absent), `selftest.py`, `selftest_blocks.py`, `examples/`. Battery
class: selftest.

## PURPOSE

Logical-level checks of three ways to avoid the compute-copy-uncompute cost of
reversible arithmetic, each with a trap that testing on computational basis
states does not catch.

1. **Measurement-based uncomputation.** A scratch value is erased by a
   Hadamard-measure-reset (HMR); the measurement leaves a phase
   (−1)^(r·f(x)) that classically controlled diagonal fix-ups must cancel.
   `verify_hmr_circuit` writes every wire as a polynomial in the inputs and
   proves, for every input and every outcome, that scratch is restored, resets
   are clean, public outputs do not depend on outcomes, every
   outcome-dependent phase cancels, and each branch is injective (the
   condition of Khattar et al.'s coherent trace-distance bound). `cuts.py` and
   `blocks.py` prove the same obligations locally for arithmetic too wide for
   one polynomial, state the arithmetic as an exact integer relation, and bind
   the program to a correction schedule.
2. **Spooky pebbling.** `validate_pebbling_schedule` checks a schedule over
   the states empty, pebble and ghost on chains and DAGs (predecessors present,
   parallel steps disjoint, weighted peak space, no ghost left at exit). The
   package also builds schedules: the A-sequence Blast/Unblast schedule of
   Kahanamoku-Meyer, Ragavan and Van Kirk (depth 2ℓ), Gidney's sequential
   schedule, and an A* search with a pebble cap.
3. **Superposition masking.** Checks that a mask map is a bijection on the
   masked register, refuses subroutines truncated for the average case when
   the masked state has full support, proves the streaming Fibonacci ladder
   identity over GL₂(Z), and gives outward-rounded total-variation bounds such
   as √(1 − φ(N)/N).

`suggest` recognizes seven rewrite rules (R1–R7) in `.kmx` circuits and
arithmetic DAGs; with `--apply` it applies the KickMix rewrites and checks the
public action of the result.

## SCOPE

Logical circuits whose branches are monomial: permutations of basis states
(X, CX, CCX, SWAP, CSWAP), diagonal phases (Z, CZ, CCZ, Z_POW, NEG), classical
logic, resets and HMR. The global symbolic check has a budget; past it,
`check` refuses (`SYMBOLIC-BUDGET-EXCEEDED`). Branch injectivity is proved by
exact rank for affine maps, by enumeration up to twelve input bits, or from
restored scratch; otherwise `check` refuses (`INJECTIVITY-NOT-PROVED`).

## USE-WHEN

- A `.kmx`, OpenQASM 3 or arithmetic-DAG circuit uses mid-circuit HMR or
  resets with conditioned phase fix-ups, spooky pebbling or masked
  accumulators, and must be checked before basis-state tests are trusted.
- Wide arithmetic (adders, comparators, multipliers) needs a proof with a
  correction schedule bound to it.
- A rewrite to save Toffolis or qubits should be applied and checked.

## NOT-FOR

- Branches that create superpositions (general H or √X inside a branch).
- Physical fault distance: lower the Clifford and HMR parts through
  `adapters/ftqc_bridge.py` to `ftqc_faults`.
- End-to-end resource estimates: the pricing functions combine logical counts
  with supplied factory, contract and latency numbers.

## INVOKE

From the BootLoops root:

```sh
python3 -m tools.ftqc_uncompute check tools/ftqc_uncompute/examples/dqi_rs_tiny_mbuc.kmx -o /tmp/dqi.receipt.json
python3 -m tools.ftqc_uncompute verify /tmp/dqi.receipt.json --target tools/ftqc_uncompute/examples/dqi_rs_tiny_mbuc.kmx
python3 -m tools.ftqc_uncompute suggest my_circuit.kmx --apply -o my_circuit_opt.kmx
python3 -m tools.ftqc_uncompute pebble --length 20 --mode kmrv
python3 -m tools.ftqc_uncompute pebble --length 12 --mode astar --pebbles 5 --window 2
python3 -m tools.ftqc_uncompute mask --kind uniform_mod_n --modulus 15 --multipliers 2,4 --runs 4
python3 -m tools.ftqc_uncompute convert tools/ftqc_uncompute/examples/gidney_adder_4bit.kmx -o /tmp/adder.qasm
python3 -m tools.ftqc_uncompute blocks tools/ftqc_uncompute/examples/block_adder_4bit.json -o /tmp/add4.receipt.json
python3 tools/ftqc_uncompute/selftest.py
```

`verify` resolves the source or takes `--target`, checks the digest and
compares every regenerated field.

## OUTPUTS

`ftqc_uncompute-receipt-1` receipts. The symbolic proof, the trajectory
sampling (a confidence bound, not a proof), and the optional dense Kraus
check (numerical, 44-digit `Decimal`, widths up to 12) are separate fields.
Refusals carry a named code and, where one exists, a witness input and
outcome.

## GATES

`python3 tools/ftqc_uncompute/selftest.py`: 22 legs in 2.8 s on a four-core
Linux virtual machine; last line `SELF-TEST: ALL PASS (rc 0)`.

- `L0` no floats or complex numbers; separate lineage of `kraus_oracle.py`.
- `L1` GF(2^m) identities, `.kmx` digest invariance, OpenQASM 3 round trip.
- `L2` the symbolic check, sampling and the Kraus oracle agree (trace distance
  0) on `del_and`, a 4-bit Gidney adder, `del_gf2_mul` and `del_lookup`.
- `L3` Blast/Unblast depth and space for ℓ = 1..30, the sequential pebbler,
  A*.
- `L4` masks against the Kraus oracle to 20 digits, the RSA-2048 bound, the
  Fibonacci ladder proof.
- `L5` all seven rules and `--apply` equivalence.
- `L6` receipt round trips on the examples and the live `ftqc_faults`
  bridges.
- `L7` the cut verifier against `hmr.py` and against evaluation on every
  input; relations against integer arithmetic; every single-operation
  corruption of small adders, comparators and multipliers refused or shown to
  compute the same function; corrupted schedule bindings refused; a 64-bit
  adder proved where the global method exceeds its budget. 39 planted defects
  in `cuts.py` and `blocks.py` each fail this leg.
- `N1–N6` planted failures: an uncompensated HMR phase, a dirty reset, a
  non-injective branch, a ghost left at exit, a non-bijective or overflowing
  mask, a float in a receipt.
- `M1–M4` four code mutations in scratch copies are caught.

## FOOTGUNS

1. **Non-injective branches void the coherent bound.** Khattar et al.'s
   bound D_tr ≤ 2√ε + 2πτ needs each fixed-outcome branch to be injective on
   basis states. Measuring an uncopied input passes a test on |0⟩ and
   collapses |+⟩ to |0⟩. `hmr.py` refuses with `NON-INJECTIVE-KICKMIX-BRANCH`
   and the sampler marks `coherent_bound_valid = False`.
2. **Dirty resets dephase.** `R q` on a qubit that can be |1⟩ traces it out;
   `DIRTY-RESET-HAZARD` comes with a witness input and outcome.
3. **Average-case truncation fails under masking.** A masked register has full
   support, so a Euclid-type loop truncated at its average length (about
   1.413n steps) overflows on the worst case of 2n − 1
   (`MASK-SUPPORT-OVERFLOW-HAZARD`).
4. **Ghosts left at exit** free their qubits but leave a branch-dependent
   phase on the output (`UNCLEANED-GHOST-AT-EXIT`).
5. **Sampling is not proof.** The trajectory bound is a Clopper–Pearson
   confidence bound at δ = 2⁻¹²⁸ for the selected oracle; only the symbolic and
   block proofs are proofs.
6. In block programs the Pauli-frame terms of the AND and ERASE gadgets, the
   logical X frame, the AND resource and the classical control are ideal; the
   relation's `reading` text is carried into the receipt unchecked.

## KNOWN GAPS

- Everything is at the logical level; physical faults enter only through the
  bridge to `ftqc_faults` for Clifford and HMR parts.
- Resource prices use supplied numbers (factory costs, contract volumes,
  decoder latency, staged service budgets).
- DAG rewriting is verified only where a concrete lowering exists; abstract
  nodes are costed and not proved.

## Wide arithmetic: cuts and block programs

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

## Code automorphisms

Pass `installed_action_per_reaction=spec` to `estimate_reaction_decoder_latency` to charge one verified physical action per reaction. Supply the allocated data width, peak ancillas and coupler budget in the usual schedule metrics. The estimator replays the route and charges its elapsed ticks once per reaction, alongside active rounds and decoder stalls. See the [installation example](../FTQC.md#code-automorphisms).

`python3 -m ftqc_uncompute price-staged STUDY --feedback-ticks F --terminal-wire-ticks T --seed-wire-ticks S` prices the reserved online wire-ticks of a staged gate program and charges the supplied feedback delay and service budgets; see the [staged logical gate guide](../ftqc_gadgets/STAGED.md).

CREDIT: Implemented from the published specifications; no third-party source is bundled. KickMix gate model, `.kmx` grammar, Fiat-Shamir trajectory verification, coherent trace-distance theorem, and 0-Toffoli $\mathrm{GF}(2^m)$ / DQI MBUC constructions follow T. Khattar, C. Gidney, M. P. Harrigan, R. Babbush and co-authors [Khattar26] (whose upstream `kickmix` package is Apache-2.0 and is never bundled here). Temporary logical-AND (`init_and`/`del_and`), ripple-carry HMR addition, and phase-gradient catalysis follow C. Gidney [Gid18]; QROM `del_lookup` follows D. W. Berry, C. Gidney, M. Motta, J. R. McClean and R. Babbush [BCB19] and M. B. Mansky et al.; sequential spooky pebbling and carry runways follow C. Gidney [Gid19a] and N. Kornerup, J. Sadun and D. Soloveichik [KSS25]; parallel spooky pebbling (`Blast`/`Unblast`, $A$-sequence, and canonical backward-$A^*$) follows G. D. Kahanamoku-Meyer, S. Ragavan and K. Van Kirk [KMRV26]; superposition masking ($\lvert M_N\rangle$, $\lvert R_N\rangle$), streaming Fibonacci traveling-eigenstate ladders, and the worst-case Euclidean Dialog support bound follow K. Van Kirk and C. Gidney [VKG26]. Bracketed keys resolve in `REFERENCES.md` at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz;
code written by Claude (Anthropic) under his supervision. Documentation is
licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository
root.
