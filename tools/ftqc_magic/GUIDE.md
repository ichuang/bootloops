# ftqc_magic — GUIDE

Overview of the FTQC tools: [`../FTQC.md`](../FTQC.md).

KIND: package (Python standard library for the default tier; no
floating-point literals; optional `cadical` for live LRAT proofs and `stim`
for cross-checks). Exact arithmetic: `rings.py` (Z[√2], Q(√2), Z[ω₈],
polynomials over Q), `ir.py` (phase-polynomial targets). Factories:
`matrix.py` (binary and PPR), `galois.py` (GF(2^s) trace codes), `limits.py`
(non-triorthogonal codes, lower bounds), `recycle.py` and `search.py`
(check-qubit recycling), `catalog.py`. Cultivation and escape: `cultivate.py`,
`caliper.py`. qLDPC factories: `qldpc_factory.py`. Composition:
`pipeline.py`, `delivery.py`, `action_program.py`, `insitu.py`, `staged.py`,
`bridges/`. `certify.py`, `cli.py`, `selftest.py`, `examples/`. Battery
class: selftest.

## PURPOSE

Exact statements about magic-state factories and the pipelines that combine
them, as `ftqc_magic-receipt-1` receipts, and the correlation rules that decide
whether their outputs may be used as independent.

- **Distillation.** For binary, PPR, GF(2^s)-trace and non-triorthogonal
  (signed stabilizer over Z[ω₈]) factories of T, CS, CCZ, CCZ^m and CCCZ: the
  target action, the distance, the exact acceptance and output error as
  polynomials in the input error rate (by dynamic programming over syndromes
  and cosets), the leading coefficient, and lower bounds on the number of
  inputs where known. Example: Bravyi–Haah 49-to-1, distance 5, 1411p⁵.
- **Check-qubit recycling.** Schedules that measure and reset check qubits in
  the middle of the circuit, with a replay of the physical qubit pool; for
  example 49-to-1 on 5 active qubits instead of 14. Optional LRAT proofs that
  fewer qubits are impossible for small cases.
- **Cultivation.** An exact simulation of the Clifford errors
  H± = (X ± Y)/√2 that arise when Pauli faults cross T checks, in time O(n³),
  which gives the |T⟩ fault distance d_T of a check schedule and the counts of
  malignant fault sets, beside the |S⟩-proxy distance d_S that a stabilizer
  simulation would report.
- **Escape stages.** Worst-case complementary gaps across open future
  boundaries of a matching graph (Caliper), cross-checked with the Blossom
  decoder of `ftqc_decoders`, and the expected volume of adaptive escape.
- **qLDPC CCZ factories.** Chain-complex checks of 3D tricycle and 4D
  quadcycle constructions, metacheck repair through `ftqc_decoders`, and the
  fan-out condition that keeps a transversal CCZ from halving the circuit
  distance.
- **Pipelines.** Stages composed with outward rounding, catalysed C2T
  conversion checked exactly over Z[ω₈], and contracts with a correlation mode
  (`INDEPENDENT_SINGLE`, `BATCH_CORRELATED`, `CATALYST_POISONED`) exported to
  `ftqc_compose`, `ftqc_extractors` and `ftqc_uncompute`.

## SCOPE

The distillation model is the standard one: independent Z-type errors of rate
p on the input T rotations, and ideal Clifford checks with postselection on
all +1 outcomes. Cultivation is analysed on its final check subcircuit with an
ideal encoded [[7,1,3]] |T⟩ input and depolarizing circuit noise on every
gate, preparation, idle and measurement. A pipeline feeds each stage's
conditional output error into the next as its input rate. Each receipt lists
its premises.

## USE-WHEN

- A factory's distance, leading error, acceptance or active-qubit count is an
  input to a resource argument and should be checked exactly.
- A cultivation schedule should be checked for Clifford-hook degradation
  (d_T < d_S) that a stabilizer simulation of |S⟩ misses.
- Factory outputs are consumed downstream and the consumer must know whether
  they are independent, batch-correlated or tied to a catalyst.

## NOT-FOR

- Circuit-level noise inside distillation Cliffords: the model treats them as
  ideal.
- Threshold-scale sampling of memories or escape stages: use
  `ftqc_decoders` with Stim or sinter and supply the resulting rate.
- Logical arithmetic cleanup: use `ftqc_uncompute`, fed by the pricing record
  this package exports.

## INVOKE

From the BootLoops root with `export PYTHONPATH="$PWD/tools"`:

    python3 -m ftqc_magic check SPEC.json -o receipt.json
    python3 -m ftqc_magic verify SPEC.json receipt.json
    python3 -m ftqc_magic recycle --factory bh_49_1_5_t [--lrat]
    python3 -m ftqc_magic cultivate --variant msc3_z_flagged --p 1/1000 [--external]
    python3 -m ftqc_magic catalog
    python3 -m ftqc_magic suggest --target CCZ --p 1/1000 --infidelity 1/1000000000
    python3 -m ftqc_magic switched-source STUDY.json

`pipeline` is an alias for `check` on a multi-stage spec. Exit codes: 0 ok;
2 refusal, which covers malformed input, a receipt that does not replay, and a
named hazard such as `UNFLAGGED-CULTIVATION-HOOK-DEGRADATION` (the refusal is
written to stderr as JSON).

## OUTPUTS

`ftqc_magic-receipt-1`: `schema`, `kind`, `spec_sha256`, `source_sha256` (a
hash of the package's Python sources), `contract` (target, distance,
active and total qubits, correlation mode, exact and decimal
`acceptance_lower`, `conditional_infidelity_upper`,
`unnormalized_failure_upper`, `expected_wire_ticks_per_output_upper`,
`premises`), `optimality_gap` where a lower bound is known, `exports` for the
sibling tools, and stage `details`.

## GATES

`python3 tools/ftqc_magic/selftest.py` followed by
`python3 -m unittest ftqc_magic.test_action_program`: 10 legs in 8 s on a
four-core Linux virtual machine without CaDiCaL or Stim.

- `L0` no floats or complex numbers; replay of the banked example receipts.
- `L1` exact rings, outward rounding, Clifford-orbit matching of targets.
- `L2` landmark binary and PPR factories, borrowed-identity transforms,
  parity lift, puncturing, GF(2^s) trace binarization, 3-twirl.
- `L3` Z[ω₈] collapse on the [[9,1,2]] 9-to-1 factory (acceptance 1/128,
  error 27/2 p²) and others, lower bounds on n.
- `L4` recycling of 15-to-1 (5 → 4) and 49-to-1 (14 → 5); banked LRAT no-go
  proofs re-checked with `ftqc_codes.lrat`.
- `L5` cultivation: unflagged MSC-3 has d_T = 2 with four malignant weight-2
  sets against d_S = 3; Z-flagged and two-round variants have d_T = 3.
- `L6` Caliper gaps and escape volume; tricycle and quadcycle checks.
- `L7` C2T, pipelines, `suggest`, bridges, CLI round trip.
- `N1–N6` the hazards below are refused; `M1–M4` four code mutations are
  caught.

## FOOTGUNS

1. **Clifford-hook degradation** (`UNFLAGGED-CULTIVATION-HOOK-DEGRADATION`).
   A Pauli fault halfway through a transversal T check becomes H±, which can
   spread into an undetected weight-2 logical error on |T⟩ even when the |S⟩
   proxy has distance 3. `cultivate` refuses a schedule whose d_T is below
   `claimed_t_distance`.
2. **Catalyst poisoning** (`CATALYST-POISONING-HAZARD`). C2T and phase
   catalysis reuse a catalyst; one error on it flips every later output.
   `CATALYST_POISONED` contracts are terminal.
3. **Batch correlation** (`MULTI-OUTPUT-CORRELATION-HAZARD`). Outputs of one
   accepted block of a multi-output factory share Z errors. Downstream use as
   independent needs `interleaving="distinct_instances"` or `"joint_3twirl"`.
4. **Depth-2 transversal CCZ** (`DEPTH2-TRANSVERSAL-CCZ-HALVING`) on tricycle
   codes without the d = 2 fan-out layer halves the circuit distance.
5. **Open boundaries** (`OPEN-BOUNDARY-NAIVE-GAP-HAZARD`). A terminal gap
   evaluated mid-circuit ignores defects on open future boundaries.
6. Receipts bind the package source hash: any code change invalidates old
   receipts, which must be reissued with `check`.

## KNOWN GAPS

- **Ideal Cliffords in distillation.** Factory error polynomials count only
  input errors. Faults in the Clifford checks are not modelled.
- **Heuristic surface-code error in pipelines.** A stage with `litinski`
  patch sizes adds 0.1·(100p)^((d+1)/2) per rotation to its input rate
  (`pipeline.LitinskiBlockSpec.topological_error_per_ppr`), a scaling formula
  that is not derived or certified here.
- **Pipeline distance is a product** of stage distances, a bookkeeping value
  rather than a proved circuit distance.
- **Supplied probabilities.** Caliper continuation and terminal acceptance
  probabilities are inputs; they enter the expected cost.
- **Cultivation acceptance.** The acceptance lower bound is 1 − Lp for L
  locations (Bernoulli). When Lp ≥ 1 the code substitutes a floor of 10⁻⁹,
  which is not a proven lower bound (the true acceptance can be smaller); do
  not use cultivation receipts in that regime.
- **Cultivation is a subcircuit.** The injection and growth stages before the
  final checks are not simulated; the encoded input is ideal.
- `suggest` searches the built-in catalog under these models; its "optimal"
  pipelines are optimal only within them.

## Delivery to encoded consumers

Use `ftqc_magic.delivery.read_delivery(spec, receipt)` for a strict replay
with normalized phases, ordered state groups, Clifford corrections,
acceptance, correlations and rational resource bounds. For binary producers
the two mask tables give the joint mass with acceptance and its conditional
law. Other supported producers can carry a state bound alone. Cost per
batch, per grouped quantum state and per native output are recorded
separately. For example, a binary CCZ can have three native output wires
but only one deliverable three-qubit state. Strict v1 receipts bind the
producer Python sources. To use a spec checked under another source revision,
reissue its receipt with `python3 -m ftqc_magic check`.

The [in-situ BB144 factory guide](../ftqc_extractors/BICYCLE_FACTORY.md)
covers injection batch replay, RM15 input assignment, signed native pivot
measurements, physical scheduling and full-attempt pricing.

## Code automorphisms

Use `check` and `verify` with a `rm15-action` specification containing a code, a signed action receipt and five ordered logical destination `slots`. Replay transports the prepare-plus axes, all 15 signed Pauli rotations, the four postselection axes and the output X and Z into that frame. The exponent applies to the signed axis. Use the installed Gross factory (`bb144-rm15` in `ftqc_compose`) for physical timing and source quality. See [code automorphisms](../FTQC.md#code-automorphisms).

`switched-source STUDY` applies the staged switched T of `ftqc_compose` to a fresh ideal encoded |+⟩ and reports a unit-yield source with the conditional infidelity bound of that model; see the [staged logical gate guide](../ftqc_gadgets/STAGED.md).

CREDIT: Implemented from the published equations and circuit specifications; no
third-party source is bundled. Punctured Reed–Muller and
triorthogonal matrix distillation follow S. Bravyi and A. Kitaev [BK05], S. Bravyi and
J. Haah [BH12], E. T. Campbell and M. Howard [CH17], J. Haah, M. B. Hastings, D. Poulin
and D. Wecker [HHPW17], and J. Haah and M. B. Hastings [HH18]; asymptotic and finite-length
$\mathrm{GF}(2^s)$ trace, Reed–Solomon, Hermitian AG, and 3-twirl concatenation protocols
follow A. Krishna and J.-P. Tillich [KT19], A. Wills, T.-C. Lin and M.-H. Hsieh [WZ25],
and A. Gong, C. A. Pattison, P. Rall and A. Wills [GPRW26]; borrowed-identity transforms,
puncturing downward closure, and free-row high-rate multi-output extensions follow S. Jain,
A. Wills and A. Singh [JWS26]; mid-circuit check-qubit recycling, canonical families
$\mathcal{F}_N^{0,1}$, and bounded-width SAT no-go formulations follow H. Jacinto,
É. Valcarce, M. Barizien, É. Gouzien and N. Sangouard [JVBS26]; non-triorthogonal
$\mathbb{Z}[\omega_8]$ stabilizer collapse, rational complement symmetry, `9T -> 1T` /
`6T -> CCCZ` / `48T -> CCZ` constructions, and $n_{\min}(\tau, d)$ lower bounds follow
M. E. Beverland, E. Campbell, M. Howard and V. Kliuchnikov [Bev20] and A. Balada Gaggioli
and J. Mareček [BGM26]; catalyzed `C2T`, phase catalysis, and asymmetric surface-code
block sizing follow C. Gidney and A. G. Fowler [GF19] and D. Litinski [Lit19]; magic-state
cultivation schedules and exact $O(n^3)$ Clifford-hook simulation follow C. Gidney,
N. Shutty and C. Jones [GJL24], T. Itogawa, Y. Takada, Y. Hirano and K. Fujii [Ito25], and
A. Chan, A. Strikis, Z. Sun and Z. Cai [CSSC26]; open-boundary Caliper mid-circuit escape
postselection follows A. Mehta, V. Menon, H. Zhou, M. D. Lukin and J. P. Bonilla Ataides
[MMZLA26]; and single-shot 3D tricycle $\to$ 4D quadcycle qLDPC `CCZ` factories and AOD
shuttle budgets follow V. Menon, A. Mehta, H. Zhou, M. D. Lukin and J. P. Bonilla Ataides
[MMZLA26q]. Bracketed keys resolve in `REFERENCES.md` at the repository root.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz;
code written by Claude (Anthropic) under his supervision. Documentation is
licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository
root.
