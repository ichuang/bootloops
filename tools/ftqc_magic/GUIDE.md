# ftqc_magic — GUIDE

KIND: package (Python, 100% stdlib for the default certified tier; `rings.py`, `ir.py`,
`matrix.py`, `galois.py`, `limits.py`, `recycle.py`, `search.py`, `cultivate.py`,
`caliper.py`, `qldpc_factory.py`, `pipeline.py`, `catalog.py`, `certify.py`, `cli.py`,
`bridges/{faults_bridge,compose_bridge,extractors_bridge}.py`, `selftest.py`, `examples/`).
Battery class: **selftest** (runs in `< 5 s` from a cold clone with zero external packages
and zero floating-point literals; optional `cadical` and `stim` backends are probed
dynamically for live LRAT proofs and stabilizer-proxy cross-checks when installed).

## PURPOSE

Certified verification, active-qubit compression, cultivation Clifford-hook analysis,
open-boundary escape postselection, single-shot qLDPC factory verification, and
multi-stage pipeline composition for magic-state factories in fault-tolerant quantum
computing (FTQC).

Between 2004 and 2026, magic-state production evolved from Bravyi–Kitaev 15-to-1
Reed–Muller distillation (`2004`) and Bravyi–Haah triorthogonal matrices (`2012`) through
catalyzed `|CCZ> -> 2|T>` conversion and asymmetric surface-code blocks (`2019`), magic
state cultivation (`2024–2025`), and a wave of `2026` advances: mid-circuit check
recycling (`49T -> 1T` in 5 active qubits), high-rate multi-output borrowed-identity
transforms, $\mathrm{GF}(2^s)$ self-dual trace binarization and 3-twirl concatenation,
non-triorthogonal signed-stabilizer collapse over $\mathbb{Z}[\omega_8]$, exact $O(n^3)$
simulation of non-Pauli Clifford hook faults in cultivation, open-boundary Caliper escape
postselection, and single-shot 3D tricycle $\to$ 4D quadcycle qLDPC `CCZ` factories.

Composing these techniques naively introduces subtle soundness hazards — unflagged
cultivation Clifford hooks drop $|T\rangle$ fault distance below $|S\rangle$-proxy
predictions, reused `C2T` catalysts correlate all output $|T\rangle$ errors across rounds,
multi-output factories emit correlated $Z$-error batches, naive depth-2 transversal `CCZ`
halves qLDPC circuit distance, and terminal gap formulas overestimate mid-circuit escape
confidence across open future boundaries. `ftqc_magic` unifies all seven frameworks
behind exact cyclotomic/rational arithmetic and fail-closed `ftqc_magic-receipt-1`
certificates:

1. **Exact zero-float cyclotomic & rational core (`rings.py`, `ir.py`)**:
   - Exact arithmetic in $\mathbb{Z}[\sqrt{2}]$, $\mathbb{Q}(\sqrt{2})$, $\mathbb{Z}[\omega_8]$
     ($\omega_8 = e^{i\pi/4}$), and `PolyQ` ($\mathbb{Q}[p]$), with integer-squaring order
     predicates (`isqrt`) and outward-rounded `Decimal` bounds (`ROUND_CEILING` for
     infidelities and spacetime wire-ticks; `ROUND_FLOOR` for acceptance probabilities).
   - `PhasePolynomialTarget` canonicalizes arbitrary degree-3 diagonal gates over
     $\mathbb{Z}_8$ (`T^k`, `CS`, `CCZ^m`, `CS_{01}CS_{02} ~ CCZ`) and degree-4 `CCCZ`,
     computing exact $4^k$ Pauli squared-modulus multisets for Clifford-orbit matching.
2. **Binary / PPR / borrowed-identity & Galois-qudit distillation (`matrix.py`, `galois.py`)**:
   - Verifies generalized triorthogonality, check decoupling, exact diagonal Clifford
     correction synthesis (`S`, `CZ`), borrowed-identity kernel coset extraction, parent
     capacity $\binom{r}{\le 2} + \binom{r}{\le 1} + 1$, parity lift (`d -> d+1`),
     puncturing downward closure, and Proposition 5 free-row extension [JWS26].
   - Computes exact all-orders `Fraction` acceptance and conditional infidelity via
     syndrome/output coset dynamic programming, cross-checked against `ftqc_codes` /
     CaDiCaL SAT+LRAT distance proofs when available.
   - Constructs $\mathrm{GF}(2^s)$ self-dual trace bases ($\mathrm{Tr}(b_i b_j) = \delta_{ij}$),
     binarizes qudit $U_3$ (`CS`) and $U_7$ (`CCZ`) factories (punctured $\mathrm{GF}(8)$
     Reed–Solomon, $\mathrm{GF}(4)$ Hermitian curve $y^2 + y = x^3$, and `64T -> 2 CCZ`
     `[[64, 6, 4]]` [GPRW26]), and evaluates exact 3-twirl correlated weight transfer.
3. **Non-triorthogonal $\mathbb{Z}[\omega_8]$ collapse & fundamental lower bounds (`limits.py`)**:
   - Implements Theorem 1 of Balada Gaggioli & Mareček [BGM26]: exact closed-form Pauli
     expectations over $\mathbb{Z}[\omega_8]$ certifying $G_{00} > 0$, rank-1
     Cauchy–Schwarz collapse $G_{00} G_{ee} = |G_{0e}|^2$ for all $|e| < d$, rational
     complement symmetry ($e \leftrightarrow e \oplus 1^n$), and target Pauli-spectrum
     matching — certifying the optimal non-triorthogonal `[[9, 1, 2]]` `9T -> 1T`
     ($p_{\mathrm{acc}}(0) = 1/128$, $\epsilon = \frac{27}{2} p^2 + O(p^3)$) and
     `[[6, 4, 1]]` `6T -> CCCZ` protocols alongside general $n_{\min}(\tau, d)$ lower bounds.
4. **Mid-circuit check recycling & LRAT spatial no-go theorems (`recycle.py`, `search.py`)**:
   - Minimizes output/check row weights over $\mathrm{GL}(r, 2)$ cosets and solves
     mid-circuit check-qubit recycling via multi-start greedy + branch-and-bound +
     optional CaDiCaL Sinz-CNF/LRAT lower bounds, compressing `15T -> 1T` from $5 \to 4$
     qubits (with `A <= 3` UNSAT LRAT proved) and Bravyi–Haah `49T -> 1T` ($d = 5$) from
     $14 \to 5$ active qubits, verified by independent physical pool replay [JVBS26].
   - Constructs canonical repetition-code families $\mathcal{F}_N^0$ (`T`) and
     $\mathcal{F}_N^1$ ($\sqrt{T}$), solves $S_N$-symmetric families (`45T -> 1T` at
     $N = 9, d = 3$), and replays banked LRAT UNSAT refutations via `ftqc_codes.lrat`
     proving $N = 4 \implies d \le 2$ for `T` and $N \le 5 \implies d \le 2$ for `CCZ`.
5. **Exact $O(n^3)$ Clifford-hook cultivation verifier (`cultivate.py`)**:
   - Propagates circuit Pauli faults across $T / T^\dagger$ checks into non-Pauli Clifford
     errors $H_{\pm} = (X \pm Y)/\sqrt{2}$ and implements Algorithm 1 + Proposition 1 of
     Chan et al. [CSSC26] over $\mathbb{Q}(\sqrt{2})$ in $O(n^3)$ binary symplectic time.
   - Reproduces exact single/double/triple malignant fault counts on `MSC-3` (`mal_T[2] = 4`,
     `mal_S[3] = 82`), proving that unflagged double-checking degrades from $d_S = 3$ on
     $|S\rangle$ to $d_T = 2$ on $|T\rangle$, while two $Z$-flag ancillas (`msc3_z_flagged`)
     or Itogawa et al.'s 2-round short check (`itogawa_2round`) restore full $d_T = 3$
     (cross-checked on $|S\rangle$ against `ftqc_faults.audit` and `stim`).
6. **Open-boundary Caliper escape & single-shot qLDPC factories (`caliper.py`, `qldpc_factory.py`)**:
   - Partitions escape DEMs into $D_{\mathrm{comm}} \sqcup D_{\mathrm{open}}$, computes
     parity-lifted shortest paths $d_b(u, v)$, and certifies worst-case soft-weight
     complementary gaps across hidden open-boundary syndromes (cross-checked against
     `ftqc_decoders.mwpm.decode_mwpm`), plus exact adaptive escape spacetime volume [MMZLA26].
   - Verifies 3D tricycle chain complexes ($H_X H_Z^T = 0$, $M_X H_X = 0$), depth-1
     $d = 2$ repetition fan-out restoration vs. depth-2 distance halving, 4D quadcycle
     mapping-cone extension ($M_X H_X = 0$, $M_Z H_Z = 0$) with live single-shot syndrome
     repair via `ftqc_decoders.ssf_metacheck.decode_with_metachecks`, homomorphic
     $|{\det}\rangle \mapsto |\mathrm{CCZ}\rangle^{\otimes 2}$ minor collapse, and
     monotone AOD wrap-class shuttle counts $\Theta^{(N)}$ [MMZLA26q].
7. **Pipeline composition, catalysis & sibling bridges (`pipeline.py`, `bridges/`)**:
   - Verifies Gidney–Fowler `C2T` (`|CCZ>|T> -> |T>^{\otimes 3}`) exactly over
     $\mathbb{Z}[\omega_8]$, sizes asymmetric Litinski $(d_X, d_Z, d_m)$ surface-code
     patches, composes multi-stage cultivation/escape/distillation/catalysis/qLDPC
     pipelines with outward rational rounding at inter-stage boundaries, and exports
     read-only contracts for `ftqc_compose`, `ftqc_extractors`, and `ftqc_uncompute`.

## USE-WHEN

- Certifying exact distance, leading error polynomial, active-qubit schedule, and
  fundamental $n_{\min}(\tau, d)$ optimality gap of a `T`, `CS`, `CCZ`, `CCZ^m`, or `CCCZ`
  magic-state factory (triorthogonal, PPR, $\mathrm{GF}(2^s)$ trace, or non-triorthogonal
  signed stabilizer).
- Compressing logical check qubits of a PPR factory (`15-to-1`, `49-to-1`, `64-to-1`,
  `47-to-CCZ`) via mid-circuit measurement and reset.
- Auditing whether a magic-state cultivation check schedule suffers non-Pauli Clifford
  hook degradation ($d_T < d_S$) that $|S\rangle$-state stabilizer simulations miss, or
  sizing mid-circuit Caliper escape postselection across open future boundaries.
- Composing hybrid multi-stage pipelines (e.g., `Z-flagged cultivation -> Caliper escape ->
  recycled 15-to-1` or `cultivation -> 8T-to-CCZ -> catalyzed C2T`) with automatic
  enforcement of correlation and catalyst-poisoning rules, or finding Pareto-optimal
  pipelines for a target `(p_phys, target_infidelity, max_active_qubits)` budget.

## NOT-FOR

- Full $10^9$-shot Monte-Carlo threshold sweeps of surface-code memory or escape stages
  — use `ftqc_decoders` (`BootLoopsSinterDecoder` / `sinter`) for large-scale sampling
  and supply the resulting physical rate to `ftqc_magic`.
- Logical reversible cleanup (`HMR` / spooky pebbling) — use `ftqc_uncompute` and feed it
  the `ftqc_uncompute_pricing` record emitted by `ftqc_magic`.

## CLI & RECEIPT SCHEMA

Run from the BootLoops repository root (`PYTHONPATH=tools python3 -m ftqc_magic` or
`python3 -m tools.ftqc_magic`):

```sh
export PYTHONPATH=tools
work=$(mktemp -d)

# Build and independently verify a source-pinned ftqc_magic-receipt-1 certificate
python3 -m ftqc_magic check tools/ftqc_magic/examples/flagged_cultivation_caliper_15to1.json \
  -o "$work/pipe.receipt.json"
python3 -m ftqc_magic verify \
  tools/ftqc_magic/examples/flagged_cultivation_caliper_15to1.json "$work/pipe.receipt.json"

# Compress active logical qubits on a catalog factory (--lrat invokes CaDiCaL if installed)
python3 -m ftqc_magic recycle --factory bh_49_1_5_t

# Run exact O(n^3) Clifford-hook cultivation analysis (--external runs ftqc_faults & Stim)
python3 -m ftqc_magic cultivate --variant msc3_z_flagged --p 1/1000

# List all 17 built-in landmark factories or suggest Pareto-optimal pipelines
python3 -m ftqc_magic catalog
python3 -m ftqc_magic suggest --target CCZ --p 1/1000 --infidelity 1/1000000000
```

Every `ftqc_magic-receipt-1` JSON document carries `schema`, `kind`, `spec_sha256`,
`source_sha256` (hash of all package `.py` sources), `contract` (target phase polynomial,
`distance`, `active_qubits`, `correlation_mode`, outward-rounded `acceptance_prob_lo`,
`conditional_infidelity_hi`, `expected_wire_ticks_per_output_hi`, and explicit
`premises`), `optimality_gap` (against proven $n_{\min}(\tau, d)$ lower bounds where
applicable), `exports` (`ftqc_compose`, `ftqc_extractors`, `ftqc_uncompute` bridge
records), and stage-specific `details`.

## KNOWN HAZARDS

1. **Unflagged cultivation Clifford-hook degradation (`N1: UNFLAGGED-CULTIVATION-HOOK-DEGRADATION`)**:
   Propagating a single Pauli fault halfway through a transversal $T / T^\dagger$ check
   creates $H_{\pm} = (X \pm Y)/\sqrt{2}$, spreading across subsequent CNOTs into a
   weight-2 undetectable logical error on $|T\rangle$ even when the Clifford $|S\rangle$
   proxy achieves $d_S = 3$ (`msc3_unflagged`: $d_T = 2$ vs $d_S = 3$). `cultivate.py`
   refuses any schedule whose exact $\mathbb{Q}(\sqrt{2})$ $|T\rangle$ fault distance falls
   below `claimed_t_distance`.
2. **Catalyst error poisoning (`N2: CATALYST-POISONING-HAZARD`)**: Gidney–Fowler `C2T` and
   phase-catalysis ladders reuse a persistent catalyst state; a single undetected error on
   the catalyst flips *every* subsequent output, violating the independent-input premise
   of downstream distillation and multi-shot `ftqc_compose` tail bounds. `CATALYST_POISONED`
   contracts are strictly terminal (`terminal_single_shot=True`).
3. **Multi-output batch correlation (`N3: MULTI-OUTPUT-CORRELATION-HAZARD`)**: Block
   factories with $k > 1$ outputs (`14T -> 2T`, `63T -> 2T`, `64T -> 2 CCZ`, qLDPC
   $|{\det}\rangle \mapsto 2\,\mathrm{CCZ}$) suffer joint $Z$-error patterns across outputs
   within the same accepted block (`BATCH_CORRELATED`). Feeding them into a downstream
   distillation stage or independent `cone-bound` study requires explicit
   `interleaving="distinct_instances"` (or `"distinct-instances"`) or `"joint_3twirl"`.
4. **Depth-2 transversal `CCZ` distance halving (`N4: DEPTH2-TRANSVERSAL-CCZ-HALVING`)**:
   On 3D tricycle qLDPC codes without a $d = 2$ repetition fan-out outer layer, physical
   `CCZ` gates share qubits across depth 2 and propagate faults to halve circuit distance
   to $\lceil d/2\rceil$.
5. **Open-boundary naive gap overconfidence (`N5: OPEN-BOUNDARY-NAIVE-GAP-HAZARD`)**:
   Evaluating terminal MWPM complementary gaps with $s_{\mathrm{open}} = 0$ during
   mid-circuit escape ignores hidden defects on open future time boundaries that erode the
   true conditional gap; `caliper.py` refuses `naive_terminal_gap` whenever
   $D_{\mathrm{open}} \ne \emptyset$.

## ACCEPTANCE BATTERY

`python3 tools/ftqc_magic/selftest.py` (registered in `tools/BATTERIES.json` under
`ftqc_magic`, `cwd: root`) executes in `< 5 s` from a cold clone with zero external
dependencies:

- `L0_purity_and_pinned_receipts`: AST scan across all `tools/ftqc_magic/**/*.py` files
  forbidding `float`/`complex` literals and `float()`/`complex()` calls, plus independent
  replay verification of all 4 banked `examples/*.receipt.json` fixtures.
- `L1_rings_and_ir`: $\mathbb{Z}[\sqrt{2}]$, $\mathbb{Q}(\sqrt{2})$, $\mathbb{Z}[\omega_8]$,
  `PolyQ`, outward `ROUND_CEILING`/`ROUND_FLOOR` soundness, and Clifford-orbit matching
  (`CS_{01}CS_{02} ~ CCZ`).
- `L2_binary_galois_and_transforms`: Landmark binary/PPR factories (`bk_15_1_3_t`,
  `bh_14_2_2_t`, `bgm_36_3_3_cs2_ccz`, `bgm_48_3_4_ccz_2level`, `jws_63_2_3_t2`),
  borrowed-identity kernel transform, parity lift, puncturing downward closure,
  $\mathrm{GF}(2^s)$ self-dual trace binarization (`gf8_rs`, `gf4_hermitian_ag`,
  `gprw_64_6_4_2ccz`), and 3-twirl global vs. local concatenation.
- `L3_collapse_and_lower_bounds`: Exact $\mathbb{Z}[\omega_8]$ stabilizer collapse on
  `bgm_9_1_2_t` (`1/128`, `27/2 p^2`), `bgm_6_4_1_cccz`, `bgm_3_2_1_cs`, and $n_{\min}$
  optimality bounds.
- `L4_recycling_and_spatial_nogo`: Active-qubit compression on `15-to-1` ($5 \to 4$) and
  `49-to-1` ($14 \to 5$), canonical $\mathcal{F}_N^0$, $S_9$ `45T -> 1T`, and banked LRAT
  spatial no-go refutations (`(4, "T", 3)` and `(5, "CCZ", 3)`) via `ftqc_codes.lrat`
  (plus live CaDiCaL LRAT when installed).
- `L5_cultivation_clifford_hooks`: Exact $O(n^3)$ `msc3_unflagged` ($d_T = 2, d_S = 3$,
  `4` vs `82`) vs. `msc3_z_flagged` ($d_T = 3$) and `itogawa_2round` ($d_T = 3$), plus
  live `ftqc_faults.audit` and `stim` cross-checks when installed.
- `L6_caliper_and_qldpc`: Caliper worst-case open-boundary gap (+ Blossom cross-check) &
  adaptive escape volume; 3D tricycle + 4D quadcycle metacheck repair +
  $|{\det}\rangle \mapsto |\mathrm{CCZ}\rangle^{\otimes 2}$ + AOD shuttle budget.
- `L7_pipeline_catalysis_bridges_and_cli`: Gidney–Fowler `C2T` $\mathbb{Z}[\omega_8]$
  proof, `compose_magic_pipeline`, `suggest_pipelines`, sibling bridges, and CLI
  `check`/`verify` subprocess round-trip.
- `N1_N6_sabotage_controls`: Deterministic fail-closed refusal checks across all six
  sabotage classes (`N1`–`N6`).
- `M1_M4_mutant_kills`: Four isolated-tempdir subprocess code-mutation self-audits across
  `limits.py`, `cultivate.py`, `pipeline.py`, and `rings.py` confirming zero surviving
  mutants.
- Expected summary line: `SELF-TEST: ALL PASS (rc 0)`.

CREDIT: Clean-room MIT implementation authored from published mathematical equations and
circuit specifications (zero third-party source lines bundled). Punctured Reed–Muller and
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

### Delivery to encoded consumers

Use `ftqc_magic.delivery.read_delivery(spec, receipt)` for a strict replay
with normalized phases, ordered state groups, Clifford corrections,
acceptance, correlations and rational resource bounds. For binary producers
the two mask tables give the joint mass with acceptance and its conditional
law. Other supported producers can carry a state bound alone. Cost per
batch, per grouped quantum state and per native output are recorded
separately. For example, a binary CCZ can have three native output wires
but only one deliverable three-qubit state. Strict v1 receipts bind the
producer Python sources. To use a historical spec at a new revision,
reissue its receipt with `python3 -m ftqc_magic check`.
