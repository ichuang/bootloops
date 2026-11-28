# ftqc_decoders — GUIDE

KIND: package (Python, 100% stdlib for the default certified tier; `decimal_bounds.py`,
`dem.py`, `bp_osd.py`, `relay_bp.py`, `union_find.py`, `mwpm.py`, `ssf_metacheck.py`,
`tanner_local.py`, `sat_mwd.py`, `effective_distance.py`, `enclosures.py`, `bridges.py`, `certify.py`,
`cli.py`, `adapters/{external_decoders,sinter_bridge}.py`, `selftest.py`, `examples/`).
Battery class: **selftest** (runs in `< 2 s` from a cold clone with zero external
packages and zero floating-point literals; optional external backends `tesseract-decoder`,
`pymatching`, `ldpc`, `chromobius`, `stim`/`sinter`, and `cadical` are probed dynamically
and audited as untrusted solvers when installed).

## PURPOSE

Certified operational syndrome decoders, quantum trapping-set / effective-distance
($d_{\mathrm{eff}}$) verifiers, Stim `.dem` & `sinter` interchange, and non-vacuous
logical-error-rate (`LER`) enclosures for fault-tolerant quantum computing (FTQC).

In a fault-tolerance proof, code distance (`ftqc_codes`) and circuit distance
(`ftqc_audit`, `ftqc_extractors`) bound the lightest *undetected* fault — what an
exponential-time minimum-weight oracle or postselection achieves. Operational FTQC
requires a polynomial-time classical decoder. On degenerate quantum LDPC (`BB72`,
`BB144`, `HGP`) and topological codes, iterative decoders can stagnate or misdecode at
strictly lower fault weight $w < \lfloor(d - 1)/2\rfloor$ on symmetric split-stabilizer
trapping sets, while naive global union bounds $\binom{L}{t+1} p^{t+1}$ saturate
vacuously at `1` on multi-round schedules ($L > 10^4$). `ftqc_decoders` closes those
gaps with replayable `ftqc_decoders-receipt-1` certificates:

1. **Zero-float deterministic stdlib decoders**:
   - `bp_osd` (`bp_osd.py`): scaled-integer Normalized Min-Sum (`NMS` with dyadic
     warm-up) and exact-rational bias Sum-Product BP (`delta = 1 - 2p in Q`, zero
     transcendentals) paired with global Ordered Statistics Decoding (`OSD-0`,
     `OSD-CS(w)`) [Fossorier & Lin 1995; Panteleev & Kalachev 2021; Roffe et al. 2020].
   - `bp_lsd` (`bp_osd.py`): **Localized Statistics Decoding** (`LSD-0`, `LSD-CS(w)`)
     performing on-the-fly `GF(2)` Gaussian elimination inside grown syndrome clusters
     rather than globally ($O(1)$ cluster rank vs $O(n)$ global rank on `BB144`)
     [Hillmann et al. 2025].
   - `relay_bp` (`relay_bp.py`): **Disordered-Memory Relay-BP** breaking symmetric
     quantum trapping sets on bivariate bicycle codes (`BB72`, `BB144`) via rational
     per-variable memory disorder $\gamma_v \in [\gamma_{\min}, \gamma_{\max}]$
     (including negative repulsive weights) relayed across legs, *without* Gaussian
     elimination [Müller et al. 2025; Bravyi et al. 2024].
   - `union_find` & `belief_find` (`union_find.py`): half-edge cluster growth +
     linear-time spanning-tree peeling on graphlike codes [Delfosse & Zémor 2020;
     Delfosse & Nickerson 2021] with cluster-local `GF(2)` elimination on hypergraphs
     [Delfosse, Londe & Beverland 2022].
   - `mwpm` (`mwpm.py`): disjoint-lineage exact graphlike Minimum-Weight Perfect
     Matching via Edmonds' primal-dual Blossom algorithm over exact `Fraction` dual
     potentials with LP complementary-slackness certification [Edmonds 1965; Dennis,
     Kitaev, Landahl & Preskill 2002].
   - `ssf` & `decode_with_metachecks` (`ssf_metacheck.py`): Small-Set-Flip for HGP /
     expander codes [Leverrier, Tillich & Zémor 2015; Fawzi, Grospellier & Leverrier
     2018] and two-stage metacheck syndrome repair ($M \tilde{s} \to \hat{s} \to \hat{e}$)
     for single-shot / mapping-cone qLDPC surgery [Campbell 2019; Cowtan, He,
     Williamson & Yoder 2025], directly consuming `ftqc_extractors/surgery_sequence.py`.
   - `tanner_local` & `simulate_tanner_single_shot_rounds` (`tanner_local.py`):
     iterative local-neighborhood mismatch-reduction decoding for Quantum Tanner codes
     (`arXiv:2206.07750`, `arXiv:2306.12470`), consuming ordered vertex neighborhoods from
     `ftqc_codes` and tracking `coset_reduced_weight` modulo stabilizers across consecutive
     noisy data + syndrome rounds.
   - `sat_mwd` (`sat_mwd.py`): inhomogeneous Tseitin-XOR + Sinz sequential-counter CNF
     exact minimum-weight decoder verified in-tree by `ftqc_codes.lrat` [Sinz 2005;
     Cruz-Filipe et al. 2017; Berent et al. 2024].
2. **Universal Stim `.dem` & `sinter.Decoder` interchange + audited external backends**:
   - Lossless pure-stdlib zero-float `.dem` parser/emitter (`dem.py`), `BootLoopsSinterDecoder`
     (`adapters/sinter_bridge.py`), and untrusted-solver wrappers (`adapters/external_decoders.py`,
     using `scipy.sparse.csc_matrix` $O(\mathrm{nnz})$ check handoffs on high-rate blocks)
     for **Tesseract** (`quantumlib/tesseract-decoder`, $A^*$ + beam-search MLE [Beni et al.
     2025, 2026]), **PyMatching v2** (Sparse Blossom [Higgott & Gidney 2025]), **ldpc v2**
     (`BpOsdDecoder`, `BpLsdDecoder`, `BeliefFindDecoder`), and **Chromobius** [Gidney &
     Jones 2023]. Every candidate mechanism mask returned by an external C++/Rust solver
     is independently audited for syndrome agreement, exact soft weight, and logical
     observable action.
3. **Operational radius, trapping-set witness search & non-vacuous `LER` enclosures**:
   - `certify_operational_radius` & `search_trapping_sets` (`effective_distance.py`):
     certifies operational decoding through weight $s$ or exhibits replayable sub-distance
     trapping-set / finite-order-OSD counterexamples via stabilizer split-support and
     logical-midpoint probes [Poulin & Chung 2008; Raveendran & Vasić 2021].
   - `enclose_exact_shells`, `enclose_peierls_cluster`, `enclose_stratified_shells`
     (`enclosures.py`): outward-rounded `Decimal` enclosures (`ROUND_FLOOR` /
     `ROUND_CEILING`), including connected-cluster Peierls percolation bounds
     [Kovalev & Pryadko 2013; Gottesman 2014] that stay $\ll 1$ where global union
     bounds saturate at `1`, and rare-event weight-shell splitting [Bravyi & Vargo 2013].
   - `bridges.py`: on-demand `Decoder.rule` hook for `ftqc_gadgets` (removing the
     $2^{n-k}$ lookup-table cap on `qr:47`, `bch`, `bb72`), operational open-probe
     decoding on `ftqc_extractors` BB144 `RecordModel`s, and classical reaction-latency
     accounting for `ftqc_uncompute` [Terhal 2015].

## USE-WHEN

- Proving that a concrete polynomial-time decoder achieves operational radius $t_{\mathrm{dec}}$
  on a CSS code or circuit-level `RecordModel`, or finding a sub-distance trapping-set
  witness where `bp_only` / `OSD-1` / `UF` fails below $\lfloor(d - 1)/2\rfloor$.
- Running `ftqc_gadgets` `ec`/`gate`/`close` analyses on codes (`qr:47`, `bch`, `bb72`)
  whose $2^{n-k}$ syndrome table exceeds memory/time budgets.
- Replacing zero-syndrome postselection with operational decoding on `ftqc_extractors`
  BB144 `gross_circuit` probes or `surgery_sequence.py` mapping-cone metachecks.
- Driving external C++ decoders (**Tesseract**, **PyMatching**, **ldpc**) on BootLoops
  circuits/codes via `.dem` while auditing every returned recovery deterministically.
- Computing non-vacuous analytical (`peierls_cluster`, `exact_shells`) or clearly tagged
  empirical (`stratified_shells`) logical-error-rate bounds and `reaction_depth` stalls.

## NOT-FOR

- Multi-core $10^9$-shot Monte-Carlo threshold sweeps in pure Python — export `.dem`
  via `convert-dem` or use `BootLoopsSinterDecoder` / `TesseractSinterDecoder` inside
  `sinter.collect` in an environment with the chosen decoder and `sinter` installed.
- Proving static or circuit *information-theoretic* minimum distance — use `ftqc_codes`,
  `ftqc_audit`, or `ftqc_extractors` for undetected-fault lower bounds.

## EXAMPLE

Run from the BootLoops repository root (supports both `PYTHONPATH=tools python3 -m ftqc_decoders`
and `python3 -m tools.ftqc_decoders`):

```sh
work=$(mktemp -d)

# 1. Discover a weight-3 symmetric split-stabilizer trapping set on BB144 (d=12)
#    under plain BP, and prove Relay-BP, BP+LSD, and BP+OSD all heal it:
python3 -m tools.ftqc_decoders trapping-set --code bb144 --decoder bp_only \
  --compare relay_bp,bp_lsd,bp_osd -o "$work/bb144_ts.receipt.json"
python3 -m tools.ftqc_decoders verify "$work/bb144_ts.receipt.json"

# 2. Check all 277 X-sector patterns of weight zero, one or two on Golay [[23,1,7]].
#    Success proves that BP+OSD-2 corrects through weight two.
python3 -m tools.ftqc_decoders certify-radius --code golay --decoder bp_osd \
  --osd-order 2 --radius 2 -o "$work/golay_r2.receipt.json"
python3 -m tools.ftqc_decoders verify "$work/golay_r2.receipt.json"

# 3. Bound BB144 logical error assuming independent decoding of disjoint clusters
#    and successful decoding of every connected cluster through size five.
python3 -m tools.ftqc_decoders enclose --code bb144 --mode peierls_cluster \
  --p 1/10000 --radius 5

# 4. With stim and tesseract-decoder installed, run audited Tesseract on BB72.
python3 -m tools.ftqc_decoders decode \
  --code bb72 --faults 2,19 --decoder external:tesseract
```

## ACCEPTANCE BATTERY

`python3 tools/ftqc_decoders/selftest.py` (registered in `tools/BATTERIES.json` under
`ftqc_decoders`) runs 8 self-contained legs in `< 2 s`:

- `L0_purity_and_lineage`: AST zero-float/zero-complex scan across all package modules
  and disjoint-lineage check on `mwpm.py`.
- `L1_dem_stim_and_mwpm`: Stim `.dem` round-trip (`repeat` blocks, `shift_detectors`)
  and `MWPM` Edmonds Blossom vs exact bitmask DP + LP strong duality on `Repetition(5)`,
  `Toric 4x4`, and 30 random odd-blossom weighted graphs.
- `L2_bp_osd_lsd_relay`: `NMS`, exact-rational `Sum-Product`, `OSD-CS(2)`, `BP+LSD`,
  and `Relay-BP` across `Steane`, `Golay`, `QR-47`, `BB72`, and `BB144`.
- `L3_uf_ssf_metachecks_sat_lrat`: half-edge `Union-Find`, `BeliefFind`, `SSF`, live
  `ftqc_extractors/surgery_sequence.py` mapping-cone metachecks, and `SAT+LRAT` MWD.
- `L4_trapping_sets_and_deff`: `BB144` split-stabilizer trapping set + `Golay` `OSD-1`
  vs `OSD-2` witness replay.
- `L5_enclosures_bridges_examples_cli`: outward `Decimal` containment, exact shells,
  Peierls cluster bound, Hoeffding shells, `ftqc_gadgets` procedural decoder bridge,
  all packaged `examples/*.json`, and optional live `Tesseract`/`PyMatching`/`ldpc` probes.
- `L6_tanner_single_shot_and_sparse_handoff`: Quantum Tanner `D_8` `[[36,8,3]]` radius-1
  certification (`X` and `Z`), consecutive-round single-shot data + syndrome error simulation
  on `D_8` `[[36,8,3]]` and `S_3` `[[288,8]]`, and sparse `scipy_csc` handoffs on `PP256` and
  `APM1152`.
- `N1_N6_sabotage_controls`: float inputs, invalid metachecks, non-graphlike MWPM,
  forged LRAT proofs, and tampered receipts all refused with expected codes.
- `M1_M4_mutant_kills`: scratch-copy mutations in `bp_osd.py`, `union_find.py`,
  `mwpm.py`, and `decimal_bounds.py` all killed.

## KNOWN-GAPS

- `peierls_cluster` bounds assume disjoint interaction-graph components decode
  independently up to `verified_component_radius`; global `OSD` couples components when
  BP fails globally, so use `peierls_cluster` with local cluster decoders (`union_find`,
  `mwpm`, `bp_lsd`).
- `sat_mwd` proves weight-0/weight-1 optimality natively in Python and invokes `cadical`
  for weight $\ge 2$ syndromes; without `cadical` on `PATH`, weight $\ge 2$ `sat_mwd`
  refuses with `CADICAL-NOT-FOUND`.
- External `chromobius` returns predicted logical observable flips directly without an
  internal mechanism vector, so its receipt marks `syndrome_verified: false`
  (`lineage: external-observable-only`), whereas `tesseract`, `pymatching`, and `ldpc`
  return full mechanism vectors audited against `H * e_hat == s`.

## CREDIT & LITERATURE REFERENCES

All code in `tools/ftqc_decoders/` is a clean-room MIT implementation from the published
literature below; external decoders are invoked out-of-tree via optional Python adapters:

- **Iterative quantum BP & trapping sets**: D. Poulin and Y. Chung, *Quantum Inf. Comput.*
  **8**(10), 987–1000 (2008), [arXiv:0801.1241](https://arxiv.org/abs/0801.1241);
  N. Raveendran and B. Vasić, *Quantum* **5**, 562 (2021),
  [arXiv:2012.15297](https://arxiv.org/abs/2012.15297); J. Chen et al., *IEEE Trans.
  Commun.* **53**(8), 1288–1299 (2005).
- **`BP+OSD`, `BP+LSD`, & `Relay-BP`**: M. P. C. Fossorier and S. Lin, *IEEE Trans. Inf.
  Theory* **41**(5), 1379–1396 (1995); P. Panteleev and G. Kalachev, *Quantum* **5**,
  585 (2021), [arXiv:1904.02703](https://arxiv.org/abs/1904.02703); J. Roffe, D. R.
  White, S. Burton, and E. Campbell, *Phys. Rev. Research* **2**, 043423 (2020),
  [arXiv:2005.07016](https://arxiv.org/abs/2005.07016); T. Hillmann, L. Berent, A. O.
  Quintavalle, J. Eisert, R. Wille, and J. Roffe, *Nature Commun.* **16**, 8214 (2025),
  [arXiv:2406.18655](https://arxiv.org/abs/2406.18655); S. Bravyi et al., *Nature*
  **627**, 778–782 (2024), [arXiv:2308.07915](https://arxiv.org/abs/2308.07915);
  T. Müller et al., [arXiv:2506.01779](https://arxiv.org/abs/2506.01779) (2025).
- **Search, Matching, & Union-Find decoders**: L. A. Beni et al., *"Tesseract: A
  Search-Based Decoder for Quantum Error Correction"*,
  [arXiv:2503.10988](https://arxiv.org/abs/2503.10988) (2025) and
  [arXiv:2602.02985](https://arxiv.org/abs/2602.02985) (2026); J. Edmonds, *Can. J.
  Math.* **17**, 449–467 (1965); E. Dennis, A. Kitaev, A. Landahl, and J. Preskill,
  *J. Math. Phys.* **43**(9), 4452–4505 (2002),
  [quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143); O. Higgott and C. Gidney,
  *Quantum* **9**, 1600 (2025), [arXiv:2303.15933](https://arxiv.org/abs/2303.15933);
  Y. Wu and L. Zhong, [arXiv:2305.08307](https://arxiv.org/abs/2305.08307); Y. Wu et al.,
  [arXiv:2508.04969](https://arxiv.org/abs/2508.04969); C. Gidney and C. Jones,
  [arXiv:2312.08813](https://arxiv.org/abs/2312.08813); N. Delfosse and G. Zémor, *Phys.
  Rev. Research* **2**, 033042 (2020), [arXiv:1703.01517](https://arxiv.org/abs/1703.01517);
  N. Delfosse and N. H. Nickerson, *Quantum* **5**, 595 (2021),
  [arXiv:1709.06218](https://arxiv.org/abs/1709.06218); N. Delfosse, V. Londe, and M. E.
  Beverland, *IEEE Trans. Inf. Theory* **68**(5), 3187–3199 (2022),
  [arXiv:2103.08049](https://arxiv.org/abs/2103.08049).
- **Single-shot `SSF`, Metachecks, SAT/LRAT, Cluster Bounds & Backlog**: A. Leverrier,
  J.-P. Tillich, and G. Zémor, *FOCS* 2015, [arXiv:1504.00822](https://arxiv.org/abs/1504.00822);
  O. Fawzi, A. Grospellier, and A. Leverrier, *FOCS* 2018 / *Commun. ACM* **64**(1)
  (2020), [arXiv:1808.03821](https://arxiv.org/abs/1808.03821); E. T. Campbell, *Quantum
  Sci. Technol.* **4**, 025006 (2019), [arXiv:1805.09271](https://arxiv.org/abs/1805.09271);
  A. Cowtan, Z. He, D. J. Williamson, and T. J. Yoder,
  [arXiv:2510.14895](https://arxiv.org/abs/2510.14895); L. Berent et al., *Quantum*
  **8**, 1506 (2024), [arXiv:2303.14237](https://arxiv.org/abs/2303.14237); C. Sinz,
  *CP* 2005; A. A. Kovalev and L. P. Pryadko, *Phys. Rev. A* **87**, 020304(R) (2013),
  [arXiv:1208.2317](https://arxiv.org/abs/1208.2317); D. Gottesman, *Quantum Inf.
  Comput.* **14**, 1338 (2014), [arXiv:1310.2984](https://arxiv.org/abs/1310.2984);
  S. Bravyi and A. Vargo, *Phys. Rev. A* **88**, 062308 (2013),
  [arXiv:1308.6270](https://arxiv.org/abs/1308.6270); B. M. Terhal, *Rev. Mod. Phys.*
  **87**, 307 (2015), [arXiv:1302.3428](https://arxiv.org/abs/1302.3428).

## Code automorphisms

Select `{"algorithm":"automorphism_ensemble","base":"bp_osd","permutations":[...]}` in a decoder study or pass it with `decode --decoder-spec FILE`. Permutations push original mechanism indices to new indices. A static CSS study can supply `actions` containing verified code-action receipts in place of `permutations`. Each action must preserve single-qubit errors in the selected sector. For DEMs, provide explicit mechanism maps. The verifier checks priors, weights, check row space and the joint check/observable row space. It inverse-transports each candidate, tests its original syndrome, and selects a valid minimum-soft-weight candidate including the identity member. See [reference results](../FTQC.md#code-automorphisms).

The [staged logical gate guide](../ftqc_gadgets/STAGED.md) describes the Steane/RM15 switch, pieceable CCZ, single-fault model and cross-tool replay.
