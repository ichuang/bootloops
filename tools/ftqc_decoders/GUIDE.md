# ftqc_decoders — GUIDE

Overview of the FTQC tools: [`../FTQC.md`](../FTQC.md).

KIND: package (Python standard library for the default tier; no
floating-point literals). Decoders: `bp_osd.py` (BP+OSD and BP+LSD),
`relay_bp.py`, `union_find.py`, `mwpm.py`, `ssf_metacheck.py`,
`tanner_local.py`, `sat_mwd.py`. Certificates and bounds:
`effective_distance.py`, `enclosures.py`, `decimal_bounds.py`. Interchange:
`dem.py`, `adapters/` (Stim/sinter and external decoders). Bridges to the
other tools: `bridges.py`, `surface_memory.py`, `surgery_signatures.py`,
`automorphisms.py`, `staged.py`. `certify.py`, `cli.py`, `selftest.py`,
`examples/`. Optional: `stim`, `sinter`, `tesseract-decoder`, `pymatching`,
`ldpc`, `chromobius`, `cadical`, each probed at run time. Battery class:
selftest.

## PURPOSE

Code distance and circuit distance describe what an exponential-time
minimum-weight decoder or postselection achieves. A real fault-tolerance
scheme uses a polynomial-time decoder, which can fail on patterns well below
half the distance. This package checks concrete decoders and writes
`ftqc_decoders-receipt-1` receipts:

- **Decoders in exact arithmetic**: BP with normalized min-sum on scaled
  integers or sum-product on rationals, followed by OSD (orders 0 and CS(w)) or
  by localized statistics decoding (LSD); disordered-memory Relay-BP;
  half-edge Union-Find and BeliefFind; Edmonds Blossom matching with exact
  rational duals and an LP optimality check; small-set-flip and two-stage
  metacheck repair for single-shot codes and surgery; local decoding for
  quantum Tanner codes; a SAT minimum-weight decoder whose optimality proof is
  re-checked with the LRAT checker of `ftqc_codes`.
- **Radius certificates** (`certify-radius`): every pattern up to weight R is
  decoded and checked; the first failure, if any, is returned as a witness.
- **Trapping-set witnesses** (`trapping-set`): targeted probes (split
  stabilizer supports, logical midpoints) that find patterns a decoder fails
  on, and comparisons showing which other decoders correct them.
- **Logical-error enclosures** (`enclose`): exact weight shells, a
  connected-cluster (Peierls) bound, and sampled weight shells, rounded
  outward and labelled analytical or empirical.
- **Interchange**: a Stim `.dem` parser and writer, a `sinter` decoder class,
  and wrappers that run Tesseract, PyMatching, ldpc and Chromobius and audit
  every returned correction against the syndrome.
- **Bridges**: procedural decoders for `ftqc_gadgets` where lookup tables are
  too large, decoding of `ftqc_extractors` gross-code probes and surgery
  metachecks, single-fault decoding of surface-code memory, decoder-latency
  estimates for `ftqc_uncompute`, automorphism ensembles from verified code
  actions, and the terminal recovery table of the staged pieceable CCZ.

## SCOPE

Static codes (one Pauli sector at a time) and detector error models, with
independent mechanisms of given priors. Radius certificates are exhaustive at
the stated radius. Enclosures hold under the assumptions they print.

## USE-WHEN

- A proof needs a polynomial-time decoder to correct every pattern up to some
  radius, or a witness that it does not.
- `ftqc_gadgets` must analyse a code whose syndrome table is too large
  (`qr:47`, BCH, BB72).
- An external decoder's output should be audited before it is trusted.
- A rough analytical logical-error bound is wanted under explicit assumptions.

## NOT-FOR

- Distance proofs: use `ftqc_codes`, `ftqc_faults` or `ftqc_extractors`.
- Monte Carlo threshold sweeps: export `.dem` (`convert-dem`) or use the
  sinter decoder class in an environment with the decoders installed.

## INVOKE

From the BootLoops root with `export PYTHONPATH="$PWD/tools"` (or
`python3 -m tools.ftqc_decoders` without it):

    python3 -m ftqc_decoders certify-radius --code golay --decoder bp_osd --osd-order 2 --radius 2 -o r.json
    python3 -m ftqc_decoders trapping-set --code bb144 --decoder bp_only --compare relay_bp,bp_lsd,bp_osd -o ts.json
    python3 -m ftqc_decoders enclose --code bb144 --mode peierls_cluster --p 1/10000 --radius 5
    python3 -m ftqc_decoders check STUDY.json -o receipt.json
    python3 -m ftqc_decoders verify receipt.json
    python3 -m ftqc_decoders decode --code bb72 --faults 2,19 --decoder external:tesseract
    python3 -m ftqc_decoders convert-dem --code bb144 -o bb144.dem
    python3 -m ftqc_decoders backends
    python3 -m ftqc_decoders.surface_memory --distance 3 --rounds 3 --basis Z

Exit codes: 0 proved or verified, and also a `trapping-set` run that finds a
witness (read `status`); 1 a radius counterexample or a decode that does not
converge; 2 refusal: malformed or floating-point input, a receipt that does
not replay, a missing optional backend (`CADICAL-NOT-FOUND`).

## OUTPUTS

Receipts with the study, the decoder configuration, `status` (`proved`,
`counterexample`, `trapping_set_found`, ...), `holds_through`, witnesses as
fault-mechanism lists with the decoder's output and the flipped observables,
and for enclosures `ler_lo`, `ler_hi`, the assumption text and the plain union
bound for comparison.

## GATES

`python3 tools/ftqc_decoders/selftest.py` followed by
`python3 -m unittest ftqc_decoders.test_automorphisms`: 10 legs in 2 s on a
four-core Linux virtual machine.

- `L0` no floats or complex numbers in the package; disjoint lineage of
  `mwpm.py`.
- `L1` `.dem` round trip (`repeat`, `shift_detectors`); Blossom against a
  bitmask dynamic program and LP duality on repetition, toric 4×4 and 30
  random odd-blossom graphs.
- `L2` the BP family, OSD-CS(2), BP+LSD and Relay-BP on Steane, Golay,
  QR-47, BB72 and BB144.
- `L3` Union-Find, BeliefFind, SSF, surgery metachecks from
  `ftqc_extractors`, SAT minimum-weight decoding.
- `L4` the BB144 trapping set; the Golay OSD-1 versus OSD-2 witness.
- `L5` outward rounding, exact shells, the Peierls bound, sampled shells, the
  gadget bridge, all examples, optional live external decoders.
- `L6` quantum Tanner radius-1 certificates and single-shot rounds; sparse
  handoff on PP256 and APM1152.
- `N1–N6` floats, invalid metachecks, non-graphlike matching, forged LRAT
  proofs and tampered receipts are refused.
- `M1–M4` mutations in `bp_osd.py`, `union_find.py`, `mwpm.py` and
  `decimal_bounds.py` are caught.

## FOOTGUNS

- **The Peierls bound takes its radius as an input.** `--radius` (the study
  field `verified_component_radius`) is the size of connected clusters the
  decoder is *assumed* to correct; the name notwithstanding, nothing verifies
  it. Establish it separately, for example with `certify-radius` on the same
  decoder, which is feasible only for small codes and radii.
- **The Peierls bound is not always better.** On the BB144 example at
  p = 10⁻⁴ and radius 5 it gives 1.6 × 10⁻¹⁴ against 1.1 × 10⁻¹⁴ for the
  plain union bound C(n, 6)p⁶ printed beside it. Both share the radius
  assumption.
- The cluster bound also assumes that disjoint clusters decode independently.
  Global OSD couples clusters; use it with local decoders (`union_find`,
  `mwpm`, `bp_lsd`).
- A trapping-set receipt is a witness, not a radius: it shows one failure and
  which decoders fix it.
- Sampled-shell enclosures are empirical and labelled so.
- Static-code studies decode one sector under independent flips; they are
  code-capacity statements, not circuit-level ones, unless the input is a
  circuit detector error model.
- **A localized decoder gives up distance.** `{"algorithm": "localized", "base": ..., "detector_range": r}`
  decodes every group of fired detectors on its own, so two faults whose detectors are
  farther apart than r are never corrected together. At a range well below the code
  distance it fails on low-weight sets that the base decoder corrects (`rep7` with
  range 1 fails on two faults). The `enclose_linked_cluster` receipt reports the size
  through which every connected set is corrected; read it before the bounds.
- **The linked-cluster radius is small.** The enclosure holds for activities below
  1 / (e (D + 1) + e² D) with D the degree of the mechanism graph at range r + 1, and a
  nonzero lower bound needs rates far below that (`toric3`: p = 10⁻⁷).

## Linked-cluster enclosure

The study kind `enclose_linked_cluster` (module `linked_cluster.py`) bounds the
probability that an observable parity is flipped from both sides, for a base decoder
localized at `detector_range`. A localized decoder acts independently on fault sets that
are separate at mechanism range r + 1; the module proves this. The flip probability is
then (1 − λ)/2 with log λ a convergent sum over connected clusters of fault sets, whose
first term is −2 times the sum of x^|C| over failing connected sets C. The receipt
enumerates and decodes every connected set through `max_cluster_size`, reports the
failing ones by size and the size through which all are corrected, and gives `ler_lo`
and `ler_hi` from that sum and the remainder 2 M θ^min(m+1, t+2). Nothing about the
decoder is assumed. Faults are independent with one probability on every mechanism.

```json
{"study": "enclose_linked_cluster", "problem": {"kind": "bb72"}, "base_decoder": "bp_osd",
 "detector_range": 2, "max_cluster_size": 2, "physical_p": "1/100000000"}
```

## KNOWN GAPS

- Radius certificates are exhaustive and so limited to small radii (277
  patterns for Golay at radius 2; Golay BP+OSD-2 fails at radius 3).
- `sat_mwd` proves optimality for syndromes of weight at least two only with
  `cadical` installed; otherwise it refuses with `CADICAL-NOT-FOUND`.
- External Chromobius returns observable flips without a mechanism vector, so
  its receipts say `syndrome_verified: false`; Tesseract, PyMatching and ldpc
  return full corrections that are audited.
- The decoder-latency estimate for `ftqc_uncompute` covers single-fault
  syndromes and supplied tick costs.

## CREDIT

The decoders here are implemented from the published literature below;
external decoders are called through optional adapters and are not bundled.

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

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz;
code written by Claude (Anthropic) under his supervision. Documentation is
licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository
root.
