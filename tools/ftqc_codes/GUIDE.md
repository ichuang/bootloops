# ftqc_codes — GUIDE

KIND: package (Python, stdlib only for the default tier; `gf2.py`, `codes.py`,
`interchange.py`, `constructions.py`, `basis.py`, `symmetry.py`,
`families_highrate.py`, `orbit_minweight.py`, `minweight.py`, `certify.py`,
`cnf.py`, `lrat.py`, `sat.py`, `bch.py`, `cuts.py`, `concat.py`, `hier.py`,
`leangen.py`, `producers/qldpc_export.py`, `selftest.py`, `examples/`; flat shim
`tools/qeccert.py` preserves legacy `import qeccert` callers). Battery class:
partial (the exact leg, symmetry-orbit leg, construction/interchange/basis
verifiers, CNF direction test, LRAT checker, and all planted controls run from
the clone; the SAT leg needs `cadical`, the Lean leg needs `lake` and a built
Lean-QEC, and both skip with a named line otherwise).

## PURPOSE

Certified code facts, replayable algebraic construction ancestry, and verified
addressable logical bases for CSS quantum error-correcting codes.

1. **Static code-fact receipts (`ftqc_codes-receipt-1`; legacy `qeccert-receipt-1`
   also accepted)**: reports `n`, `k`, logical distances `d_x`, `d_z`,
   `d = min(d_x, d_z)`, correctable weight `t = (d-1)//2`, lightest stabilizer
   weights `wmin_x`, `wmin_z`, and structural flags (`self_dual`,
   `doubly_even_x/z`, `all_ones_logical`) as intervals `[lower, upper]` backed
   by dot-product-verified witnesses (`upper`) and named proof legs (`lower`).
2. **Sparse interchange & construction receipts (`interchange.py`,
   `constructions.py`, `families_highrate.py`)**: imports sparse JSON and
   MatrixMarket (`.mtx`) CSS matrices (`ftqc_codes-import-1`) with explicit push
   coordinate maps (`perm_q`, `perm_hx`, `perm_hz`), expands HGP, cyclic and
   general finite-group lifted products (`2BGA`, balanced products, Quantum
   Tanner left-right Cayley complexes, radial codes, and pair-partition CPM/APM
   high-rate codes `PP184`, `PP256`, `PP576`, `APM1152`), and emits
   `ftqc_codes-construction-1` receipts keeping `paper_claims`, `witness_upper`,
   `proved_lower`, and `simulation_estimates` strictly separate.
3. **Addressable paired logical bases & lift symmetries (`basis.py`,
   `symmetry.py`)**: derives or verifies named symplectic logical slot pairs
   (`ftqc_codes-basis-1`) bound to both `code_sha256` and `basis_sha256`,
   supports exact mod-4 signed Pauli (`Y = +iXZ`) embeddings, stabilizer
   dressings and symplectic basis conversions (`verify_basis_transform`), Lee's
   odd-order Abelian character-packet spectroscopy (`arXiv:2607.05386`),
   canonical cyclic LP seed-orbit coverage (`arXiv:2607.28605`), and
   Hirasaki–Lee cyclic extension pullback/pushforward chain maps
   (`arXiv:2607.28621`).

Seven legs prove lower bounds, with separate lineages:

- **A, exact enumeration** (`minweight.py`): Brouwer-Zimmermann information-set
  enumeration over disjoint column sets, integer arithmetic only.
- **B, SAT** (`cnf.py`, `sat.py`, `lrat.py`): a CNF stating "a word of weight at
  most d-1 exists"; CaDiCaL refutes it and emits an LRAT proof, and the
  in-package checker `lrat.py` re-checks that proof against the CNF this
  package wrote. The solver is not trusted.
- **C, Lean kernel** (`leangen.py`): generates a Lean 4 file in the form of the
  Lean-QEC examples (Ehatamm, Lee, Wu, Tao, arXiv:2605.16523) and builds it
  against a Lean-QEC checkout, ending in
  `theorem <name>_dist : d ≤ <css>.toBSM.distance`. The receipt records the
  axioms the theorem depends on, as printed by `#print axioms`.
- **D, BCH bound** (`bch.py`), for cyclic codes of odd length: the generator
  polynomial and its zero set are computed exactly in GF(2^m), and the longest
  run of zeros in arithmetic progression gives the bound.
- **E, cut certificate** (`cuts.py`), for codes with k = 1: m pairwise disjoint
  supports of nontrivial Z-type logical operators force every X-type logical to
  have weight at least m (and the mirror statement).
- **F, concatenation** (`concat.py`): checks that a code's stabilizer spaces
  are exactly those of "outer code with each position replaced by an inner
  k = 1 code", then applies the product bound d >= d(inner) * d(outer).
- **G, symmetry-reduced orbit enumeration** (`orbit_minweight.py`): verifies
  supplied or inferred code automorphism generators (`H_X` and `H_Z` rowspace
  invariance), certifies a complete qubit-orbit partition with explicit
  generator-word witnesses to each canonical orbit root, and runs progressive
  orbit-restricted syndrome branch-and-bound (`arXiv:2607.14091`). Cheap
  `verify` re-checks every automorphism and orbit-word witness; `verify --full`
  also re-runs the branch-and-bound search.
- **Hierarchical decoding** (`hier.py`, not a distance leg): blockwise inner +
  outer decoder (`hier.Hier`) and radius certificate (`hier.certificate`).

## USE-WHEN

- A fault-tolerance argument takes a code's distance, `t`, minimum stabilizer
  weight, construction ancestry, or named logical basis as an input and you want
  that input to carry its own replayable receipt.
- You want to pass an addressable high-rate CSS block (`LP`, `BB`, `radial90`,
  `QuantumTanner`, `PP256`, `APM1152`) into `ftqc_gadgets`, `ftqc_extractors`,
  `ftqc_decoders`, or `ftqc_compose` with verified `code_sha256` and
  `basis_sha256` digests.

## NOT-FOR

- **Circuit-level fault equivalence, ZX gadget synthesis, or stochastic gadget
  contracts**: static code facts do not count hook errors, measurement flips, or
  coupler schedules. Use **`ftqc_faults`**, **`ftqc_decoders`**,
  **`ftqc_gadgets`**, **`ftqc_extractors`**, and **`ftqc_compose`**.
- Non-CSS stabilizer codes, subsystem codes, qudit codes (`interchange.py` and
  `cpm_pair_partition` refuse non-CSS/symplectic inputs fail-closed).

## INVOKE

With `tools/` on `PYTHONPATH`:

    python3 -m ftqc_codes certify golay -o golay.json
    python3 -m ftqc_codes certify bb18 --legs G -o bb18_orbit.json
    python3 -m ftqc_codes construct pp256 -o pp256.construct.json
    python3 -m ftqc_codes basis lp6_2 --method spectral -o lp6_2.basis.json
    python3 -m ftqc_codes verify lp6_2.basis.json --code lp6_2
    python3 -m ftqc_codes verify golay.json [--full]
    python3 -m ftqc_codes spec radial90 -o radial90.json
    python3 -m ftqc_codes families

`SPEC` is a JSON file (`CSSCode` matrix JSON, `ftqc_codes-import-1` manifest, or
`ftqc_codes-construction-spec-1` spec) or `family[:args]`, including named
instances `bb18 bb54 gb54 bb70 bb72 bb90 bb108 bb144 radial90 lp6_2 pp184 pp256
pp576 apm1152 qt_d8_36 qt_s3_288`.

Exit codes: 0 ok; 1 refusal (not a CSS code, receipt fails verification, two
legs disagree); 2 usage or malformed spec.

## GATES

`python3 tools/ftqc_codes/selftest.py` (default tier, ~5 s, >3,200 checks):

- GF(2) identities, Leg A vs full enumeration, hierarchical decoding, CNF/LRAT,
  BCH, cuts, concatenation, SAT (`cadical`), and CLI exit codes.
- Sparse/MatrixMarket interchange (`qLDPC` `BB144`, mapped `LP(6,2)`,
  Okada–Kasai `PP256` `.mtx`, `QuantumExpanders.jl` `D_8` `[[36,8,3]]`),
  coordinate-map validation, pre-allocation expansion budgets, and
  `ftqc_codes-construction-1` replay.
- Paired basis verification (`ftqc_codes-basis-1`), slot-swap & stabilizer
  dressing transforms, mixed-Pauli (`Y`) sign checks, Lee odd-order Abelian
  spectroscopy (including even-order refusal), cyclic LP seed-orbit coverage,
  and Hirasaki–Lee extension chain maps.
- High-rate families (`L=1` HGP equivalence, Raveendran finite LP, nonabelian
  `2BGA` & BB reconciliation, nonabelian balanced product, `D_8`/`S_3` Quantum
  Tanner with `"checked"`/`"assumed"` hypothesis tags, `PP184`/`PP256`/`PP576`/
  `APM1152` with `paper_claims` separated from `proved_lower`) and Leg G orbit
  search cross-checked against Leg A.

## FOOTGUNS

- `lower` (`proved_lower`) is what a proof leg certified, never a paper claim.
  Published distances that exceed exact search budgets live only under
  `claims.paper_claims` until a leg proves them.
- Changing logical slot ordering or dressing representatives by stabilizers
  preserves `code_sha256` while changing `basis_sha256` and the physical Pauli
  operators measured by downstream surgery; consumers must bind both digests.
- Lee's character-packet spectral basis requires odd group orders (`gcd(|G|,2)=1`);
  for even lifts (`LP(8,4)`, `PP256`), use `derive_paired_basis` (`--method binary`).

CREDIT: Lean-QEC SAT/BitVec distance reduction: M. Ehatamm, Y. Lee, X. Wu,
R. Tao (arXiv:2605.16523, Apache-2.0). Information-set enumeration: A. Brouwer
and K.-H. Zimmermann (Grassl 2006). Sinz sequential counter (CP 2005); LRAT:
Cruz-Filipe et al. (CADE 2017); CaDiCaL: A. Biere et al. High-rate QEC
constructions and symmetries: J.-P. Tillich & G. Zémor (arXiv:0903.0566),
P. Panteleev & G. Kalachev (arXiv:2012.04068, arXiv:2111.03654), N. P.
Breuckmann & J. N. Eberhardt (arXiv:2012.09271), A. Leverrier & G. Zémor
(arXiv:2202.13641), H.-K. Lin & L. P. Pryadko (arXiv:2306.16400), S. Bravyi et
al. (Nature 627, 2024), T. R. Scruby, T. Hillmann & J. Roffe (arXiv:2406.14445),
N. Raveendran, D. Declercq & B. Vasić (arXiv:2503.07567), J. Y. Lee
(arXiv:2607.05386), K. Okada & K. Kasai (arXiv:2607.14091), H. Zheng et al.
(arXiv:2607.28605), Y. Hirasaki & J. Y. Lee (arXiv:2607.28621), L. Berent,
L. Z. Cohen & A. O. Quintavalle (arXiv:2609.11723), J. Y. Lee et al.
(arXiv:2609.30069), and N. Maskara et al. (arXiv:2610.06749).

## Code automorphisms

Use `action SPEC steps.json -o action.json` to certify an ordered push permutation, Clifford layers, or a phase fold on a named logical basis. Use `verify action.json --code SPEC` to regenerate the signed images. The action command also accepts signed stabilizer codes with explicit logical slots. The [worked examples](../FTQC.md#code-automorphisms) include a non-CSS code, a signed GB fold and optional external discovery.

The [staged logical gate guide](../ftqc_gadgets/STAGED.md) describes the Steane/RM15 switch, pieceable CCZ, single-fault model and cross-tool replay.
