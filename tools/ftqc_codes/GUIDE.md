# ftqc_codes — GUIDE

Tool page: none yet. Overview of the FTQC tools: [`../FTQC.md`](../FTQC.md).

KIND: package (Python, standard library for the default tier). Construction
and GF(2) algebra: `gf2.py`, `codes.py`, `interchange.py`, `constructions.py`,
`families_highrate.py`. Logical bases, symmetries and actions: `basis.py`,
`symmetry.py`, `automorphisms.py`, `action_domain.py`. Distance legs:
`minweight.py` (A), `cnf.py` + `sat.py` + `lrat.py` (B), `leangen.py` (C),
`bch.py` (D), `cuts.py` (E), `concat.py` (F), `orbit_minweight.py` (G).
`hier.py` (hierarchical decoding), `certify.py` (receipts), `producers/`
(optional exporters run in separate environments), `selftest.py`,
`examples/`. The flat shim `tools/qeccert.py` keeps `import qeccert` working.
Battery class: partial (the SAT leg needs `cadical`, the Lean leg needs `lake`
and a built Lean-QEC; both skip with a named line otherwise).

## PURPOSE

Certified facts about CSS quantum error-correcting codes, as receipts that the
other FTQC tools take as inputs.

- **Code facts** (`ftqc_codes-receipt-1`; archived `qeccert-receipt-1` is
  still read): n, k, d_X, d_Z, d = min(d_X, d_Z), t = ⌊(d−1)/2⌋, the lightest
  stabilizer weights `wmin_x` and `wmin_z`, and the flags `self_dual`,
  `doubly_even_x/z` and `all_ones_logical`. Every number is an interval. The
  upper end has a witness checked with dot products; the lower end names the
  leg that proved it.
- **Constructions** (`ftqc_codes-construction-1`): hypergraph products,
  cyclic and finite-group lifted products, two-block group algebra codes,
  balanced products, quantum Tanner codes, radial codes and pair-partition
  codes (`pp184`, `pp256`, `pp576`, `apm1152`), and imports of sparse JSON and
  MatrixMarket matrices with explicit coordinate maps (`ftqc_codes-import-1`).
  Published claims, upper witnesses, proved lower bounds and simulation
  estimates are kept in separate fields.
- **Logical bases** (`ftqc_codes-basis-1`): named, paired symplectic logical
  slots bound to both the code digest and the basis digest, with signs
  (Y = +iXZ), stabilizer dressings and basis changes.
- **Signed actions**: an ordered word of push permutations, local Clifford
  layers and phase folds, with the signed image of every logical X and Z on a
  named basis. Signed (non-CSS) stabilizer codes are accepted here.

The seven legs for lower bounds have separate lineages:

- **A, exact enumeration** (`minweight.py`): Brouwer–Zimmermann
  information-set enumeration, integer arithmetic only. The default leg.
- **B, SAT** (`cnf.py`, `sat.py`, `lrat.py`): a CNF stating that a logical of
  weight at most d − 1 exists. CaDiCaL refutes it and writes an LRAT proof,
  which the checker in `lrat.py` re-checks against the CNF this package wrote.
  The solver is not trusted.
- **C, Lean** (`leangen.py`): a Lean 4 file in the form of the Lean-QEC
  examples (Ehatamm, Lee, Wu, Tao), built against a Lean-QEC checkout and
  ending in `theorem <name>_dist : d ≤ <css>.toBSM.distance`. The receipt
  records the axioms printed by `#print axioms`.
- **D, BCH bound** (`bch.py`), cyclic codes of odd length: the zeros of the
  generator polynomial are computed exactly in GF(2^m), and the longest run of
  zeros in arithmetic progression gives the bound.
- **E, cuts** (`cuts.py`), k = 1: m pairwise disjoint supports of nontrivial
  Z logicals force every X logical to have weight at least m, and the mirror
  statement.
- **F, concatenation** (`concat.py`, library only): checks that the
  stabilizer spaces are exactly those of the outer code with each position
  replaced by an inner k = 1 code, then applies d ≥ d(inner)·d(outer).
- **G, symmetry orbits** (`orbit_minweight.py`): verifies supplied or
  inferred automorphism generators, certifies a complete partition of the
  qubits into orbits with generator words to each root, and runs the
  orbit-restricted branch-and-bound of Okada and Kasai.

`hier.py` is not a distance leg. It gives the blockwise inner-then-outer
decoder of a concatenated code (`hier.Hier`) and its radius certificate
(`hier.certificate`), which `ftqc_gadgets` uses.

## SCOPE

Binary CSS codes on qubits, given by check matrices or a named family, for
code facts, constructions and bases. Signed stabilizer codes with explicit
logical slots for the `action` command only. Lower bounds come only from the
legs above. A published distance that no leg reaches stays in
`claims.paper_claims`.

## USE-WHEN

- A fault-tolerance argument takes a code's distance, t, lightest stabilizer
  weight, construction or named logical basis as an input, and that input
  should carry its own replayable evidence.
- A high-rate block (LP, BB, radial, quantum Tanner, pair-partition) is passed
  to `ftqc_extractors`, `ftqc_gadgets`, `ftqc_decoders` or `ftqc_compose`,
  which must bind both its code digest and its basis digest.
- A claimed automorphism or transversal gate should be checked with signs
  before it is installed or used in a decoder ensemble.

## NOT-FOR

- Anything about circuits: hook errors, measurement faults, schedules and
  couplers are invisible to static distance. Use `ftqc_faults`,
  `ftqc_extractors` and `ftqc_gadgets`.
- Decoder performance: use `ftqc_decoders`.
- Distance proofs for non-CSS, subsystem or qudit codes. `interchange.py` and
  the pair-partition constructors refuse non-CSS input.

## INVOKE

From the BootLoops root with `export PYTHONPATH="$PWD/tools"`:

    python3 -m ftqc_codes certify golay -o golay.json             # leg A by default
    python3 -m ftqc_codes certify bb18 --legs G -o bb18_orbit.json
    python3 -m ftqc_codes certify qr:47 --legs A,D -o qr47.json
    python3 -m ftqc_codes verify golay.json [--full]                # --full re-runs the searches
    python3 -m ftqc_codes construct pp256 -o pp256.construct.json
    python3 -m ftqc_codes basis lp6_2 --method spectral -o lp6_2.basis.json
    python3 -m ftqc_codes verify lp6_2.basis.json --code lp6_2
    python3 -m ftqc_codes action SPEC steps.json -o action.json
    python3 -m ftqc_codes verify action.json --code SPEC
    python3 -m ftqc_codes spec radial90 -o radial90.json
    python3 -m ftqc_codes families

`SPEC` is a JSON file (`CSSCode` matrices, an `ftqc_codes-import-1` manifest,
an `ftqc_codes-construction-spec-1` spec, or an `ftqc_codes-stabilizer-1`
signed code for `action`) or `family[:args]`. Families: `steane`, `golay`,
`rm15`, `qr:P`, `cyclic:N,l1,...`, `rotated_surface:D`, `planar:L`, `toric:L`,
`bb:L,M,A,B`; named instances `bb18 bb54 gb54 bb70 bb72 bb90 bb108 bb144
radial90 lp6_2 pp184 pp256 pp576 apm1152 qt_d8_36 qt_s3_288`.

`certify` options: `--legs A,B,C,D,E,G|none`; `--target q=v`, which stops the
search for quantity q once q ≥ v is proved; `--max-combos N` (leg A budget,
default 5,000,000); `--workdir` (CNF, LRAT and Lean files); `--sat-timeout`;
`--lean-timeout`. `basis --method binary|candidate|spectral`.

Leg B: put `cadical` on `PATH` or set `FTQC_CODES_CADICAL`. Leg C: set
`FTQC_CODES_LEANQEC` to a Lean-QEC checkout built with
`lake exe cache get && lake build` (the generator follows commit
`3ecc50283caa`), and make `lake` findable (`PATH`, `FTQC_CODES_LAKE` or
`~/.elan/bin/lake`). `python3 tools/ftqc_codes/selftest.py --lean` adds the
Lean tier.

Library: `ftqc_codes.codes.CSSCode`, `ftqc_codes.certify`,
`ftqc_codes.concat.run_leg` (leg F), `ftqc_codes.hier.Hier`,
`ftqc_codes.basis`, `ftqc_codes.automorphisms`.

Exit codes: 0 ok; 1 refusal (not a CSS code, receipt fails verification, two
legs disagree); 2 usage or malformed spec.

## INPUTS

Check matrices as rows of bits, an import manifest with `perm_q`, `perm_hx`
and `perm_hz`, or a construction spec. Bit 0 is coordinate 0. For `action`, a
steps file: a list of `{"perm": [...]}`, `{"fold": [...], "sdg": [...]}` and
`{"layer": [[gate, site, ...], ...]}` items, or
`{"actions": [...], "word": [...]}` to compose verified actions. Input files
are trusted JSON; do not run files from sources you do not trust.

## OUTPUTS

JSON receipts: `ftqc_codes-receipt-1` (code facts with per-quantity `lower`
and `upper`, witnesses and evidence per leg), `ftqc_codes-construction-1`,
`ftqc_codes-basis-1`, and action receipts (`basis`, `steps`, `code_sha256`,
`basis_sha256`, signed `logical_x` and `logical_z` images). `certify` prints
one summary line, for example `QR47: n=47 k=1 d=11 (d_x=11, d_z=11) t=5
wmin_x=12 wmin_z=12 self_dual=True legs=A`.

## GATES

The registered battery runs `python3 tools/ftqc_codes/selftest.py` and then
`python3 -m unittest ftqc_codes.test_discovered_actions`: 3,188 checks and
two named skips (SAT, Lean) in 8.5 s on a four-core Linux virtual machine
without CaDiCaL or Lean.

- GF(2) identities; leg A against full enumeration; hierarchical decoding;
  the CNF direction test; the LRAT checker on banked proofs; BCH; cuts;
  concatenation; CLI exit codes.
- Interchange: qLDPC BB144, a mapped LP(6,2), the Okada–Kasai PP256 `.mtx`,
  the QuantumExpanders.jl D₈ [[36,8,3]]; coordinate-map validation; expansion
  budgets checked before allocation; construction replay.
- Bases: paired-basis verification, slot swaps and stabilizer dressings,
  signed Y checks, Lee's odd-order spectral bases with the even-order
  refusal, cyclic LP seed orbits, Hirasaki–Lee extension maps.
- High-rate families: L = 1 HGP equivalence, finite LP, nonabelian two-block
  group algebra and BB reconciliation, balanced products, D₈ and S₃ quantum
  Tanner codes with `checked`/`assumed` hypothesis tags, pair-partition codes
  with `paper_claims` kept apart from `proved_lower`; leg G against leg A.
- Actions: signed replay of externally discovered generators and of the GB
  phase fold.
- Planted controls, each refused: a lower bound above the witness or without
  evidence, a forged lighter witness, a stabilizer offered as a logical, a
  matrix edited under a stale digest, a wrong k or t, edited cut
  certificates, and a forged enumeration count (caught by `--full`).

## FOOTGUNS

- `lower` is what a leg proved, never a paper claim. Read `claims` before
  quoting a distance for a large family.
- Leg A is the default. For large codes it stops at its combination budget
  and says so; a lower end below the upper end is an open interval, not a
  refutation.
- Leg D rests on the BCH theorem (cited, not kernel-checked) and an exact
  computation of the zeros; leg E on a five-line lemma; leg F on a ten-line
  product lemma. Leg C rests on the Lean kernel and the compiled-evaluation
  axioms it prints.
- Reordering logical slots or dressing representatives by stabilizers keeps
  `code_sha256` and changes `basis_sha256` and the physical Paulis that a
  surgery measures. Consumers must bind both digests.
- Lee's character-packet spectral basis needs an odd group order. For even
  lifts (`LP(8,4)`, `PP256`) use `--method binary`.
- An action receipt is a statement about the code and its basis. It says
  nothing about couplers, timing or faults; `ftqc_gadgets installed-action`
  checks a physical route.
- The optional producers under `producers/` need packages that shadow
  BootLoops names (the PyPI `galois`). Run them with `PYTHONPATH` unset.

## KNOWN GAPS

- Exact distances of the large high-rate families are out of reach; their
  receipts carry upper witnesses and partial lower bounds.
- Lean proofs exist only where leg C was run against a built Lean-QEC; the
  default battery does not run it.
- No distance legs for non-CSS codes.

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
(arXiv:2609.30069), and N. Maskara et al. (arXiv:2610.06749). The optional
producers drive qLDPC and AutQEC.

## Code automorphisms

`action SPEC steps.json -o action.json` certifies an ordered push
permutation, Clifford layers or a phase fold on a named logical basis, and
`verify action.json --code SPEC` regenerates the signed images. The
[worked examples](../FTQC.md#code-automorphisms) include a non-CSS code, a
signed GB fold and optional external discovery. `certify rm15` gives the
punctured Reed–Muller code of the [staged logical gates](../ftqc_gadgets/STAGED.md),
with sites 0 to 6 in the Steane order.

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz;
code written by Claude (Anthropic) under his supervision. Documentation is
licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository
root; third-party components keep their own licenses (THIRD_PARTY.md).
