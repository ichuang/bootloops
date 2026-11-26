# qeccert — GUIDE

KIND: package (Python, stdlib only for the default tier; `gf2.py`, `codes.py`,
`minweight.py`, `certify.py`, `cnf.py`, `lrat.py`, `sat.py`, `leangen.py`,
`selftest.py`, `examples/`). Battery class: partial (the exact leg, the CNF
direction test, the LRAT checker and all planted controls run from the clone;
the SAT leg needs `cadical`, the Lean leg needs `lake` and a built Lean-QEC,
and both skip with a named line otherwise).

## PURPOSE

Certified code facts for CSS quantum error-correcting codes. Input: a pair of
binary check matrices (H_X, H_Z), or a family name. Output: a JSON receipt with
`n`, `k`, the logical distances `d_x`, `d_z`, `d = min(d_x, d_z)`, the
correctable weight `t = (d-1)/2`, the lightest stabilizer weights `wmin_x`,
`wmin_z`, and the structural flags `self_dual`, `doubly_even_x/z`,
`all_ones_logical`. Each number is carried as an interval `[lower, upper]`:
the upper end has a witness word that is checked by dot products alone, the
lower end names the leg that proved it. A value is reported as exact only when
a proven lower bound meets a witnessed upper bound.

Three legs prove lower bounds, with separate lineages:

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

## USE-WHEN

- A fault-tolerance argument takes a code's distance, `t`, or minimum stabilizer
  weight as an input and you want that input to carry its own proof instead of
  a number computed by the same program that consumes it.
- You have a new CSS code (check matrices) and want its parameters with a
  receipt another person can re-verify (`verify`, `verify --full`).
- You want a Lean-QEC distance theorem for a code that Lean-QEC does not ship
  (the generator needs only the matrices).

## NOT-FOR

- Circuit-level fault distance (hook errors, syndrome-extraction schedules,
  measurement faults). The facts here are about the code, not about any circuit.
- Non-CSS stabilizer codes, subsystem codes, qudit codes.
- Large codes. Leg A is exponential in the weight searched; leg B took 186 s for
  QR-47 (d = 11) and is out of reach for the [[889,k,d]] BCH codes; leg C is
  bounded by Lean's `bv_decide` (Lean-QEC reports 90 qubits in the kernel).
  When no leg reaches the claimed value the receipt shows the interval that was
  proven and nothing more.
- Decoder performance. `t` is the bounded-distance radius implied by `d`; the
  statement "weight at most t is corrected exactly" is Lean-QEC's
  `gfidelity_boundedWeightDecoder_eq_one`, which this package cites but does
  not instantiate per code.

## INVOKE

With `tools/` on `PYTHONPATH` (or from `tools/`):

    python3 -m qeccert certify golay -o golay.json
    python3 -m qeccert certify qr:47 --legs A,B --workdir out/ -o qr47.json
    python3 -m qeccert certify code.json --legs A,B,C --target d=6 --workdir out/
    python3 -m qeccert verify golay.json [--full]
    python3 -m qeccert spec bb72 -o bb72.json
    python3 -m qeccert families

`SPEC` is a JSON file (`{"name", "n", "hx": ["0101..."], "hz": [...]}`, rows as
0/1 strings of length n, coordinate 0 first) or `family[:args]`: `steane`,
`golay`, `qr:P`, `toric:L`, `planar:L`, `rotated_surface:D`,
`bb:L,M,x3+y1+y2,y3+x1+x2`, and the named instances `bb18 bb54 gb54 bb70 bb72
bb90 bb108 bb144` (parameters as in the Lean-QEC examples).

`--target q=v` (q in `d`, `d_x`, `d_z`, `wmin_x`, `wmin_z`) stops leg A once
`q >= v` is proven and makes legs B and C prove exactly that bound. Without a
target, legs B and C prove the witnessed upper bound already found.

Library: `from qeccert import load, build, verify, summary`;
`certify.add_lower_bound`, `sat.run_leg`, `leangen.run_leg`.

Exit codes: 0 ok; 1 refusal (not a CSS code, receipt fails verification, two
legs disagree); 2 usage or malformed spec.

ENV: `QECCERT_CADICAL` (solver binary; default `cadical` on PATH),
`QECCERT_LEANQEC` (path to a Lean-QEC checkout with `lake build` done; no
default), `QECCERT_LAKE` (default `lake` on PATH, then `~/.elan/bin/lake`).

## GATES

`python3 tools/qeccert/selftest.py` (default tier, a few seconds):

- GF(2) kernel/rank identities on 200 random matrices.
- Leg A against full enumeration on 300 random (code, subspace) pairs,
  including target stops and capped runs (a capped bound must never exceed the
  truth).
- Family table: Steane, Golay, QR-47, toric 3 and 4, planar 3, rotated surface
  3 and 5, BB18, BB54, GB54 against their known `(n, k, d, wmin)`.
- Receipt controls, each of which must be refused: non-commuting checks, a
  lower bound above a witness, a lower bound with no evidence entry, a forged
  lighter witness, a stabilizer passed off as a logical, a matrix edited under
  a stale hash, wrong `k`, wrong `t`, a forged enumeration claim (caught by
  `--full`), a second leg contradicting a witness.
- CNF direction: for every quantity of five codes the witness extends to a
  model of the formula at its own weight and not below it; on Steane the set
  of error vectors whose induced assignment satisfies the formula equals the
  set satisfying the statement, for every weight cap.
- LRAT checker: accepts a hand-written refutation in text and binary form;
  rejects a missing hint, a non-unit hint, a RAT hint, a deleted or unknown
  hint, a non-increasing id, a fresh variable, a satisfiable CNF, a truncated
  proof.
- SAT leg (needs cadical): all four bounds of five codes proved with a checked
  proof and equal to leg A; a target one above the truth comes back refuted
  with a verified word; a proof with a used lemma removed, and a proof replayed
  against a different CNF, are rejected.
- CLI exit codes.

`--heavy` adds BB54, GB54, BB72 and QR-47 to the SAT leg (about 4 minutes).
`--lean` builds generated files for Steane and BB18 and checks that a false
claim (Steane d >= 4) fails to build.

## FOOTGUNS

- `lower` is what was proven, not what is true. A capped leg A run on a large
  code returns a small `lower` with `budget` set in its evidence entry; read
  the evidence before quoting the number.
- `d_x` is the least weight of an X-type logical (a vector of ker H_Z outside
  rowspace H_X). Papers differ on which of the two is called d_X.
- Row order and redundant rows of H do not change any fact, but they change
  the code hash. Leg C drops dependent rows before generating Lean (the
  theorem is then about the matrices printed in the generated file, and the
  receipt records their hash separately).
- Leg B's soundness rests on one direction of the CNF encoding (every genuine
  word extends to a model). That direction is tested, not proven; leg C is the
  leg with a verified reduction.
- Leg C's theorems do not rest on the three standard axioms alone. Each
  `bv_decide` and `native_decide` call adds an axiom asserting the result of
  compiled evaluation, in the generated file and also inside Lean-QEC library
  lemmas the reduction uses (`dot_product_correct`, `zmod2_dot_bool_step`,
  `zmod2_or_val_ge_left/right` at the pinned commit). The receipt lists all of
  them under `axioms` and groups them under `axiom_classes`; an entry in
  `other` (for example `sorryAx`) means the theorem must not be relied on.
- Recent CaDiCaL prints partial models unless `--modelalllits=true` is given.
  `sat.py` passes the flag when the solver knows it and refuses a partial model
  otherwise.
- `verify` without `--full` re-checks witnesses and consistency only; it does
  not re-prove lower bounds. LRAT and Lean artifacts are re-checked only by
  re-running their legs.

## What runs publicly

Everything in the default tier. The SAT leg runs wherever `cadical` is
installed. The Lean leg needs elan, a Lean-QEC checkout at the commit named in
`leangen.LEANQEC_COMMIT` and its Mathlib cache (about 8 GB).

CREDIT: the SAT-to-distance reduction, the BitVec flattening and the Lean
development that leg C drives are Lean-QEC (M. Ehatamm, Y. Lee, X. Wu, R. Tao,
arXiv:2605.16523, Apache-2.0); generated files follow its example files. The
enumeration bound is that of A. Brouwer and K.-H. Zimmermann as described in
M. Grassl, "Searching for linear codes with large minimum distance" (2006).
The cardinality encoding is C. Sinz (CP 2005). LRAT: L. Cruz-Filipe, M. Heule,
W. Hunt, M. Kaufmann, P. Schneider-Kamp (CADE 2017). Solver: CaDiCaL (A. Biere
et al.). Bivariate bicycle codes: S. Bravyi et al., Nature 627 (2024).
