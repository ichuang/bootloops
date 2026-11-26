# faultcert — GUIDE

KIND: package (Python, stdlib only; imports `qeccert` from this repository for
the SAT leg and for code constructions; `circuit.py`, `stab.py`, `model.py`,
`analysis.py`, `certify.py`, `library.py`, `stimleg.py`, `selftest.py`).
Battery class: partial (everything exact runs from the clone; the SAT leg needs
`cadical`, the Stim leg needs the `stim` Python package; both skip with a named
line otherwise).

## PURPOSE

Circuit-level fault properties of Clifford gadgets, as receipts. `qeccert`
says what a code's distance is; this package says what faults inside a circuit
do. Input: a circuit of preparations, H, CX, idles and measurements in ticks,
optionally with declared detectors, logical observables and output qubits.
Fault model: one fault is one location misbehaving (any non-identity Pauli
after a preparation, idle, H or CX; a flipped outcome at a measurement), and
the weight of a fault set is its number of locations.

Outputs:

- **fault model**: for every location and every Pauli, which detectors fire,
  which observables flip, and whether the output state differs from the one the
  reported outcomes imply (residual error and wrong Pauli frame together);
- **circuit distance**: least number of faults that fire no detector and flip a
  logical observable, as an interval with a replayable witness;
- **fault distance**: the same with "harm" widened to "the output state is
  wrong", for gadgets whose product is a state;
- **strict fault tolerance through order s**: every set of j <= s faults that
  fires no detector leaves the output within Pauli weight j of the believed
  state, or the first counterexample;
- **heavier than w**: least number of faults that fire no detector and leave
  the output further than Pauli weight w from the believed state, with the
  fault set (for a purification step: how many faults an accepted heavy error
  needs);
- **census**: how many locations are harmless, how many can leave only a
  single-qubit error (per output qubit), how many can leave a heavier one.

Why exact answers are possible: for Pauli faults in a Clifford circuit the
effect of a fault set is the XOR of the effects of its members, and the set of
effects one location can produce is closed under XOR. "At most j faults"
therefore means "a sum of at most j columns" of a GF(2) matrix whose columns
are the distinct single-location effects. Circuit distance is a minimum-weight
problem of the same shape as code distance, which is why `qeccert`'s CNF
encoding and LRAT checker are reused unchanged.

Legs:

- **A, enumeration** (`analysis.py`): meet in the middle on the detector
  signature; exact up to the weight the cap allows.
- **B, SAT** (`certify._leg_b`): the statement as a `qeccert.cnf.Formula`;
  CaDiCaL's LRAT proof is re-checked by `qeccert.lrat`; a SAT answer is turned
  into a fault set and replayed.
- **S, Stim** (`stimleg.py`): the circuit is exported with the same noise and
  detectors; Stim must accept the detectors as deterministic, its detector
  error model must have the same set of symptoms as this package's columns,
  and any logical error it finds must not be lighter than a bound proven here.
  Stim is a cross-check from a separate lineage, not a source of bounds.

## USE-WHEN

- An argument uses a count or a property of fault locations in a concrete
  gadget (how many locations can leave a multi-qubit error, whether a
  verification step catches every single fault) and you want it recomputed
  from the circuit by a second program.
- You changed a CX schedule and want to know whether hook errors now line up
  with a logical operator (`surface:D,R` against `surface:D,R,z,bad`).
- You have a preparation-with-verification circuit and want the lowest order
  at which an accepted fault set leaves too heavy an error, with the fault set.

## NOT-FOR

- All-orders statements. Everything here is "for all fault sets of at most j
  locations" with j small; a probability bound still needs an argument such as
  a cone or counting theorem on top.
- Non-Clifford gates, and classically controlled operations (decoders and
  adaptive corrections inside the circuit). Postselection is supported as
  "no detector fires".
- Leakage, coherent errors, correlated faults across locations.
- Large circuits at high weight: leg A enumerates subsets of the distinct
  fault effects (a few thousand of them for a distance-5 surface patch), so it
  proves "no harm below weight w" for w up to about 5 there; leg B goes
  somewhat further.

## INVOKE

With `tools/` on `PYTHONPATH` (or from `tools/`):

    python3 -m faultcert analyze surface:3,2 --legs A,B,S -o r.json
    python3 -m faultcert analyze surface:3,2,z,bad
    python3 -m faultcert analyze cat:4,verify --strict 2 --heavier 1,2
    python3 -m faultcert analyze mycircuit.json --max-weight 5 --strict 2
    python3 -m faultcert verify r.json --full
    python3 -m faultcert export css_prep:golay --stim

Circuit JSON: `{"name", "n", "ticks": [[["RX",0],["R",1]], [["CX",0,1]], ...],
"detectors": [[measurement indices], ...] | null, "observables": [[...]],
"outputs": [qubits]}`. Operations: `R` (|0>), `RX` (|+>), `I`, `H`, `CX c t`,
`M` (Z), `MX` (X); prefix a name with `!` to make that operation noiseless.
Measurements are numbered in order of appearance. `detectors: null` means
"discover them": every outcome parity that is deterministic in the noiseless
circuit. A circuit that declares observables must declare its detectors.

Library gadgets (`library.py`): `cat:N[,verify]`, `css_prep:CODE` (CODE a
qeccert spec), `surface:D,R[,z|x[,bad]]`.

Library use: `from faultcert import load, build, verify, Model`;
`Model(circuit).replay([(location, "XZ"), ...])` gives the (detectors,
observables, syndrome) signature of any fault set.

Exit codes: 0 ok; 1 refusal (receipt fails, legs disagree, Stim disagrees);
2 usage or malformed circuit.

## GATES

`python3 tools/faultcert/selftest.py` (about 2 s; 5126 checks without Stim,
5179 with Stim and the heavy tier):

- Stabilizer simulator on GHZ, Bell, and reset-and-reuse circuits.
- **Model against a dense state-vector simulation**: 60 random circuits on 2 to
  4 qubits, 8 random fault sets each. In every measurement branch the detector
  pattern equals the model's, the reported outcomes are possible in a noiseless
  run exactly when no detector fires, and the reduced output state equals the
  believed one exactly when the model's syndrome is zero. Coverage counters
  require at least 20 cases each of detected, harmless, harmful-undetected and
  frame-only (outcome flips alone) harm.
- Linearity: the effect of a fault set is the XOR of its members; each
  location's effects are closed under addition.
- Search against brute force on 200 random column sets, including capped runs.
- Strict fault tolerance against the definition (all location subsets, all
  Paulis) on five small gadgets.
- Known answers: the cat-4 chain has a single-fault weight-2 hook; one
  verification ancilla removes it; rotated surface d = 3 has circuit distance 3
  with the standard CX order and 2 with the hook-aligned order, in both bases
  (d = 5: 5 and 3 in the heavy tier).
- Receipt controls, each refused: a witness with a fault removed or moved, a
  lower bound above the witness or without evidence, a circuit edited under a
  stale hash, an edited census, a forged enumeration claim (`--full`), a
  strict-FT counterexample that is not one; malformed circuits; a declared
  detector that is random; observables without declared detectors.
- SAT leg: surface d = 3 distances with checked proofs, equal to leg A.
- Stim leg: identical mechanism sets on eight library circuits and 40 random
  circuits; a model with mechanisms removed, a wrong detector and an
  over-claimed bound are each reported as disagreements.

Four deliberate defects were planted in `model.py` during development (frame
term dropped, CX Z-propagation dropped, H not swapping, MX flipping on the wrong
component); the dense test failed on each (41 to 142 failing checks).

## FOOTGUNS

- **Detectors decide everything.** With discovered detectors, every
  deterministic parity counts as a check, including a logical value that happens
  to be fixed by the preparation. That is right for preparation-and-verify
  gadgets and wrong for memory experiments, which is why observables require
  declared detectors.
- **Weight counts locations, not qubits.** A two-qubit Pauli after one CX is
  one fault. A fault on an idle is a fault; `library.with_idles` inserts idles,
  a hand-written circuit has only the idles you write.
- **"Strict" is Pauli weight on the output**, not per-sector weight. A CSS
  argument that counts X and Z dirt separately is weaker than this property;
  a gadget can fail strict FT at order j and still satisfy the sector version.
- **Wrong Pauli frame is counted as harm.** A flipped outcome that an
  experimenter would use to fix the frame of the output makes the output state
  wrong even though no physical error is on it.
- Discovered detectors need the symbolic simulation of the whole circuit, whose
  cost grows with qubits squared per measurement: 368 qubits and 3451 locations
  (a 16-copy Golay purification tree) take 12 s, a 64-copy tree is out of reach.
- `lower` is what was proven. A capped leg A run says so in its evidence note.
- Leg S compares fault mechanisms, not probabilities, and Stim's logical-error
  search is heuristic: it gives an upper bound only.
- Preparations assume a fresh qubit; mid-circuit `R`/`RX` discards whatever was
  there, including errors.

## What runs publicly

Everything except the two engine-dependent legs, which run wherever `cadical`
and `stim` are installed.

CREDIT: symbolic stabilizer states with sign forms over measurement outcomes,
and the formulation of gadget fault tolerance as a property over fault counts,
follow K. Chen, Y. Liu, W. Fang, J. Paykin, X.-C. Wu, A. Schmitz, S. Zdancewic,
G. Li, "Verifying Fault-Tolerance of Quantum Error Correction Codes" (CAV 2025,
arXiv:2501.14380). The view of a circuit under noise as the set of its
undetectable faults is that of B. Rodatz, B. Poor, A. Kissinger, "Fault
Tolerance by Construction" (arXiv:2506.17181); their fault equivalence is not
implemented here. Tableau updates follow S. Aaronson and D. Gottesman (2004).
Detector error models and the cross-check leg: Stim (C. Gidney, Quantum 5, 497,
2021). Strict fault tolerance of preparations: P. Aliferis, D. Gottesman,
J. Preskill (2006). SAT leg: qeccert, CaDiCaL (A. Biere et al.), LRAT
(Cruz-Filipe, Heule, Hunt, Kaufmann, Schneider-Kamp 2017).
