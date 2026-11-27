# ftqc_gadgets — GUIDE

KIND: package (Python, stdlib only; imports `ftqc_audit` and `ftqc_codes` from this
repository; `gadget.py`, `effects.py`, `props.py`, `library.py`, `certify.py`,
`selftest.py`). Battery class: selftest (everything runs from the clone).

## PURPOSE

Fault properties of Clifford gadgets that act on unknown inputs, and comparison of
two gadgets fault by fault. `ftqc_audit` handles circuits that start from prepared
qubits; this package adds what a fault-tolerance argument is actually made of:
gates and error correction applied to an arbitrary encoded state, decoders, the
conditions under which a gadget may be replaced by a simpler description of it.

How unknown inputs become a finite check. Each input block enters as one half of a
noiseless encoded Bell pair; the other half is a set of reference qubits that
nothing touches. The joint state of references and outputs then determines the
gadget's action on every input state, and it is a stabilizer state, so "the gadget
did the right thing" is "no stabilizer of that state is violated". The errors an
input arrives with are faults on one extra idle per input qubit, counted as input
weight and never as gadget faults.

Decoders. A decoder is a lookup table from check values to a lightest error, applied
either to the measurement record (transversal measurement, then classical decoding,
as in teleportation-based error correction) or as a Pauli frame on output qubits
(syndrome extraction). Everything before the lookup is linear over GF(2), so fault
sets are still sums of columns; the lookup is applied to each sum.

Outputs (as receipts, every failure with a replayable input error and fault set):

- **check prep | gate | ec, t**: for all input weight r and s faults with
  r + s <= t, every accepted run leaves each output block within weight `bound` of
  the ideal state, logical action included. `prep`: bound s. `gate`: bound r + s.
  `ec`: bound s (what came in is gone).
- **close s**: for *any* input, j <= s faults leave each output block within weight
  j of its code space. `ec` and `close` together are the two error-correction
  properties of Aliferis, Gottesman and Preskill.
- **refine IMPL SPEC, order w**: every accepted fault set of j faults in IMPL (with
  any input error of weight r' <= r, r' + j <= w) has exactly the effect of an
  accepted set of at most j faults in SPEC with the same input error. Two gadgets
  that refine each other are fault-equivalent through order w. A specification is
  an ordinary gadget, usually the noiseless gadget followed by one noisy wire per
  output qubit (`,ideal`, `prep:...,spec`, `library.idealize`).
- **compose**: a gadget analysed with preparations given as specifications, plus a
  refinement receipt for each real preparation, gives a composition receipt: the
  gadget with the real preparations has the property through
  min(order of the property, orders of the refinements). The statement and its
  conditions are printed in the receipt (`certify.THEOREM`). An implementation may
  itself be a composed gadget (pass `dict(composition=, target=, parts=)` as the
  part), which is how a 64-copy tree is reached from 16-copy trees.

- **local** (`--local`): the coordinate-locality certificate. Every fault and every
  input error is shown to act within one coordinate (position) of each block, by a
  check that is linear in the size of the gadget. A lemma printed in the receipt
  (`props.LEMMA`) then gives `gate`, `ec` and `close` whenever faults and input
  errors together touch at most T coordinates, T the decoders' radius, however many
  faults lie inside those coordinates. This is the route for codes whose t is too
  large to enumerate (QR-47, BCH-127): a candidate radius can come from the code's
  `ftqc_codes` receipt (`--code-receipt`, `library.radii`) and decoder tables are
  skipped (`--no-tables`). For a lightest-pattern decoder the verifier then
  proves from the actual check matrix that its kernel contains no nonzero word
  of weight at most twice the requested radius. A code receipt or an explicit
  number supplies the target. An insufficient search bound yields a refusal.
  A hierarchical decoder instead uses its verified inner unique-decoding radius.
  Wire order fixes the physical coordinate for every fault mode. All gates and
  classical controls in the body must respect those coordinates. The exception
  for a fresh specified ancilla is checked against its complete preparation
  circuit and its noisy boundary, before the ancilla couples to live data.
- **levels** (`levels.py`): the outer level of a concatenated scheme. A qubit of the
  outer gadget stands for an inner block and a fault for a failed inner gadget;
  because the locality certificate counts coordinates, it already says "at most T2
  bad columns, hence correct". `levels.reduction` joins that certificate to one
  inner receipt per kind of operation the outer circuit uses and records what is
  quoted (the level-reduction step) and what is assumed (outer ancillas given as
  specifications). `levels.two_level` builds the same gadget on the concatenated
  code with hierarchical decoders (`ftqc_codes.hier`), and `levels.Direct` tests the
  statement on it: fault sets with at most t2 columns holding more than t1 faults
  must act ideally, exhaustively at low weight and by sampling above, with a search
  for the first failure one bad column further.
- **tree** (`--tree`, `factories.py`): the factory certificate of a purification tree,
  for the statement a level of concatenation needs and enumeration cannot reach: any
  number of faults that stay on their own coordinate are harmless, and only faults
  that can spread are counted. The circuit is recognised operation by operation as
  raw preparations (leaves) joined by checks (transversal CX, the partner measured,
  acceptance on the fixed parities of the state); every leaf location is classed per
  sector as harmless, local at one coordinate or junk-capable; the detectors of each
  check must span exactly the state's parities; the minimum weight of the patterns
  that act trivially is proved. The filter rule printed in the receipt then gives the
  contract: the output is its specification with wire faults only on coordinates of
  faulty local and check locations, unless K leaves hold a junk-capable fault, with
  K reported per sector. `--strict` keeps only the last idle round of a leaf local
  (`prep:CODE,STATE,tree=SEQ,rounds=K`): one level up an idle is an inner error
  correction, and a raw outer preparation is not transversal. `--light` skips the
  propagation of every fault through the whole circuit, which the certificate does
  not need. `factories.test` compares the rule with propagation set by set.
  `levels.reduction(..., factories={role: receipt})` replaces an assumed outer ancilla
  by its factory and records K: the runs with K or more junk leaves are the junk
  event, which the receipt leaves to a probability bound. `Tree.bounds` is such a
  bound, under one hypothesis (faults of different leaves and of different checks
  are independent): from per-leaf and per-check failure probabilities it gives, by
  the recursion printed as `factories.RECURSION`, the junk probability given
  acceptance, a lower bound on each check's acceptance and the expected number of
  leaves consumed. It is tested against sampling and against the formulas written
  out; it credits a check with accepting only when both inputs are clean, so it is
  loose, and for wide blocks with frequent leaf failures it gives nothing.
  `Tree.profile()` (and `factories.profile(receipt)`, which verifies the receipt first) gives
  the tree as plain data for a probability calculus: junk-capable and local locations of
  every kind of leaf, and the gate, measurement and idle locations of every check; the
  `filter-contract` study of `ftqc_compose` reads it. `factories.accepted_law(g, p)` sums
  exactly over all fault sets and returns the law of the accepted output's error class under
  independent faults; it is feasible when the state has few classes (Steane) and is what
  that calculus is tested against.
- **cones** (`cones` command, `cones.py`): fault cones of a whole program. A program
  is a list of layers of transversal operations on code blocks, with error
  correction placed wherever the schedule puts it (`ECa`, `ECb`), classically
  controlled gates, decoded measurements and delivered resources. For every decoding
  event the receipt gives its cone: the location classes and ancilla atoms whose
  faults can be part of the pattern its decoder sees, with the weight per coordinate
  that a probability bound needs. With an outer level it gives the outer cones too,
  whose sources are the inner events. Three routes are compared. `Flow` applies the
  propagation rules (unions along gates, a correction as a sink). `exact` re-derives
  the cones by linear algebra over GF(2), modelling a correction as the measurement
  and feedback it is, for several settings of the classical controls, and must stay
  inside the rule's cones. `Direct` builds the program qubit by qubit on a code and
  runs real lookup decoders with corrections carried through later gates: the first
  decoding that goes wrong must have more than t dirty coordinates in its cone.
  `cones.ledger(flow)` reports what a schedule spends: block-ticks, the largest number of live
  blocks, inner and outer corrections and resources delivered by kind; `cones.count(program)`
  gives the same counts from the layers alone, for sweeps over many schedules.
  `cones.correction_pattern(program)` measures how gates and idle ticks are followed by
  corrections (runs of steps without one). The `volume` study of `ftqc_compose` prices the
  ledger.
- **independent_of_input**: which decoders' check values no input error can change.
  For repeated error correction (`steane_ec:CODE,rounds=2`) the later rounds have
  this property and the first does not: after one round, what a block does next is
  driven by that round's faults and ancillas, not by what the block arrived with.

Weight is Pauli weight per block by default; `--sector` uses the larger of X weight
and Z weight (the convention of arguments that treat the two sectors of a CSS code
separately). The two give different verdicts and the receipt says which was used.

## USE-WHEN

- A proof uses "this gadget satisfies the gate property" or "the error-correction
  properties" for a concrete circuit and you want that checked on the circuit, for
  all inputs, at the orders the proof needs.
- A large gadget contains sub-circuits (ancilla factories) too big to analyse in
  place: prove each one refines its specification, analyse the gadget with the
  specifications, compose.
- Two circuits are claimed interchangeable under noise (a rescheduled or optimised
  gadget against the reference one): `equiv`.
- A decoder is part of the claim (a hook that a lookup decoder turns into a
  weight-2 error): `seq_syndrome:steane` is the worked negative example.

## NOT-FOR

- All-orders statements and probabilities. Everything is "for all fault sets of at
  most j locations", j small. Composition transfers such statements; it does not
  produce a threshold.
- Non-Clifford gates and adaptive circuits. Decoders may only feed a Pauli frame or
  correct the record; they cannot choose later gates.
- Output weights beyond the radius in which the lightest error of a block is unique
  (above t for a distance 2t + 1 code): the package refuses and says so.
- Fault equivalence in the full sense of Rodatz, Poor and Kissinger (all weights, on
  ZX diagrams with their rewrite rules). This is the truncated circuit-level
  relation, decided by enumeration.
- Large gadgets at high order: the cost is the number of accepted fault sets. A
  transversal gate or an error-correction step on the Golay code at r + s <= 3 takes
  one to two minutes; a 16-copy purification tree (3451 locations) against its
  specification at order 3 takes about a minute. Sector weights on two Golay blocks
  at input weight 3 are 2.7e8 input cases and are refused by the budget.

## INVOKE

With `tools/` on `PYTHONPATH` (or from `tools/`):

    python3 -m ftqc_gadgets analyze cnot:golay --check gate:3 -o r.json
    python3 -m ftqc_gadgets analyze knill:steane,anc=spec --check ec:1 --close 1 -o ec.json
    python3 -m ftqc_gadgets analyze prep:steane,zero,tree=XZ --refines prep:steane,zero,spec:1 -o zero.json
    python3 -m ftqc_gadgets analyze prep:steane,plus,tree=ZX --refines prep:steane,plus,spec:1 -o plus.json
    python3 -m ftqc_gadgets compose ec.json zero=zero.json plus=plus.json -o composed.json
    python3 -m ftqc_gadgets equiv cnot:steane cnot:steane,ideal --order 2 --input-weight 1
    python3 -m ftqc_gadgets verify r.json --full
    python3 -m ftqc_codes certify qr:47 -o qr47.json
    python3 -m ftqc_gadgets analyze steane_ec:qr:47,anc=spec --local --code-receipt qr47.json --no-tables
    python3 -m ftqc_gadgets analyze prep:golay,zero,tree=XZXZ --tree -o tree.json
    python3 -m ftqc_gadgets analyze prep:qr:47,zero,tree=XZXZ,rounds=2 --tree --strict --light --wmin 11
    python3 -m ftqc_gadgets cones program.json -o cones.json     # {"layers": [[["CX","A","B"]], ...], "inputs": [...], "outputs": [...]}
    python3 -m ftqc_gadgets verify cones.json

Library gadgets (`library.py`): `prep:CODE,STATE[,raw|spec|tree=SEQ][,rounds=K]`,
`cnot:CODE[,leaky]`, `idle:CODE`, `knill:CODE[,anc=METHOD]`,
`steane_ec:CODE[,anc=METHOD][,rounds=K]`, `measure:CODE[,x|z]`, `seq_syndrome:CODE`,
`fc:ftqc_audit_SPEC`; suffixes
`,ideal` and `,idles`. CODE is a ftqc_codes spec without commas, or `rep3`.

Your own gadget, in Python:

    from ftqc_gadgets import Builder, Analysis, check
    b = Builder("my_gadget")
    A = b.input(code)                       # block with its references and input-error idles
    anc = b.alloc(code.n)
    b.tick([("R", q) for q in anc])
    b.tick([("CX", a, q) for a, q in zip(A, anc)])
    recs = b.measure(anc, "Z")
    b.decoder("x_errors", recs, code.hz, "X", A)     # or b.checks(recs, masks) to postselect
    b.output(A, code)
    print(check(Analysis(b.build()), "ec", 1))

Exit codes: 0 ok; 1 a property fails or a receipt is refused; 2 usage or malformed
gadget.

## Fault cone receipts

Run the [two-level correction schedule](examples/two-level-cones.json) from the BootLoops root.

```sh
export PYTHONPATH="$PWD/tools"
out=$(mktemp -d)
python3 -m ftqc_gadgets cones tools/ftqc_gadgets/examples/two-level-cones.json -o "$out/cones.json"
python3 -m ftqc_gadgets verify "$out/cones.json"
```

Read `inner_events` for per-coordinate location weights and ancilla contributions, and `outer_events` for the inner events and ancilla errors that can reach each outer decoder. The receipt embeds the schedule and its digest. For a program with classical controls, `--samples N` selects N seeded random settings in addition to all on and all off for the inner algebra check. The one-level physical decoder checks run in the gadget self-test. [Interpret the example](../FTQC.md#trace-a-fault-through-a-correction-schedule)

## GATES

`python3 tools/ftqc_gadgets/selftest.py` (about 35 s, 6943 checks; `--heavy` adds the
Golay gadgets at t = 3 and larger two-level tests, about 5 minutes):

- **Against a state-vector simulation**: seven small gadgets (teleportation and
  syndrome-extraction error correction on the 3-qubit code, a transversal CX,
  sequential syndrome bits, two purification steps, a bare input encoded and
  decoded), random faults and input errors, every measurement branch. The detector
  pattern and the decoder check values equal the model's, and the state of
  references and outputs equals the believed one (noiseless run forced to the
  corrected record, then the Pauli frame) exactly when the effect is zero.
- **Against the definition**: the three properties on ten gadgets by replay over
  every location subset, every Pauli label and every input Pauli, with the bounds
  written out separately in the test.
- Refinement against brute-force effect tables; refinement of a gadget against its
  idealized form equals `ftqc_audit`'s strict fault tolerance on five circuits.
- Coset-leader tables and unique radius against brute force (Golay: 2048 syndromes,
  unique through weight 3).
- Known answers: transversal CX on Steane and Golay holds through t and is
  fault-equivalent to its specification; a CX pair inside a block breaks it at one
  fault; teleportation error correction holds with specified or purified ancillas
  and fails at one fault with raw ones; sequential syndrome extraction fails at one
  fault.
- Composition: Steane teleportation with XZ and ZX trees; the direct analysis of
  the composed circuit agrees; a two-level composition agrees with the direct
  analysis of the larger tree; a part that refines only through order 0 lowers the
  conclusion; missing, mismatched and inflated compositions are refused.
- Receipt controls, each refused: witness with the fault removed or moved,
  inflated order, a positive claim without cases, edited circuit or decoder under a
  stale hash, edited statistics, forged case count (`--full`), a refinement witness
  the specification matches.

- Locality: on eleven library gadgets the certificate's radius and properties are the
  expected ones and enumeration agrees through that radius; random sets of many
  faults inside few coordinates stay within the bound; a CX inside a block,
  sequential extraction, raw ancillas, a correction applied one position over, and a
  decoder whose positions are numbered differently from the block's are each refused
  with the offending fault or input error; without decoder tables and with radii
  from a code receipt the certificate is the same.
- Two levels: on Steane over Steane (49 qubits per block), decoded measurement,
  transversal CX and Steane-type correction built physically with hierarchical
  decoders: every set of at most one (measurement: two) faults or input errors and
  thousands of sampled sets with one bad column act ideally; with two bad columns a
  failing set is found; the outer certificate tolerates one coordinate. Reduction
  receipts: missing, mismatched, weaker and inflated hypotheses are refused or
  lower the conclusion.
- Repeated correction: second-round decoders are independent of the input, first-round
  ones are not, with a replayable input error.

- Factories: on nine Steane trees the rule agrees with propagation for every set of
  at most two locations (one on the larger trees) and thousands of sampled sets,
  among them sets of up to twenty faults on local locations, the same fault mirrored
  in neighbouring leaves, and single-coordinate errors on the two sides of a check
  that add up to a stabilizer. K is the least number of junk leaves that raises the
  flag, by search over leaves; enumeration against the specification confirms the
  order and, on the Golay code, stops exactly where the contract does. Location
  classes equal the definition by replay. Eight things that are not a tree (missing
  or foreign detector, reversed or skewed or partial gate, a gate inside a block,
  an idling partner, leaves in different states) are refused. Error correction with
  real trees: faults confined to at most T coordinates, however many, with fewer
  than K junk leaves, leave the output within that weight, and two mirrored junk
  faults defeat it. Levels receipts with factories; refusals for a loose, foreign or
  wrong-state factory and for missing inner receipts.

- Cones: known answers on small programs (what a measurement after a gate sees; a
  full correction leaves only its own sources, which is also what a standard block
  carries; a controlled gate has weight two). On forty random programs linear
  propagation never leaves a cone, and two equal gates show the rule overestimating.
  Physical tests on the Steane code (Golay with `--heavy`), controls on, off and
  random: no decoding goes wrong while every cone has at most t dirty coordinates,
  with sets of up to a dozen faults, and a failure is found with t + 1. A cone with
  one gate removed is caught. Two levels: outer correction is a sink for earlier
  inner events. Malformed programs and edited receipts are refused.

Seventeen defects were planted in the cone code (a gate not carrying a sector, a
correction that is not a sink or forgets an ancilla's errors or its measurement,
idles dropped, sectors exchanged, the algebra without its feedback or without the
ancilla, the physical test without carried corrections or without recognising a
wrong decoding, outer correction not a sink, inner events not entering the outer
level, a non-Clifford location pooling one sector, a controlled gate of weight
one). The battery fails on each; two needed a new check first.

Six more were planted in the probability recursion (sum and product exchanged,
fresh dirt or the partner's pattern dropped, acceptance not depending on the inputs,
retries not counted, a leaf's rate counting one location). Three needed exact
formula checks before the battery failed on them; sampling alone did not see them.

Twenty-one defects were planted in the factory code (the intersection taken without
the weight condition, step faults dropped on either side, the two sectors' rules
exchanged or weakened, junk read as harmless, several coordinates accepted as local,
detectors not compared with the state, the count summed in the wrong sector, the
measured block not tied to the gate, strict or physical not enforced, leaves not
compared, an idling partner allowed, four omissions in the levels receipt, the
certificate not recomputed on verification). The battery fails on each.

Ten further defects were planted in the locality and two-level code (hierarchical
decoder ignoring the inner step, an error allowed to leave its coordinate, an input
error allowed at another coordinate, the largest decoder radius taken, `ec` claimed
with input errors left, correction not applied, bad columns miscounted, a given
radius unchecked, later rounds reading only their own outcomes, the any-input
property dropped from a reduction). The battery fails on each; four needed a new
control first.

Eleven defects were planted during development (decoder correction applied as the
wrong kind, decoders ignored, syndrome left in a non-canonical basis, input X and Z
exchanged, decoder table not lightest-first, subset enumeration double counting,
gate bound without the input weight, refinement accepting a costlier match, output
frame ignored, logical remainder ignored, weight bound not enforced). The battery
fails on each. The last was a real bug found this way.

## FOOTGUNS

- **The specification decides what refinement means.** One fault of a transversal
  two-block gate can leave an error on both blocks, so its specification has one
  noisy two-qubit location per qubit pair (`cnot:CODE,ideal`), not one wire per
  qubit. Against single wires the gate does not refine, correctly.
- **Refinement is exact matching of effects, the properties are weight bounds.**
  A gadget can satisfy `gate` and not refine a given specification.
- **Pauli weight and sector weight differ.** One CX fault can leave X on one qubit
  and Z on another: weight 2 as a Pauli, weight 1 per sector.
- **Decoders with ties.** The lookup takes the first lightest pattern in position
  order. `unique_through` in the receipt's statistics says up to which weight that
  choice is forced; beyond it the verdict depends on the tie-break.
- **Detectors are postselection.** In a gadget with decoders every postselected
  parity must be declared (`Builder.checks`); deterministic parities that are
  neither declared nor fed to a decoder are simply not used.
- **Idles exist only where written**, or everywhere with `,idles`. Composition
  says nothing about the time an output of a preparation waits before it is used.
- `close` for any input is decided either by showing that the output's code
  syndrome cannot depend on the input (teleportation), or by enumerating the span
  of input errors (2^8 signatures for Steane); otherwise it reports "not decided".
- A property that stops on its budget reports `holds_through` for the orders it
  finished and a `note`. Read the note before reading a low `holds_through` as a
  failure: only a `violation` is one.
- **The locality lemma, like composition, is an argument.** It is printed in the
  receipt and tested against enumeration and against direct two-level analysis.
  What the tool computes is that its hypotheses hold for the circuit.
- **A levels receipt is not a threshold proof.** It certifies "few bad columns,
  hence correct" for the outer gadget and the inner properties level reduction
  needs. The step between them is quoted from the literature for the fixed pattern
  of error correction, and outer ancillas are listed as assumed unless their
  factories are supplied.
- **Radii from code receipts assume lightest-pattern decoding and need the lightest
  stabilizer to be heavier than the distance** (true for the quadratic-residue and
  BCH codes here, false for surface codes, where `library.radii` returns a small
  radius by design).
- **Cones are sets of sources, not probabilities.** The receipt gives each event's
  weight and atoms; the tail bound over coordinates is a separate step. A cone is an
  overestimate wherever paths cancel.
- **The cone theorem's induction is an argument.** The tool checks its hypotheses on
  the program (rule against algebra) and tests its conclusion on physical circuits
  of small programs. Non-Clifford locations and outer operations are not built
  physically; they enter by rule only.
- **A program's cones depend on its schedule.** The receipt is bound to the layers by
  hash; moving one correction changes it.
- **A factory contract is not a refinement order.** It says nothing about a run in
  which K leaves hold a junk-capable fault, and it counts leaves, not columns: one
  bad column crosses every leaf. In the levels receipt the columns of failed local
  and check locations count toward the radius; the junk event is separate.
- **Which sector a check protects.** A check of one sector copies the partner's
  errors of the other sector onto the survivor. A candidate verified against fresh
  unverified copies therefore has K = 1 in the sector whose verifiers come first:
  one fault in a later verifier's preparation is accepted. Trees avoid this because
  the partner of every check is itself a verified subtree.
- **Minimum weights and radii of large cyclic codes** rest on the BCH bound when the
  search stops short (zeros computed exactly from the matrix at hand, theorem cited);
  the receipt says which. In a non-cyclic qubit order there is no such fallback.
- Without decoder tables (`--no-tables`) only the locality certificate is available.
- The composition theorem is an argument printed in the receipt and tested against
  direct analysis on small gadgets. It is not machine-checked.

## What runs publicly

Everything.

CREDIT: the properties are those of P. Aliferis, D. Gottesman, J. Preskill,
"Quantum accuracy threshold for concatenated distance-3 codes" (2006). Fault
equivalence and the idea of proving a circuit fault tolerant by relating it to an
idealised specification: B. Rodatz, B. Poor, A. Kissinger, "Fault Tolerance by
Construction" (arXiv:2506.17181). Symbolic stabilizer states and fault tolerance as
a property over fault counts: K. Chen et al., "Verifying Fault-Tolerance of Quantum
Error Correction Codes" (CAV 2025, arXiv:2501.14380). Teleportation-based error
correction: E. Knill (2005). Syndrome extraction with encoded ancillas: A. Steane
(1997). Fault model and simulator: ftqc_audit. Codes: ftqc_codes.

New receipts use `ftqc_gadgets-receipt-1` and `ftqc_gadgets-composition-1`. The verifier also reads archived `gadgetcert-receipt-1` and `gadgetcert-composition-1` receipts.

Two-level reduction binds each inner operation to its signed logical action on
canonical encoded Bell inputs. Input and output order identifies the CNOT
control and target. The receipt gives the measurement parity and any tracked
Pauli-frame parities as raw record indices. Idles, Hadamards, preparations and
destructive X/Z measurements have distinct roles. Unimplemented operations are
refused even if marked ideal in the outer body. Input-error placeholders and
matched fresh ancilla specifications have their own boundary contracts.

## Signed surgery instruments

See [the surgery interface and replay examples](../ftqc_gadgets/SURGERY.md) for physical measurement receipts, installed feed-forward schedules and conditional residual laws.

### Magic transfer certificates

Run `python3 -m ftqc_gadgets magic-transfer spec.json -o receipt.json` and
`python3 -m ftqc_gadgets verify receipt.json`. A `fresh` input specifies
`n`, CSS check masks `hx` and `hz`, and a port with the block name, ordered
logical X/Z masks and frame/correction. Its equal-amplitude cosets and erased
input wire certify the open isometry. A `swap` input specifies a block with
signed checks, a complete ordered logical basis, and the destination slot.
Checking the Pauli images certifies transfer into that slot and identity on
every spectator, including its reference. The physical circuit and its wire
order are part of the receipt.
