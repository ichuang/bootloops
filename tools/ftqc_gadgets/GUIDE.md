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

Library gadgets (`library.py`): `prep:CODE,STATE[,raw|spec|tree=SEQ]`,
`cnot:CODE[,leaky]`, `idle:CODE`, `knill:CODE[,anc=METHOD]`,
`steane_ec:CODE[,anc=METHOD]`, `seq_syndrome:CODE`, `fc:ftqc_audit_SPEC`; suffixes
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

## GATES

`python3 tools/ftqc_gadgets/selftest.py` (about 12 s, 6691 checks; `--heavy` adds the
Golay gadgets at t = 3, about 2 minutes):

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
