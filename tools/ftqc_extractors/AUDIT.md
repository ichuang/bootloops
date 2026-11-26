# Comparing ftqc_extractors and ftqc_audit

Use ftqc_extractors to prove which logical measurement or channel an HGP construction implements and to establish its distance under the stated boundaries. Use ftqc_audit to ask which faults in a concrete Clifford gadget escape its flags and how much physical error they can leave on its outputs. The bridge lets both programs answer the same logical-distance question on supported schedules and lets ftqc_audit inspect the exact ancilla ftqc_extractors will consume.

| Question | ftqc_extractors | ftqc_audit |
|---|---|---|
| What is supplied at the input? | HGP seeds, logical ports and an encoded input that may be entangled with a reference. The conditional T construction also assumes an encoded magic state. | A tick circuit with explicit preparations. A symbolic stabilizer simulation tracks the ideal state and its dependence on measurement records. |
| What counts as detection? | Declared parity or flag checks and the ideal terminal HGP syndrome. | Declared measurement-record checks, or deterministic parities discovered in a preparation. An output-state stabilizer violation can count as harm. |
| What constitutes failure? | A silent change to the requested readout, protected output or a spectator logical. Static surgery and the canonical extractor have separate mathematical distance statements. | Circuit distance uses a flipped declared deterministic observable. Fault distance also counts an incorrect prepared state. Strict FT through order s requires every accepted j-fault set, for j ≤ s, to leave an error equivalent on the output to a Pauli of weight ≤ j. |
| What evidence can I replay? | Finite algebraic checks, signed ideal flows, explicit fault witnesses, bounded searches, structural bounds and supported SAT receipts with embedded LRAT. | Native witnesses, exact meet-in-the-middle enumeration, checked SAT with proof digests, and optional Stim detector/observable comparisons. Native full verification reconstructs the model and checks finite bounds by enumeration. |
| How is a logical T handled? | The branch identities, magic input and conditional frame establish the ideal logical action. Each compiled Clifford arm then has a Pauli fault calculation. ftqc_extractors also has a separate source-inclusive model. | The Clifford transport can check those arm calculations. A native stabilizer preparation cannot certify a magic state or the adaptive non-Clifford channel. |

Both count distinct faulty locations. An arbitrary correlated Pauli on the two CX operands costs one fault. Measurement faults flip the reported bit. For a preparation, ftqc_extractors represents the eigenstate flip once. Native ftqc_audit includes all three one-qubit Paulis. The extra stabilizing Pauli has zero effect on the shared tests and the remaining Pauli is equivalent to the flip. Both use the sibling ftqc_codes SAT/LRAT implementation. Its checker is a shared dependency.

## Cross-check a logical schedule

Run these commands from the BootLoops root with Python 3.

```sh
export PYTHONPATH=tools
python3 -m ftqc_extractors.audit_bridge compare tools/ftqc_extractors/examples/edge-extractor.json
python3 -m ftqc_extractors.audit_bridge compare tools/ftqc_extractors/examples/t-gate.json --arm Y
python3 -m ftqc_extractors.audit_bridge replay tools/ftqc_extractors/examples/edge-hook.json --fault 19:8
```

`compare` constructs the ftqc_extractors schedule, checks its ideal circuit conditions, and retains its locations, acceptance checks and logical tests. For every elementary fault mode it compares ftqc_extractors's adjoint signature with ftqc_audit's forward propagation of the Pauli and flipped records. ftqc_audit's meet-in-the-middle search then computes a distance interval on the common signatures. A finite upper endpoint includes a physical witness replayed in both programs. A matching lower and upper endpoint establishes exact distance for that specification. An exhausted search budget leaves the upper endpoint open. `--max-weight` and `--cap` limit the subset search. Constructing all columns can be expensive for a long schedule.

`replay` checks a named witness directly. Zero-based location numbers and numeric modes have the same meaning as in an ftqc_extractors receipt. The report lists fired detector indices and changed logical-test indices in ftqc_extractors order. The schedule is serial and includes exactly its listed idle locations. ftqc_audit's forward calculation acts on error differences, so it works on an arbitrary encoded input. An initial noiseless tick declares the data wires for its parser and contributes zero fault locations. Using that tick as a physical initialization in a native state-preparation analysis would change the problem.

For a T specification, `--arm gate` ends at the gate output. `--arm X` and `--arm Y` continue through the corresponding signed readout. The ftqc_extractors acceptance and harm predicates are shared inputs to the comparison. Agreement tests the physical transport and search under those predicates. The branch identity, resource assumption, spectator predicate and any structural lower bound still require their ftqc_extractors checks. Use `ftqc_extractors verify` for the original logical receipt.

## Analyze an ancilla before it touches data

```sh
export PYTHONPATH=tools
python3 -m ftqc_extractors.audit_bridge export-cat tools/ftqc_extractors/examples/t-gate.json --round 0 -o /tmp/ftqc_extractors-cat.json
python3 -m ftqc_audit analyze /tmp/ftqc_extractors-cat.json --legs A,B,S --strict 2 --heavier 1 -o /tmp/cat-receipt.json
python3 -m ftqc_audit verify /tmp/cat-receipt.json --full
```

The export stops after the selected round's cat preparation and verification, before coupling to the HGP data. Each flag is an explicit acceptance detector. Cat qubits are the outputs. Add `ftqc_origin.first_location` to a local fault index to recover its location in the full ftqc_extractors schedule. The origin also binds the input digest, selected arm, round and qubit offset. Stim and CaDiCaL are optional for generation. The native full replay uses Python.

For the four-qubit cat above, a Z fault just after the initial RX flips the GHZ phase and passes all three flags. Native fault distance is one. Every accepted one- or two-fault set leaves residual weight at most its fault count. The complete eight-qubit-block T protocol has logical distance two in all three cases. The physical cat error and the logical gate error are measured at different boundaries. In this example the same Z fault flips the first joint readout and fires detector 6 when the next round disagrees. Conversely, an accepted heavy cat error is a candidate for a hook. Replay it through the consumer to establish logical harm.

See the [ftqc_extractors model and examples](GUIDE.md), [logical T proof and assumptions](T_PROOFS.md) and [native ftqc_audit guide](../ftqc_audit/GUIDE.md).
