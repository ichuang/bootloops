# ftqc_compose

Build and replay finite bad-support calculations from He, Nguyen and Pattison, §2. A support is bad when it contains a listed generator. The receipt counts inclusion-minimal generators by size. This polynomial gives a union bound when the event family covers physical failures and the supplied noise law obeys the local-stochastic subset condition.

With `tools` on PYTHONPATH, run `python3 -m ftqc_compose check SPEC -o RECEIPT` and `python3 -m ftqc_compose verify SPEC RECEIPT`. The `family` input takes an `expression` whose operations are `family` (sites and minimal supports), `union`, `product` (left and right), or `substitute` (outer expression and an expression for each outer site under `at`). Products and substituted children require disjoint physical domains. Every generated family receipt checks every subset against an independent Boolean interpreter, up to 16 physical sites. Larger domains are refused. [Paper](https://arxiv.org/abs/2508.08246)

## Typed contracts

A `contract` study supplies an ordered tree under `program` using `sequence` and `parallel`. The `accepted-diagonal-pauli` scope tracks correctness of accepted branches under stochastic Pauli faults. It permits correlated fault supports. Every join matches the ordered wire coordinates, classical or quantum kind, code, logical basis and decoder. Avoidance of the output's bad residual supports must imply avoidance of the next input's bad supports. Sequential and parallel composition take the union of the physical bad families on disjoint location domains.

The built-in classical leaves are three fallible preparations of zero (`prepare0`), an ideal classical `majority` decoder and an ideal `bit-wire`. A `hypothesis` leaf can name other simulation obligations and one factor of a friendly filter per input bundle, or null if friendliness is unproved. These declarations yield a conditional receipt with its unresolved premises. Set `require_closed` to require checked leaves, and `require_friendly` to require a filter derivation. Parallel composition needs both filters. Sequential composition inherits the first gadget's filter under a valid join (Prop. 4.9). A friendly conditional receipt depends on its declared hypotheses. Fully coherent faults and unconditional postselection probabilities require separate proofs.

## One-edge HGP contract

Use `{"kind":"hgp-repeat","variant":"edge-extractor"}` to check two nondestructive measurements of logical Z1 on the same eight-qubit memory. Use `edge-hook` to replay the negative control. The underlying `hgp-z` leaf takes an instance name and a variant. It can also appear in an ordered contract tree. Its output bundles are the updated quantum memory and one classical measurement bit. A second measurement in parallel with a classical wire preserves the first reported bit.

The promise assumes an encoded input with no physical error, the schedule's noisy gates and classical feedback operation, ideal syndrome tests at each output boundary, and ideal storage of the reported bits. The check conditions on acceptance. It covers diagonal Pauli alternatives at the listed locations. The bundle filter for arbitrary damaged quantum inputs is open, so requesting friendliness is refused. Basis masks (X1,Z1)=(3,5) and (X2,Z2)=(80,48) use qubit zero as the low bit. The signed Heisenberg flows check the ideal instrument on arbitrary encoded states. After projection onto a measurement branch, a stabilizer or the measured Z contributes only a phase. The output test reduces the data Pauli modulo these branch stabilizers and checks the reported bit and spectator. This certifies a clean conditional output at each successful join.

Each good leaf rejects 121 single-location modes and accepts 12 with the correct readout and residual. Taking every two-location subset within either leaf gives 342 sufficient minimal generators for the composed cover, with enumerator 342 x². Some pairs in that cover are harmless. A separate flat simulator rebases the second measurement's record controls and compares both boundary states with independent audit transport for every zero-, one- and two-location assignment. Its witness establishes distance two. In the hook circuit a Z fault on the target of the first extra CNOT (location 19, mode 8) escapes all checks and flips the spectator X test in both measurements.

A checked counterexample makes `check` exit 1 and still writes its receipt. A successful `verify` exits 0 for a faithfully replayed proof or counterexample. Invalid inputs, failed joins and receipt mismatches exit 2. A conditional receipt lists the hypotheses on which its conclusion depends.

## Define and replay a contract

This closed classical example recovers a zero bit after up to one faulty preparation. The three physical location IDs must be unique in the enclosing circuit.

```json
{
  "kind": "contract",
  "scope": "accepted-diagonal-pauli",
  "require_closed": true,
  "require_friendly": true,
  "program": {
    "op": "sequence",
    "left": {"op": "prepare0", "sites": ["p0", "p1", "p2"]},
    "right": {"op": "majority"}
  }
}
```

An external hypothesis has fields `op` (hypothesis), `id`, `ideal`, `inputs`, `outputs`, `faults` (a family leaf) and `filters`. Each ordered port supplies `kind` (Q or C), `wires` (ordered local coordinates), `code`, `basis`, `decoder`, and `bad` (minimal forbidden residual supports in these coordinates). For a nontrivial code those identifiers refer to the mathematical data of the premise. Checking that two identifiers match does not discharge the premise. A closed proof uses checked leaves only.

The canonical battery is registered as `selftest`. It requires Python 3.10+ and the sibling packages in this checkout. It covers empty and overlapping families, physical aliasing, all 512 two-level majority inputs, incompatible joins, conditional filters, tampered evidence, and the HGP proof and hook control. Replaying a receipt rebuilds the derivation against the input SHA-256. A scope or resource refusal produces no new receipt.

Credit: Zhiyang He, Quynh T. Nguyen and Christopher A. Pattison, *Composable Quantum Fault-Tolerance*, §§2 and 4. The support operations and typed composition rules follow their framework. The classical fixtures and finite diagonal-Pauli HGP checks are BootLoops implementations. The HGP ideal-flow checks come from ftqc_extractors and its forward transporter from ftqc_audit.

## Level reduction

A `level` study constructs the physical circuit by expanding an ordered DAG of gadgets. Name its external input bundles, output wires, physical sites and vertices. Each vertex binds a unique ID to a gadget and ordered input/output wire names. Quantum and classical wires are consumed once. All vertex location domains, including ideal sites, must partition the stated physical universe. The `outer` family is over all contracted vertex IDs. Each of its generators is replaced with the disjoint product of the local bad families. The result records the minimal physical support cover.

Theorem 4.11 needs trivial external codes and output-adjacent gadgets with empty bad families. The implementation supports individual unencoded bit and qubit bundles at the boundaries. Native prepare0, majority, vote3 and bit-wire discharge unconditional local obligations. For another gadget, an `assume-4.8` node explicitly wraps a hypothesis for unconditional simulation with environment, non-diagonal Pauli terms and bundle-wise friendly filters. The result remains conditional. Accepted-run contracts cannot serve as this premise.

The theorem gives an avoiding logical fault, potentially a general operation with environment. To conclude correctness one must also establish that every such logical fault preserves the intended output. Set `output_check` to `classical` for exhaustive logical truth-table verification on the native classical network, or `none` for the simulation conclusion alone. Use `require_closed` to refuse declared premises. The classical checker limits external bits plus vertices to 16; the explicit family lift limits generated terms to 200,000.
