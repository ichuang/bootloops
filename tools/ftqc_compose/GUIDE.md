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

## Local-stochastic bounds

A `local-stochastic` study embeds a checked `source` study and an exact rate `p`, for example `"0.01"` or `"1/100"`. The checker regenerates the support cover and evaluates its nonnegative weight enumerator using rational arithmetic. It reports the raw sum and the probability bound capped at one. The premise bounds the joint inclusion probability for **every** subset of locations. It allows correlations. Per-site marginals alone are insufficient.

Optionally provide a finite `law`, a list of objects with fields `support` and `mass`. The checker verifies normalization and every subset inequality on up to 16 physical sites, then calculates the actual covered-event probability. Without a law, the distribution condition is an explicit assumption. A family source gives an event bound. A level source with a checked logical output proof also bounds the unconditional output distance. An accepted-run contract bounds the joint probability of acceptance and harm under its input promise. A conditional error rate on accepted runs would also require an acceptance lower bound.

## Correcting 15-to-1 distillation

Use `output_check: "distill15"` on a level-reduction circuit of 15 distinct `magic-source` preparations and one `rm15-ideal` processing vertex. The latter denotes a fault-free encoded-plus preparation, 15 coordinate-wise T injections, syndrome extraction with minimum-weight correction, and logical decoding. The output is unconditional. The sources are the fallible physical locations. An explicit support law must assign zero mass to faults on ideal sites.

The checker constructs the punctured Reed–Muller CSS code on columns 1 through 15, with the four coordinate functions as X checks. The odd coset encodes logical zero and the even coset logical one. In this basis transversal T gives phases ω⁷ and 1, which differ by the required logical T phase (ω = exp(iπ/4)). It enumerates every weight-one syndrome, checks syndrome correction for left/right Pauli pairs and checks the corrected logical Kraus amplitudes for injection on every coordinate with integer arithmetic in the eighth roots of unity. The equality of the two logical amplitudes proves that a single faulty resource changes only the scalar and syndrome on the corrected branches. Linearity covers a general one-resource error.

The logical cover must include every pair of sources and a fault at the distillation vertex, or a stronger cover. Substitution removes the ideal vertex's singleton. The supplied study gives 105 physical pairs and the bound min(1, 105p²) on output trace distance. It also checks that two distinct resource Z faults fool the decoder into a third Z correction and leave logical Z. In the banked correlated example at p = 1/100 each pair has mass 1/10000, with the remaining mass on no faults. These disjoint atoms attain the bound if the listed sources suffer Z errors, since Z|T⟩ is orthogonal to |T⟩. Independent noise generally gives a different result.

[He, Nguyen and Pattison, Theorem 4.11, Proposition 5.2 and Proposition 7.9](https://arxiv.org/abs/2508.08246)

## Gross-code residual data and frame

A `gross-residual` input embeds a native `gross_circuit` specification under `extractor` and its complete `distance_receipt`. Verification replays the distance computation and binds both complementary probes by their circuit hashes. It then transports every location and every allowed Pauli or measurement flip to the 144 output data qubits and the 18 X split records. The counterfactual ideal data readouts reveal the raw X and Z displacements. Their outcomes can be individually random. The fault-free checks and logical tests must still be deterministic.

For split displacement s, port q gets a Z update equal to the parity of s on its spanning-tree path. On each fault the checker evaluates all terminal checks and logicals on (raw X, raw Z XOR frame Z) and compares these signs and the split cycles to the audit model. It checks that the kernel of all terminal tests is generated by the X stabilizers in one sector and the Z stabilizers plus the measured Z in the other. Multiplying by the measured Z only contributes a phase on a projected branch. All results use clean encoded input and the native stochastic Pauli location model.

History and split-cycle acceptance leave an open data block. The receipt counts its accepted single-location modes with a nontrivial residual and gives forward-replayed examples. The reported representative weights are constructive upper bounds. A zero representative is exact by row reduction. For t locations, XOR of representatives gives the stated linear spread bound even when multiple detector flips cancel. Under the additional ideal terminal syndrome test the native distance lower bound rules out an incorrect bit or a nontrivial residual below that weight. A downstream circuit that only observes the history must carry the open residual.

## From a distilled qubit to logical T

A `distilled-logical-t` study embeds a `local-stochastic` factory with the `distill15` logical output check, a native extractor `t_gate` and its complete `gate_receipt`. Set `handoff` to `fresh-product-css`. The logical input may be entangled with a reference. The factory supplies a fresh resource independent of that joint input and subsequent faults. Factory processing uses the fault-free operations in the RM15 contract. The gate uses the native flags, repeated parities, and ideal terminal syndrome checks.

The encoder takes the consumer's actual HGP checks and positive logical X and Z representatives. It prepares the uniform positive superposition over the X stabilizers, copies the raw basis bit b into the logical X string, and erases the raw wire with the logical Z parity. Binary propagation through its CNOT list verifies both cosets and disentanglement for arbitrary amplitudes. The receipt binds the ordered schedule by hash. The other resource logicals start in a fixed state. Data-block spectators may be entangled with the addressed data qubit.

For parity p and X result m the ideal consumer measures Z⊗Z and resource X, then uses the frame Sᵖ Zᵐ. The checker constructs all four Kraus operators on the joint data and resource input, verifies their completeness with exact eighth-root phases, and matches their action to the physical consumer's verified branches. With |T⟩ the action is T. With Z|T⟩ it is ZT. Encoding by an isometry preserves trace distance. Appending a fresh resource and applying this trace-preserving instrument contracts it on data plus any reference, including when the transcript is retained. Consequently half the diamond distance on logical data with fault-free encoding and consumption is at most the factory bound δ.

The two exact rates `encoder_p` and `online_p` specify separate noise assumptions. For N encoder preparations and CNOTs at per-location rate e, the probability of any encoder fault is at most Ne. A Z on the raw control after its first fanout CNOT propagates to the encoded logical Z. This single fault passes the otherwise ideal consumer, so the deterministic encoder has accepted fault distance one. The receipt gives its physical location. The encoder model includes the listed operations, with ideal discard and classical frame processing.

For a native consumer case with L locations and verified lower distance d, a joint wrong-and-accepted test requires at least d online faults with the ideal resource. Under a local-stochastic online rate q, the union bound for this event is C(L,d)qᵈ. Adding δ and Ne bounds the same test with the supplied factory and encoder. All probabilities are capped at one. A test here is a verifier that would pass a correct output with certainty. For the complete adaptive signed-X readout the checker adds the two possible arms' fault tails, counting the resource contribution once. Noise bounds apply to potential locations before selection of an arm. A rate conditioned on acceptance would additionally need an acceptance lower bound.

## Two-block reading residuals

A `block-reading-residual` input embeds the `gross_block_reading` specification as `extractor`, its replayable `distance_receipt`, and a nonnegative integer `faults`. The result sums the two component spread bounds after the split frame and caps the resulting physical support by the 288 output data qubits. It applies to the open output under the declared clean-input and gate-location model. For one fault the X-reading fixture gives an upper bound of 74 and the Z-reading fixture gives 90. Both count residual support on the 288 data outputs after the frame update.

## Simultaneous Y source and ideal consumers

The [BB144 consumer example](examples/bb-injection-consumers.json) takes logical slots 2 and 3 of the four-output Y injection. In each ideal S gadget the data controls a CNOT into the resource, its Z measurement selects a Z correction on the data, and a wrong Y sign induces a Z after S. Replay checks both measurement branches on both input basis columns, which determines the linear map on an arbitrary data state. The interface assumes ideal terminal decoding and transfer. Its frame gives one detector-parity bit mask per resource in circuit detector order.

The `injection-batch` verifier rebuilds the actual location responses. Within one physical channel its Pauli alternatives are exclusive. Distinct listed channels are independent and unlisted channels are ideal. A finite convolution tracks the complete postselection syndrome and all output signs, then conditions on their shared acceptance. Set locations to the literal string `all` to request every noisy channel. Exceeding the state budget leaves that guarantee unavailable. The example enables three channels including the noisy CX on physical sites 38 and 254. Conditional on acceptance the two consumers both acquire a Z with probability 1/14988. The receipt supplies the full joint distribution. The empirical BP+OSD histograms in the sampling study are a separate evidence type.
