# Fault-tolerant quantum computing in BootLoops

The eight tools on `ftqc_integration` cover different parts of a fault-tolerance and logical-compilation argument. Start with the object you want to certify and identify what the next step needs from its receipt.

| Tool | Input and question | Evidence |
|---|---|---|
| [ftqc_codes](ftqc_codes/GUIDE.md) | Given CSS check matrices, which logical errors can the memory detect? | Code parameters, bounds on X and Z distance, witnesses and hierarchical decoding guarantees. |
| [ftqc_extractors](ftqc_extractors/GUIDE.md) | Which logical measurements can a fixed apparatus perform, and what fault distance does a surgery sequence preserve? | Fixed-coupler graphs and bridges, transversal block readings, padded hypergraph surgery, physical fault witnesses, and resource-conditional logical T receipts. |
| [ftqc_audit](ftqc_audit/GUIDE.md) | Which physical faults can evade a Clifford circuit’s detectors? | Detector and observable flips per location, replayable fault sets, circuit-distance intervals and preparation checks. |
| [ftqc_decoders](ftqc_decoders/GUIDE.md) | Does a polynomial-time decoder achieve the code or circuit distance without trapping sets, and what operational error rate does it certify? | Operational radius proofs, split-stabilizer trapping-set witnesses, two-stage metacheck repair, Peierls cluster LER bounds, Stim `.dem` export and audited external (`Tesseract`, `PyMatching`, `ldpc`) runs. |
| [ftqc_gadgets](ftqc_gadgets/GUIDE.md) | What does a decoded Clifford gadget do to an unknown encoded input with incoming errors and internal faults? | Residual-error bounds, component refinement, two-level checks on concatenated codes, and the faults that can reach each decoding event in a scheduled program. |
| [ftqc_zx](ftqc_zx/GUIDE.md) | Can a candidate replace a reference circuit under the specified noise model? | Pauli-web equivalence checks, ZX rewrite proofs, flag decoder tables and stochastic gadget contracts. |
| [ftqc_compose](ftqc_compose/GUIDE.md) | Do the output and fault guarantees of one operation meet the next operation’s input requirements? | Typed joins, residual data and frame checks, level reduction and local-stochastic error bounds for a composed study. |
| [ftqc_uncompute](ftqc_uncompute/GUIDE.md) | Does measurement-based cleanup preserve the logical quantum action? Is a pebbling schedule legal, and what error does a declared mask allow? | Exact symbolic HMR checks, sampled trajectory bounds, mask and schedule checks, and source-bound replay of `ftqc_uncompute-receipt-1`. |

## Where the tools connect

- **From a code to a noisy circuit.** The audit, gadget and extractor packages reuse code constructions or binary linear algebra from `ftqc_codes`. A static distance counts data qubits supporting a logical Pauli. To count faulty circuit locations, the study also specifies preparations, gates, idles, measurements, detectors and protected outputs.
- **From extractor to audit to composition.** The extractor defines the measured logical and its spectators, and exports physical X and Z probes as native audit circuits. The compose gross-code checker replays the extractor distance receipt and uses the audit fault model to calculate the residual data Pauli after applying the split-record frame. The HGP bridge also supports replay in the enclosing circuit and export of verified-cat preparations. [HGP bridge](ftqc_extractors/AUDIT.md) · [Gross-code probes](ftqc_extractors/GROSS_CIRCUIT.md) · [Residual checker](ftqc_compose/GUIDE.md#gross-code-residual-data-and-frame)
- **From code or circuit distance to operational decoding.** Static and circuit distances count undetected faults under an ideal minimum-weight oracle or postselection. `ftqc_decoders` checks concrete polynomial-time decoders (`BP+OSD`, `BP+LSD`, `Relay-BP`, `Union-Find`, `MWPM`, `SSF`, and audited `.dem` solvers such as `Tesseract`), plugs procedural decoding into `ftqc_gadgets` beyond lookup-table limits, repairs noisy surgery syndromes via metachecks, and bounds classical feedback stall latency for `ftqc_uncompute`. [Decoder guide](ftqc_decoders/GUIDE.md)
- **From a Clifford schedule to a concatenated gadget.** `ftqc_gadgets` uses the audit fault propagator and adds encoded Bell inputs, reference qubits and decoders. Its two-level checker binds inner receipts to signed logical operations and outer faults to physical columns. It verifies the decoder radius and records any specified ancillas and the rectangle argument assumed by the reduction. A purification tree that makes such an ancilla gets a factory certificate, which the reduction accepts in place of the assumption. [Gadget guide](ftqc_gadgets/GUIDE.md)
- **From gadgets to a circuit-wide noise bound.** `ftqc_compose` checks ordered wires, residual filters and bad fault supports. Its level-reduction checker substitutes physical supports into a logical circuit. A local-stochastic bound also needs a joint fault-support model. A two-level gadget receipt can inform the local claim; a Definition 4.8 hypothesis supplied to compose remains an explicit premise. [Composition guide](ftqc_compose/GUIDE.md#level-reduction)
- **Comparing a ZX design with a location-fault study.** Run the ZX verifier on the proposed circuit and its reference. To use its result in another argument, match the gates and idles, input codespace, flag or detector handling, output frame and fault-counting convention. The ZX contract has its own verifier and noise model. [ZX guide](ftqc_zx/GUIDE.md)
- **From logical cleanup to physical audit.** The Clifford exporter lowers Hadamard-measure-reset (HMR) to X measurement and reset. A classically controlled Z follows the value of its bit, including a stored constant, an inversion or a later measurement. Conditional quantum gates beyond this direct Z feedback and non-Clifford operations require another implementation before export. The audit bridge builds a weight-one physical fault model for the chosen preparation and outputs. [Audit bridge](ftqc_uncompute/adapters/ftqc_bridge.py)
- **Pricing a logical circuit.** The resource estimator combines logical gate counts and reaction depth with the supplied ZX contract enclosures and numerical factory cost and escape assumptions. A composed logical-T receipt can justify a factory input separately. The caller supplies that number to the estimator. [Pricing implementation](ftqc_uncompute/adapters/ftqc_bridge.py) · [Logical-T composition](ftqc_compose/GUIDE.md#from-a-distilled-qubit-to-logical-t)

## Follow a gross-code measurement into the audit tool

Run these commands from the BootLoops root with Python 3.10 or later.

```sh
export PYTHONPATH="$PWD/tools"
out=$(mktemp -d)
spec=tools/ftqc_extractors/examples/gross-circuit.json
python3 -m ftqc_extractors schedule "$spec" --probe Z -o "$out/z-probe.json"
python3 -m ftqc_audit records "$out/z-probe.json" -o "$out/records.json"
python3 -m ftqc_extractors check "$spec" -o "$out/distance.json"
python3 -m ftqc_extractors verify "$spec" "$out/distance.json"
```

The audit report counts 25,048 noisy locations in this probe and checks faults through weight two. The extractor receipt combines both CSS probes, the time history of the logical readout, and a physical witness. It proves exact distance **7** for the seven-round measurement. Flipping one fixed port measurement in every round changes the reported logical Z and passes all detectors.

## Cross-block measurements

For two BB144 memories the full reading measures twelve independent corresponding logical parities. The new apparatus uses 72 data qubits and 144 check meters. With twelve syndrome rounds before and after the reading, Proposition 3.9 gives phenomenological distance at least twelve. A separate one-merge-round circuit with ideal encoded inputs and terminal measurements has circuit distance between 3 and 12 in the two complementary audit probes. Compose replays that evidence before using the split-frame residual bound. [BB144 benchmark](ftqc_extractors/GUIDE.md#two-block-bb144-benchmark)

A general fast-surgery study supplies chain maps, auxiliary cycle checks, metachecks and the ordered temporal boundaries. A block-map study instead supplies hardware ports, matching bridges and logical Clifford+T gates. Its output lists the signed measurements available on the installed couplers for each adaptive branch. [Sequence example](ftqc_extractors/examples/three-block-sequence.json) · [Compiler example](ftqc_extractors/examples/bridged-pbc.json) · [ZX comparison](ftqc_extractors/examples/homomorphic-surgery.json)

## Three questions about Steane-code operations

These commands use the same shell setup. The first pair checks the static code parameters. The second checks error correction and output closeness through order one for a Knill gadget whose ancillas are supplied as specifications. The last compares a five-ancilla X-syndrome extractor to its reference through three faults in the ZX tool's CSS circuit noise model.

```sh
python3 -m ftqc_codes certify steane -o "$out/steane.json"
python3 -m ftqc_codes verify "$out/steane.json" --full
python3 -m ftqc_gadgets analyze knill:steane,anc=spec --check ec:1 --close 1 -o "$out/gadget.json"
python3 -m ftqc_gadgets verify "$out/gadget.json" --full
python3 -m ftqc_zx equiv steane_5anc steane_spec --w 3 --noise-model circuit_css --code steane_7_1_3
```

## QR-47 as an outer code

The [[47,1,11]] quadratic-residue code supplies a candidate decoding radius of five. With lookup tables disabled, the gadget verifier proves that radius from the actual decoder check matrix and checks that each fault stays in its assigned physical column.

```sh
python3 -m ftqc_codes certify qr:47 -o "$out/qr47.json"
python3 -m ftqc_gadgets analyze steane_ec:qr:47,anc=spec --local --code-receipt "$out/qr47.json" --no-tables -o "$out/outer.json"
python3 -m ftqc_gadgets verify "$out/outer.json"
```

The report states that the gate, error-correction and any-input properties hold when faults and input errors touch at most five coordinates. Read with a coordinate as a column of inner blocks, that is the outer half of a two-level argument; `ftqc_gadgets/levels.py` joins it to the inner receipts.

That gadget takes its two ancillas as specifications. A purification tree that makes one is certified separately. One level up an idle is an inner correction, so a raw outer preparation ends with idle rounds, and only the last round counts as local.

```sh
python3 -m ftqc_gadgets analyze prep:qr:47,zero,tree=XZXZ,rounds=2 --tree --strict --light --wmin 11 -o "$out/zero.json"
python3 -m ftqc_gadgets verify "$out/zero.json"
```

The report says the sixteen-copy tree acts as its specification, with wire faults only in columns holding a failed local or check location, unless four leaves each hold a failed junk-capable location. `levels.reduction(..., factories=...)` then lists no assumed ancilla and records that count. Given independent faults across leaves and checks, the separate `Tree.bounds` calculation bounds the probability of junk conditional on acceptance, each check’s acceptance probability from below, and expected leaf consumption. [Factory guide](ftqc_gadgets/GUIDE.md)

## Trace a fault through a correction schedule

A decoder sees faults that can reach its measured block and Pauli sector by that point in the schedule. The `cones` command in `ftqc_gadgets` follows these sources through transversal gates, idles, measurements and corrections. An inner cone counts location classes per coordinate and residual errors from delivered ancillas. At the outer level the receipt lists the contributing inner events and outer ancilla sources.

Run the seven-layer example from the BootLoops root with the shell setup above.

```sh
python3 -m ftqc_gadgets cones tools/ftqc_gadgets/examples/two-level-cones.json -o "$out/cones.json"
python3 -m ftqc_gadgets verify "$out/cones.json"
```

The receipt lists eight inner and four outer decoding events. In the X sector, the outer correction sees inner events 0 and 2. After that correction the output cone sees events 2, 4 and 6. This shows exactly where the earlier X error stops contributing, assuming successful decoding. Verification recomputes the cones and checks inner propagation over GF(2) for the tested control settings. The gadget suite also exercises physical one-level Clifford examples. Converting a cone into a failure probability requires a decoder radius and a joint noise model for its sources. [Example schedule](ftqc_gadgets/examples/two-level-cones.json) · [Gadget guide](ftqc_gadgets/GUIDE.md#fault-cone-receipts)

## Carry a distillation bound into a logical T gate

The eight-qubit HGP example feeds a correcting 15-to-1 factory through a CSS encoder into the extractor’s logical T gate. Compose checks all four adaptive measurement branches. With ideal encoding and consumption, the example bounds the logical-channel half-diamond distance by 21/2000. A single encoder fault can produce an accepted logical Z. The test also exercises separate encoder and gate fault rates. [Logical T example and checks](ftqc_compose/test_logical_t.py)

```sh
PYTHONPATH=tools python3 -m unittest ftqc_compose.test_logical_t
```

## Certify and simplify measurement-based uncomputation

The bilinear circuit over GF(2³) uses HMR to erase a scratch product. The global symbolic check tests scratch restoration, clean resets, branch injectivity, outcome-independent public basis outputs, and cancellation of measurement-dependent phases. If the proof exceeds its budget or leaves injectivity unresolved, `check` refuses the certificate. Trajectories compare against the supplied oracle (identity by default). Their confidence bound and the optional small-width numerical Kraus checks have separate fields in the receipt.

```sh
dqi=tools/ftqc_uncompute/examples/dqi_rs_tiny_mbuc.kmx
python3 -m tools.ftqc_uncompute check "$dqi" -o "$out/dqi.receipt.json"
python3 -m tools.ftqc_uncompute verify "$out/dqi.receipt.json" --target "$dqi"
python3 -m tools.ftqc_uncompute pebble --length 12 --mode kmrv -o "$out/pebble12.json"
python3 -m tools.ftqc_uncompute check "$out/pebble12.json" -o "$out/pebble12.receipt.json"
python3 -m tools.ftqc_uncompute verify "$out/pebble12.receipt.json" --target "$out/pebble12.json"
```

Receipt verification needs the matching input and compares the complete regenerated evidence. The rewrite advisor recognizes seven rule families. For KickMix, `--apply` implements temporary-AND cleanup, bilinear cleanup and scratch reuse, then compares the symbolic public action. For arithmetic DAGs the verification checks the lowered circuit’s HMR obligations. Abstract nodes can still be costed and inspected. [Uncomputation guide](ftqc_uncompute/GUIDE.md) · [Packaged examples](ftqc_uncompute/examples)

## Expose a qLDPC trapping set and certify operational decoding

On the `[[144,12,12]]` bivariate bicycle gross code, a weight-3 half of a weight-6 stabilizer lies strictly inside the unique-decoding radius five yet stagnates plain iterative belief propagation. `ftqc_decoders` finds that trapping-set witness deterministically, verifies that `Relay-BP`, `BP+LSD` (with cluster rank 1) and `BP+OSD` all heal it, and computes a non-vacuous connected-cluster error bound. When the optional `uv` decoder environment is active, the same `.dem` can be solved and audited with `Tesseract`. [Decoder guide](ftqc_decoders/GUIDE.md)

```sh
python3 -m ftqc_decoders check tools/ftqc_decoders/examples/bb72_trapping_set.json -o "$out/bb_ts.receipt.json"
python3 -m ftqc_decoders verify "$out/bb_ts.receipt.json"
python3 -m ftqc_decoders check tools/ftqc_decoders/examples/bb144_peierls_cluster.json -o "$out/bb144_ler.receipt.json"
```

## Reproduce the integration checks

From BootLoops, run the eight registered suites together. Install Stim and put CaDiCaL on `PATH` (or set `FTQC_CODES_CADICAL`) to exercise both external audit checks. The optional Lean build has separate prerequisites and an opt-in command in the [code guide](ftqc_codes/GUIDE.md).

```sh
python3 run_selftests.py ftqc_codes ftqc_audit ftqc_extractors ftqc_decoders ftqc_gadgets ftqc_zx ftqc_compose ftqc_uncompute
```

## Supply states from inside a CSS memory

The simultaneous injection checker verifies a fresh CSS source from the signed logical action on its carriers. Its BB144 and HGP225 inputs include the circuit-specific first-round selection mask. Four closed Y probes use native fault propagation, Stim cross-checks and optional BP+OSD sampling. The two-state S consumer conditions a joint error distribution on batch acceptance. A T consumption study requires suitable evidence for a T source. [Injection study and replay commands](ftqc_extractors/INJECTION.md) · [Exact batch contract](ftqc_compose/GUIDE.md#simultaneous-y-source-and-ideal-consumers)
