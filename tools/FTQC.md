# Fault-tolerant quantum computing in BootLoops

Eight tools are available on `ftqc_integration`.

The arrows show what one tool supplies to another in selected studies.

![FTQC tool relationships across the eight packages](figures/ftqc-tools.svg)

| Tool | Input and question | Evidence |
|---|---|---|
| [ftqc_codes](ftqc_codes/GUIDE.md) | Given CSS check matrices or a high-rate algebraic seed, which logical errors can the memory detect and how are its logical slots addressed? | Code parameters, construction ancestry (`ftqc_codes-construction-1`), paired logical bases (`ftqc_codes-basis-1`), lift symmetries, bounds on X and Z distance, witnesses and hierarchical decoding guarantees. |
| [ftqc_extractors](ftqc_extractors/GUIDE.md) | Which logical measurements can a fixed apparatus perform, and what fault distance does a surgery sequence preserve? | Fixed-coupler graphs and bridges, transversal block readings, padded hypergraph surgery, physical fault witnesses, and resource-conditional logical T receipts. |
| [ftqc_faults](ftqc_faults/GUIDE.md) | How do faults propagate in this circuit, and can an implementation replace a reference under a specified fault model? | Physical fault signatures, distance intervals and witnesses from the audit engine. ZX equivalence, synthesis, flag decoder tables, and stochastic contracts from the ZX engine. |
| [ftqc_decoders](ftqc_decoders/GUIDE.md) | Does a polynomial-time decoder achieve the code or circuit distance without trapping sets, and what operational error rate does it certify? | Operational radius proofs, split-stabilizer trapping-set witnesses, two-stage metacheck repair, Peierls cluster LER bounds, Stim `.dem` export and audited external (`Tesseract`, `PyMatching`, `ldpc`) runs. |
| [ftqc_gadgets](ftqc_gadgets/GUIDE.md) | What does a decoded Clifford gadget do to an unknown encoded input with incoming errors and internal faults? | Residual-error bounds, component refinement, two-level checks on concatenated codes, and the faults that can reach each decoding event in a scheduled program. |
| [ftqc_compose](ftqc_compose/GUIDE.md) | Do the output and fault guarantees of one operation meet the next operation’s input requirements? | Typed joins, residual data and frame checks, level reduction and local-stochastic error bounds for a composed study. |
| [ftqc_uncompute](ftqc_uncompute/GUIDE.md) | Does measurement-based cleanup preserve the logical quantum action? Is a pebbling schedule legal, and what error does a declared mask allow? | Exact symbolic HMR checks, sampled trajectory bounds, mask and schedule checks, and source-bound replay of `ftqc_uncompute-receipt-1`. |
| [ftqc_magic](ftqc_magic/GUIDE.md) | What output quality, yield, and active-qubit footprint can a magic-state factory or cultivation pipeline deliver, and are its output correlations safe to compose? | Exact $\mathbb{Z}[\omega_8]$ / $\mathbb{Q}(\sqrt{2})$ / rational distillation & stabilizer-collapse proofs, recycled check schedules, $O(n^3)$ Clifford-hook cultivation distances, Caliper open-boundary escape gaps, single-shot qLDPC `CCZ` contracts, and `ftqc_magic-receipt-1` exports. |

## Where the tools connect

- **From a code to a noisy circuit.** The audit, gadget and extractor packages reuse code constructions or binary linear algebra from `ftqc_codes`. A static distance counts data qubits supporting a logical Pauli. To count faulty circuit locations, the study also specifies preparations, gates, idles, measurements, detectors and protected outputs.
- **From extractor to audit to composition.** The extractor defines the measured logical and its spectators, and exports physical X and Z probes as native audit circuits. The compose gross-code checker replays the extractor distance receipt and uses the audit fault model to calculate the residual data Pauli after applying the split-record frame. The HGP bridge also supports replay in the enclosing circuit and export of verified-cat preparations. [HGP bridge](ftqc_extractors/AUDIT.md) · [Gross-code probes](ftqc_extractors/GROSS_CIRCUIT.md) · [Residual checker](ftqc_compose/GUIDE.md#gross-code-residual-data-and-frame)
- **From code or circuit distance to operational decoding.** Static and circuit distances count undetected faults under an ideal minimum-weight oracle or postselection. `ftqc_decoders` checks concrete polynomial-time decoders (`BP+OSD`, `BP+LSD`, `Relay-BP`, `Union-Find`, `MWPM`, `SSF`, and audited `.dem` solvers such as `Tesseract`), plugs procedural decoding into `ftqc_gadgets` beyond lookup-table limits, and repairs noisy surgery syndromes via metachecks. To estimate feedback stalls for `ftqc_uncompute`, it measures decoder work on every single-fault-mechanism syndrome and applies the supplied tick costs and reaction depth. The estimate covers that set of syndromes. [Decoder guide](ftqc_decoders/GUIDE.md)
- **From a Clifford schedule to a concatenated gadget.** `ftqc_gadgets` uses the audit fault propagator and adds encoded Bell inputs, reference qubits and decoders. Its two-level checker binds inner receipts to signed logical operations and outer faults to physical columns. It verifies the decoder radius and records any specified ancillas and the rectangle argument assumed by the reduction. A purification tree that makes such an ancilla gets a factory certificate, which the reduction accepts in place of the assumption. [Gadget guide](ftqc_gadgets/GUIDE.md)
- **From gadgets to a circuit-wide noise bound.** `ftqc_compose` checks ordered wires, residual filters and bad fault supports. Its level-reduction checker substitutes physical supports into a logical circuit. A local-stochastic bound also needs a joint fault-support model. A two-level gadget receipt can inform the local claim; a Definition 4.8 hypothesis supplied to compose remains an explicit premise. [Composition guide](ftqc_compose/GUIDE.md#level-reduction)
- **Comparing a ZX design with a location-fault study.** Run the ZX verifier on the proposed circuit and its reference. To use its result in another argument, match the gates and idles, input codespace, flag or detector handling, output frame and fault-counting convention. The ZX contract has its own verifier and noise model. [ZX guide](ftqc_zx/GUIDE.md)
- **From logical cleanup to physical audit.** The Clifford exporter lowers Hadamard-measure-reset (HMR) to X measurement and reset. A classically controlled Z follows the value of its bit, including a stored constant, an inversion or a later measurement. Conditional quantum gates beyond this direct Z feedback and non-Clifford operations require another implementation before export. The audit bridge builds a weight-one physical fault model for the chosen preparation and outputs. [Audit bridge](ftqc_uncompute/adapters/ftqc_bridge.py)
- **Pricing a logical circuit.** The resource estimator combines logical gate counts and reaction depth with the supplied ZX contract enclosures and numerical factory cost and escape assumptions. A composed logical-T receipt can justify a factory input separately. The caller supplies that number to the estimator. [Pricing implementation](ftqc_uncompute/adapters/ftqc_bridge.py) · [Logical-T composition](ftqc_compose/GUIDE.md#from-a-distilled-qubit-to-logical-t)
- **From cultivation and distillation to composed magic contracts.** `ftqc_magic` verifies binary/PPR, $\mathrm{GF}(2^s)$ trace, and non-triorthogonal $\mathbb{Z}[\omega_8]$ stabilizer factories, mid-circuit check recycling, $O(n^3)$ Clifford-hook cultivation (`d_T` vs `d_S`), Caliper open-boundary escape gaps (cross-checked with `ftqc_decoders` Blossom), and single-shot 3D tricycle $\to$ 4D quadcycle `CCZ` factories (repaired via `ftqc_decoders` metachecks). Its pipeline composer enforces `INDEPENDENT_SINGLE`, `BATCH_CORRELATED`, and `CATALYST_POISONED` rules before exporting `ftqc_compose` ancilla contracts and `ftqc_uncompute` pricing records. [Magic factory guide](ftqc_magic/GUIDE.md)

## Compare a surface-code memory or lattice-surgery measurement

The rotated patch uses one data qubit per grid point, with one meter per check
during a memory round. One such round comprises preparation, four CNOT layers
and measurement. For a distance-three memory, the standard and hook-aligned
CNOT orders have circuit distances three and two in the two-round X and Z
probes. The sector decoder projects each complete physical fault onto the
measured CSS sector and applies equal-weight MWPM. Its report counts failures
among all physical single-location alternatives. The exact audit measures
undetected distance separately.

```sh
export PYTHONPATH="$PWD/tools"
out=$(mktemp -d)
python3 -m ftqc_codes certify rotated_surface:3 -o "$out/code.json"
python3 -m ftqc_codes verify "$out/code.json" --full
python3 -m ftqc_faults audit analyze surface:3,2,z,bad -o "$out/hook.json"
python3 -m ftqc_decoders.surface_memory --distance 3 --rounds 3 --basis Z
for axis in xx zz; do
    python3 -m ftqc_extractors.surgery_bench "surface-$axis" -o "$out/$axis.json"
    python3 -m ftqc_extractors check "$out/$axis.json" -o "$out/$axis-proof.json"
    python3 -m ftqc_extractors verify "$out/$axis.json" "$out/$axis-proof.json"
done
```

The installation accepts odd distances of at least three. It begins with one
split-code round, measures the merged code for R rounds, then remeasures the
split code for R rounds. The shared seam meters measure opposing half checks
in successive subrounds. Their complementary seam checks encode the signed
parity. The last split record controls the recovery on live data. Gates,
feedback and idles carry physical faults. The initial encoded state and terminal
syndrome are ideal. Static square-grid placements and a gate range of one
check every coupler. Every ideal measurement branch preserves the named
surviving logical pair. The distance receipt includes single-location and
two-location exclusions and a forward-replayed physical witness, or a smaller
counterexample. At d = R = 3 the exact circuit distance is three in both bases.
A larger patch is opt-in with `--distance 5`; `--rounds` overrides the default
R = d.

To price a decoded measurement, carry the installed Gadgets proof into Compose:

```sh
python3 - "$out" <<'PY'
import json, pathlib, sys
directory = pathlib.Path(sys.argv[1])
installed = json.loads((directory / "zz-proof.json").read_text())
study = dict(kind="decoded-surgery", receipt=installed["result"]["gadget"],
             incoming_weight=0, p="1/100000")
(directory / "channel.json").write_text(json.dumps(study))
PY
python3 -m ftqc_compose check "$out/channel.json" -o "$out/bound.json"
python3 -m ftqc_compose verify "$out/channel.json" "$out/bound.json"
```

The decoder uses one joint syndrome table for public bits and the residual
Pauli frame. It scores every physical single-location alternative with its
within-location probability, including correlated CNOT faults. For independent
sites with common fault probability p, Compose adds a pair union bound for
all higher orders. The resulting public, residual and joint failure bounds
assume a clean encoded input. Residual failure means the corrected data Pauli
lies outside the stabilizers and measured branch subgroup. The research
checkout's `studies/surface-baseline` contains the matched HGP comparison,
seeded sampling and twelve-parity resource accounting.

## Follow a gross-code measurement into the audit tool

Run these commands from the BootLoops root with Python 3.10 or later.

```sh
export PYTHONPATH="$PWD/tools"
out=$(mktemp -d)
spec=tools/ftqc_extractors/examples/gross-circuit.json
python3 -m ftqc_extractors schedule "$spec" --probe Z -o "$out/z-probe.json"
python3 -m ftqc_faults audit records "$out/z-probe.json" -o "$out/records.json"
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
python3 -m ftqc_faults zx equiv steane_5anc steane_spec --w 3 --noise-model circuit_css --code steane_7_1_3
```

The five-ancilla bridge verifies both sector-specific ZX derivations and cross-checks all 455 single-location fault modes in each physical schedule. The two complementary terminal probes recover both components of the raw data Pauli and compare every measurement-record displacement. Replay the study with the same shell setup. [Model, boundaries and fault weights](ftqc_faults/GUIDE.md#five-ancilla-steane-extractor)

```sh
python3 -m ftqc_faults steane -o "$out/steane-bridge.json"
python3 -m ftqc_faults steane --verify "$out/steane-bridge.json"
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

The receipt lists eight inner and four outer decoding events. In the X sector, the outer correction sees inner events 0 and 2. After that correction the output cone sees events 2, 4 and 6. This shows exactly where the earlier X error stops contributing, assuming successful decoding. Verification recomputes the cones and checks inner propagation over GF(2) for the tested control settings. The gadget suite also exercises physical one-level Clifford examples. A `cone-bound` study in `ftqc_compose` turns the cone into a conditional upper bound on decoding failure using a declared decoder radius, noise model and ancilla contracts. It supports independent faults at one or two levels and local-stochastic faults at one level. The two-level calculation also assumes independent columns. A `filter-contract` study derives the ancilla contracts themselves from a certified purification tree, with conditioning on acceptance, at the physical level and for a tree of encoded blocks. A `volume` study prices the same schedule in expected physical locations from a ledger of its block-ticks, corrections and delivered resources, with ancilla costs taken from the filter studies. [Volume](ftqc_compose/GUIDE.md#volume-of-a-schedule) · [Filter calculus](ftqc_compose/GUIDE.md#filter-calculus-for-purification-trees) · [Tail bounds and premises](ftqc_compose/GUIDE.md#tail-bounds-over-fault-cones) [Example schedule](ftqc_gadgets/examples/two-level-cones.json) · [Gadget guide](ftqc_gadgets/GUIDE.md#fault-cone-receipts)

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

Wide arithmetic exceeds the global symbolic budget: a 16-bit ripple adder is already refused. The `blocks` command certifies such programs compositionally. It names intermediate values, proves every deletion locally, proves the arithmetic statement as one exact integer relation, and binds a correction schedule in the format of `ftqc_gadgets cones` to the program, so the schedule whose failure probability `cone-bound` bounds is shown to compute the stated function. [Cuts and block programs](ftqc_uncompute/GUIDE.md#wide-arithmetic-cuts-and-block-programs)

Receipt verification needs the matching input and compares the complete regenerated evidence. The rewrite advisor recognizes seven rule families. For KickMix, `--apply` implements temporary-AND cleanup, bilinear cleanup and scratch reuse, then compares the symbolic public action. For arithmetic DAGs the verification checks the lowered circuit’s HMR obligations. Abstract nodes can still be costed and inspected. [Uncomputation guide](ftqc_uncompute/GUIDE.md) · [Packaged examples](ftqc_uncompute/examples)

## Expose a qLDPC trapping set and certify operational decoding

On the `[[144,12,12]]` BB144 gross code, plain belief propagation stalls on a weight-three half of a weight-six stabilizer. The first receipt replays this witness and verifies that Relay-BP, BP+LSD (maximum cluster rank one), and BP+OSD each correct it. [Trapping-set input](ftqc_decoders/examples/bb72_trapping_set.json)

The second receipt computes a Peierls upper bound on the X-sector logical error rate at physical error probability 1/10,000. It assumes that disjoint fault clusters in the interaction graph decode independently and that the decoder corrects every connected cluster of size at most five. [Peierls input](ftqc_decoders/examples/bb144_peierls_cluster.json) · [Decoder guide](ftqc_decoders/GUIDE.md#known-gaps)

```sh
python3 -m ftqc_decoders check tools/ftqc_decoders/examples/bb72_trapping_set.json -o "$out/bb_ts.receipt.json"
python3 -m ftqc_decoders verify "$out/bb_ts.receipt.json"
python3 -m ftqc_decoders check tools/ftqc_decoders/examples/bb144_peierls_cluster.json -o "$out/bb144_ler.receipt.json"
python3 -m ftqc_decoders verify "$out/bb144_ler.receipt.json"
```

## Certify a recycled factory and a flagged cultivation pipeline

On the Bravyi–Haah `49T -> 1T` ($d = 5$, $1411 p^5$) factory, mid-circuit check measurement and reset compresses 14 logical rows down to **5** active qubits. On `MSC-3` magic-state cultivation, exact $O(n^3)$ simulation of propagated Clifford errors $H_{\pm} = (X \pm Y)/\sqrt{2}$ exposes 4 weight-two malignant faults on $|T\rangle$ ($d_T = 2$ despite $d_S = 3$ on $|S\rangle$) when unflagged, and certifies $d_T = 3$ once two $Z$-flag ancillas are added and chained through Caliper open-boundary escape into a recycled Litinski 15-to-1 stage. [Factory examples](ftqc_magic/examples) · [Magic guide](ftqc_magic/GUIDE.md)

```sh
python3 -m ftqc_magic check tools/ftqc_magic/examples/bh49_recycled_5qubits.json -o "$out/bh49.receipt.json"
python3 -m ftqc_magic verify tools/ftqc_magic/examples/bh49_recycled_5qubits.json "$out/bh49.receipt.json"
python3 -m ftqc_magic check tools/ftqc_magic/examples/flagged_cultivation_caliper_15to1.json -o "$out/cult_pipe.receipt.json"
python3 -m ftqc_magic verify tools/ftqc_magic/examples/flagged_cultivation_caliper_15to1.json "$out/cult_pipe.receipt.json"
```

## Reproduce the integration checks

From BootLoops, run the nine suites, including the separate audit and ZX engine batteries, then run the cross-engine checks. Install Stim and put CaDiCaL on `PATH` (or set `FTQC_CODES_CADICAL`) to exercise external SAT+LRAT and stabilizer checks. The optional Lean build has separate prerequisites and an opt-in command in the [code guide](ftqc_codes/GUIDE.md).

```sh
python3 run_selftests.py ftqc_codes ftqc_audit ftqc_extractors ftqc_decoders ftqc_gadgets ftqc_zx ftqc_compose ftqc_uncompute ftqc_magic
python3 run_selftests.py ftqc_faults
```

## Addressable high-rate blocks and single-slot Magic delivery

`ftqc_codes` emits replayable `ftqc_codes-construction-1` and `ftqc_codes-basis-1` receipts binding ordered matrices, named symplectic logical slots, stabilizer dressings, and lift symmetries (`LP`, `2BGA`, Quantum Tanner, `radial90`, and pair-partition `PP256` / `APM1152`). `ftqc_extractors.addressable_surgery` installs signed slot measurements and classical feed-forward on explicit coupler inventories, flattens schedules through `ftqc_faults` and `ftqc_decoders`, verifies typed residual contracts in `ftqc_compose`, and prices reaction stalls in `ftqc_uncompute`. `ftqc_compose.addressable_magic` delivers verified single-output `ftqc_magic` $T$ states via `occupied_swap` into a named resource slot while preserving untouched data slots and Bell references. [High-rate QEC roadmap](../PLANS/high-rate-qec.md)

```sh
python3 -m ftqc_codes construct pp256 -o "$out/pp256.construct.json"
python3 -m ftqc_codes basis lp6_2 --method spectral -o "$out/lp6_2.basis.json"
PYTHONPATH=tools python3 -m unittest ftqc_extractors.test_highrate_qec ftqc_compose.test_addressable_magic
```

## Supply states from inside a CSS memory

The simultaneous injection checker verifies a fresh CSS source from the signed logical action on its carriers. Its BB144 and HGP225 inputs include the circuit-specific first-round selection mask. Four closed Y probes use native fault propagation, Stim cross-checks and optional BP+OSD sampling. The two-state S consumer conditions a joint error distribution on batch acceptance. A T consumption study requires suitable evidence for a T source. [Injection study and replay commands](ftqc_extractors/INJECTION.md) · [Exact batch contract](ftqc_compose/GUIDE.md#simultaneous-y-source-and-ideal-consumers)
