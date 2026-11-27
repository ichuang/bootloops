# Composable signed Pauli measurements

Gadgets certifies a physical projective measurement on an arbitrary encoded input. Extractors binds these circuits to a shared apparatus, allocates record bits, inserts feedback and checks the complete resulting circuit. Compose can consume the accepted diagonal-Pauli contract or calculate a joint accepted error law. For the delivered features, outstanding work and original paper links, see the [milestone status for the surgery round](../../PLANS/extractor-surgery.md#implementation-status-at-6613a17).

## Run a measurement

From the BootLoops root, set `PYTHONPATH=tools`. The commands below generate a finite linear-PCPP example with an explicit proof coordinate and extract its gate-level single-fault witness. `surgery` exits 0 for an exclusion through the requested order, 1 for a counterexample, 2 for malformed input and 3 for an unfinished search. `verify` replays the supplied receipt, including a witnessed failure.

```sh
python3 -m ftqc_extractors.surgery_bench pcpp -o design.json
python3 -m ftqc_gadgets build-surgery design.json --geometry geometry.json -o instrument.json
python3 -m ftqc_gadgets surgery instrument.json --order 1 -o receipt.json
python3 -m ftqc_gadgets verify receipt.json
python3 -m ftqc_gadgets.surgery_operational residual instrument.json --incoming 1 --order 1 --cap 1000000
python3 -m ftqc_gadgets.surgery_operational sample instrument.json --p 1/1000 --seed 42 --shots 2048 --export-prefix physical
```

A design can use `homomorphic` for an encoded ancilla, `gauging` for an explicit connected graph, `enhanced` for a paired Tanner restriction with finite cut augmentation, `bridge` for a matching between two connected components, `adapter` for a SkipTree repetition adapter, `parsimonious` for a cellulated cone, `branch` for coherent repetition-tree fanout, `branch-code` for a CSS branching tree of mapping cylinders, or `linear` for explicit binary maps. The graph and linear merge constructors return deformed CSS matrices, logical images and split corrections. Independent mixed-Pauli subgroups sharing compatible single-site bases can use the branch constructor. The small rotated-patch tests compare bridge and enhanced constructions.

Qubit zero is the first letter after the sign in a Pauli string and the least significant bit in binary matrix rows. Data wires occupy `0..n-1`. A standalone spec gives `n`, signed `checks` and independent signed `targets`, disjoint gate `ticks`, one `readouts` record-index parity and sign per target, zero-parity `detectors`, and the boolean `terminal_syndrome`. Measurement records start at zero. A leading exclamation point on an operation declares a noiseless location. Optionally bind an ordered `basis.input` list of named signed X/Z pairs and a complete `basis.surviving` list of their names. Verification checks symplectic pairing, independence and preservation on all branches. The example HGP export is available through `ftqc_extractors export-surgery`.

The reference is an encoded Bell state. After the instrument, each public bit must agree with both the input-reference target and the post-split target. The unmeasured centralizer must retain its input correlations even after private records are discarded. The checker transports detector, readout and residual signatures for all Pauli alternatives at each physical gate. A measurement fault flips its recorded outcome. Physical residuals are reduced modulo the checks and measured targets, so a nonzero representative weight is an upper bound on the minimum representative weight. The scope of a finite bound includes the stated ideal terminal syndrome and the explicitly scheduled idles. Two-location evidence can also be obtained with `surgery_faults.pair_search`.

The optional sampler exports a Stim circuit and detector-error model, uses an independent-location depolarizing model and reports a Wilson interval. Its default decoder matches complete single-mechanism signatures by syndrome and predicts zero on unseen syndromes. Pass `--decoder bp_osd` to use the registered joint decoder on small cases. The deterministic enumeration and the exact categorical law keep alternatives at a single location mutually exclusive. `law` uses rational arithmetic and returns an inconclusive status when its state cap is exceeded.

## Install and compose

```sh
python3 -m ftqc_extractors.surgery_bench two-block -o installed.json
python3 -m ftqc_extractors check installed.json -o joint.json
python3 -m ftqc_extractors verify installed.json joint.json
```

An `installed_surgery` setting includes the complete Gadgets receipt and an injective placement. An optional control `[earlier_bit, [global_sites...]]` wraps a local instrument in a bit-controlled Z frame and updates its public parity. Batches can share a time slot when physical supports are disjoint. The compiler binds a coupler inventory, replays the flattened schedule on all data, checks the ordered joint target subgroup and reports peak live wires. Optional `transport` supplies two-dimensional integer coordinates, synchronous straight-line movements, per-layer durations and an interaction radius. The timing model tests moving collisions and two-qubit range. Its durations are input assumptions. See `shuttle-packed` and `shuttle-serial` for matched finite cases.

For a PBC input, optional `surgery_installation` supplies `couplers`, `idles` and a list of placed single-output Gadgets settings. The compiler matches each signed request modulo memory stabilizers, chooses the shortest verified setting, records any public sign correction and serializes its physical measurement ticks. Coverage applies to the requested program branches. Resource-state preparations and the joint fault tolerance of noncommuting PBC steps retain their existing premises. For a commuting adaptive workload, `installed_surgery` can instead certify the enclosing physical handoff as one instrument.

Compose's `surgery` leaf takes `instance`, ordered data `wires` and a Gadgets `receipt`; the placed receipt is at `result.gadget` of an Extractors check. Optional `incoming` requires a fresh joint incoming/internal enumeration and zero accepted damage everywhere admitted. The separate Compose input `{"kind":"surgery-channel","receipt":...,"incoming_weight":1,"p":"1/1000","cap":100000}` calculates the exact correlated output law for clean inputs and a worst-case conditional damage envelope over admitted input weights. It charges each physical location once on the enclosing circuit.

## Reproduce the larger cases

`ftqc_extractors.surgery_bench --help` lists finite HGP, cyclic, radial, PCPP, product-connection and architecture cases. Padded HGP emits a sequence-checker input. Radial, PCPP and product-connection cases emit Gadgets designs. Installed cases emit Extractors check inputs. Set `noise_locations` to `phenomenological` on a linear or graph merge design to charge one data/merged-qubit error per qubit at each round boundary and noisy syndrome/split measurements, with ideal coupling gates. Match the noise locations and target subgroup before comparing costs.

The companion research checkout banks results and full replay commands under `studies/surgery/`, including checksums and a downloader for the authors' trapped-ion QIR and shot files, raw-histogram analysis and a separate simulated and decoded comparison. The optional Stim backend accelerates large Choi tests. The finite radial example uses the paper's protographs with a dense rank-one equivariant map. The cellulated cone and coherent fanout have finite witnesses for their generated geometries. Published optimal-distance, asymptotic-overhead and native scheduling claims require additional hypotheses.

