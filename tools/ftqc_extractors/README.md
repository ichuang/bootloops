# ftqc_extractors

ftqc_extractors checks fault-tolerance claims about logical operations on encoded quantum memory. Give it a code and a measurement construction, a circuit schedule, or a logical T protocol. It produces a JSON receipt that can be verified against the same input. A failed distance target can be just as informative as a proof. Its witness identifies the faults that evade detection and damage the logical operation.

An **extractor system** couples auxiliary qubits to a quantum memory so that a choice of checks measures a logical Pauli on the data. A full extractor can address every logical Pauli on the block. A partial extractor covers a specified set. ftqc_extractors verifies finite fixed-coupler graphs for signed stabilizer codes and auxiliary-hypergraph surgeries for CSS memories. It also has specialized constructions for hypergraph product (HGP) codes and the bivariate bicycle gross code. Each construction specifies its auxiliary graph, port connections and measurement procedure. [He et al., §§4.2–4.5](https://arxiv.org/abs/2503.10390v2)

## Try a circuit proof

From the BootLoops repository root, run the following with Python 3.

```sh
export PYTHONPATH=tools
work=$(mktemp -d)
spec=tools/ftqc_extractors/examples/edge-extractor.json
python3 -m ftqc_extractors check "$spec" -o "$work/edge.json"
python3 -m ftqc_extractors verify "$spec" "$work/edge.json"
python3 -m ftqc_extractors.audit_bridge compare "$spec"
```

The example measures a weight-two logical Z on an eight-qubit HGP block. Two rounds should agree. The first command proves that every single fault either triggers detection or preserves the protected logical information. The verifier prints `verified`. The comparison also finds an undetected two-fault witness, establishing exact circuit distance two. Adding a pair of cancelling CNOTs leaves the ideal measurement intact but gives a single fault a path to a spectator logical. [Hook control](examples/edge-hook.json)

## Choose the claim

| Input kind | What the receipt establishes |
|---|---|
| `block_reading`, `gross_block_reading` | Full transversal reading of a CSS pattern, the measured subgroup, spectators and split frame. The BB144 benchmark separately audits a physical schedule. |
| `fixed_graph`, `pbc` | Routing, relative expansion and couplers of a fixed full or partial extractor, plus signed measurement settings for a bridged Clifford+T program. |
| `surgery_sequence`, `intermediate_reading` | Sequence-level phenomenological bounds with checked auxiliary complexes, distances and temporal boundaries. Intermediate time uses a separate soundness test. |
| `modular_expansion`, `surgery_zx` | Finite hypergraph expansion and thickening, or a ZX comparison of a homomorphic measurement with a surgery instrument. |
| `gross_circuit` | Circuit-level fault distance of the explicit noisy BB144 logical Z measurement, including the split frame and surviving logicals. |
| `gross_z` | Ideal unprimed logical Z measurement and exact static distances for the BB144 memory and its merged code. |
| `hgp` | Block length, logical dimension and base Pauli distance from the two classical seed matrices. |
| `column` | Whether a specified single-column Z surgery on HGP memory meets a target static code distance. |
| `canonical` | A phenomenological fault-distance lower bound for the full cyclic HGP extractor, by checking the finite hypotheses of the published construction. [Blue et al., Appendices E–F](https://arxiv.org/abs/2606.03507v2) |
| `circuit` | Correct logical Z readout, preservation of surviving logicals and a bound on undetected faulty locations in the listed schedule. |
| `t_gate`, `t_bounds` | The conditional logical T protocol with an ideal encoded magic state. The tool treats the gate output and each signed follow-on readout arm separately. |
| `t_source` | The specified in-situ preparation and all-plus acceptance rule, including an accepted raw-state fault that propagates to the logical output. |

Static distance counts physical qubits in a nontrivial logical Pauli. Circuit distance counts faulty **locations** that evade all checks and change a protected result. One faulty CNOT can affect both operands. The circuit proofs assume an ideal input codespace and terminal syndrome. Include any routing gates or idle locations whose faults the claim should cover. An exact distance needs both a lower bound and a matching witness. A search that exhausts its budget leaves the obligation open.

The default exact backend and receipt replay run with the packages in this checkout. Generating SAT/LRAT evidence for `circuit` or `t_gate` additionally requires CaDiCaL. Run `python3 run_selftests.py ftqc_extractors` for the package tests. The guides cover the [input schema and exit codes](GUIDE.md), [logical T and resource assumptions](T_PROOFS.md), and [ftqc_audit cross-checks](AUDIT.md).

## Measure a gross-code logical Z

The 12 port checks of the unprimed Z setting select one logical observable of the [[144,12,12]] bivariate bicycle memory. ftqc_extractors exports the ideal measurements, split-frame rule and merged CSS checks. The distance receipt proves that the resulting [[162,11,12]] code preserves distance 12 in both Pauli sectors. Its four lower bounds and four weight-12 witnesses can be replayed without an external solver. [Run the gross-code example](GROSS.md)

With noisy syndrome extraction, the seven-round schedule has exact circuit distance 7 and the twelve-round schedule has exact distance 10. ftqc_extractors exports both CSS probes to ftqc_audit and verifies the lower bounds and matching physical faults. [Noisy schedule and replay instructions](GROSS_CIRCUIT.md)
