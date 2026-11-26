# FTQC

FTQC checks fault-tolerance claims about logical operations on encoded quantum memory. Give it a code and a measurement construction, a circuit schedule, or a logical T protocol. It produces a JSON receipt that can be verified against the same input. A failed distance target can be just as informative as a proof. Its witness identifies the faults that evade detection and damage the logical operation.

An **extractor system** couples auxiliary qubits to a quantum memory so that a choice of checks measures a logical Pauli on the data. A full extractor can address every logical Pauli on the block. A partial extractor covers a specified set. Hypergraph product (HGP) describes the **memory code family** used by the extractor constructions currently implemented in FTQC. Other memories can use different auxiliary graphs, port maps and measurement schedules. [He et al., §§4.2–4.5](https://arxiv.org/abs/2503.10390v2)

## Try a circuit proof

From the BootLoops repository root, run the following with Python 3.

```sh
export PYTHONPATH=tools
work=$(mktemp -d)
spec=tools/ftqc/examples/edge-extractor.json
python3 -m ftqc check "$spec" -o "$work/edge.json"
python3 -m ftqc verify "$spec" "$work/edge.json"
python3 -m ftqc.faultcert_bridge compare "$spec"
```

The example measures a weight-two logical Z on an eight-qubit HGP block. Two rounds should agree. The first command proves that every single fault either triggers detection or preserves the protected logical information. The verifier prints `verified`. The comparison also finds an undetected two-fault witness, establishing exact circuit distance two. Adding a pair of cancelling CNOTs leaves the ideal measurement intact but gives a single fault a path to a spectator logical. [Hook control](examples/edge-hook.json)

## Choose the claim

| Input kind | What the receipt establishes |
|---|---|
| `hgp` | Block length, logical dimension and base Pauli distance from the two classical seed matrices. |
| `column` | Whether a specified single-column Z surgery on HGP memory meets a target static code distance. |
| `canonical` | A phenomenological fault-distance lower bound for the full cyclic HGP extractor, by checking the finite hypotheses of the published construction. [Blue et al., Appendices E–F](https://arxiv.org/abs/2606.03507v2) |
| `circuit` | Correct logical Z readout, preservation of surviving logicals and a bound on undetected faulty locations in the listed schedule. |
| `t_gate`, `t_bounds` | The conditional logical T protocol with an ideal encoded magic state. The tool treats the gate output and each signed follow-on readout arm separately. |
| `t_source` | The specified in-situ preparation and all-plus acceptance rule, including an accepted raw-state fault that propagates to the logical output. |

Static distance counts physical qubits in a nontrivial logical Pauli. Circuit distance counts faulty **locations** that evade all checks and change a protected result. One faulty CNOT can affect both operands. The circuit proofs assume an ideal input codespace and terminal syndrome. Include any routing gates or idle locations whose faults the claim should cover. An exact distance needs both a lower bound and a matching witness. A search that exhausts its budget leaves the obligation open.

The default exact backend and receipt replay run with the packages in this checkout. Generating SAT/LRAT evidence for `circuit` or `t_gate` additionally requires CaDiCaL. Run `python3 run_selftests.py ftqc` for the package tests. The guides cover the [input schema and exit codes](GUIDE.md), [logical T and resource assumptions](T_PROOFS.md), and [faultcert cross-checks](FAULTCERT.md).

## A bivariate bicycle partial extractor as the next system

A concrete second target is the **103-ancilla gauge-fixed surgery system on the [[144,12,12]] gross code**. Its X and Z auxiliary graphs, selectable connections to the data, and bridge support eight named Pauli measurements, one at a time. He et al. identify this hardware as a partial extractor. Cross et al. give its logical supports, gauge checks and wiring, and report distance 12 for the merged codes. Their circuit-noise simulations are a separate comparison. [Extractor interpretation, §4.4](https://arxiv.org/abs/2503.10390v2) · [Gross-code construction, Table 1 and §4](https://arxiv.org/abs/2407.18393v4)

For a **first certificate**, choose the unprimed logical Z from Table 1. The sibling `qeccert` package already exports the memory via `python3 -m qeccert spec bb144`. Map its qubit ordering to the paper's left/right lattice coordinates, then build the Z measurement's merged checks and port connections. Verify check commutation, the measured logical class and preservation of the other eleven logical qubits. Target a replayable exclusion of surviving logicals below weight 12 and a physical weight-12 witness. A lower-weight witness would pinpoint a discrepancy in the proposed reconstruction. [Existing memory constructor](../qeccert/codes.py) · [Logical representatives, Eqs. 61–64, and gauge checks, Figs. 8–9](https://arxiv.org/abs/2407.18393v4)

That first target is a **static** certificate. A circuit-level extension would take a fully scheduled repeated measurement with explicit detectors and fault locations, check the signed readout, and replay candidate hooks with `faultcert`. Each of the other seven port settings would need its own receipt or a verified symmetry relating it to one already checked. Bicycle extractor support is proposed work.
