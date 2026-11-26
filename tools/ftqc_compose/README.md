# ftqc_compose

Check bad-support families, join typed gadgets in series or parallel, and replay a diagonal-Pauli proof for two one-edge HGP measurements on the same encoded memory. The receipt records which inputs, decoders, physical locations, ideal boundaries and unresolved premises its conclusion uses.

From the BootLoops root:

```sh
export PYTHONPATH=tools
out=$(mktemp -d)
echo '{"kind":"hgp-repeat","variant":"edge-extractor"}' > "$out/spec.json"
python3 -m ftqc_compose check "$out/spec.json" -o "$out/receipt.json"
python3 -m ftqc_compose verify "$out/spec.json" "$out/receipt.json"
python3 run_selftests.py ftqc_compose
```

See the [package guide](GUIDE.md) for the contract language, evidence and limits, and the [FTQC overview](../FTQC.md) for related tools.

The [level-reduction and noise rules](GUIDE.md#level-reduction) also check unconditional circuits with typed interfaces and calculate exact rational local-stochastic bounds. The 15-to-1 example verifies all single-resource corrections and the 105 two-resource bad supports under fault-free processing.

The [gross-code residual checker](GUIDE.md#gross-code-residual-data-and-frame) transports data and split-frame errors to the BB144 output. The [distilled logical T study](GUIDE.md#from-a-distilled-qubit-to-logical-t) checks the encoding interface, all four adaptive branches, and separate bounds for factory, encoder and online faults.
