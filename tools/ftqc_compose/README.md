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
