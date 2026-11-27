# `ftqc_uncompute`

Certified verification and simplification toolkit for **measurement-based uncomputation (MBUC)**, **spooky & parallel spooky pebbling**, and **superposition-masked quantum arithmetic**.

See [`GUIDE.md`](GUIDE.md) for full documentation, refusal codes, known hazards, and citations, and [`../FTQC.md`](../FTQC.md) for how `ftqc_uncompute` connects to [`ftqc_codes`](../ftqc_codes/GUIDE.md), [`ftqc_audit`](../ftqc_audit/GUIDE.md), [`ftqc_gadgets`](../ftqc_gadgets/GUIDE.md), [`ftqc_extractors`](../ftqc_extractors/README.md), [`ftqc_zx`](../ftqc_zx/GUIDE.md), and [`ftqc_compose`](../ftqc_compose/GUIDE.md).

## Quickstart

```sh
# Run the acceptance battery, including verifier regressions
python3 tools/ftqc_uncompute/selftest.py

# Certify a .kmx circuit and re-verify its zero-float receipt
python3 -m tools.ftqc_uncompute check tools/ftqc_uncompute/examples/gidney_adder_4bit.kmx -o /tmp/adder.receipt.json
python3 -m tools.ftqc_uncompute verify /tmp/adder.receipt.json --target tools/ftqc_uncompute/examples/gidney_adder_4bit.kmx

# Apply supported KickMix rewrites
python3 -m tools.ftqc_uncompute suggest my_circuit.kmx --apply -o my_circuit_opt.kmx
```
