# `ftqc_uncompute`

Logical-level checks for measurement-based uncomputation (Hadamard-measure-
reset with phase fix-ups), spooky and parallel spooky pebbling, and
superposition-masked arithmetic. The symbolic check proves, for every input and
measurement outcome, that scratch is restored, resets are clean, outputs do not
depend on outcomes, phases cancel and branches are injective; the `blocks`
command does the same locally for arithmetic too wide for one polynomial. A
rewrite advisor applies seven simplification rules and checks the result.

```sh
python3 tools/ftqc_uncompute/selftest.py
python3 -m tools.ftqc_uncompute check tools/ftqc_uncompute/examples/gidney_adder_4bit.kmx -o /tmp/adder.receipt.json
python3 -m tools.ftqc_uncompute verify /tmp/adder.receipt.json --target tools/ftqc_uncompute/examples/gidney_adder_4bit.kmx
python3 -m tools.ftqc_uncompute suggest my_circuit.kmx --apply -o my_circuit_opt.kmx
```

Physical faults are out of scope here; Clifford and HMR parts lower to
`ftqc_faults`. See the [guide](GUIDE.md) and [`../FTQC.md`](../FTQC.md).
