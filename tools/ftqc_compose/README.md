# ftqc_compose

Joins operations certified by the other FTQC tools and computes what follows
for the whole: typed joins of wires, codes, logical bases and decoders;
bad-support families and level reduction after He, Nguyen and Pattison;
probability bounds under declared local-stochastic or independent noise; tail
bounds over the fault cones of a correction schedule; contracts for
purification trees, under independent faults or under local-stochastic noise
alone; the closure of a two-level program, a bound whose only probabilistic
premise is the noise condition; regions of several noise rates; schedule
volumes; and factory bounds carried into a logical T. Every result is a receipt that lists the premises it rests on, and
is marked conditional when it rests on a quoted theorem, a supplied contract or
a declared noise law.

From the BootLoops root:

```sh
export PYTHONPATH=tools
out=$(mktemp -d)
echo '{"kind":"hgp-repeat","variant":"edge-extractor"}' > "$out/spec.json"
python3 -m ftqc_compose check "$out/spec.json" -o "$out/receipt.json"
python3 -m ftqc_compose verify "$out/spec.json" "$out/receipt.json"
python3 run_selftests.py ftqc_compose
```

The example checks two nondestructive measurements of logical Z on the same
eight-qubit HGP memory and their composition; `edge-hook` replays the negative
control. See the [guide](GUIDE.md) for every study kind, its premises and its
limits, and [`../FTQC.md`](../FTQC.md) for how the tools fit together.
