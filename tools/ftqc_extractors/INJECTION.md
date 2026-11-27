# Simultaneous in-situ injection

The source implements Sections III–IV of Liu et al., arXiv:2604.05126v1. A signed logical-intertwining check applies to arbitrary input states on the specified carriers. The noisy benchmarks use logical Y eigenstates and ideal terminal Pauli-product measurements. The general carrier Clifford is emitted explicitly. Each source prepares a new memory block. Other logical states in that block have no preservation guarantee.

## Reproduce the source and physical audit

From BootLoops, set `PYTHONPATH=tools`. Run `python -m ftqc_extractors check tools/ftqc_extractors/examples/bb144-injection-4.json`. For the noisy study run `python -m ftqc_audit injection bb-uniform`. Add `--faults modes.jsonl --decode` to inspect accepted elementary faults after BP+OSD. The standard finite checker uses Stim for the physical probes. Decoder and sampling support additionally require NumPy, stimbposd and ldpc.

The reference checkout is kunliu7/qldpc-msi at 09b005ef4f2038baaa5fea565468ec18044567e7. To regenerate the examples, install that checkout in an isolated Python 3.12 environment with its declared dependencies. Run the following from BootLoops using that interpreter.

```sh
env -u PYTHONPATH python tools/ftqc_extractors/import_injection.py /path/to/qldpc-msi tools/ftqc_extractors/examples --circuits tools/ftqc_extractors/examples/injection-circuits
PYTHONPATH=tools python -m ftqc_audit.injection_sampling bb-uniform --seed 260405126 --shots 2000000 -o bb-samples.json
```

The BB uniform correlation driver uses two syndrome rounds. The scalability and asymmetric studies use three. The first two rounds supply the postselection decisions in each case. Additional detectors, including the terminal syndrome, feed BP+OSD. The conditional histogram includes the zero mask. A marginal sums every bin containing that output; total batch error divided by the output count is a different statistic. Seeds, software versions, batch size and shot counts travel with each banked sample. Exact repeatability of a Stim sample also depends on the processor's SIMD implementation. Wilson intervals describe Monte Carlo uncertainty for each displayed proportion.

| Experiment | Physical rate | First-order total | First-order correlated |
|---|---:|---:|---:|
| BB144 four Y states, uniform | 0.001 | 94p/15 | p/15 |
| BB144 four Y states, asymmetric | 0.001 on two-qubit operations | 8p/3 | p/15 |
| HGP225 six Y states, released product-coloration circuit | 0.0001 | 48p | p |

Coefficients count single faults silent on all detectors. Reported discard proxies use the product of independent detector-error-mechanism no-event probabilities. Resource estimates count the physical footprint and the two-qubit layers for a full attempt, with no early-abort credit. The HGP export retains the upstream single-qubit channels and omits explicit raw preparation flips and depolarizing two-qubit channels. Transfer to a different hardware schedule requires a separate physical audit. [Source inputs and counterexample](GUIDE.md#in-situ-simultaneous-injection)

For larger fixed carrier and Bell layouts, `python -m ftqc_extractors.injection_optimize SOURCE --objective protected -o proposed.json` uses SciPy's MILP solver for the off-support preparation bases. Use `--objective fixed` for the competing objective. The resulting source undergoes the same direct isometry and fixed-check tests; its postselection includes every fixed check. The solver optimum refers only to that layout. The separately replayable finite `injection_search` also varies carrier sites and pairings.

The [general BB144 source](examples/bb144-injection-general.json) injects all twelve carriers and the [general HGP225 source](examples/hgp225-injection-general.json) injects all nine. Both replay the carrier Clifford and projection frame. The published four- and six-output noisy benchmarks use the selected preparations instead.

## Seeded study on bish.local

With the published two-round BB circuit at uniform p = 0.001, 158,809 of 2,000,000 attempts passed. There were 1,029 wrong batches including 15 with multiple wrong outputs. The conditional batch error was 0.006479 (individual 95% Wilson interval 0.006097–0.006886), compared with Fig. 2's 0.00647. Dividing the batch rate by four gives 0.001620. The two-round discard fraction was 0.92060. The separate three-round run discarded 0.93194, consistent with the published roughly 93% figure. Under asymmetric noise, 95,181 of 500,000 three-round attempts passed and 251 had a wrong output. The batch rate divided by four was 0.000659. For six HGP outputs at p = 0.0001, 35,929 of 100,000 attempts passed and 184 of those batches contained an error. The full histograms include each marginal and the correlated error mass; empirical estimates do not certify a noisy T source.
