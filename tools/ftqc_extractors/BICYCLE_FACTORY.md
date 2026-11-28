# In-situ T sources and the BB144 pivot factory

The `injection_export` Extractors study has three models. `closed-Y` computes
the exact joint sign law for selected channels of a released stabilizer
experiment. `sampled-Y` recomputes statistics from supplied counts, binding
the circuit and postselection. `raw-T` starts from explicit T carriers and
replays a signed CSS injection isometry. Its independent input events may
act jointly on any listed carriers and may herald rejection. Encoding,
syndrome projection and frame correction are ideal in this model. A raw T
law for noisy open-output encoding remains an open extension.

An `insitu-rm15` input embeds each source spec, receipt and list of
`[member, column]` assignments. The Magic and Compose CLIs replay the same
joined receipt. All 15 zero-based columns must be assigned once. Source
masses include the original acceptance event. Distinct source attempts can
be multiplied under the explicit independent-noise premise. A zero-yield
factory produces no conditional output law. `insitu-retries` prices a
synchronized group with explicit lane durations and a realized transcript.
`magic-logical-t` accepts the replayed output under ideal input decoding
and transfer into the reference factory and charges its HGP output encoder.

The Extractors `bicycle_rm15_map` study gives the signed BB basis and all 15
physical joint parities. In this basis eight use the spare in logical zero.
Each parity is a coordinate translation of the one-round Gross Z
merge/split, including its public sign and physical split correction. The
translation preserves both stabilizer groups. The six-coordinate layout is
given in `bicycle.layout`; these are our port choices. Rotation signs use
positive T phase on odd parity, or exponent −iπP/8 up to a global phase.
The final S brings the surviving output to the positive T convention.

For a `bb144-rm15` Compose study, embed the source study in `source` and
specify nonnegative integer feed-forward `latency`. The compiler places
each CSS input block outside the live factory. It checks the source frame,
prepares the pivot, measures the source-pivot ZZ and Gross joint parity,
reads the source in X, chooses pivot X or Y, applies the physical byproduct
and returns the translated memory to its original coordinates. Each Tanner
generator is two simultaneous layers of SWAPs expanded into six CNOT
ticks. Public readout bits index physical conditional arms in little-endian
order. Each arm reserves the longest arm's window. Supplied `couplers` and
`schedule` can be checked against the installation.

The expanded schedule contains a template dictionary, per-step coordinates,
branch conditions, record names, readout parities and detector definitions.
The preparation primitive T initializes a raw carrier. All other generated
operations are Clifford operations, resets or measurements. Peak qubits
are reserved physical addresses. The installation uses a physical parity
meter with additional fanout for external input, pivot and final checks.
This implementation prices only the declared geometry and timing. The
paper's 378-qubit gross pivot uses a different LPU and adapter. Reproducing
that inventory and its logical operation error model is separate work.
[Xu et al., Fig. 4 and Tables III–IV](https://arxiv.org/html/2602.20546v2)

A `bb144-rm15-price` study adds independent X faults immediately before
specific physical source-X meter readouts. Replay traces each bit flip
through the pivot byproduct, final acceptance and output. Optional seeded
sampling reports accepted and wrong counts and individual Wilson intervals.
Other physical locations, idles and classical decoding are ideal. The
arbitrary-fault circuit distance is open. Complete physical Clifford probes
replace every raw T by +i and compare the final signed output against a
separate amplitude calculation while preserving a Bell spectator.

The tests `ftqc_extractors.test_bicycle`, `ftqc_compose.test_insitu` and
`ftqc_compose.test_bicycle_install` build and replay these studies:

```sh
PYTHONPATH=tools python3 -m unittest ftqc_extractors.test_bicycle ftqc_compose.test_insitu ftqc_compose.test_bicycle_install
```
