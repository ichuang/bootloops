# Fault-tolerant quantum computing in BootLoops

Choose a tool by the claim you need to prove. A distance can count qubits in a code or faulty locations in a circuit. A gadget proof can instead ask how much error survives correction on an unknown input, or whether an implementation can replace a reference circuit. The five packages on the `ftqc_integration` branch make these questions executable in one checkout.

| Question | Tool | Evidence |
|---|---|---|
| What can the memory code protect? | [`ftqc_codes`](ftqc_codes/GUIDE.md) | CSS code parameters, logical Pauli witnesses, and bounds on code distance from the check matrices. |
| Which faults can a particular Clifford circuit hide? | [`ftqc_audit`](ftqc_audit/GUIDE.md) | Detector and observable flips per physical location, replayable fault sets, circuit-distance bounds and strict fault-tolerance checks for preparations. |
| What does an encoded gate or correction gadget guarantee for an unknown input? | [`ftqc_gadgets`](ftqc_gadgets/GUIDE.md) | Bounds on residual errors with input errors and internal faults counted separately. Refinement and composition receipts justify replacing a specified sub-gadget with an implementation. |
| Does an extractor measure the intended logical observable and preserve the other logical qubits? | [`ftqc_extractors`](ftqc_extractors/README.md) | HGP and gross-code measurement proofs, static and scheduled distance bounds, and conditional logical T proofs with explicit magic-state assumptions. |
| Can a proposed Clifford circuit replace a reference under a given fault model? | [`ftqc_zx`](ftqc_zx/GUIDE.md) | Pauli-web equivalence checks, ZX rewrite proofs, flag decoder tables and stochastic bounds for a specified noise model. |

## Where the tools connect

- **Code to circuit.** The audit, gadget and extractor packages import code constructions or binary linear algebra from `ftqc_codes`. The audit and extractor SAT backends also use its CNF and LRAT routines. A static distance counts the minimum weight of a logical Pauli. To obtain a circuit distance, specify preparations, gates, idle ticks, measurements, detectors and protected outputs.
- **Circuit to encoded gadget.** `ftqc_gadgets` uses the audit circuit and fault propagator. The gadget layer adds reference qubits for unknown encoded inputs, decoder rules, and tests of the ideal logical action. Its `fc:SPEC` input runs a circuit from the audit library as a preparation gadget. Composition requires a refinement receipt for each substituted part.
- **Extractor to audit.** The extractor package exports the physical X and Z probes of a gross-code measurement as native audit circuits. For HGP schedules its audit bridge can replay faults in the enclosing circuit or export the verified-cat preparation. The extractor defines which logical results and spectators to protect. [HGP bridge](ftqc_extractors/AUDIT.md) · [Gross-code probes and proof](ftqc_extractors/GROSS_CIRCUIT.md)
- **ZX design and verification.** `ftqc_zx` has its own circuit representation and verification engine. It can check a candidate against an ideal specification under a selected noise model. To compare a ZX result with a location-fault study, match the gates and idles, input codespace, detector or flag handling, output frame and fault-counting convention. Each tool's verifier consumes its own evidence format.

## Follow a gross-code measurement into the audit tool

Run these commands from the BootLoops root with Python 3.10 or later.

```sh
export PYTHONPATH="$PWD/tools"
out=$(mktemp -d)
spec=tools/ftqc_extractors/examples/gross-circuit.json
python3 -m ftqc_extractors schedule "$spec" --probe Z -o "$out/z-probe.json"
python3 -m ftqc_audit records "$out/z-probe.json" -o "$out/records.json"
python3 -m ftqc_extractors check "$spec" -o "$out/distance.json"
python3 -m ftqc_extractors verify "$spec" "$out/distance.json"
```

The audit report counts 25,048 noisy locations in this probe and checks faults through weight two. The extractor receipt combines both CSS probes, the time history of the logical readout, and a physical witness. It proves exact distance **7** for the seven-round measurement. Flipping one fixed port measurement in every round changes the reported logical Z and passes all detectors. The research checkout also supplies a twelve-round input and a replayable distance-**10** result for the same 14-tick cycle.

## Three questions about Steane-code operations

These commands use the same shell setup. The first pair checks the static code parameters. The second checks error correction and output closeness through order one for a Knill gadget whose ancillas are supplied as specifications. The last compares a five-ancilla X-syndrome extractor to its reference through three faults in the ZX tool's CSS circuit noise model.

```sh
python3 -m ftqc_codes certify steane -o "$out/steane.json"
python3 -m ftqc_codes verify "$out/steane.json" --full
python3 -m ftqc_gadgets analyze knill:steane,anc=spec --check ec:1 --close 1 -o "$out/gadget.json"
python3 -m ftqc_gadgets verify "$out/gadget.json" --full
python3 -m ftqc_zx equiv steane_5anc steane_spec --w 3 --noise-model circuit_css --code steane_7_1_3
```

## Reproduce the integration checks

From BootLoops, run the five registered suites together. Install Stim and put CaDiCaL on `PATH` (or set `FTQC_CODES_CADICAL`) to exercise both external audit checks. The optional Lean build has separate prerequisites and an opt-in command in the [code guide](ftqc_codes/GUIDE.md).

```sh
python3 run_selftests.py ftqc_codes ftqc_audit ftqc_extractors ftqc_gadgets ftqc_zx
```

From the research repository root with its pinned BootLoops submodule initialized, run the command below to verify the saved inputs, evidence hashes, bounds and witnesses. The twelve-round search also requires a C++17 compiler with unsigned 128-bit integers. [Research studies and replay instructions](https://github.mit.edu/ichuang/RESEARCH-extractors-proof-system-2026/tree/main/studies)

```sh
sh studies/replay.sh
```
