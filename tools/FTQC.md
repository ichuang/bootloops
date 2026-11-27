# Fault-tolerant quantum computing in BootLoops

The six tools on `ftqc_integration` cover different parts of a fault-tolerance argument. Start with the object you want to certify and identify what the next step needs from its receipt.

| Tool | Input and question | Evidence |
|---|---|---|
| [ftqc_codes](ftqc_codes/GUIDE.md) | Given CSS check matrices, which logical errors can the memory detect? | Code parameters, bounds on X and Z distance, witnesses and hierarchical decoding guarantees. |
| [ftqc_extractors](ftqc_extractors/GUIDE.md) | Does an HGP or gross-code extractor measure the intended logical observable and preserve the spectators? | Static and scheduled distance bounds, logical readout and frame checks, and logical T gate receipts with a specified resource input. |
| [ftqc_audit](ftqc_audit/GUIDE.md) | Which physical faults can evade a Clifford circuit’s detectors? | Detector and observable flips per location, replayable fault sets, circuit-distance intervals and preparation checks. |
| [ftqc_gadgets](ftqc_gadgets/GUIDE.md) | What does a decoded Clifford gadget do to an unknown encoded input with incoming errors and internal faults? | Residual-error bounds, refinement receipts for component substitution and two-level checks on concatenated codes. |
| [ftqc_zx](ftqc_zx/GUIDE.md) | Can a candidate replace a reference circuit under the specified noise model? | Pauli-web equivalence checks, ZX rewrite proofs, flag decoder tables and stochastic gadget contracts. |
| [ftqc_compose](ftqc_compose/GUIDE.md) | Do the output and fault guarantees of one operation meet the next operation’s input requirements? | Typed joins, residual data and frame checks, level reduction and local-stochastic error bounds for a composed study. |

## Where the tools connect

- **From a code to a noisy circuit.** The audit, gadget and extractor packages reuse code constructions or binary linear algebra from `ftqc_codes`. A static distance counts data qubits supporting a logical Pauli. To count faulty circuit locations, the study also specifies preparations, gates, idles, measurements, detectors and protected outputs.
- **From extractor to audit to composition.** The extractor defines the measured logical and its spectators, and exports physical X and Z probes as native audit circuits. The compose gross-code checker replays the extractor distance receipt and uses the audit fault model to calculate the residual data Pauli after applying the split-record frame. The HGP bridge also supports replay in the enclosing circuit and export of verified-cat preparations. [HGP bridge](ftqc_extractors/AUDIT.md) · [Gross-code probes](ftqc_extractors/GROSS_CIRCUIT.md) · [Residual checker](ftqc_compose/GUIDE.md#gross-code-residual-data-and-frame)
- **From a Clifford schedule to a concatenated gadget.** `ftqc_gadgets` uses the audit fault propagator and adds encoded Bell inputs, reference qubits and decoders. Its two-level checker binds inner receipts to signed logical operations and outer faults to physical columns. It verifies the decoder radius and records any specified ancillas and the rectangle argument assumed by the reduction. A purification tree that makes such an ancilla gets a factory certificate, which the reduction accepts in place of the assumption. For a program with sparse corrections it computes the fault cone of every decoding event, at one or two levels. [Gadget guide](ftqc_gadgets/GUIDE.md)
- **From gadgets to a circuit-wide noise bound.** `ftqc_compose` checks ordered wires, residual filters and bad fault supports. Its level-reduction checker substitutes physical supports into a logical circuit. A local-stochastic bound also needs a joint fault-support model. A two-level gadget receipt can inform the local claim; a Definition 4.8 hypothesis supplied to compose remains an explicit premise. [Composition guide](ftqc_compose/GUIDE.md#level-reduction)
- **Comparing a ZX design with a location-fault study.** Run the ZX verifier on the proposed circuit and its reference. To use its result in another argument, match the gates and idles, input codespace, flag or detector handling, output frame and fault-counting convention. The ZX contract has its own verifier and noise model. [ZX guide](ftqc_zx/GUIDE.md)

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

## QR-47 as an outer code

The [[47,1,11]] quadratic-residue code supplies a candidate decoding radius of five. With lookup tables disabled, the gadget verifier proves that radius from the actual decoder check matrix and checks that each fault stays in its assigned physical column.

```sh
python3 -m ftqc_codes certify qr:47 -o "$out/qr47.json"
python3 -m ftqc_gadgets analyze steane_ec:qr:47,anc=spec --local --code-receipt "$out/qr47.json" --no-tables -o "$out/outer.json"
python3 -m ftqc_gadgets verify "$out/outer.json"
```

The report states that the gate, error-correction and any-input properties hold when faults and input errors touch at most five coordinates. Read with a coordinate as a column of inner blocks, that is the outer half of a two-level argument; `ftqc_gadgets/levels.py` joins it to the inner receipts.

That gadget takes its two ancillas as specifications. A purification tree that makes one is certified separately. One level up an idle is an inner correction, so a raw outer preparation ends with idle rounds, and only the last round counts as local.

```sh
python3 -m ftqc_gadgets analyze prep:qr:47,zero,tree=XZXZ,rounds=2 --tree --strict --light --wmin 11 -o "$out/zero.json"
python3 -m ftqc_gadgets verify "$out/zero.json"
```

The report says the sixteen-copy tree acts as its specification, with wire faults only in columns holding a failed local or check location, unless four leaves each hold a failed junk-capable location. `levels.reduction(..., factories=...)` then lists no assumed ancilla and records that count. The probability of the junk event is not computed.

A whole program, with corrections placed where its schedule puts them, gets a cone receipt: for every decoding event, the sources whose faults its decoder can see.

```sh
cat > "$out/prog.json" <<'JSON'
{"inputs": ["A", "B"], "outputs": ["A"],
 "layers": [[["CX", "A", "B"]], [["ECa", "A"]], [["ECb", "A"]], [["MZ", "B", "m"]]]}
JSON
python3 -m ftqc_gadgets cones "$out/prog.json" -o "$out/cones.json"
python3 -m ftqc_gadgets verify "$out/cones.json"
```

The receipt lists each event's cone weight and ancilla atoms, and reports that linear propagation under several settings of the classical controls stays inside every cone. It gives no probability; the weights are what a tail bound takes as input.

## Carry a distillation bound into a logical T gate

The research study feeds the correcting 15-to-1 factory into an explicit CSS encoder and the extractor’s HGP logical T gate. Compose replays the gate receipt and verifies the four adaptive measurement branches. For the banked correlated factory model, ideal encoding and gate execution give a logical-channel half-diamond distance of at most 21/2000. A single encoder fault can produce an accepted logical Z, so the noisy study accounts for factory, encoder and gate faults separately.

Replay the eight-qubit example from the **research repository root**. The companion 450-qubit example and both gross-code residual studies run as part of the full replay. [Composition studies](https://github.mit.edu/ichuang/RESEARCH-extractors-proof-system-2026/tree/main/studies/compose)

```sh
PYTHONPATH=bootloops/tools python3 -m ftqc_compose verify studies/compose/distilled-t-hgp8.json studies/compose/distilled-t-hgp8-receipt.json
```

## Reproduce the integration checks

From BootLoops, run the six registered suites together. Install Stim and put CaDiCaL on `PATH` (or set `FTQC_CODES_CADICAL`) to exercise both external audit checks. The optional Lean build has separate prerequisites and an opt-in command in the [code guide](ftqc_codes/GUIDE.md).

```sh
python3 run_selftests.py ftqc_codes ftqc_audit ftqc_extractors ftqc_gadgets ftqc_zx ftqc_compose
```

From the research repository root with its pinned BootLoops submodule initialized, run the command below to verify the saved inputs, evidence hashes, bounds and witnesses. The twelve-round search also requires a C++17 compiler with unsigned 128-bit integers. [Research studies and replay instructions](https://github.mit.edu/ichuang/RESEARCH-extractors-proof-system-2026/tree/main/studies)

```sh
sh studies/replay.sh
```
