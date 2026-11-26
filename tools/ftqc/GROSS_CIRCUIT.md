# Circuit distance of the gross-code logical Z measurement

For the unprimed logical Z measurement of BB144, seven noisy rounds have exact circuit distance **7** and twelve noisy rounds have exact circuit distance **10** under the schedule below. Distance counts the minimum number of physical locations that can change the reported measurement or the protected logical output while all detectors remain silent. The merged code has static distance **12**. [Static construction and split frame](GROSS.md)

## Run the seven-round example

From the BootLoops root, use Python 3.10 or later.

```sh
export PYTHONPATH=tools
spec=tools/ftqc/examples/gross-circuit.json
work=$(mktemp -d)
python3 -m ftqc check "$spec" -o "$work/receipt.json"
python3 -m ftqc verify "$spec" "$work/receipt.json"
python3 -m ftqc schedule "$spec" --probe Z -o "$work/z.json"
python3 -m faultcert records "$work/z.json"
```

The input requires `schema: 1`, `kind: "gross_circuit"`, an integer `rounds` from 1 to 32, and a positive `target`. The example sets both integers to 7. Export the other basis with `--probe X`. Add `--stim` to the faultcert command to compare every distinct detector and observable signature with Stim. Add one `--fault LOCATION:LABEL` per physical fault to replay an explicit witness. Locations index the noisy operations in export order, starting at zero. Pauli letters follow operand order. A faulty measurement takes the label `flip`.

For the twelve-round result, the research study supplies an input with `spacetime: true` and two ten-location witnesses in `fault_witnesses`. Each witness is a list of [location, label] pairs. FTQC checks every witness by forward propagation through both probes. Generate that receipt with `--budget 2000000000`, or use `verify` to recover the budget from the saved receipt. The time-resolved search compiles a temporary C++17 executable, requiring GCC or Clang with unsigned 128-bit integers. The budget limits projected search nodes per probe. An unfinished exclusion gives a lower bound and can leave a target inconclusive.

## Which operations can fail?

The 319 wires comprise 144 memory qubits, 18 edge qubits and 157 measurement qubits. Prepare the edges in the X basis, repeat the following 14-tick cycle, then measure the edges in X to compute the split Pauli frame.

| Tick | Operations |
|---|---|
| 1 | Prepare 73 X-check meters in the X basis and 84 Z-check meters in the Z basis. |
| 2–5 | Couple the edge qubits to their port and gauge meters in the early, red, green and blue layers. |
| 6 | Couple each of the 12 data ports to its Z meter. |
| 7–13 | Apply the seven gross-memory CNOT layers. The first layer also couples the 18 incident X checks to their edges. |
| 14 | Read all 157 meters in their preparation bases. |

The compiler uses the memory order of Fig. 13(a), fixes the early gates and local coloring in Fig. 14 (bottom), and completes the remaining coloring of Fig. 9 deterministically. The preparation and readout barriers, placement of the data-port tick and resulting idle locations are choices of this implementation. The paper studies a pipelined schedule and estimates circuit distance from simulations. The bounds here apply to the exported ticks. The compiler checks disjointness within each tick. It pulls each meter observable backwards through the cycle to check the measured Pauli, and checks that every generator commuting with the measured checks is preserved. [Cross et al., §4.2 and Figs. 9–10, 13–14](https://arxiv.org/abs/2407.18393v4) · [Compiler](gross_cycle.py)

A noisy CNOT admits all 15 nonidentity two-qubit Paulis at unit cost. A preparation or idle admits any single-qubit Pauli, and a measurement fault flips its reported bit. Every waiting data qubit and every live waiting meter has an explicit idle. Edge preparation and split readout are noisy, with one idle per memory qubit in each of those ticks. The encoded input and the terminal measurements used to test the output memory are ideal. The split correction is a classical update computed from the recorded bits. All listed couplings are directly available.

## What the two probes protect

Both experiments use the same physical fault locations. The **Z probe** starts with encoded zeros, detects changes in the memory Z syndromes and repeated port records, and tests the logical Z readout and twelve final logical Z values. The **X probe** starts with encoded plus states, detects changes in the X-check records, enforces seven cycle constraints at the split, and tests the corrected final X syndrome and eleven surviving logical X values. The measured qubit's conjugate X is consumed by the intended measurement. A Pauli with both X and Z components has sector projections on subsets of its faulty locations. For this CSS circuit, a harmful undetected Pauli must therefore evade the checks in a sector where it changes a protected observable. [Probe definitions](gross_probes.py)

## How to audit the lower and upper bounds

1. Faultcert propagates every detector and logical parity backwards to each fault location. Each physical mode gets a binary vector of the bits it flips. Vector addition is XOR. Replaying a witness uses a separate forward propagator. Stim can independently compare the full set of single-fault vectors. [Record model](../faultcert/record_model.py)
2. Project the vectors onto the final syndrome, split consistency checks and protected output tests. Coalescing equal vectors and dropping earlier checks relaxes the search. An exhaustive search rules out harmful projected sums through weight six for seven rounds. At twelve rounds, the first harmful projected supports have size eight in the X probe and nine in the Z probe. [Projected search](../faultcert/projected.py)
3. For twelve rounds, enumerate eligible projected supports through weight nine. Each nonzero column has a finite set of possible time-resolved detector vectors from actual fault modes. The X search finds 288 eight-column and 7,780 nine-column candidates. The Z search finds 360 nine-column candidates. At this width, duplicate nonzero columns would cancel in pairs to give a harmful support below the projected minimum. An eight-column candidate can also take one fault with zero projection. The checker tests this possibility. An affine-span filter leaves 435 X and 308 Z candidates, and exact selection from the remaining finite domains excludes every one. Relaxing distinct-location constraints at this stage makes each exclusion valid for the physical model. Replay reruns this finite enumeration. [Support enumeration](../faultcert/projected_supports.cpp) · [Time-resolved checks](../faultcert/spacetime.py)
4. Bound a wrong readout with correct terminal data. Fewer faults than rounds leave at least one clean round. The checker verifies on each elementary fault signature that, in a later clean round, the 72 old Z syndromes and the port-readout parity equal its contributions to the terminal syndrome and target Z. A future fault cannot change an earlier record. Silence requires agreement between rounds. A wrong readout in the clean round therefore requires a harmful, syndrome-free prefix. Returning all output tests and the syndrome to their correct values then requires a second such contribution in the disjoint suffix. The readout-only lower bound is the smaller of the round count and twice the projected Z-probe distance bound. [Readout check and bound assembly](gross_noisy.py)
5. For seven rounds, flipping the same port meter's readout each round changes only the logical measurement. For twelve rounds the study supplies two distinct ten-location witnesses, one changing surviving X tests and the other changing surviving Z tests. Forward replay checks silence and the flipped observables. Matching upper and lower bounds give the exact distances above.
