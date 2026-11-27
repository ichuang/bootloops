# From in-situ injection to a 15-to-1 bicycle factory

**First deliverable.** feed a replayed, correlated batch of encoded T inputs from `ftqc_extractors` into the postselected `bk_15_1_3_t` evaluator in `ftqc_magic`. **Second deliverable.** execute that same factory on named logical ports of a BB144 gross block and verify its acceptance, output frame and cost. The decisive end-to-end test is a shared input fault that survives postselection. Its identity and probability must survive both integrations.

At BootLoops `77de822`, the [existing handoff](magic-handoff.md) delivers a Magic output to an encoded consumer. The reverse path has two useful starting points. Extractors already verifies signed CSS injection isometries, and Compose already enumerates a restricted BB144 noise model with a shared acceptance decision. Its noisy observables are logical Y probes. Magic already knows the RM15 factory matrix, diagonal Clifford corrections and check-recycling schedules. The milestones below are proposed work. [Liu et al., in-situ injection](https://arxiv.org/abs/2604.05126v1) · [Xu et al., bicycle factories](https://arxiv.org/abs/2602.20546v2)

## Integration 1. In-situ injection as a Magic input

### 1. Export an accepted, ordered batch from Extractors

Factor the physical-channel enumeration in `ftqc_compose/injection_batch.py` into a producer-side routine alongside `ftqc_extractors/injection_circuit.py`. Preserve the existing Y-to-S consumer as a regression. Give the export a versioned schema and a verifier that rebuilds it from the injection spec, circuit, noise locations, decoder or frame rule, and acceptance predicate. Identify the preparation attempt, block, ordered logical slots, signed basis, state target, and any spectator initialization assumptions.

For an exact Z-mask model export rational masses $w(m)=P(A,m)$, including the clean mask. Their sum is source acceptance. Specify which faults are shared between outputs and which locations are ideal. Give sampled histograms their own evidence type with attempted and accepted shot counts, decoder configuration and uncertainty. A state bound can support a conservative bound when a consumer has a valid propagation rule. Only an exact law enters the exact evaluator.

**Exit tests.** Replay the restricted three-channel BB144 Y study and recover the conditional two-consumer Z probability `1/14988`. Replay the two-round BB144 sample summary with 158,809 accepted attempts out of 2,000,000 and 15 accepted batches with multiple wrong outputs. Mutate the frame, slot order or source digest and require replay to detect the change. Exercise zero acceptance and exhaustion of the enumeration budget. Preserve a zero-yield or inconclusive result as appropriate.

### 2. Establish a T source and propagate its joint law through RM15

Start with a small CSS isometry already proved for arbitrary carrier inputs. Prepare named $T|+\rangle$ carriers, declare ideal encoding and introduce a finite set of known pre-encoding Z faults, including one mechanism acting on multiple carriers. The isometry transports these faults to named logical Z operators. This supplies an exact T fixture with explicit noise assumptions. Extending it to noisy syndrome extraction requires a replay of the accepted open-output channel or a verified twirl and its noise assumptions. The banked Y measurements continue to test the decoder and correlations.

In `ftqc_magic`, add a factory-input alternative that binds ordered batch members to rotation columns. Sum their **joint** masses into the existing syndrome and logical-label calculation. Permit multiplication across attempts only with distinct provenance and an explicit independent-noise premise. For a finite attempt group, report $P(A)$, $P(A\cap F)$, $P(A\cap F\cap E)$ and their supported conditional ratios, where $A$ means all required source batches accepted, $F$ means the factory accepted, and $E$ means the output is wrong. Report retries under a named policy separately. Compose owns cross-tool replay so the Magic evaluator can take a verified law without importing the Extractors CLI.

**Exit tests.** Independently enumerate all 32,768 masks at input rate 1/1000 and recover the existing postselected RM15 acceptance and conditional output error. For the correlated control, let one source event of conditional probability $q$ flip precisely the three columns with syndrome labels 1, 2 and 3. Their syndrome XOR is zero and the output parity is odd. With the other inputs and factory operations ideal, the factory must accept every delivered batch and give conditional output error exactly $q$. Check that replacing those columns by independent inputs of the same marginal produces a different result. Reject reused members, a Y receipt in a T port, a sampled law requested as exact, and division by zero final acceptance.

### 3. Bank a replayable injection-to-factory study

Use `ftqc_compose` to join the native injection evidence, column assignment, factory spec and output receipt. Record how 15 distinct inputs are supplied, where accepted batches wait, which members are unused, and when a block may be reinitialized. Extend the Magic delivery reader to verify the new output evidence. Give the delivered output a named encoded port and frame. Count a shared batch preparation or retry once.

**Exit tests.** Replay a rejected batch followed by its replacement. Compare delivered-state yield and cost against a hand-calculated small retry example. Use the existing HGP8 logical-T consumer as an output check, including its Z-flipped resource branches, with any transfer explicitly charged. Bank the exact T example and the empirical Y comparison as separately labeled studies. This is the release point for integration 1.

## Integration 2. Compile postselected RM15 onto a bicycle block

### 4. Map signed factory operations to BB144

Have Magic export the ordered Pauli-product rotations, their signs and Clifford corrections, preparations, final checks, acceptance rule and surviving output. Begin with the unrecycled five-wire form. In Extractors, pin the BB144 check matrices, signed logical basis, logical processing unit (LPU) measurement menu, and allowed code automorphisms. Map each factory wire to a slot. To implement a requested Pauli $P$ by measuring $PZ_s$, prove that spare slot $s$ is in the specified Z eigenstate at that point and propagate its sign. Use Gadgets to verify the signed measurement and the action on surviving logicals. Record the automorphism, byproduct and physical setting for each rotation.

**Exit tests.** Reproduce the paper’s target of 15 native rotations after masking. A native rotation still consumes an automorphism sequence and an LPU measurement. Replay their signed action. Test a wrong automorphism, a reversed sign and masking on an unknown or occupied spectator. Check that untouched workspace logicals preserve a reference. A missing physical setting should identify exactly which rotation has only a logical implementation. [Xu et al., §§IV-A–IV-C and Fig. 5](https://arxiv.org/html/2602.20546v2)

### 5. Install and verify one complete factory attempt

Choose the paper’s **pivot-injection** variant for the first benchmark. Bind the 15 incoming resources, pivot, masks, data slots, measurements, resets, routing, decoder latency and feed-forward to a physical timeline. For an in-situ source in a separate BB block, include the transfer to the pivot. A colocated source needs a proof that preparation and reset preserve every live factory slot. Track the active logical permutation between rotations. Reuse the installed-surgery machinery in Extractors and the signed instrument verifier in Gadgets. Check the output after the final frame correction, including postselection. Add check recycling or a shorter automorphism order as a second schedule only after the baseline passes.

**Exit tests.** With ideal hardware, compare the compiled instrument’s acceptance and output Z action with the Magic oracle for all 32,768 input masks. Use exact branch checks or verified stabilizer reductions for adaptive outcomes, with an independent phase-sensitive output check for clean and wrong states. Replay the shared triple from milestone 2 on the mapped slots. Detect a missing byproduct, an early reset, a measurement before its source is ready and a coupler collision. Any optimized schedule must give the same accepted instrument as the baseline.

### 6. Price the installed factory and identify its physical error claim

Compose should report peak physical qubits, elapsed cycles per attempt, acceptance, error conditional on final acceptance, and expected qubit-cycles per accepted output. Include LPU and pivot ancillas, automorphisms, idles, buffering and source retries. Publish both a factory-only boundary and an inclusive source-to-output boundary. For variable attempt cost, accumulate costs across actual rejects and successful attempts under the stated retry policy.

**Exit tests.** With zero hardware noise, recover integration 1’s accepted output law. Inject a named physical fault and replay its detectors, acceptance and output frame. For a finite restricted noise model cross-check the whole-attempt result by direct enumeration. For a larger simulation report counts and intervals. Keep whole-circuit fault-distance claims open until an enclosing circuit proof covers them. As an external comparison, Xu et al.’s gross pivot row at $p_{\rm phys}=p_{\rm in}=10^{-4}$ reports 378 physical qubits and 5,999 timesteps. Match their input, timing and retry conventions before comparing volume or explaining a discrepancy. A distinct in-situ source run should use its own measured or proved inputs. [Xu et al., Tables III and IV](https://arxiv.org/html/2602.20546v2)

## Rollout and replay

Land the milestones in order on a feature branch, with one tested commit per milestone. Register tests in the owning package batteries. At each integration release, run the following from BootLoops, and run the existing injection and Magic handoff replay scripts from the research checkout. Bank the new specs, native receipts, source and slot maps, schedules, noise and retry assumptions, expected results and a single replay command under `studies/`. Update the research pin after replay against the merged integration revision.

```sh
python3 run_selftests.py ftqc_magic ftqc_extractors ftqc_gadgets ftqc_compose ftqc_faults
```

## Implementation and evidence

Milestones 1–6 have runnable studies. The first integration preserves the
shared triple's exact error probability through rejection, RM15 and an HGP8
T consumer. The second realizes all 15 mapped pivot parities in the Gross Z
orbit, checks all 32,768 source masks, and prices the installed adaptive
factory under a restricted physical readout-noise model. The signed basis,
adapter, timing and evidence boundaries are documented in the
[implementation guide](../tools/ftqc_extractors/BICYCLE_FACTORY.md).
The research archive includes a compressed expanded pulse schedule,
reference-entangled Clifford controls, exact results and seeded counts.
Whole-circuit noisy distance, a noisy T law from the paper's empirical Y
probes, and a reproduction of the 378-qubit paper hardware remain open.
