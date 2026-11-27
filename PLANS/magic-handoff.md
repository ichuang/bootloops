# Delivering verified magic states to encoded operations

**First deliverable:** consume one accepted T state from Magic's 15-to-1 factory, encode it into a fresh resource block of the eight-qubit HGP example, and replay its use on a logical data qubit. The resulting Compose receipt should identify the factory attempt, the physical transfer and the online measurement branches, and bound the error under their stated assumptions. Follow with an occupied high-rate block and then correlated outputs.

## Starting point at the Magic merge, `827e0b9`

| Owner | What exists | Work at the handoff |
|---|---|---|
| Magic | `ftqc_magic-receipt-1`, `MagicSourceContract`, factory replay, pipeline ledgers, correlation modes and bridge summaries. | Expose the verified quantum target, output grouping, error model and conditioning in a form the consumer can use. [Receipt verifier](../tools/ftqc_magic/certify.py) · [Source contract](../tools/ftqc_magic/ir.py) |
| Gadgets | Clifford fault propagation, encoded Bell-input checks, signed measurements, named surviving pairs, and incoming/internal residual profiles. | Certify the **transfer** of an arbitrary input qubit to a specified encoded slot, including spectators, frame and noisy locations. Magic retains ownership of non-Clifford preparation. [Gadgets guide](../tools/ftqc_gadgets/GUIDE.md) · [Surgery instruments](../tools/ftqc_gadgets/SURGERY.md) |
| Extractors | An HGP logical-T instrument with four feed-forward branches, installed measurement schedules and PBC compilation. PBC currently requests a fresh independent encoded T at a cache block. | Bind that request to a verified delivery, a named slot, a readiness event and a single use. [PBC compiler](../tools/ftqc_extractors/pbc.py) · [Installed surgery](../tools/ftqc_extractors/installed_surgery.py) |
| Compose | The existing `distilled-logical-t` checker replays a correcting `distill15` proof, checks a CSS encoder and verifies the T branches. An injection-batch checker also follows shared faults into two ideal S consumers. | Add a Magic-backed source alternative, reuse the existing encoder and branch checks, and join the three native receipts. [Logical-T checker](../tools/ftqc_compose/logical_t.py) · [Batch consumer](../tools/ftqc_compose/injection_batch.py) |

## 1. Make a Magic receipt consumable

Add one Magic-owned reader taking the factory spec and receipt. Call `verify_receipt` with its normal strict source check. Return typed evidence derived from the replay. Today the verifier compares the kind, contract, optimality gap and exports. The optional `details` payload is outside that comparison. Any newly consumed phase correction, joint mask law or physical schedule from that payload needs its own replay or inclusion among verified fields. Keep existing v1 inputs usable and explicitly version new serialized evidence. [Verifier](../tools/ftqc_magic/certify.py) · [Current exports](../tools/ftqc_magic/bridges/extractors_bridge.py)

Bind the normalized phase target, ordered wires per state, number of states per batch, sign or Clifford-frame convention, acceptance rule, denominator, noise premises and cost units. A concrete regression is the catalog 15-to-1 T factory: it returns the polynomial label `0`, meaning T on wire zero, with phases `(0,1)` modulo 8. Comparing a display string to `"T"` would reject it. Conversely, two CCZ states have six resource wires even though the qLDPC factory reports `k_outputs=2`. [Target representation](../tools/ftqc_magic/ir.py) · [qLDPC producer](../tools/ftqc_magic/qldpc_factory.py)

Admit the binary 15-to-1 factory first. Its exact evaluator supplies joint masses for **acceptance and** a nonzero Z mask. Recover the clean mass from exact acceptance and normalize once. Mark exact laws, upper bounds and unavailable laws separately. For a fresh T state with a proved Z-error mixture the conditional Z-error probability is its trace distance from the ideal state. For a general pure-target infidelity bound ε, use the conservative trace-distance bound √ε with outward rounding. A general state bound alone supplies no Pauli fault pattern for Gadgets. Preserve independence from the data, reference and subsequent noise as an explicit premise. [Binary evaluator](../tools/ftqc_magic/matrix.py) · [Current scalar bridge](../tools/ftqc_magic/bridges/compose_bridge.py) · [Fuchs–van de Graaf inequalities](https://arxiv.org/abs/quant-ph/9712042)

## 2. Land one end-to-end T and its negative controls

Factor the existing CSS encoder and its check into a reusable Gadgets transfer verifier. Start with a fresh resource block, ordered logical X/Z representatives, specified spectator initializations, and a discarded input wire. Check the logical isometry with an entangled reference. Bind the actual preparations, couplings, measurements, locations and frame to its receipt. For a source already encoded in the requested block, allow direct delivery only after checking code, slot, wire order and frame. Charging each transfer location once also matters when factory and consumer share hardware. [Existing encoder and single-fault witness](../tools/ftqc_compose/logical_t.py)

Add an optional source binding to the existing T consumer. Compose should replay the source, transfer and extractor receipts and reuse the existing corrected Kraus-branch check. Report factory acceptance per attempt, conditional resource error, transfer fault contribution, and online joint acceptance-and-failure separately. For a conditional final error claim, require a positive lower bound on the **same final acceptance event** before dividing. Keep the current legacy source path as a regression control.

**Exit test.** Run `bk_15_1_3_t` in correcting and postselected modes at input Z-error rate 1/1000 with ideal factory Clifford checks. Correcting mode accepts every attempt and has 105 harmful input pairs. Postselection has a distinct acceptance probability and distance three. Compare correcting mode to the existing Compose upper bound under their shared assumptions, and independently enumerate the 32,768 input masks to check both factory laws. Replay all four online branches with clean and Z-flipped resources. Turn on encoder noise and reproduce the single encoder fault that causes an accepted logical Z. Reject a T† or CCZ substitution, swapped logical pair, altered source receipt, and unproved frame correction. A known frame difference can pass with its verified correction.

## 3. Install the handoff on a high-rate memory

At the PBC request, bind the source instance, attempt, batch member, resource wires, destination block/slot, acceptance bit and ready tick. A rejected attempt frees its wires for retry and produces no consumable slot. A successful output can be consumed once along any realized execution path. Track retries, buffer occupancy, shared couplers, live memory idles and decoder/feed-forward latency. Magic already reports expected wire-ticks per output with factory retries included where applicable. Record expected cost, peak physical occupancy and elapsed ticks separately, with explicit units and retry assumptions.

Use a tiny addressable lifted-product block as the next fixture. Deliver by a checked swap or teleportation into a free slot while a different logical qubit holds unknown data. The fresh-block encoder initializes its spectators and is unsuitable for that operation. Verify the spectator using an untouched reference, then execute a logical T and a follow-on measurement on a known test input with a phase-sensitive readout. Test a rejected attempt followed by success, early consumption, double consumption, and collision with an installed surgery. Graduate to BB144 or HGP450 after the tiny program passes. A per-setting ideal boundary leaves the corresponding whole-program fault claim open until the enclosing physical schedule is replayed. [High-rate roadmap](high-rate-qec.md) · [Surgery milestone status](extractor-surgery.md#implementation-status-at-6613a17)

## 4. Extend to batches and CCZ

Start with two T outputs and two consumers. Preserve the shared acceptance event, ordered error masks and common fault identities. Prove that a deliberately shared fault can affect both consumers and compare the composed distribution with direct enumeration. For sources offering bounds only, propagate a justified joint bound or leave the numerical claim open. Interleaving requires evidence for distinct factory instances and the stated independence of their noise. For the existing `joint_3twirl` pipeline option, verify transfer of the full downstream distribution before using it to price correlated failures. A normalization check alone does not establish that transfer. [Pipeline](../tools/ftqc_magic/pipeline.py) · [Existing correlated-consumer control](../tools/ftqc_compose/test_injection_batch.py)

Next implement the three-wire CCZ consumption instrument and its branch-dependent Clifford corrections. Exercise the paired qLDPC CCZ outputs with separate resource ports and one batch identity. Catalyzed streams require an explicit catalyst lifetime and shared history before reuse. Keep the first catalysis experiment terminal. A physical factory-to-memory claim additionally needs the actual preparation and transfer schedules: the current Magic fault bridge exercises a Clifford frame probe. [qLDPC factory](../tools/ftqc_magic/qldpc_factory.py) · [Fault bridge](../tools/ftqc_magic/bridges/faults_bridge.py)

## Acceptance and rollout

Land stages 1 and 2 together as the minimum usable integration. Place positive, negative and inconclusive cases in the owning package suites. Bank source specs, native receipts, installed schedules, noise and retry premises, and a replay command in the research study. Run `python3 run_selftests.py ftqc_magic ftqc_gadgets ftqc_extractors ftqc_compose ftqc_faults` from BootLoops, then replay the existing research studies at the candidate pin. Stages 3 and 4 can then extend the same consumer without changing how legacy factories or conditional T fixtures are verified.

## Implementation status

The `magic-logical-t`, `magic-batch` and `magic-interleaved` Compose studies
and the PBC `magic_delivery` installation implement the four milestones
above. Gadgets exposes `magic-transfer` and native receipt replay. The
correcting and postselected RM15 laws, BH14 shared-output law, paired qLDPC
CCZ ports, terminal C2T, the installed four-qubit lifted product, and the
HGP450 reference-preserving SWAP have executable controls. The research
repository banks inputs, receipts and replay instructions in
`studies/magic-handoff/`.

The largest block check covers transfer. The installed tiny program has
per-setting evidence at ideal boundaries. Enclosing whole-program noisy
distance, noisy factory-to-memory scheduling, reusable catalyst histories,
arbitrary correlated distillation cascades, and a conditional final logical
error bound remain open. The joint CCZ twirl is verified for one or two CCZ
states followed by the ideal consumption instrument. See the
[Compose guide](../tools/ftqc_compose/GUIDE.md) and
[Extractors guide](../tools/ftqc_extractors/GUIDE.md).
