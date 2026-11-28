# ftqc_extractors — GUIDE

Overview of the FTQC tools: [`../FTQC.md`](../FTQC.md).

KIND: package (Python standard library for generation and replay; CaDiCaL
optional for SAT evidence; Stim optional for some injection studies; a C++17
compiler for the twelve-round gross-code search). Core: `algebra.py`,
`extractor.py`, `circuit.py`, `proof.py`. Gross code: `gross.py`,
`gross_cycle.py`, `gross_distance.py`, `gross_noisy.py`, `gross_probes.py`,
`css_cover.py`. Surgery and compilation: `fixed_graph.py`, `block_reading.py`,
`surgery_*.py`, `installed_surgery.py`, `addressable_surgery.py`, `pbc.py`,
`pbc_surgery.py`, `bicycle.py`, `staged.py`. Logical T and injection:
`t_*.py`, `magic*.py`, `injection*.py`, `import_injection.py`. Bridges:
`audit.py`, `audit_bridge.py`, `action_orbit.py`. Tests `test_*.py`,
`examples/`. Battery class: selftest.

## PURPOSE

Fault-tolerance claims about logical operations performed on encoded memory
by an apparatus attached to it. An *extractor system* couples auxiliary
qubits to a memory so that a choice of checks measures a logical Pauli; a full
extractor reaches every logical Pauli of the block, a partial one a stated
set (He, Cowtan, Williamson and Yoder, §§4.2–4.5). Given a code and a
measurement construction, a circuit schedule, a logical T protocol or an
injection source, the package writes a receipt with one of these statements:

- **static**: the deformed (merged) code keeps a target distance, and the
  ideal instrument measures the requested signed logical and preserves the
  surviving logicals;
- **phenomenological**: a lower bound from the checked finite hypotheses of a
  published theorem (Blue et al., Appendices E–F; Cowtan, He, Williamson and
  Yoder, Proposition 3.9; intermediate-time soundness);
- **circuit**: the number of faulty locations in an explicit schedule needed
  to change the readout or a protected logical while every detector stays
  silent, with an ideal input and terminal syndrome;
- **channel**: the conditional logical T of an ideal encoded magic state, or
  a signed injection isometry from raw carriers.

A failed target is reported with the witness that breaks it.

## SCOPE

CSS memories given by check matrices (HGP seeds, BB144, general CSS); signed
stabilizer codes for fixed-coupler graphs. Fault model: any nonidentity Pauli
after a gate (one location per gate), a flipped preparation or record, idles
only where scheduled. Boundaries: ideal encoded input and ideal terminal
syndrome unless the input says otherwise. Published theorems are applied
after their finite hypotheses are checked; they are not re-derived.

## USE-WHEN

- A lattice-surgery or extractor construction should be shown to measure the
  intended signed logical, keep its spectators and keep a distance.
- A noisy schedule for such a measurement needs an exact or bounded circuit
  distance with a physical witness.
- A logical T protocol, a Pauli-based compilation onto installed couplers, or
  an in-situ injection source needs its logical action checked branch by
  branch.

## NOT-FOR

- Code facts alone (`ftqc_codes`), decoders (`ftqc_decoders`), probability
  bounds (`ftqc_compose`).
- Generic circuits without a logical measurement: use `ftqc_faults`.
- Threshold estimates: the samplers here report empirical histograms only.

## INVOKE

From the BootLoops root with `export PYTHONPATH="$PWD/tools"`:

    python3 -m ftqc_extractors check SPEC.json -o receipt.json [--backend exact|sat] [--budget N] [--seconds S]
    python3 -m ftqc_extractors verify SPEC.json receipt.json
    python3 -m ftqc_extractors schedule SPEC.json --probe X|Z -o probe.json     # export to the audit engine
    python3 -m ftqc_extractors measurement gross-z.json -o measurement.json
    python3 -m ftqc_extractors export-surgery SPEC.json -o instrument.json
    python3 -m ftqc_extractors injection bb-uniform [--faults modes.jsonl --decode]
    python3 -m ftqc_extractors.audit_bridge compare|replay|export-cat SPEC.json ...
    python3 -m ftqc_extractors.surgery_bench --help

Exit codes of `check`: 0 proved, 1 counterexample, 2 invalid input, 3
inconclusive search. `verify` returns 0 for a faithfully replayed positive or
negative receipt.

## INPUTS

JSON with `schema: 1`, a `kind` and kind-specific fields; bit 0 is coordinate
0; HGP qubits are row-major in the two product blocks. Kinds:

| Kind | What the receipt establishes |
|---|---|
| `hgp` | n, k and base distance from two classical seed matrices |
| `column` | a single-column Z surgery on HGP memory meets a static distance target |
| `canonical` | phenomenological lower bound min(d, k², kλ) for the full cyclic HGP extractor, from checked hypotheses |
| `circuit` | logical Z readout, surviving logicals and undetected-fault distance of a listed schedule |
| `t_gate`, `t_bounds`, `t_source` | the conditional logical T with an ideal encoded magic state, per branch and readout arm; the in-situ preparation and its accepted counterexample |
| `gross_z` | the ideal unprimed logical Z measurement of BB144 and exact static distance 12 of memory and merged code |
| `gross_circuit` | circuit distance of the noisy BB144 logical Z measurement (exact 7 at seven rounds) |
| `block_reading`, `gross_block_reading` | transversal reading of logical parities across blocks, split frame, phenomenological bound |
| `fixed_graph`, `pbc` | routing, expansion and couplers of a fixed extractor; signed settings for a bridged Clifford+T program |
| `surgery_sequence`, `intermediate_reading`, `modular_expansion`, `surgery_zx` | sequence-level phenomenological bounds, intermediate-time soundness, hypergraph expansion, ZX comparison of measurement and surgery |
| `installed_surgery`, `staged_install` | placed Gadgets instruments on a coupler inventory, replayed on all data |
| `action_orbit` | reachability of signed logical settings under verified code actions |
| `injection`, `injection_search`, `injection_circuit`, `injection_export`, `bicycle_rm15_map` | signed injection sources, finite source search, closed Y probes, raw-T sources and the BB144 RM15 parity map |

Inputs are trusted JSON; do not run files from sources you do not trust.

## OUTPUTS

A receipt bound to the SHA-256 of the canonical input, with the status, the
lower and upper ends and their evidence (search counts, LRAT proof, structural
certificate), witnesses as `[location, mode]` pairs, and the stated model and
boundaries.

## GATES

The registered battery (sixteen `unittest` modules, listed in
`tools/BATTERIES.json`): 72 tests, 6 skipped, in 52 s on a four-core Linux
virtual machine without Stim or CaDiCaL. It exceeds the toolkit's 20 s
standard. The tests cover the algebra, receipts and their tampering, the
edge-extractor proof and its hook control, logical T branches, the audit
bridge, the gross-code static and circuit receipts, architecture and PBC
compilation, injection, installed surgery, layouts, Magic delivery, the
bicycle factory map, and the LP(8,4) translated-seed comparison. A banked
LRAT proof is re-checked even without CaDiCaL.

## FOOTGUNS

- Static distance counts qubits; circuit distance counts locations. One
  faulty CNOT can damage both operands at weight one.
- Exact distance needs both a lower bound and a matching witness. A search
  that exhausts `--budget` or `--seconds` leaves the claim open (exit 3).
- Include every routing gate and idle whose faults the claim should cover;
  unlisted locations are noiseless.
- Phenomenological and circuit receipts have different boundaries and are not
  interchangeable; the two-block BB144 reading is ≥ 12 phenomenologically and
  only in [3, 12] at circuit level.
- The basis label of a dual (X) experiment names the probe before the global
  Hadamard conjugation.

## KNOWN GAPS

- The cited theorems (Blue et al., Cowtan–He–Williamson–Yoder, Liu et al.) and
  the counting lemma of the verified-cat T protocol are dependencies, not
  re-derived or machine-checked.
- Whole-program noisy fault distance of a compiled PBC program, of parallel
  surgery compositions and of the BB144 RM15 factory is open.
- Hardware routing, storage and timing are inputs; the T protocol is not
  compiled onto a fixed coupler network.
- Decoding is postselection or a lookup over single faults unless
  `ftqc_decoders` is used.

The sections below document each input kind in detail.

## Claims
* Base distance counts nontrivial Pauli normalizer elements modulo stabilizers.
* Static surgery compares the deformed CSS code with the input code (Blue et al., arXiv:2606.03507v2, Definition 11). Definition 16 and Theorem 17 concern single-column Z surgery. Cap comparisons by the base distance when applying the reduction.
* The canonical full cyclic extractor in Appendix F has phenomenological lower bound min(d,k²,k λ). The inputs must establish every nonzero single-column case and use the canonical assembly and d-round protocol of Appendix E.1. k is the seed dimension.
* A circuit claim counts faulty scheduled locations with ideal input code and final syndrome boundary. It requires a checked logical measurement and preservation of the remaining logical algebra. The weight of a correlated two-qubit Pauli at one gate is one.
* A lower bound is obtained by ruling out every damaging undetected fault below the threshold. Exact distance also requires a witness at the threshold. A resource limit leaves an open obligation.

## Running a verification

From the BootLoops root:

    export PYTHONPATH=tools
    python3 -m ftqc_extractors check tools/ftqc_extractors/examples/column.json -o /tmp/column.json
    python3 -m ftqc_extractors verify tools/ftqc_extractors/examples/column.json /tmp/column.json
    python3 -m ftqc_extractors check tools/ftqc_extractors/examples/edge-extractor.json -o /tmp/edge.json
    python3 -m ftqc_extractors check tools/ftqc_extractors/examples/edge-hook.json -o /tmp/hook.json
    python3 -m ftqc_extractors verify tools/ftqc_extractors/examples/edge-hook.json /tmp/hook.json
    python3 run_selftests.py ftqc_extractors

Check returns 0 for a proved claim, 1 for a counterexample, 2 for invalid input, and 3 for an inconclusive search. Verify returns 0 for a valid positive or negative receipt. An exact search checks increasing weights. A failure with matching lower and upper endpoints establishes the exact distance. Increase the target by one to search for an upper witness.

Use --backend sat with a circuit to generate a proof or fault witness. Install a CaDiCaL with text LRAT support and set FTQC_CODES_CADICAL to the executable. Verification uses only Python and the sibling ftqc_codes package. The verifier reconstructs the CNF using backward propagation, then checks the LRAT proof, or replays the witness on the original schedule. The same site's Pauli alternatives are mutually exclusive. --seconds limits the external solver. --budget limits exhaustive combinations. Inconclusive receipts retain the reason. The default battery verifies a banked LRAT even when the solver is absent.

## Inputs and model

Each example is a complete JSON input. 'a' and 'b' are lists of binary strings, coordinate zero first. Both HGP blocks are retained, including the transpose-code sector.

* hgp takes a and b. Its d is the minimum classical distance over active seed and transpose sectors. For an empty logical quotient d is null.
* column adds column (zero-based), word (integer bit mask), vertices, edges, ports (one vertex per row coordinate), optional paths (one edge mask per seed check), and target. The active ports must be injective. The measured Z has support word in the indicated left-block column. Vertex, path and cycle checks are constructed explicitly. The subcode Z test applies up to the original HGP distance. For a larger target the tool searches the full merged CSS code. A failed reduction is lifted to a physical logical witness.
* canonical takes n, polynomial (cyclic check polynomial bit mask) and undirected edges on vertices 0 through n-1. Input labels define the identity single-column port map. It denotes the full extractor of Appendix F.2–F.3, including its canonical port priority (Eq. 64), path deformations and Appendix E.1 protocol. An optional assembly is the ordered edge list from full_graph and must match exactly. The result applies to that mathematical construction. Replay checks the reduced CSS matrices for every nonzero seed word, verifies that the examined cosets exhaust the surviving quotient, and independently computes the min-cut. The published theorem and HGP distance formula are mathematical dependencies. Replays do not formalize their proofs in Lean.
* circuit takes a, b, measured_z (physical bit mask), ops, readout, detectors and target. readout is a list of measurement-record indices whose XOR reports the requested logical Z. Each detector is a list of record indices with deterministic even parity. The ordered ops are a serial schedule. Qubits 0 through n-1 are arbitrary encoded input data. Larger labels are ancillas. RZ/RX prepare +Z/+X eigenstates, MZ/MX measure and retire an ancilla, CX takes control and target, H and I act on one live qubit. Z_IF takes a data/ancilla qubit and an earlier zero-based measurement index for classically controlled Z. Every ancilla lifetime needs a preparation and a terminating measurement. Include explicit I locations wherever idle faults should count.
* Every preparation can flip its eigenstate. Every measurement can flip its record. At each other location any nonidentity Pauli on its operands can follow the ideal operation, including on Z_IF's quantum operand. Witnesses are pairs [zero-based gate index, mode]. On a unitary, bits 0 and 1 select X and Z on its first operand, bits 2 and 3 on its second. A preparation or measurement uses mode 1. The input codespace and the terminal ideal data stabilizer syndrome define the boundaries. The checker also validates the ideal measured value and preservation of a spanning set of the surviving logical algebra. A damaging pattern evades all specified detectors and the terminal syndrome while changing the readout or protected data. This is a finite-block undetected-fault distance under those boundaries. Decoder thresholds, unlisted idle locations, and noisy terminal readout require their own models.

## Worked controls and limits

edge-extractor measures the two vertex checks Z(data)Z(edge) of the one-edge graph twice, retires the edge in X, and updates the Pauli frame. The product of the first two records measures logical Z on data 0 and 2. Its two parity detectors compare the vertex readings across rounds. edge-hook appends two identical CX gates on data 4 and 5. Their ideal action cancels. A Z on the target after the first gate spreads back through the second and evades the terminal syndrome. The smaller repeat-z example measures parity with a single meter; its banked LRAT exercises the checker on a cold clone.

The canonical benchmark is an implicit full architecture at phenomenological noise level. Gate-level receipts cover explicit schedules and their declared boundaries. Each receipt binds to a SHA-256 of the canonical input JSON. The exact checker redoes finite searches during replay. The SAT checker accepts the LRAT fragment supported by ftqc_codes. Search time and proof size grow rapidly with threshold and schedule size. Resource limits yield open obligations.

## Addressable logical T

See [T_PROOFS.md](T_PROOFS.md) for the exact channel identity, fault contract, hardware assumptions and structural proof. Use the same check and verify commands with examples/t-gate.json and examples/t-hook.json. The exact and SAT backends both handle t_gate. A gate receipt has three cases: the T output and the two arms of a subsequent signed logical X measurement. A counterexample means at least one case failed the requested lower bound.

Inputs specify the HGP seeds, the anticommuting logical_x and logical_z masks on each block, rounds, target and resource. The contract for t_gate and t_bounds is ideal-encoded-T. verification_rounds defaults to one. representatives optionally supplies a pair of equivalent physical masks per measurement round. An experimental trailing_identity CNOT sequence must be checked as identity on every involved qubit. The t_bounds kind omits target and returns checked lower and upper bounds per case. The t_source kind takes resource=in-situ-all-plus and verifies the accepted preparation and resulting raw-Z counterexample. The 225-qubit benchmark uses the published classical seed with our selected single-site port. BootLoops works from the input matrices without requiring the source repository.

    export PYTHONPATH=tools
    python3 -m ftqc_extractors check tools/ftqc_extractors/examples/t-gate.json -o /tmp/t.json
    python3 -m ftqc_extractors verify tools/ftqc_extractors/examples/t-gate.json /tmp/t.json

The physical adapter is a verified-cat parity meter with direct access to the listed sites. Implementing the protocol on a particular fixed coupler network requires compiling and checking its routing and storage locations.

## Cross-checks with ftqc_faults.audit

Compare the two fault propagators on an extractor or fixed T arm, replay a logical witness, or export a verified cat for a native strict-FT analysis. The [comparison guide](AUDIT.md) gives the commands and the input, acceptance and output assumptions for each claim.

## Bivariate bicycle partial extractor

Use input kind `gross_z` for the published BB144 unprimed Z setting. The exact backend computes and certifies both CSS distances. The measurement command exports the merged matrices, ideal readout and split frame. See the [gross-code guide](GROSS.md) for the input, coordinate conventions, proof format and scope.

Input kind `gross_circuit` compiles the noisy BB144 operation. Use `schedule --probe X` or `schedule --probe Z` to export a native ftqc_faults.audit circuit, and `check` or `verify` for its distance. [Noise model, optional input fields and proof search](GROSS_CIRCUIT.md)

## Full block reading

Run `check` and `verify` on [the two-block fixture](examples/block-reading.json). Matrices are arrays of nonnegative integer row masks, least significant bit first. Each independent pattern row selects blocks. Its port parities measure one logical per member of the memory logical basis. The receipt gives the data-only projections, surviving dimension, deformed matrices and metachecks. In each pattern row, decompose the measured auxiliary split syndrome in the recorded syndrome basis and apply the corresponding correction on every selected block (Z for Z reading, X for its dual). The static bound applies Lemma 3.8 after checking the input distance.

## Fixed-coupler measurement graphs

[Four-port example](examples/fixed-graph.json). Signed strings use wire zero first. A port of -1 excludes that physical site. The theorem covers all logical representatives supported on the included sites. The verifier checks relative expansion on all cuts, every even subset via generating pairs within each check’s allowed routing pool, a complete cycle basis, the input distance and the maximum degree on both sides of the physical check-data wiring. These are finite-instance degree counts. The selected logical produces signed deformed checks, the parity of all vertex outcomes with its sign, and split paths in support order. Edge qubits use the vertex-Z convention: prepare plus, measure X on split, and apply the local factor of the measured Pauli for each odd path. The phenomenological protocol uses target rounds before, during and after the merge.

## Padded hypergraph surgery

[Three-block sequence](examples/three-block-sequence.json) measures two independent patterns in order. Each auxiliary CSS code has a fourth term supplied by its complete sparse cycle basis. The three chain maps are **column** masks into the original memory. Metacheck masks refer to original Z checks followed by the operation’s vertex checks. The verifier builds the simultaneous mapping cone for the distance proof; the physical protocol performs the listed operations sequentially with one merge round apiece and the stated boundary rounds. It checks both chain squares, vertex non-overlap, sparsity, injectivity and independence of measured logical classes. The auxiliary 1-cosystole search checks X cocycles modulo X checks. A failed sufficient condition leaves the theorem unavailable; it is not a physical fault witness. The distance model is phenomenological.

The [partial transversal fixture](examples/partial-block-sequence.json) embeds the auxiliary [[4,2,2]] code into one factor of each of two eight-qubit memories. It measures two parities among eight logical qubits, leaves six encoded, and replays a distance-two lower bound with shared boundary padding.

## Two-block BB144 benchmark

Run check and verify on [Z reading](examples/gross-block-Z.json) or [X reading](examples/gross-block-X.json). Each measures twelve independent corresponding parities and leaves twelve logical qubits. Each adds 72 data qubits and 144 check meters to the two memories. Proposition 3.9 gives phenomenological distance at least 12 with twelve syndrome rounds at each boundary. The underlying BB144 distance exclusion runs independently in both sectors.

The circuit benchmark uses ideal encoded inputs, one noisy merge, noisy split and ideal terminal readout. Greedy disjoint CNOT layers have explicit idles on live wires. Separate X and Z probes protect the readout and surviving algebra. Use `schedule SPEC --probe X` or `--probe Z` to export them for ftqc_faults.audit. For either reading the two probes give circuit distance in [3,12]: all faults on one or two distinct locations are excluded, and a weight-12 fault on data is replayed. The basis label on the dual experiment denotes the probe before global Hadamard conjugation.

## Expansion and intermediate time

[Modular expansion](examples/modular-expansion.json) checks the exact minimum in Definition 7.2 before and after Cartesian path thickening, with ports on the first layer. It returns both minimizing cuts. [Intermediate reading](examples/intermediate-reading.json) applies Theorem 6.5. Breadth-first search on attainable Pauli syndromes computes quantum soundness from minimum physical Pauli weights, including Y. Each individual cone must preserve the memory target and its auxiliary 1-cosystolic distance times its declared merge rounds must reach that target. Both temporal boundaries are checked. These enumerations are intended for finite studies and return an open obligation at their budget.

## Homomorphic measurement and surgery

[ZX comparison](examples/homomorphic-surgery.json) compiles two physical circuits from the same chain map: an encoded ancilla with data-to-ancilla CNOTs and destructive readout, and a vertex-check measurement with X split on the auxiliary edges. The ideal ZX boundary and outcome map are compared. Vertex strings label cosets modulo the auxiliary X checks; the extra split outcomes drive the returned Z-frame columns. Each unencoded open-wire implementation has its own ZX circuit-noise search through weight one. These circuit bounds have different boundaries from the padded phenomenological theorem.

## Bridges, adapters and Pauli-based compilation

[Two-block Clifford+T example](examples/bridged-pbc.json) installs a matching between ports, its connecting cycle checks and a signed logical X/Z basis on each block. The compiler checks the combined graph and the degree of the union of installed couplers. Physical site and meter labels in the hardware inventory identify shared wires across settings. The same construction accepts different memories at the two ends of an adapter when their joined apparatus passes the checks.

Each block reserves its last logical slot. Workspace gates use [block, slot] operands. In-block Clifford gates propagate into later Pauli axes. A cross-block CX uses a minus-Y ancilla, a joint parity, an X read and record-dependent Pauli corrections. A T uses a fresh independent encoded magic state on its designated cache slot and a Clifford correction. Resource supply is an input assumption. The compiler tracks corrections into final signed Z measurements. Settings include the physical checks, split paths and sign of the reported parity. Cycle numbers serialize shared blocks and permit disjoint blocks to measure in parallel. By default all two-bit rotation outcomes are compiled. A supplied outcomes list selects a single trace; b is its low bit and m its high bit. Logical rotation identities and physical graph hypotheses underlie the receipt. There is no gate-noise estimate for this entire compiled computation.

The compiler waits for the resource X readout before a dependent block uses its updated frame. Rotations using disjoint data and cache blocks can share a cycle. A partial layout can use any nontrivial logical supported on its installed ports. Each compiled setting also checks its particular signed representative.

## In-situ simultaneous injection

The `injection` input supplies a CSS memory and ordered signed logical pairs [X mask, Z mask, X sign bit, Z sign bit]. General mode cleans a complete basis onto carriers and synthesizes the physical carrier Clifford. Selected mode specifies raw carriers, plus and zero sites, and positive Bell pairs. Replay checks both signed logical generators on arbitrary inputs, including inputs entangled with a reference. Each listed frame column corrects a minus syndrome on the corresponding independent check while commuting with the cleaned logicals. Nonselected logical qubits are unspecified. The [passing](examples/bell-injection.json) and [failing](examples/split-bell-injection.json) ten-qubit controls differ only in their Bell pairing. The failed pairing satisfies the overlap-class rule of Liu et al. Eq. (15), but logical X on the first carrier splits two pairs. [General signed example](examples/general-injection.json)

Use `injection_search` with memory, logicals, count and objective (`protected` or `fixed`) for a finite exhaustive search. The search traverses carrier choices, Bell orientations and free X/Z initializations up to the work budget. It groups Bell endpoints by their full logical-support signatures. Search returns a source that can be replayed independently. Larger published inputs supply a configuration directly. `fixed_checks` lists every deterministic row, while `postselect` selects the verified subset actually used by the protocol. BB144 selects 32 of 37. HGP225 sizes three and four leave site 2 at the simulator's initial zero. Reference inputs come from qldpc-msi commit 09b005ef4f2038baaa5fea565468ec18044567e7, BB fixture 491 and HGP fixtures 32, 45, 38, 35, 4 and 5.

### Closed physical Y probes

Use an `injection_circuit` input with case `bb-uniform` (two rounds, correlation driver), `bb-uniform-three` (three rounds), `bb-asymmetric` (three rounds) or `hgp-release` (three rounds). Recompile the packaged Stim schedules with the pinned reference importer. ftqc_faults.audit binds the first detector mask to the verified fixed checks, the next mask to all first-to-second round parities, and subsequent decoder checks to later rounds and the ideal terminal Pauli measurements. Ideal logical Y signs include the phase of iXZ. Fault responses and coefficients are evaluated on the exported physical channels. The reference HGP builder includes single-qubit depolarization and omits explicit raw reset flips and two-qubit depolarizing instructions. Its reported benchmark applies to that released circuit.

`python3 -m ftqc_extractors injection bb-uniform --faults modes.jsonl` exports every accepted physical location/mode, its exact probability, detector syndrome and pre-decoding logical mask. Add `--decode` in the study environment to apply BP+OSD and obtain the residual error mask. Operation indices count flattened non-annotation instructions other than TICK. The coefficient in the standard receipt counts first-order faults silent on all detectors, matching Appendix F. First-order coefficients are not finite-noise upper bounds.

## Signed surgery instruments

See [the surgery interface and replay examples](../ftqc_gadgets/SURGERY.md) for physical measurement receipts, installed feed-forward schedules and conditional residual laws.

### Installing an accepted Magic output in a PBC program

An optional `magic_delivery` on a PBC input embeds a factory spec and strict
receipt, realized attempts, physical transfer couplers, decoder/feed-forward
latency in ticks and one ordered request per T. Each attempt names its wire
inventory, output wire and half-open execution interval. A successful member
occupies its raw wire until transfer. Rejected attempts release their wires.
The same member may be used on mutually exclusive branches, once per branch.

Gadgets checks a three-parity-CNOT SWAP into the reserved scratch slot and
the identity on every other logical Pauli pair, including entangled
spectators. The ledger serializes the chosen measurements and charges
transfer, waiting, buffering and latency. Physical settings use their checked
durations, placements and memory idles. With abstract settings a measurement
uses one model tick. Peak reserved qubits and elapsed time describe the
realized transcript. Factory time windows and independent retry noise are
study assumptions. Whole-program noisy fault distance is open. See
`test_magic_delivery.py` for the installed [[4,2,2]] lifted product, and
`ftqc_gadgets/test_occupied_transfer.py` for the 32-logical HGP450 port check.

The [in-situ BB144 factory guide](../ftqc_extractors/BICYCLE_FACTORY.md)
covers injection batch replay, RM15 input assignment, signed native pivot
measurements, physical scheduling and full-attempt pricing.

## Code automorphisms

The `action_orbit` study takes `code`, a paired `basis`, a list of verified `actions`, and signed logical `seed` and `targets`. It returns a shortest generator word for each reachable target and physically replays its sign and stabilizer coset. A completed orbit certifies unreachability under the supplied generators. Exceeding the state budget leaves the result inconclusive. The [LP and Gross examples](../FTQC.md#code-automorphisms) compare measurement settings with installed implementations.

The [staged logical gate guide](../ftqc_gadgets/STAGED.md) describes the Steane/RM15 switch, pieceable CCZ, single-fault model and cross-tool replay.

CREDIT: extractor systems: Z. He, A. Cowtan, D. J. Williamson and T. J. Yoder
(arXiv:2503.10390); the canonical HGP extractor: Blue, He, Zhou and Chuang
(arXiv:2606.03507v2); gross-code partial extractor: Cross, He, Rall and Yoder
(arXiv:2407.18393v4); fast surgery and Proposition 3.9: Cowtan, He, Williamson
and Yoder (arXiv:2510.14895); simultaneous injection: Liu et al.
(arXiv:2604.05126); in-situ bicycle factory: Xu et al. (arXiv:2602.20546).

## License

MIT License. Copyright (c) 2026 Anthropic, PBC. Created by Matthew D. Schwartz;
code written by Claude (Anthropic) under his supervision. Documentation is
licensed CC BY 4.0. See LICENSE, LICENSE-CONTENT and NOTICE at the repository
root.
