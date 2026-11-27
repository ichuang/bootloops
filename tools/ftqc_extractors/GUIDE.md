# ftqc_extractors
Version 1 input is a JSON object with schema integer 1, kind, and model-specific fields. Bit zero is coordinate zero. HGP qubits use row-major ordering in the two product blocks. Checks and supports are binary.

## Claims
* Base distance counts nontrivial Pauli normalizer elements modulo stabilizers.
* Static surgery compares the deformed CSS code with the input code (Blue et al., arXiv:2606.03507v2, Definition 11). Definition 16 and Theorem 17 concern single-column Z surgery. Cap comparisons by the base distance when applying the reduction.
* The canonical full cyclic extractor in Appendix F has phenomenological lower bound min(d,k²,k λ). The inputs must establish every nonzero single-column case and use the canonical assembly and d-round protocol of Appendix E.1. k is the seed dimension.
* A circuit claim counts faulty scheduled locations with ideal input code and final syndrome boundary. It requires a checked logical measurement and preservation of the remaining logical algebra. The weight of a correlated two-qubit Pauli at one gate is one.
* A lower bound is obtained by ruling out every damaging undetected fault below the threshold. Exact distance also requires a witness at the threshold. A resource limit leaves an open obligation.

Credit: He, Cowtan, Williamson and Yoder, arXiv:2503.10390; Blue, He, Zhou and Chuang, arXiv:2606.03507v2.

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

## Cross-checks with ftqc_audit

Compare the two fault propagators on an extractor or fixed T arm, replay a logical witness, or export a verified cat for a native strict-FT analysis. The [comparison guide](AUDIT.md) gives the commands and the input, acceptance and output assumptions for each claim.

## Bivariate bicycle partial extractor

Use input kind `gross_z` for the published BB144 unprimed Z setting. The exact backend computes and certifies both CSS distances. The measurement command exports the merged matrices, ideal readout and split frame. See the [gross-code guide](GROSS.md) for the input, coordinate conventions, proof format and scope.

Input kind `gross_circuit` compiles the noisy BB144 operation. Use `schedule --probe X` or `schedule --probe Z` to export a native ftqc_audit circuit, and `check` or `verify` for its distance. [Noise model, optional input fields and proof search](GROSS_CIRCUIT.md)

## Full block reading

Run `check` and `verify` on [the two-block fixture](examples/block-reading.json). Matrices are arrays of nonnegative integer row masks, least significant bit first. Each independent pattern row selects blocks. Its port parities measure one logical per member of the memory logical basis. The receipt gives the data-only projections, surviving dimension, deformed matrices and metachecks. In each pattern row, decompose the measured auxiliary split syndrome in the recorded syndrome basis and apply the corresponding correction on every selected block (Z for Z reading, X for its dual). The static bound applies Lemma 3.8 after checking the input distance.

## Fixed-coupler measurement graphs

[Four-port example](examples/fixed-graph.json). Signed strings use wire zero first. A port of -1 excludes that physical site. The theorem covers all logical representatives supported on the included sites. The verifier checks relative expansion on all cuts, every even subset via generating pairs within each check’s allowed routing pool, a complete cycle basis, the input distance and the maximum degree on both sides of the physical check-data wiring. These are finite-instance degree counts. The selected logical produces signed deformed checks, the parity of all vertex outcomes with its sign, and split paths in support order. Edge qubits use the vertex-Z convention: prepare plus, measure X on split, and apply the local factor of the measured Pauli for each odd path. The phenomenological protocol uses target rounds before, during and after the merge.

## Padded hypergraph surgery

[Three-block sequence](examples/three-block-sequence.json) measures two independent patterns in order. Each auxiliary CSS code has a fourth term supplied by its complete sparse cycle basis. The three chain maps are **column** masks into the original memory. Metacheck masks refer to original Z checks followed by the operation’s vertex checks. The verifier builds the simultaneous mapping cone for the distance proof; the physical protocol performs the listed operations sequentially with one merge round apiece and the stated boundary rounds. It checks both chain squares, vertex non-overlap, sparsity, injectivity and independence of measured logical classes. This implementation uses the stronger sufficient condition of full auxiliary CSS distance. A failed sufficient condition leaves the theorem unavailable; it is not a physical fault witness. The distance model is phenomenological.

## Two-block BB144 benchmark

Run check and verify on [Z reading](examples/gross-block-Z.json) or [X reading](examples/gross-block-X.json). Each measures twelve independent corresponding parities and leaves twelve logical qubits. Each adds 72 data qubits and 144 check meters to the two memories. Proposition 3.9 gives phenomenological distance at least 12 with twelve syndrome rounds at each boundary. The underlying BB144 distance exclusion runs independently in both sectors.

The circuit benchmark uses ideal encoded inputs, one noisy merge, noisy split and ideal terminal readout. Greedy disjoint CNOT layers have explicit idles on live wires. Separate X and Z probes protect the readout and surviving algebra. Use `schedule SPEC --probe X` or `--probe Z` to export them for ftqc_audit. For either reading the two probes give circuit distance in [3,12]: all faults on one or two distinct locations are excluded, and a weight-12 fault on data is replayed. The basis label on the dual experiment denotes the probe before global Hadamard conjugation.

## Expansion and intermediate time

[Modular expansion](examples/modular-expansion.json) checks the exact minimum in Definition 7.2 before and after Cartesian path thickening, with ports on the first layer. It returns both minimizing cuts. [Intermediate reading](examples/intermediate-reading.json) applies Theorem 6.5. Breadth-first search on attainable Pauli syndromes computes quantum soundness from minimum physical Pauli weights, including Y. Each individual cone must preserve the memory target and its auxiliary distance times its declared merge rounds must reach that target. Both temporal boundaries are checked. These enumerations are intended for finite studies and return an open obligation at their budget.
