# FTQC
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
    python3 -m ftqc check tools/ftqc/examples/column.json -o /tmp/column.json
    python3 -m ftqc verify tools/ftqc/examples/column.json /tmp/column.json
    python3 -m ftqc check tools/ftqc/examples/edge-extractor.json -o /tmp/edge.json
    python3 -m ftqc check tools/ftqc/examples/edge-hook.json -o /tmp/hook.json
    python3 -m ftqc verify tools/ftqc/examples/edge-hook.json /tmp/hook.json
    python3 run_selftests.py ftqc

Check returns 0 for a proved claim, 1 for a counterexample, 2 for invalid input, and 3 for an inconclusive search. Verify returns 0 for a valid positive or negative receipt. An exact search checks increasing weights. A failure with matching lower and upper endpoints establishes the exact distance. Increase the target by one to search for an upper witness.

Use --backend sat with a circuit to generate a proof or fault witness. Install a CaDiCaL with text LRAT support and set QECCERT_CADICAL to the executable. Verification uses only Python and the sibling qeccert package. The verifier reconstructs the CNF using backward propagation, then checks the LRAT proof, or replays the witness on the original schedule. The same site's Pauli alternatives are mutually exclusive. --seconds limits the external solver. --budget limits exhaustive combinations. Inconclusive receipts retain the reason. The default battery verifies a banked LRAT even when the solver is absent.

## Inputs and model

Each example is a complete JSON input. 'a' and 'b' are lists of binary strings, coordinate zero first. Both HGP blocks are retained, including the transpose-code sector.

* hgp takes a and b. Its d is the minimum classical distance over active seed and transpose sectors. For an empty logical quotient d is null.
* column adds column (zero-based), word (integer bit mask), vertices, edges, ports (one vertex per row coordinate), optional paths (one edge mask per seed check), and target. The active ports must be injective. The measured Z has support word in the indicated left-block column. Vertex, path and cycle checks are constructed explicitly. The subcode Z test applies up to the original HGP distance. For a larger target the tool searches the full merged CSS code. A failed reduction is lifted to a physical logical witness.
* canonical takes n, polynomial (cyclic check polynomial bit mask) and undirected edges on vertices 0 through n-1. Input labels define the identity single-column port map. It denotes the full extractor of Appendix F.2–F.3, including its canonical port priority (Eq. 64), path deformations and Appendix E.1 protocol. An optional assembly is the ordered edge list from full_graph and must match exactly. The result applies to that mathematical construction. Replay checks the reduced CSS matrices for every nonzero seed word, verifies that the examined cosets exhaust the surviving quotient, and independently computes the min-cut. The published theorem and HGP distance formula are mathematical dependencies. Replays do not formalize their proofs in Lean.
* circuit takes a, b, measured_z (physical bit mask), ops, readout, detectors and target. readout is a list of measurement-record indices whose XOR reports the requested logical Z. Each detector is a list of record indices with deterministic even parity. The ordered ops are a serial schedule. Qubits 0 through n-1 are arbitrary encoded input data. Larger labels are ancillas. RZ/RX prepare +Z/+X eigenstates, MZ/MX measure and retire an ancilla, CX takes control and target, H and I act on one live qubit. Z_IF takes a data/ancilla qubit and an earlier zero-based measurement index for classically controlled Z. Every ancilla lifetime needs a preparation and a terminating measurement. Include explicit I locations wherever idle faults should count.
* Every preparation can flip its eigenstate. Every measurement can flip its record. At each other location any nonidentity Pauli on its operands can follow the ideal operation, including on Z_IF's quantum operand. Witnesses are pairs [zero-based gate index, mode]. On a unitary, bits 0 and 1 select X and Z on its first operand, bits 2 and 3 on its second. A preparation or measurement uses mode 1. The input codespace and the terminal ideal data stabilizer syndrome define the boundaries. The checker also validates the ideal measured value and preservation of a spanning set of the surviving logical algebra. A damaging pattern evades all specified detectors and the terminal syndrome while changing the readout or protected data. This is a finite-block undetected-fault distance under those boundaries. Decoder thresholds, unlisted idle locations, and noisy terminal readout require their own models.

## Worked controls and limits

edge-extractor measures the two vertex checks Z(data)Z(edge) of the one-edge graph twice, retires the edge in X, and updates the Pauli frame. The product of the first two records measures logical Z on data 0 and 2. Its two parity detectors compare the vertex readings across rounds. edge-hook appends two identical CX gates on data 4 and 5. Their ideal action cancels. A Z on the target after the first gate spreads back through the second and evades the terminal syndrome. The smaller repeat-z example measures parity with a single meter; its banked LRAT exercises the checker on a cold clone.

The canonical benchmark is an implicit full architecture at phenomenological noise level. Gate-level receipts cover explicit schedules and their declared boundaries. Each receipt binds to a SHA-256 of the canonical input JSON. The exact checker redoes finite searches during replay. The SAT checker accepts the LRAT fragment supported by qeccert. Search time and proof size grow rapidly with threshold and schedule size. Resource limits yield open obligations.

## Addressable logical T

See [T_PROOFS.md](T_PROOFS.md) for the exact channel identity, fault contract, hardware assumptions and structural proof. Use the same check and verify commands with examples/t-gate.json and examples/t-hook.json. The exact and SAT backends both handle t_gate. A gate receipt has three cases: the T output and the two arms of a subsequent signed logical X measurement. A counterexample means at least one case failed the requested lower bound.

Inputs specify the HGP seeds, the anticommuting logical_x and logical_z masks on each block, rounds, target and resource. The contract for t_gate and t_bounds is ideal-encoded-T. verification_rounds defaults to one. representatives optionally supplies a pair of equivalent physical masks per measurement round. An experimental trailing_identity CNOT sequence must be checked as identity on every involved qubit. The t_bounds kind omits target and returns checked lower and upper bounds per case. The t_source kind takes resource=in-situ-all-plus and verifies the accepted preparation and resulting raw-Z counterexample. The 225-qubit benchmark uses the published classical seed with our selected single-site port. BootLoops works from the input matrices without requiring the source repository.

    export PYTHONPATH=tools
    python3 -m ftqc check tools/ftqc/examples/t-gate.json -o /tmp/t.json
    python3 -m ftqc verify tools/ftqc/examples/t-gate.json /tmp/t.json

The physical adapter is a verified-cat parity meter with direct access to the listed sites. Implementing the protocol on a particular fixed coupler network requires compiling and checking its routing and storage locations.

## Cross-checks with faultcert

Compare the two fault propagators on an extractor or fixed T arm, replay a logical witness, or export a verified cat for a native strict-FT analysis. The [comparison guide](FAULTCERT.md) gives the commands and the input, acceptance and output assumptions for each claim.

## Bivariate bicycle partial extractor

Use input kind `gross_z` for the published BB144 unprimed Z setting. The exact backend computes and certifies both CSS distances. The measurement command exports the merged matrices, ideal readout and split frame. See the [gross-code guide](GROSS.md) for the input, coordinate conventions, proof format and scope.

Input kind `gross_circuit` compiles the noisy BB144 operation. Use `schedule --probe X` or `schedule --probe Z` to export a native faultcert circuit, and `check` or `verify` for its distance. [Noise model, optional input fields and proof search](GROSS_CIRCUIT.md)
