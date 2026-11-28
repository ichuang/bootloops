# Gross-code unprimed Z measurement

`gross_z` implements the unprimed logical Z setting of the 103-ancilla partial extractor of Cross, He, Rall and Yoder. The full apparatus can select eight measurements. This setting uses 18 additional code coordinates. Both the original and the merged CSS codes have exact static distance 12. The measurement export gives ideal commuting Pauli operations and a Pauli frame for the split. The compiled noisy measurement has exact circuit distances seven and ten at seven and twelve rounds. [Noisy schedule and receipts](GROSS_CIRCUIT.md) [Construction, §4.1 and Fig. 9](https://arxiv.org/abs/2407.18393v4) · [Partial extractor interpretation, §4.4](https://arxiv.org/abs/2503.10390v2)

```sh
# From the BootLoops root
export PYTHONPATH=tools
work=$(mktemp -d)
spec=tools/ftqc_extractors/examples/gross-z.json
python3 -m ftqc_extractors measurement "$spec" -o "$work/measurement.json"
python3 -m ftqc_extractors check "$spec" -o "$work/receipt.json"
python3 -m ftqc_extractors verify "$spec" "$work/receipt.json"
```

The input has `schema: 1`, `kind: "gross_z"` and an integer `target` for the merged distance. Generation uses the exact backend and one shared search-node budget across both codes and axes. The example asks for 12. A target of 13 produces a verified counterexample at weight 12. A small budget can yield an inconclusive result. Generation and replay require only Python and the sibling packages in this checkout.

## Reconstruct the measured code

Use the `ftqc_codes` BB144 ordering. Cell (i,j) is coordinate 6i+j on the left and 72+6i+j on the right, with periods 12 and 6. The X checks have polynomial blocks A = x³+y+y² and B = y³+x+x². The support S from Eqs. (63)–(64) has six sites per half. Exactly 18 X checks meet S, each at two sites. Label each such check by its cell and assign it one ancilla data coordinate in increasing label order. They form edges of a connected graph whose 12 vertices are the port sites.

* Extend each incident memory X check onto its own ancilla. Keep the other memory checks. Add the X gauge check numbered 6 in Fig. 9 on the ancillas labelled y⁴, xy⁴, x³y² and x³y. Its physical support is [146,150,155,156].
* For each port v, add a Z check on v and the ancillas for edges incident to v. Multiplying all 12 port checks cancels every edge twice and leaves the published logical Z on the data.
* The merged code has 162 coordinates and X and Z ranks 73 and 78. Its stabilizers acting only as Z on memory data span precisely the old Z stabilizers and the measured Z. Projecting the merged X centralizer back to memory gives all old X centralizers commuting with the measured Z. The merged X stabilizers project to the original X stabilizers. Old Z logicals lift on data alone, so their commutation with these lifted X representatives is preserved. This identifies the surviving eleven logical pairs.

As a negative control, deleting gauge check 6 leaves its four-site support as a nontrivial X logical. It commutes with all Z checks and is independent of the remaining X rows. The package tests verify both facts.

## Interpret the ideal records

Prepare the 18 edge ancillas in |+⟩, measure the 12 Z ports, then measure every ancilla in X. In each measurement list record 0 denotes eigenvalue +1 and record 1 denotes −1. The XOR of the port bits gives the logical Z outcome. Each exported `data_z_frame` entry names the ancilla X bits whose parity calls for a Z correction on one data site. The `x_record_constraints` must have even parity on ideal records.

To derive this rule, write D for the vertex-by-edge binary incidence matrix. On a computational basis state with port bits u, the first measurements with record b select the edge strings a satisfying D a = b + u. Connectedness makes this solvable precisely when the parity of u matches the parity of b. An X measurement record m contributes a Fourier sign (−1)^(m·a). Summing over the edge strings cancels unless m = Dᵀc for a vertex potential c. The remaining sign is (−1)^(c·b+c·u). Its first term is a branch phase, and its second is the data Z frame. With c fixed to zero at the root, its value at each port is the parity along the exported tree path. The code checks every path boundary and the rank of the cycle constraints. The corrected branch therefore projects onto the reported Z eigenspace for arbitrary encoded states, including those entangled with spectators.

## Check the distance proof

For each code and each Pauli axis the receipt supplies a weight-12 support, the check-matrix digest and a compressed tree proving that all smaller logical supports are excluded. The verifier reconstructs the matrices and checks each witness for zero syndrome and nonmembership in the same-axis stabilizer span. For a CSS code a nontrivial mixed Pauli has a nontrivial X or Z component, so the two sector bounds suffice.

Suppose a minimum nontrivial support of weight at most 11 exists. Start at its smallest site. A node records the sites chosen so far, the available sites, the syndrome, and the number of remaining picks. The proof uses these rules.

1. An odd check needs an odd number of additional sites on that check. Branch over the first available site on the check belonging to the putative support, excluding earlier alternatives from later branches. This covers every completion. An odd check with no available site closes the branch.
2. If r picks remain, r+1 odd checks with pairwise disjoint available supports need too many sites. The proof lists the check indices, and the verifier checks their parities, distinctness and disjointness.
3. A nonzero syndrome with no picks left closes the branch. At zero syndrome the chosen support must be a stabilizer. A minimum nontrivial logical cannot contain such a nonempty subset because removing it gives a smaller representative of the same logical. A zero-syndrome nonstabilizer refutes the lower bound.

The checker traverses every root and every required child using the split or packing supplied at each internal node. It checks stabilizer leaves with a separate row-space implementation and requires exact consumption of the proof bytes and node count. Bounds and witnesses are checked for both the memory and the surgery. [Generator and proof checker](css_cover.py) · [Certificate assembly](gross_distance.py)
