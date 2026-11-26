# FTQC
Version 1 input is a JSON object with schema integer 1, kind, and model-specific fields. Bit zero is coordinate zero. HGP qubits use row-major ordering in the two product blocks. Checks and supports are binary.

## Claims
* Base distance counts nontrivial Pauli normalizer elements modulo stabilizers.
* Static surgery compares the deformed CSS code with the input code (Blue et al., arXiv:2606.03507v2, Definition 11). Definition 16 and Theorem 17 concern single-column Z surgery. Cap comparisons by the base distance when applying the reduction.
* The canonical full cyclic extractor in Appendix F has phenomenological lower bound min(d,k²,k λ). The inputs must establish every nonzero single-column case and use the canonical assembly and d-round protocol of Appendix E.1. k is the seed dimension.
* A circuit claim counts faulty scheduled locations with ideal input code and final syndrome boundary. It requires a checked logical measurement and preservation of the remaining logical algebra. The weight of a correlated two-qubit Pauli at one gate is one.
* A lower bound is obtained by ruling out every damaging undetected fault below the threshold. Exact distance also requires a witness at the threshold. A resource limit leaves an open obligation.

Credit: He, Cowtan, Williamson and Yoder, arXiv:2503.10390; Blue, He, Zhou and Chuang, arXiv:2606.03507v2.
