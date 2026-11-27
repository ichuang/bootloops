# ftqc_gadgets

Fault properties of Clifford gadgets that act on unknown inputs (gates, error
correction with decoders), fault equivalence and refinement against a
specification, and composition of the resulting receipts. Built on `ftqc_audit`
and `ftqc_codes`.

    python3 -m ftqc_gadgets analyze knill:steane,anc=spec --check ec:1 --close 1   # from tools/
    python3 tools/ftqc_gadgets/selftest.py                                         # battery, from the repository root

The coordinate-locality certificate (`--local`) and `levels.py` extend the same
questions to large codes and to the outer level of a concatenated scheme; the
factory certificate (`--tree`, `factories.py`) does the same for the purification
trees that make the ancillas. `cones.py` (`cones` command) gives the fault cone of
every decoding event of a whole program with corrections placed sparsely.

[Run and interpret the cone example](GUIDE.md#fault-cone-receipts) · [Full gadget guide](GUIDE.md)
