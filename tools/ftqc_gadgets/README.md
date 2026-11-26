# ftqc_gadgets

Fault properties of Clifford gadgets that act on unknown inputs (gates, error
correction with decoders), fault equivalence and refinement against a
specification, and composition of the resulting receipts. Built on `ftqc_audit`
and `ftqc_codes`.

    python3 -m ftqc_gadgets analyze knill:steane,anc=spec --check ec:1 --close 1   # from tools/
    python3 tools/ftqc_gadgets/selftest.py                                         # battery, from the repository root

The coordinate-locality certificate (`--local`) and `levels.py` extend the same
questions to large codes and to the outer level of a concatenated scheme.

See `GUIDE.md`.
