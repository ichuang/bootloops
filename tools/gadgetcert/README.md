# gadgetcert

Fault properties of Clifford gadgets that act on unknown inputs (gates, error
correction with decoders), fault equivalence and refinement against a
specification, and composition of the resulting receipts. Built on `faultcert`
and `qeccert`.

    python3 -m gadgetcert analyze knill:steane,anc=spec --check ec:1 --close 1   # from tools/
    python3 tools/gadgetcert/selftest.py                                         # battery, from the repository root

See `GUIDE.md`.
