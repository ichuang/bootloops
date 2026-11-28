# ftqc_gadgets

What a gadget with its decoder does to an arbitrary encoded input, under
incoming errors and internal faults: the gate and error-correction properties
of Aliferis, Gottesman and Preskill at low order, refinement against a
specification, and composition of the resulting receipts. Each input enters as
half of an encoded Bell pair, so "for every input" is a property of one
stabilizer state. Built on the audit engine of `ftqc_faults` and on
`ftqc_codes`.

    export PYTHONPATH="$PWD/tools"                       # from the BootLoops root
    python3 -m ftqc_gadgets analyze knill:steane,anc=spec --check ec:1 --close 1
    python3 tools/ftqc_gadgets/selftest.py               # battery

For codes too large to enumerate, the coordinate-locality certificate
(`--local`) and `levels.py` carry the same properties to the outer level of a
concatenated scheme, and the factory certificate (`--tree`) does so for the
purification trees that make ancillas. `cones` gives the fault cone of every
decoding event of a whole correction schedule. Signed surgery instruments
([SURGERY.md](SURGERY.md)) and staged switched-T and pieceable-CCZ gates
([STAGED.md](STAGED.md)) have their own guides.

The statements hold for fault sets up to a small order; the lemmas that extend
them are printed in the receipts and tested, not machine-checked. See the
[guide](GUIDE.md) and [`../FTQC.md`](../FTQC.md).
