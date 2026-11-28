# ftqc_codes

Certified facts about CSS quantum error-correcting codes: n, k, d_X, d_Z, t,
the lightest stabilizer weights and structural flags, each bound carried as an
interval with a replayable witness and the proof leg behind it (exhaustive
enumeration, SAT with a re-checked LRAT proof, a Lean 4 theorem on Lean-QEC,
the BCH bound, logical cuts, concatenation, symmetry orbits). Also
construction receipts for high-rate families, named logical bases, and signed
receipts for code automorphisms and transversal Cliffords.

    export PYTHONPATH="$PWD/tools"                       # from the BootLoops root
    python3 -m ftqc_codes certify golay -o golay.json
    python3 -m ftqc_codes verify golay.json --full
    python3 -m ftqc_codes construct pp256 -o pp256.construct.json
    python3 -m ftqc_codes basis lp6_2 -o lp6_2.basis.json
    python3 tools/ftqc_codes/selftest.py                 # battery

Static distance says nothing about circuits or decoders; see the
[guide](GUIDE.md) for scope, gates and hazards and [`../FTQC.md`](../FTQC.md)
for how the FTQC tools fit together. `tools/qeccert.py` keeps the old
`import qeccert` working.
