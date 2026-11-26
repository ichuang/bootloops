# qeccert

Certified code facts (n, k, distances, lightest stabilizer, structural flags)
for CSS quantum error-correcting codes, as receipts with re-checkable evidence
from up to three independent legs: exact enumeration, SAT with a checked LRAT
proof, and a Lean 4 theorem built on Lean-QEC (arXiv:2605.16523).

    python3 -m qeccert certify golay -o golay.json     # from tools/, or with tools/ on PYTHONPATH
    python3 -m qeccert verify golay.json --full
    python3 tools/qeccert/selftest.py                  # battery, from the repository root

See `GUIDE.md` for scope, gates and hazards. `examples/` holds two code specs
and their receipts.
