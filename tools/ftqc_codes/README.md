# ftqc_codes

Certified code facts (`n`, `k`, distances `d_x`/`d_z`/`d`, `t`, lightest
stabilizer weights, and structural flags) for CSS quantum error-correcting
codes, as receipts (`ftqc_codes-receipt-1`) with re-checkable evidence from up
to six independent legs: exact enumeration, SAT with an in-package checked LRAT
proof, a Lean 4 theorem built on Lean-QEC (arXiv:2605.16523), the BCH bound,
disjoint logical cuts, and checked code concatenation.

    python3 -m ftqc_codes certify golay -o golay.json     # from tools/, or with tools/ on PYTHONPATH
    python3 -m ftqc_codes verify golay.json --full
    python3 tools/ftqc_codes/selftest.py                  # battery, from the repository root

See [`GUIDE.md`](GUIDE.md) for scope, gates, hazards, and when to reach for
`ftqc_codes` vs sibling circuit-level tools (`ftqc_zx`, `ftqc_audit`,
`ftqc_gadgets`, `ftqc_extractors`). `examples/` holds two banked code specs and their
receipts; `tools/qeccert.py` preserves legacy `import qeccert` compatibility.
