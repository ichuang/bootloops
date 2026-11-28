# ftqc_codes

Certified code facts (`n`, `k`, distances `d_x`/`d_z`/`d`, `t`, lightest
stabilizer weights, and structural flags), construction receipts
(`ftqc_codes-construction-1`), sparse/MatrixMarket interchange
(`ftqc_codes-import-1`), and addressable paired logical basis receipts
(`ftqc_codes-basis-1`) for CSS quantum error-correcting codes, with
re-checkable evidence from up to seven independent legs: exact enumeration, SAT
with an in-package checked LRAT proof, a Lean 4 theorem built on Lean-QEC
(arXiv:2605.16523), the BCH bound, disjoint logical cuts, checked code
concatenation, and symmetry-reduced automorphism orbit exclusion.

    python3 -m ftqc_codes certify golay -o golay.json      # from tools/, or with tools/ on PYTHONPATH
    python3 -m ftqc_codes construct pp256 -o pp256.construct.json
    python3 -m ftqc_codes basis lp6_2 -o lp6_2.basis.json
    python3 -m ftqc_codes verify golay.json --full
    python3 tools/ftqc_codes/selftest.py                   # battery, from the repository root

See [`GUIDE.md`](GUIDE.md) for scope, gates, hazards, and when to reach for
`ftqc_codes` vs sibling circuit-level tools (`ftqc_faults`, `ftqc_decoders`,
`ftqc_gadgets`, `ftqc_extractors`, `ftqc_compose`). `examples/` holds banked
code specs, external import manifests, and receipts; `tools/qeccert.py`
preserves legacy `import qeccert` compatibility.
