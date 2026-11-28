# ftqc_decoders

Code and circuit distance describe an ideal minimum-weight decoder.
`ftqc_decoders` checks concrete polynomial-time decoders instead. Given a CSS
code from `ftqc_codes`, a circuit record model from `ftqc_faults` or
`ftqc_extractors`, a surgery cone with metachecks, or a Stim `.dem` file, it
runs decoders in exact arithmetic, certifies that every pattern up to a small
radius is corrected or returns the pattern that is not, finds trapping-set
witnesses, computes logical-error bounds under printed assumptions and a
two-sided linked-cluster enclosure for localized decoders, and audits
the output of external decoders such as Tesseract. Results are
`ftqc_decoders-receipt-1` receipts.

## Quickstart

From the BootLoops repository root:

```sh
export PYTHONPATH=tools
work=$(mktemp -d)

# 1. Run the battery (about 2 s)
python3 tools/ftqc_decoders/selftest.py

# 2. Expose a weight-3 symmetric split-stabilizer trapping set on BB144 [[144,12,12]]
#    under plain BP, and check that Relay-BP, BP+LSD and BP+OSD each correct it.
#    The fixture retains its original filename. Its problem.kind field selects bb144.
python3 -m ftqc_decoders check tools/ftqc_decoders/examples/bb72_trapping_set.json -o "$work/bb144-trap.json"
python3 -m ftqc_decoders verify "$work/bb144-trap.json"

# 3. Export BB144 [[144,12,12]] to Stim .dem and decode with audited Tesseract
#    Run this optional example with stim and tesseract-decoder installed.
python3 -m ftqc_decoders convert-dem --code bb144 -o "$work/bb144.dem"
python3 -m ftqc_decoders decode \
  --code "$work/bb144.dem" --faults 0,18,41 --decoder external:tesseract
```

## Included Decoders & Verifiers

| Algorithm / Study | Module | Primary References |
|---|---|---|
| `bp_osd` (`NMS` & exact-$\mathbb{Q}$ Sum-Product + `OSD-0`/`OSD-CS`) | `bp_osd.py` | [Panteleev & Kalachev 2021](https://arxiv.org/abs/1904.02703) · [Roffe et al. 2020](https://arxiv.org/abs/2005.07016) |
| `bp_lsd` (Localized Statistics Decoding, on-the-fly cluster `GF(2)`) | `bp_osd.py` | [Hillmann et al. 2025](https://arxiv.org/abs/2406.18655) |
| `relay_bp` (Disordered-memory message passing without `OSD`) | `relay_bp.py` | [Müller et al. 2025](https://arxiv.org/abs/2506.01779) · [Bravyi et al. 2024](https://arxiv.org/abs/2308.07915) |
| `union_find` & `belief_find` (Half-edge peeling + hypergraph `GF(2)`) | `union_find.py` | [Delfosse & Nickerson 2021](https://arxiv.org/abs/1709.06218) · [Delfosse, Londe & Beverland 2022](https://arxiv.org/abs/2103.08049) |
| `mwpm` (Disjoint-lineage Blossom with LP strong-duality check) | `mwpm.py` | [Dennis et al. 2002](https://arxiv.org/abs/quant-ph/0110143) · Edmonds (1965) |
| `ssf` & `metacheck_decode` (Single-shot HGP + 2-stage surgery repair) | `ssf_metacheck.py` | [Leverrier, Tillich & Zémor 2015](https://arxiv.org/abs/1504.00822) · [Cowtan et al. 2025](https://arxiv.org/abs/2510.14895) |
| `sat_mwd` (Inhomogeneous Sinz-CNF + `ftqc_codes.lrat` optimality proof) | `sat_mwd.py` | [Berent et al. 2024](https://arxiv.org/abs/2303.14237) · Sinz (2005) |
| `external:tesseract`, `pymatching`, `ldpc_*`, `chromobius` + `Sinter` | `adapters/` | [Beni et al. 2025](https://arxiv.org/abs/2503.10988) · [Higgott & Gidney 2025](https://arxiv.org/abs/2303.15933) |
| `certify_radius`, `trapping_set`, `enclose_*`, & sibling bridges | `effective_distance.py`, `enclosures.py`, `bridges.py` | [Poulin & Chung 2008](https://arxiv.org/abs/0801.1241) · [Raveendran & Vasić 2021](https://arxiv.org/abs/2012.15297) · [Gottesman 2014](https://arxiv.org/abs/1310.2984) |

See [GUIDE.md](GUIDE.md) for full CLI options, receipt fields, operational-vs-postselected semantics, and the complete bibliography, and [../FTQC.md](../FTQC.md) for cross-tool workflows with `ftqc_codes`, `ftqc_faults.audit`, `ftqc_extractors`, `ftqc_gadgets`, `ftqc_compose`, and `ftqc_uncompute`.
