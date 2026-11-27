#!/usr/bin/env python3
"""Optional standalone exporter from Python `qLDPC` into `ftqc_codes-import-1` JSON.

IMPORTANT: Run in an isolated environment where `qLDPC` (e.g. `qLDPC==0.2.9`,
`numpy<2.5`, `scipy`) is installed, WITHOUT `bootloops/tools` on `PYTHONPATH`
(putting `bootloops/tools` on `PYTHONPATH` shadows PyPI `galois` with
`tools/galois`; see `PLANS/qec-library-survey.md`).

Examples:
    /path/to/qec_producers/bin/python tools/ftqc_codes/producers/qldpc_export.py \
        --preset bb144 -o tools/ftqc_codes/examples/qldpc_bb144.import.json
    /path/to/qec_producers/bin/python tools/ftqc_codes/producers/qldpc_export.py \
        --preset lp6_2 --map-to-native -o tools/ftqc_codes/examples/qldpc_lp6_2.import.json
"""
import argparse
import hashlib
import json
import sys


def _check_environment():
    try:
        import galois
    except ImportError as exc:
        raise SystemExit(
            f"PyPI 'galois' not installed in this interpreter ({exc}). "
            "Run qldpc_export.py in an isolated qLDPC venv."
        ) from None
    if not hasattr(galois, "GF"):
        raise SystemExit(
            "Import name 'galois' resolved to a package without galois.GF "
            "(likely BootLoops tools/galois on PYTHONPATH). "
            "Unset PYTHONPATH before running qldpc_export.py."
        )
    try:
        import numpy as np
        import qldpc
        import scipy
    except Exception as exc:
        raise SystemExit(f"Failed to import qLDPC producer stack: {exc}") from None
    return qldpc, np, scipy, galois


def _css_sha256_from_supports(n, hx_supp, hz_supp):
    def row_str(supp):
        s = set(supp)
        return "".join("1" if j in s else "0" for j in range(n))

    txt = (
        f"n={n}\nHX\n"
        + "\n".join(row_str(r) for r in hx_supp)
        + "\nHZ\n"
        + "\n".join(row_str(r) for r in hz_supp)
        + "\n"
    )
    return hashlib.sha256(txt.encode("utf-8")).hexdigest()


def _supports_sha256(supports):
    blob = json.dumps([sorted(int(x) for x in r) for r in supports], separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def _matrix_to_supports(mat):
    out = []
    for r in range(mat.shape[0]):
        out.append([int(c) for c in range(mat.shape[1]) if int(mat[r, c]) & 1])
    return out


def build_preset(preset, *, map_to_native=False):
    qldpc, np, scipy, galois = _check_environment()
    from qldpc.abstract import CyclicGroup, RingMember
    from qldpc.objects import Pauli
    import sympy

    cmap = None
    if preset == "bb144":
        x, y = sympy.symbols("x y")
        code = qldpc.codes.BBCode({x: 12, y: 6}, x**3 + y + y**2, y**3 + x + x**2)
        name = "BB144"
        prov = {
            "paper": "Bravyi et al., High-threshold and low-overhead fault-tolerant quantum memory",
            "version": "arXiv:2308.07915 / Nature 627, 778 (2024)",
            "locator": "[[144, 12, 12]] Gross code: l=12, m=6, A=x^3+y+y^2, B=y^3+x+x^2",
        }
        p_claims = {"d": 12, "d_x": 12, "d_z": 12}
    elif preset == "lp6_2":
        G = CyclicGroup(3)
        elems = list(G.generate())
        el = RingMember(G, elems[0], elems[1])
        code = qldpc.codes.LPCode([[el]], [[el]])
        name = "LP(6,2)"
        prov = {
            "paper": "Panteleev-Kalachev lifted product control",
            "version": "PLANS/high-rate-qec.md & PLANS/qec-library-survey.md",
            "locator": "LPCode([[1+x]], [[1+x]]) over GF(2)[C_3]",
        }
        p_claims = {"d": 2, "d_x": 2, "d_z": 2}
        if map_to_native:
            cmap = {
                "convention": "push_external_to_native",
                "perm_q": [0, 1, 2, 5, 3, 4],
                "perm_hx": [0, 1, 2],
                "perm_hz": [2, 0, 1],
            }
    else:
        raise SystemExit(f"unknown preset {preset!r}; choose bb144 or lp6_2")

    hx_mat = np.asarray(code.matrix_x, dtype=int)
    hz_mat = np.asarray(code.matrix_z, dtype=int)
    lx_mat = np.asarray(code.get_logical_ops(Pauli.X), dtype=int)
    lz_mat = np.asarray(code.get_logical_ops(Pauli.Z), dtype=int)

    n = int(hx_mat.shape[1])
    k = int(lx_mat.shape[0])
    hx_supp = _matrix_to_supports(hx_mat)
    hz_supp = _matrix_to_supports(hz_mat)
    lx_supp = _matrix_to_supports(lx_mat)
    lz_supp = _matrix_to_supports(lz_mat)
    raw_sha = _css_sha256_from_supports(n, hx_supp, hz_supp)

    if cmap is not None:
        pq, px, pz = cmap["perm_q"], cmap["perm_hx"], cmap["perm_hz"]
        hx_t = [None] * len(hx_supp)
        for i, r in enumerate(hx_supp):
            hx_t[px[i]] = sorted(pq[c] for c in r)
        hz_t = [None] * len(hz_supp)
        for i, r in enumerate(hz_supp):
            hz_t[pz[i]] = sorted(pq[c] for c in r)
        final_sha = _css_sha256_from_supports(n, hx_t, hz_t)
    else:
        final_sha = raw_sha

    return {
        "schema": "ftqc_codes-import-1",
        "name": name,
        "code_type": "css",
        "field": "GF(2)",
        "index_base": 0,
        "producer": {
            "name": "qLDPC",
            "version": str(qldpc.__version__),
            "python": sys.version.split()[0],
            "numpy": str(np.__version__),
            "scipy": str(scipy.__version__),
            "galois": str(galois.__version__),
        },
        "parameters": {
            "n": n,
            "k": k,
            "rows_hx": len(hx_supp),
            "rows_hz": len(hz_supp),
        },
        "raw_sha256": raw_sha,
        "sha256": final_sha,
        "coordinate_map": cmap,
        "sources": {
            "hx": {
                "format": "sparse_rows",
                "rows": len(hx_supp),
                "n": n,
                "index_base": 0,
                "sha256": _supports_sha256(hx_supp),
                "supports": hx_supp,
            },
            "hz": {
                "format": "sparse_rows",
                "rows": len(hz_supp),
                "n": n,
                "index_base": 0,
                "sha256": _supports_sha256(hz_supp),
                "supports": hz_supp,
            },
        },
        "candidate_logicals": {
            "convention": "positive_css",
            "coordinates": "external",
            "slot_names": [f"q{i}" for i in range(k)],
            "signs_x": [0] * k,
            "signs_z": [0] * k,
            "lx": {
                "format": "sparse_rows",
                "rows": k,
                "n": n,
                "index_base": 0,
                "sha256": _supports_sha256(lx_supp),
                "supports": lx_supp,
            },
            "lz": {
                "format": "sparse_rows",
                "rows": k,
                "n": n,
                "index_base": 0,
                "sha256": _supports_sha256(lz_supp),
                "supports": lz_supp,
            },
        },
        "provenance": prov,
        "paper_claims": p_claims,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", choices=("bb144", "lp6_2"), required=True)
    ap.add_argument("--map-to-native", action="store_true")
    ap.add_argument("-o", "--out")
    args = ap.parse_args(argv)
    doc = build_preset(args.preset, map_to_native=args.map_to_native)
    txt = json.dumps(doc, indent=2) + "\n"
    if args.out:
        with open(args.out, "w") as f:
            f.write(txt)
    else:
        sys.stdout.write(txt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
