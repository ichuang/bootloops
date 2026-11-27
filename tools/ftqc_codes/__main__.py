"""Command line:  python3 -m ftqc_codes <command> ...   (with tools/ on PYTHONPATH)

  certify   SPEC [-o receipt.json] [--legs A,B,C,D,E,G|none] [--target q=v ...] [--max-combos N]
  construct SPEC [-o construction.receipt.json] [--static-receipt receipt.json]
  basis     SPEC [-o basis.receipt.json] [--method binary|candidate|spectral]
  verify    RECEIPT.json [--full] [--static-receipt receipt.json] [--code SPEC]
  spec      SPEC [-o code.json]          write the code's matrices as JSON
  families                               list family constructors
SPEC is a code/import/construction JSON file or family[:args], e.g. golay, qr:47, toric:5, bb72, radial90, pp256.

Exit codes: 0 ok, 1 refusal (not a CSS code, receipt fails, legs disagree),
2 usage or malformed spec.
"""
import argparse
import json
import os
import sys

from . import basis, certify, codes, constructions, symmetry


def _targets(items):
    out = {}
    for it in items or []:
        q, _, v = it.partition("=")
        if q not in ("d_x", "d_z", "wmin_x", "wmin_z", "d") or not v.isdigit():
            raise codes.SpecError(f"bad --target {it!r}; use d=7, d_x=7, wmin_x=8, ...")
        if q == "d":
            out["d_x"] = out["d_z"] = int(v)
        else:
            out[q] = int(v)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ftqc_codes", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("certify")
    c.add_argument("spec")
    c.add_argument("-o", "--out")
    c.add_argument("--legs", default="A")
    c.add_argument("--target", action="append")
    c.add_argument("--max-combos", type=int, default=5_000_000)
    c.add_argument("--workdir", help="directory for CNF, LRAT and Lean files (legs B, C)")
    c.add_argument("--sat-timeout", type=int, default=600)
    c.add_argument("--lean-timeout", type=int, default=3600)
    cs = sub.add_parser("construct")
    cs.add_argument("spec")
    cs.add_argument("-o", "--out")
    cs.add_argument("--static-receipt")
    bs = sub.add_parser("basis")
    bs.add_argument("spec")
    bs.add_argument("-o", "--out")
    bs.add_argument("--method", choices=("binary", "candidate", "spectral"), default="binary")
    v = sub.add_parser("verify")
    v.add_argument("receipt")
    v.add_argument("--full", action="store_true")
    v.add_argument("--static-receipt")
    v.add_argument("--code")
    s = sub.add_parser("spec")
    s.add_argument("spec")
    s.add_argument("-o", "--out")
    sub.add_parser("families")
    a = ap.parse_args(argv)
    try:
        if a.cmd == "families":
            for k, (_, doc) in sorted(codes.FAMILIES.items()):
                print(f"{k:16s} {doc}")
            print("named:", " ".join(sorted(codes.NAMED)))
            return 0
        if a.cmd == "spec":
            txt = json.dumps(codes.load(a.spec).to_json(), indent=1)
            if a.out:
                with open(a.out, "w") as f:
                    f.write(txt + "\n")
            else:
                print(txt)
            return 0
        if a.cmd == "construct":
            code = codes.load(a.spec)
            srec = None
            if a.static_receipt:
                with open(a.static_receipt) as f:
                    srec = json.load(f)
            base_dir = os.path.dirname(os.path.abspath(a.spec)) if a.spec.endswith(".json") else None
            crec = constructions.build_construction_receipt(code, static_receipt=srec, base_dir=base_dir)
            constructions.verify_construction_receipt(crec, static_receipt=srec, base_dir=base_dir)
            m = crec["metrics"]
            print(f"{crec['payload']['name']}: n={m['n']} k={m['k']} rate={m['rate'][0]}/{m['rate'][1]} "
                  f"ancestry={crec['payload']['ancestry']} sha256={crec['code_sha256'][:12]}")
            if a.out:
                with open(a.out, "w") as f:
                    json.dump(crec, f, indent=1)
                    f.write("\n")
            return 0
        if a.cmd == "basis":
            code = codes.load(a.spec)
            if a.method == "binary":
                brec = basis.derive_paired_basis(code)
            elif a.method == "candidate":
                brec = basis.build_basis_receipt(code, source="imported_candidate")
            else:
                cmeta = code.meta.get("construction") or {}
                if cmeta.get("kind") != "lifted_product":
                    raise certify.Refusal("spectral basis CLI requires a cyclic lifted_product spec")
                si = cmeta["seed_inputs"]
                gc = cmeta["generator_convention"]
                _, brec, _ = symmetry.spectral_abelian_basis(
                    si["lift"],
                    si["a_reduced_exponents"],
                    si["b_reduced_exponents"],
                    coeff_encoding="exponents",
                    coordinate_order=gc["coordinate_order"],
                    name=code.name,
                )
            basis.verify_basis_receipt(code, brec)
            print(f"{code.name}: k={brec['k']} source={brec['source']} "
                  f"code_sha256={brec['code_sha256'][:12]} basis_sha256={brec['basis_sha256'][:12]}")
            if a.out:
                with open(a.out, "w") as f:
                    json.dump(brec, f, indent=1)
                    f.write("\n")
            return 0
        if a.cmd == "verify":
            with open(a.receipt) as f:
                rec = json.load(f)
            if rec.get("schema") == constructions.RECEIPT_SCHEMA:
                srec = None
                if a.static_receipt:
                    with open(a.static_receipt) as sf:
                        srec = json.load(sf)
                base_dir = os.path.dirname(os.path.abspath(a.receipt))
                for line in constructions.verify_construction_receipt(rec, static_receipt=srec, base_dir=base_dir):
                    print("ok  ", line)
                m = rec["metrics"]
                print(f"{rec['payload']['name']}: n={m['n']} k={m['k']} rate={m['rate'][0]}/{m['rate'][1]} "
                      f"ancestry={rec['payload']['ancestry']}")
                return 0
            if rec.get("schema") == basis.SCHEMA:
                if not a.code:
                    raise codes.SpecError("verifying a basis receipt requires --code SPEC")
                code = codes.load(a.code)
                for line in basis.verify_basis_receipt(code, rec):
                    print("ok  ", line)
                print(f"{rec['code_name']}: k={rec['k']} basis_sha256={rec['basis_sha256'][:12]}")
                return 0
            for line in certify.verify(rec, full=a.full):
                print("ok  ", line)
            print(certify.summary(rec))
            return 0
        raw_legs = [x.strip().upper() for x in a.legs.split(",") if x.strip()]
        legs = [] if raw_legs == ["NONE"] else raw_legs
        if set(legs) - {"A", "B", "C", "D", "E", "G"}:
            raise codes.SpecError("--legs takes a subset of A,B,C,D,E,G or none")
        code = codes.load(a.spec)
        targets = _targets(a.target)
        rec = certify.build(code, targets, a.max_combos, legs)
        if "D" in legs:
            from . import bch
            bch.run_leg(rec, code)
        if "E" in legs:
            from . import cuts
            cuts.run_leg(rec, code)
        if "G" in legs:
            from . import orbit_minweight
            orbit_minweight.run_leg(rec, code, targets=targets, max_nodes=a.max_combos)
        if "B" in legs:
            from . import sat
            sat.run_leg(rec, code, targets, a.workdir, a.sat_timeout)
        if "C" in legs:
            from . import leangen
            leangen.run_leg(rec, code, targets, a.workdir, a.lean_timeout)
        certify.verify(rec)
        for ev in rec["evidence"]:
            print(f"leg {ev['leg']} {ev['quantity']:7s} {ev['status']:8s} "
                  f"{ev.get('claim') or ev.get('reason', '')}")
        print(certify.summary(rec))
        if a.out:
            with open(a.out, "w") as f:
                json.dump(rec, f, indent=1)
                f.write("\n")
        return 0
    except certify.Refusal as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 1
    except (codes.SpecError, FileNotFoundError, json.JSONDecodeError) as e:
        print(f"usage/spec error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
