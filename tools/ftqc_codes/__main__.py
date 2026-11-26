"""Command line:  python3 -m ftqc_codes <command> ...   (with tools/ on PYTHONPATH)

  certify SPEC [-o receipt.json] [--legs A,B,C,D,E] [--target q=v ...] [--max-combos N]
  verify  RECEIPT.json [--full]
  spec    SPEC [-o code.json]          write the code's matrices as JSON
  families                             list family constructors
SPEC is a code JSON file or family[:args], e.g. golay, qr:47, toric:5, bb72.

Exit codes: 0 ok, 1 refusal (not a CSS code, receipt fails, legs disagree),
2 usage or malformed spec.
"""
import argparse
import json
import sys

from . import certify, codes


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
    v = sub.add_parser("verify")
    v.add_argument("receipt")
    v.add_argument("--full", action="store_true")
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
        if a.cmd == "verify":
            with open(a.receipt) as f:
                rec = json.load(f)
            for line in certify.verify(rec, full=a.full):
                print("ok  ", line)
            print(certify.summary(rec))
            return 0
        legs = [x.strip().upper() for x in a.legs.split(",") if x.strip()]
        if set(legs) - {"A", "B", "C", "D", "E"}:
            raise codes.SpecError("--legs takes a subset of A,B,C,D,E")
        code = codes.load(a.spec)
        targets = _targets(a.target)
        rec = certify.build(code, targets, a.max_combos, legs)
        if "D" in legs:
            from . import bch
            bch.run_leg(rec, code)
        if "E" in legs:
            from . import cuts
            cuts.run_leg(rec, code)
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
