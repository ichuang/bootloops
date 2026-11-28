"""Command line:  python3 -m ftqc_gadgets <command> ...   (with tools/ on PYTHONPATH)

  analyze SPEC [--check KIND:T]... [--close S] [--refines SPEC2:ORDER[:R]]... [--sector]
               [--local [--radius R | --code-receipt CODE.json] [--no-tables]]
               [-o receipt.json] [--cap N]
  equiv   SPEC1 SPEC2 --order W [--input-weight R]
  verify  RECEIPT.json [--full]
  compose TARGET.json ROLE=PART.json... [-o composition.json]
  export  SPEC [-o gadget.json]

SPEC is a gadget JSON file or a library gadget (see ftqc_gadgets.library), for example
cnot:steane, knill:golay,anc=spec, prep:golay,zero,tree=XZXZ, fc:cat:4,verify.
KIND is prep, gate or ec.

Exit codes: 0 ok, 1 refusal or a property that fails, 2 usage or malformed gadget.
"""
import argparse
import json
import sys

from ftqc_codes.certify import Refusal as CodeRefusal
from . import certify, library, props
from .effects import Analysis
from .gadget import GadgetError


def _split_ref(text):
    parts = text.rsplit(":", 2)
    if len(parts) == 3 and parts[1].isdigit() and parts[2].isdigit():
        return parts[0], int(parts[1]), int(parts[2])
    spec, _, order = text.rpartition(":")
    if not spec or not order.isdigit():
        raise GadgetError(f"--refines takes SPEC:ORDER[:R], got {text!r}")
    return spec, int(order), 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ftqc_gadgets", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    staged = sub.add_parser("staged", help="certify a staged encoded gate or conversion")
    staged.add_argument("spec")
    staged.add_argument("-o", "--out")
    installed = sub.add_parser("installed-action", help="verify and audit a physical Clifford route")
    installed.add_argument("spec")
    installed.add_argument("-o", "--out")
    build = sub.add_parser("build-surgery", help="compile a finite constructor to a standalone instrument")
    build.add_argument("spec")
    build.add_argument("-o","--out")
    build.add_argument("--geometry", help="write the finite geometry report")
    surgery = sub.add_parser("surgery", help="certify a signed projective measurement")
    surgery.add_argument("spec")
    surgery.add_argument("--order", type=int, default=1)
    surgery.add_argument("--cap", type=int, default=1000000)
    surgery.add_argument("-o", "--out")
    transfer = sub.add_parser("magic-transfer", help="certify a fresh CSS encoding or occupied-slot SWAP")
    transfer.add_argument("spec")
    transfer.add_argument("-o","--out")
    a = sub.add_parser("analyze")
    a.add_argument("spec")
    a.add_argument("--check", action="append", default=[])
    a.add_argument("--close", type=int)
    a.add_argument("--refines", action="append", default=[])
    a.add_argument("--sector", action="store_true")
    a.add_argument("--local", action="store_true", help="coordinate-locality certificate")
    a.add_argument("--radius", type=int, help="decoder radius for --local (default: from the tables)")
    a.add_argument("--code-receipt", help="ftqc_codes receipt of the block code: decoder radii for --local")
    a.add_argument("--tree", action="store_true", help="factory certificate of a purification tree")
    a.add_argument("--strict", action="store_true",
                   help="for --tree: only the last idle round of a leaf is local (outer-level factories)")
    a.add_argument("--light", action="store_true",
                   help="for --tree alone: do not propagate every fault through the whole circuit (large trees)")
    a.add_argument("--wmin", type=int, help="for --tree: minimum stabilizer weight to prove (large codes)")
    a.add_argument("--no-tables", action="store_true", help="skip decoder tables (large codes)")
    a.add_argument("-o", "--out")
    a.add_argument("--cap", type=int, default=3_000_000)
    q = sub.add_parser("equiv")
    q.add_argument("spec1")
    q.add_argument("spec2")
    q.add_argument("--order", type=int, required=True)
    q.add_argument("--input-weight", type=int, default=0)
    v = sub.add_parser("verify")
    v.add_argument("receipt")
    v.add_argument("--full", action="store_true")
    c = sub.add_parser("compose")
    c.add_argument("target")
    c.add_argument("parts", nargs="+")
    c.add_argument("-o", "--out")
    k = sub.add_parser("cones", help="fault cones of a block-level program (JSON: layers, inputs, outputs)")
    k.add_argument("program")
    k.add_argument("--no-events", action="store_true", help="statistics and a hash instead of the event list")
    k.add_argument("--samples", type=int, default=8, help="random settings of the classical controls")
    k.add_argument("-o", "--out")
    e = sub.add_parser("export")
    e.add_argument("spec")
    e.add_argument("-o", "--out")
    args = ap.parse_args(argv)

    def emit(obj, path):
        txt = json.dumps(obj, indent=1) + "\n"
        if path:
            with open(path, "w") as f:
                f.write(txt)
        return txt

    try:
        if args.cmd == "staged":
            from .staged import analyze
            with open(args.spec) as stream:
                rec = analyze(json.load(stream))
            txt = emit(rec, args.out)
            if not args.out:
                print(txt, end="")
            return int(rec["result"]["status"] != "proved")
        if args.cmd == "installed-action":
            from .installed_action import analyze
            with open(args.spec) as stream:
                result = analyze(json.load(stream))
            txt = emit(result, args.out)
            if not args.out:
                print(txt)
            return 0
        if args.cmd == "magic-transfer":
            from .magic_transfer import transfer_receipt
            with open(args.spec) as f:
                rec = transfer_receipt(json.load(f))
            txt = emit(rec,args.out)
            if not args.out:
                print(txt,end="")
            return 0
        if args.cmd == "build-surgery":
            from .surgery_design import build
            with open(args.spec) as f:
                design = build(json.load(f))
            txt = emit(design["spec"],args.out)
            if not args.out:
                print(txt,end="")
            if args.geometry:
                emit(design["geometry"],args.geometry)
            return 0
        if args.cmd == "surgery":
            from .surgery import receipt
            with open(args.spec) as f:
                r = receipt(json.load(f), args.order, args.cap)
            txt = emit(r, args.out)
            if not args.out:
                print(txt, end="")
            return {"proved":0,"counterexample":1,"inconclusive":3}[r["faults"]["status"]]
        if args.cmd == "cones":
            from . import cones
            with open(args.program) as f:
                d = json.load(f)
            if "layers" not in d or "inputs" not in d or "outputs" not in d:
                raise GadgetError("a program file has layers, inputs and outputs")
            rec = cones.build(cones.Program.from_json(d), samples=args.samples, events=not args.no_events)
            emit(rec, args.out)
            print(cones.summary(rec))
            return 1 if rec["algebra"]["outside"] else 0
        if args.cmd == "verify":
            with open(args.receipt) as f:
                rec = json.load(f)
            if rec.get("schema") == "ftqc_gadgets-staged-1":
                from .staged import verify
                verify(rec)
                print("verified staged gadget")
                return 0
            if rec.get("schema") in ("ftqc-gadgets-fresh-transfer-1","ftqc-gadgets-occupied-swap-1"):
                from .magic_transfer import verify_transfer
                verify_transfer(rec)
                print("verified magic transfer")
                return 0
            if rec.get("schema") == "ftqc_gadgets-surgery-1":
                from .surgery import verify
                verify(rec)
                print("verified signed surgery instrument")
                return 0
            if rec.get("schema") == "ftqc_gadgets-cones-1":
                from . import cones
                for line in cones.verify(rec):
                    print("ok  ", line)
                print(cones.summary(rec))
                return 0
            for line in certify.verify(rec, full=args.full):
                print("ok  ", line)
            print(certify.summary(rec))
            return 0
        if args.cmd == "compose":
            with open(args.target) as f:
                target = json.load(f)
            parts = {}
            for item in args.parts:
                role, _, path = item.partition("=")
                if not path:
                    raise GadgetError("parts are given as ROLE=receipt.json")
                with open(path) as f:
                    parts[role] = json.load(f)
            comp = certify.compose(target, parts)
            emit(comp, args.out)
            print(certify.summary(comp))
            return 0
        if args.cmd == "export":
            txt = emit(library.load(args.spec).to_json(), args.out)
            if not args.out:
                sys.stdout.write(txt)
            return 0
        if args.cmd == "equiv":
            a1, a2 = Analysis(library.load(args.spec1)), Analysis(library.load(args.spec2))
            r12 = props.refine(a1, a2, args.order, args.input_weight)
            r21 = props.refine(a2, a1, args.order, args.input_weight)
            both = min(r12["holds_through"], r21["holds_through"])
            print(f"{a1.g.name} refines {a2.g.name} through order {r12['holds_through']}")
            print(f"{a2.g.name} refines {a1.g.name} through order {r21['holds_through']}")
            for r in (r12, r21):
                if r["violation"]:
                    print("  unmatched:", json.dumps(r["violation"]))
            print(f"fault-equivalent through order {both}" if both >= 0 else "not fault-equivalent at order 0")
            return 0 if both >= args.order else 1
        g = library.load(args.spec)
        checks = []
        for item in args.check:
            kind, _, t = item.partition(":")
            if not t.isdigit():
                raise GadgetError(f"--check takes KIND:T, got {item!r}")
            checks.append((kind, int(t)))
        refs = [(library.load(s), o, r) for s, o, r in map(_split_ref, args.refines)]
        loc = None
        if args.local:
            loc = True if args.radius is None else args.radius
            if args.code_receipt:
                with open(args.code_receipt) as f:
                    loc = library.radii(g, json.load(f))
        tree = dict(strict=args.strict, wmin=args.wmin, light=args.light) if args.tree else None
        rec = certify.build(g, checks, args.close, refs, args.sector, args.cap, local=loc,
                            tables=not args.no_tables, tree=tree)
        emit(rec, args.out)
        print(certify.summary(rec))
        failed = any(r.get("violation") or r.get("local") is False or r.get("tree") is False
                     for r in rec["results"])
        return 1 if failed else 0
    except (certify.Refusal, CodeRefusal) as e:
        print("REFUSED:", e, file=sys.stderr)
        return 1
    except (ValueError, KeyError, TypeError, IndexError, OSError) as e:
        print("error:", e, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
