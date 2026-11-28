"""Command line:  python3 -m ftqc_faults audit <command> ...   (with tools/ on PYTHONPATH)

  analyze SPEC [-o receipt.json] [--max-weight W] [--strict S] [--heavier 1,2] [--legs A,B,S]
               [--cap N] [--workdir DIR] [--sat-timeout T]
  verify  RECEIPT.json [--full]
  export  SPEC [-o circuit.json] [--stim]
SPEC is a circuit JSON file or a library gadget: cat:4, cat:4,verify, css_prep:golay,
surface:3,2, surface:3,2,z,bad.

Exit codes: 0 ok, 1 refusal (receipt fails, legs disagree), 2 usage or malformed circuit.
"""
import argparse
import json
import sys

from . import certify, library, stimleg
from .circuit import CircuitError


def main(argv=None, *, prog="ftqc_faults audit"):
    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser(prog=prog, description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("analyze")
    a.add_argument("spec")
    a.add_argument("-o", "--out")
    a.add_argument("--max-weight", type=int, default=4)
    a.add_argument("--strict", type=int, default=0)
    a.add_argument("--heavier", default="", help="comma-separated weights w")
    a.add_argument("--legs", default="A")
    a.add_argument("--cap", type=int, default=3_000_000)
    a.add_argument("--workdir")
    a.add_argument("--sat-timeout", type=int, default=900)
    v = sub.add_parser("verify")
    v.add_argument("receipt")
    v.add_argument("--full", action="store_true")
    e = sub.add_parser("export")
    e.add_argument("spec")
    e.add_argument("-o", "--out")
    e.add_argument("--stim", action="store_true")
    r = sub.add_parser(
        "records", help="check fault signatures of a CSS record-only circuit"
    )
    r.add_argument("spec")
    r.add_argument(
        "--stim", action="store_true", help="compare all detector/observable mechanisms"
    )
    r.add_argument("--fault", action="append", default=[], metavar="LOCATION:LABEL")
    r.add_argument("-o", "--out")
    injection = sub.add_parser("injection", help="forward a published injection study to ftqc_extractors")
    if argv[:1] == ["injection"]:
        from ftqc_extractors.injection_cli import configure_parser
        configure_parser(injection)
    args = ap.parse_args(argv)
    try:
        if args.cmd=="injection":
            from ftqc_extractors.injection_cli import command
            return command(args)
        if args.cmd == "verify":
            with open(args.receipt) as f:
                rec = json.load(f)
            for line in certify.verify(rec, full=args.full):
                print("ok  ", line)
            print(certify.summary(rec))
            return 0
        circ = library.load(args.spec)
        if args.cmd == "records":
            from .record_model import record_report

            result = record_report(circ, args.fault, args.stim)
            out = json.dumps(result, indent=2, sort_keys=True) + "\n"
            if args.out:
                with open(args.out, "w") as f:
                    f.write(out)
            else:
                sys.stdout.write(out)
            return 0
        if args.cmd == "export":
            if args.stim:
                from .model import Model

                txt = stimleg.to_stim(circ, Model(circ))
            else:
                txt = json.dumps(circ.to_json()) + "\n"
            if args.out:
                with open(args.out, "w") as f:
                    f.write(txt)
            else:
                sys.stdout.write(txt)
            return 0
        legs = [x.strip().upper() for x in args.legs.split(",") if x.strip()]
        if set(legs) - {"A", "B", "S"}:
            raise CircuitError("--legs takes a subset of A,B,S")
        heavier = tuple(int(x) for x in args.heavier.split(",") if x.strip())
        rec, model = certify.build(circ, args.max_weight, args.strict, legs, args.cap,
                                   args.workdir, args.sat_timeout, heavier=heavier)
        if "S" in legs:
            lo = rec["results"].get("circuit_distance", {}).get("lower")
            ev = stimleg.compare(circ, model, lo)
            rec["results"]["stim"] = ev
            if ev["status"] == "DISAGREES":
                raise certify.Refusal("Stim disagrees: " + ev["reason"])
        certify.verify(rec)
        for kind in ("circuit_distance", "fault_distance"):
            for ev in rec["results"].get(kind, {}).get("evidence", []):
                print(f"leg {ev['leg']} {kind:17s} {ev['status']:8s} "
                      f"{ev.get('claim') or ev.get('reason') or ''} {ev.get('note') or ''}")
        if "stim" in rec["results"]:
            ev = rec["results"]["stim"]
            print(f"leg S {'fault model':17s} {ev['status']:8s} {ev.get('reason', '')}")
        for row in rec["results"].get("heavier_than", []):
            val = row["lower"] if row["exact"] else f">= {row['lower']}"
            print(f"accepted error heavier than weight {row['weight']} needs {val} faults")
        print(certify.summary(rec))
        if args.out:
            with open(args.out, "w") as f:
                json.dump(rec, f)
                f.write("\n")
        return 0
    except certify.Refusal as ex:
        print(f"REFUSED: {ex}", file=sys.stderr)
        return 1
    except (CircuitError, FileNotFoundError, json.JSONDecodeError, ValueError) as ex:
        print(f"usage/circuit error: {ex}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
