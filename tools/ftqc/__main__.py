"""Generate or replay a versioned FTQC receipt."""

import argparse
import hashlib
import json
import sys
from pathlib import Path
from .algebra import Limit, seed, hgp_parameters, distance_below
from .extractor import canonical, check_column
from .audit import audit_canonical, audit_hgp
from .circuit import Schedule, search
from .proof import solve, replay_sat

from .t_verify import analyze as analyze_t

FIELDS = {
    "t_gate": ("a b logical_x logical_z rounds target resource","trailing_identity"),
    "hgp": ("a b", ""),
    "column": ("a b column word vertices ports edges target", "paths"),
    "canonical": ("n polynomial edges", "assembly"),
    "circuit": ("a b ops measured_z readout detectors target", ""),
}


def digest(s):
    return hashlib.sha256(
        json.dumps(s, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    ).hexdigest()


def validate(s):
    if (
        not isinstance(s, dict)
        or type(s.get("schema")) != int
        or s["schema"] != 1
        or s.get("kind") not in FIELDS
    ):
        raise ValueError("unknown input schema or kind")
    required, optional = FIELDS[s["kind"]]
    keys = {"kind", "schema", *required.split()}
    if not keys <= s.keys() or s.keys() - keys - set(optional.split()):
        raise ValueError("missing or unsupported input fields")
    if "target" in s and (type(s["target"]) != int or not 1 <= s["target"] <= 10000):
        raise ValueError("target must be an integer in [1,10000]")


def evaluate(s, backend, budget, seconds, reverse=False):
    validate(s)
    if (
        type(budget) != int
        or budget < 0
        or seconds <= 0
        or backend not in ("exact", "sat")
    ):
        raise ValueError("invalid search options")
    kind = s["kind"]
    if backend == "sat" and kind not in ("circuit","t_gate"):
        raise ValueError("SAT backend requires circuit input")
    try:
        if kind == "t_gate":
            return analyze_t(s,backend,budget,seconds,reverse)
        if kind == "hgp":
            result = dict(
                status="proved", **hgp_parameters(*seed(s["a"]), *seed(s["b"]))
            )
            if reverse:
                audit_hgp(s, result)
            return result
        if kind == "canonical":
            result = canonical(s)
            if reverse and result["status"] == "proved":
                audit_canonical(s, result)
            return result
        if kind == "column":
            return check_column(s, budget, reverse)
        circuit = Schedule(s)
        return (
            solve(circuit, s["target"], seconds)
            if backend == "sat"
            else search(circuit, s["target"], budget, reverse)
        )
    except Limit as err:
        return dict(status="inconclusive", reason=str(err))


def receipt(s, backend="exact", budget=2000000, seconds=60):
    return dict(
        version=1,
        input_sha256=digest(s),
        backend=backend,
        budget=budget,
        result=evaluate(s, backend, budget, seconds),
    )


def verify(s, r):
    validate(s)
    if (
        set(r) != {"version", "input_sha256", "backend", "budget", "result"}
        or r["version"] != 1
        or r["input_sha256"] != digest(s)
    ):
        raise ValueError("receipt header or input digest mismatch")
    result = r["result"]
    if result["status"] not in ("proved", "counterexample"):
        raise ValueError("receipt contains an open obligation")
    if r["backend"] == "sat" and s['kind']=='t_gate':
        if digest(result)!=digest(analyze_t(s,'sat',r['budget'],60,True,result)):
            raise ValueError("T-gate SAT receipt mismatch")
        return
    if r["backend"] == "sat":
        if s["kind"] != "circuit":
            raise ValueError("wrong backend")
        replay_sat(Schedule(s), s["target"], result)
    elif r["backend"] != "exact" or digest(
        evaluate(s, "exact", r["budget"], 60, True)
    ) != digest(result):
        raise ValueError("independent replay mismatch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("check")
    run.add_argument("spec", type=Path)
    run.add_argument("-o", type=Path)
    run.add_argument("--backend", choices=["exact", "sat"], default="exact")
    run.add_argument("--budget", type=int, default=2000000)
    run.add_argument("--seconds", type=int, default=60)
    check = sub.add_parser("verify")
    check.add_argument("spec", type=Path)
    check.add_argument("receipt", type=Path)
    args = parser.parse_args()
    try:
        s = json.loads(args.spec.read_text())
        if args.command == "verify":
            verify(s, json.loads(args.receipt.read_text()))
            print("verified")
            return 0
        r = receipt(s, args.backend, args.budget, args.seconds)
        text = json.dumps(r, sort_keys=True, indent=2) + "\n"
        if args.o:
            args.o.write_text(text)
        else:
            print(text, end="")
        return {"proved": 0, "counterexample": 1, "inconclusive": 3}[
            r["result"]["status"]
        ]
    except (OSError, ValueError, KeyError, TypeError, IndexError, OverflowError) as err:
        print(f"FTQC refusal: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
