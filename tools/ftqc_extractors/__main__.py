"""Generate or replay a versioned ftqc_extractors receipt."""

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

from .t_source import source_record

from .t_bounds import structural

from .gross_distance import gross_certificate

FIELDS = {
    "injection": ("method memory logicals", "preparation"),
    "pbc": ("blocks bridges magic_supplies gates target max_degree", "outcomes"),
    "surgery_zx": ("memory operation max_degree", ""),
    "modular_expansion": ("vertices edges components layers target", ""),
    "intermediate_reading": ("memory operations padding max_degree target rounds", ""),
    "gross_block_reading": ("axis padding target", ""),
    "surgery_sequence": ("memory operations padding max_degree target", ""),
    "fixed_graph": (
        "n checks vertices edges ports cycles routes max_degree logical target",
        "couplers",
    ),
    "block_reading": ("memory blocks pattern axis target", "expected"),
    "gross_circuit": ("rounds target", "spacetime fault_witnesses"),
    "gross_z": ("target", ""),
    "t_bounds": (
        "a b logical_x logical_z rounds resource",
        "verification_rounds representatives",
    ),
    "t_source": ("a b logical_x logical_z rounds target resource", ""),
    "t_gate": (
        "a b logical_x logical_z rounds target resource",
        "trailing_identity verification_rounds representatives",
    ),
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
    if backend == "sat" and kind not in ("circuit", "t_gate"):
        raise ValueError("SAT backend requires circuit input")
    try:
        if kind == "injection":
            from .injection import analyze

            return analyze(s)
        if kind == "fixed_graph":
            from .fixed_graph import analyze

            return analyze(s, budget)
        if kind == "surgery_sequence":
            from .surgery_sequence import analyze

            return analyze(s, budget)
        if kind == "gross_block_reading":
            from .reading_audit import analyze

            return analyze(s, budget)
        if kind == "intermediate_reading":
            from .surgery_bounds import intermediate

            return intermediate(s, budget)
        if kind == "modular_expansion":
            from .surgery_bounds import expansion

            return expansion(s, budget)
        if kind == "surgery_zx":
            from .surgery_zx import analyze

            return analyze(s, budget)
        if kind == "pbc":
            from .pbc import analyze

            return analyze(s, budget)
        if kind == "block_reading":
            from .block_reading import analyze

            return analyze(s, budget)
        if kind == "gross_circuit":
            from .gross_noisy import report

            return report(s, budget)
        if kind == "gross_z":
            return gross_certificate(s, budget)
        if kind == "t_bounds":
            return structural(s)
        if kind == "t_source":
            return source_record(s, reverse)
        if kind == "t_gate":
            return analyze_t(s, backend, budget, seconds, reverse)
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
    if s["kind"] == "gross_z":
        if r["backend"] != "exact" or type(r["budget"]) is not int or r["budget"] < 0:
            raise ValueError("invalid gross-code receipt backend or budget")
        if digest(result) != digest(gross_certificate(s, r["budget"], result)):
            raise ValueError("gross-code distance receipt mismatch")
        return
    if r["backend"] == "sat" and s["kind"] == "t_gate":
        if digest(result) != digest(analyze_t(s, "sat", r["budget"], 60, True, result)):
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
    export = sub.add_parser(
        "schedule", help="expand a conditional T gadget to physical fault locations"
    )
    export.add_argument("spec", type=Path)
    export.add_argument(
        "--probe", choices=["X", "Z"], help="gross circuit readout basis"
    )
    export.add_argument("--arm", choices=("gate", "X", "Y"), default="gate")
    export.add_argument("-o", type=Path)
    meas = sub.add_parser(
        "measurement",
        help="export the ideal gross-code merge/split protocol and merged checks",
    )
    meas.add_argument("spec", type=Path)
    meas.add_argument("-o", type=Path)
    check = sub.add_parser("verify")
    check.add_argument("spec", type=Path)
    check.add_argument("receipt", type=Path)
    args = parser.parse_args()
    try:
        s = json.loads(args.spec.read_text())
        if args.command == "measurement":
            from .gross import measurement

            validate(s)
            if s["kind"] != "gross_z":
                raise ValueError(
                    "measurement export requires an unprimed gross Z input"
                )
            data = dict(input_sha256=digest(s), **measurement())
            text = json.dumps(data, indent=2) + "\n"
            if args.o:
                args.o.write_text(text)
            else:
                print(text, end="")
            return 0
        if args.command == "schedule":
            from .t_schedule import TSchedule
            from .algebra import bits

            validate(s)
            if s["kind"] == "gross_block_reading":
                from ftqc_codes.codes import CSSCode, load
                from .reading_audit import probe

                if (
                    args.probe is None
                    or args.arm != "gate"
                    or s["axis"] not in ("X", "Z")
                    or s["target"] != 12
                    or s["padding"] != [12, 12]
                ):
                    raise ValueError("select a canonical BB144 reading probe X or Z")
                mem = load("bb144")
                if s["axis"] == "X":
                    mem = CSSCode("dual", mem.n, mem.hz, mem.hx)
                text = json.dumps(probe(mem, s["axis"], args.probe)[0].to_json()) + "\n"
                if args.o:
                    args.o.write_text(text)
                else:
                    print(text, end="")
                return 0
            if s["kind"] == "gross_circuit":
                from .gross_probes import experiment

                if args.probe is None or args.arm != "gate":
                    raise ValueError("gross circuit export requires --probe X or Z")
                data = experiment(s["rounds"], args.probe).to_json()
                text = json.dumps(data, separators=(",", ":")) + "\n"
                if args.o:
                    args.o.write_text(text)
                else:
                    print(text, end="")
                return 0
            if args.probe is not None:
                raise ValueError("--probe applies to gross_circuit schedules")
            if (
                s["kind"] not in ("t_gate", "t_bounds")
                or s["resource"] != "ideal-encoded-T"
            ):
                raise ValueError("schedule export requires a conditional T input")
            circuit = TSchedule(s, {"gate": None, "X": 0, "Y": 1}[args.arm])
            data = dict(
                input_sha256=digest(s),
                arm=args.arm,
                operations=circuit.ops,
                sign_by_reported_x=circuit.read_sign,
                meters=[
                    dict(
                        start=start,
                        stop=end,
                        x_support=list(bits(x)),
                        z_support=list(bits(z)),
                        readout=list(bits(r)),
                        flags=flags,
                    )
                    for start, end, x, z, r, flags in circuit.parts
                ],
            )
            text = json.dumps(data, indent=2) + "\n"
            if args.o:
                args.o.write_text(text)
            else:
                print(text, end="")
            return 0
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
        print(f"ftqc_extractors refusal: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
