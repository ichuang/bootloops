"""Command-line interface for `ftqc_magic`."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Sequence

from .catalog import get_catalog_factory, list_catalog
from .certify import build_receipt, verify_receipt
from .cultivate import (
    analyze_cultivation_double_check,
    build_itogawa_short_check_schedule,
    build_msc3_double_check_schedule,
)
from .ir import MagicRefusalError
from .matrix import BinaryFactory
from .pipeline import suggest_pipelines
from .recycle import compress_active_qubits


def _write_or_print(obj: object, out_path: str | None) -> None:
    txt = json.dumps(obj, indent=2, sort_keys=True) + "\n"
    if out_path:
        Path(out_path).write_text(txt, encoding="utf-8")
    else:
        sys.stdout.write(txt)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="ftqc_magic",
        description="Certified magic-state distillation factories, cultivation, and pipelines.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)
    switched = sub.add_parser("switched-source", help="certify fresh encoded plus through a switched T")
    switched.add_argument("study", type=Path)
    switched.add_argument("-o", "--out")

    p_check = sub.add_parser("check", help="Build a certified ftqc_magic-receipt-1 from a JSON spec")
    p_check.add_argument("spec", help="Path to input JSON specification")
    p_check.add_argument("-o", "--output", help="Optional path to write receipt JSON")
    p_check.add_argument("--external", action="store_true", help="Run live CaDiCaL / Stim / ftqc_faults checks")

    p_ver = sub.add_parser("verify", help="Independently replay and verify an ftqc_magic-receipt-1")
    p_ver.add_argument("spec", help="Path to input JSON specification")
    p_ver.add_argument("receipt", help="Path to receipt JSON")
    p_ver.add_argument("--ignore-source-sha", action="store_true", help="Allow verification across source revisions")

    p_rec = sub.add_parser("recycle", help="Minimize active logical qubits via mid-circuit check recycling")
    p_rec.add_argument("--factory", required=True, help="Catalog factory name (e.g. bk_15_1_3_t, bh_49_1_5_t)")
    p_rec.add_argument("--lrat", action="store_true", help="Certify column-ordering lower bound via CaDiCaL LRAT")
    p_rec.add_argument("-o", "--output", help="Optional output JSON path")

    p_cult = sub.add_parser("cultivate", help="Run exact O(n^3) Clifford-hook cultivation fault-distance analysis")
    p_cult.add_argument(
        "--variant",
        default="msc3_z_flagged",
        choices=["msc3_unflagged", "msc3_z_flagged", "itogawa_2round"],
    )
    p_cult.add_argument("--p", default="1/1000", help="Exact rational physical error rate")
    p_cult.add_argument("--external", action="store_true", help="Cross-check S-proxy against ftqc_faults and Stim")
    p_cult.add_argument("-o", "--output", help="Optional output JSON path")

    p_pipe = sub.add_parser("pipeline", help="Alias for `check` on a multi-stage pipeline JSON spec")
    p_pipe.add_argument("spec", help="Path to pipeline JSON specification")
    p_pipe.add_argument("-o", "--output", help="Optional path to write receipt JSON")
    p_pipe.add_argument("--external", action="store_true", help="Run live external cross-checks")

    p_sug = sub.add_parser("suggest", help="Return Pareto-optimal hybrid magic-state pipelines")
    p_sug.add_argument("--target", default="T", choices=["T", "CCZ"])
    p_sug.add_argument("--p", default="1/1000", help="Exact rational physical error rate")
    p_sug.add_argument("--infidelity", default="1/10000000000", help="Maximum conditional output infidelity")
    p_sug.add_argument("--max-qubits", type=int, default=None, help="Maximum active logical/physical qubits")
    p_sug.add_argument("-o", "--output", help="Optional output JSON path")

    sub.add_parser("catalog", help="List built-in landmark magic-state factories")

    args = parser.parse_args(argv)
    try:
        if args.cmd == "switched-source":
            from .staged import prepare
            _write_or_print(prepare(json.loads(args.study.read_text())), args.out)
            return 0
        if args.cmd in ("check", "verify"):
            spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
            if spec.get("kind") == "switched-source":
                from .staged import source_receipt
                from .certify import canonical_json_sha256
                result = source_receipt(spec)
                if args.cmd == "check":
                    _write_or_print(result, args.output)
                elif canonical_json_sha256(result) != canonical_json_sha256(json.loads(Path(args.receipt).read_text())):
                    raise ValueError("switched source receipt differs from replay")
                else:
                    print("verified switched source")
                return 0
            if spec.get("kind") in ("insitu-rm15", "rm15-action"):
                # Native Extractors owns source replay; share the same joined
                # receipt with Compose's logical gate and hardware consumers.
                from ftqc_compose.__main__ import receipt as replay_inputs
                from .certify import canonical_json_sha256
                from ftqc_codes.certify import Refusal as CodeRefusal
                try:
                    joined = replay_inputs(spec)
                except CodeRefusal as exc:
                    raise ValueError(str(exc)) from exc
                if args.cmd == "check":
                    _write_or_print(joined, args.output)
                elif canonical_json_sha256(joined) != canonical_json_sha256(json.loads(Path(args.receipt).read_text())):
                    raise ValueError("RM15 receipt differs from replay")
                else:
                    print("verified", spec["kind"])
                return 0
        if args.cmd in ("check", "pipeline"):
            spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
            rec = build_receipt(spec, crosscheck_external=bool(args.external))
            _write_or_print(rec, args.output)
            return 0
        if args.cmd == "verify":
            spec = json.loads(Path(args.spec).read_text(encoding="utf-8"))
            rec = json.loads(Path(args.receipt).read_text(encoding="utf-8"))
            res = verify_receipt(spec, rec, strict_source_sha=not bool(args.ignore_source_sha))
            _write_or_print(res, None)
            return 0
        if args.cmd == "recycle":
            fac = get_catalog_factory(args.factory)
            if not isinstance(fac, BinaryFactory):
                raise ValueError(f"{args.factory!r} is not a binary PPR factory")
            sched = compress_active_qubits(fac, certify_lrat=bool(args.lrat))
            _write_or_print(sched.to_dict(), args.output)
            return 0
        if args.cmd == "cultivate":
            if args.variant == "msc3_unflagged":
                s = build_msc3_double_check_schedule(z_flagged=False)
            elif args.variant == "msc3_z_flagged":
                s = build_msc3_double_check_schedule(z_flagged=True)
            else:
                s = build_itogawa_short_check_schedule()
            c, meta = analyze_cultivation_double_check(
                s, p_phys=args.p, crosscheck_faults_and_stim=bool(args.external)
            )
            _write_or_print({"contract": c.to_dict(), "analysis": meta}, args.output)
            return 0
        if args.cmd == "suggest":
            front = suggest_pipelines(
                target_gate=args.target,
                p_phys=args.p,
                target_infidelity=args.infidelity,
                max_active_qubits=args.max_qubits,
            )
            _write_or_print({"target": args.target, "pareto_frontier": front}, args.output)
            return 0
        if args.cmd == "catalog":
            _write_or_print({"catalog": list_catalog()}, None)
            return 0
    except MagicRefusalError as exc:
        sys.stderr.write(json.dumps(exc.to_dict(), sort_keys=True) + "\n")
        return 2
    except (ValueError, KeyError, TypeError, OSError) as exc:
        sys.stderr.write(f"ftqc_magic: {exc}"+chr(10))
        return 2
    return 0
