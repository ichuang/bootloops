"""Command-line interface for `ftqc_uncompute`.

Exit codes:
  - `0` : command / verification / selftest succeeded
  - `1` : soundness refusal, failed verification, or failed selftest
  - `2` : CLI usage or malformed input error
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

from .adapters.qiskit_qasm3 import kmx_from_qasm3, kmx_to_qasm3
from .certify import (
    ReceiptRefusalError,
    build_receipt_for_kmx,
    build_receipt_for_spec_dict,
    parse_strict_json,
    verify_receipt_file,
)
from .ir_dag import ArithDAG
from .ir_kickmix import KickMixCircuit
from .masking import verify_mask_spec
from .pebbling import (
    PebbleDAG,
    compile_gidney_sequential_schedule,
    compile_kmrv_blast_schedule,
    solve_optimal_pebbling_astar,
    validate_pebbling_schedule,
)
from .suggest import (
    apply_dag_rewrites,
    apply_kmx_rewrites,
    suggest_for_dag,
    suggest_for_kmx,
)


def _write_or_print(text: str, out_path: str | None) -> None:
    if out_path:
        Path(out_path).write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
        if not text.endswith("\n"):
            sys.stdout.write("\n")


def _cmd_check(ns: argparse.Namespace) -> int:
    in_path = Path(ns.target)
    if not in_path.exists():
        sys.stderr.write(f"error: target file not found: {in_path}\n")
        return 2
    try:
        if in_path.suffix.lower() == ".kmx":
            circ = KickMixCircuit.from_kmx(in_path.read_text(encoding="utf-8"))
            receipt = build_receipt_for_kmx(
                circ,
                source_name=str(in_path),
                shots=ns.shots,
                run_kraus=True if ns.kraus else None,
            )
        else:
            raw = in_path.read_bytes()
            spec = parse_strict_json(raw.decode("utf-8"))
            receipt = build_receipt_for_spec_dict(spec, raw_bytes=raw, source_name=str(in_path))
    except (ValueError, KeyError, ReceiptRefusalError) as exc:
        sys.stderr.write(f"error: {exc}\n")
        return 2

    out_json = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    _write_or_print(out_json, ns.output)
    return 0 if receipt["verdict"] == "CERTIFIED" else 1


def _cmd_verify(ns: argparse.Namespace) -> int:
    rpath = Path(ns.receipt)
    if not rpath.exists():
        sys.stderr.write(f"error: receipt file not found: {rpath}\n")
        return 2
    try:
        ok, errs = verify_receipt_file(rpath, target_path=ns.target)
    except ReceiptRefusalError as exc:
        sys.stderr.write(f"REFUSED ({exc.code}): {exc.message}\n")
        return 1
    except Exception as exc:
        sys.stderr.write(f"error: {exc}\n")
        return 2

    report = {"verified": ok, "errors": errs, "receipt": str(rpath)}
    sys.stdout.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return 0 if ok else 1


def _cmd_suggest(ns: argparse.Namespace) -> int:
    in_path = Path(ns.target)
    if not in_path.exists():
        sys.stderr.write(f"error: target file not found: {in_path}\n")
        return 2
    zx_contract: dict[str, Any] | None = None
    if ns.zx_contract:
        zx_contract = parse_strict_json(Path(ns.zx_contract).read_text(encoding="utf-8"))

    if in_path.suffix.lower() == ".kmx":
        circ = KickMixCircuit.from_kmx(in_path.read_text(encoding="utf-8"))
        if ns.apply:
            rewritten, sugs = apply_kmx_rewrites(circ, verify=True)
            if ns.output:
                Path(ns.output).write_text(rewritten.to_kmx(), encoding="utf-8")
            payload = {
                "applied": bool(sugs),
                "before_toffolis": circ.operation_statistics().num_toffoli_like,
                "after_toffolis": rewritten.operation_statistics().num_toffoli_like,
                "before_qubits": circ.num_qubits,
                "after_qubits": rewritten.num_qubits,
                "suggestions": [s.to_dict() for s in sugs],
            }
            sys.stdout.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
            return 0
        sugs = suggest_for_kmx(circ, zx_contract=zx_contract)
        sys.stdout.write(json.dumps([s.to_dict() for s in sugs], indent=2, sort_keys=True) + "\n")
        return 0

    dag = ArithDAG.from_json(in_path.read_text(encoding="utf-8"))
    if ns.apply:
        rew_dag, sugs = apply_dag_rewrites(dag, allow_approximate_mask=ns.allow_approx_mask)
        if ns.output:
            Path(ns.output).write_text(rew_dag.to_json() + "\n", encoding="utf-8")
        payload = {
            "applied": bool(sugs),
            "before_costs": dag.estimate_costs(),
            "after_costs": rew_dag.estimate_costs(),
            "suggestions": [s.to_dict() for s in sugs],
        }
        sys.stdout.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return 0
    sugs = suggest_for_dag(dag, zx_contract=zx_contract)
    sys.stdout.write(json.dumps([s.to_dict() for s in sugs], indent=2, sort_keys=True) + "\n")
    return 0


def _cmd_pebble(ns: argparse.Namespace) -> int:
    length = ns.length
    mode = ns.mode
    if mode == "kmrv":
        steps = compile_kmrv_blast_schedule(length)
    elif mode == "gidney":
        steps = compile_gidney_sequential_schedule(length)
    elif mode == "astar":
        if ns.pebbles is None:
            sys.stderr.write("error: --pebbles required for --mode astar\n")
            return 2
        res = solve_optimal_pebbling_astar(length, ns.pebbles, window=ns.window)
        if res is None:
            sys.stderr.write(f"INFEASIBLE: cannot pebble length={length} within space={ns.pebbles}\n")
            return 1
        steps = res
    else:
        return 2

    dag = PebbleDAG.line(length, window=ns.window)
    val = validate_pebbling_schedule(dag, steps)
    out_obj = {
        "kind": "pebbling_schedule",
        "mode": mode,
        "length": length,
        "window": ns.window,
        "validation": val.to_dict(),
        "steps": [s.to_dict() for s in steps],
    }
    _write_or_print(json.dumps(out_obj, indent=2, sort_keys=True) + "\n", ns.output)
    return 0 if val.valid else 1


def _cmd_mask(ns: argparse.Namespace) -> int:
    mults = [int(x) for x in ns.multipliers.split(",") if x.strip()] if ns.multipliers else []
    res = verify_mask_spec(
        ns.kind,
        modulus=ns.modulus,
        width=ns.width,
        multipliers=mults,
        runs=ns.runs,
        min_prime_factor=ns.min_prime_factor,
        runway_bits=ns.runway_bits,
        max_added_delta=ns.max_delta,
        subroutine_kind=ns.subroutine,
        subroutine_steps=ns.steps,
    )
    _write_or_print(json.dumps(res.to_dict(), indent=2, sort_keys=True) + "\n", ns.output)
    return 0 if res.valid else 1


def _cmd_convert(ns: argparse.Namespace) -> int:
    in_path = Path(ns.input)
    if not in_path.exists():
        sys.stderr.write(f"error: input file not found: {in_path}\n")
        return 2
    txt = in_path.read_text(encoding="utf-8")
    in_suf = in_path.suffix.lower()
    out_suf = Path(ns.output).suffix.lower() if ns.output else ".qasm"

    if in_suf == ".kmx":
        circ = KickMixCircuit.from_kmx(txt)
    elif in_suf in (".qasm", ".qasm3"):
        circ = kmx_from_qasm3(txt)
    elif in_suf == ".json":
        circ = ArithDAG.from_json(txt).lower_to_kmx()
    else:
        sys.stderr.write(f"error: unsupported input extension {in_suf!r}\n")
        return 2

    if out_suf in (".qasm", ".qasm3"):
        rendered = kmx_to_qasm3(circ)
    elif out_suf == ".kmx":
        rendered = circ.to_kmx()
    else:
        sys.stderr.write(f"error: unsupported output extension {out_suf!r}\n")
        return 2

    _write_or_print(rendered, ns.output)
    return 0


def _cmd_blocks(ns: argparse.Namespace) -> int:
    from . import blocks

    in_path = Path(ns.spec)
    if not in_path.exists():
        sys.stderr.write(f"error: spec file not found: {in_path}\n")
        return 2
    try:
        spec = json.loads(in_path.read_text(encoding="utf-8"))
        if ns.verify:
            blocks.verify_receipt(spec, json.loads(Path(ns.verify).read_text(encoding="utf-8")))
            print("verified receipt")
            return 0
        rec = blocks.receipt(spec)
    except (ValueError, KeyError, TypeError) as exc:
        sys.stderr.write(f"error: {exc}\n")
        return 2 if not ns.verify else 1
    _write_or_print(json.dumps(rec, indent=2, sort_keys=True), ns.output)
    return 0 if rec["result"]["status"] == "proved" else 1


def _cmd_selftest(_ns: argparse.Namespace) -> int:
    from .selftest import run_selftest

    return run_selftest()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="ftqc_uncompute", description="MBUC, spooky pebbling & superposition masking toolkit")
    sub = p.add_subparsers(dest="command", required=True)

    p_chk = sub.add_parser("check", help="Validate a .kmx circuit or JSON spec and emit a receipt")
    p_chk.add_argument("target")
    p_chk.add_argument("--shots", type=int, default=512)
    p_chk.add_argument("--kraus", action="store_true")
    p_chk.add_argument("-o", "--output", default=None)
    p_chk.set_defaults(func=_cmd_check)

    p_ver = sub.add_parser("verify", help="Verify an ftqc_uncompute-receipt-1 JSON receipt")
    p_ver.add_argument("receipt")
    p_ver.add_argument("--target", default=None)
    p_ver.set_defaults(func=_cmd_verify)

    p_sug = sub.add_parser("suggest", help="Suggest (and optionally --apply) R1-R7 simplifications")
    p_sug.add_argument("target")
    p_sug.add_argument("--apply", action="store_true")
    p_sug.add_argument("--allow-approx-mask", action="store_true")
    p_sug.add_argument("--zx-contract", default=None)
    p_sug.add_argument("-o", "--output", default=None)
    p_sug.set_defaults(func=_cmd_suggest)

    p_peb = sub.add_parser("pebble", help="Synthesize and validate spooky pebbling schedules")
    p_peb.add_argument("--length", type=int, required=True)
    p_peb.add_argument("--mode", choices=["kmrv", "gidney", "astar"], default="kmrv")
    p_peb.add_argument("--pebbles", type=int, default=None)
    p_peb.add_argument("--window", type=int, default=None)
    p_peb.add_argument("-o", "--output", default=None)
    p_peb.set_defaults(func=_cmd_pebble)

    p_msk = sub.add_parser("mask", help="Validate superposition masks and compute Decimal error bounds")
    p_msk.add_argument("--kind", required=True)
    p_msk.add_argument("--modulus", type=int, default=None)
    p_msk.add_argument("--width", type=int, default=None)
    p_msk.add_argument("--multipliers", default="")
    p_msk.add_argument("--runs", type=int, default=1)
    p_msk.add_argument("--min-prime-factor", type=int, default=None)
    p_msk.add_argument("--runway-bits", type=int, default=None)
    p_msk.add_argument("--max-delta", type=int, default=1)
    p_msk.add_argument("--subroutine", default=None)
    p_msk.add_argument("--steps", type=int, default=None)
    p_msk.add_argument("-o", "--output", default=None)
    p_msk.set_defaults(func=_cmd_mask)

    p_cnv = sub.add_parser("convert", help="Convert between .kmx, .qasm (OpenQASM 3), and ArithDAG .json")
    p_cnv.add_argument("input")
    p_cnv.add_argument("-o", "--output", default=None)
    p_cnv.set_defaults(func=_cmd_convert)

    p_blk = sub.add_parser("blocks", help="Certify a block program: deletions, integer relation, binding to layers")
    p_blk.add_argument("spec")
    p_blk.add_argument("--verify", default=None, help="a receipt to recompute and compare")
    p_blk.add_argument("-o", "--output", default=None)
    p_blk.set_defaults(func=_cmd_blocks)

    p_st = sub.add_parser("selftest", help="Run cold-clone positive, sabotage, and mutant battery")
    p_st.set_defaults(func=_cmd_selftest)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    try:
        ns = parser.parse_args(argv)
    except SystemExit as exc:
        return int(exc.code) if isinstance(exc.code, int) else 2
    return int(ns.func(ns))
