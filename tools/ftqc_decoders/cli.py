"""Command-line interface for `ftqc_decoders`.

Exit codes (BootLoops convention):
  * `0`: PASS / certified / receipt verified
  * `1`: Counterexample / trapping-set failure witness found / verification mismatch
  * `2`: Refusal (`RefusalError`, invalid input, float input, missing optional backend)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

from .adapters.external_decoders import probe_external_backends
from .certify import build_receipt, load_problem_from_spec, verify_receipt
from .decimal_bounds import RefusalError


def _parse_code_shorthand(code_arg: str, prior: str = "1/1000", sector: str = "X") -> dict[str, Any]:
    if code_arg == "bb72":
        return {"kind": "bb72", "sector": sector, "prior": prior}
    if code_arg == "bb144":
        return {"kind": "bb144", "sector": sector, "prior": prior}
    if code_arg.startswith("toric:"):
        return {"kind": "toric", "length": int(code_arg.split(":", 1)[1]), "prior": prior}
    if code_arg.startswith("repetition:"):
        return {"kind": "repetition", "distance": int(code_arg.split(":", 1)[1]), "prior": prior}
    if code_arg.endswith(".dem") or Path(code_arg).is_file():
        return {"kind": "stim_dem", "path": code_arg, "prior": prior}
    return {"kind": "css_code", "code": code_arg, "sector": sector, "prior": prior}


def _emit_json(data: Any, out_path: Path | None) -> None:
    text = json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if out_path is not None:
        out_path.write_text(text)
    else:
        sys.stdout.write(text)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="ftqc_decoders",
        description="Certified zero-float QEC decoders, trapping-set / d_eff verifiers, Stim .dem IO, and audited backends.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # 1. check SPEC.json [-o RECEIPT.json]
    p_check = sub.add_parser("check", help="Run a declarative JSON study and emit an ftqc_decoders-receipt-1")
    p_check.add_argument("spec", type=Path, help="Input JSON study specification")
    p_check.add_argument("-o", "--output", type=Path, default=None, help="Output receipt path")

    # 2. verify RECEIPT.json [--target SPEC.json]
    p_ver = sub.add_parser("verify", help="Deterministically replay and verify an ftqc_decoders-receipt-1")
    p_ver.add_argument("receipt", type=Path, help="Receipt JSON file")
    p_ver.add_argument("--target", type=Path, default=None, help="Optional original study JSON to verify against")

    # 2b. decode
    p_dec = sub.add_parser("decode", help="Decode a syndrome or injected fault list with any stdlib or external decoder")
    p_dec.add_argument("--code", required=True, help="Code spec (steane, golay, qr:47, bb72, bb144, toric:5, FILE.dem)")
    p_dec.add_argument("--sector", default="X", choices=["X", "Z"])
    p_dec.add_argument("--decoder", default="bp_osd", help="Decoder (bp_osd, bp_lsd, relay_bp, union_find, mwpm, ssf, sat_mwd, external:tesseract, ...)")
    p_dec.add_argument("--faults", default="", help="Comma-separated injected fault mechanism indices (e.g. 0,18,41)")
    p_dec.add_argument("--syndrome", type=int, default=None, help="Integer syndrome bitmask (if --faults omitted)")
    p_dec.add_argument("-o", "--output", type=Path, default=None)

    # 3. certify-radius
    p_rad = sub.add_parser("certify-radius", help="Certify operational decoding radius through weight R")
    p_rad.add_argument("--code", required=True, help="Code spec (e.g. steane, golay, qr:47, bb72, bb144, toric:5, FILE.dem)")
    p_rad.add_argument("--sector", default="X", choices=["X", "Z"])
    p_rad.add_argument("--decoder", default="bp_osd", help="Algorithm (bp_osd, bp_lsd, relay_bp, union_find, belief_find, mwpm, ssf, sat_mwd, external:tesseract, ...)")
    p_rad.add_argument("--osd-order", type=int, default=2)
    p_rad.add_argument("--radius", type=int, default=1)
    p_rad.add_argument("--budget", type=int, default=200_000)
    p_rad.add_argument("-o", "--output", type=Path, default=None)

    # 4. trapping-set
    p_ts = sub.add_parser("trapping-set", help="Search for split-stabilizer & midpoint trapping-set witnesses")
    p_ts.add_argument("--code", required=True)
    p_ts.add_argument("--sector", default="X", choices=["X", "Z"])
    p_ts.add_argument("--decoder", default="bp_only")
    p_ts.add_argument("--compare", default="relay_bp,bp_lsd,bp_osd", help="Comma-separated decoders to test on found witness")
    p_ts.add_argument("--max-exhaustive-weight", type=int, default=0)
    p_ts.add_argument("-o", "--output", type=Path, default=None)

    # 5. enclose
    p_enc = sub.add_parser("enclose", help="Compute zero-float outward-rounded logical error rate enclosures")
    p_enc.add_argument("--code", required=True)
    p_enc.add_argument("--sector", default="X", choices=["X", "Z"])
    p_enc.add_argument("--mode", default="peierls_cluster", choices=["exact_shells", "peierls_cluster", "stratified_shells"])
    p_enc.add_argument("--decoder", default="bp_osd")
    p_enc.add_argument("--p", default="1/1000", help="Exact rational or decimal physical error probability string")
    p_enc.add_argument("--radius", type=int, default=1, help="Verified operational / component radius")
    p_enc.add_argument("--max-weight", type=int, default=2, help="Max shell weight for exact_shells / stratified_shells")
    p_enc.add_argument("-o", "--output", type=Path, default=None)

    # 6. convert-dem
    p_dem = sub.add_parser("convert-dem", help="Export any code / DEM spec to canonical Stim .dem text")
    p_dem.add_argument("--code", required=True)
    p_dem.add_argument("--sector", default="X", choices=["X", "Z"])
    p_dem.add_argument("--prior", default="1/1000")
    p_dem.add_argument("-o", "--output", type=Path, default=None)

    # 7. backends
    sub.add_parser("backends", help="List available stdlib and optional external QEC decoder backends")

    # 8. selftest
    sub.add_parser("selftest", help="Run the self-contained <2s cold-clone verification & mutant battery")

    args = parser.parse_args(argv)

    try:
        if args.command == "selftest":
            from .selftest import run_all

            return run_all()

        if args.command == "backends":
            info = {
                "stdlib_zero_float_decoders": [
                    "bp_only",
                    "bp_osd",
                    "bp_lsd",
                    "relay_bp",
                    "union_find",
                    "belief_find",
                    "mwpm",
                    "ssf",
                    "sat_mwd",
                ],
                "external_optional_backends": probe_external_backends(),
            }
            _emit_json(info, None)
            return 0

        if args.command == "convert-dem":
            prob = load_problem_from_spec(_parse_code_shorthand(args.code, prior=args.prior, sector=args.sector))
            dem_txt = prob.to_stim_dem_text()
            if args.output is not None:
                args.output.write_text(dem_txt)
            else:
                sys.stdout.write(dem_txt)
            return 0

        if args.command == "check":
            raw_spec = json.loads(args.spec.read_text(), parse_float=lambda s: (_ for _ in ()).throw(RefusalError("FLOAT-INPUT-REFUSED", f"Float literal {s} forbidden in JSON spec; quote as string")))
            rec = build_receipt(raw_spec, base_dir=args.spec.parent)
            _emit_json(rec, args.output)
            st = rec["result"].get("status", "proved")
            return 1 if st in ("counterexample", "failed", "non_converged") else 0

        if args.command == "verify":
            raw_rec = json.loads(args.receipt.read_text(), parse_float=lambda s: (_ for _ in ()).throw(RefusalError("FLOAT-INPUT-REFUSED", f"Float literal {s} forbidden in receipt")))
            target_obj = None
            base_d = args.receipt.parent
            if args.target is not None:
                target_obj = json.loads(args.target.read_text(), parse_float=lambda s: (_ for _ in ()).throw(RefusalError("FLOAT-INPUT-REFUSED", f"Float literal {s} forbidden in target")))
                base_d = args.target.parent
            rep = verify_receipt(raw_rec, target_study=target_obj, base_dir=base_d)
            _emit_json(rep, None)
            return 0

        if args.command == "decode":
            study_d: dict[str, Any] = {
                "study": "decode",
                "problem": _parse_code_shorthand(args.code, sector=args.sector),
                "decoder": args.decoder,
            }
            if args.faults.strip():
                study_d["fault_mechanisms"] = [int(x.strip()) for x in args.faults.split(",") if x.strip()]
            else:
                study_d["syndrome"] = int(args.syndrome or 0)
            rec = build_receipt(study_d)
            _emit_json(rec, args.output)
            return 0 if rec["result"]["converged"] else 1

        if args.command == "certify-radius":
            dec_spec: str | dict[str, Any] = (
                {"algorithm": args.decoder, "osd_order": args.osd_order}
                if args.decoder in ("bp_osd", "bp_lsd")
                else args.decoder
            )
            study = {
                "study": "certify_radius",
                "problem": _parse_code_shorthand(args.code, sector=args.sector),
                "decoder": dec_spec,
                "target_radius": args.radius,
                "budget_cap": args.budget,
            }
            rec = build_receipt(study)
            _emit_json(rec, args.output)
            return 1 if rec["result"]["status"] == "counterexample" else 0

        if args.command == "trapping-set":
            comps = [c.strip() for c in args.compare.split(",") if c.strip()]
            study = {
                "study": "trapping_set",
                "problem": _parse_code_shorthand(args.code, sector=args.sector),
                "decoder": args.decoder,
                "max_exhaustive_weight": args.max_exhaustive_weight,
                "compare_decoders": comps,
            }
            rec = build_receipt(study)
            _emit_json(rec, args.output)
            return 0

        if args.command == "enclose":
            mode_map = {
                "exact_shells": "enclose_exact_shells",
                "peierls_cluster": "enclose_peierls_cluster",
                "stratified_shells": "enclose_stratified_shells",
            }
            study = {
                "study": mode_map[args.mode],
                "problem": _parse_code_shorthand(args.code, sector=args.sector),
                "decoder": args.decoder,
                "physical_p": args.p,
            }
            if args.mode == "peierls_cluster":
                study["verified_component_radius"] = args.radius
            elif args.mode == "exact_shells":
                study["verified_radius"] = args.radius
                study["max_exact_weight"] = args.max_weight
            else:
                study["verified_radius"] = args.radius
                study["max_shell_weight"] = args.max_weight
            rec = build_receipt(study)
            _emit_json(rec, args.output)
            return 0

    except RefusalError as exc:
        sys.stderr.write(f"REFUSED: {exc}\n")
        return 2
    except (KeyError, ValueError, json.JSONDecodeError) as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        return 2

    return 0
