"""Command-Line Interface for `ftqc_faults.zx`.

Exit codes (BootLoops convention):
  0 : PASS / certified / clean
  1 : Mathematical finding (`FAULT-EQUIV-FAIL`, `MAP-MISMATCH`, `PROOF-INVALID`, `CONTRACT-UNSAFE-SINGLE-FAULT`)
  2 : Input / precondition refusal (`RefusalError`) or CLI usage error

Zero float literals or float builtin references appear anywhere in this module.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

if __package__ in (None, ""):
    _PKG_DIR = Path(__file__).resolve().parent
    sys.path = [p for p in sys.path if Path(p or ".").resolve() != _PKG_DIR]
    sys.path.insert(0, str(_PKG_DIR.parents[1]))
    from ftqc_faults.zx.adapters import circuit_to_qasm2, circuit_to_stim_text
    from ftqc_faults.zx.contracts import build_gadget_contract
    from ftqc_faults.zx.diagram import CliffordCircuit, RefusalError, ZXDiagram, circuit_to_zx_fe, load_spec_file
    from ftqc_faults.zx.equivalence import check_fault_equiv, circuit_distance, resolve_code
    from ftqc_faults.zx.pauliweb import analyze_pauli_webs, compile_outcome_map
    from ftqc_faults.zx.rewrites import DerivationCertificate, verify_derivation
    from ftqc_faults.zx.synthesize import (
        make_cat_spec,
        make_full_syndrome_spec,
        make_steane_7_1_3,
        make_syndrome_meas_spec,
        synth_cat_state,
        synth_sequential_cat4_steane_bad,
        synth_shor_syndrome,
        synth_steane_goto_prep,
        synth_steane_syndrome,
        synth_steane_zero_like,
        synth_unflagged_cat_state,
    )
else:
    from .adapters import circuit_to_qasm2, circuit_to_stim_text
    from .contracts import build_gadget_contract
    from .diagram import CliffordCircuit, RefusalError, ZXDiagram, circuit_to_zx_fe, load_spec_file
    from .equivalence import check_fault_equiv, circuit_distance, resolve_code
    from .pauliweb import analyze_pauli_webs, compile_outcome_map
    from .rewrites import DerivationCertificate, verify_derivation
    from .synthesize import (
        make_cat_spec,
        make_full_syndrome_spec,
        make_steane_7_1_3,
        make_syndrome_meas_spec,
        synth_cat_state,
        synth_sequential_cat4_steane_bad,
        synth_shor_syndrome,
        synth_steane_goto_prep,
        synth_steane_syndrome,
        synth_steane_zero_like,
        synth_unflagged_cat_state,
    )

PRESET_NAMES = (
    "cat4",
    "cat4_unflagged",
    "cat",
    "cat_spec",
    "shor_standard",
    "shor_opt",
    "shor_teleport",
    "shor_spec",
    "steane_goto",
    "steane_zero_like",
    "steane_8anc",
    "steane_5anc",
    "steane_spec",
    "steane_bad_unflagged_cat4",
)


def _resolve_preset_or_file(
    target: str,
    n_qubits: int = 4,
    w_param: int | None = None,
    basis: str = "Z",
    pauli: str = "ZZZZ",
    sector: str = "X",
) -> tuple[CliffordCircuit | ZXDiagram, DerivationCertificate | None]:
    p_path = Path(target)
    if p_path.is_file():
        return load_spec_file(str(p_path)), None

    if target == "cat4":
        c, cert = synth_cat_state(4, w=None, basis=basis)
        return c, cert
    if target == "cat4_unflagged":
        return synth_unflagged_cat_state(4, basis=basis), None
    if target == "cat":
        c, cert = synth_cat_state(n_qubits, w=w_param, basis=basis)
        return c, cert
    if target == "cat_spec":
        return make_cat_spec(n_qubits, basis=basis), None
    if target in ("shor_standard", "shor_opt", "shor_teleport"):
        vmap = {"shor_standard": "standard", "shor_opt": "optimised", "shor_teleport": "teleport"}
        c, cert = synth_shor_syndrome(pauli, variant=vmap[target])
        return c, cert
    if target == "shor_spec":
        return make_syndrome_meas_spec(pauli, mode="ft"), None
    if target == "steane_goto":
        return synth_steane_goto_prep(), None
    if target == "steane_zero_like":
        return synth_steane_zero_like(), None
    if target == "steane_8anc":
        c, cert = synth_steane_syndrome(make_steane_7_1_3(), mode="standard_8anc", sector=sector)
        return c, cert
    if target == "steane_5anc":
        c, cert = synth_steane_syndrome(make_steane_7_1_3(), mode="optimised_5anc", sector=sector)
        return c, cert
    if target == "steane_spec":
        return make_full_syndrome_spec(make_steane_7_1_3(), sector=sector), None
    if target == "steane_bad_unflagged_cat4":
        return synth_sequential_cat4_steane_bad(), None

    raise RefusalError(
        "UNKNOWN-TARGET",
        f"{target!r} is neither an existing JSON file nor a known preset {PRESET_NAMES}",
    )


def _emit(payload: dict[str, Any], as_json: bool, summary: str) -> None:
    if as_json:
        print(json.dumps(payload, indent=2, sort_keys=False))
    else:
        print(summary)


def cmd_audit(args: argparse.Namespace) -> int:
    obj, _ = _resolve_preset_or_file(args.target, n_qubits=args.n, basis=args.basis, pauli=args.pauli, sector=args.sector)
    code = resolve_code(args.code) if args.code else None
    zx = circuit_to_zx_fe(obj, noise_model=args.noise_model) if isinstance(obj, CliffordCircuit) else obj
    ana = analyze_pauli_webs(zx)
    dist = circuit_distance(obj, w_max=args.w_max, noise_model=args.noise_model, codespace=code)
    omap = compile_outcome_map(zx) if ana.is_nonzero else None

    payload = {
        "verdict": "AUDIT-OK" if ana.is_nonzero else "AUDIT-ZERO-MAP",
        "name": obj.name,
        "noise_model": args.noise_model,
        "is_nonzero": ana.is_nonzero,
        "n_zx_nodes": len(zx.nodes),
        "n_zx_edges": len(zx.edges),
        "n_detecting_regions": len(ana.detecting_regions),
        "n_boundary_stabilizers": len(ana.signed_boundary_stabs),
        "n_distinct_branches": omap.n_distinct_branches if omap else 0,
        "distance": dist.to_dict(),
    }
    summary = (
        f"[{payload['verdict']}] {obj.name} ({args.noise_model}): "
        f"nodes={len(zx.nodes)} edges={len(zx.edges)} "
        f"detectors={len(ana.detecting_regions)} stabs={len(ana.signed_boundary_stabs)} "
        f"dist={'==' + str(dist.exact_distance) if dist.exact_distance is not None else '>=' + str(dist.lower_bound)}"
    )
    _emit(payload, args.json, summary)
    return 0 if ana.is_nonzero else 1


def cmd_equiv(args: argparse.Namespace) -> int:
    lhs, _ = _resolve_preset_or_file(args.lhs, n_qubits=args.n, basis=args.basis, pauli=args.pauli, sector=args.sector)
    rhs, _ = _resolve_preset_or_file(args.rhs, n_qubits=args.n, basis=args.basis, pauli=args.pauli, sector=args.sector)
    code = resolve_code(args.code) if args.code else None
    res = check_fault_equiv(
        lhs,
        rhs,
        w=args.w,
        noise_model_lhs=args.noise_model,
        noise_model_rhs=args.noise_model_rhs,
        codespace=code,
    )
    payload = res.to_dict()
    summary = (
        f"[{res.verdict}] w={res.w} lhs={lhs.name} ({res.noise_model_lhs}) "
        f"vs rhs={rhs.name} ({res.noise_model_rhs}) "
        f"codespace={res.codespace_name} spackle={res.spackle_crosschecked}"
    )
    if res.counterexample:
        summary += f" cex={res.counterexample}"
    _emit(payload, args.json, summary)
    return 0 if res.equivalent else 1


def cmd_verify_proof(args: argparse.Namespace) -> int:
    p_path = Path(args.cert)
    if p_path.is_file():
        data = json.loads(p_path.read_text(encoding="utf-8"))
        cert = DerivationCertificate.from_dict(data)
    else:
        _, cert_opt = _resolve_preset_or_file(
            args.cert, n_qubits=args.n, w_param=args.w, basis=args.basis, pauli=args.pauli, sector=args.sector
        )
        if cert_opt is None:
            raise RefusalError("NO-PROOF-CERTIFICATE", f"preset {args.cert!r} has no constructive proof certificate")
        cert = cert_opt

    res = verify_derivation(cert)
    summary = (
        f"[{res['verdict']}] {cert.name}: certified_w={res.get('certified_w', 'inf')} "
        f"steps={res.get('n_steps', 0)} subproofs={res.get('n_subproofs', 0)}"
    )
    _emit(res, args.json, summary)
    return 0 if res["verdict"] == "PROOF-VALID" else 1


def cmd_synth(args: argparse.Namespace) -> int:
    obj, cert = _resolve_preset_or_file(
        args.preset, n_qubits=args.n, w_param=args.w, basis=args.basis, pauli=args.pauli, sector=args.sector
    )
    if not isinstance(obj, CliffordCircuit):
        raise RefusalError("SYNTH-EXPECTED-CIRCUIT", "synth preset did not produce a CliffordCircuit")

    if args.out:
        Path(args.out).write_text(json.dumps(obj.to_dict(), indent=2) + "\n", encoding="utf-8")
    if args.proof_out:
        if cert is None:
            raise RefusalError("NO-PROOF-CERTIFICATE", f"preset {args.preset!r} has no derivation certificate")
        Path(args.proof_out).write_text(json.dumps(cert.to_dict(), indent=2) + "\n", encoding="utf-8")
    if args.stim_out:
        stim_txt = circuit_to_stim_text(obj, p_phys=args.p, noise_model=args.noise_model)
        Path(args.stim_out).write_text(stim_txt, encoding="utf-8")
    if args.qasm_out:
        qasm_txt = circuit_to_qasm2(obj)
        Path(args.qasm_out).write_text(qasm_txt, encoding="utf-8")

    payload = {
        "verdict": "SYNTH-OK",
        "preset": args.preset,
        "circuit": obj.to_dict(),
        "has_proof_certificate": cert is not None,
    }
    summary = (
        f"[SYNTH-OK] {obj.name}: n_qubits={obj.n_qubits} ops={len(obj.ops)} "
        f"cert={'yes' if cert is not None else 'no'}"
    )
    _emit(payload, args.json, summary)
    return 0


def cmd_contract(args: argparse.Namespace) -> int:
    obj, _ = _resolve_preset_or_file(
        args.target, n_qubits=args.n, w_param=args.w, basis=args.basis, pauli=args.pauli, sector=args.sector
    )
    if not isinstance(obj, CliffordCircuit):
        raise RefusalError("CONTRACT-REQUIRES-CIRCUIT", "stochastic contracts require a physical CliffordCircuit")
    code = resolve_code(args.code) if args.code else None
    postselect: bool | None = True if args.postselect else (False if args.inline else None)
    contract = build_gadget_contract(
        obj,
        p_phys=args.p,
        idle_p_ratio=args.idle_ratio,
        noise_model=args.noise_model,
        codespace=code,
        postselect_detectors=postselect,
        extracted_sector=args.sector,
        certified_w=args.certified_w,
    )
    if args.out:
        Path(args.out).write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    sc = contract["schedule"]
    au = contract["single_fault_audit"]
    en = contract["enclosures"]
    summary = (
        f"[{contract['verdict']}] {obj.name} ({contract['noise_model']}, mode={contract['operating_mode']}): "
        f"ticks={sc['ticks']} locs={sc['total_locations']} W={sc['wire_tick_volume']} "
        f"bad_labels={au['bad_unheralded_labels']} esc_hi={en['conditional_escape_hi']}"
    )
    _emit(contract, args.json, summary)
    return 0 if contract["verdict"] == "CONTRACT-CERTIFIED" else 1


def cmd_selftest(args: argparse.Namespace) -> int:
    if __package__ in (None, ""):
        from ftqc_faults.zx.battery import run_battery
    else:
        from .battery import run_battery

    report = run_battery(include_mutations=args.mutation_controls)
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        for leg in report["legs"]:
            print(f"  [{'PASS' if leg['passed'] else 'FAIL'}] {leg['id']}: {leg['summary']}")
        print(f"[{'PASS' if report['passed'] else 'FAIL'}] ftqc_faults.zx battery ({report['passed_count']}/{report['total_count']} legs)")
    return 0 if report["passed"] else 1


def build_parser(prog="ftqc_faults zx") -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog=prog,
        description="ZX-calculus w-fault equivalence, Pauli-web verifier, FTQC synthesizer, and zero-float contract compiler.",
    )
    sub = ap.add_subparsers(dest="subcmd", required=True)

    def add_common_preset_args(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--n", type=int, default=4, help="Qubit count for cat presets (default: 4)")
        sp.add_argument("--basis", choices=("Z", "X"), default="Z", help="Cat basis (default: Z)")
        sp.add_argument("--pauli", default="ZZZZ", help="Pauli word for Shor presets (default: ZZZZ)")
        sp.add_argument("--sector", choices=("X", "Z"), default="X", help="CSS sector for Steane presets (default: X)")

    p_aud = sub.add_parser("audit", help="Pauli-web audit, outcome push-out map, and circuit distance")
    p_aud.add_argument("target", help="JSON spec file or preset name")
    p_aud.add_argument("--noise-model", default="circuit_depolarizing")
    p_aud.add_argument("--w-max", type=int, default=3)
    p_aud.add_argument("--code", default=None, help="Optional CSS codespace preset (e.g. steane_7_1_3)")
    p_aud.add_argument("--json", action="store_true")
    add_common_preset_args(p_aud)

    p_eq = sub.add_parser("equiv", help="Bidirectional semantic w-fault equivalence check (+ spackle cross-check)")
    p_eq.add_argument("lhs", help="LHS JSON file or preset")
    p_eq.add_argument("rhs", help="RHS JSON file or preset")
    p_eq.add_argument("--w", type=int, required=True, help="Adversarial fault-equivalence bound w >= 1")
    p_eq.add_argument("--noise-model", default="circuit_depolarizing")
    p_eq.add_argument("--noise-model-rhs", default=None)
    p_eq.add_argument("--code", default=None, help="Optional CSS codespace preset (e.g. steane_7_1_3)")
    p_eq.add_argument("--json", action="store_true")
    add_common_preset_args(p_eq)

    p_vp = sub.add_parser("verify-proof", help="Polynomial-time ZX derivation certificate verifier")
    p_vp.add_argument("cert", help="DerivationCertificate JSON file or preset name")
    p_vp.add_argument("--w", type=int, default=None, help="Optional w parameter for cat preset")
    p_vp.add_argument("--json", action="store_true")
    add_common_preset_args(p_vp)

    p_syn = sub.add_parser("synth", help="Synthesize an FTQC gadget circuit + ZX proof certificate")
    p_syn.add_argument("preset", choices=PRESET_NAMES)
    p_syn.add_argument("--w", type=int, default=None, help="Optional w truncation for cat preset")
    p_syn.add_argument("--out", default=None, help="Output circuit JSON path")
    p_syn.add_argument("--proof-out", default=None, help="Output DerivationCertificate JSON path")
    p_syn.add_argument("--stim-out", default=None, help="Output .stim path (with Pauli-web DETECTORs)")
    p_syn.add_argument("--qasm-out", default=None, help="Output OpenQASM 2.0 path")
    p_syn.add_argument("--p", default="0", help="Physical fault rate for .stim export (default: 0)")
    p_syn.add_argument("--noise-model", default="circuit_depolarizing")
    p_syn.add_argument("--json", action="store_true")
    add_common_preset_args(p_syn)

    p_con = sub.add_parser("contract", help="Emit zero-float outward-rounded ftqc_faults-zx-contract-1 JSON")
    p_con.add_argument("target", help="JSON circuit file or preset name")
    p_con.add_argument("--p", default="1e-4", help="Exact physical fault rate string (e.g. 1e-4 or 1/10000)")
    p_con.add_argument("--idle-ratio", default="1", help="Exact idle-to-gate fault ratio (e.g. 1 or 1/10)")
    p_con.add_argument("--noise-model", default="circuit_depolarizing")
    p_con.add_argument("--code", default=None, help="Optional CSS codespace preset (e.g. steane_7_1_3)")
    p_con.add_argument("--postselect", action="store_true", help="Force postselected factory operating mode")
    p_con.add_argument("--inline", action="store_true", help="Force deterministic inline operating mode")
    p_con.add_argument("--certified-w", type=int, default=3)
    p_con.add_argument("--w", type=int, default=None)
    p_con.add_argument("--out", default=None, help="Output contract JSON path")
    p_con.add_argument("--json", action="store_true")
    add_common_preset_args(p_con)

    p_st = sub.add_parser("selftest", help="Run self-test and negative-control battery")
    p_st.add_argument("--mutation-controls", action="store_true", help="Also run 4-mutant scratch sabotage checks")
    p_st.add_argument("--json", action="store_true")

    return ap


def main(argv: list[str] | None = None, *, prog="ftqc_faults zx") -> int:
    parser = build_parser(prog)
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        code = exc.code if isinstance(exc.code, int) else 2
        return 0 if code == 0 else 2

    try:
        if args.subcmd == "audit":
            return cmd_audit(args)
        if args.subcmd == "equiv":
            return cmd_equiv(args)
        if args.subcmd == "verify-proof":
            return cmd_verify_proof(args)
        if args.subcmd == "synth":
            return cmd_synth(args)
        if args.subcmd == "contract":
            return cmd_contract(args)
        if args.subcmd == "selftest":
            return cmd_selftest(args)
        return 2
    except RefusalError as exc:
        if getattr(args, "json", False):
            print(json.dumps(exc.to_dict(), indent=2))
        else:
            print(f"REFUSAL {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
