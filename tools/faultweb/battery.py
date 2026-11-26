"""In-tree Synthetic Verification & 4-Mutant Sabotage Battery for `faultweb`.

Runs in < 5 s with zero external dependencies, verifying:
  - L1: AST zero-float literal & builtin purity scan across the entire package,
  - L2: Cat-state family (Props. 7.1-7.3), dual-lineage (`pauliweb` == `spackle`), and 64-qubit proof verifier,
  - L3: Shor syndrome variants (`standard`, `optimised`, `teleport`; Props. 7.5-7.7) + Appendix A outcome push-out,
  - L4: Steane 8-ancilla & 5-ancilla (Props. 7.8-7.12, X & Z sectors), closed circuit distance, + Appendix C.2 flag decoder,
  - L5: Zero-float outward-rounded `faultweb-contract-v1` enclosures, `examples/*.json`, QASM/Stim bridges, and CLI exit codes,
  - N1-N5: Planted negative controls and fail-closed refusals,
  - M1-M4 (when `include_mutations=True`): Scratch-copy single-line sabotages proving every verification gate is live!

Zero float literals or float builtin references appear anywhere in this module.
"""

from __future__ import annotations

import ast
import contextlib
from fractions import Fraction
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any, Callable

from .adapters import (
    circuit_from_qasm2,
    circuit_to_qasm2,
    circuit_to_stim_text,
    from_pyzx_graph,
    pyzx_available,
    sample_stim_detectors,
    stim_available,
    to_pyzx_graph,
)
from .cli import main as cli_main
from .contracts import build_gadget_contract, fraction_to_outward_interval_strs, parse_exact_rate
from .decoder_synth import synthesize_flag_decoder
from .diagram import (
    CircuitOp,
    CliffordCircuit,
    RefusalError,
    ZXDiagram,
    circuit_to_zx_fe,
    extract_circuit_fe,
    load_spec_file,
)
from .equivalence import (
    check_fault_equiv,
    circuit_distance,
    make_golay_23_1_7,
    make_shor_9_1_3,
    make_steane_7_1_3,
)
from .pauliweb import compile_outcome_map
from .rewrites import RewriteStep, apply_zx_rewrite, verify_derivation
from .synthesize import (
    make_cat_spec,
    make_full_syndrome_spec,
    make_syndrome_meas_spec,
    synth_cat_state,
    synth_sequential_cat4_steane_bad,
    synth_shor_syndrome,
    synth_steane_goto_prep,
    synth_steane_syndrome,
    synth_steane_zero_like,
    synth_unflagged_cat_state,
)

PKG_DIR = Path(__file__).resolve().parent


def _assert_json_zero_floats(payload: Any) -> None:
    raw = json.dumps(payload)

    def _forbid_float(token: str) -> None:
        raise AssertionError(f"JSON output contains forbidden floating-point number {token!r}")

    json.loads(raw, parse_float=_forbid_float)


def leg_l1_purity_scan() -> str:
    count = 0
    for py_path in sorted(PKG_DIR.rglob("*.py")):
        count += 1
        tree = ast.parse(py_path.read_text(encoding="utf-8"), filename=str(py_path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and type(node.value).__name__ == "float":
                raise AssertionError(f"Float literal {node.value!r} in {py_path.name}:{node.lineno}")
            if isinstance(node, ast.Name) and node.id == "float":
                raise AssertionError(f"Builtin `float` reference in {py_path.name}:{node.lineno}")
    return f"scanned {count} Python modules; 0 float literals or float builtin references"


def leg_l2_cat_and_64q_proof() -> str:
    spec4 = make_cat_spec(4, basis="Z")
    cat4, cert4 = synth_cat_state(4, w=None, basis="Z")
    r4 = check_fault_equiv(cat4, spec4, w=4, noise_model_lhs="circuit_depolarizing")
    assert r4.equivalent and r4.spackle_crosschecked
    assert verify_derivation(cert4)["verdict"] == "PROOF-VALID"

    spec5 = make_cat_spec(5, basis="Z")
    cat5, cert5 = synth_cat_state(5, w=None, basis="Z")
    assert check_fault_equiv(cat5, spec5, w=3, noise_model_lhs="circuit_depolarizing").equivalent
    assert verify_derivation(cert5)["verdict"] == "PROOF-VALID"

    spec8 = make_cat_spec(8, basis="Z")
    cat8_full, cert8_full = synth_cat_state(8, w=None, basis="Z")
    assert check_fault_equiv(cat8_full, spec8, w=4, noise_model_lhs="circuit_depolarizing").equivalent
    assert verify_derivation(cert8_full)["certified_w"] is None

    for wt in (2, 3):
        cat8_w, cert8_w = synth_cat_state(8, w=wt, basis="Z")
        assert check_fault_equiv(cat8_w, spec8, w=wt, noise_model_lhs="circuit_depolarizing").equivalent
        assert not check_fault_equiv(cat8_w, spec8, w=wt + 1, noise_model_lhs="circuit_depolarizing").equivalent
        pv = verify_derivation(cert8_w)
        assert pv["verdict"] == "PROOF-VALID" and pv["certified_w"] == wt

    t0 = time.perf_counter_ns()
    _, cert64 = synth_cat_state(64, w=3, basis="Z")
    pv64 = verify_derivation(cert64)
    dt_ms = (time.perf_counter_ns() - t0) // 1_000_000
    assert pv64["verdict"] == "PROOF-VALID" and pv64["certified_w"] == 3
    return f"Props 7.1-7.3 verified with spackle cross-check; 64-cat w=3 proof verified in {dt_ms} ms"


def leg_l3_shor_and_outcome_map() -> str:
    for pstr in ("ZZZZ", "XXXX", "XZYZ"):
        spec = make_syndrome_meas_spec(pstr, mode="ft")
        for var in ("standard", "optimised", "teleport"):
            circ, cert = synth_shor_syndrome(pstr, variant=var)
            res = check_fault_equiv(circ, spec, w=3, noise_model_lhs="circuit_depolarizing")
            assert res.equivalent and res.spackle_crosschecked, (pstr, var, res)
            assert verify_derivation(cert)["verdict"] == "PROOF-VALID"

    sh_opt, _ = synth_shor_syndrome("ZZZZ", variant="optimised")
    omap = compile_outcome_map(circuit_to_zx_fe(sh_opt, "fault_free"))
    assert omap.n_distinct_branches == 1 and len(omap.parity_checks) == 1
    return "Props 7.5-7.7 (standard/optimised/teleport) verified on ZZZZ, XXXX, XZYZ + Appendix A push-out"


def leg_l4_steane_and_flag_decoder() -> str:
    steane = make_steane_7_1_3()
    assert len(make_shor_9_1_3().hz_rows) == 6 and len(make_golay_23_1_7().hx_rows) == 11

    goto_p = synth_steane_goto_prep()
    ideal_z = CliffordCircuit(
        7,
        (),
        tuple(range(7)),
        [
            CircuitOp(
                "ENC_PREP",
                tuple(range(7)),
                mode="ff",
                params={
                    "input_qubits": [],
                    "stabilizers": steane.stabilizers_pauli_strings(),
                    "logical_z": steane.logical_z_strings(),
                },
            )
        ],
        name="ideal_steane_zero",
    )
    assert check_fault_equiv(goto_p, ideal_z, w=3, noise_model_lhs="circuit_css").equivalent

    # Closed Goto + ideal verifier has exact distance 3 (Theorem 3.11)
    ver_spec = make_full_syndrome_spec(steane, sector="BOTH")
    seg_cnt = {q: 0 for q in range(8)}
    v_ops = list(goto_p.ops)
    v_ff: set[tuple[int, int]] = set()
    for op in goto_p.ops:
        if op.kind not in ("PREP_0", "PREP_1", "PREP_PLUS", "PREP_MINUS", "MEAS_Z", "MEAS_X", "MEAS_Y"):
            for q in op.qubits:
                seg_cnt[q] += 1
    for m_i, op in enumerate(ver_spec.ops):
        if m_i > 0:
            for q in range(7):
                v_ff.add((q, seg_cnt[q]))
        v_ops.append(op)
        for q in range(7):
            seg_cnt[q] += 1
    for q in range(7):
        v_ff.add((q, seg_cnt[q]))
    goto_closed = CliffordCircuit(8, (), tuple(range(7)), v_ops, fault_free_wires=frozenset(v_ff), name="goto_closed")
    assert circuit_distance(goto_closed, w_max=3, noise_model="circuit_css").exact_distance == 3

    zl = synth_steane_zero_like()
    zl_spec = CliffordCircuit(
        zl.n_qubits,
        zl.input_qubits,
        zl.output_qubits,
        [CircuitOp(op.kind, op.qubits, mode="ff", label=op.label, params=dict(op.params)) for op in zl.ops],
        name="ideal_zero_like",
    )
    assert check_fault_equiv(zl, zl_spec, w=3, noise_model_lhs="circuit_css").equivalent

    for sec in ("X", "Z"):
        spec_s = make_full_syndrome_spec(steane, sector=sec)
        for mode in ("standard_8anc", "optimised_5anc"):
            c_st, cert_st = synth_steane_syndrome(steane, mode=mode, sector=sec)
            res = check_fault_equiv(c_st, spec_s, w=3, noise_model_lhs="circuit_css", codespace=steane)
            assert res.equivalent and res.spackle_crosschecked
            assert verify_derivation(cert_st)["verdict"] == "PROOF-VALID"

    c5x, _ = synth_steane_syndrome(steane, mode="optimised_5anc", sector="X")
    dec = synthesize_flag_decoder(c5x, steane, extracted_sector="X")
    assert dec.h_prime_rows == ([1, 1, 1, 0], [1, 1, 0, 1], [0, 1, 1, 0])
    assert dec.unflagged_max_coset_weight <= 1 and dec.unambiguous_weight1
    return "Props 7.8-7.12 (8-anc & 5-anc X/Z) + closed dist==3 + Appendix C.2 flag decoder certified"


def leg_l5_contracts_examples_and_cli() -> str:
    assert parse_exact_rate("1e-4") == Fraction(1, 10000)
    for tf in (Fraction(1, 3), Fraction(1, 7), Fraction(9999, 10000) ** 60):
        lo_s, hi_s = fraction_to_outward_interval_strs(tf)
        assert Fraction(lo_s) < tf < Fraction(hi_s)

    steane = make_steane_7_1_3()
    c8, _ = synth_steane_syndrome(steane, mode="standard_8anc", sector="X")
    c5, _ = synth_steane_syndrome(steane, mode="optimised_5anc", sector="X")
    con8 = build_gadget_contract(c8, p_phys="1e-4", idle_p_ratio="1/10", noise_model="circuit_css", codespace=steane)
    con5 = build_gadget_contract(c5, p_phys="1e-4", idle_p_ratio="1/10", noise_model="circuit_css", codespace=steane)
    _assert_json_zero_floats(con8)
    _assert_json_zero_floats(con5)
    assert con5["verdict"] == "CONTRACT-CERTIFIED" and con8["verdict"] == "CONTRACT-CERTIFIED"
    assert Fraction(con5["enclosures"]["conditional_escape_hi"]) < Fraction(
        con8["enclosures"]["conditional_escape_hi"]
    )

    ex_dir = PKG_DIR / "examples"
    for ex_name in ("cat4_flagged.json", "shor_opt_zzzz.json", "steane_5anc_x.json"):
        loaded = load_spec_file(str(ex_dir / ex_name))
        assert isinstance(loaded, CliffordCircuit)

    cat4, _ = synth_cat_state(4)
    qasm = circuit_to_qasm2(cat4)
    cat4_rt = circuit_from_qasm2(qasm, input_qubits=(), output_qubits=(0, 1, 2, 3))
    assert check_fault_equiv(cat4_rt, make_cat_spec(4), w=4).equivalent

    sh_opt, _ = synth_shor_syndrome("ZZZZ", variant="optimised")
    stim_txt = circuit_to_stim_text(sh_opt, p_phys="1e-4")
    assert "DETECTOR rec[" in stim_txt

    opt_notes: list[str] = []
    if pyzx_available():
        zx4 = circuit_to_zx_fe(cat4, "circuit_css")
        assert len(from_pyzx_graph(to_pyzx_graph(zx4)).outputs) == 4
        opt_notes.append("pyzx=live")
    else:
        opt_notes.append("pyzx=skip")
    if stim_available():
        assert sample_stim_detectors(sh_opt, shots=64, p_phys="0")["triggered_shots"] == 0
        opt_notes.append("stim=live")
    else:
        opt_notes.append("stim=skip")

    with tempfile.TemporaryDirectory() as tmp:
        out_con = Path(tmp) / "c.json"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            assert cli_main([
                "contract", "steane_5anc", "--p", "1e-4", "--noise-model", "circuit_css",
                "--code", "steane_7_1_3", "--out", str(out_con), "--json",
            ]) == 0
            assert out_con.is_file()
            assert cli_main(["equiv", "cat4_unflagged", "cat_spec", "--w", "2", "--json"]) == 1
            assert cli_main(["equiv", "cat4", "cat_spec", "--w", "0", "--json"]) == 2

    return f"outward intervals, zero-float JSON contracts, examples, QASM/Stim, CLI (0/1/2) ok ({', '.join(opt_notes)})"


def leg_negatives_and_refusals() -> str:
    steane = make_steane_7_1_3()
    spec_sx = make_full_syndrome_spec(steane, sector="X")

    # N1: Unflagged 4-cat fails w=2
    assert not check_fault_equiv(synth_unflagged_cat_state(4), make_cat_spec(4), w=2).equivalent

    # N2: Bare CNOT extraction under circuit_depolarizing -> UNGUARDED-CNOT-EXTRACTION
    bare_cx = CliffordCircuit(2, (0, 1), (0, 1), [CircuitOp("CX", (0, 1))], name="bare_cx")
    bare_zx = circuit_to_zx_fe(bare_cx, noise_model="circuit_css")
    try:
        extract_circuit_fe(bare_zx, bare_cx, target_noise_model="circuit_depolarizing")
        raise AssertionError("Expected UNGUARDED-CNOT-EXTRACTION")
    except RefusalError as exc:
        assert exc.code == "UNGUARDED-CNOT-EXTRACTION"

    # N3: Steane 5-anc without codespace fails w=2 (Appendix C.1); unflagged-cat4 Steane fails w=2 with codespace
    c5, _ = synth_steane_syndrome(steane, mode="optimised_5anc", sector="X")
    assert not check_fault_equiv(c5, spec_sx, w=2, noise_model_lhs="circuit_css", codespace=None).equivalent
    assert not check_fault_equiv(
        synth_sequential_cat4_steane_bad(), spec_sx, w=2, noise_model_lhs="circuit_depolarizing", codespace=steane
    ).equivalent

    # N4: fuse_1_fe side conditions & ideal_fuse_fe purple-edge side condition
    zx_sc = ZXDiagram("scalar")
    e_sc = zx_sc.add_edge(zx_sc.add_node("Z", 0), zx_sc.add_node("Z", 0), fault_free=False)
    try:
        apply_zx_rewrite(zx_sc, RewriteStep("fuse_1_fe", {"edge_id": e_sc}))
        raise AssertionError("Expected FUSE1-ISOLATED-SCALAR-REFUSED")
    except RefusalError as exc:
        assert exc.code == "FUSE1-ISOLATED-SCALAR-REFUSED"

    zx_hub = ZXDiagram("hubs")
    hu, hv = zx_hub.add_node("Z", 0), zx_hub.add_node("Z", 0)
    for _ in range(2):
        zx_hub.add_edge(zx_hub.add_boundary(True), hu, fault_free=False)
        zx_hub.add_edge(hv, zx_hub.add_boundary(False), fault_free=False)
    e_hub = zx_hub.add_edge(hu, hv, fault_free=False)
    try:
        apply_zx_rewrite(zx_hub, RewriteStep("fuse_1_fe", {"edge_id": e_hub}))
        raise AssertionError("Expected FUSE1-DEGREE-REFUSED")
    except RefusalError as exc:
        assert exc.code == "FUSE1-DEGREE-REFUSED"
    try:
        apply_zx_rewrite(zx_hub, RewriteStep("ideal_fuse_fe", {"edge_id": e_hub}))
        raise AssertionError("Expected NON-IDEAL-EDGE-IN-IDEAL-RULE")
    except RefusalError as exc:
        assert exc.code == "NON-IDEAL-EDGE-IN-IDEAL-RULE"

    # N5: Non-Clifford QASM T gate refused
    try:
        circuit_from_qasm2("OPENQASM 2.0;\nqreg q[1];\nt q[0];\n")
        raise AssertionError("Expected NON-CLIFFORD-OR-UNSUPPORTED-QASM")
    except RefusalError as exc:
        assert exc.code == "NON-CLIFFORD-OR-UNSUPPORTED-QASM"

    return "all 5 planted negative/refusal families (unflagged cat, unguarded CNOT, missing encoder, bad rewrites, non-Clifford QASM) caught"


def _run_mutant_sabotage(
    mutant_id: str,
    rel_file: str,
    target_snippet: str,
    replacement_snippet: str,
) -> str:
    """Copy `faultweb/` to a scratch dir, apply 1-line sabotage, and require battery to FAIL!"""
    with tempfile.TemporaryDirectory() as tmp:
        stage_pkg = Path(tmp) / "faultweb"
        shutil.copytree(PKG_DIR, stage_pkg, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        f_path = stage_pkg / rel_file
        src = f_path.read_text(encoding="utf-8")
        if target_snippet not in src:
            raise AssertionError(f"Mutant {mutant_id} target snippet missing in {rel_file}")
        f_path.write_text(src.replace(target_snippet, replacement_snippet, 1), encoding="utf-8")

        sub_env = dict(os.environ)
        sub_env["PYTHONPATH"] = str(tmp)
        proc = subprocess.run(
            [
                sys.executable,
                "-c",
                "import sys; from faultweb.battery import run_battery; "
                "rep = run_battery(include_mutations=False); sys.exit(0 if rep['passed'] else 1)",
            ],
            cwd=tmp,
            env=sub_env,
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode == 0:
            raise AssertionError(f"Mutant {mutant_id} SURVIVED (expected failure!)")
        return f"sabotage in {rel_file} killed by core battery (rc={proc.returncode})"


def leg_m1_sabotage_grouped_y_weight() -> str:
    # Sabotage grouped-fault subspace expansion so Y = XZ (sub == 3) is omitted from weight-1!
    return _run_mutant_sabotage(
        "M1_grouped_Y_weight",
        "equivalence.py",
        "for sub in range(1, 1 << width):",
        "for sub in (1,):",
    )


def leg_m2_sabotage_pi_half_web() -> str:
    # Sabotage +/-pi/2 spider opposite-colour parity coupling in `pauliweb.py`!
    return _run_mutant_sabotage(
        "M2_pi_half_spider",
        "pauliweb.py",
        "if nd.phase_k & 1:",
        "if False:",
    )


def leg_m3_sabotage_cnot_extraction_guard() -> str:
    # Sabotage Proposition 6.8 `UNGUARDED-CNOT-EXTRACTION` gate in `diagram.py`!
    return _run_mutant_sabotage(
        "M3_unguarded_cnot",
        "diagram.py",
        "if n_mixed_gadgets < 4 * len(unguarded_cxs):",
        "if False:",
    )


def leg_m4_sabotage_outward_rounding() -> str:
    # Sabotage outward upper-bound rounding from `ROUND_CEILING` to `ROUND_FLOOR` in `contracts.py`!
    return _run_mutant_sabotage(
        "M4_outward_rounding",
        "contracts.py",
        "_OUT_HI = Context(prec=OUT_PREC, rounding=ROUND_CEILING)",
        "_OUT_HI = Context(prec=OUT_PREC, rounding=ROUND_FLOOR)",
    )


def run_battery(include_mutations: bool = True) -> dict[str, Any]:
    legs_spec: list[tuple[str, Callable[[], str]]] = [
        ("L1_purity_zero_float", leg_l1_purity_scan),
        ("L2_cat_and_64q_proof", leg_l2_cat_and_64q_proof),
        ("L3_shor_and_outcome_map", leg_l3_shor_and_outcome_map),
        ("L4_steane_and_flag_decoder", leg_l4_steane_and_flag_decoder),
        ("L5_contracts_examples_cli", leg_l5_contracts_examples_and_cli),
        ("N1_N5_negatives_and_refusals", leg_negatives_and_refusals),
    ]
    if include_mutations:
        legs_spec.extend(
            [
                ("M1_mutant_grouped_weight", leg_m1_sabotage_grouped_y_weight),
                ("M2_mutant_pi_half_web", leg_m2_sabotage_pi_half_web),
                ("M3_mutant_cnot_guard", leg_m3_sabotage_cnot_extraction_guard),
                ("M4_mutant_outward_rounding", leg_m4_sabotage_outward_rounding),
            ]
        )

    results: list[dict[str, Any]] = []
    passed_cnt = 0
    for leg_id, fn in legs_spec:
        t0 = time.perf_counter_ns()
        try:
            summary = fn()
            ok = True
        except Exception as exc:  # noqa: BLE001
            summary = f"{type(exc).__name__}: {exc}"
            ok = False
        dt_ms = (time.perf_counter_ns() - t0) // 1_000_000
        if ok:
            passed_cnt += 1
        results.append(
            {
                "id": leg_id,
                "passed": ok,
                "elapsed_ms": dt_ms,
                "summary": summary,
            }
        )

    return {
        "package": "faultweb",
        "passed": passed_cnt == len(results),
        "passed_count": passed_cnt,
        "total_count": len(results),
        "legs": results,
    }
