#!/usr/bin/env python3
"""Cold-clone self-certification, sabotage-refusal, and mutation battery for `ftqc_magic`."""

from __future__ import annotations

import ast
from decimal import Decimal
from fractions import Fraction
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Callable

PKG_DIR = Path(__file__).resolve().parent
TOOLS_DIR = PKG_DIR.parent
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

from ftqc_magic import (  # noqa: E402
    CATALOG_NAMES,
    BinaryFactory,
    CliffordTableau,
    CorrelationMode,
    GF2Ext,
    LitinskiBlockSpec,
    MagicRefusalError,
    PhasePolynomialTarget,
    PolyQ,
    QSqrt2,
    SignedStabilizerFactory,
    ZOmega8,
    ZSqrt2,
    adaptive_escape_spacetime_volume,
    analytical_lower_bound_n,
    analyze_cultivation_double_check,
    apply_catalyzed_ccz_to_2t,
    apply_phase_catalysis_ladder,
    audit_ppr_check_frame,
    build_canonical_t_factory,
    build_gf8_reed_solomon_ccz_factory,
    build_itogawa_short_check_schedule,
    build_msc3_double_check_schedule,
    build_receipt,
    build_surface_escape_dem,
    build_tricycle_code,
    caliper_worst_case_gap,
    certify_active_qubits_sat,
    colour_code_t_conditional_fidelity,
    compare_global_vs_local_concatenation,
    compose_magic_pipeline,
    compute_aod_shuttle_budget,
    compress_active_qubits,
    from_borrowed_identity,
    get_catalog_factory,
    parent_capacity,
    parity_lift,
    puncture_base_code,
    round_down_str,
    round_up_str,
    search_symmetric_sn_family,
    suggest_pipelines,
    to_compose_ancilla_contract,
    to_compose_logical_t_summary,
    to_extractors_magic_record,
    to_uncompute_pricing_kwargs,
    verify_banked_spatial_nogo,
    verify_determinant_to_2ccz_collapse,
    verify_gidney_fowler_c2t_identity,
    verify_puncture_downward_closure,
    verify_quadcycle_jump_and_singleshot,
    verify_receipt,
    verify_tricycle_ccz_schedule,
    verify_zx_clifford_frame,
)


def _assert_raises_refusal(expected_code: str, fn: Callable[[], object]) -> MagicRefusalError:
    try:
        fn()
    except MagicRefusalError as exc:
        if exc.code != expected_code:
            raise AssertionError(f"Expected refusal {expected_code!r}, got {exc.code!r}: {exc}") from exc
        return exc
    raise AssertionError(f"Expected MagicRefusalError({expected_code!r}), but call succeeded")


def _cadical_available() -> bool:
    env_bin = os.environ.get("FTQC_CODES_CADICAL")
    if env_bin and Path(env_bin).is_file():
        return True
    return shutil.which("cadical") is not None


def leg_l0_purity_and_pinned_receipts() -> str:
    py_files = sorted(p for p in PKG_DIR.rglob("*.py") if "__pycache__" not in p.parts)
    if not py_files:
        raise AssertionError("No Python files discovered in ftqc_magic")
    for path in py_files:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, (float, complex)):
                raise AssertionError(f"Forbidden float/complex literal in {path.name}:{getattr(node, 'lineno', 0)}")
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in ("float", "complex"):
                raise AssertionError(f"Forbidden {node.func.id}() call in {path.name}:{getattr(node, 'lineno', 0)}")

    ex_dir = PKG_DIR / "examples"
    spec_paths = sorted(p for p in ex_dir.glob("*.json") if not p.name.endswith(".receipt.json"))
    if len(spec_paths) < 4:
        raise AssertionError(f"Expected >= 4 packaged example specs, found {len(spec_paths)}")
    for sp in spec_paths:
        rp = ex_dir / f"{sp.stem}.receipt.json"
        if not rp.is_file():
            raise AssertionError(f"Missing banked receipt {rp.name} for {sp.name}")
        spec = json.loads(sp.read_text(encoding="utf-8"))
        rec = json.loads(rp.read_text(encoding="utf-8"))
        vres = verify_receipt(spec, rec, strict_source_sha=True)
        assert vres["verified"] is True
    return f"scanned={len(py_files)}_py, verified={len(spec_paths)}_receipts"


def leg_l1_rings_and_ir() -> str:
    s2 = ZSqrt2(0, 1)
    assert s2 * s2 == ZSqrt2(2, 0)
    assert ZSqrt2(3, -2) > 0 and ZSqrt2(2, -2) < 0

    q = QSqrt2(Fraction(1, 2), Fraction(1, 4))
    assert q * q.inv() == QSqrt2(Fraction(1), Fraction(0))
    lo, hi = q.enclosure_fractions(30)
    assert lo < hi and (QSqrt2(lo) <= q <= QSqrt2(hi))

    # Outward Decimal rounding soundness (tested by M4!):
    up_third = Decimal(round_up_str(Fraction(1, 3), 12))
    dn_third = Decimal(round_down_str(Fraction(1, 3), 12))
    assert up_third * 3 > Decimal(1)
    assert dn_third * 3 < Decimal(1)

    w = ZOmega8.omega_pow(1)
    assert w * ZOmega8.omega_pow(7) == ZOmega8.one()
    assert (ZOmega8.one() + w).abs_sq() == ZSqrt2(2, 1)

    poly = PolyQ.from_weight_enumerator(3, [1, 0, 3, 1])
    assert poly.eval_exact(Fraction(1, 5)) == Fraction(64 + 12 + 1, 125)

    t_ccz = PhasePolynomialTarget.named("CCZ")
    t_cs2 = PhasePolynomialTarget.from_jws_polynomial(3, "0 1+0 2")
    assert t_ccz.pauli_abs_sq_spectrum() == t_cs2.pauli_abs_sq_spectrum()
    t_ccz_cz = PhasePolynomialTarget(
        k=3, phases=tuple(((4 if x == 7 else 0) + (4 if (x & 3) == 3 else 0)) % 8 for x in range(8))
    )
    assert t_ccz.equivalent_up_to_diagonal_clifford(t_ccz_cz)
    return "ZSqrt2+QSqrt2+ZOmega8+PolyQ+orbit_ok"


def leg_l2_binary_galois_and_transforms() -> str:
    assert len(CATALOG_NAMES) >= 17
    bk15 = get_catalog_factory("bk_15_1_3_t")
    assert isinstance(bk15, BinaryFactory)
    c_ps = bk15.to_source_contract("1/1000", mode="postselect")
    c_ec = bk15.to_source_contract("1/1000", mode="ec")
    assert c_ps.distance == 3 and c_ps.leading_order_coeff == 35
    assert c_ec.distance == 2 and c_ec.leading_order_coeff == 105

    bh14 = get_catalog_factory("bh_14_2_2_t")
    assert isinstance(bh14, BinaryFactory)
    c14 = bh14.to_source_contract("1/1000")
    assert c14.distance == 2 and c14.leading_order_coeff == 7 and c14.correlation_mode == CorrelationMode.BATCH_CORRELATED

    bgm36 = get_catalog_factory("bgm_36_3_3_cs2_ccz")
    c36 = bgm36.to_source_contract("1/1000")
    assert c36.distance == 3 and c36.leading_order_coeff == 120

    bgm48 = get_catalog_factory("bgm_48_3_4_ccz_2level")
    c48 = bgm48.to_source_contract("1/1000")
    assert c48.distance == 4 and c48.leading_order_coeff == 864

    jws63 = get_catalog_factory("jws_63_2_3_t2")
    c63 = jws63.to_source_contract("1/1000")
    assert c63.distance == 3 and c63.k_outputs == 2

    # Borrowed-identity & transforms:
    assert parent_capacity(bk15.check_rows, bk15.n)["capacity_kappa"] == 1
    id16_cols = [1 | (s << 1) for s in range(16)]
    borrowed = from_borrowed_identity(1, 4, id16_cols, [1], name="borrowed_bk15")
    assert borrowed.n == 15 and borrowed.certify_distance().distance == 3

    lifted = parity_lift(get_catalog_factory("jvbs_47_3_3_ccz"))
    assert lifted.n == 48 and lifted.certify_distance().distance == 4

    rm16_rows = [
        (1 << 16) - 1,
        *[sum((1 << x) for x in range(16) if (x >> b) & 1) for b in range(4)],
    ]
    punc = puncture_base_code(rm16_rows, 16, [0], name="punctured_rm15")
    assert punc.n == 15 and punc.certify_distance().distance == 3
    dc = verify_puncture_downward_closure(rm16_rows, 16, [0], [0, 1])
    assert dc["downward_closed"] is True

    # Galois self-dual trace binarization + 3-twirl concatenation:
    for s_deg in (2, 3, 4):
        fld = GF2Ext(s_deg)
        basis = fld.self_dual_basis()
        for i in range(s_deg):
            for j in range(s_deg):
                assert fld.trace(fld.mul(basis[i], basis[j])) == (1 if i == j else 0)

    gf8 = get_catalog_factory("gf8_rs_8_1_2_ccz")
    gf4 = get_catalog_factory("gf4_hermitian_ag_8_1_2_cs")
    gprw64 = get_catalog_factory("gprw_64_6_4_2ccz")
    assert gf8.certify_distance().distance >= 2
    assert gf4.certify_distance().distance >= 2
    c_gprw = gprw64.to_source_contract("1/1000")
    assert c_gprw.distance == 4 and c_gprw.k_outputs == 6

    cmp_twirl = compare_global_vs_local_concatenation(
        gprw64, build_gf8_reed_solomon_ccz_factory(), "1/1000"
    )
    assert Fraction(cmp_twirl["local_saves_factor"]) > 1
    return "binary+borrowed+galois+3twirl_ok"


def leg_l3_collapse_and_lower_bounds() -> str:
    bgm9 = get_catalog_factory("bgm_9_1_2_t")
    assert isinstance(bgm9, SignedStabilizerFactory)
    c9, rep9 = bgm9.certify_collapse(PhasePolynomialTarget.named("T"), claimed_distance=2, p_input="1/1000")
    assert c9.distance == 2
    assert rep9["rational_polynomials"]["ideal_acceptance"] == "1/128"
    assert rep9["rational_polynomials"]["leading_error_ratio"] == "27/2"
    assert c9.optimality_lower_bound_n == 9

    bgm6 = get_catalog_factory("bgm_6_4_1_cccz")
    c6, _ = bgm6.certify_collapse(PhasePolynomialTarget.named("CCCZ"), claimed_distance=1, p_input="1/1000")
    assert c6.distance >= 1 and c6.optimality_lower_bound_n == 6

    bgm3 = get_catalog_factory("bgm_3_2_1_cs")
    c3, _ = bgm3.certify_collapse(PhasePolynomialTarget.named("CS"), claimed_distance=1, p_input="1/1000")
    assert c3.distance >= 1 and c3.optimality_lower_bound_n == 3

    assert analytical_lower_bound_n(PhasePolynomialTarget.named("CCZ"), 2)["lower_bound_n"] == 8
    assert analytical_lower_bound_n(PhasePolynomialTarget.named("CCZ"), 3)["lower_bound_n"] == 10
    return "9T->1T(1/128,27/2p^2)+6T->CCCZ+lower_bounds_ok"


def leg_l4_recycling_and_spatial_nogo() -> str:
    bk15 = get_catalog_factory("bk_15_1_3_t")
    r15 = compress_active_qubits(bk15)
    assert r15.original_wires == 5 and r15.active_qubits == 4

    bh49 = get_catalog_factory("bh_49_1_5_t")
    r49 = compress_active_qubits(bh49)
    assert r49.original_wires == 14 and r49.active_qubits <= 5

    f5 = build_canonical_t_factory(5)
    assert f5.n == 15 and f5.certify_distance().distance == 3
    sn9 = search_symmetric_sn_family(9, 3)
    assert sn9 and sn9[0]["n"] == 45 and sn9[0]["d"] == 3

    nogo_t4 = verify_banked_spatial_nogo(4, "T", 3)
    nogo_ccz5 = verify_banked_spatial_nogo(5, "CCZ", 3)
    assert nogo_t4["verdict"] == "UNSAT-NOGO-CERTIFIED"
    assert nogo_ccz5["verdict"] == "UNSAT-NOGO-CERTIFIED"

    ext_tag = "banked_lrat"
    if _cadical_available():
        sat_opt = certify_active_qubits_sat(bk15, 3)
        assert sat_opt["verdict"] == "UNSAT-OPTIMALITY-PROVED" and sat_opt["certified_min_active_qubits"] == 4
        ext_tag = "banked+live_cadical_lrat"
    return f"15->1(4q),49->1({r49.active_qubits}q),S9(45T),{ext_tag}"


def leg_l5_cultivation_clifford_hooks() -> str:
    assert colour_code_t_conditional_fidelity(CliffordTableau.identity(7)) == 1
    assert (
        colour_code_t_conditional_fidelity(
            CliffordTableau.from_gates(7, [("Z", q) for q in range(7)])
        )
        == 0
    )

    c_unflag, meta_unflag = analyze_cultivation_double_check(
        build_msc3_double_check_schedule(z_flagged=False), p_phys="1/1000"
    )
    assert c_unflag.distance == 2
    assert meta_unflag["d_S_proxy"] == 3
    assert meta_unflag["malignant_counts_T"][2] == 4
    assert meta_unflag["malignant_counts_S"][2] == 0
    assert meta_unflag["malignant_counts_S"][3] == 82

    c_flag, meta_flag = analyze_cultivation_double_check(
        build_msc3_double_check_schedule(z_flagged=True),
        claimed_t_distance=3,
        p_phys="1/1000",
        crosscheck_faults_and_stim=True,
    )
    assert c_flag.distance == 3 and meta_flag["malignant_counts_T"][2] == 0

    c_ito, _ = analyze_cultivation_double_check(
        build_itogawa_short_check_schedule(), claimed_t_distance=3, p_phys="1/1000"
    )
    assert c_ito.distance == 3

    s_xchk = meta_flag.get("external_crosschecks") or {}
    aud_ent = s_xchk.get("ftqc_faults_audit")
    stim_ent = s_xchk.get("stim")
    audit_ok = bool(aud_ent.get("agrees_d_ge_3", False)) if isinstance(aud_ent, dict) else False
    stim_ok = bool(stim_ent.get("agrees_d_eq_3", False)) if isinstance(stim_ent, dict) else False
    return f"unflag(d_T=2,mal=4),z_flag(d_T=3),itogawa(d_T=3),audit={audit_ok},stim={stim_ok}"


def leg_l6_caliper_and_qldpc() -> str:
    edem = build_surface_escape_dem(3, 4, 2, prior="1/1000", use_stim_if_available=True)
    gap = caliper_worst_case_gap(edem, committed_syndrome=1, mode="caliper", crosscheck_blossom=True)
    assert gap["caliper_worst_soft_gap"] <= gap["naive_soft_gap"]
    assert gap["gap_erosion"] > 0

    vol = adaptive_escape_spacetime_volume([15, 35, 65, 98], ["4/5", "4/5", "4/5"], "1/2")
    assert Decimal(vol["volume_reduction_factor_lo"]) > Decimal(1)

    tri = build_tricycle_code(3, 3, 3)
    cc_rep = tri.verify_chain_complex()
    assert cc_rep["css_commutes"] is True and cc_rep["x_metacheck_valid"] is True
    fo_rep = verify_tricycle_ccz_schedule(tri, ccz_schedule="fanout_depth1", claimed_circuit_distance=3)
    assert fo_rep["effective_circuit_distance"] == 3 and fo_rep["physical_ccz_depth"] == 1

    qj = verify_quadcycle_jump_and_singleshot(L=2, prior="1/1000")
    assert qj["css_commutes"] and qj["meta_x_valid"] and qj["meta_z_valid"]

    det = verify_determinant_to_2ccz_collapse()
    assert det["verified_all_64_basis_states_both_branches"] is True
    aod = compute_aod_shuttle_budget([(1, 0, 0), (1, 1, 0)], [4, 4, 4])
    assert aod["theta_aod_shuttle_moves"] == 2 + 4
    return f"caliper_erosion={gap['gap_erosion']},blossom={gap['blossom_crosschecks']},quadcycle_ok"


def leg_l7_pipeline_catalysis_bridges_and_cli() -> str:
    c2t_id = verify_gidney_fowler_c2t_identity()
    assert c2t_id["verified_exact_z_omega8"] is True

    cube_c = get_catalog_factory("cube_8_3_2_ccz").to_source_contract("1/1000")
    c2t_c, _ = apply_catalyzed_ccz_to_2t(cube_c)
    assert c2t_c.correlation_mode == CorrelationMode.CATALYST_POISONED and c2t_c.k_outputs == 2

    t_c = get_catalog_factory("bk_15_1_3_t").to_source_contract("1/1000")
    lad_c, _ = apply_phase_catalysis_ladder(cube_c, t_c, level_ell=4)
    assert lad_c.correlation_mode == CorrelationMode.CATALYST_POISONED

    lit = LitinskiBlockSpec(d_X=7, d_Z=13, d_m=7)
    assert lit.physical_qubits_for_patches(1, 4) > 0

    pareto_t = suggest_pipelines(target_gate="T", p_phys="1/1000", target_infidelity="1/100000000")
    pareto_ccz = suggest_pipelines(target_gate="CCZ", p_phys="1/1000", target_infidelity="1/100000")
    assert len(pareto_t) >= 2 and len(pareto_ccz) >= 2

    # Exercise explicit interleaving & joint_3twirl multi-output pipeline modes:
    c_j3t, _ = compose_magic_pipeline(
        [
            {"kind": "raw"},
            {"kind": "distillation", "factory": "bh_14_2_2_t"},
            {"kind": "distillation", "factory": "bk_15_1_3_t", "correlation_propagation": "joint_3twirl"},
        ]
    )
    assert c_j3t.distance == 6 and "joint_3twirl_kernel_norm" in c_j3t.metadata

    # Bridges:
    f_aud = audit_ppr_check_frame(get_catalog_factory("bk_15_1_3_t"))
    zx_aud = verify_zx_clifford_frame()
    if f_aud["available"]:
        assert f_aud["verified"] is True
    if zx_aud["available"]:
        assert zx_aud["cnot_cancel_equiv"] is True

    comp_anc = to_compose_ancilla_contract(t_c, atom_tag="t1")
    comp_sum = to_compose_logical_t_summary(t_c)
    ext_rec = to_extractors_magic_record(t_c)
    unc_kw = to_uncompute_pricing_kwargs(t_c, ccz_contract=cube_c)
    assert "independent_contract" in comp_anc and comp_sum["distance"] == 3
    assert ext_rec["fault_distance"] == 3 and unc_kw["t_factory_wire_ticks"] >= 1

    # Subprocess CLI round-trip (`check` + `verify`):
    spec_p = PKG_DIR / "examples" / "bgm_9_to_1_nontriorthogonal.json"
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(TOOLS_DIR), *([env["PYTHONPATH"]] if env.get("PYTHONPATH") else [])])
    with tempfile.TemporaryDirectory() as tmp:
        out_p = Path(tmp) / "cli.receipt.json"
        p1 = subprocess.run(
            [sys.executable, "-m", "ftqc_magic", "check", str(spec_p), "-o", str(out_p)],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        assert p1.returncode == 0 and out_p.is_file()
        p2 = subprocess.run(
            [sys.executable, "-m", "ftqc_magic", "verify", str(spec_p), str(out_p)],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        assert p2.returncode == 0 and json.loads(p2.stdout)["verified"] is True
    return f"pareto_T={len(pareto_t)},pareto_CCZ={len(pareto_ccz)},cli_verified"


def leg_n1_n6_sabotage_controls() -> str:
    # N1: Unflagged cultivation claiming d_T = 3
    _assert_raises_refusal(
        "UNFLAGGED-CULTIVATION-HOOK-DEGRADATION",
        lambda: analyze_cultivation_double_check(
            build_msc3_double_check_schedule(z_flagged=False), claimed_t_distance=3
        ),
    )

    # N2: Catalyst poisoning fed into downstream distillation or independent compose tail bound
    _assert_raises_refusal(
        "CATALYST-POISONING-HAZARD",
        lambda: compose_magic_pipeline(
            [
                {"kind": "raw"},
                {"kind": "distillation", "factory": "cube_8_3_2_ccz"},
                {"kind": "catalyzed_c2t"},
                {"kind": "distillation", "factory": "bk_15_1_3_t"},
            ]
        ),
    )
    cat_stream, _ = apply_catalyzed_ccz_to_2t(get_catalog_factory("cube_8_3_2_ccz").to_source_contract("1/1000"))
    _assert_raises_refusal("CATALYST-POISONING-HAZARD", lambda: to_compose_ancilla_contract(cat_stream))

    # N3: Multi-output correlated batch fed into distillation without interleaving or joint 3-twirl
    _assert_raises_refusal(
        "MULTI-OUTPUT-CORRELATION-HAZARD",
        lambda: compose_magic_pipeline(
            [
                {"kind": "raw"},
                {"kind": "distillation", "factory": "bh_14_2_2_t"},
                {"kind": "distillation", "factory": "bk_15_1_3_t"},
            ]
        ),
    )
    batch_c = get_catalog_factory("bh_14_2_2_t").to_source_contract("1/1000")
    _assert_raises_refusal("MULTI-OUTPUT-CORRELATION-HAZARD", lambda: to_compose_ancilla_contract(batch_c))

    # N4: Naive depth-2 transversal CCZ on 3D tricycle claiming full distance d = 3
    _assert_raises_refusal(
        "DEPTH2-TRANSVERSAL-CCZ-HALVING",
        lambda: verify_tricycle_ccz_schedule(
            build_tricycle_code(3, 3, 3), ccz_schedule="naive_depth2", claimed_circuit_distance=3
        ),
    )

    # N5: Naive terminal gap on open-boundary mid-circuit escape DEM
    _assert_raises_refusal(
        "OPEN-BOUNDARY-NAIVE-GAP-HAZARD",
        lambda: caliper_worst_case_gap(
            build_surface_escape_dem(3, 4, 2), committed_syndrome=1, mode="naive_terminal_gap"
        ),
    )

    # N6: Non-collapsing stabilizer code + non-Clifford check leak + tampered receipt
    bad_stab = SignedStabilizerFactory(n=3, k=1, stabilizers=("ZZI", "IZZ"), name="bad_rep3")
    _assert_raises_refusal(
        "NON-COLLAPSE-OR-PARITY-VIOLATION",
        lambda: bad_stab.certify_collapse(PhasePolynomialTarget.named("T"), claimed_distance=2),
    )
    bk15 = get_catalog_factory("bk_15_1_3_t")
    bad_bin = BinaryFactory(
        k=1,
        r=bk15.r,
        n=bk15.n,
        output_rows=bk15.output_rows,
        check_rows=(bk15.check_rows[0] ^ 1, *bk15.check_rows[1:]),
        name="bad_bk15",
    )
    _assert_raises_refusal("NON-CLIFFORD-CHECK-LEAK", lambda: bad_bin.verify_algebraic_structure())

    spec = {"kind": "signed_stabilizer", "catalog": "bgm_9_1_2_t", "target": "T", "claimed_distance": 2}
    rec = build_receipt(spec)
    tampered = json.loads(json.dumps(rec))
    tampered["contract"]["distance"] = 99
    _assert_raises_refusal("RECEIPT-VERIFICATION-MISMATCH", lambda: verify_receipt(spec, tampered))
    return "N1..N6_all_refused_cleanly"


def _run_mutant_probe() -> None:
    """Fast (<0.15s) probe executed inside isolated mutated copies for M1..M4."""
    # Covers M4 (ROUND_CEILING soundness):
    assert Decimal(round_up_str(Fraction(1, 3), 12)) * 3 > Decimal(1)
    # Covers M1 (limits.py non-CSS Pauli phase + Cauchy-Schwarz collapse):
    c9, rep9 = get_catalog_factory("bgm_9_1_2_t").certify_collapse(
        PhasePolynomialTarget.named("T"), claimed_distance=2, p_input="1/1000"
    )
    assert c9.distance == 2 and rep9["rational_polynomials"]["ideal_acceptance"] == "1/128"
    # Covers M2 (cultivate.py Lemma 4 quadratic sign chi(u)):
    c_unflag, meta_u = analyze_cultivation_double_check(
        build_msc3_double_check_schedule(z_flagged=False), p_phys="1/1000"
    )
    assert c_unflag.distance == 2 and meta_u["malignant_counts_T"][2] == 4
    # Covers M3 (pipeline.py CATALYST_POISONED guard):
    _assert_raises_refusal(
        "CATALYST-POISONING-HAZARD",
        lambda: compose_magic_pipeline(
            [
                {"kind": "raw"},
                {"kind": "distillation", "factory": "cube_8_3_2_ccz"},
                {"kind": "catalyzed_c2t"},
                {"kind": "distillation", "factory": "bk_15_1_3_t"},
            ]
        ),
    )


def leg_m1_m4_mutant_kills() -> str:
    # First prove the baseline unmutated probe succeeds cleanly!
    _run_mutant_probe()
    mutants = [
        (
            "M1_limits_pauli_phase",
            "limits.py",
            "return (x1 ^ x2, z1 ^ z2, (p1 + p2 + 2 * comm_flips) % 4)",
            "return (x1 ^ x2, z1 ^ z2, (p1 + p2) % 4)",
        ),
        (
            "M2_cultivate_lemma4_chi",
            "cultivate.py",
            "return (((num // 2) % 2) ^ cross_zx) & 1",
            "return 0",
        ),
        (
            "M3_pipeline_catalyst_guard",
            "pipeline.py",
            "if current.correlation_mode == CorrelationMode.CATALYST_POISONED:",
            "if False and current.correlation_mode == CorrelationMode.CATALYST_POISONED:",
        ),
        (
            "M4_rings_outward_rounding",
            "rings.py",
            "return _fraction_to_decimal_exact_bound(Fraction(x), sig, ROUND_CEILING)",
            "return _fraction_to_decimal_exact_bound(Fraction(x), sig, ROUND_FLOOR)",
        ),
    ]
    killed = 0
    for label, rel_file, target, replacement in mutants:
        with tempfile.TemporaryDirectory() as tmp:
            dst_pkg = Path(tmp) / "ftqc_magic"
            shutil.copytree(
                PKG_DIR,
                dst_pkg,
                ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "examples"),
            )
            fpath = dst_pkg / rel_file
            src = fpath.read_text(encoding="utf-8")
            if target not in src:
                raise AssertionError(f"Mutant {label} target string missing in {rel_file}")
            fpath.write_text(src.replace(target, replacement, 1), encoding="utf-8")

            env = dict(os.environ)
            # Prepend `tmp` BEFORE `TOOLS_DIR` so `import ftqc_magic` loads the mutated copy!
            env["PYTHONPATH"] = os.pathsep.join([str(tmp), str(TOOLS_DIR)])
            proc = subprocess.run(
                [sys.executable, str(dst_pkg / "selftest.py"), "--mutant-probe"],
                env=env,
                capture_output=True,
                text=True,
                timeout=20,
            )
            if proc.returncode == 0:
                raise AssertionError(f"SURVIVING MUTANT {label}: probe exited with rc=0!")
            killed += 1
    return f"killed={killed}/4_mutants"


def main(argv: list[str] | None = None) -> int:
    args = argv if argv is not None else sys.argv[1:]
    if "--mutant-probe" in args:
        _run_mutant_probe()
        return 0

    legs: list[tuple[str, Callable[[], str]]] = [
        ("L0_purity_and_pinned_receipts", leg_l0_purity_and_pinned_receipts),
        ("L1_rings_and_ir", leg_l1_rings_and_ir),
        ("L2_binary_galois_and_transforms", leg_l2_binary_galois_and_transforms),
        ("L3_collapse_and_lower_bounds", leg_l3_collapse_and_lower_bounds),
        ("L4_recycling_and_spatial_nogo", leg_l4_recycling_and_spatial_nogo),
        ("L5_cultivation_clifford_hooks", leg_l5_cultivation_clifford_hooks),
        ("L6_caliper_and_qldpc", leg_l6_caliper_and_qldpc),
        ("L7_pipeline_catalysis_bridges_and_cli", leg_l7_pipeline_catalysis_bridges_and_cli),
        ("N1_N6_sabotage_controls", leg_n1_n6_sabotage_controls),
        ("M1_M4_mutant_kills", leg_m1_m4_mutant_kills),
    ]
    t0 = time.monotonic_ns()
    for name, fn in legs:
        lt0 = time.monotonic_ns()
        detail = fn()
        elapsed_ms = (time.monotonic_ns() - lt0) // 1_000_000
        print(f"[PASS] {name} ({elapsed_ms} ms): {detail}")
    total_ms = (time.monotonic_ns() - t0) // 1_000_000
    print(f"SELF-TEST: ALL PASS (rc 0) [{len(legs)} legs in {total_ms} ms]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
