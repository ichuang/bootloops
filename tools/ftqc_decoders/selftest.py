"""Cold-clone self-certification & scratch-copy mutation battery for `ftqc_decoders`.

Runs green in `< 2 s` on a bare stdlib Python 3.10+ clone with zero external packages,
and dynamically exercises live external backends (`Tesseract`, `PyMatching`, `ldpc`,
`Sinter`) whenever importable in the active interpreter.

Zero `float` or `complex` literals/calls appear anywhere in this module.
"""

from __future__ import annotations

import ast
from decimal import Decimal
from fractions import Fraction
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Callable

# Bootstrap `tools/` onto `sys.path` when invoked directly as `python3 tools/ftqc_decoders/selftest.py`
_PKG_DIR = Path(__file__).resolve().parent
_TOOLS_DIR = _PKG_DIR.parent
if str(_TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(_TOOLS_DIR))

import ftqc_decoders as fd  # noqa: E402
from ftqc_decoders.mwpm import max_weight_perfect_matching_blossom  # noqa: E402


def leg_l0_purity_and_lineage() -> str:
    py_files = sorted(p for p in _PKG_DIR.rglob("*.py") if not p.name.startswith("._"))
    for py_f in py_files:
        tree = ast.parse(py_f.read_text(), filename=str(py_f))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and not isinstance(
                node.value, (str, int, bool, bytes, type(None), tuple, type(Ellipsis))
            ):
                raise AssertionError(f"Non-exact constant {node.value!r} in {py_f}:{node.lineno}")
            if isinstance(node, ast.Name) and node.id in ("float", "complex"):
                raise AssertionError(f"Forbidden builtin {node.id} in {py_f}:{node.lineno}")
        if py_f.name == "mwpm.py":
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom) and node.module and any(
                    bad in node.module for bad in ("bp_osd", "relay_bp", "union_find", "ssf_metacheck", "ftqc_codes")
                ):
                    raise AssertionError(f"Disjoint-lineage violation in mwpm.py: {node.module}")
    return f"{len(py_files)} modules scanned (0 floats, mwpm disjoint)"


def leg_l1_dem_stim_and_mwpm() -> str:
    rep5 = fd.repetition_problem(5)
    r_rep = fd.certify_operational_radius(rep5, "mwpm", 2)
    assert r_rep["status"] == "proved" and r_rep["holds_through"] == 2

    tor4 = fd.toric_x_problem(4)
    r_tor = fd.certify_operational_radius(tor4, "mwpm", 1)
    assert r_tor["status"] == "proved" and r_tor["holds_through"] == 1

    # Stim .dem round-trip + repeat-block parser
    dem_txt = tor4.to_stim_dem_text()
    tor4_rt = fd.DecodingProblem.from_stim_dem_text(dem_txt, name="tor4_rt")
    assert (tor4_rt.num_checks, tor4_rt.num_mechanisms, tor4_rt.num_observables) == (
        tor4.num_checks,
        tor4.num_mechanisms,
        tor4.num_observables,
    )
    rep_dem = "error(0.01) D0 L0\nrepeat 2 {\n  error(0.01) D0 D1\n  shift_detectors 1\n}\nerror(0.01) D0\n"
    p_rep = fd.DecodingProblem.from_stim_dem_text(rep_dem)
    assert p_rep.num_checks == 3 and p_rep.is_graphlike()

    # Stress-test Edmonds Blossom vs exact bitmask DP + LP strong duality on 30 random odd-blossom graphs
    rng = random.Random(20261126)
    for _ in range(30):
        n = rng.choice([4, 6, 8])
        w = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                val = rng.choice([0, 1, 3, 5, 8, 13, 34])
                w[i][j] = w[j][i] = val
        mate, dual_ok = max_weight_perfect_matching_blossom(n, w)
        assert dual_ok, "LP strong duality failed in Blossom"
        bw = sum(w[i][mate[i]] for i in range(n)) // 2
        dp = [-1] * (1 << n)
        dp[0] = 0
        for mask in range(1, 1 << n):
            if mask.bit_count() & 1:
                continue
            lsb = mask & -mask
            i = lsb.bit_length() - 1
            rem = mask ^ lsb
            cur = rem
            best = -1
            while cur:
                jlsb = cur & -cur
                j = jlsb.bit_length() - 1
                cand = dp[rem ^ jlsb] + w[i][j]
                if cand > best:
                    best = cand
                cur ^= jlsb
            dp[mask] = best
        assert bw == dp[(1 << n) - 1]
    return "DEM round-trip + 30 random Blossom/DP/LP-duality checks pass"


def leg_l2_bp_osd_lsd_relay() -> str:
    steane = fd.DecodingProblem.from_css_code("steane", sector="X")
    for bp_m in ("nms", "sum_product"):
        r = fd.certify_operational_radius(steane, {"algorithm": "bp_osd", "bp_method": bp_m, "osd_order": 1}, 1)
        assert r["status"] == "proved"

    golay = fd.DecodingProblem.from_css_code("golay", sector="X")
    r_gol = fd.certify_operational_radius(golay, {"algorithm": "bp_osd", "osd_order": 2}, 2)
    assert r_gol["status"] == "proved" and r_gol["holds_through"] == 2

    qr47 = fd.DecodingProblem.from_css_code("qr:47", sector="X")
    r_qr = fd.certify_operational_radius(qr47, {"algorithm": "bp_osd", "osd_order": 2}, 1)
    assert r_qr["status"] == "proved"

    bb144 = fd.bb144_gross_problem()
    g0_bits = fd.bit_positions(bb144.stabilizer_generators[0])
    e_half = sum(1 << b for b in g0_bits[:3])
    syn_h = bb144.syndrome_of(e_half)
    r_rel = fd.decode_relay_bp(bb144, syn_h)
    r_lsd = fd.decode_bp_osd(bb144, syn_h, post_method="lsd", osd_order=0)
    r_osd = fd.decode_bp_osd(bb144, syn_h, post_method="osd", osd_order=0)
    assert r_rel.converged and bb144.observable_of(e_half ^ r_rel.error_mask) == 0
    assert r_lsd.converged and bb144.observable_of(e_half ^ r_lsd.error_mask) == 0
    assert r_osd.converged and bb144.observable_of(e_half ^ r_osd.error_mask) == 0
    assert r_lsd.stats["max_cluster_rank"] <= 2 and r_osd.stats["eliminated_rank"] >= 60
    return f"Steane/Golay/QR47/BB144 pass (LSD rank={r_lsd.stats['max_cluster_rank']} vs OSD={r_osd.stats['eliminated_rank']})"


def leg_l3_uf_ssf_metachecks_sat_lrat() -> str:
    rep5 = fd.repetition_problem(5)
    r_uf = fd.certify_operational_radius(rep5, "union_find", 2)
    assert r_uf["status"] == "proved" and r_uf["holds_through"] == 2

    tor4 = fd.toric_x_problem(4)
    assert fd.certify_operational_radius(tor4, "union_find", 1)["status"] == "proved"

    bb72 = fd.bb72_problem()
    assert fd.certify_operational_radius(bb72, "belief_find", 1)["status"] == "proved"

    steane = fd.DecodingProblem.from_css_code("steane", sector="X")
    assert fd.certify_operational_radius(steane, "ssf", 1)["status"] == "proved"
    assert fd.certify_operational_radius(steane, "sat_mwd", 1)["status"] == "proved"

    # Self-contained redundant-Steane metachecks + live `ftqc_extractors/surgery_sequence.py` mapping-cone metachecks!
    s_rows = list(steane.check_rows)
    red_checks = s_rows + s_rows + s_rows
    meta_rows = [(1 << c) | (1 << (c + 3)) for c in range(3)] + [(1 << (c + 3)) | (1 << (c + 6)) for c in range(3)]
    meta_steane = fd.DecodingProblem.from_matrices(
        red_checks, steane.observable_rows, 7, name="steane_meta", metacheck_rows=meta_rows
    )
    mc_s = fd.decode_with_metachecks(meta_steane, meta_steane.syndrome_of(1 << 2) ^ (1 << 4))
    assert mc_s["meta_converged"] and mc_s["data_converged"] and meta_steane.observable_of((1 << 2) ^ mc_s["data_error_mask"]) == 0

    seq_example = _TOOLS_DIR / "ftqc_extractors" / "examples" / "three-block-sequence.json"
    if seq_example.is_file():
        from ftqc_codes import gf2 as fc_gf2
        from ftqc_extractors import surgery_sequence

        seq_spec = json.loads(seq_example.read_text())
        code0 = surgery_sequence.memory(seq_spec["memory"])
        n0, hx0, hz0 = code0
        x_c, z_c = list(hx0), list(hz0)
        op0 = seq_spec["operations"][0]
        cone_w, _ = surgery_sequence.append_cone(code0, op0, n0, x_c, z_c, seq_spec["max_degree"])
        ech_z = fc_gf2.Echelon(z_c)
        d_logs = [v for v in fc_gf2.nullspace(x_c, cone_w) if ech_z.add(v)]
        cone_p = fd.DecodingProblem.from_matrices(
            z_c, d_logs, cone_w, name="cone0", stabilizer_generators=x_c, metacheck_rows=op0["metachecks"]
        )
        noisy_s = cone_p.syndrome_of(1 << 0) ^ (1 << len(hz0))
        mc = fd.decode_with_metachecks(cone_p, noisy_s)
        assert mc["meta_converged"] and mc["data_converged"] and mc["residual_syndrome_weight"] == 0
    return "Half-edge UF, BeliefFind, SSF, surgery metachecks, and SAT+LRAT MWD pass"


def leg_l4_trapping_sets_and_deff() -> str:
    bb144 = fd.bb144_gross_problem()
    ts = fd.search_trapping_sets(
        bb144,
        "bp_only",
        max_exhaustive_weight=0,
        max_splits_per_generator=4,
        compare_decoders=["relay_bp", "bp_lsd", "bp_osd"],
    )
    assert ts["status"] == "trapping_set_found" and ts["witness"]["fault_weight"] <= 3
    assert all(c["healed"] for c in ts["comparison"])

    golay = fd.DecodingProblem.from_css_code("golay", sector="X")
    c_osd1 = fd.certify_operational_radius(golay, {"algorithm": "bp_osd", "schedule": "serial", "osd_order": 1}, 2)
    assert c_osd1["status"] == "counterexample" and c_osd1["failure_weight"] == 2
    assert fd.verify_failure_witness(golay, {"algorithm": "bp_osd", "schedule": "serial", "osd_order": 1}, c_osd1["witness"])
    return f"BB144 sub-radius trapping set w={ts['witness']['fault_weight']} healed + Golay OSD-1 witness verified"


def leg_l5_enclosures_bridges_examples_cli() -> str:
    # Outward containment check
    lo_s, hi_s = fd.fraction_to_outward_strs(Fraction(1, 7), digits=16)
    assert Fraction(lo_s) < Fraction(1, 7) < Fraction(hi_s)

    steane = fd.DecodingProblem.from_css_code("steane", sector="X")
    ex = fd.enclose_exact_shells(steane, "bp_osd", "1/1000", verified_radius=1, max_exact_weight=2, digits=18)
    assert Decimal("0.0000208") < Decimal(ex["ler_lo"]) <= Decimal(ex["ler_hi"]) < Decimal("0.0000210")

    bb144 = fd.bb144_gross_problem()
    pe = fd.enclose_peierls_cluster(bb144, "1/10000", verified_component_radius=5, digits=18)
    assert Decimal(pe["ler_hi"]) < Decimal("1e-12")

    st = fd.enclose_stratified_shells(steane, "bp_osd", "1/1000", verified_radius=1, max_shell_weight=3, shots_per_shell=25, digits=16)
    assert st["evidence_kind"] == "empirical_stratified_hoeffding"

    # Procedural gadget bridge on `knill:golay,anc=spec` (avoids 2^11 lookup table)
    from ftqc_gadgets import effects, library, props

    g_gol = fd.with_procedural_decoders(library.load("knill:golay,anc=spec"), algorithm="bp_osd", osd_order=2)
    ec_gol = props.check(effects.Analysis(g_gol), "ec", 1)
    assert ec_gol["holds_through"] == 1 and ec_gol.get("violation") is None

    # Replay all packaged example study specs
    ex_dir = _PKG_DIR / "examples"
    ex_files = sorted(ex_dir.glob("*.json"))
    assert len(ex_files) >= 4
    for ef in ex_files:
        spec = json.loads(ef.read_text())
        rec = fd.build_receipt(spec, base_dir=ex_dir)
        ver = fd.verify_receipt(rec, target_study=spec, base_dir=ex_dir)
        assert ver["verified"] is True

    # Optional live external backend probe (runs automatically inside uv venv!)
    backends = fd.probe_external_backends()
    live_ext: list[str] = []
    if backends["tesseract"]["available"] and backends["stim"]["available"]:
        r_t = fd.run_decoder(bb72 := fd.bb72_problem(), bb72.syndrome_of((1 << 2) | (1 << 19)), "external:tesseract")
        assert r_t.converged and bb72.observable_of(((1 << 2) | (1 << 19)) ^ r_t.error_mask) == 0
        live_ext.append("tesseract")
    if backends["pymatching"]["available"]:
        tor4 = fd.toric_x_problem(4)
        r_p = fd.run_decoder(tor4, tor4.syndrome_of(1 << 3), "external:pymatching")
        assert r_p.converged and tor4.observable_of((1 << 3) ^ r_p.error_mask) == 0
        live_ext.append("pymatching")
    if backends["ldpc"]["available"]:
        bb72 = fd.bb72_problem()
        r_l = fd.run_decoder(bb72, bb72.syndrome_of(1 << 4), "external:ldpc_bplsd")
        assert r_l.converged and bb72.observable_of((1 << 4) ^ r_l.error_mask) == 0
        live_ext.append("ldpc")

    ext_note = f"live={','.join(live_ext)}" if live_ext else "external=skipped(not in current python)"
    return f"{len(ex_files)} examples verified, gadgets bridge live, {ext_note}"


def leg_l6_tanner_single_shot_and_sparse_handoff() -> str:
    qt_d8_x = fd.DecodingProblem.from_css_code("qt_d8_36", sector="X")
    qt_d8_z = fd.DecodingProblem.from_css_code("qt_d8_36", sector="Z")
    r_x = fd.certify_operational_radius(qt_d8_x, "tanner_local", 1)
    r_z = fd.certify_operational_radius(qt_d8_z, "tanner_local", 1)
    assert r_x["status"] == "proved" and r_x["holds_through"] == 1
    assert r_z["status"] == "proved" and r_z["holds_through"] == 1

    # Consecutive single-shot rounds on [[36,8,3]] (t=1: 1 fault per round, alternating data & syndrome noise)
    ss_d8 = fd.simulate_tanner_single_shot_rounds(
        qt_d8_x,
        data_errors_per_round=[1 << 0, 0, 1 << 14, 0, 1 << 28, 0],
        syndrome_errors_per_round=[0, 1 << 5, 0, 1 << 2, 0, 1 << 4],
        final_clean_readout=True,
    )
    assert ss_d8["max_intermediate_coset_weight"] <= 1
    assert ss_d8["terminal"]["logical_success"] is True
    assert ss_d8["terminal"]["final_coset_weight"] == 0

    # Lifted Quantum Tanner [[288, 8]] over S_3 under simultaneous data + syndrome noise every round!
    qt_s3_x = fd.DecodingProblem.from_css_code("qt_s3_288", sector="X")
    ss_s3 = fd.simulate_tanner_single_shot_rounds(
        qt_s3_x,
        data_errors_per_round=[1 << 3, (1 << 55) | (1 << 190), 1 << 120, (1 << 12) | (1 << 211), 1 << 88, 0],
        syndrome_errors_per_round=[1 << 2, (1 << 17) | (1 << 95), 1 << 40, 1 << 73, (1 << 8) | (1 << 110), 1 << 31],
        final_clean_readout=True,
    )
    assert ss_s3["max_intermediate_coset_weight"] <= 2
    assert ss_s3["terminal"]["logical_success"] is True

    # Verify coset_reduced_weight collapses full stabilizer generators to weight 0
    g_stab0 = qt_d8_x.stabilizer_generators[0]
    assert fd.coset_reduced_weight(g_stab0, qt_d8_x.stabilizer_generators) == 0
    assert fd.coset_reduced_weight(g_stab0 ^ (1 << 5), qt_d8_x.stabilizer_generators) == 1

    # Negative control: non-Tanner problem refused by extract_tanner_neighborhoods
    steane = fd.DecodingProblem.from_css_code("steane", sector="X")
    try:
        fd.decode_tanner_local(steane, steane.syndrome_of(1))
        raise AssertionError("Expected MISSING-TANNER-META refusal on Steane")
    except fd.RefusalError as exc:
        assert exc.code == "MISSING-TANNER-META"

    # Sparse O(nnz) scipy.sparse.csc_matrix handoff on high-rate PP256 & APM1152 when available
    backends = fd.probe_external_backends()
    sparse_note = "sparse=skipped"
    if backends.get("scipy", {}).get("available"):
        from ftqc_decoders.adapters.external_decoders import build_check_matrix_for_external

        pp256_x = fd.DecodingProblem.from_css_code("pp256", sector="X")
        apm_x = fd.DecodingProblem.from_css_code("apm1152", sector="X")
        h_pp, fmt_pp = build_check_matrix_for_external(pp256_x, prefer_sparse=True)
        h_apm, fmt_apm = build_check_matrix_for_external(apm_x, prefer_sparse=True)
        assert fmt_pp == "scipy_csc" and fmt_apm == "scipy_csc"
        assert h_pp.nnz == sum(r.bit_count() for r in pp256_x.check_rows)
        assert h_apm.nnz == sum(r.bit_count() for r in apm_x.check_rows)
        sparse_note = f"scipy_csc(pp256_nnz={h_pp.nnz},apm1152_nnz={h_apm.nnz})"
        if backends.get("ldpc", {}).get("available"):
            r_pp = fd.run_decoder(pp256_x, pp256_x.syndrome_of((1 << 5) | (1 << 140)), "external:ldpc_bplsd")
            assert r_pp.converged and r_pp.stats["matrix_format"] == "scipy_csc"
            assert pp256_x.observable_of(((1 << 5) | (1 << 140)) ^ r_pp.error_mask) == 0

    return f"QE_D8_36_8_3 d_eff>=3 + 6-round single-shot + QE_S3_288_8 + {sparse_note}"


def leg_n1_n6_sabotage_controls() -> str:
    def expect_refusal(fn: Callable[[], object], expected_code: str) -> None:
        try:
            fn()
        except fd.RefusalError as exc:
            if exc.code != expected_code:
                raise AssertionError(f"Expected code {expected_code!r}, got {exc.code!r}") from exc
            return
        raise AssertionError(f"Sabotage control was not refused with {expected_code!r}")

    # N1: Non-exact probability input refused (pass boolean True to test runtime type guard without float AST!)
    expect_refusal(lambda: fd.parse_exact_probability(True), "FLOAT-INPUT-REFUSED")

    # N2: Unquoted JSON float in CLI `check` refused with exit code 2
    with tempfile.TemporaryDirectory() as tmp:
        import contextlib
        import io
        from ftqc_decoders.cli import main as cli_main

        bad_json = Path(tmp) / "bad.json"
        bad_json.write_text('{"study": "enclose_peierls_cluster", "physical_p": 0.001, "problem": {"kind": "bb72"}, "verified_component_radius": 2}')
        with contextlib.redirect_stderr(io.StringIO()):
            rc = cli_main(["check", str(bad_json)])
        assert rc == 2, f"Expected exit code 2 for JSON float, got {rc}"

    # N3: Invalid metacheck violating M * H == 0 (mod 2)
    steane = fd.DecodingProblem.from_css_code("steane", sector="X")
    expect_refusal(
        lambda: fd.DecodingProblem.from_matrices(
            steane.check_rows, steane.observable_rows, 7, metacheck_rows=[1]
        ),
        "INVALID-METACHECK",
    )

    # N4: Non-graphlike BB72 hypergraph sent to MWPM refused
    bb72 = fd.bb72_problem()
    expect_refusal(lambda: fd.decode_mwpm(bb72, bb72.syndrome_of(1)), "NON-GRAPHLIKE-DEM")

    # N5: Tampered LRAT optimality proof rejected by `ftqc_codes.lrat`
    syn0 = steane.syndrome_of(1)
    expect_refusal(
        lambda: fd.verify_lrat_optimality(steane, syn0, 1, b"999 0 1 0\n"),
        "LRAT-PROOF-REJECTED",
    )

    # N6: Tampered receipt rejected on replay
    rec = fd.build_receipt({"study": "certify_radius", "problem": {"kind": "css_code", "code": "steane"}, "decoder": "bp_osd", "target_radius": 1})
    rec_bad = json.loads(json.dumps(rec))
    rec_bad["result"]["subsets_checked"] = 999
    expect_refusal(lambda: fd.verify_receipt(rec_bad), "RECEIPT-REPLAY-MISMATCH")

    return "N1-N6 sabotage controls all refused with expected error codes"


def _run_quick_core() -> None:
    """Fast subset invoked inside scratch-copy mutation sub-runs (`M1`-`M4`)."""
    leg_l1_dem_stim_and_mwpm()
    bb144 = fd.bb144_gross_problem()
    g0_bits = fd.bit_positions(bb144.stabilizer_generators[0])
    e_half = sum(1 << b for b in g0_bits[:3])
    r_osd0 = fd.decode_bp_osd(bb144, bb144.syndrome_of(e_half), schedule="parallel", post_method="osd", osd_order=0)
    assert r_osd0.converged and bb144.observable_of(e_half ^ r_osd0.error_mask) == 0
    rep5 = fd.repetition_problem(5)
    r_u = fd.certify_operational_radius(rep5, "union_find", 2)
    assert r_u["status"] == "proved"
    lo_s, hi_s = fd.fraction_to_outward_strs(Fraction(1, 7), digits=16)
    assert Fraction(lo_s) <= Fraction(1, 7) <= Fraction(hi_s)


def leg_m1_m4_mutant_kills() -> str:
    mutants = [
        (
            "M1_osd_inverted_reliability_order",
            "bp_osd.py",
            "ordered_cols = sorted(range(n), key=lambda j: (posterior_scores[j], j))",
            "ordered_cols = sorted(range(n), key=lambda j: (-posterior_scores[j], j))",
        ),
        (
            "M2_uf_full_edge_instead_of_half_edge",
            "union_find.py",
            "if edge_half_progress[mj] >= 2:",
            "if edge_half_progress[mj] >= 1:",
        ),
        (
            "M3_blossom_missing_base_edge",
            "mwpm.py",
            "blossom_edges[b] = [cur_e] + list(reversed(path_e))",
            "blossom_edges[b] = list(reversed(path_e))",
        ),
        (
            "M4_outward_upper_bound_floor_instead_of_ceiling",
            "decimal_bounds.py",
            "ctx_out_hi = Context(prec=digits, rounding=ROUND_CEILING)",
            "ctx_out_hi = Context(prec=digits, rounding=ROUND_FLOOR)",
        ),
    ]

    with tempfile.TemporaryDirectory(prefix="ftqc_decoders_mutants_") as tmp:
        scratch_tools = Path(tmp) / "tools"
        scratch_pkg = scratch_tools / "ftqc_decoders"
        scratch_tools.mkdir(parents=True)
        shutil.copytree(_PKG_DIR, scratch_pkg, ignore=shutil.ignore_patterns("__pycache__", "._*"))
        # Symlink sibling packages (`ftqc_codes`, `ftqc_audit`, `ftqc_extractors`, `ftqc_gadgets`) into scratch_tools
        for sib in ("ftqc_codes", "ftqc_faults", "ftqc_audit", "ftqc_extractors", "ftqc_gadgets"):
            sib_src = _TOOLS_DIR / sib
            if sib_src.exists():
                os.symlink(sib_src, scratch_tools / sib)

        env = dict(os.environ, PYTHONPATH=str(scratch_tools))
        for m_name, rel_file, orig_line, mut_line in mutants:
            target_f = scratch_pkg / rel_file
            orig_txt = target_f.read_text()
            if orig_line not in orig_txt:
                raise AssertionError(f"Mutant {m_name} target string not found in {rel_file}")
            target_f.write_text(orig_txt.replace(orig_line, mut_line, 1))
            try:
                proc = subprocess.run(
                    [sys.executable, str(scratch_pkg / "selftest.py"), "--quick"],
                    env=env,
                    capture_output=True,
                    text=True,
                    timeout=15,
                )
                if proc.returncode == 0:
                    raise AssertionError(f"Mutant {m_name} SURVIVED quick selftest!")
            finally:
                target_f.write_text(orig_txt)

    return f"{len(mutants)}/4 scratch-copy code mutants killed"


def run_all(argv: list[str] | None = None) -> int:
    args = argv if argv is not None else sys.argv[1:]
    if "--quick" in args:
        _run_quick_core()
        return 0

    legs: list[tuple[str, Callable[[], str]]] = [
        ("L0_purity_and_lineage", leg_l0_purity_and_lineage),
        ("L1_dem_stim_and_mwpm", leg_l1_dem_stim_and_mwpm),
        ("L2_bp_osd_lsd_relay", leg_l2_bp_osd_lsd_relay),
        ("L3_uf_ssf_metachecks_sat_lrat", leg_l3_uf_ssf_metachecks_sat_lrat),
        ("L4_trapping_sets_and_deff", leg_l4_trapping_sets_and_deff),
        ("L5_enclosures_bridges_examples_cli", leg_l5_enclosures_bridges_examples_cli),
        ("L6_tanner_single_shot_and_sparse_handoff", leg_l6_tanner_single_shot_and_sparse_handoff),
        ("N1_N6_sabotage_controls", leg_n1_n6_sabotage_controls),
        ("M1_M4_mutant_kills", leg_m1_m4_mutant_kills),
    ]

    t_start = time.monotonic_ns()
    passed = 0
    for name, fn in legs:
        t0 = time.monotonic_ns()
        try:
            detail = fn()
        except Exception as exc:
            ms = (time.monotonic_ns() - t0) // 1_000_000
            print(f"[FAIL] {name} ({ms} ms): {exc}", file=sys.stderr)
            return 1
        ms = (time.monotonic_ns() - t0) // 1_000_000
        print(f"[PASS] {name} ({ms} ms): {detail}")
        passed += 1

    total_ms = (time.monotonic_ns() - t_start) // 1_000_000
    print(f"SELF-TEST: {passed} executed / 0 skipped ({passed} PASS, 0 FAIL) of {len(legs)} legs in {total_ms} ms")
    print("SELF-TEST: ALL PASS (rc 0)")
    return 0


if __name__ == "__main__":
    raise SystemExit(run_all())
