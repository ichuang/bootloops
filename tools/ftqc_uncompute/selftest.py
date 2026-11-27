"""Cold-clone acceptance battery for `ftqc_uncompute` (`< 10 s`, 100% Python stdlib).

Runs:
  - Positive verification legs `L0`–`L6` (AST zero-float & disjoint-lineage purity,
    `.kmx`/`GF(2^m)`/QASM3, HMR + Fiat-Shamir + Kraus dual-lineage, KMRV26/Gidney/A*
    spooky pebbling, superposition masking & Fibonacci ladders, `R1`–`R7` advisor &
    `--apply` equivalence, receipts & live sibling FTQC bridges).
  - Sabotage refusal controls `N1`–`N6`.
  - Scratch-copy code-mutation self-audits `M1`–`M4` (mutating isolated temp copies
    of `hmr.py`, `kickmix_sim.py`, `pebbling.py`, and `masking.py` in subprocesses and
    confirming the battery catches every single mutant).

100% Python standard library, MIT-licensed, zero float literals or float() calls.
"""

from __future__ import annotations

import ast
from decimal import Decimal
from fractions import Fraction
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Callable

# Support both `python3 tools/ftqc_uncompute/selftest.py` (from repo root) and `-m tools.ftqc_uncompute selftest`
_PKG_DIR = Path(__file__).resolve().parent
_TOOLS_DIR = _PKG_DIR.parent
_REPO_DIR = _TOOLS_DIR.parent
for _p in (str(_REPO_DIR), str(_TOOLS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from tools.ftqc_uncompute import (  # noqa: E402
    ArithDAG,
    CircuitBuilder,
    DAGNode,
    GF2Field,
    KickMixCircuit,
    LadderStep,
    PebbleDAG,
    PebbleStep,
    ReceiptRefusalError,
    RegDecl,
    adapters,
    apply_dag_rewrites,
    apply_kmx_rewrites,
    build_receipt_for_kmx,
    carry_runway_mask_bound,
    clopper_pearson_upper_bound,
    compile_gidney_sequential_schedule,
    compile_kmrv_blast_schedule,
    exact_ensemble_trace_distance,
    kmrv_theoretical_pebble_bound,
    main as cli_main,
    parse_strict_json,
    simulate_trajectories,
    solve_optimal_pebbling_astar,
    suggest_for_dag,
    uniform_residue_mask_bound,
    validate_pebbling_schedule,
    verify_channel_kraus,
    verify_hmr_circuit,
    verify_mask_spec,
    verify_receipt_file,
    verify_streaming_fibonacci_ladder,
)
from tools.ftqc_uncompute.kraus_oracle import CTX, ONE_D, ZERO_D  # noqa: E402


def _leg_l0_purity() -> str:
    py_files = sorted(_PKG_DIR.rglob("*.py"))
    assert len(py_files) >= 15
    for py_path in py_files:
        src = py_path.read_text(encoding="utf-8")
        tree = ast.parse(src, filename=str(py_path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, (float, complex)):
                raise AssertionError(f"Forbidden float/complex literal in {py_path.name}:{node.lineno}")
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in ("float", "complex"):
                raise AssertionError(f"Forbidden {node.func.id}() call in {py_path.name}:{node.lineno}")

    # Enforce disjoint lineage: kraus_oracle.py must never import hmr, kickmix_sim, or gf2
    ko_tree = ast.parse((_PKG_DIR / "kraus_oracle.py").read_text(encoding="utf-8"))
    for node in ast.walk(ko_tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            for banned in ("hmr", "kickmix_sim", "gf2"):
                if banned in node.module:
                    raise AssertionError(f"Disjoint-lineage violation in kraus_oracle.py: imports {node.module}")
    return f"{len(py_files)} modules zero-float & disjoint-lineage clean"


def _leg_l1_ir_and_algebra() -> str:
    for m in range(1, 7):
        fld = GF2Field(m)
        for a in range(1, min(fld.order, 12)):
            assert fld.mul(a, fld.invert(a)) == 1

    ex_adder = KickMixCircuit.from_kmx((_PKG_DIR / "examples" / "gidney_adder_4bit.kmx").read_text(encoding="utf-8"))
    assert ex_adder.operation_statistics().num_toffoli_like == 3
    assert ex_adder.reaction_depth() == 3
    assert KickMixCircuit.from_kmx("# comment\n" + ex_adder.to_kmx()).sha256() == ex_adder.sha256()

    qasm = adapters.kmx_to_qasm3(ex_adder)
    rt = adapters.kmx_from_qasm3(qasm)
    assert rt.operation_statistics().num_hmr == 3 and verify_hmr_circuit(rt).valid
    return "GF(2^m), .kmx SHA-256, and OpenQASM 3 round-trip verified"


def _leg_l2_hmr_sim_kraus() -> str:
    # 1. Gidney 4-bit adder dual-lineage check (HMR + 256-shot Fiat-Shamir + 11-qubit Kraus!)
    ex_add = KickMixCircuit.from_kmx((_PKG_DIR / "examples" / "gidney_adder_4bit.kmx").read_text(encoding="utf-8"))
    assert verify_hmr_circuit(ex_add).valid
    oracle_add = lambda d: ({"a": d["a"], "b": (d["a"] + d["b"]) & 0xF}, Fraction(0, 1))
    assert simulate_trajectories(ex_add, shots=256, oracle=oracle_add).passed
    kr_add = verify_channel_kraus(ex_add, oracle=oracle_add)
    assert kr_add.exact_match and kr_add.trace_distance == Decimal(0)

    # 2. DQI GF(2^3) MBUC fixture
    fld3 = GF2Field(3)
    ex_dqi = KickMixCircuit.from_kmx((_PKG_DIR / "examples" / "dqi_rs_tiny_mbuc.kmx").read_text(encoding="utf-8"))
    assert verify_hmr_circuit(ex_dqi).valid
    oracle_dqi = lambda d: (
        d,
        Fraction(fld3.mul(d["x"], d["y"]) & 1, 1) + Fraction((fld3.mul(d["x"], d["y"]) >> 1) & 1, 2),
    )
    assert simulate_trajectories(ex_dqi, shots=256, oracle=oracle_dqi).passed
    kr_dqi = verify_channel_kraus(ex_dqi, oracle=oracle_dqi)
    assert kr_dqi.exact_match and kr_dqi.trace_distance == Decimal(0)

    # 3. QROM del_lookup fixture
    tbl = [5, 3, 7, 2, 6]
    ex_lut = KickMixCircuit.from_kmx((_PKG_DIR / "examples" / "qrom_stream_deltas.kmx").read_text(encoding="utf-8"))
    assert verify_hmr_circuit(ex_lut).valid
    oracle_lut = lambda d: (
        d,
        Fraction((tbl[d["addr"]] if d["addr"] < 5 else 0) & 1, 1)
        + Fraction(1 if (((tbl[d["addr"]] if d["addr"] < 5 else 0) >> 1) & 1) else 0, 4),
    )
    assert simulate_trajectories(ex_lut, shots=256, oracle=oracle_lut).passed
    kr_lut = verify_channel_kraus(ex_lut, oracle=oracle_lut)
    assert kr_lut.exact_match and kr_lut.trace_distance == Decimal(0)

    eps0 = clopper_pearson_upper_bound(0, 4096, delta_neg_log2=128)
    assert Decimal("0.02") < eps0 < Decimal("0.023")
    return "HMR + Fiat-Shamir + Kraus exact (D_tr == 0) on iadd, del_gf2_mul, del_lookup"


def _leg_l3_spooky_pebbling() -> str:
    for ell in range(1, 31):
        sched = compile_kmrv_blast_schedule(ell)
        val = validate_pebbling_schedule(PebbleDAG.line(ell), sched)
        assert val.valid and val.depth == 2 * ell
        assert val.peak_pebbles <= kmrv_theoretical_pebble_bound(ell)

    g_val = validate_pebbling_schedule(PebbleDAG.line(16), compile_gidney_sequential_schedule(16))
    assert g_val.valid and g_val.is_sequential and g_val.peak_pebbles <= 6

    astar_s = solve_optimal_pebbling_astar(7, 4)
    assert astar_s is not None and validate_pebbling_schedule(PebbleDAG.line(7), astar_s, max_pebbles=4).valid
    win_s = solve_optimal_pebbling_astar(7, 5, window=2, window_ancilla_weight=1)
    assert win_s is not None and validate_pebbling_schedule(PebbleDAG.line(7, window=2), win_s, max_weighted_space=5).valid
    return "KMRV26 depth-2*ell (ell=1..30), Gidney sequential, and A* windowed schedules certified"


def _leg_l4_masking_and_ladders() -> str:
    dim16 = 16
    vec_m15 = [(ZERO_D, ZERO_D)] * dim16
    vec_r15 = [(ZERO_D, ZERO_D)] * dim16
    inv8 = CTX.divide(ONE_D, CTX.sqrt(Decimal(8)))
    inv15 = CTX.divide(ONE_D, CTX.sqrt(Decimal(15)))
    for r in range(15):
        vec_r15[r] = (inv15, ZERO_D)
        if math.gcd(r, 15) == 1:
            vec_m15[r] = (inv8, ZERO_D)
    dtr_exact = exact_ensemble_trace_distance([(ONE_D, vec_r15)], [(ONE_D, vec_m15)])
    dtr_bound = Decimal(uniform_residue_mask_bound(15)["single_mask_trace_distance_upper"])
    assert dtr_exact <= dtr_bound and (dtr_bound - dtr_exact) < Decimal("1e-20")

    rw_b = Decimal(carry_runway_mask_bound(4, 1)["single_step_trace_distance_upper"])
    assert Decimal("0.34") < rw_b < Decimal("0.36")

    rsa_ub = uniform_residue_mask_bound((1 << 2047) + 99, runs=45, min_prime_factor=(1 << 1023))
    assert Decimal(rsa_ub["total_tv_distance_upper"]) < Decimal("1e-150")

    fib_spec = parse_strict_json((_PKG_DIR / "examples" / "fib_masked_iqpe_n15.json").read_text(encoding="utf-8"))
    steps = [LadderStep.from_dict(s) for s in fib_spec["steps"]]
    l_res = verify_streaming_fibonacci_ladder(
        15, 2, steps,
        expected_control_weights={"c0": 1, "c1": 1, "c2": 2, "c3": 3},
        subroutine_kind="euclidean_dialog",
        subroutine_steps=7,
    )
    assert l_res.valid and l_res.reverse_uncompute_eliminated
    return "Mask D_tr bounds match Kraus to 20+ digits; Fibonacci ladder certified"


def _leg_l5_advisor_and_rewrites() -> str:
    cb = CircuitBuilder()
    qx = cb.create_quantum_register(3, name="x")
    s0, s1 = cb.alloc_qubits(2)
    cb.ccx(qx[0], qx[1], s0)
    cb.z_pow(s0, Fraction(1, 4))
    cb.ccx(qx[0], qx[1], s0)
    cb.ccx(qx[1], qx[2], s1)
    cb.z(s1)
    cb.ccx(qx[1], qx[2], s1)
    orig = cb.finish_circuit()
    rew, sugs = apply_kmx_rewrites(orig, verify=True)
    assert rew.operation_statistics().num_toffoli_like == 2 and rew.num_qubits == 4
    assert {"R1_TEMP_AND_HMR", "R5_STREAM_MBUC"} <= {s.rule_id for s in sugs}

    dag = ArithDAG(
        name="r1_to_r7",
        registers=(
            RegDecl("a", "quint", 4, role="io"),
            RegDecl("b", "quint", 4, role="io"),
            RegDecl("gx", "qgf2", 3, role="input"),
            RegDecl("gy", "qgf2", 3, role="input"),
            RegDecl("and_t", "quint", 1, role="scratch"),
            RegDecl("gf1", "qgf2", 3, role="scratch"),
            RegDecl("gf2", "qgf2", 3, role="scratch"),
            RegDecl("lut_t", "quint", 2, role="scratch"),
            RegDecl("m1", "qmodint", 4, role="scratch", modulus=15),
            RegDecl("m2", "qmodint", 4, role="scratch", modulus=15),
            RegDecl("m3", "qmodint", 4, role="scratch", modulus=15),
            RegDecl("m4", "qmodint", 4, role="scratch", modulus=15),
        ),
        nodes=(
            DAGNode("ca", "and_gate", ("a", "b"), ("and_t",), chirality="compute"),
            DAGNode("ua", "and_gate", ("a", "b"), ("and_t",), chirality="uncompute", uncompute_mode="reverse", pair_with="ca"),
            DAGNode("cg1", "gf2_mul", ("gx", "gy"), ("gf1",), chirality="compute"),
            DAGNode("ug1", "gf2_mul", ("gx", "gy"), ("gf1",), chirality="uncompute", uncompute_mode="reverse", pair_with="cg1"),
            DAGNode("cg2", "gf2_mul", ("gx", "gy"), ("gf2",), chirality="compute"),
            DAGNode("ug2", "gf2_mul", ("gx", "gy"), ("gf2",), chirality="uncompute", uncompute_mode="hmr", pair_with="cg2"),
            DAGNode("cl", "qrom_lookup", ("gx",), ("lut_t",), chirality="compute", params={"table": [1, 2, 3, 1, 2, 3, 0, 1]}),
            DAGNode("ul", "qrom_lookup", ("gx",), ("lut_t",), chirality="uncompute", uncompute_mode="reverse", pair_with="cl", params={"table": [1, 2, 3, 1, 2, 3, 0, 1]}),
            DAGNode("add0", "int_add", ("a",), ("b",), chirality="inplace", params={"use_hmr": False}),
            DAGNode("ph0", "phase_by_bits", ("a",), (), chirality="phase", uncompute_mode="none"),
            DAGNode("s1", "mod_sqr", ("a",), ("m1",), chirality="compute"),
            DAGNode("s2", "mod_sqr", ("m1",), ("m2",), chirality="compute"),
            DAGNode("s3", "mod_sqr", ("m2",), ("m3",), chirality="compute"),
            DAGNode("s4", "mod_sqr", ("m3",), ("m4",), chirality="compute"),
            DAGNode("u4", "mod_sqr", ("m3",), ("m4",), chirality="uncompute", uncompute_mode="reverse", pair_with="s4"),
            DAGNode("u3", "mod_sqr", ("m2",), ("m3",), chirality="uncompute", uncompute_mode="reverse", pair_with="s3"),
            DAGNode("u2", "mod_sqr", ("m1",), ("m2",), chirality="uncompute", uncompute_mode="reverse", pair_with="s2"),
            DAGNode("u1", "mod_sqr", ("a",), ("m1",), chirality="uncompute", uncompute_mode="reverse", pair_with="s1"),
        ),
    )
    found = {s.rule_id for s in suggest_for_dag(dag)}
    assert len(found) == 7, f"Expected all 7 rules R1-R7, got {found}"
    return "All 7 rules (R1-R7) detected and `--apply` equivalence verified"


def _leg_l6_receipts_and_bridges() -> str:
    ex_files = [
        _PKG_DIR / "examples" / "gidney_adder_4bit.kmx",
        _PKG_DIR / "examples" / "dqi_rs_tiny_mbuc.kmx",
        _PKG_DIR / "examples" / "qrom_stream_deltas.kmx",
        _PKG_DIR / "examples" / "regev_spooky_chain.json",
        _PKG_DIR / "examples" / "fib_masked_iqpe_n15.json",
    ]
    with tempfile.TemporaryDirectory() as tmpd:
        td = Path(tmpd)
        for idx, ef in enumerate(ex_files):
            rf = td / f"rec_{idx}.json"
            assert cli_main(["check", str(ef), "--shots", "128", "-o", str(rf)]) == 0
            ok, errs = verify_receipt_file(rf, target_path=ef)
            assert ok and not errs, f"Receipt verify failed for {ef.name}: {errs}"

    bridge_notes: list[str] = ["5 example receipts verified"]
    if adapters.ftqc_audit_available():
        hmr_c = KickMixCircuit.from_kmx(
            'REGISTER r0 "x"\nAPPEND_TO_REGISTER q0 r0\nCX q0 q1\nHMR q1 b0\nZ q0 if b0\nBIT_STORE0 b0\n'
        )
        rep = adapters.audit_clifford_hmr_with_ftqc_audit(hmr_c)
        assert rep["verified"]
        bridge_notes.append("ftqc_audit live")
    if adapters.ftqc_zx_available():
        try:
            import ftqc_zx  # type: ignore[import-not-found]

            cat4, _ = ftqc_zx.synth_cat_state(4)
            zxc = ftqc_zx.build_gadget_contract(cat4, p_phys="1/1000")
            ex_dqi = KickMixCircuit.from_kmx((_PKG_DIR / "examples" / "dqi_rs_tiny_mbuc.kmx").read_text(encoding="utf-8"))
            priced = adapters.price_circuit_with_zx_contract(ex_dqi, zxc, t_factory_escape_hi="21/2000")
            assert Decimal(priced["total_wire_ticks_upper"]) > Decimal(0)
            bridge_notes.append("ftqc_zx live")
        except ModuleNotFoundError:
            bridge_notes.append("ftqc_zx partial-local")
    if adapters.ftqc_compose_available():
        import ftqc_compose  # type: ignore[import-not-found]

        _ = ftqc_compose
        bridge_notes.append("ftqc_compose live")
    return ", ".join(bridge_notes)


def _leg_sabotage_controls() -> str:
    # N1: Missing CZ fixup in del_and
    c_n1 = KickMixCircuit.from_kmx(
        'REGISTER r0 "x"\nAPPEND_TO_REGISTER q0 r0\nAPPEND_TO_REGISTER q1 r0\nCCX q0 q1 q2\nHMR q2 b0\nBIT_STORE0 b0\n'
    )
    r_n1 = verify_hmr_circuit(c_n1)
    assert not r_n1.valid and any(f.code == "UNCOMPENSATED-HMR-PHASE" for f in r_n1.findings)
    assert simulate_trajectories(c_n1, shots=128).failures > 0
    assert not verify_channel_kraus(c_n1).exact_match

    # N2: Dirty Reset `R` on active qubit
    c_n2 = KickMixCircuit.from_kmx(
        'REGISTER r0 "x"\nAPPEND_TO_REGISTER q0 r0\nAPPEND_TO_REGISTER q1 r0\nCCX q0 q1 q2\nR q2\n'
    )
    r_n2 = verify_hmr_circuit(c_n2)
    assert not r_n2.valid and any(f.code == "DIRTY-RESET-HAZARD" for f in r_n2.findings)
    s_n2 = simulate_trajectories(c_n2, shots=128)
    assert s_n2.dirty_resets > 0 and not s_n2.coherent_bound_valid

    # N3: Non-injective branch (`HMR q0 b0` erasing input)
    c_n3 = KickMixCircuit.from_kmx('REGISTER r0 "x"\nAPPEND_TO_REGISTER q0 r0\nHMR q0 b0\nBIT_STORE0 b0\n')
    r_n3 = verify_hmr_circuit(c_n3)
    assert not r_n3.valid and any(f.code == "NON-INJECTIVE-KICKMIX-BRANCH" for f in r_n3.findings)

    # N4: Uncleaned ghost at exit & parallel active-set conflict
    v_n4a = validate_pebbling_schedule(
        PebbleDAG.line(2),
        [PebbleStep(pebble=(1,)), PebbleStep(pebble=(2,), ghost=(1,)), PebbleStep(pebble=(1,)), PebbleStep(unpebble=(2,), ghost=(1,))],
    )
    assert not v_n4a.valid and any(f.code == "UNCLEANED-GHOST-AT-EXIT" for f in v_n4a.findings)
    v_n4b = validate_pebbling_schedule(
        PebbleDAG.line(2),
        [PebbleStep(pebble=(1,)), PebbleStep(pebble=(2,), unpebble=(1,))],
    )
    assert not v_n4b.valid and any(f.code == "PARALLEL-ACTIVE-SET-CONFLICT" for f in v_n4b.findings)

    # N5: Non-coprime mask multiplier + truncated Euclidean Dialog overflow
    m_n5a = verify_mask_spec("mult_mod_n", modulus=15, multipliers=[3])
    assert not m_n5a.valid and any(f.code == "MASK-NON-BIJECTIVE-GARBAGE" for f in m_n5a.findings)
    m_n5b = verify_mask_spec(
        "mult_mod_n", modulus=(1 << 2048) - 159, multipliers=[2], subroutine_kind="euclidean_dialog", subroutine_steps=2900
    )
    assert not m_n5b.valid and any(f.code == "MASK-SUPPORT-OVERFLOW-HAZARD" for f in m_n5b.findings)

    # N6: Float literal injected into certified JSON
    try:
        parse_strict_json('{"schema": "ftqc_uncompute-receipt-1", "eps": 0.01}')
        raise AssertionError("N6 failed to reject float literal!")
    except ReceiptRefusalError as exc:
        assert exc.code == "FLOAT-LITERAL-IN-RECEIPT"

    return "N1-N6 sabotage controls all refused with expected error codes"


def _leg_mutants() -> str:
    """M1-M4: Mutate isolated temporary scratch copies of the verifier modules and assert failure!"""
    mutants = [
        (
            "M1_hmr_skip_kickback",
            "hmr.py",
            "_phase_add_inplace(phase_poly, _anf_to_phase_poly(kick, Fraction(1, 1)))",
            "pass  # MUTANT M1",
        ),
        (
            "M2_sim_ignore_dirty_reset",
            "kickmix_sim.py",
            "dirty_resets += dirty.bit_count()",
            "pass  # MUTANT M2",
        ),
        (
            "M3_pebble_ignore_ghost_exit",
            "pebbling.py",
            'if state[v] == NODE_GHOST:',
            'if False and state[v] == NODE_GHOST:',
        ),
        (
            "M4_mask_ignore_gcd",
            "masking.py",
            "if g_mod == 0 or g_gcd != 1:",
            "if False:",
        ),
    ]

    with tempfile.TemporaryDirectory() as tmpd:
        base = Path(tmpd) / "tools"
        base.mkdir(parents=True)
        scratch_pkg = base / "ftqc_uncompute"
        shutil.copytree(_PKG_DIR, scratch_pkg, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))

        probe_code = (
            "import sys; sys.path.insert(0, sys.argv[1]); "
            "from tools.ftqc_uncompute.selftest import _leg_sabotage_controls; "
            "_leg_sabotage_controls()"
        )

        for m_name, rel_file, target_str, repl_str in mutants:
            target_file = scratch_pkg / rel_file
            orig_src = target_file.read_text(encoding="utf-8")
            if target_str not in orig_src:
                raise AssertionError(f"Mutant {m_name} target string missing in {rel_file}")
            try:
                target_file.write_text(orig_src.replace(target_str, repl_str, 1), encoding="utf-8")
                proc = subprocess.run(
                    [sys.executable, "-c", probe_code, str(Path(tmpd))],
                    capture_output=True,
                    text=True,
                )
                if proc.returncode == 0:
                    raise AssertionError(f"SURVIVING MUTANT {m_name}: sabotage gate did not catch mutation!")
            finally:
                target_file.write_text(orig_src, encoding="utf-8")

    return "M1-M4 scratch-copy code mutants all killed"


def _regress_global_budget() -> str:
    n = 14
    header = ['REGISTER r0 "x"'] + [f"APPEND_TO_REGISTER q{i} r0" for i in range(n)]
    invert = [f"X q{i}" for i in range(n)]
    chain = [f"CCX q0 q1 q{n}"] + [
        f"CCX q{n+j-2} q{j} q{n+j-1}" for j in range(2, n)
    ]
    bad = [f"CCX q0 q1 q{n}", f"HMR q{n} b0"]
    bare = KickMixCircuit.from_kmx("\n".join(header + bad))
    assert any(f.code == "UNCOMPENSATED-HMR-PHASE" for f in verify_hmr_circuit(bare).findings)
    padded = KickMixCircuit.from_kmx("\n".join(header + invert + chain + chain[::-1] + invert + bad))
    receipt = build_receipt_for_kmx(padded, shots=8)
    assert receipt["verdict"] == "REFUSED"
    assert receipt["hmr_verification"]["findings"][0]["code"] == "SYMBOLIC-BUDGET-EXCEEDED"
    assert not receipt["fiat_shamir_trajectories"]["coherent_bound_valid"]
    return "identity padding cannot certify a missing phase fixup"


def _regress_receipt_replay() -> str:
    import json

    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp)
        target = directory / "relocated.kmx"
        target.write_text("X q0\nX q0\n")
        path = directory / "receipt.json"
        path.write_text(json.dumps({"schema": "ftqc_uncompute-receipt-1", "verdict": "CERTIFIED"}))
        assert not verify_receipt_file(path)[0]
        circuit = KickMixCircuit.from_kmx(target.read_text())
        receipt = build_receipt_for_kmx(circuit, source_name="<memory>", shots=7,
                                        salt="custom", run_kraus=False)
        path.write_text(json.dumps(receipt))
        assert not verify_receipt_file(path)[0]
        assert verify_receipt_file(path, target_path=target) == (True, [])
        for section, key, value in (("hmr_verification", "phase_cancelled", False),
                                    ("summary", "num_qubits", True),
                                    ("fiat_shamir_trajectories", "epsilon_upper", "0")):
            edited = dict(receipt)
            edited[section] = {**receipt[section], key: value}
            path.write_text(json.dumps(edited))
            assert not verify_receipt_file(path, target_path=target)[0]
        legacy = build_receipt_for_kmx(circuit, shots=9)
        del legacy["verification_options"]
        path.write_text(json.dumps(legacy))
        assert verify_receipt_file(path, target_path=target) == (True, [])
        for malformed in ('{"x":1,"x":2}', '{"x":Infinity}', '{"x":NaN}'):
            try:
                parse_strict_json(malformed)
            except ReceiptRefusalError:
                continue
            raise AssertionError(malformed)
    return "source binding, replay options, strict JSON and evidence checked"


def run_selftest() -> int:
    t0 = time.perf_counter_ns()
    legs: list[tuple[str, Callable[[], str]]] = [
        ("_regress_global_budget", _regress_global_budget),
        ("_regress_receipt_replay", _regress_receipt_replay),
        ("L0_purity_and_lineage", _leg_l0_purity),
        ("L1_ir_gf2_qasm3", _leg_l1_ir_and_algebra),
        ("L2_hmr_sim_kraus_dual", _leg_l2_hmr_sim_kraus),
        ("L3_spooky_pebbling", _leg_l3_spooky_pebbling),
        ("L4_superposition_masking", _leg_l4_masking_and_ladders),
        ("L5_advisor_r1_r7", _leg_l5_advisor_and_rewrites),
        ("L6_receipts_and_bridges", _leg_l6_receipts_and_bridges),
        ("N1_N6_sabotage_controls", _leg_sabotage_controls),
        ("M1_M4_mutant_kills", _leg_mutants),
    ]

    passed = 0
    failed = 0
    for name, fn in legs:
        t_leg = time.perf_counter_ns()
        try:
            detail = fn()
            dt_ms = (time.perf_counter_ns() - t_leg) // 1_000_000
            print(f"[PASS] {name} ({dt_ms} ms): {detail}")
            passed += 1
        except Exception as exc:
            dt_ms = (time.perf_counter_ns() - t_leg) // 1_000_000
            print(f"[FAIL] {name} ({dt_ms} ms): {exc}", file=sys.stderr)
            failed += 1

    total_ms = (time.perf_counter_ns() - t0) // 1_000_000
    print(f"SELF-TEST: {passed} executed / 0 skipped ({passed} PASS, {failed} FAIL) of {len(legs)} legs in {total_ms} ms")
    if failed == 0:
        print("SELF-TEST: ALL PASS (rc 0)")
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(run_selftest())
