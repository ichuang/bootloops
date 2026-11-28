"""Addressable high-rate QEC surgery adapter, translated-seed substitution, and scale-up studies.

Connects verified `ftqc_codes` receipts (`ftqc_codes-receipt-1` / `ftqc_codes-construction-1`),
paired logical basis receipts (`ftqc_codes-basis-1`), basis conversions, and cyclic symmetry
certificates to `ftqc_gadgets.surgery`, `ftqc_extractors.installed_surgery`, `ftqc_faults`,
`ftqc_decoders`, `ftqc_compose`, and `ftqc_uncompute`.
"""

import copy
import hashlib
import json
from types import SimpleNamespace

from ftqc_codes import (
    basis as fc_basis,
    certify as fc_certify,
    codes as fc_codes,
    constructions as fc_constructions,
    families_highrate as fc_highrate,
    gf2 as fc_gf2,
    symmetry as fc_symmetry,
)
from ftqc_gadgets.surgery import Group, Instrument, commute, paired, pauli, receipt as surgery_receipt, verify as verify_surgery
from ftqc_gadgets.surgery_build import graph as build_gauging_graph
from ftqc_gadgets.surgery_design import build as build_surgery_design
from ftqc_gadgets.surgery_geometry import restricted_cone
from ftqc_gadgets.surgery_merge import word as css_word
from ftqc_faults.audit import build as build_audit_receipt, verify as verify_audit_receipt
from ftqc_faults.audit.circuit import Circuit as AuditCircuit
from ftqc_faults.audit.stab import Gen, Sim
from .installed_surgery import analyze as analyze_installed, compile_schedule
from .surgery_layout import translated_settings
from .surgery_products import lifted_product, radial90, radial_design, radial_homotopy_design


def _canonical_sha256(obj):
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def code_checks_pauli(code):
    """Return CSS stabilizer checks as `+X...` and `+Z...` strings (qubit 0 first)."""
    return [fc_basis.pauli_from_xz(code.n, r, 0, 0) for r in code.hx] + [
        fc_basis.pauli_from_xz(code.n, 0, r, 0) for r in code.hz
    ]


def verify_bound_code_and_basis(
    code_or_receipt,
    basis_receipt,
    *,
    old_basis_receipt=None,
    transform_receipt=None,
):
    """Re-verify code + `ftqc_codes-basis-1` receipts and enforce digest bindings."""
    if isinstance(code_or_receipt, fc_codes.CSSCode):
        code = code_or_receipt
        if not code.commutes():
            raise ValueError("code checks do not commute")
        code_rec_sha = None
    elif isinstance(code_or_receipt, dict):
        schema = code_or_receipt.get("schema")
        if schema == fc_certify.SCHEMA or schema in fc_certify.LEGACY_SCHEMAS:
            fc_certify.verify(code_or_receipt)
            code = fc_codes.CSSCode.from_json(code_or_receipt["code"])
            code_rec_sha = _canonical_sha256(code_or_receipt)
        elif schema == fc_constructions.RECEIPT_SCHEMA:
            fc_constructions.verify_construction_receipt(code_or_receipt)
            if code_or_receipt["payload"]["ancestry"] == "reproduced":
                code = fc_constructions.from_spec(code_or_receipt["payload"]["replay_spec"])
            else:
                si = code_or_receipt["payload"]["seed_inputs"]
                n_m = code_or_receipt["metrics"]["n"]
                hx = [fc_gf2.from_supp(s) for s in si["hx_supports"]]
                hz = [fc_gf2.from_supp(s) for s in si["hz_supports"]]
                code = fc_codes.CSSCode(code_or_receipt["payload"]["name"], n_m, hx, hz)
            code_rec_sha = code_or_receipt["construction_sha256"]
        else:
            raise ValueError(f"unsupported code receipt schema {schema!r}")
    else:
        raise ValueError("expected CSSCode or verified ftqc_codes receipt")

    try:
        fc_basis.verify_basis_receipt(code, basis_receipt)
        if transform_receipt is not None:
            if old_basis_receipt is None:
                raise ValueError("transform_receipt requires old_basis_receipt")
            fc_basis.verify_basis_transform(code, old_basis_receipt, basis_receipt, transform_receipt)
    except (fc_certify.Refusal, fc_codes.SpecError) as exc:
        raise ValueError(f"basis binding check failed: {exc}") from exc

    return code, {
        "code_sha256": code.sha256(),
        "code_receipt_sha256": code_rec_sha,
        "basis_sha256": basis_receipt["basis_sha256"],
        "transform_sha256": _canonical_sha256(transform_receipt) if transform_receipt else None,
    }


def build_slot_surgery_instrument(
    code_or_receipt,
    basis_receipt,
    slot_name,
    axis,
    *,
    sign=0,
    method="cone",
    rounds=2,
    edges=None,
    surviving_slots=None,
    audit_order=0,
    audit_cap=100_000,
):
    """Build a verified signed `ftqc_gadgets` surgery instrument for a named logical slot."""
    code, digests = verify_bound_code_and_basis(code_or_receipt, basis_receipt)
    try:
        target_word = fc_basis.embed_named_pauli(basis_receipt, [(slot_name, axis)], sign=sign)
        fc_basis.verify_logical_target(code, basis_receipt, [(slot_name, axis)], target_word, sign=sign)
        b_export = fc_basis.export_surgery_basis(
            basis_receipt,
            measured_slots=[slot_name],
            surviving_slots=surviving_slots,
        )
    except (fc_certify.Refusal, fc_codes.SpecError) as exc:
        raise ValueError(f"slot target verification failed: {exc}") from exc

    px, pz, _ = fc_basis.parse_pauli(target_word, code.n)
    supp_sites = fc_gf2.supp(px | pz)
    if edges is None:
        edges = [(i, i + 1) for i in range(len(supp_sites) - 1)]

    checks = code_checks_pauli(code)
    spec = build_gauging_graph(code.n, checks, target_word, edges, method=method, rounds=rounds)
    spec["basis"] = b_export["surgery_basis"]
    rec = surgery_receipt(spec, order=audit_order, cap=audit_cap)
    return spec, rec, {
        **digests,
        "slot": slot_name,
        "axis": axis,
        "sign": sign,
        "target_word": target_word,
        "support_sites": supp_sites,
        "surviving_slots": list(b_export["surgery_basis"]["surviving"]),
    }


def _derive_z_control_sites_for_slot(code, basis_receipt, target_slot, target_axis):
    """Find physical Z-operator support on `code` flipping `target_slot`'s `X`/`Y` sign cleanly.

    Because `installed_surgery.compile_schedule` conjugates the controlled gadget by
    `prod_{q in sites} Z_q` before and after execution, choosing `sites = supp(Z_{target_slot})`
    when `Z_{target_slot}` is pure Z (or solving in `ker(H_X)` for general dressed bases)
    guarantees commutation with all code checks and all other logical slots while flipping
    the sign of `X_{target_slot}` or `Y_{target_slot}`!
    """
    if target_axis not in ("X", "Y"):
        raise ValueError("Z_IF classical control flips X or Y logical reads; got Z")
    t_idx = fc_basis.slot_index(basis_receipt, target_slot)
    z_word = basis_receipt["slots"][t_idx]["z"]
    zx, zz, _ = fc_basis.parse_pauli(z_word, code.n)
    if zx == 0 and zz != 0:
        return fc_gf2.supp(zz)

    # General symplectic fallback when slot Z was dressed by X-stabilizers:
    # strip X-stabilizer dressing by solving in ker(H_X) for a pure-Z vector pairing only with X_t!
    z_cands = fc_gf2.nullspace(code.hx, code.n)
    k = basis_receipt["k"]
    slot_xs = [fc_basis.parse_pauli(basis_receipt["slots"][i]["x"], code.n)[0] for i in range(k)]
    # Build linear system over z_cands: dot(z_cand, slot_xs[i]) == (1 if i == t_idx else 0)
    rows = []
    for c_mask in z_cands:
        sig = sum(fc_gf2.dot(c_mask, slot_xs[i]) << i for i in range(k))
        rows.append((sig, c_mask))
    pivs = {}
    for sig, c_mask in rows:
        s_cur, m_cur = sig, c_mask
        for p in sorted(pivs, reverse=True):
            if (s_cur >> p) & 1:
                s_cur ^= pivs[p][0]
                m_cur ^= pivs[p][1]
        if s_cur:
            pivs[s_cur.bit_length() - 1] = (s_cur, m_cur)
    want_sig = 1 << t_idx
    sol_mask = 0
    for p in sorted(pivs, reverse=True):
        if (want_sig >> p) & 1:
            want_sig ^= pivs[p][0]
            sol_mask ^= pivs[p][1]
    if want_sig != 0 or sol_mask == 0:
        raise ValueError(f"could not synthesize pure-Z control support for slot {target_slot!r}")
    return fc_gf2.supp(sol_mask)


def install_addressable_program(
    code_or_receipt,
    basis_receipt,
    steps,
    *,
    old_basis_receipt=None,
    transform_receipt=None,
    couplers=None,
    reuse_ancillas=True,
    idles=True,
    require_explicit_idles=True,
    terminal_syndrome=True,
    target_distance=2,
    budget=200_000,
):
    """Compile and verify a multi-step addressable surgery program on an explicit coupler graph.

    Each entry in `steps` is a dict `{"slot": str, "axis": "X"|"Y"|"Z", "sign": 0|1,
    "rounds": int, "method": "cone"|"homomorphic", "control_from_step": int | None,
    "control_sites": list[int] | None}`.
    """
    code, digests = verify_bound_code_and_basis(
        code_or_receipt,
        basis_receipt,
        old_basis_receipt=old_basis_receipt,
        transform_receipt=transform_receipt,
    )
    if require_explicit_idles and not idles:
        raise ValueError("boundary/schedule policy requires explicit idle locations (idles=True)")

    all_slot_names = [s["name"] for s in basis_receipt["slots"]]
    surviving = list(all_slot_names)
    measured_so_far = []
    batches = []
    expected_targets = []
    step_ledger = []
    anc_offset = code.n
    auto_couplers = set()

    for step_idx, st in enumerate(steps):
        slot_nm = st["slot"]
        axis = st["axis"]
        sign = int(st.get("sign", 0))
        rounds = int(st.get("rounds", 2))
        method = st.get("method", "cone")
        if slot_nm not in surviving:
            raise ValueError(f"step {step_idx}: slot {slot_nm!r} is not among surviving slots {surviving}")
        surviving.remove(slot_nm)
        measured_so_far.append(slot_nm)

        local_spec, local_rec, slot_meta = build_slot_surgery_instrument(
            code,
            basis_receipt,
            slot_nm,
            axis,
            sign=sign,
            method=method,
            rounds=rounds,
            surviving_slots=None,
            audit_order=0,
        )
        local_width = 1 + max(
            [code.n - 1]
            + [
                q
                for row in local_spec["ticks"]
                for op in row
                for q in (op[1:2] if op[0].lstrip("!") == "Z_IF" else op[1:])
            ]
        )
        n_anc_local = local_width - code.n
        start_anc = code.n if reuse_ancillas else anc_offset
        if not reuse_ancillas:
            anc_offset += n_anc_local
        placement = list(range(code.n)) + list(range(start_anc, start_anc + n_anc_local))

        ctrl_spec = None
        ctrl_sites = None
        if st.get("control_from_step") is not None:
            c_idx = int(st["control_from_step"])
            ctrl_sites = (
                list(st["control_sites"])
                if st.get("control_sites") is not None
                else _derive_z_control_sites_for_slot(code, basis_receipt, slot_nm, axis)
            )
            ctrl_spec = [c_idx, ctrl_sites]

        setting = {"gadget": local_rec, "placement": placement}
        if ctrl_spec is not None:
            setting["control"] = ctrl_spec
        batches.append([setting])
        expected_targets.append(slot_meta["target_word"])

        for row in local_spec["ticks"]:
            for op in row:
                if op[0].lstrip("!") == "CX":
                    auto_couplers.add(tuple(sorted((placement[op[1]], placement[op[2]]))))

        step_ledger.append(
            {
                "step": step_idx,
                "slot": slot_nm,
                "axis": axis,
                "sign": sign,
                "target_word": slot_meta["target_word"],
                "support_sites": slot_meta["support_sites"],
                "surviving_after_step": list(surviving),
                "ancilla_wires": list(range(start_anc, start_anc + n_anc_local)),
                "control_from_step": st.get("control_from_step"),
                "control_z_sites": ctrl_sites,
            }
        )

    coupler_list = (
        [list(p) for p in couplers]
        if couplers is not None
        else [list(p) for p in sorted(auto_couplers)]
    )
    sched_input = {
        "schema": 1,
        "kind": "installed_surgery",
        "n": code.n,
        "checks": code_checks_pauli(code),
        "targets": expected_targets,
        "couplers": coupler_list,
        "batches": batches,
        "idles": bool(idles),
        "terminal_syndrome": bool(terminal_syndrome),
        "target": int(target_distance),
    }

    compiled_spec, peak_live = compile_schedule(sched_input)
    if require_explicit_idles:
        idle_ops = sum(1 for row in compiled_spec["ticks"] for op in row if op[0] == "I")
        if idle_ops == 0:
            raise ValueError("compiled schedule contains zero explicit idle operations")

    joint_basis_export = fc_basis.export_surgery_basis(
        basis_receipt,
        measured_slots=measured_so_far,
        surviving_slots=surviving,
    )
    compiled_spec["basis"] = joint_basis_export["surgery_basis"]
    joint_receipt = surgery_receipt(compiled_spec, order=max(0, target_distance - 1), cap=max(1, budget))

    # Extract split-frame Z_IF updates, measurement IDs, and per-wire ancilla lifetimes
    zif_updates = []
    anc_lifetimes = {}
    meas_counter = 0
    raw_step_records = []
    for tick_i, row in enumerate(compiled_spec["ticks"]):
        for op in row:
            gname = op[0].lstrip("!")
            if gname == "Z_IF":
                zif_updates.append({"tick": tick_i, "qubit": op[1], "record": op[2]})
            if gname in ("R", "RX") and op[1] >= code.n:
                anc_lifetimes.setdefault(str(op[1]), []).append({"alloc_tick": tick_i, "release_tick": None})
            if gname in ("M", "MZ", "MX"):
                if op[1] >= code.n and str(op[1]) in anc_lifetimes:
                    anc_lifetimes[str(op[1])][-1]["release_tick"] = tick_i
                meas_counter += 1

    for idx_s, entry in enumerate(step_ledger):
        entry["compiled_readout"] = copy.deepcopy(compiled_spec["readouts"][idx_s])
        if entry["control_from_step"] is not None:
            ctl_ro = compiled_spec["readouts"][entry["control_from_step"]]
            raw_recs = sorted(set(compiled_spec["readouts"][idx_s]["records"]) ^ set(ctl_ro["records"]))
            raw_sign = compiled_spec["readouts"][idx_s]["sign"] ^ ctl_ro["sign"]
        else:
            raw_recs = list(compiled_spec["readouts"][idx_s]["records"])
            raw_sign = compiled_spec["readouts"][idx_s]["sign"]
        entry["raw_physical_readout"] = {"records": raw_recs, "sign": raw_sign}
        raw_step_records.append(entry["raw_physical_readout"])

    return {
        "schema": "ftqc_extractors-addressable-surgery-1",
        "code_name": code.name,
        "n": code.n,
        "k": basis_receipt["k"],
        **digests,
        "reuse_ancillas": bool(reuse_ancillas),
        "installed_couplers": len({tuple(sorted(p)) for p in coupler_list}),
        "couplers": coupler_list,
        "peak_live_qubits": peak_live,
        "peak_ancilla_qubits": peak_live - code.n,
        "ticks": len(compiled_spec["ticks"]),
        "total_measurements": meas_counter,
        "surviving_slots": surviving,
        "step_ledger": step_ledger,
        "split_frame_z_if_updates": zif_updates,
        "ancilla_lifetimes": anc_lifetimes,
        "schedule_input": sched_input,
        "compiled_spec": compiled_spec,
        "surgery_receipt": joint_receipt,
    }


def _eval_linear_form(form, var_assignment):
    """Evaluate a `Sim` symbolic phase/outcome linear form (`bit 0` = const, `bit v>=1` = var `v`)."""
    out = form & 1
    vars_mask = form >> 1
    while vars_mask:
        lsb = vars_mask & -vars_mask
        v_id = lsb.bit_length()
        out ^= int(var_assignment.get(v_id, 0)) & 1
        vars_mask ^= lsb
    return out


def _solve_branch_assignment(forms_and_targets):
    """Find a concrete 0/1 variable assignment satisfying `[(form_i, bit_i), ...]` over GF(2)."""
    pivs = {}
    all_vars = set()
    for form, target in forms_and_targets:
        vmask = form >> 1
        rhs = (form & 1) ^ (target & 1)
        cur_v = vmask
        while cur_v:
            lsb = cur_v & -cur_v
            all_vars.add(lsb.bit_length())
            cur_v ^= lsb
        for p in sorted(pivs, reverse=True):
            if (vmask >> (p - 1)) & 1:
                vmask ^= pivs[p][0]
                rhs ^= pivs[p][1]
        if vmask == 0:
            if rhs != 0:
                raise ValueError("inconsistent branch outcome assignment")
            continue
        lead = vmask.bit_length()
        for p in list(pivs):
            if (pivs[p][0] >> (lead - 1)) & 1:
                pivs[p] = (pivs[p][0] ^ vmask, pivs[p][1] ^ rhs)
        pivs[lead] = (vmask, rhs)

    assign = {v: 0 for v in all_vars}
    for lead, (vmask, rhs) in sorted(pivs.items()):
        rem = vmask ^ (1 << (lead - 1))
        val = rhs
        while rem:
            lsb = rem & -rem
            v_id = lsb.bit_length()
            val ^= assign.get(v_id, 0)
            rem ^= lsb
        assign[lead] = val
    return assign


def _replay_encoded_vs_physical_control(code, basis_receipt, installed_bundle):
    """Replay encoded Bell-reference tableau across all 16 input Pauli frames & all 4 `(m0, m1_raw)` branches,
    checking branch-by-branch agreement against a direct 2-qubit physical reference circuit!
    """
    spec = installed_bundle["compiled_spec"]
    inst = Instrument(spec)
    shift = inst.c.n
    width = shift + code.n

    s0_x = pauli(basis_receipt["slots"][0]["x"], code.n)
    s0_z = pauli(basis_receipt["slots"][0]["z"], code.n)
    s1_x = pauli(basis_receipt["slots"][1]["x"], code.n)
    s1_z = pauli(basis_receipt["slots"][1]["z"], code.n)

    ro0 = installed_bundle["step_ledger"][0]["raw_physical_readout"]
    ro1_raw = installed_bundle["step_ledger"][1]["raw_physical_readout"]
    ro1_comp = spec["readouts"][1]

    checked_branches = 0
    # Test all 16 logical Pauli input frames `(fx0, fz0, fx1, fz1)` (including Y = iXZ)!
    for frame_bits in range(16):
        fx0 = (frame_bits >> 0) & 1
        fz0 = (frame_bits >> 1) & 1
        fx1 = (frame_bits >> 2) & 1
        fz1 = (frame_bits >> 3) & 1

        # Build physical Pauli frame P(f_x, f_z) applying logical X0^fx0 Z0^fz0 X1^fx1 Z1^fz1 on data!
        frame_pauli_str = fc_basis.embed_logical_pauli(
            basis_receipt,
            fx0 | (fx1 << 1),
            fz0 | (fz1 << 1),
            sign=0,
        )
        f_phys = pauli(frame_pauli_str, code.n)

        # Direct 2-qubit physical reference circuit:
        # Data qubits {0, 1}, Bell partners {2, 3}, initial Bell stabilizers X0X2, Z0Z2, X1X3, Z1Z3.
        # Input frame X0^fx0 Z0^fz0 X1^fx1 Z1^fz1 flips:
        #   Z0 Z2 -> (-1)^fx0,  X0 X2 -> (-1)^fz0,  Z1 Z3 -> (-1)^fx1,  X1 X3 -> (-1)^fz1.
        # Step 0 measures +Z0 with outcome m0 in {0, 1} -> post-state has (-1)^m0 Z0 = +1 and (-1)^{m0 ^ fx0} Z2 = +1.
        # Step 1 measures controlled (-1)^m0 X1 with raw outcome m1_raw in {0, 1}:
        #   (-1)^{m1_raw} ((-1)^m0 X1) = +1  <=>  (-1)^{m1_raw ^ m0} X1 = +1
        #   and on Bell partner 3: (-1)^{m1_raw ^ m0 ^ fz1} X3 = +1!

        sim = Sim(width)
        sim.gens = [
            Gen(
                g.x,
                g.z,
                g.s ^ fc_gf2.dot(g.x, f_phys.z) ^ fc_gf2.dot(g.z, f_phys.x),
            )
            for g in inst.initial
        ]
        sim.run(SimpleNamespace(ticks=inst.c.ticks[1:]))

        # All detectors must remain deterministic 0 under any logical Pauli input frame!
        for d_mask in inst.detectors:
            if sim.parity_form(fc_gf2.supp(d_mask)) != 0:
                raise ValueError("detector fired on clean logical Pauli input frame")

        form_m0 = sim.parity_form(ro0["records"]) ^ ro0["sign"]
        form_m1_raw = sim.parity_form(ro1_raw["records"]) ^ ro1_raw["sign"]
        form_m1_comp = sim.parity_form(ro1_comp["records"]) ^ ro1_comp["sign"]

        # Verify compiled readout algebraically equals raw_readout_1 XOR m0!
        if (form_m1_raw ^ form_m0) != form_m1_comp:
            raise ValueError("compiled feed-forward readout form disagrees with m1_raw XOR m0")

        for m0 in (0, 1):
            for m1_raw in (0, 1):
                assign = _solve_branch_assignment([(form_m0, m0), (form_m1_raw, m1_raw)])
                branch_gens = [
                    Gen(g.x, g.z, _eval_linear_form(g.s, assign)) for g in sim.gens
                ]
                b_grp = Group(width, branch_gens)

                # Expected physical-reference branch stabilizers embedded on encoded data + reference:
                want_z0_sign = s0_z.s ^ m0
                want_ref_z0_sign = s0_z.s ^ fc_gf2.dot(s0_z.x, s0_z.z) ^ m0 ^ fx0
                want_x1_sign = s1_x.s ^ m1_raw ^ m0
                want_ref_x1_sign = s1_x.s ^ fc_gf2.dot(s1_x.x, s1_x.z) ^ m1_raw ^ m0 ^ fz1

                for req in (
                    [Gen(chk.x, chk.z, chk.s) for chk in inst.checks]
                    + [Gen(chk.x << shift, chk.z << shift, chk.s ^ fc_gf2.dot(chk.x, chk.z)) for chk in inst.checks]
                    + [
                        Gen(s0_z.x, s0_z.z, want_z0_sign),
                        Gen(s0_z.x << shift, s0_z.z << shift, want_ref_z0_sign),
                        Gen(s1_x.x, s1_x.z, want_x1_sign),
                        Gen(s1_x.x << shift, s1_x.z << shift, want_ref_x1_sign),
                    ]
                ):
                    if any(not commute(req, h) for h in branch_gens):
                        raise ValueError("branch stabilizer fails to commute with final state")
                    red = b_grp.reduce(req)
                    if red.x or red.z or red.s:
                        raise ValueError(
                            f"branch mismatch at frame={frame_bits}, m0={m0}, m1_raw={m1_raw}"
                        )
                checked_branches += 1

    return {
        "logical_frames_tested": 16,
        "outcome_branches_per_frame": 4,
        "total_verified_branches": checked_branches,
        "first_read_outcomes_verified": [0, 1],
    }


def replay_lp6_2_addressable_experiment():
    """End-to-end Milestone 4 experiment on `LP(6,2)` (`c6860dd462fe`).

    1. Installs `LP(6,2)` (`[[6, 2, 2]]`) with generic verified `ftqc_codes-basis-1`,
       measures `+Z` on `q0`, and uses its recorded public bit to control the sign of
       `+X` on surviving `q1`.
    2. Replays all 16 input Pauli frames x 4 `(m0, m1_raw)` branches against a direct
       2-qubit physical reference circuit.
    3. Repeats with a stabilizer-dressed basis + verified `verify_basis_transform`.
    4. Verifies all 4 required negative controls fail closed:
       - slot swap without conversion,
       - stale code or basis digest,
       - reversed classical record sign,
       - omitted boundary/schedule idle.
    """
    lp6 = fc_codes.load("lp6_2")
    if not lp6.sha256().startswith("c6860dd462fe"):
        raise ValueError(f"unexpected LP(6,2) sha256 {lp6.sha256()}")
    static_rec = fc_certify.build(lp6)
    base_basis = fc_basis.derive_paired_basis(lp6)

    steps = [
        {"slot": "q0", "axis": "Z", "sign": 0, "rounds": 2, "method": "cone"},
        {"slot": "q1", "axis": "X", "sign": 0, "rounds": 2, "method": "cone", "control_from_step": 0},
    ]
    installed_base = install_addressable_program(
        static_rec,
        base_basis,
        steps,
        idles=True,
        require_explicit_idles=True,
        terminal_syndrome=True,
        target_distance=2,
        budget=200_000,
    )
    branch_base = _replay_encoded_vs_physical_control(lp6, base_basis, installed_base)

    # Repeat with stabilizer-dressed basis + verified basis conversion!
    dressed_basis, transform = fc_basis.dress_basis(
        lp6,
        base_basis,
        dress_x_hx=[[0], [1]],
        dress_z_hz=[[1], [0]],
    )
    installed_dressed = install_addressable_program(
        static_rec,
        dressed_basis,
        steps,
        old_basis_receipt=base_basis,
        transform_receipt=transform,
        idles=True,
        require_explicit_idles=True,
        terminal_syndrome=True,
        target_distance=2,
        budget=200_000,
    )
    branch_dressed = _replay_encoded_vs_physical_control(lp6, dressed_basis, installed_dressed)

    # Negative controls (all 4 must fail!)
    neg_controls = {}

    # (a) Slot swap without conversion: claiming `base_basis` slot `q0` target for swapped basis!
    swapped_basis, _ = fc_basis.dress_basis(lp6, base_basis, permute_slots=[1, 0])
    try:
        fc_basis.verify_logical_target(
            lp6,
            base_basis,
            [("q0", "Z")],
            fc_basis.embed_named_pauli(swapped_basis, [("q0", "Z")]),
        )
        raise AssertionError("slot swap without conversion was not refused")
    except fc_certify.Refusal as exc:
        neg_controls["slot_swap_without_conversion"] = f"refused: {exc}"

    # (b) Stale code_sha256 and stale basis_sha256
    stale_basis = copy.deepcopy(base_basis)
    stale_basis["code_sha256"] = "0" * 64
    try:
        verify_bound_code_and_basis(static_rec, stale_basis)
        raise AssertionError("stale code_sha256 was not refused")
    except ValueError as exc:
        neg_controls["stale_code_digest"] = str(exc)

    stale_bsha = copy.deepcopy(base_basis)
    stale_bsha["basis_sha256"] = "f" * 64
    try:
        verify_bound_code_and_basis(static_rec, stale_bsha)
        raise AssertionError("stale basis_sha256 was not refused")
    except ValueError as exc:
        neg_controls["stale_basis_digest"] = str(exc)

    # (c) Reversed classical record sign in compiled instrument
    bad_sign_spec = copy.deepcopy(installed_base["compiled_spec"])
    bad_sign_spec["readouts"][0]["sign"] ^= 1
    try:
        Instrument(bad_sign_spec)
        raise AssertionError("reversed record sign was not refused")
    except ValueError as exc:
        neg_controls["reversed_record_sign"] = str(exc)

    # (d) Omitted boundary/schedule idle
    try:
        install_addressable_program(
            static_rec,
            base_basis,
            steps,
            idles=False,
            require_explicit_idles=True,
        )
        raise AssertionError("omitted boundary idle was not refused")
    except ValueError as exc:
        neg_controls["omitted_boundary_idle"] = str(exc)

    return {
        "code_sha256": lp6.sha256(),
        "base_basis_sha256": base_basis["basis_sha256"],
        "dressed_basis_sha256": dressed_basis["basis_sha256"],
        "transform_sha256": installed_dressed["transform_sha256"],
        "base_schedule": {
            "installed_couplers": installed_base["installed_couplers"],
            "peak_live_qubits": installed_base["peak_live_qubits"],
            "ticks": installed_base["ticks"],
            "faults": installed_base["surgery_receipt"]["faults"],
            "branch_verification": branch_base,
        },
        "dressed_schedule": {
            "installed_couplers": installed_dressed["installed_couplers"],
            "peak_live_qubits": installed_dressed["peak_live_qubits"],
            "ticks": installed_dressed["ticks"],
            "faults": installed_dressed["surgery_receipt"]["faults"],
            "branch_verification": branch_dressed,
        },
        "negative_controls": neg_controls,
        "installed_base_bundle": installed_base,
    }


def compare_translated_seed_vs_ordinary_surgery():
    """Substitute a verified translated-seed gadget on canonical cyclic `LP(8,4)` (`arXiv:2607.28605`)
    and compare couplers, rounds, ancillas, ticks, and fault-distance intervals against
    repeated ordinary basis-derived gauging surgery on the same 2D logical Z subspace.
    """
    lp8 = fc_constructions.lifted_product([[5]], [[5]], 4, name="LP(8,4)")
    static_rec = fc_certify.build(lp8)
    basis_rec = fc_basis.derive_paired_basis(lp8)
    k = basis_rec["k"]

    # Find canonical weight-minimal Z seed whose shifts {0, 1} give 2 independent logical Z reads!
    chosen_seed = None
    orbit_cert = None
    for cand in sorted(fc_gf2.nullspace(lp8.hx, lp8.n), key=lambda w: (w.bit_count(), w)):
        if cand == 0:
            continue
        try:
            translated_settings(lp8.n, lp8.hx, lp8.hz, 4, cand, [0, 1], 2)
        except ValueError:
            continue
        cert = fc_symmetry.verify_seed_orbit_coverage(lp8, basis_rec, [cand], 4, sector="Z")
        if cert["covered_rank"] == 2:
            chosen_seed = cand
            orbit_cert = cert
            break
    if chosen_seed is None:
        raise ValueError("no rank-2 translated Z seed found on LP(8,4)")

    # 1. Translated-seed schedule (equivariant cyclic port translation of one seed!)
    t_designs = translated_settings(lp8.n, lp8.hx, lp8.hz, 4, chosen_seed, [0, 1], 2)
    from ftqc_codes import automorphisms as aut
    shift = fc_symmetry.cyclic_cell_shift_perm(lp8.n, 4)
    generator = aut.build(lp8, basis_rec, [{"perm": shift}])
    logical = (0, orbit_cert["orbits"][0]["translates"][0]["logical_mask"], 0)
    group = fc_basis._SignedStabilizerGroup(lp8.n, fc_basis._code_signed_checks(lp8))
    target_slots = []
    for design in t_designs:
        if logical[0] or logical[1].bit_count() != 1 or logical[2]:
            raise ValueError("LP comparison requires a positive single-slot Z setting")
        target_slots.append(basis_rec["slots"][logical[1].bit_length() - 1]["name"])
        physical = fc_basis.parse_pauli(design["targets"][0], lp8.n)
        expected = fc_basis.parse_pauli(fc_basis.embed_logical_pauli(basis_rec, *logical), lp8.n)
        delta = fc_basis.mul_pauli(physical, expected)
        if delta[2] % 2 or not group.contains_exact((delta[0], delta[1], delta[2] // 2)):
            raise ValueError("translated LP instrument measures the wrong signed logical")
        logical = aut.logical_image(generator, logical)
    t_parts = [build_surgery_design(d)["spec"] for d in t_designs]
    t_batches = []
    t_couplers = set()
    for p_spec in t_parts:
        w_loc = Instrument(p_spec).c.n
        pos = list(range(w_loc))
        t_batches.append([{"gadget": surgery_receipt(p_spec, 0), "placement": pos}])
        for row in p_spec["ticks"]:
            for op in row:
                if op[0].lstrip("!") == "CX":
                    t_couplers.add(tuple(sorted((pos[op[1]], pos[op[2]]))))
    t_sched = {
        "schema": 1,
        "kind": "installed_surgery",
        "n": lp8.n,
        "checks": t_parts[0]["checks"],
        "targets": [p_spec["targets"][0] for p_spec in t_parts],
        "couplers": [list(c) for c in sorted(t_couplers)],
        "batches": t_batches,
        "idles": True,
        "terminal_syndrome": True,
        "target": 2,
    }
    t_res = analyze_installed(t_sched, 200_000)

    # 2. Ordinary surgery on 2 slots using stabilizer-dressed basis representatives
    dz_dress = [[i % len(lp8.hz)] for i in range(k)]
    dressed_b, _ = fc_basis.dress_basis(lp8, basis_rec, dress_z_hz=dz_dress)
    ord_bundle = install_addressable_program(
        static_rec,
        dressed_b,
        [
            {"slot": target_slots[0], "axis": "Z", "rounds": 2, "method": "cone"},
            {"slot": target_slots[1], "axis": "Z", "rounds": 2, "method": "cone"},
        ],
        reuse_ancillas=True,
        idles=True,
        target_distance=2,
        budget=200_000,
    )

    return {
        "code": {"name": lp8.name, "n": lp8.n, "k": k, "d": static_rec["facts"]["d"], "sha256": lp8.sha256()},
        "seed_orbit_certificate": {
            "seed_support": fc_gf2.supp(chosen_seed),
            "signed_z_slots_in_order": target_slots,
            "seed_weight": chosen_seed.bit_count(),
            "covered_rank": orbit_cert["covered_rank"],
            "full_coverage": orbit_cert["full_coverage"],
            "covered_unit_slots": orbit_cert["covered_unit_slots"],
        },
        "translated_seed_gadget": {
            "installed_couplers": t_res["installed_couplers"],
            "rounds": 4,
            "peak_live_qubits": t_res["peak_live"],
            "peak_ancillas": t_res["peak_live"] - lp8.n,
            "ticks": t_res["ticks"],
            "faults": t_res["gadget"]["faults"],
        },
        "ordinary_basis_surgery": {
            "installed_couplers": ord_bundle["installed_couplers"],
            "rounds": 4,
            "peak_live_qubits": ord_bundle["peak_live_qubits"],
            "peak_ancillas": ord_bundle["peak_ancilla_qubits"],
            "ticks": ord_bundle["ticks"],
            "faults": ord_bundle["surgery_receipt"]["faults"],
        },
    }


def run_full_composed_chain():
    """Execute the complete `Codes -> Extractors/Gadgets -> Faults -> Decoders -> Compose -> Uncompute`
    chain on `LP(6,2)` (`c6860dd462fe`), checking detector/observable alignment, single-charged
    shared boundary faults across sequential stages, and 3-layout reaction-stall pricing.
    """
    from ftqc_compose.contracts import derive as derive_contract
    from ftqc_compose.surgery import channel as compose_surgery_channel
    import ftqc_decoders as fd
    from ftqc_faults.zx.diagram import CircuitOp, CliffordCircuit
    from ftqc_faults.zx.equivalence import circuit_distance
    from ftqc_uncompute import KickMixCircuit, estimate_reaction_decoder_latency

    lp6 = fc_codes.load("lp6_2")
    static_rec = fc_certify.build(lp6)
    basis_rec = fc_basis.derive_paired_basis(lp6)

    # Build two composable single-slot stages: Stage 1 measures Z on q0 (`surviving=['q1']`),
    # Stage 2 measures X on surviving q1 (`surviving=[]`)!
    stage1 = install_addressable_program(
        static_rec,
        basis_rec,
        [{"slot": "q0", "axis": "Z", "rounds": 2, "method": "cone"}],
        idles=True,
        target_distance=2,
        budget=200_000,
    )
    # After Stage 1 projects +Z_{q0}, Stage 2 runs on the post-Stage-1 stabilizer code (checks + +Z_{q0})
    # OR as a second standalone measurement on the code with `incoming=1` check!
    # Let us also test `op="seq"` in `ftqc_compose` by chaining Stage 1 (`LP(6,2)` -> `LP(6,2)` surviving `q1`)
    # with an output-matching Stage 2 whose input basis matches Stage 1's surviving output port!
    stage1_spec = stage1["compiled_spec"]
    stage1_rec = stage1["surgery_receipt"]

    # 1. Faults (`ftqc_faults.audit` & `ftqc_faults.zx`) on flattened Stage 1 circuit
    from ftqc_faults.audit.record_model import RecordModel

    audit_circ = AuditCircuit(
        [[("!R", q) for q in range(lp6.n)], *stage1_spec["ticks"]],
        n=Instrument(stage1_spec).c.n,
        detectors=[tuple(d) for d in stage1_spec["detectors"]],
        observables=[tuple(stage1_spec["readouts"][0]["records"])],
        outputs=list(range(lp6.n)),
        name="lp6_2_slot_q0_z",
    )
    audit_rec, _ = build_audit_receipt(audit_circ, max_weight=1)
    verify_audit_receipt(audit_rec, full=True)

    # Record-endpoint pure-CSS homomorphic counterpart (`R, RX, CX, I, M`) for `RecordModel -> DecodingProblem`!
    stage1_hom = install_addressable_program(
        static_rec,
        basis_rec,
        [{"slot": "q0", "axis": "Z", "rounds": 2, "method": "homomorphic"}],
        idles=True,
        target_distance=2,
        budget=200_000,
    )
    hom_spec = stage1_hom["compiled_spec"]
    hom_audit = AuditCircuit(
        [[("!R", q) for q in range(lp6.n)], *hom_spec["ticks"]],
        n=Instrument(hom_spec).c.n,
        detectors=[tuple(d) for d in hom_spec["detectors"]],
        outputs=list(range(lp6.n)),
    )
    n_mid_meas = len(hom_audit.meas)
    term_dets = [tuple(d) for d in hom_spec["detectors"]] + [
        tuple(n_mid_meas + q for q in fc_gf2.supp(z_chk)) for z_chk in lp6.hz
    ]
    z_q0_mask = fc_basis.parse_pauli(basis_rec["slots"][0]["z"], lp6.n)[1]
    term_obs = [
        tuple(hom_spec["readouts"][0]["records"]),
        tuple(n_mid_meas + q for q in fc_gf2.supp(z_q0_mask)),
    ]
    rec_circ = AuditCircuit(
        [
            [("!R", q) for q in range(lp6.n)],
            *hom_spec["ticks"],
            [("!M", q) for q in range(lp6.n)],
        ],
        n=hom_audit.n,
        detectors=term_dets,
        observables=term_obs,
        outputs=[],
        name="lp6_2_slot_q0_z_record_endpoint",
    )
    record_model = RecordModel(rec_circ)

    # Also build open-wire Clifford ZX diagram of the core CNOT/measurement layer for ZX distance sanity
    zx_ops = []
    for row in stage1_spec["ticks"]:
        for op in row:
            nm = op[0].lstrip("!")
            if nm == "R":
                zx_ops.append(CircuitOp("PREP_0", (op[1],)))
            elif nm == "RX":
                zx_ops.append(CircuitOp("PREP_PLUS", (op[1],)))
            elif nm == "H":
                zx_ops.append(CircuitOp("H", (op[1],)))
            elif nm == "CX":
                zx_ops.append(CircuitOp("CX", (op[1], op[2])))
            elif nm == "M":
                zx_ops.append(CircuitOp("MEAS_Z", (op[1],)))
            elif nm == "MX":
                zx_ops.append(CircuitOp("MEAS_X", (op[1],)))
    zx_circ = CliffordCircuit(
        audit_circ.n,
        tuple(range(lp6.n)),
        tuple(range(lp6.n)),
        zx_ops,
        name="lp6_2_q0_zx",
    )
    zx_dist = circuit_distance(zx_circ, w_max=1).to_dict()

    # 2. Decoders (`ftqc_decoders`): feed identical detector/observable orders from `record_model`
    # AND static `LP(6,2)` CSS problem across `bp_osd`, `bp_lsd`, `relay_bp`!
    dem_prob = fd.DecodingProblem.from_record_model(
        record_model,
        name="lp6_2_stage1_record_model",
    )
    if dem_prob.num_checks != len(term_dets) or dem_prob.num_observables != len(term_obs):
        raise ValueError("detector/observable count mismatch between ftqc_faults and ftqc_decoders")

    lp6_x_prob = fd.DecodingProblem.from_css_code(lp6, sector="X")
    dec_results = {}
    for alg in ("bp_osd", "bp_lsd", "relay_bp"):
        r_dec = fd.run_decoder(lp6_x_prob, lp6_x_prob.syndrome_of(0), alg)
        r_dem = fd.run_decoder(dem_prob, dem_prob.syndrome_of(1) if dem_prob.num_mechanisms else 0, alg)
        dec_results[alg] = {
            "static_converged": r_dec.converged,
            "circuit_dem_converged": r_dem.converged,
        }

    # 3. Compose (`ftqc_compose`): reverify producer receipt, named basis, residual filter,
    #    and single-charged boundary wires across a 2-stage `seq` composition!
    from ftqc_compose.contracts import Violation

    wires = [f"lp6_q{i}" for i in range(lp6.n)]
    leaf1 = derive_contract(
        {
            "op": "surgery",
            "instance": "stage1_q0_z",
            "wires": wires,
            "receipt": stage1_rec,
            "incoming": 0,
        }
    )
    # On d=2 LP(6,2), claiming incoming=1 simultaneously with 1 internal fault (1+1=2=d) must raise Violation!
    try:
        derive_contract(
            {
                "op": "surgery",
                "instance": "stage1_overclaim",
                "wires": wires,
                "receipt": stage1_rec,
                "incoming": 1,
            }
        )
        raise AssertionError("expected Violation when claiming incoming=1 on d=2 LP(6,2)")
    except Violation:
        overclaim_refused = True

    chan1 = compose_surgery_channel(
        {
            "kind": "surgery_channel",
            "receipt": stage1_rec,
            "incoming_weight": 1,
            "p": "1/1000",
            "cap": 200_000,
        }
    )

    # Second sequential leaf on the same 6 boundary data wires (re-measuring idempotent +Z_{q0}
    # with basis stripped to physical signed Pauli so input/output port schemas match directly in `seq`,
    # verifying boundary wires `lp6_q0..lp6_q5` are single-charged at the interface rather than double-counted!)
    raw_s1 = copy.deepcopy(stage1_spec)
    raw_s1.pop("basis", None)
    raw_r1 = surgery_receipt(raw_s1, order=1, cap=200_000)
    s1_bit0_port = {
        "kind": "C",
        "wires": ["s1/bit/0"],
        "code": "bit",
        "basis": "Z",
        "decoder": "identity",
        "bad": [["s1/bit/0"]],
    }
    seq_thm = derive_contract(
        {
            "op": "sequence",
            "left": {"op": "surgery", "instance": "s1", "wires": wires, "receipt": raw_r1, "incoming": 0},
            "right": {
                "op": "parallel",
                "left": {"op": "surgery", "instance": "s2", "wires": wires, "receipt": raw_r1, "incoming": 0},
                "right": {
                    "op": "hypothesis",
                    "id": "carry_s1_bit0",
                    "inputs": [s1_bit0_port],
                    "outputs": [s1_bit0_port],
                    "ideal": "identity-bit",
                    "faults": {"op": "family", "sites": [], "minimal": []},
                    "filters": ["identity"],
                },
            },
        }
    )
    s1_locs = sum(not op[0].startswith("!") for row in raw_s1["ticks"] for op in row)

    # 4. Uncompute (`ftqc_uncompute`): price reaction depth + decoder stalls across 3 layouts!
    kmx = KickMixCircuit.from_kmx("HMR q0 b0\nZ q1 if b0\nHMR q1 b1\nZ q0 if b1")
    ticks_per_round = max(1, stage1["ticks"] // 2)
    sched_metrics = {
        "num_data_qubits": lp6.n,
        "peak_ancilla_qubits": stage1["peak_ancilla_qubits"],
        "rounds_per_reaction": 2,
        "ticks_per_round": ticks_per_round,
        "installed_couplers": stage1["installed_couplers"],
    }
    layouts = {}
    for lname, stall, shuttle in (
        ("planar_local", "5/1", 0),
        ("multitier_coupler", "3/1", 1),
        ("reconfigurable_shuttle", "7/2", 4),
    ):
        layouts[lname] = estimate_reaction_decoder_latency(
            kmx,
            sched_metrics,
            decoder_stall_ticks_per_reaction=stall,
            shuttle_ticks_per_reaction=shuttle,
            surviving_logicals=1,
            bare_code_k=2,
            layout_name=lname,
        )

    return {
        "audit": {
            "num_ticks": len(audit_circ.ticks),
            "num_measurements": len(audit_circ.meas),
            "num_detectors": len(term_dets),
            "num_observables": len(term_obs),
            "verified": True,
        },
        "zx": {"open_wire_check_w1": zx_dist},
        "decoders": {
            "dem_checks": dem_prob.num_checks,
            "dem_mechanisms": dem_prob.num_mechanisms,
            "algorithms": dec_results,
        },
        "compose": {
            "leaf1_outputs": len(leaf1["outputs"]),
            "leaf1_surviving_basis": json.loads(leaf1["outputs"][0]["basis"]),
            "d2_incoming1_overclaim_refused": overclaim_refused,
            "channel_status": chan1["status"],
            "seq_single_stage_locations": s1_locs,
            "seq_total_fault_universe_size": len(seq_thm["family"].sites),
            "boundary_wires_single_charged": len(seq_thm["family"].sites) == 2 * s1_locs,
        },
        "uncompute_pricing": layouts,
    }


def build_bb144_decoder_study(*, max_splits=6):
    """Scale-up Study 1: `BB144` (`[[144, 12, 12]]`, `0ae98ac8a6e5`) static memory & Gross
    circuit-level syndrome extraction compared across `BP+OSD`, `BP+LSD`, and `Relay-BP`.
    """
    import ftqc_decoders as fd
    from pathlib import Path

    bb144 = fc_codes.load("bb144")
    crec = fc_constructions.build_construction_receipt(
        bb144,
        provenance={"paper": "Bravyi et al., Nature 627, 778 (2024)", "arxiv": "2308.07915", "locator": "Table 1 [[144,12,12]]"},
        paper_claims={"d": 12, "d_x": 12, "d_z": 12},
    )
    fc_constructions.verify_construction_receipt(crec)

    static_prob = fd.bb144_gross_problem()
    # Circuit-level BB144 Gross Z syndrome-extraction DEM (`gross_probes.experiment` -> `RecordModel`)
    from ftqc_faults.audit.record_model import RecordModel
    from .gross_probes import experiment as gross_probe_experiment

    gross_circ_prob = fd.DecodingProblem.from_record_model(
        RecordModel(gross_probe_experiment(2, "Z")),
        name="BB144_GrossZ_2R_RecordModel",
    )

    # Compare decoders on adversarial split-stabilizer trapping sets + operational radius w=1 + circuit faults
    ts = fd.search_trapping_sets(
        static_prob,
        "bp_only",
        max_exhaustive_weight=0,
        max_splits_per_generator=max_splits,
        compare_decoders=["bp_osd", "bp_lsd", "relay_bp"],
    )

    static_bench = {}
    circ_bench = {}
    # Deterministic weight-2 static errors and weight-1/weight-2 circuit-level faults
    static_probes = [
        (1 << 0) | (1 << 73),
        (1 << 11) | (1 << 58) | (1 << 119),
        sum(1 << b for b in fd.bit_positions(static_prob.stabilizer_generators[0])[:3]),
    ]
    circ_probes = [
        1 << 0,
        (1 << 5) | (1 << 42),
        (1 << 17) | (1 << 88),
    ]

    specs = [("bp_osd", {"algorithm": "bp_osd", "osd_order": 0}),
             ("bp_lsd", {"algorithm": "bp_lsd", "osd_order": 0}),
             ("relay_bp", {"algorithm": "relay_bp", "num_relays": 16})]
    backends = fd.probe_external_backends()
    if backends.get("ldpc", {}).get("available"):
        specs.append(("external:ldpc_bplsd", {"algorithm": "external:ldpc_bplsd", "osd_order": 0}))

    for label, dspec in specs:
        rad1 = fd.certify_operational_radius(static_prob, dspec, 1, budget_cap=500)
        s_healed = 0
        max_elim = 0
        max_clust = 0
        for em in static_probes:
            res = fd.run_decoder(static_prob, static_prob.syndrome_of(em), dspec)
            if res.converged and static_prob.observable_of(em ^ res.error_mask) == 0:
                s_healed += 1
            max_elim = max(max_elim, int(res.stats.get("eliminated_rank", 0)))
            max_clust = max(max_clust, int(res.stats.get("max_cluster_rank", 0)))
        static_bench[label] = {
            "radius_w1_status": rad1["status"],
            "holds_through": rad1["holds_through"],
            "probes_succeeded": f"{s_healed}/{len(static_probes)}",
            "max_eliminated_rank": max_elim,
            "max_cluster_rank": max_clust,
        }

        c_ok = 0
        for cm in circ_probes:
            if cm >> gross_circ_prob.num_mechanisms:
                continue
            rc = fd.run_decoder(gross_circ_prob, gross_circ_prob.syndrome_of(cm), dspec)
            if rc.converged and gross_circ_prob.observable_of(cm ^ rc.error_mask) == 0:
                c_ok += 1
        circ_bench[label] = {
            "circuit_probes_succeeded": f"{c_ok}/{len(circ_probes)}",
        }

    peierls = fd.enclose_peierls_cluster(static_prob, "1/10000", verified_component_radius=5, digits=18)
    return {
        "study": "bb144_memory_and_extractor_decoders",
        "code_sha256": bb144.sha256(),
        "construction_sha256": crec["construction_sha256"],
        "budget": {"max_splits_per_generator": max_splits, "radius_budget_cap": 500, "seed": 20260601},
        "static_problem": {"checks": static_prob.num_checks, "mechanisms": static_prob.num_mechanisms, "observables": static_prob.num_observables},
        "circuit_problem": {"checks": gross_circ_prob.num_checks, "mechanisms": gross_circ_prob.num_mechanisms, "observables": gross_circ_prob.num_observables},
        "trapping_set_probe": {
            "bp_only_status": ts["status"],
            "fault_weight": ts["witness"]["fault_weight"],
            "healed_by": {c["decoder"]: c["healed"] for c in ts["comparison"]},
        },
        "static_decoder_comparison": static_bench,
        "circuit_decoder_comparison": circ_bench,
        "peierls_upper_bound_p1e4": peierls["ler_hi"],
    }


def build_radial90_map_comparison_study():
    """Scale-up Study 2: `radial90` (`[[90, 8, 10]]`, `arXiv:2609.11723`) Z_5-homotopy-sparsified
    lifted measurement map vs dense rank-one `radial_design()` at matched `rounds=1` and boundaries.
    """
    import ftqc_decoders as fd

    r_code = fc_codes.load("radial90")
    dense_spec = radial_design()
    hom_spec = radial_homotopy_design(steps=240, seed=20260910)

    dense_built = build_surgery_design(dense_spec)
    hom_built = build_surgery_design(hom_spec)

    z_ech = fc_gf2.Echelon(r_code.hz)
    same_coset = (
        z_ech.reduce(dense_built["geometry"]["cone"]["measured"][0])
        == z_ech.reduce(hom_built["geometry"]["cone"]["measured"][0])
    )
    if not same_coset:
        raise ValueError("homotopy-sparsified radial90 map left the target logical coset")

    def summarize_variant(label, d_in, built):
        sp = built["spec"]
        cone = built["geometry"]["cone"]
        inst = Instrument(sp)
        rec = surgery_receipt(sp, order=1, cap=500_000)
        cx_count = sum(1 for row in sp["ticks"] for op in row if op[0].lstrip("!") == "CX")
        # Decode merged static cone X-errors under Z-checks (`cone["hz"]` has 90 checks on 135 merged qubits)
        merged_css = fc_codes.CSSCode(f"Radial90_{label}", cone["n"], cone["hx"], cone["hz"])
        m_prob = fd.DecodingProblem.from_css_code(merged_css, sector="X")
        dec_ok = {}
        probe_err = (1 << 0) | (1 << 47)
        for alg in ("bp_osd", "bp_lsd", "relay_bp"):
            r_d = fd.run_decoder(m_prob, m_prob.syndrome_of(probe_err), alg)
            dec_ok[alg] = bool(r_d.converged and m_prob.observable_of(probe_err ^ r_d.error_mask) == 0)
        return {
            "rounds": d_in["rounds"],
            "terminal_syndrome": sp["terminal_syndrome"],
            "merged_qubits": cone["n"],
            "merged_logicals": cone["logical_qubits"],
            "max_check_weight": cone["max_check"],
            "max_qubit_degree": cone["max_qubit_degree"],
            "port_nonzeros": sum(p.bit_count() for p in d_in["ports"]),
            "cx_gates": cx_count,
            "physical_wires": inst.c.n,
            "ticks": len(sp["ticks"]),
            "fault_locations": sum(1 for _ in inst.c.locations()),
            "fault_audit_order1": rec["faults"],
            "merged_cone_decoded_w2": dec_ok,
            "circuit_sha256": rec["circuit_sha256"],
        }

    dense_sum = summarize_variant("dense_rank1", dense_spec, dense_built)
    hom_sum = summarize_variant("z5_homotopy_sparse", hom_spec, hom_built)

    return {
        "study": "radial90_lifted_map_vs_dense_reference",
        "provenance": {"paper": "Berent, Cohen & Quintavalle, Lifted surgery (2026)", "arxiv": "2609.11723", "locator": "Section IV.E, [[90,8,10]]"},
        "code_sha256": r_code.sha256(),
        "matched_rounds": 1,
        "matched_boundary_terminal_syndrome": True,
        "same_logical_coset_verified": same_coset,
        "dense_rank_one_reference": dense_sum,
        "z5_homotopy_sparsified": hom_sum,
        "degree_reduction": {
            "max_check_delta": dense_sum["max_check_weight"] - hom_sum["max_check_weight"],
            "max_qubit_degree_delta": dense_sum["max_qubit_degree"] - hom_sum["max_qubit_degree"],
            "cx_gate_delta": dense_sum["cx_gates"] - hom_sum["cx_gates"],
        },
    }


def _edge_color_bipartite_cnot_depth(check_rows):
    """Compute a deterministic conflict-free CNOT tick coloring for check extraction."""
    # Greedy Vizing-style coloring of bipartite check-to-data edges: each tick uses disjoint wires!
    edges = []
    for c_i, r in enumerate(check_rows):
        for q in fc_gf2.supp(r):
            edges.append((c_i, q))
    used_c = {}
    used_q = {}
    max_color = 0
    for c_i, q in edges:
        bc = used_c.get(c_i, 0) | used_q.get(q, 0)
        col = 0
        while (bc >> col) & 1:
            col += 1
        used_c[c_i] = used_c.get(c_i, 0) | (1 << col)
        used_q[q] = used_q.get(q, 0) | (1 << col)
        if col + 1 > max_color:
            max_color = col + 1
    return max_color, len(edges)


def build_pp256_orbit_surgery_study():
    """Scale-up Study 3: Maskara et al. (`arXiv:2610.06749v1` Table V & Eq. C1-C12) `PP256`
    (`[[256, 68, 12]]`) parallel 34-logical orbit surgery (`[[576, 34, 12]]`) vs 34-step serial
    measurement of the *same signed X subgroup*, strictly separating:
      - bare `k/n` vs physical qubits per surviving logical,
      - `static`, `phenomenological`, `circuit`, and `decoded` claims.
    """
    import ftqc_decoders as fd

    pp = fc_highrate.pp256_68()
    crec = fc_constructions.build_construction_receipt(pp)
    fc_constructions.verify_construction_receipt(crec)

    n, hx, hz = pp.n, pp.hx, pp.hz
    P = 32
    # Attachment gamma_1: Z-pivot blocks S_Z={3,4,5} (96 cols) + pair program (6, 7, shift=13) (32 cols) = 128 cols!
    gamma1 = [1 << (b * P + j) for b in (3, 4, 5) for j in range(P)] + [
        (1 << (6 * P + j)) | (1 << (7 * P + ((j + 13) % P))) for j in range(P)
    ]
    dP_cols = [sum(fc_gf2.dot(hz[r], gamma1[c]) << r for r in range(96)) for c in range(128)]

    # Verify measured 34-logical X subgroup via `restricted_cone` on starter m=1!
    aux_m1 = [sum(((dP_cols[c] >> a) & 1) << c for c in range(128)) for a in range(96)]
    cone_m1 = restricted_cone(n, hz, hx, gamma1, aux_m1)
    measured_x_basis = list(cone_m1["measured"])
    if len(measured_x_basis) != 34 or cone_m1["logical_qubits"] != 34:
        raise ValueError("unexpected PP256 measured subgroup dimension")

    # Build exact Table V m=2 spatial-convolved (`f = 1 + x + x^3`, `g = 1`) merged CSS code `[[576, 34]]`!
    f_exps, g_exps = (0, 1, 3), (0,)
    mx_rows = list(hx)
    for c in range(128):
        blk, pos = divmod(c, P)
        r1 = gamma1[c] | (dP_cols[c] << 256)
        for e in f_exps:
            r1 |= 1 << (448 + blk * P + ((pos + e) % P))
        mx_rows.append(r1)
    for c in range(128):
        blk, pos = divmod(c, P)
        r2 = dP_cols[c] << 352
        for e in g_exps:
            r2 |= 1 << (448 + blk * P + ((pos + e) % P))
        mx_rows.append(r2)

    mz_rows = [hz[r] | (1 << (256 + r)) for r in range(96)]
    for r in range(96):
        blk, pos = divmod(r, P)
        zr = 0
        for e in f_exps:
            zr |= 1 << (256 + blk * P + ((pos - e) % P))
        for e in g_exps:
            zr |= 1 << (352 + blk * P + ((pos - e) % P))
        for c in range(128):
            if (dP_cols[c] >> r) & 1:
                zr |= 1 << (448 + c)
        mz_rows.append(zr)

    merged576 = fc_codes.CSSCode("PP256_Merged576", 576, mx_rows, mz_rows)
    if not merged576.commutes():
        raise ValueError("PP256_Merged576 checks do not commute")
    k_surv = 576 - fc_gf2.rank(mx_rows) - fc_gf2.rank(mz_rows)
    if k_surv != 34 or len(mx_rows) != 352 or len(mz_rows) != 192:
        raise ValueError("PP256_Merged576 disagrees with Maskara Table V parameters")

    # Verify algebraically that EVERY measured X-operator in `measured_x_basis` (34 independent logicals!)
    # becomes an X-stabilizer on data when projected from `rowspace(mx_rows)` with zero ancilla support!
    # Specifically, vectors in `rowspace(mx_rows)` supported purely on `0..255` modulo `rowspace(hx)` must equal `span(measured_x_basis)`!
    ech_hx = fc_gf2.Echelon(hx)
    # In `mx_rows`, rows 96..223 have `z1`, rows 224..351 have `z2`: any combination `(z1, z2)` supported on `0..255`
    # has `f z1 + g z2 = 0 => z2 = f z1`, `dP z2 = 0`, `dP z1 = 0 => dP z1 = 0` (`ker dP`), producing `gamma1(z1)`!
    ker_dP = fc_gf2.nullspace([sum(((dP_cols[c] >> r) & 1) << c for c in range(128)) for r in range(96)], 128)
    m2_data_stabilized = []
    seen_m2 = fc_gf2.Echelon(hx)
    for z1_mask in ker_dP:
        # Build `z2_mask = f * z1_mask` and verify `mx_rows` combination has ZERO ancilla bits on `256..575`!
        comb = 0
        for c in fc_gf2.supp(z1_mask):
            comb ^= mx_rows[96 + c]
            blk, pos = divmod(c, P)
            for e in f_exps:
                c2 = blk * P + ((pos + e) % P)
                comb ^= mx_rows[224 + c2]
        if (comb >> 256) != 0:
            raise ValueError("m=2 starter cycle left nonzero ancilla support")
        if seen_m2.add(comb):
            m2_data_stabilized.append(comb)
    # Verify exact equality of 34-dim measured X-subspace mod `hx`!
    sub_m1 = fc_gf2.Echelon(hx + measured_x_basis)
    sub_m2 = fc_gf2.Echelon(hx + m2_data_stabilized)
    same_subgroup = (
        len(m2_data_stabilized) == 34
        and sub_m1.rank == sub_m2.rank
        and all(v in sub_m1 for v in m2_data_stabilized)
    )
    if not same_subgroup:
        raise ValueError("m=2 spatial convolution measured a different X-subgroup than m=1")

    # Serial baseline measuring the EXACT SAME 34 signed `+X` generators one by one via gauging graphs!
    serial_weights = [w.bit_count() for w in measured_x_basis]
    max_serial_wt = max(serial_weights)
    # For weight `w` cone gauging: `w` vertex ancillas + `w-1` edge ancillas = `2w - 1` ancillas
    # plus `192` bare-code check meters if maintaining full syndrome extraction (`256 + 192 + 2*w_max - 1`)!
    serial_peak_anc = 2 * max_serial_wt - 1
    parallel_total_with_meters = 576 + 352 + 192  # Eq. (F12): n_merged + m_X + m_Z = 1120

    # 2. Phenomenological 3-round spacetime Z-error / X-check DEM for `PP256_Merged576` (`rounds = 3`)
    rounds_par = 3
    rounds_ser = 34 * rounds_par
    # Build 3-round repeated phenomenological Z-fault + X-syndrome-measurement DEM:
    # Mechanisms: `3 * 576` data Z-faults + `2 * 352` noisy intermediate X-check measurement flips (`2432` mechs),
    # Detectors: `3 * 352 = 1056` difference syndromes (with clean terminal round)!
    phen_checks = []
    # Let mech indices be: round r in 0..2 data q in 0..575 -> r*576 + q; meas err round r in 0..1 check c in 0..351 -> 1728 + r*352 + c
    for r_i in range(rounds_par):
        for c_i, xrow in enumerate(mx_rows):
            det_row = xrow << (r_i * 576)
            if r_i < rounds_par - 1:
                det_row |= 1 << (1728 + r_i * 352 + c_i)
            if r_i > 0:
                det_row |= 1 << (1728 + (r_i - 1) * 352 + c_i)
            phen_checks.append(det_row)
    # Surviving X-logical observables accumulate across all 3 data rounds!
    ech_mx = fc_gf2.Echelon(mx_rows)
    surv_x_logs = [v for v in fc_gf2.nullspace(mz_rows, 576) if ech_mx.add(v)]
    phen_obs = [v | (v << 576) | (v << 1152) for v in surv_x_logs]
    phen_prob = fd.DecodingProblem.from_matrices(
        phen_checks,
        phen_obs,
        1728 + 2 * 352,
        name="PP256_Merged576_Phenom3R",
    )

    # 3. Circuit schedule metrics (conflict-free bipartite CNOT edge-coloring)
    cnot_depth_x, cx_x = _edge_color_bipartite_cnot_depth(mx_rows)
    cnot_depth_z, cx_z = _edge_color_bipartite_cnot_depth(mz_rows)
    ticks_per_merged_round = 2 + cnot_depth_x + cnot_depth_z  # prep + X-CNOTs + Z-CNOTs + meas
    parallel_circuit_ticks = rounds_par * ticks_per_merged_round

    base_depth_x, base_cx_x = _edge_color_bipartite_cnot_depth(hx)
    base_depth_z, base_cx_z = _edge_color_bipartite_cnot_depth(hz)
    ticks_per_serial_round = 2 + base_depth_x + base_depth_z + 4  # + gauging vertex/edge CNOTs
    serial_circuit_ticks = rounds_ser * ticks_per_serial_round

    # 4. Decoded claims (static merged Z-sector & 3-round phenomenological data+syndrome recovery!)
    static_merged_z = fd.DecodingProblem.from_css_code(merged576, sector="Z")
    decoded_summary = {}
    s_err = (1 << 10) | (1 << 277) | (1 << 460)
    p_err = (1 << 10) | (1 << (576 + 205)) | (1 << (1728 + 44))  # data faults in rounds 0,1 + meas flip!

    dec_list = ["bp_osd", "bp_lsd", "relay_bp"]
    if fd.probe_external_backends().get("ldpc", {}).get("available"):
        dec_list.append("external:ldpc_bplsd")

    for alg in dec_list:
        rs = fd.run_decoder(static_merged_z, static_merged_z.syndrome_of(s_err), alg)
        rp = fd.run_decoder(phen_prob, phen_prob.syndrome_of(p_err), alg)
        decoded_summary[alg] = {
            "static_w3_converged": rs.converged,
            "static_w3_logical_ok": bool(rs.converged and static_merged_z.observable_of(s_err ^ rs.error_mask) == 0),
            "phenom_3round_converged": rp.converged,
            "phenom_3round_logical_ok": bool(rp.converged and phen_prob.observable_of(p_err ^ rp.error_mask) == 0),
        }

    return {
        "study": "pp256_parallel_orbit_vs_serial_subgroup",
        "provenance": pp.meta["provenance"],
        "base_code_sha256": pp.sha256(),
        "merged576_sha256": merged576.sha256(),
        "same_34_generator_signed_x_subgroup_verified": same_subgroup,
        "space_accounting": {
            "bare_code_parameters": [256, 68],
            "bare_code_rate_k_over_n": "68/256",
            "bare_data_qubits_per_logical": "256/68",
            "measured_subgroup_rank": 34,
            "surviving_logicals_during_surgery": 34,
            "parallel_merged_code_qubits": 576,
            "parallel_merged_qubits_per_surviving_logical": "576/34",
            "parallel_total_physical_qubits_with_meters": parallel_total_with_meters,
            "parallel_physical_qubits_per_surviving_logical": f"{parallel_total_with_meters}/34",
            "serial_peak_code_plus_gauging_qubits": 256 + serial_peak_anc,
            "serial_total_physical_qubits_with_meters": 256 + 192 + serial_peak_anc,
            "serial_physical_qubits_per_surviving_logical": f"{256 + 192 + serial_peak_anc}/34",
        },
        "claims_by_model": {
            "static": {
                "base_code": {"n": 256, "k": 68, "max_check_weight": 8, "paper_claims": {"d": 12}, "proved_lower": None},
                "parallel_starter_m1": {"n": 352, "k_surviving": 34, "m_x": 224, "m_z": 96},
                "parallel_spatial_m2": {
                    "n": 576,
                    "k_surviving": 34,
                    "m_x": 352,
                    "m_z": 192,
                    "max_check_weight_x": max(r.bit_count() for r in mx_rows),
                    "max_check_weight_z": max(r.bit_count() for r in mz_rows),
                    "paper_claims": {"d": 12},
                    "proved_lower": None,
                },
                "serial_34_steps": {
                    "steps": 34,
                    "generator_weight_min": min(serial_weights),
                    "generator_weight_max": max_serial_wt,
                },
            },
            "phenomenological": {
                "parallel_rounds": rounds_par,
                "serial_rounds_matched_per_step": rounds_ser,
                "round_reduction_factor": "34/1",
                "parallel_dem_checks": phen_prob.num_checks,
                "parallel_dem_mechanisms": phen_prob.num_mechanisms,
                "parallel_surviving_observables": phen_prob.num_observables,
            },
            "circuit": {
                "parallel_cnot_depth_per_round": cnot_depth_x + cnot_depth_z,
                "parallel_cnots_per_round": cx_x + cx_z,
                "parallel_total_ticks_3rounds": parallel_circuit_ticks,
                "parallel_qubit_ticks_3rounds": parallel_circuit_ticks * parallel_total_with_meters,
                "serial_total_ticks_102rounds": serial_circuit_ticks,
                "serial_qubit_ticks_102rounds": serial_circuit_ticks * (256 + 192 + serial_peak_anc),
            },
            "decoded": decoded_summary,
        },
    }
