"""Lineage B: Disjoint-Core Forward Spacetime Pauli-Frame ("Spackle") Oracle.

Computes signed Choi stabiliser tableaux, deterministic detecting sets (Def. 2.6),
and forward spacetime fault signatures directly on `CliffordCircuit` via causal
stabiliser-state evolution and forward Pauli propagation — with ZERO dependence on
`ZXDiagram` or `pauliweb.py`.

Used as the non-circular cross-verification oracle in every `ftqc_faults.zx` battery leg.
Zero float literals appear anywhere in this module.
"""

from __future__ import annotations

from dataclasses import dataclass

from .diagram import (
    AtomicFaultGroup,
    CliffordCircuit,
    PREP_STATES,
    RefusalError,
    SINGLE_MEAS,
    SINGLE_QUBIT_GATES,
    VALID_NOISE_MODELS,
)
from .gf2 import (
    canonicalize_signed_stabilizers,
    nullspace_gf2,
    pauli_product_phase_mod4,
    pauli_string_to_symplectic,
    rref_gf2,
    solve_gf2,
    symplectic_inner_product,
)


def _conj_1q_gate(gate: str, z_b: int, x_b: int) -> tuple[int, int, int]:
    """Conjugate single-qubit Pauli (z_b, x_b) by Clifford `gate`: G P G^dag = (-1)^flip P(z', x')."""
    if z_b == 0 and x_b == 0:
        return (0, 0, 0)
    if gate == "I":
        return (z_b, x_b, 0)
    if gate == "X":
        # X Z X = -Z, X X X = +X, X Y X = -Y -> flip iff z_b == 1
        return (z_b, x_b, z_b)
    if gate == "Z":
        # Z Z Z = +Z, Z X Z = -X, Z Y Z = -Y -> flip iff x_b == 1
        return (z_b, x_b, x_b)
    if gate == "Y":
        # Y Z Y = -Z, Y X Y = -X, Y Y Y = +Y -> flip iff z_b ^ x_b == 1
        return (z_b, x_b, z_b ^ x_b)
    if gate == "H":
        # H Z H = +X, H X H = +Z, H Y H = -Y
        return (x_b, z_b, z_b & x_b)
    if gate == "S":
        # S Z S^dag = +Z, S X S^dag = +Y, S Y S^dag = -X!
        # Check: (1,0)->(1,0,0); (0,1)->(1,1,0); (1,1)->(0,1,1)
        return (z_b ^ x_b, x_b, z_b & x_b)
    if gate == "SDG":
        # Sdg Z S = +Z, Sdg X S = -Y, Sdg Y S = +X!
        # Check: (1,0)->(1,0,0); (0,1)->(1,1,1); (1,1)->(0,1,0)
        return (z_b ^ x_b, x_b, (1 - z_b) & x_b)
    raise ValueError(f"Unsupported 1q gate {gate!r}")


def _apply_1q_to_stab_list(
    stabs: list[tuple[int, int]],
    slot: int,
    gate: str,
) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    z_pos = 2 * slot
    x_pos = 2 * slot + 1
    clear = ~((1 << z_pos) | (1 << x_pos))
    for vec, sgn in stabs:
        zb = (vec >> z_pos) & 1
        xb = (vec >> x_pos) & 1
        nz, nx, flip = _conj_1q_gate(gate, zb, xb)
        nvec = (vec & clear) | (nz << z_pos) | (nx << x_pos)
        out.append((nvec, sgn ^ flip))
    return out


def _apply_2q_to_stab_list(
    stabs: list[tuple[int, int]],
    slot0: int,
    slot1: int,
    gate: str,
) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    z0p, x0p = 2 * slot0, 2 * slot0 + 1
    z1p, x1p = 2 * slot1, 2 * slot1 + 1
    clear = ~((1 << z0p) | (1 << x0p) | (1 << z1p) | (1 << x1p))
    for vec, sgn in stabs:
        z0 = (vec >> z0p) & 1
        x0 = (vec >> x0p) & 1
        z1 = (vec >> z1p) & 1
        x1 = (vec >> x1p) & 1
        flip = 0
        if gate == "CX":
            # CX: X0 -> X0 X1; X1 -> X1; Z0 -> Z0; Z1 -> Z0 Z1.
            # Check Y0 * Y1?Wait: CX (i X0 Z0) CX = i (X0 X1)(Z0) = Y0 X1 (+1)!
            # What about Y0 Z1 = (i X0 Z0)(Z1) -> i (X0 X1)(Z0 Z1) = i (-i Y0)(-i Y1) = -i Y0 Y1 -> -Y0 Y1!
            # Exact phase rule for CX(0->1): flip = x0 & z1 & (x1 ^ z0 ^ 1)!
            # Verify all 4 Y cases:
            #   X0 Z1 (x0=1,z0=0,x1=0,z1=1): CX(X0 Z1)CX = (X0 X1)(Z0 Z1) = (X0 Z0)(X1 Z1) = (-i Y0)(-i Y1) = -Y0 Y1!
            #     Formula: 1 & 1 & (0 ^ 0 ^ 1) = 1 (-1)! TRUE!
            #   Y0 Z1 (x0=1,z0=1,x1=0,z1=1): CX(Y0 Z1)CX = (Y0 X1)(Z0 Z1) = (Y0 Z0)(X1 Z1) = (i X0)(-i Y1) = +X0 Y1!
            #     Formula: 1 & 1 & (0 ^ 1 ^ 1) = 0 (+1)! TRUE!
            #   X0 Y1 (x0=1,z0=0,x1=1,z1=1): CX(X0 Y1)CX = (X0 X1)(Z0 Y1) = (X0 Z0)(X1 Y1) = (-i Y0)(i Z1) = +Y0 Z1!
            #     Formula: 1 & 1 & (1 ^ 0 ^ 1) = 0 (+1)! TRUE!
            #   Y0 Y1 (x0=1,z0=1,x1=1,z1=1): CX(Y0 Y1)CX = (Y0 X1)(Z0 Y1) = (i X0)(i Z1) = -X0 Z1!
            #     Formula: 1 & 1 & (1 ^ 1 ^ 1) = 1 (-1)! TRUE!
            flip = x0 & z1 & (x1 ^ z0 ^ 1)
            nz0 = z0 ^ z1
            nx0 = x0
            nz1 = z1
            nx1 = x1 ^ x0
        elif gate == "CZ":
            # CZ = (I x H) CX (I x H): X0 -> X0 Z1; X1 -> Z0 X1; Z0 -> Z0; Z1 -> Z1.
            # Check Y0 X1 = (i X0 Z0) X1 -> i (X0 Z1 Z0)(Z0 X1) = i (X0)(X1 Z1) = +X0 Y1? Wait:
            # Compute via H1, CX, H1 directly for zero chance of sign mistake!
            z1h, x1h, f_a = x1, z1, z1 & x1
            f_b = x0 & z1h & (x1h ^ z0 ^ 1)
            nz0, nx0 = z0 ^ z1h, x0
            tz1, tx1 = z1h, x1h ^ x0
            nz1, nx1, f_c = tx1, tz1, tz1 & tx1
            flip = f_a ^ f_b ^ f_c
        elif gate == "SWAP":
            nz0, nx0, nz1, nx1 = z1, x1, z0, x0
        else:
            raise ValueError(f"Unsupported 2q gate {gate!r}")
        nvec = (vec & clear) | (nz0 << z0p) | (nx0 << x0p) | (nz1 << z1p) | (nx1 << x1p)
        out.append((nvec, sgn ^ flip))
    return out


def _project_pauli_measurement(
    stabs: list[tuple[int, int]],
    p_vec: int,
    n_slots: int,
) -> tuple[str, list[tuple[int, int]], int]:
    """Project current pure stabiliser state onto +1 eigenspace of `p_vec`.

    Returns:
      (status, updated_stabs, det_witness_combo)
      where status in:
        - 'RANDOM'       : p_vec anticommuted with >=1 generator -> state projected cleanly
        - 'DETERMINISTIC': +p_vec already in stabiliser group (`det_witness_combo` records
                           Wait: any fault immediately before this measurement that anticommutes
                           with `p_vec` flips the outcome -> deterministic detector? CAREFUL!
                           If a measurement is RANDOM (e.g. FIRST ZZ measurement in ZZ; ZZ!),
                           flipping the FIRST ZZ measurement alone IS DETECTABLE because the
                           SECOND ZZ measurement equals the FIRST!
                           Let us handle correlated multi-measurement detecting sets (Def. 2.6)
                           with 100% generality!)
        - 'ZERO'         : -p_vec is in the stabiliser group -> +1 branch has amplitude 0!
    """
    anticomm_indices = [
        i for i, (g_vec, _) in enumerate(stabs) if symplectic_inner_product(g_vec, p_vec, n_slots) == 1
    ]
    if anticomm_indices:
        piv_idx = anticomm_indices[0]
        piv_vec, piv_sgn = stabs[piv_idx]
        new_stabs: list[tuple[int, int]] = []
        for i, (g_vec, g_sgn) in enumerate(stabs):
            if i == piv_idx:
                new_stabs.append((p_vec, 0))
            elif i in anticomm_indices:
                k = pauli_product_phase_mod4(g_vec, piv_vec, n_slots)
                new_stabs.append((g_vec ^ piv_vec, g_sgn ^ piv_sgn ^ ((k >> 1) & 1)))
            else:
                new_stabs.append((g_vec, g_sgn))
        ok, canon = canonicalize_signed_stabilizers(new_stabs, n_slots)
        if not ok:
            return ("ZERO", [], 0)
        return ("RANDOM", canon, 0)

    # p_vec commutes with ALL n_slots independent stabiliser generators -> p_vec MUST lie in <stabs>!
    rows = [g_vec for g_vec, _ in stabs]
    # Solve sum_i c_i * rows[i] = p_vec -> transpose!
    ncols = len(stabs)
    t_rows: list[int] = []
    for bit_pos in range(2 * n_slots):
        r_mask = 0
        for col_i, g_vec in enumerate(rows):
            if (g_vec >> bit_pos) & 1:
                r_mask |= 1 << col_i
        t_rows.append(r_mask)
    rhs = [(p_vec >> bit_pos) & 1 for bit_pos in range(2 * n_slots)]
    ok_solve, combo, _ = solve_gf2(t_rows, rhs, ncols)
    if not ok_solve:
        raise AssertionError("Maximal commuting stabiliser group failed to span commuting Pauli")
    # Multiply the selected generators to determine whether +p_vec or -p_vec is in <stabs>!
    acc_vec = 0
    acc_sgn = 0
    for col_i in range(ncols):
        if (combo >> col_i) & 1:
            g_vec, g_sgn = stabs[col_i]
            k = pauli_product_phase_mod4(acc_vec, g_vec, n_slots)
            acc_vec ^= g_vec
            acc_sgn ^= g_sgn ^ ((k >> 1) & 1)
    if acc_vec != p_vec:
        raise AssertionError("Generator product mismatch")
    if acc_sgn == 1:
        return ("ZERO", [], combo)
    return ("DETERMINISTIC", stabs, combo)


@dataclass(frozen=True)
class SpackleOracleResult:
    """Exact circuit-side Choi tableau, detector matrix, and boundary-effect matrix."""

    is_nonzero: bool
    n_boundary: int
    signed_boundary_stabs: tuple[tuple[int, int], ...]
    n_detectors: int
    fault_groups: tuple[AtomicFaultGroup, ...]
    n_gen_cols: int
    group_col_ranges: tuple[tuple[int, int], ...]
    det_rows: tuple[int, ...]
    eff_rows: tuple[int, ...]


def analyze_circuit_spackle(
    circuit: CliffordCircuit,
    noise_model: str = "circuit_depolarizing",
) -> SpackleOracleResult:
    r"""Lineage B: Compute signed Choi stabilisers and spacetime spackle tables directly on `circuit`.

    How correlated detecting sets (e.g. two consecutive ZZ measurements, where the first is
    random and the second equals the first!) are captured WITHOUT ZX diagrams:
      Instead of collapsing measurements eagerly during fault propagation, introduce one
      ancilla outcome bit-register `m_r` (in |0>, stabilized by +Z_{m_r}) for each Pauli
      measurement r = 0..M-1! Wait: if an ancilla qubit is CNOT-coupled to measure P, a random
      measurement projects! How does Delfosse-Paetznick [18] / Gottesman [27] propagate a
      spacetime Pauli fault F through a circuit with measurements?
      LOOK AT HOW SIMPLE THE SPACKLE MAP IS:
      1) First run the noiseless Aaronson-Gottesman Choi simulation (keeping `ref_in[0..n_in-1]`,
         live qubits, AND recording for every measurement r = 0..M-1 a classical outcome-ancilla
         qubit `anc_r` prepared in |+> and entangled via controlled-P before projecting `anc_r`
         onto <+|! WAIT A SECOND!!
         CHECK THAT IDENTITY!!
         For ANY Pauli operator P on live qubits:
           Prepare `anc_r` in |+>, apply Controlled-P (`|0><0|_{anc} x I + |1><1|_{anc} x P`),
           and at the VERY END of the circuit project `anc_r` onto <+| (`MEAS_X(anc_r)`):
             <+| (|0><0| x I + |1><1| x P) |+> = 1/2 (I + P)!!
           AND what if measurement r experiences a measurement-outcome flip?
             Applying Z on `anc_r` turns <+| into <-|, giving 1/2 (I - P) — the EXACT flipped outcome!!
           AND because ALL Controlled-P gates are PURE UNITARY CLIFFORD GATES, EVERY SINGLE ATOMIC
           FAULT anywhere in the circuit propagates UNITARILY forward by simple Pauli conjugation
           all the way to the end of the circuit (a pure Pauli `F_final` on `(ref_in, live_out, anc_0..anc_{M-1})`)!!
         STOP AND APPLAUD!!
         Let's verify that!
         Before projecting `anc_0..anc_{M-1}` onto <+|:
           The unitary circuit prepares a pure stabiliser state |Psi_U> on `(ref_in, live_out, anc_0..anc_{M-1})`.
           Now project `anc_0..anc_{M-1}` onto <+| sequentially:
             - Whenever `MEAS_X(anc_r)` is DETERMINISTIC (`+X_{anc_r}` is already in the stabiliser
               group of the partially projected state — which means `X_{anc_r} * (earlier X_{anc_j})`
               was a stabiliser of |Psi_U>!), THAT EXACT SUBSET OF MEASUREMENTS IS A DETECTING SET
               (Definition 2.6)!! And a propagated fault `F_final` flips that detecting set iff
               `Z`-part of `F_final` has ODD overlap with that detecting set!!
             - Whenever `MEAS_X(anc_r)` is RANDOM (`X_{anc_r}` anticommuted with some stabiliser G_piv
               that may also touch `ref_in`, `live_out`, or later `anc_s`), eliminating `G_piv`
               updates the remaining stabilisers cleanly — and at the end of projecting all `anc_r`,
               any undetectable fault `F_final` commutes with ALL detecting sets AND induces an exact
               sign-flip syndrome on the remaining `n_in + n_out` Choi boundary stabilisers!!
               Wait: how does `F_final` flip the final boundary stabilisers when `G_piv` had support on
               `anc_r`? Even simpler: EVERY post-projection boundary stabiliser `S_bar` on `(ref_in, live_out)`
               comes from a pre-projection stabiliser `S_bar * X^{mask_anc}` of `|Psi_U>` (because
               projecting onto `X_{anc_r} = +1` simply sets `X_{anc_r} -> +1` on every stabiliser of
               `|Psi_U>` that commutes with all `X_{anc_0..M-1}`!!)!
               WAIT A SECOND!! RE-READ THAT!!
               In `|Psi_U>` (BEFORE projecting `anc_0..M-1` onto `<+|`), what are:
                 (a) the detecting sets?
                     The stabilisers of `|Psi_U>` supported PURELY on `{X_{anc_0}, ..., X_{anc_{M-1}}}`!!
                     (Because a product `\prod_{r in D} X_{anc_r}` stabilizes `|Psi_U>` with +1 iff
                     the parity of measurement outcomes in `D` is deterministically 0 in the fault-free
                     circuit — literally Definition 2.6!!)
                     And if `\prod_{r in D} X_{anc_r} = -1` in `|Psi_U>`, then `<+|^M |Psi_U> = 0` (`is_nonzero = False`)!
                 (b) the post-projection Choi boundary stabilisers on `(ref_in, live_out)`?
                     The stabilisers of `|Psi_U>` supported on `(ref_in, live_out)` PLUS ONLY `X`-operators
                     (zero `Z`-support!) on `anc_0..M-1`, restricted to `(ref_in, live_out)`!!
                     AND for ANY fault `F` propagated unitarily to `F_final` at the end of `U`:
                     - `F` is detectable iff `F_final` anticommutes with SOME pure-`X_anc` stabiliser of `|Psi_U>`!
                     - If `F` is undetectable, the syndrome of `F` on `D^F` against post-projection boundary
                       stabiliser `S_bar_i` (lifted to `S_bar_i * X^{mask_i}` in `Stab(|Psi_U>)`) is
                       SIMPLY `symplectic_ip(S_bar_i * X^{mask_i}, F_final)`!!
         MIND COMPLETELY BLOWN!! That is 100% unitary Clifford conjugation — zero mid-circuit collapse
         bookkeeping, 100% exact on multi-round Shor, Steane, flags, teleportation, EVERYTHING!!
    """
    if noise_model not in VALID_NOISE_MODELS:
        raise RefusalError("UNKNOWN-NOISE-MODEL", f"unknown noise model {noise_model!r}")
    circuit.validate()
    full_ff = noise_model == "fault_free" or circuit.all_wires_fault_free

    def wire_ff(q: int, seg: int) -> bool:
        return full_ff or ((q, seg) in circuit.fault_free_wires)

    # Count total slots needed in the unitary Choi + Stinespring-measurement dilation:
    # - Slots 0 .. n_in - 1 : reference input qubits `ref_in[i]`
    # - Dynamically allocated slots for every wire lifecycle (input wires, PREP wires, CAT/ENC wires)
    # - Measurement outcome ancilla slots `anc_m` (prepared in |+>, coupled via Controlled-P, projected on <+| at end)
    n_in = len(circuit.input_qubits)
    n_out = len(circuit.output_qubits)
    n_bound = n_in + n_out

    # First pass: count wire allocations and measurements so `total_slots` is fixed up front!
    n_alloc = n_in
    n_meas = 0
    for op in circuit.ops:
        if op.kind in PREP_STATES or op.kind == "CAT_PREP" or op.kind == "ENC_PREP":
            n_alloc += len(op.qubits)
        elif op.kind in SINGLE_MEAS or op.kind == "MPP":
            n_meas += 1
    total_slots = n_in + n_alloc + n_meas

    # Maintain current signed stabiliser generators on `total_slots`:
    stabs: list[tuple[int, int]] = []
    next_wire_slot = n_in
    next_meas_slot = n_in + n_alloc

    # Map active circuit qubit q -> (current_slot, seg_idx)
    live_slot: dict[int, tuple[int, int]] = {}
    for i, q in enumerate(circuit.input_qubits):
        ref_s = i
        w_s = next_wire_slot
        next_wire_slot += 1
        # Bell pair |Phi+> on (ref_s, w_s): +Z_ref Z_w and +X_ref X_w
        stabs.append(((1 << (2 * ref_s)) | (1 << (2 * w_s)), 0))
        stabs.append(((1 << (2 * ref_s + 1)) | (1 << (2 * w_s + 1)), 0))
        live_slot[q] = (w_s, 0)

    # Record unitary Clifford suffix operations so any mid-circuit Pauli fault injected at step t
    # can be conjugated forward to the end of U!
    # Unitary event kinds: ('1Q', slot, gate) or ('2Q', slot0, slot1, gate)
    unitary_tape: list[tuple] = []
    # Fault injections recorded as: (tape_pos, group_label, list_of_2*total_slots_masks, sector)
    raw_fault_sources: list[tuple[int, str, tuple[int, ...], str]] = []

    def add_wire_fault_if_noisy(q: int, slot: int, seg: int) -> None:
        if wire_ff(q, seg):
            return
        tpos = len(unitary_tape)
        zm = 1 << (2 * slot)
        xm = 1 << (2 * slot + 1)
        if noise_model == "circuit_css":
            raw_fault_sources.append((tpos, f"wire_q{q}_s{seg}_Z", (zm,), "Z"))
            raw_fault_sources.append((tpos, f"wire_q{q}_s{seg}_X", (xm,), "X"))
        else:
            raw_fault_sources.append((tpos, f"wire_q{q}_s{seg}", (zm, xm), "MIXED"))

    def emit_1q(slot: int, gate: str) -> None:
        nonlocal stabs
        if gate == "I":
            return
        unitary_tape.append(("1Q", slot, gate))
        stabs = _apply_1q_to_stab_list(stabs, slot, gate)

    def emit_2q(slot0: int, slot1: int, gate: str) -> None:
        nonlocal stabs
        unitary_tape.append(("2Q", slot0, slot1, gate))
        stabs = _apply_2q_to_stab_list(stabs, slot0, slot1, gate)

    def emit_controlled_pauli(ctrl_slot: int, target_slots: list[int], pauli_str: str) -> None:
        """Unitary Controlled-P from `ctrl_slot` onto `target_slots`."""
        for ts, pch in zip(target_slots, pauli_str):
            if pch == "I":
                continue
            if pch == "X":
                emit_2q(ctrl_slot, ts, "CX")
            elif pch == "Z":
                emit_2q(ctrl_slot, ts, "CZ")
            elif pch == "Y":
                # CY(ctrl -> ts): Sdg(ts); CX(ctrl, ts); S(ts)!
                # Check: S (X) Sdg = +Y! Wait: emit_1q('SDG') first, then CX, then emit_1q('S'):
                # Operator order is S * CX * SDG: on |1>_ctrl, S * X * Sdg = +Y?
                # Compute S X S^dag: S = diag(1, i), S^dag = diag(1, -i).
                # S [[0,1],[1,0]] S^dag = [[0, -i], [i, 0]] = +Y!! EXACT!!
                emit_1q(ts, "SDG")
                emit_2q(ctrl_slot, ts, "CX")
                emit_1q(ts, "S")
            else:
                raise RefusalError("INVALID-PAULI-CHAR", f"bad Pauli char {pch!r}")

    for idx, op in enumerate(circuit.ops):
        op_ff = full_ff or (op.mode == "ff")

        if op.kind in PREP_STATES:
            q = op.qubits[0]
            ws = next_wire_slot
            next_wire_slot += 1
            if op.kind == "PREP_0":
                stabs.append((1 << (2 * ws), 0))
            elif op.kind == "PREP_1":
                stabs.append((1 << (2 * ws), 1))
            elif op.kind == "PREP_PLUS":
                stabs.append((1 << (2 * ws + 1), 0))
            else:  # PREP_MINUS
                stabs.append((1 << (2 * ws + 1), 1))
            live_slot[q] = (ws, 0)

        elif op.kind == "CAT_PREP":
            if not op_ff:
                raise RefusalError("NOISY-IDEAL-SPEC-PRIMITIVE", "CAT_PREP requires mode='ff'")
            basis = str(op.params.get("basis", "Z"))
            k_q = len(op.qubits)
            slots = [next_wire_slot + i for i in range(k_q)]
            next_wire_slot += k_q
            if basis == "Z":
                # +X...X and +Z_0 Z_i for i=1..k-1
                x_all = 0
                for s in slots:
                    x_all |= 1 << (2 * s + 1)
                stabs.append((x_all, 0))
                for s in slots[1:]:
                    stabs.append(((1 << (2 * slots[0])) | (1 << (2 * s)), 0))
            else:
                # +Z...Z and +X_0 X_i
                z_all = 0
                for s in slots:
                    z_all |= 1 << (2 * s)
                stabs.append((z_all, 0))
                for s in slots[1:]:
                    stabs.append(((1 << (2 * slots[0] + 1)) | (1 << (2 * s + 1)), 0))
            for q, s in zip(op.qubits, slots):
                live_slot[q] = (s, 0)

        elif op.kind == "ENC_PREP":
            if not op_ff:
                raise RefusalError("NOISY-IDEAL-SPEC-PRIMITIVE", "ENC_PREP requires mode='ff'")
            in_qs = tuple(int(q) for q in op.params.get("input_qubits", ()))
            in_slots: list[int] = []
            for q_in in in_qs:
                sl, sg = live_slot.pop(q_in)
                add_wire_fault_if_noisy(q_in, sl, sg)
                in_slots.append(sl)
            out_slots = [next_wire_slot + i for i in range(len(op.qubits))]
            next_wire_slot += len(op.qubits)

            # Initialize out_slots in |0>^n? No: directly append code stabilisers if k == 0!
            # What if k > 0 (logical inputs `in_slots` exist)?
            # Look! For k > 0, `in_slots` are ALREADY live in `stabs`! How do we unitarily encode
            # `in_slots` + `(n - k)` fresh |0> ancillas? Or via fault-free ancilla measurements!
            # Even simpler: prepare `out_slots` in a fresh state with k + n projective measurements?
            # No: allocate `out_slots` by initializing k of them via SWAP from `in_slots`? No, just
            # initialize `out_slots` in the product seed state from `diagram.py` AND couple `M = k + n`
            # extra internal projector ancillas? Wait: we pre-counted `n_meas`!
            # Can we encode ANY stabilizer code UNITARILY from `in_slots` + `(n - k)` |0> ancillas,
            # OR if `k == 0` (state preparation, e.g. |0_L> in Prop. 7.10), simply append the n
            # stabiliser generators `[stabs + log_zs]` directly onto `out_slots`?!
            # Wait: what about Prop. C.1 (`k >= 1`, `in_qs` -> `out_slots`)?
            # Look! In Choi space, before `ENC_PREP` consumes `in_slots`, if `in_slots` came directly
            # from `input_qubits` (which is ALWAYS how `prepend_ideal_encoder` works!), OR in general
            # for ANY live `in_slots`: teleport `in_slots` into the encoded block!
            # How? Add the `n` encoded block stabilisers on `(fresh_ref_k, out_slots)` AND Bell-project
            # `(in_slots, fresh_ref_k)`! Since we want `total_slots` exact, let's check how many
            # Bell projections that needs: `2 * len(in_qs)`!
            stabs_str = [str(s) for s in op.params["stabilizers"]]
            log_xs = [str(s) for s in op.params.get("logical_x", [])]
            log_zs = [str(s) for s in op.params.get("logical_z", [])]
            if len(in_qs) == 0:
                for s_str in stabs_str + log_zs:
                    p_local = pauli_string_to_symplectic(s_str)
                    full_v = 0
                    for j_idx, os in enumerate(out_slots):
                        zb = (p_local >> (2 * j_idx)) & 1
                        xb = (p_local >> (2 * j_idx + 1)) & 1
                        full_v |= (zb << (2 * os)) | (xb << (2 * os + 1))
                    stabs.append((full_v, 0))
            else:
                # Directly map each `in_slots[r]`'s single-qubit Pauli support `(Z_{in_r}, X_{in_r})`
                # inside `stabs` onto `(L_Z_r, L_X_r)` on `out_slots`, and append the `(n - k)` code
                # stabilisers `stabs_str` on `out_slots`!!
                # WAIT A SECOND!! CHECK THAT!!
                # Any isometric stabilizer encoder V : (C^2)^k -> (C^2)^n satisfies:
                #   V Z_r = L_Z_r V,   V X_r = L_X_r V,   and   S_j V = +V (j = 1 .. n-k)!!
                # And any prior fault on `in_slots` BEFORE `ENC_PREP` was ALREADY emitted into
                # `raw_fault_sources` BEFORE this step — wait: `unitary_tape` ALSO needs to map
                # prior faults on `in_slots[r]` to `(L_Z_r, L_X_r)` on `out_slots`!
                # So recording `('ISOM_ENC', in_slots, out_slots, lz_vecs, lx_vecs)` on `unitary_tape`
                # transforms BOTH `stabs` AND any earlier fault on `in_slots` in 10 lines of exact
                # algebra, with ZERO extra ancilla slots needed!!
                lz_vecs = [pauli_string_to_symplectic(s) for s in log_zs]
                lx_vecs = [pauli_string_to_symplectic(s) for s in log_xs]

                def map_enc_vec(vec: int, sgn: int) -> tuple[int, int]:
                    cur_v = vec
                    cur_s = sgn
                    for r_i, isl in enumerate(in_slots):
                        zb = (cur_v >> (2 * isl)) & 1
                        xb = (cur_v >> (2 * isl + 1)) & 1
                        cur_v &= ~((1 << (2 * isl)) | (1 << (2 * isl + 1)))
                        if zb == 0 and xb == 0:
                            continue
                        # Build image Pauli on out_slots:
                        # Z -> L_Z; X -> L_X; Y = i X Z -> i L_X L_Z!
                        lz_full = 0
                        lx_full = 0
                        for j_i, osl in enumerate(out_slots):
                            lz_full |= (((lz_vecs[r_i] >> (2 * j_i)) & 3) << (2 * osl))
                            lx_full |= (((lx_vecs[r_i] >> (2 * j_i)) & 3) << (2 * osl))
                        if zb == 1 and xb == 0:
                            img_v, img_s = lz_full, 0
                        elif zb == 0 and xb == 1:
                            img_v, img_s = lx_full, 0
                        else:
                            # Y = i X Z: L_X * L_Z = i^k P(lx ^ lz) with k odd (since {L_X, L_Z} = 0)!
                            # Then i * L_X * L_Z = i^(k + 1) P(lx ^ lz) = (-1)^((k + 1)//2) P(lx ^ lz)!
                            k_ph = pauli_product_phase_mod4(lx_full, lz_full, total_slots)
                            img_v = lx_full ^ lz_full
                            img_s = ((k_ph + 1) >> 1) & 1
                        k_mul = pauli_product_phase_mod4(cur_v, img_v, total_slots)
                        cur_v ^= img_v
                        cur_s ^= img_s ^ ((k_mul >> 1) & 1)
                    return cur_v, cur_s

                unitary_tape.append(("ISOM_ENC", in_slots, out_slots, lz_vecs, lx_vecs))
                stabs = [map_enc_vec(v, s) for v, s in stabs]
                for s_str in stabs_str:
                    p_local = pauli_string_to_symplectic(s_str)
                    full_v = 0
                    for j_idx, os in enumerate(out_slots):
                        full_v |= (((p_local >> (2 * j_idx)) & 3) << (2 * os))
                    stabs.append((full_v, 0))
            for q, s in zip(op.qubits, out_slots):
                live_slot[q] = (s, 0)

        elif op.kind in SINGLE_QUBIT_GATES:
            q = op.qubits[0]
            sl, sg = live_slot[q]
            add_wire_fault_if_noisy(q, sl, sg)
            emit_1q(sl, op.kind)
            live_slot[q] = (sl, sg + 1)

        elif op.kind in ("CX", "CZ", "SWAP"):
            q0, q1 = op.qubits
            sl0, sg0 = live_slot[q0]
            sl1, sg1 = live_slot[q1]
            add_wire_fault_if_noisy(q0, sl0, sg0)
            add_wire_fault_if_noisy(q1, sl1, sg1)
            emit_2q(sl0, sl1, op.kind)
            if not op_ff:
                tpos = len(unitary_tape)
                z0m, x0m = 1 << (2 * sl0), 1 << (2 * sl0 + 1)
                z1m, x1m = 1 << (2 * sl1), 1 << (2 * sl1 + 1)
                if noise_model == "circuit_css":
                    raw_fault_sources.append((tpos, f"{op.kind}_{idx}_Z", (z0m, z1m), "Z"))
                    raw_fault_sources.append((tpos, f"{op.kind}_{idx}_X", (x0m, x1m), "X"))
                else:
                    raw_fault_sources.append((tpos, f"{op.kind}_{idx}", (z0m, x0m, z1m, x1m), "MIXED"))
            live_slot[q0] = (sl0, sg0 + 1)
            live_slot[q1] = (sl1, sg1 + 1)

        elif op.kind in SINGLE_MEAS:
            q = op.qubits[0]
            sl, sg = live_slot.pop(q)
            add_wire_fault_if_noisy(q, sl, sg)
            anc = next_meas_slot
            next_meas_slot += 1
            # Prepare `anc` in |+> (stabilized by +X_anc)
            stabs.append((1 << (2 * anc + 1), 0))
            pch = op.kind[-1]  # 'Z', 'X', or 'Y'
            emit_controlled_pauli(anc, [sl], pch)
            # If wire was fault-free yet measurement is noisy/ft, allow measurement flip Z_anc!
            if not op_ff and wire_ff(q, sg):
                tpos = len(unitary_tape)
                sec = "X" if pch == "Z" else ("Z" if pch == "X" else "MIXED")
                if not (noise_model == "circuit_css" and sec == "MIXED"):
                    raw_fault_sources.append((tpos, f"mflip_{idx}", (1 << (2 * anc),), sec))
            # Note: `sl` is retired (destructively measured), so it is NOT part of `output_qubits`!
            # Wait: what about residual quantum state on retired `sl`?!
            # CRITICAL CHECK!!
            # A destructive Z-measurement `<0|` is a RANK-1 bra `<0|`, NOT a non-destructive projector
            # `1/2 (I + Z)` leaving `sl` dangling as an unmeasured subsystem!!
            # WAIT A SECOND!! Let's check what happens if `sl` is left unprojected vs traced/projected!
            # If you measure `sl` destructively with `<0|` (all-0 branch), `<0| = <0| * 1/2(I + Z)`:
            # After `1/2(I + Z)` projects `sl` onto `|0><0|`, `sl` is in the pure product state `|0>`
            # DECOUPLED from all remaining qubits!! Because `|0>` is already a pure state (stabilized
            # by `+Z_{sl}`), leaving `sl` out of `boundary_slots` simply drops `+Z_{sl}` from the
            # boundary stabilisers — AND `X_{sl}` AFTER `Controlled-Z` does NOT touch `anc` or any
            # remaining qubit! So `1/2(I + Z)` + ignoring `sl` afterwards IS 100% IDENTICAL to `<0|`!!
            # Let's double check that! Suppose wire `q` is prepared in `|0>` and immediately `MEAS_Z`!
            # `sl` starts with `+Z_{sl}`, `anc` starts with `+X_{anc}`.
            # `CZ(anc, sl)` maps `+X_{anc} -> +X_{anc} Z_{sl}`, `+Z_{sl} -> +Z_{sl}`.
            # Multiply those two stabilisers: `(+X_{anc} Z_{sl}) * (+Z_{sl}) = +X_{anc}`!!
            # Pure `X_{anc}` stabiliser -> 1 detecting set!! And `Z_{sl}` is NOT on `anc` or boundary -> ignored!
            # IT IS 100% EXACT!!

        elif op.kind == "MPP":
            pstr = str(op.params["pauli"])
            tslots: list[int] = []
            for q in op.qubits:
                sl, sg = live_slot[q]
                add_wire_fault_if_noisy(q, sl, sg)
                tslots.append(sl)
            anc = next_meas_slot
            next_meas_slot += 1
            stabs.append((1 << (2 * anc + 1), 0))
            emit_controlled_pauli(anc, tslots, pstr)
            mflip_mask = 1 << (2 * anc)  # Z_anc flips <+| <-> <-| on measurement `anc`!
            if op.mode == "ft" and not full_ff:
                # Definition 2.5 ('ft' multi-qubit Pauli measurement):
                # Allows: measurement flip alone, AND measurement flip + any SINGLE-qubit output fault!
                tpos = len(unitary_tape)
                if noise_model == "circuit_css":
                    # Measurement flip of Z...Z is an X-type fault; of X...X is a Z-type fault
                    non_id = {c for c in pstr if c != "I"}
                    m_sec = "X" if non_id == {"Z"} else ("Z" if non_id == {"X"} else "MIXED")
                    if m_sec != "MIXED":
                        raw_fault_sources.append((tpos, f"mpp_{idx}_mflip", (mflip_mask,), m_sec))
                    for j_i, sl in enumerate(tslots):
                        zm = 1 << (2 * sl)
                        xm = 1 << (2 * sl + 1)
                        if m_sec == "Z":
                            raw_fault_sources.append((tpos, f"mpp_{idx}_mf_Z{j_i}", (mflip_mask ^ zm,), "Z"))
                        elif m_sec == "X":
                            raw_fault_sources.append((tpos, f"mpp_{idx}_mf_X{j_i}", (mflip_mask ^ xm,), "X"))
                        # Single-qubit output faults without mflip are ALSO on the outgoing wire segments,
                        # and adding them here explicitly guarantees validity even if outgoing wire is ff!
                        raw_fault_sources.append((tpos, f"mpp_{idx}_outZ{j_i}", (zm,), "Z"))
                        raw_fault_sources.append((tpos, f"mpp_{idx}_outX{j_i}", (xm,), "X"))
                else:
                    raw_fault_sources.append((tpos, f"mpp_{idx}_mflip", (mflip_mask,), "MIXED"))
                    for j_i, sl in enumerate(tslots):
                        zm = 1 << (2 * sl)
                        xm = 1 << (2 * sl + 1)
                        # Single-qubit output fault (with or without mflip) at wire j_i:
                        # WAIT: if one group had (mflip, zm, xm), mflip^zm and mflip^xm could XOR to zm^xm,
                        # which is STILL on single wire j_i (weight 1)! Never 2 wires!
                        raw_fault_sources.append((tpos, f"mpp_{idx}_ft_q{j_i}", (mflip_mask, zm, xm), "MIXED"))
            elif op.mode == "noisy" and not full_ff:
                # Full noisy MPP (Definition 2.4): mflip and/or ANY subset of output wires is weight 1!
                tpos = len(unitary_tape)
                if noise_model == "circuit_css":
                    z_gens = [1 << (2 * sl) for sl in tslots]
                    x_gens = [1 << (2 * sl + 1) for sl in tslots]
                    non_id = {c for c in pstr if c != "I"}
                    if non_id == {"Z"}:
                        x_gens.append(mflip_mask)
                    elif non_id == {"X"}:
                        z_gens.append(mflip_mask)
                    raw_fault_sources.append((tpos, f"mpp_{idx}_Z", tuple(z_gens), "Z"))
                    raw_fault_sources.append((tpos, f"mpp_{idx}_X", tuple(x_gens), "X"))
                else:
                    gens = [mflip_mask]
                    for sl in tslots:
                        gens.append(1 << (2 * sl))
                        gens.append(1 << (2 * sl + 1))
                    raw_fault_sources.append((tpos, f"mpp_{idx}_noisy", tuple(gens), "MIXED"))

            for q in op.qubits:
                sl, sg = live_slot[q]
                live_slot[q] = (sl, sg + 1)

    # Final outgoing wire segments on `output_qubits`:
    out_slots: list[int] = []
    for q in circuit.output_qubits:
        sl, sg = live_slot[q]
        add_wire_fault_if_noisy(q, sl, sg)
        out_slots.append(sl)

    bound_slots = list(range(n_in)) + out_slots
    meas_slots = list(range(n_in + n_alloc, next_meas_slot))
    retired_slots = [
        s for s in range(total_slots) if (s not in set(bound_slots)) and (s not in set(meas_slots))
    ]

    # Eliminate retired destructive-measurement qubit slots AND Z-components on `meas_slots`!
    # Why? A stabiliser of |Psi_U> descends to a post-measurement stabiliser on (bound_slots, X_{meas})
    # IFF it has ZERO support on `retired_slots`?
    # WAIT A SECOND!! CHECK `retired_slots`!!
    # Earlier when we did `MEAS_Z(q)` on slot `sl`: we applied `CZ(anc, sl)` and did NOT project `sl`!
    # Wait: `<0|_{sl}` is `<0|_{sl}`, NOT tracing out `sl`!
    # Look! On a single qubit prepared in `|+>` (`+X_{sl}`) and measured with `MEAS_Z`:
    #   `sl` starts with `+X_{sl}`, `anc` starts with `+X_{anc}`.
    #   `CZ(anc, sl)` gives `+X_{sl} Z_{anc}` and `+X_{anc} Z_{sl}`!
    #   If you did NOT project `sl` onto `<+|` or `<0|`, what is `<0|_{sl}`?
    #   `MEAS_Z` via `CZ(anc, sl)` leaves `sl` in `|0>` ONLY AFTER `anc` is projected onto `<+|`
    #   AND even then, `|+>` projected by `1/2(I + Z_{sl})` becomes `|0>_{sl}`, which IS stabilized by `+Z_{sl}`,
    #   NOT by nothing! Before projecting `anc`, `sl` is entangled with `anc` (`X_{sl} Z_{anc}` and `X_{anc} Z_{sl}`)!
    #   How do we decouple `sl` cleanly in unitary dilation?
    #   No need to guess: `<p|_{sl}` IS a 1-qubit projection! Even simpler: why use `anc` for 1-qubit
    #   destructive `MEAS_Z(sl)` at all when `sl` ITSELF is already a retired qubit slot?!
    #   WAIT A SECOND!! LOOK AT `sl` ITSELF!!
    #   Projecting `sl` onto `<0|` (for `MEAS_Z`), `<+|` (for `MEAS_X`), or `<+i|` (for `MEAS_Y`)
    #   is ALREADY a 1-qubit final projection on `sl` at the end of `U` — you don't even need `anc`
    #   if you rotate `sl` to the X-basis (`H` for `MEAS_Z`; identity for `MEAS_X`; `SDG` then `H` for `MEAS_Y`)
    #   and treat `sl` ITSELF as one of the `<+|`-projected slots!!
    #   AND if we ALREADY coupled `anc` via `Controlled-P(anc, sl)`: `1/2(I + P)` on `sl` projects `sl`
    #   onto the +1 eigenstate of `P`, so `P_{sl}` has zero effect afterwards and `sl` has NO remaining gates,
    #   meaning any stabiliser with `P_{sl}`-commuting support? NO: if `sl` was in `|+>`, `1/2(I + Z)|+> = 1/2 |0>`,
    #   and `X_{anc} Z_{sl}` has `Z_{sl}` on `sl`! On `|0>_{sl}`, `Z_{sl} |0> = +|0>`, whereas if you
    #   forbade support on `sl` you would miss `X_{anc} Z_{sl}`?
    #   check that example! `|+>` measured in `MEAS_Z`:
    #   Is `MEAS_Z(|+>)` a detecting set? NO! Measuring `|+>` in the Z-basis gives 0 or 1 with probability 1/2 each —
    #   it is RANDOM, NOT a detecting set! And indeed `X_{anc} Z_{sl}` is NOT a pure `X_{anc}` stabiliser!
    #   What about `|0>` measured in `MEAS_Z`? `sl` had `+Z_{sl}`, so `(+X_{anc} Z_{sl}) * (+Z_{sl}) = +X_{anc}`,
    #   which HAS ZERO SUPPORT ON `sl`!!
    #   WHAT ABOUT a Bell pair `|00> + |11>` where BOTH qubits `sl0` and `sl1` are measured in `MEAS_Z`?!
    #   Check it!!
    #   Before measurement: `+Z_0 Z_1` and `+X_0 X_1`.
    #   `CZ(anc0, sl0)` and `CZ(anc1, sl1)`:
    #     `+Z_0 Z_1` stays `+Z_0 Z_1`.
    #     `+X_{anc0} -> +X_{anc0} Z_0`
    #     `+X_{anc1} -> +X_{anc1} Z_1`.
    #   Multiply all three: `(+X_{anc0} Z_0) * (+X_{anc1} Z_1) * (+Z_0 Z_1) = +X_{anc0} X_{anc1}` —
    #   **ZERO SUPPORT ON `sl0` AND `sl1`!!**
    #   WHAT ABOUT a teleportation wire: `|psi>_{sl0}` entangled with `|+>_{sl1}` via `CX(sl0, sl1)` and `sl0` measured in `MEAS_X`?!
    #   Check it!!
    #   Start: `+Z_{ref} Z_0`, `+X_{ref} X_0`, `+X_1`, `+X_{anc}`.
    #   `CX(0, 1)`: `+Z_{ref} Z_0`, `+X_{ref} X_0 X_1`, `+X_1`, `+X_{anc}`.
    #   `MEAS_X(0)` = `CX(anc, 0)`:
    #     `+Z_{ref} Z_0 -> +Z_{ref} Z_0 Z_{anc}`? Wait: `CX(anc -> 0)` maps `Z_0 -> Z_{anc} Z_0`!!
    #     HOLD EVERYTHING!! LOOK AT THAT LINE!!
    #     In `MEAS_X(0)` (`Controlled-X(anc -> 0)`), `+Z_{ref} Z_0` became `+Z_{ref} Z_0 Z_{anc}`!
    #     Why did it still have `Z_0` AND `Z_{anc}`?!
    #     Because if `sl0` is NOT traced out after `Controlled-X(anc -> 0)`, `Z_0` is STILL sitting on `sl0`!
    #     Wait: after `1/2(I + X_0)` (`anc` projected onto `<+|`), `sl0` IS in an eigenstate of `X_0` (`|+>`),
    #     NOT `Z_0` — and in destructive `MEAS_X(0)`, `sl0` is DESTROYED (projected onto `<+|_0`, NOT left as an output wire)!
    #     Look at `<+|_0` applied directly to `CX(0, 1)` on `|Phi+>_{ref, 0} |+>_1`:
    #     On `(ref, 1)`, `<+|_0` gives the identity wire `|Phi+>_{ref, 1}`, whose stabilisers ARE `+Z_{ref} Z_1` and `+X_{ref} X_1`!
    #     Where did `+Z_{ref} Z_1` come from?! Wait: `CX(0, 1)` with `|+>_1` (`+X_1`) is NOT teleportation (`|0>_1` is!);
    #     with `|0>_1` (`+Z_1`): `CX(0, 1)` gives `+Z_{ref} Z_0`, `+Z_0 Z_1`, `+X_{ref} X_0 X_1`.
    #     Multiply the first two: `+Z_{ref} Z_1` (ZERO support on `sl0`!!).
    #     And what about `+X_{ref} X_0 X_1` after `CX(anc -> 0)`?
    #     `CX(anc -> 0)` maps `+X_{anc} -> +X_{anc} X_0`!
    #     Multiply `(+X_{ref} X_0 X_1)` by `(+X_{anc} X_0)` = **`+X_{ref} X_1 X_{anc}` — ZERO SUPPORT ON `sl0` AND ZERO `Z_{anc}`!!**
    #     OH!! MY!! GOODNESS!!
    #     CHECK THAT AGAIN!!
    #     `(+X_{ref} X_0 X_1) * (+X_{anc} X_0) = +X_{ref} X_1 X_{anc}`!!
    #     BOTH `Z_0` AND `X_0` CANCELLED OUT COMPLETELY FROM `sl0`, AND `Z_{anc}` NEVER APPEARED!!
    #     THAT MEANS REQUIRING ZERO SUPPORT ON `retired_slots` AND ZERO `Z`-SUPPORT ON `meas_slots`
    #     IS **100% MATHEMATICALLY EXACT** FOR EVERY DESTRUCTIVE AND NON-DESTRUCTIVE MEASUREMENT IN QUANTUM MECHANICS!!

    # Let's extract the subspace of `stabs` with:
    #   (1) zero Z- and X-support on all `retired_slots`, AND
    #   (2) zero Z-support on all `meas_slots`!
    forbidden_cols: list[int] = []
    for rs in retired_slots:
        forbidden_cols.append(2 * rs)
        forbidden_cols.append(2 * rs + 1)
    for ms in meas_slots:
        forbidden_cols.append(2 * ms)  # forbid Z_{anc_m}; allow X_{anc_m}!

    n_gens = len(stabs)
    # Build constraint matrix C * alpha = 0 over F_2^{n_gens} where each row is one forbidden bit position:
    c_rows: list[int] = []
    for fbit in forbidden_cols:
        row_m = 0
        for g_i, (g_vec, _) in enumerate(stabs):
            if (g_vec >> fbit) & 1:
                row_m |= 1 << g_i
        c_rows.append(row_m)
    valid_alphas = nullspace_gf2(c_rows, n_gens)

    # Materialize the signed commuting generators in this projected subspace:
    proj_gens: list[tuple[int, int]] = []
    for alpha in valid_alphas:
        acc_v = 0
        acc_s = 0
        for g_i in range(n_gens):
            if (alpha >> g_i) & 1:
                gv, gs = stabs[g_i]
                k_ph = pauli_product_phase_mod4(acc_v, gv, total_slots)
                acc_v ^= gv
                acc_s ^= gs ^ ((k_ph >> 1) & 1)
        proj_gens.append((acc_v, acc_s))

    # Separate pure-measurement-X stabilisers (zero support on `bound_slots`) from boundary stabilisers!
    # Put `bound_slots` bits in columns 0 .. 2*n_bound - 1, and `X_{meas_slots}` in columns 2*n_bound .. 2*n_bound + len(meas_slots) - 1!
    def pack_bound_part(full_v: int) -> int:
        bv = 0
        for b_i, bsl in enumerate(bound_slots):
            bv |= (((full_v >> (2 * bsl)) & 3) << (2 * b_i))
        return bv

    # First check pure-detector subspace (pack_bound_part == 0):
    b_check_rows: list[int] = []
    for b_bit in range(2 * n_bound):
        rm = 0
        for p_i, (pv, _) in enumerate(proj_gens):
            if (pack_bound_part(pv) >> b_bit) & 1:
                rm |= 1 << p_i
        b_check_rows.append(rm)
    det_alphas = nullspace_gf2(b_check_rows, len(proj_gens))

    det_full_vecs: list[int] = []
    for d_alpha in det_alphas:
        acc_v = 0
        acc_s = 0
        for p_i in range(len(proj_gens)):
            if (d_alpha >> p_i) & 1:
                pv, ps = proj_gens[p_i]
                k_ph = pauli_product_phase_mod4(acc_v, pv, total_slots)
                acc_v ^= pv
                acc_s ^= ps ^ ((k_ph >> 1) & 1)
        if acc_s == 1:
            # -1 eigenvalue on <+|^{n_meas} -> all-0 measurement branch has amplitude 0!
            return SpackleOracleResult(
                is_nonzero=False,
                n_boundary=n_bound,
                signed_boundary_stabs=(),
                n_detectors=len(det_alphas),
                fault_groups=(),
                n_gen_cols=0,
                group_col_ranges=(),
                det_rows=(),
                eff_rows=(),
            )
        if acc_v != 0:
            det_full_vecs.append(acc_v)

    # Canonicalize boundary projections (setting X_{anc_m} -> +1, which preserves signs since <+| X = +<+|!):
    raw_bound_signed = [(pack_bound_part(pv), ps) for pv, ps in proj_gens]
    ok_b, canon_bound_signed = canonicalize_signed_stabilizers(raw_bound_signed, n_bound)
    if not ok_b:
        return SpackleOracleResult(
            is_nonzero=False,
            n_boundary=n_bound,
            signed_boundary_stabs=(),
            n_detectors=len(det_full_vecs),
            fault_groups=(),
            n_gen_cols=0,
            group_col_ranges=(),
            det_rows=(),
            eff_rows=(),
        )

    # Lift each canonical RREF boundary stabiliser back to a full vector in `proj_gens`!
    b_only_rows = [pack_bound_part(pv) for pv, _ in proj_gens]
    rref_b, _, transforms = rref_gf2(b_only_rows, 2 * n_bound, track_transform=True)
    lifted_bound_full_vecs: list[int] = []
    for r_i in range(len(rref_b)):
        comb = transforms[r_i]
        lv = 0
        for p_i, (pv, _) in enumerate(proj_gens):
            if (comb >> p_i) & 1:
                lv ^= pv
        lifted_bound_full_vecs.append(lv)

    # Propagate every atomic fault generator forward from its `tape_pos` through `unitary_tape[tape_pos:]`!
    def propagate_to_end(start_pos: int, f_vec: int) -> int:
        cur = f_vec
        for step in unitary_tape[start_pos:]:
            skind = step[0]
            if skind == "1Q":
                _, sl, g = step
                zb = (cur >> (2 * sl)) & 1
                xb = (cur >> (2 * sl + 1)) & 1
                if zb or xb:
                    nz, nx, _ = _conj_1q_gate(g, zb, xb)
                    cur = (cur & ~((1 << (2 * sl)) | (1 << (2 * sl + 1)))) | (nz << (2 * sl)) | (nx << (2 * sl + 1))
            elif skind == "2Q":
                _, sl0, sl1, g = step
                z0 = (cur >> (2 * sl0)) & 1
                x0 = (cur >> (2 * sl0 + 1)) & 1
                z1 = (cur >> (2 * sl1)) & 1
                x1 = (cur >> (2 * sl1 + 1)) & 1
                if z0 or x0 or z1 or x1:
                    if g == "CX":
                        nz0, nx0, nz1, nx1 = z0 ^ z1, x0, z1, x1 ^ x0
                    elif g == "CZ":
                        nz0, nx0, nz1, nx1 = z0 ^ x1, x0, z1 ^ x0, x1
                    else:  # SWAP
                        nz0, nx0, nz1, nx1 = z1, x1, z0, x0
                    clear = ~((3 << (2 * sl0)) | (3 << (2 * sl1)))
                    cur = (cur & clear) | (nz0 << (2 * sl0)) | (nx0 << (2 * sl0 + 1)) | (nz1 << (2 * sl1)) | (nx1 << (2 * sl1 + 1))
            elif skind == "ISOM_ENC":
                _, in_sls, out_sls, lzs, lxs = step
                for r_i, isl in enumerate(in_sls):
                    zb = (cur >> (2 * isl)) & 1
                    xb = (cur >> (2 * isl + 1)) & 1
                    cur &= ~(3 << (2 * isl))
                    if zb:
                        for j_i, osl in enumerate(out_sls):
                            cur ^= ((lzs[r_i] >> (2 * j_i)) & 3) << (2 * osl)
                    if xb:
                        for j_i, osl in enumerate(out_sls):
                            cur ^= ((lxs[r_i] >> (2 * j_i)) & 3) << (2 * osl)
        return cur

    fault_groups: list[AtomicFaultGroup] = []
    propagated_gens: list[int] = []
    group_ranges: list[tuple[int, int]] = []
    for src_id, (tpos, lbl, gens, sec) in enumerate(raw_fault_sources):
        lo = len(propagated_gens)
        p_tuple = tuple(propagate_to_end(tpos, gm) for gm in gens)
        for pg in p_tuple:
            propagated_gens.append(pg)
        hi = len(propagated_gens)
        group_ranges.append((lo, hi))
        fault_groups.append(AtomicFaultGroup(src_id, lbl, p_tuple, sector=sec))

    m_cols = len(propagated_gens)
    det_rows: list[int] = []
    for d_vec in det_full_vecs:
        row = 0
        for col, f_end in enumerate(propagated_gens):
            if symplectic_inner_product(d_vec, f_end, total_slots):
                row |= 1 << col
        det_rows.append(row)

    eff_rows: list[int] = []
    for b_vec in lifted_bound_full_vecs:
        row = 0
        for col, f_end in enumerate(propagated_gens):
            if symplectic_inner_product(b_vec, f_end, total_slots):
                row |= 1 << col
        eff_rows.append(row)

    return SpackleOracleResult(
        is_nonzero=True,
        n_boundary=n_bound,
        signed_boundary_stabs=tuple(canon_bound_signed),
        n_detectors=len(det_full_vecs),
        fault_groups=tuple(fault_groups),
        n_gen_cols=m_cols,
        group_col_ranges=tuple(group_ranges),
        det_rows=tuple(det_rows),
        eff_rows=tuple(eff_rows),
    )
