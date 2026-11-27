"""Zero-Float Outward-Rounded Stochastic Gadget Contracts (`ftqc_zx-contract-v1`).

Bridges adversarial w-fault equivalence (Rodatz, Poor, and Kissinger, arXiv:2506.17181v4)
to rigorous stochastic FTQC space-time ledgers (`ftqc_masking_v2c_certified_carbon` and
`ftqc_masking_v3a_fable` Filter Calculus):
  1) Packs any physical `CliffordCircuit` into a discrete ASAP tick schedule, counting
     concurrent active qubits, idle wire locations, gate counts, total locations `L`,
     and exact wire-tick volume `W`.
  2) Exhaustively audits all single-location faults under `circuit_depolarizing` or
     `circuit_css`, computing:
       - `v2c` per-coordinate active single-fault location counts `(v_X[c], v_Z[c])` for
         cone-reset `S_ec` coin envelopes, AND
       - `v3a` Filter-Calculus partition (`benign_zero`, `local_weight1`, `flagged_hook`,
         `bad_unheralded`) plus flag-conditioned decoder verification when inline.
  3) Evaluates exact all-orders fault-count polynomials in `fractions.Fraction` and emits
     outward-rounded `decimal.Decimal` string intervals (`ROUND_FLOOR` for lower bounds,
     `ROUND_CEILING` for upper bounds) with runtime containment verification!

Purity invariant: ZERO `float` literals or `float` builtin references appear in this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
from typing import Any

from .decoder_synth import FlagConditionedDecoder, _expand_rowspace, _min_coset_rep, _syndrome_bits, synthesize_flag_decoder
from .diagram import (
    CircuitOp,
    CliffordCircuit,
    PREP_STATES,
    RefusalError,
    SINGLE_MEAS,
    SINGLE_QUBIT_GATES,
    TWO_QUBIT_GATES,
    circuit_to_zx_fe,
)
from .equivalence import CSSCodeSpec
from .gf2 import rref_gf2, solve_gf2, swap_zx_bits
from .pauliweb import analyze_pauli_webs, compile_outcome_map

WORK_PREC = 100
OUT_PREC = 30

_CTX_LO = Context(prec=WORK_PREC, rounding=ROUND_FLOOR)
_CTX_HI = Context(prec=WORK_PREC, rounding=ROUND_CEILING)
_OUT_LO = Context(prec=OUT_PREC, rounding=ROUND_FLOOR)
_OUT_HI = Context(prec=OUT_PREC, rounding=ROUND_CEILING)


def parse_exact_rate(text: str | Fraction) -> Fraction:
    """Parse a rational or scientific-notation rate string (e.g. '1e-4', '1/10000', '0.0001') into Fraction."""
    if isinstance(text, Fraction):
        frac = text
    elif isinstance(text, str):
        s = text.strip()
        if not s:
            raise RefusalError("INVALID-RATE-STRING", "empty error rate string")
        try:
            if "/" in s:
                num_s, den_s = s.split("/", 1)
                frac = Fraction(int(num_s.strip()), int(den_s.strip()))
            else:
                # Fraction(Decimal_str) in Python stdlib parses '1e-4' and '0.0001' with 0 floats!
                frac = Fraction(s)
        except (ValueError, ZeroDivisionError) as exc:
            raise RefusalError("INVALID-RATE-STRING", f"cannot parse exact rate {text!r}: {exc}") from exc
    else:
        raise RefusalError("NON-EXACT-RATE-REFUSED", f"rate must be str or Fraction, got {type(text).__name__}")

    if frac < 0 or frac >= 1:
        raise RefusalError("RATE-OUT-OF-RANGE", f"physical fault rate must satisfy 0 <= p < 1, got {frac}")
    return frac


def fraction_to_outward_interval_strs(val: Fraction) -> tuple[str, str]:
    """Convert exact nonnegative `Fraction` to outward-rounded `(lo_str, hi_str)` Decimal strings."""
    if val < 0:
        raise ValueError(f"Expected nonnegative Fraction, got {val}")
    num_d = Decimal(val.numerator)
    den_d = Decimal(val.denominator)
    lo = _OUT_LO.plus(_CTX_LO.divide(num_d, den_d))
    hi = _OUT_HI.plus(_CTX_HI.divide(num_d, den_d))
    # Runtime containment assertion: exact Fraction MUST lie inside [Fraction(lo), Fraction(hi)]!
    if not (Fraction(str(lo)) <= val <= Fraction(str(hi))):
        raise AssertionError(f"Outward interval containment failed for {val}: [{lo}, {hi}]")
    return (format(lo, "f") if lo == 0 else str(lo), format(hi, "f") if hi == 0 else str(hi))


@dataclass(frozen=True)
class ScheduledLocation:
    """Single physical spacetime location `(tick, op)` in a packed circuit schedule."""

    loc_id: int
    tick: int
    kind: str                  # op.kind or 'IDLE'
    qubits: tuple[int, ...]
    label: str
    is_idle: bool
    meas_index: int = -1       # 0..n_meas-1 when `kind in SINGLE_MEAS`, else -1


@dataclass(frozen=True)
class TickSchedule:
    """Greedy ASAP discrete tick schedule and exact spacetime volume accounting."""

    n_qubits_total: int
    peak_concurrent_qubits: int
    peak_concurrent_ancillas: int
    ticks: int
    cnot_count: int
    meas_count: int
    prep_count: int
    single_gate_count: int
    idle_location_count: int
    total_locations: int
    wire_tick_volume: int
    locations: tuple[ScheduledLocation, ...]
    tick_ops: tuple[tuple[ScheduledLocation, ...], ...]


def schedule_circuit_ticks(circuit: CliffordCircuit) -> TickSchedule:
    """Pack `circuit` into greedy ASAP ticks and insert explicit `IDLE` locations on live wires."""
    circuit.validate()
    for op in circuit.ops:
        if op.kind in ("CAT_PREP", "ENC_PREP", "MPP") or op.mode != "noisy":
            raise RefusalError(
                "IDEALISED-OP-IN-CONTRACT",
                f"stochastic hardware contracts require physical noisy primitives; got {op.kind} (mode={op.mode!r})",
            )

    # Compute earliest tick per qubit; input qubits are live starting at tick 0!
    next_free_tick: dict[int, int] = {q: 0 for q in circuit.input_qubits}
    # Track segment lifetimes `(birth_tick, death_tick_inclusive)` for every qubit usage!
    placed_ops: list[tuple[int, int, CircuitOp]] = []  # (tick, orig_idx, op)

    for idx, op in enumerate(circuit.ops):
        if op.kind in PREP_STATES:
            q = op.qubits[0]
            t_use = next_free_tick.get(q, 0)
            placed_ops.append((t_use, idx, op))
            next_free_tick[q] = t_use + 1
        elif op.kind in SINGLE_MEAS:
            q = op.qubits[0]
            t_use = next_free_tick[q]
            placed_ops.append((t_use, idx, op))
            next_free_tick[q] = t_use + 1
        else:
            t_use = max(next_free_tick[q] for q in op.qubits)
            placed_ops.append((t_use, idx, op))
            for q in op.qubits:
                next_free_tick[q] = t_use + 1

    # Just-in-time `PREP_*` tightening: advance each `PREP_*` to `t_next_on_qubit - 1` so freshly
    # prepared ancillas do not sit idling before their first gate!
    per_qubit_timeline: dict[int, list[int]] = {}
    for p_pos, (_, _, op) in enumerate(placed_ops):
        for q in op.qubits:
            per_qubit_timeline.setdefault(q, []).append(p_pos)
    for q, p_indices in per_qubit_timeline.items():
        for k_i in range(len(p_indices) - 1):
            pos_cur = p_indices[k_i]
            pos_nxt = p_indices[k_i + 1]
            t_cur, orig_i, op_cur = placed_ops[pos_cur]
            t_nxt, _, _ = placed_ops[pos_nxt]
            if op_cur.kind in PREP_STATES and t_nxt - 1 > t_cur:
                placed_ops[pos_cur] = (t_nxt - 1, orig_i, op_cur)

    total_ticks = max(next_free_tick.values(), default=0)
    if total_ticks == 0 and circuit.output_qubits:
        total_ticks = 1

    # Reconstruct live intervals on each qubit across `0 .. total_ticks - 1`:
    #   - Input qubits start live at tick 0.
    #   - Prepared qubits become live at their `PREP_*` tick, and cease being live after `MEAS_*` tick!
    #   - Output qubits remain live through `total_ticks - 1`.
    # Map original op index of each `SINGLE_MEAS` to its 0..n_meas-1 order in `circuit.ops`
    orig_op_to_meas_idx: dict[int, int] = {}
    m_run = 0
    for idx, op in enumerate(circuit.ops):
        if op.kind in SINGLE_MEAS:
            orig_op_to_meas_idx[idx] = m_run
            m_run += 1

    ops_by_tick: list[list[tuple[int, CircuitOp]]] = [[] for _ in range(total_ticks)]
    for t_use, idx, op in placed_ops:
        ops_by_tick[t_use].append((idx, op))

    live_now: set[int] = set(circuit.input_qubits)
    all_locs: list[ScheduledLocation] = []
    tick_loc_groups: list[tuple[ScheduledLocation, ...]] = []

    loc_id = 0
    wire_ticks = 0
    peak_q = len(live_now)
    n_data_base = len(circuit.input_qubits)
    peak_anc = max(0, len(live_now) - n_data_base)

    cnot_cnt = 0
    meas_cnt = 0
    prep_cnt = 0
    sq_cnt = 0
    idle_cnt = 0

    for t in range(total_ticks):
        # Births at tick t (`PREP_*`) enter `live_now` at start of tick t
        for _, op in ops_by_tick[t]:
            if op.kind in PREP_STATES:
                live_now.add(op.qubits[0])

        active_snapshot = sorted(live_now)
        wire_ticks += len(active_snapshot)
        if len(active_snapshot) > peak_q:
            peak_q = len(active_snapshot)
        cur_anc = max(0, len(active_snapshot) - n_data_base)
        if cur_anc > peak_anc:
            peak_anc = cur_anc

        used_in_tick: set[int] = set()
        this_tick_locs: list[ScheduledLocation] = []
        deaths_after_tick: list[int] = []

        for idx, op in ops_by_tick[t]:
            for q in op.qubits:
                used_in_tick.add(q)
            lbl = op.label or f"{op.kind}_t{t}_op{idx}"
            sloc = ScheduledLocation(
                loc_id, t, op.kind, op.qubits, lbl, is_idle=False, meas_index=orig_op_to_meas_idx.get(idx, -1)
            )
            loc_id += 1
            all_locs.append(sloc)
            this_tick_locs.append(sloc)

            if op.kind in TWO_QUBIT_GATES:
                cnot_cnt += 1
            elif op.kind in SINGLE_MEAS:
                meas_cnt += 1
                deaths_after_tick.append(op.qubits[0])
            elif op.kind in PREP_STATES:
                prep_cnt += 1
            elif op.kind in SINGLE_QUBIT_GATES:
                sq_cnt += 1

        for q in active_snapshot:
            if q not in used_in_tick:
                sloc = ScheduledLocation(loc_id, t, "IDLE", (q,), f"IDLE_t{t}_q{q}", is_idle=True)
                loc_id += 1
                idle_cnt += 1
                all_locs.append(sloc)
                this_tick_locs.append(sloc)

        tick_loc_groups.append(tuple(this_tick_locs))
        for q_dead in deaths_after_tick:
            live_now.discard(q_dead)

    return TickSchedule(
        n_qubits_total=circuit.n_qubits,
        peak_concurrent_qubits=peak_q,
        peak_concurrent_ancillas=peak_anc,
        ticks=total_ticks,
        cnot_count=cnot_cnt,
        meas_count=meas_cnt,
        prep_count=prep_cnt,
        single_gate_count=sq_cnt,
        idle_location_count=idle_cnt,
        total_locations=len(all_locs),
        wire_tick_volume=wire_ticks,
        locations=tuple(all_locs),
        tick_ops=tuple(tick_loc_groups),
    )


class _StabilizerSubspace:
    """Polynomial-time O(n_out) stabiliser subspace membership and low-weight coset reducer."""

    def __init__(self, gen_rows: tuple[int, ...] | list[int], n_bits: int) -> None:
        self.n_bits = n_bits
        rref_rows, piv_cols, _ = rref_gf2(gen_rows, n_bits)
        self.rref_rows = tuple(rref_rows)
        self.piv_cols = tuple(piv_cols)
        if len(self.rref_rows) <= 10:
            self.elements: list[int] | None = _expand_rowspace(self.rref_rows, n_bits)
        else:
            self.elements = None

    def reduce_vec(self, vec: int) -> int:
        cur = vec
        for row, pcol in zip(self.rref_rows, self.piv_cols):
            if (cur >> pcol) & 1:
                cur ^= row
        return cur

    def contains(self, vec: int) -> bool:
        return self.reduce_vec(vec) == 0

    def min_coset_rep(self, err_mask: int) -> tuple[int, int]:
        if self.contains(err_mask):
            return 0, 0
        for q in range(self.n_bits):
            bit = 1 << q
            if self.contains(err_mask ^ bit):
                return 1, bit
        if self.elements is not None:
            return _min_coset_rep(err_mask, self.elements, self.n_bits)
        # Check weight 2 in O(n_bits^2) even when subspace has dimension > 10!
        for q1 in range(self.n_bits):
            b1 = 1 << q1
            for q2 in range(q1 + 1, self.n_bits):
                cand = b1 | (1 << q2)
                if self.contains(err_mask ^ cand):
                    return 2, cand
        red = self.reduce_vec(err_mask)
        return red.bit_count(), red


def _joint_min_coset_rep(
    out_x: int,
    out_z: int,
    x_space: _StabilizerSubspace,
    z_space: _StabilizerSubspace,
    n_out: int,
) -> tuple[int, int, int]:
    """Minimize joint Pauli support weight `((out_x ^ sx) | (out_z ^ sz)).bit_count()` in polynomial time."""
    if x_space.contains(out_x) and z_space.contains(out_z):
        return 0, 0, 0
    for q in range(n_out):
        bit = 1 << q
        if x_space.contains(out_x ^ bit) and z_space.contains(out_z):
            return 1, bit, 0
        if x_space.contains(out_x) and z_space.contains(out_z ^ bit):
            return 1, 0, bit
        if x_space.contains(out_x ^ bit) and z_space.contains(out_z ^ bit):
            return 1, bit, bit
    if (
        x_space.elements is not None
        and z_space.elements is not None
        and len(x_space.elements) * len(z_space.elements) <= 4096
    ):
        best_w = n_out + 1
        best_rx = out_x
        best_rz = out_z
        for sx in x_space.elements:
            cx = out_x ^ sx
            for sz in z_space.elements:
                cz = out_z ^ sz
                w = (cx | cz).bit_count()
                if w < best_w or (w == best_w and (cx, cz) < (best_rx, best_rz)):
                    best_w = w
                    best_rx = cx
                    best_rz = cz
        return best_w, best_rx, best_rz
    _, rx = x_space.min_coset_rep(out_x)
    _, rz = z_space.min_coset_rep(out_z)
    return (rx | rz).bit_count(), rx, rz


def _stabilizers_for_output(
    circuit: CliffordCircuit,
    codespace: CSSCodeSpec | None,
) -> tuple[_StabilizerSubspace, _StabilizerSubspace, int]:
    """Return `(x_space, z_space, n_out)` on `0..len(output_qubits)-1` for coset weight reduction."""
    n_out = len(circuit.output_qubits)
    if codespace is not None and n_out == codespace.n:
        if len(circuit.input_qubits) == 0:
            return (
                _StabilizerSubspace(codespace.hx_rows, n_out),
                _StabilizerSubspace(codespace.hz_rows + codespace.lz_rows, n_out),
                n_out,
            )
        return (
            _StabilizerSubspace(codespace.hx_rows, n_out),
            _StabilizerSubspace(codespace.hz_rows, n_out),
            n_out,
        )

    if n_out >= 2:
        n_in = len(circuit.input_qubits)
        in_mask = (1 << (2 * n_in)) - 1
        pw = analyze_pauli_webs(circuit_to_zx_fe(circuit, "fault_free"))
        x_rows: list[int] = []
        z_rows: list[int] = []
        for sym, sgn in pw.signed_boundary_stabs:
            if sgn != 0 or (sym & in_mask) != 0:
                continue
            out_sym = sym >> (2 * n_in)
            z_part = 0
            x_part = 0
            for j in range(n_out):
                if (out_sym >> (2 * j)) & 1:
                    z_part |= 1 << j
                if (out_sym >> (2 * j + 1)) & 1:
                    x_part |= 1 << j
            if z_part == 0 and x_part != 0:
                x_rows.append(x_part)
            elif x_part == 0 and z_part != 0:
                z_rows.append(z_part)
        return _StabilizerSubspace(x_rows, n_out), _StabilizerSubspace(z_rows, n_out), n_out

    return _StabilizerSubspace((), n_out), _StabilizerSubspace((), n_out), n_out


def propagate_single_fault(
    sched: TickSchedule, output_qubits: tuple[int, ...], target_loc_id: int,
    fx_qmask: int, fz_qmask: int,
) -> tuple[int, int, int]:
    """Raw output Pauli and measurement displacements for a scheduled physical fault."""
    xm = 0
    zm = 0
    meas_bits = 0
    for tick_locs in sched.tick_ops:
        for sloc in tick_locs:
            isf = sloc.loc_id == target_loc_id
            k = sloc.kind
            if k in SINGLE_MEAS and isf:
                xm ^= fx_qmask
                zm ^= fz_qmask
            if k in PREP_STATES:
                q = sloc.qubits[0]
                xm &= ~(1 << q)
                zm &= ~(1 << q)
            elif k == "H":
                q = sloc.qubits[0]
                xb, zb = (xm >> q) & 1, (zm >> q) & 1
                xm = (xm & ~(1 << q)) | (zb << q)
                zm = (zm & ~(1 << q)) | (xb << q)
            elif k in ("S", "SDG"):
                q = sloc.qubits[0]
                if (xm >> q) & 1:
                    zm ^= 1 << q
            elif k == "Z":
                pass
            elif k == "X":
                pass
            elif k == "Y":
                pass
            elif k == "CX":
                c, t = sloc.qubits
                if (xm >> c) & 1:
                    xm ^= 1 << t
                if (zm >> t) & 1:
                    zm ^= 1 << c
            elif k == "CZ":
                q0, q1 = sloc.qubits
                if (xm >> q0) & 1:
                    zm ^= 1 << q1
                if (xm >> q1) & 1:
                    zm ^= 1 << q0
            elif k == "SWAP":
                q0, q1 = sloc.qubits
                x0, x1 = (xm >> q0) & 1, (xm >> q1) & 1
                z0, z1 = (zm >> q0) & 1, (zm >> q1) & 1
                xm = (xm & ~((1 << q0) | (1 << q1))) | (x1 << q0) | (x0 << q1)
                zm = (zm & ~((1 << q0) | (1 << q1))) | (z1 << q0) | (z0 << q1)
            elif k in SINGLE_MEAS:
                q = sloc.qubits[0]
                m_idx = sloc.meas_index
                if k == "MEAS_Z":
                    bit = (xm >> q) & 1
                elif k == "MEAS_X":
                    bit = (zm >> q) & 1
                else:
                    bit = ((xm >> q) ^ (zm >> q)) & 1
                if bit:
                    meas_bits |= 1 << m_idx
                xm &= ~(1 << q)
                zm &= ~(1 << q)
            if (k not in SINGLE_MEAS) and isf:
                xm ^= fx_qmask
                zm ^= fz_qmask

    out_x = 0
    out_z = 0
    for o_i, q in enumerate(output_qubits):
        if (xm >> q) & 1:
            out_x |= 1 << o_i
        if (zm >> q) & 1:
            out_z |= 1 << o_i
    return out_x, out_z, meas_bits

def build_gadget_contract(
    circuit: CliffordCircuit,
    p_phys: str | Fraction = "1e-4",
    idle_p_ratio: str | Fraction = "1",
    noise_model: str = "circuit_depolarizing",
    codespace: CSSCodeSpec | None = None,
    postselect_detectors: bool | None = None,
    extracted_sector: str = "X",
    certified_w: int = 3,
    t_allow: int = 1,
) -> dict[str, Any]:
    """Audit `circuit` and emit a zero-float outward-rounded `ftqc_zx-contract-v1` dictionary."""
    if noise_model not in ("circuit_depolarizing", "circuit_css"):
        raise RefusalError("CONTRACT-NOISE-MODEL", f"unsupported contract noise model {noise_model!r}")
    p = parse_exact_rate(p_phys)
    idle_ratio = parse_exact_rate(idle_p_ratio) if str(idle_p_ratio).strip() != "1" else Fraction(1)
    p_idle = p * idle_ratio
    if p_idle >= 1:
        raise RefusalError("RATE-OUT-OF-RANGE", f"p_idle = {p_idle} >= 1")

    sched = schedule_circuit_ticks(circuit)
    n_out = len(circuit.output_qubits)

    # Default operating mode: state-preparation circuits (`len(input_qubits) == 0`) postselect flags,
    # whereas inline syndrome-refresh circuits (`len(input_qubits) > 0`) run deterministically!
    if postselect_detectors is None:
        postselect = len(circuit.input_qubits) == 0
    else:
        postselect = postselect_detectors

    # Compute noiseless detector parity checks over measurements via Appendix A `compile_outcome_map`!
    # Wait: if `codespace` is provided on an inline refresh circuit (`len(input_qubits) == codespace.n`),
    # `compile_outcome_map` on `circuit` ALONE (without prepending encoder) finds the INTERNAL ancilla
    # detector checks (such as flags `steane5_flag == 0` or Shor `k3 ^ k4 == 0`), WITHOUT treating the
    # syndrome measurements themselves as postselected zero checks! That is EXACTLY what we want!
    zx_ff = circuit_to_zx_fe(circuit, "fault_free")
    ana_ff = analyze_pauli_webs(zx_ff)
    omap = compile_outcome_map(zx_ff)
    meas_label_to_bit = {lbl: i for i, lbl in enumerate(omap.meas_labels)}
    det_parity_checks = omap.parity_checks  # bitmasks over `0..len(omap.meas_labels)-1`

    # Precompute Appendix A (Proposition A.1) classical Pauli-frame compensation for valid
    # measurement outcomes `mbits in ker(det_parity_checks)` by solving directly over pure
    # OUTPUT-wire Pauli corrections `D^{P(k)} = C_out * D^I` (Kissinger Eq. 9, p. 50)!
    n_in = len(circuit.input_qubits)
    n_meas = len(omap.meas_labels)
    n_b_tot = n_in + n_out
    out_dual_rows = [
        swap_zx_bits(s_vec >> (2 * n_in), n_out) for s_vec, _ in ana_ff.signed_boundary_stabs
    ]
    outcome_frame_pivots: list[tuple[int, int, int]] = []  # (bit_idx_in_mbits, frame_x_mask, frame_z_mask)
    if n_meas > 0 and n_out > 0 and omap.valid_outcome_basis:
        basis_rref, piv_cols_m, transforms_m = rref_gf2(
            list(omap.valid_outcome_basis), n_meas, track_transform=True
        )
        n_stab_rows = len(ana_ff.signed_boundary_stabs)
        for r_i, pcol in enumerate(piv_cols_m):
            comb = transforms_m[r_i]
            syn_piv = 0
            for j_b, b_syn in enumerate(omap.basis_boundary_syndromes):
                if (comb >> j_b) & 1:
                    syn_piv ^= b_syn
            rhs_bits = [(syn_piv >> s_i) & 1 for s_i in range(n_stab_rows)]
            ok_pure, c_out_symp, _ = solve_gf2(out_dual_rows, rhs_bits, 2 * n_out)
            if (not ok_pure) and n_in == n_out:
                # Non-trivial measurement branch (Kissinger Appendix A, Eq. 7 p. 48):
                # Check single-qubit symmetric branch flip `B_sym = P_{in,q} P_{out,q}` + pure output frame `C_out`!
                in_dual_rows = [
                    swap_zx_bits(s_vec & ((1 << (2 * n_in)) - 1), n_in)
                    for s_vec, _ in ana_ff.signed_boundary_stabs
                ]
                for q_br in range(n_in):
                    for sym_pauli in (2, 1, 3):  # try X, Z, Y symmetric branch flip on wire q_br
                        sym_mask = sym_pauli << (2 * q_br)
                        adj_rhs = [
                            rhs_bits[s_i]
                            ^ (((in_dual_rows[s_i] & sym_mask).bit_count() ^ (out_dual_rows[s_i] & sym_mask).bit_count()) & 1)
                            for s_i in range(n_stab_rows)
                        ]
                        ok_br, cand_c_out, _ = solve_gf2(out_dual_rows, adj_rhs, 2 * n_out)
                        if ok_br:
                            ok_pure = True
                            c_out_symp = cand_c_out
                            break
                    if ok_pure:
                        break
            if ok_pure and c_out_symp != 0:
                fx_piv = 0
                fz_piv = 0
                for o_i in range(n_out):
                    if (c_out_symp >> (2 * o_i)) & 1:
                        fz_piv |= 1 << o_i
                    if (c_out_symp >> (2 * o_i + 1)) & 1:
                        fx_piv |= 1 << o_i
                outcome_frame_pivots.append((pcol, fx_piv, fz_piv))

    # Optional Appendix C.2 flag-conditioned decoder when inline on a CSS code block:
    flag_decoder: FlagConditionedDecoder | None = None
    if (
        (not postselect)
        and codespace is not None
        and tuple(circuit.input_qubits) == tuple(range(codespace.n))
        and tuple(circuit.output_qubits) == tuple(range(codespace.n))
        and len(omap.meas_labels) > 0
    ):
        try:
            flag_decoder = synthesize_flag_decoder(circuit, codespace, extracted_sector=extracted_sector)
        except RefusalError:
            flag_decoder = None

    x_stabs, z_stabs, _ = _stabilizers_for_output(circuit, codespace)

    # Audit every location and fault label!
    total_labels = 0
    accepted_labels = 0
    rejected_labels = 0
    bad_unheralded_labels = 0
    bad_css_sector_labels = 0
    bad_post_decoder_labels = 0

    active_x_locs_sets: list[set[int]] = [set() for _ in range(n_out)]
    active_z_locs_sets: list[set[int]] = [set() for _ in range(n_out)]

    n_benign_zero_locs = 0
    n_local_w1_locs = 0
    n_flagged_hook_locs = 0
    n_bad_unheralded_locs = 0

    # Per-location exact label fractions for stochastic polynomial bounds:
    loc_probs: list[Fraction] = []
    loc_reject_fracs: list[Fraction] = []
    loc_x_coord_fracs: list[list[Fraction]] = [[] for _ in range(n_out)]
    loc_z_coord_fracs: list[list[Fraction]] = [[] for _ in range(n_out)]
    loc_bad_fracs: list[Fraction] = []

    for sloc in sched.locations:
        p_loc = p_idle if sloc.is_idle else p
        loc_probs.append(p_loc)
        qs = sloc.qubits
        # Enumerate fault labels `(fx_qmask, fz_qmask, combined_coset_check)`:
        labels: list[tuple[int, int]] = []
        if len(qs) == 1:
            q0 = qs[0]
            if noise_model == "circuit_css":
                labels = [(1 << q0, 0), (0, 1 << q0)]
            else:
                labels = [(1 << q0, 0), (0, 1 << q0), (1 << q0, 1 << q0)]
        else:
            q0, q1 = qs
            if noise_model == "circuit_css":
                for sub in (1, 2, 3):
                    m = ((1 << q0) if (sub & 1) else 0) | ((1 << q1) if (sub & 2) else 0)
                    labels.append((m, 0))
                    labels.append((0, m))
            else:
                for px in range(4):
                    for pz in range(4):
                        if (px, pz) == (0, 0):
                            continue
                        mx = ((1 << q0) if (px & 1) else 0) | ((1 << q1) if (px & 2) else 0)
                        mz = ((1 << q0) if (pz & 1) else 0) | ((1 << q1) if (pz & 2) else 0)
                        labels.append((mx, mz))

        n_labs = len(labels)
        total_labels += n_labs
        rej_here = 0
        bad_here = 0
        max_unflagged_wt_here = 0
        saw_flagged_hook_here = False
        x_hits_here = [0] * n_out
        z_hits_here = [0] * n_out

        for fx_m, fz_m in labels:
            out_x, out_z, mbits = propagate_single_fault(sched, circuit.output_qubits, sloc.loc_id, fx_m, fz_m)
            det_fired = any(((chk & mbits).bit_count() & 1) == 1 for chk in det_parity_checks)

            if (not det_fired) and mbits != 0 and outcome_frame_pivots:
                for bit_idx, fx_piv, fz_piv in outcome_frame_pivots:
                    if (mbits >> bit_idx) & 1:
                        out_x ^= fx_piv
                        out_z ^= fz_piv

            wx, repx = x_stabs.min_coset_rep(out_x)
            wz, repz = z_stabs.min_coset_rep(out_z)
            if noise_model == "circuit_depolarizing":
                w_joint, repx, repz = _joint_min_coset_rep(out_x, out_z, x_stabs, z_stabs, n_out)
            else:
                w_joint = max(wx, wz)

            if det_fired:
                rejected_labels += 1
                rej_here += 1
                if w_joint > t_allow:
                    saw_flagged_hook_here = True
                # If inline deterministic (`not postselect`) with a flag decoder, verify recovery!
                if (not postselect) and flag_decoder is not None:
                    if flag_decoder.propagated_error_type == "X":
                        f_mask = 0
                        for f_i, flbl in enumerate(flag_decoder.flag_meas_labels):
                            if (mbits >> meas_label_to_bit[flbl]) & 1:
                                f_mask |= 1 << f_i
                        if f_mask != 0 and codespace is not None:
                            nsyn = _syndrome_bits(out_x, codespace.hz_rows)
                            rec = flag_decoder.lookup_flagged.get((f_mask, nsyn), 0)
                            res_wt, _ = x_stabs.min_coset_rep(out_x ^ rec)
                            if res_wt > 0:
                                bad_post_decoder_labels += 1
                    else:
                        f_mask = 0
                        for f_i, flbl in enumerate(flag_decoder.flag_meas_labels):
                            if (mbits >> meas_label_to_bit[flbl]) & 1:
                                f_mask |= 1 << f_i
                        if f_mask != 0 and codespace is not None:
                            nsyn = _syndrome_bits(out_z, codespace.hx_rows)
                            rec = flag_decoder.lookup_flagged.get((f_mask, nsyn), 0)
                            res_wt, _ = z_stabs.min_coset_rep(out_z ^ rec)
                            if res_wt > 0:
                                bad_post_decoder_labels += 1
            else:
                accepted_labels += 1
                if w_joint > max_unflagged_wt_here:
                    max_unflagged_wt_here = w_joint
                if w_joint > t_allow:
                    bad_unheralded_labels += 1
                    bad_here += 1
                if max(wx, wz) > t_allow:
                    bad_css_sector_labels += 1
                for c in range(n_out):
                    if (repx >> c) & 1:
                        active_x_locs_sets[c].add(sloc.loc_id)
                        x_hits_here[c] += 1
                    if (repz >> c) & 1:
                        active_z_locs_sets[c].add(sloc.loc_id)
                        z_hits_here[c] += 1

        # Classify location in `v3a` Filter-Calculus partition:
        if bad_here > 0:
            n_bad_unheralded_locs += 1
        elif saw_flagged_hook_here:
            n_flagged_hook_locs += 1
        elif max_unflagged_wt_here > 0:
            n_local_w1_locs += 1
        else:
            n_benign_zero_locs += 1

        loc_reject_fracs.append(Fraction(rej_here, n_labs))
        loc_bad_fracs.append(Fraction(bad_here, n_labs))
        for c in range(n_out):
            loc_x_coord_fracs[c].append(Fraction(x_hits_here[c], n_labs))
            loc_z_coord_fracs[c].append(Fraction(z_hits_here[c], n_labs))

    # Exact all-orders fault-count polynomial `dp[0 .. k_keep]` in exact `Fraction`!
    k_keep = max(2, certified_w)
    dp: list[Fraction] = [Fraction(0)] * k_keep
    dp[0] = Fraction(1)
    for p_l in loc_probs:
        q_l = Fraction(1) - p_l
        for j in range(k_keep - 1, 0, -1):
            dp[j] = dp[j] * q_l + dp[j - 1] * p_l
        dp[0] = dp[0] * q_l

    p_zero_exact = dp[0]
    p_one_exact = dp[1] if k_keep > 1 else Fraction(0)
    p_ge2_exact = Fraction(1) - p_zero_exact - p_one_exact
    p_ge_w_exact = Fraction(1) - sum(dp[:certified_w])

    # Exact 1-fault accepted probability contribution:
    # sum_l p_l * (1 - reject_frac_l) * prod_{m != l} (1 - p_m)
    p_one_accepted = Fraction(0)
    p_one_bad = Fraction(0)
    p_one_x_coord = [Fraction(0) for _ in range(n_out)]
    p_one_z_coord = [Fraction(0) for _ in range(n_out)]

    for idx_l, p_l in enumerate(loc_probs):
        if p_l == 0:
            continue
        # prod_{m != l} (1 - p_m) = dp[0] / (1 - p_l)
        others_clean = p_zero_exact / (Fraction(1) - p_l)
        w_l = p_l * others_clean
        p_one_accepted += w_l * (Fraction(1) - loc_reject_fracs[idx_l])
        p_one_bad += w_l * loc_bad_fracs[idx_l]
        for c in range(n_out):
            p_one_x_coord[c] += w_l * loc_x_coord_fracs[c][idx_l]
            p_one_z_coord[c] += w_l * loc_z_coord_fracs[c][idx_l]

    if postselect:
        # Certified lower bound on acceptance probability: 0-fault + accepted 1-fault trajectories!
        p_acc_lo_exact = p_zero_exact + p_one_accepted
        # Certified upper bound on acceptance probability: 1 - rejected 1-fault trajectories!
        p_one_rejected = p_one_exact - p_one_accepted
        p_acc_hi_exact = Fraction(1) - p_one_rejected
    else:
        p_acc_lo_exact = Fraction(1)
        p_acc_hi_exact = Fraction(1)

    # Certified conditional escape upper bound (unheralded coset weight > t_allow):
    # Any bad event requires EITHER a bad weight-1 fault (`p_one_bad`, which is 0 when certified!)
    # OR >= 2 faults (`p_ge2_exact`)!
    raw_escape_hi_exact = p_one_bad + p_ge2_exact
    cond_escape_hi_exact = min(Fraction(1), raw_escape_hi_exact / p_acc_lo_exact)

    # Certified per-coordinate outgoing error upper bounds (`S_ec` cone-reset coins in `v2c`/`v3a`):
    qx_hi_strs: list[str] = []
    qz_hi_strs: list[str] = []
    for c in range(n_out):
        qx_exact_ub = min(Fraction(1), (p_one_x_coord[c] + p_ge2_exact) / p_acc_lo_exact)
        qz_exact_ub = min(Fraction(1), (p_one_z_coord[c] + p_ge2_exact) / p_acc_lo_exact)
        _, qx_hi = fraction_to_outward_interval_strs(qx_exact_ub)
        _, qz_hi = fraction_to_outward_interval_strs(qz_exact_ub)
        qx_hi_strs.append(qx_hi)
        qz_hi_strs.append(qz_hi)

    exp_wire_ticks_hi_exact = Fraction(sched.wire_tick_volume, 1) / p_acc_lo_exact

    p_zero_lo_s, p_zero_hi_s = fraction_to_outward_interval_strs(p_zero_exact)
    p_acc_lo_s, _ = fraction_to_outward_interval_strs(p_acc_lo_exact)
    _, p_acc_hi_s = fraction_to_outward_interval_strs(p_acc_hi_exact)
    _, p_ge2_hi_s = fraction_to_outward_interval_strs(p_ge2_exact)
    _, p_gew_hi_s = fraction_to_outward_interval_strs(p_ge_w_exact)
    _, esc_hi_s = fraction_to_outward_interval_strs(cond_escape_hi_exact)
    _, exp_wt_hi_s = fraction_to_outward_interval_strs(exp_wire_ticks_hi_exact)

    vx_counts = [len(s) for s in active_x_locs_sets]
    vz_counts = [len(s) for s in active_z_locs_sets]

    verdict = (
        "CONTRACT-CERTIFIED"
        if (bad_unheralded_labels == 0 and bad_post_decoder_labels == 0)
        else "CONTRACT-UNSAFE-SINGLE-FAULT"
    )

    return {
        "schema": "ftqc_zx-contract-v1",
        "verdict": verdict,
        "circuit_name": circuit.name,
        "noise_model": noise_model,
        "codespace": codespace.name if codespace else None,
        "operating_mode": "postselected_factory" if postselect else "deterministic_inline",
        "certified_w": certified_w,
        "t_allow": t_allow,
        "rates": {
            "p_phys_rational": f"{p.numerator}/{p.denominator}",
            "idle_p_ratio_rational": f"{idle_ratio.numerator}/{idle_ratio.denominator}",
            "p_idle_rational": f"{p_idle.numerator}/{p_idle.denominator}",
        },
        "schedule": {
            "n_qubits_total": sched.n_qubits_total,
            "n_data_in": len(circuit.input_qubits),
            "n_data_out": n_out,
            "peak_concurrent_qubits": sched.peak_concurrent_qubits,
            "peak_concurrent_ancillas": sched.peak_concurrent_ancillas,
            "ticks": sched.ticks,
            "cnot_count": sched.cnot_count,
            "meas_count": sched.meas_count,
            "prep_count": sched.prep_count,
            "single_gate_count": sched.single_gate_count,
            "idle_location_count": sched.idle_location_count,
            "total_locations": sched.total_locations,
            "wire_tick_volume": sched.wire_tick_volume,
        },
        "single_fault_audit": {
            "total_labels": total_labels,
            "accepted_labels": accepted_labels,
            "rejected_labels": rejected_labels,
            "bad_unheralded_labels": bad_unheralded_labels,
            "bad_css_sector_labels": bad_css_sector_labels,
            "bad_post_decoder_labels": bad_post_decoder_labels,
            "active_x_locs_per_output": vx_counts,
            "active_z_locs_per_output": vz_counts,
            "max_active_x_locs": max(vx_counts, default=0),
            "max_active_z_locs": max(vz_counts, default=0),
        },
        "filter_calculus_v3a": {
            "benign_zero_locs": n_benign_zero_locs,
            "local_weight1_locs": n_local_w1_locs,
            "flagged_hook_locs": n_flagged_hook_locs,
            "bad_unheralded_locs": n_bad_unheralded_locs,
            "is_1_filtered": bad_unheralded_labels == 0,
        },
        "enclosures": {
            "p_zero_fault": [p_zero_lo_s, p_zero_hi_s],
            "p_accept": [p_acc_lo_s, p_acc_hi_s],
            "p_ge2_faults_hi": p_ge2_hi_s,
            "p_ge_w_faults_hi": p_gew_hi_s,
            "conditional_escape_hi": esc_hi_s,
            "per_output_x_dirt_hi": qx_hi_strs,
            "per_output_z_dirt_hi": qz_hi_strs,
            "expected_wire_ticks_per_shot_hi": exp_wt_hi_s,
        },
        "flag_decoder": flag_decoder.to_dict() if flag_decoder is not None else None,
    }
