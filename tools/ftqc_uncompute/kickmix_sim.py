"""64-lane bit-packed Fiat-Shamir KickMix trajectory simulator and zero-float bounds.

Implements randomized trajectory verification and certified error enclosures from:
  - Khattar et al., "Verifiable Quantum Circuit Compilation" (arXiv:2610.23832v2, Section II.B-C)

100% Python standard library, MIT-licensed, zero float literals or float() calls.
All real-valued statistical and trace-distance bounds use outward-rounded `decimal.Decimal`.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
import hashlib
import math
from typing import Any, Callable

from .decimal_bounds import sqrt_upper

from .ir_kickmix import (
    ANNOTATION_OPS,
    KickMixCircuit,
    Qubit,
    z_pow_hierarchy_level,
)

# Strictly greater than pi (pi = 3.14159265358979323846264338327950288419716939937510...)
PI_UPPER = Decimal("3.1415926535897932384626433832795028841972")
LANE_WIDTH = 64
LANE_MASK = (1 << LANE_WIDTH) - 1


class FiatShamirPRNG:
    """Platform-independent deterministic CSPRNG driven by counter-mode SHA-256."""

    def __init__(self, seed_material: bytes) -> None:
        self._seed = hashlib.sha256(seed_material).digest()
        self.seed_hex = self._seed.hex()
        self._counter = 0
        self._buf: list[int] = []

    def next_u64(self) -> int:
        if not self._buf:
            ctr_bytes = self._counter.to_bytes(8, "little")
            self._counter += 1
            digest = hashlib.sha256(self._seed + ctr_bytes).digest()
            self._buf = [
                int.from_bytes(digest[0:8], "little"),
                int.from_bytes(digest[8:16], "little"),
                int.from_bytes(digest[16:24], "little"),
                int.from_bytes(digest[24:32], "little"),
            ]
        return self._buf.pop()

    def randbelow(self, modulus: int) -> int:
        if modulus <= 1:
            return 0
        bit_len = (modulus - 1).bit_length()
        words = (bit_len + 63) // 64
        mask = (1 << bit_len) - 1
        while True:
            val = 0
            for w in range(words):
                val |= self.next_u64() << (64 * w)
            val &= mask
            if val < modulus:
                return val


def fiat_shamir_seed(circuit: KickMixCircuit, salt: str = "ftqc_uncompute-fs-v1") -> bytes:
    """Canonical comment-independent Fiat-Shamir seed material."""
    return circuit.canonical_bytes() + b"\x00" + salt.encode("utf-8")


def _binom_cdf_upper(k: int, n: int, p: Decimal, ctx_ceil: Context) -> Decimal:
    """Outward-rounded (`ROUND_CEILING`) upper bound on `P[Binomial(n, p) <= k]`."""
    if p <= Decimal(0):
        return Decimal(1)
    if p >= Decimal(1):
        return Decimal(1) if k >= n else Decimal(0)
    # To upper-bound sum_{j=0}^k C(n, j) p^j (1 - p)^(n - j), round (1 - p) UP and p UP.
    q_up = ctx_ceil.subtract(Decimal(1), p)
    total = Decimal(0)
    # For k == 0 (dominant zero-failure case), exact 1-term power:
    if k == 0:
        return ctx_ceil.power(q_up, Decimal(n))
    for j in range(0, k + 1):
        comb = Decimal(math.comb(n, j))
        pj = Decimal(1) if j == 0 else ctx_ceil.power(p, Decimal(j))
        qnj = Decimal(1) if j == n else ctx_ceil.power(q_up, Decimal(n - j))
        term = ctx_ceil.multiply(comb, ctx_ceil.multiply(pj, qnj))
        total = ctx_ceil.add(total, term)
    return total


def clopper_pearson_upper_bound(
    failures: int,
    trials: int,
    delta_neg_log2: int = 128,
    decimal_digits: int = 24,
) -> Decimal:
    """Certified outward-rounded (`ROUND_CEILING`) Clopper-Pearson upper bound `eps^+`
    at failure significance `delta = 2 ** (-delta_neg_log2)`.

    Guarantees mathematically that `P[Binomial(trials, eps^+) <= failures] <= delta`,
    using zero floating-point operations.
    """
    if trials < 1:
        raise ValueError("trials must be >= 1")
    if not (0 <= failures <= trials):
        raise ValueError(f"failures ({failures}) must lie in [0, {trials}]")
    if delta_neg_log2 < 1:
        raise ValueError("delta_neg_log2 must be >= 1")
    if failures == trials:
        return Decimal(1)

    work_prec = max(64, decimal_digits + 32 + (delta_neg_log2 // 3))
    ctx_floor = Context(prec=work_prec, rounding=ROUND_FLOOR)
    ctx_ceil = Context(prec=work_prec, rounding=ROUND_CEILING)
    ctx_out = Context(prec=decimal_digits, rounding=ROUND_CEILING)

    # Lower-bound delta so `cdf_upper(p_hi) <= delta_lo <= 2^(-s)` is strictly conservative.
    two_pow_s = Decimal(1 << delta_neg_log2)
    delta_lo = ctx_floor.divide(Decimal(1), two_pow_s)

    if failures == 0:
        # Solve (1 - p)^n <= delta_lo  <=>  p >= 1 - delta_lo^(1/n)
        # Compute initial estimate via ln/exp, then self-certify via bisection if needed.
        ln_delta = ctx_floor.ln(delta_lo)  # negative, rounded toward -inf -> smaller exp -> larger p!
        exponent = ctx_floor.divide(ln_delta, Decimal(trials))
        q_est = ctx_floor.exp(exponent)
        p_cand = ctx_ceil.subtract(Decimal(1), q_est)
        if p_cand > Decimal(1):
            p_cand = Decimal(1)
        if _binom_cdf_upper(0, trials, p_cand, ctx_ceil) <= delta_lo:
            return ctx_out.create_decimal(p_cand)

    # Certified bisection on [0, 1] using outward-rounded `_binom_cdf_upper`.
    lo = ctx_floor.divide(Decimal(failures), Decimal(trials)) if failures > 0 else Decimal(0)
    hi = Decimal(1)
    steps = max(120, decimal_digits * 4)
    half = Decimal("0.5")
    for _ in range(steps):
        mid = ctx_ceil.multiply(ctx_ceil.add(lo, hi), half)
        if mid >= hi or mid <= lo:
            break
        if _binom_cdf_upper(failures, trials, mid, ctx_ceil) <= delta_lo:
            hi = mid
        else:
            lo = mid
    return ctx_out.create_decimal(hi)


def coherent_trace_distance_bound(
    epsilon_upper: Decimal,
    phase_tol_turns: Fraction | int | str | Decimal = 0,
    decimal_digits: int = 24,
) -> Decimal:
    """Outward-rounded KickMix coherent trace-distance bound (Khattar et al. Eq. 3):

    `D_tr <= min(1, 2 * sqrt(epsilon_upper) + 2 * pi * tau)`,
    valid when branch maps `V_r` are injective on the input support.
    """
    if isinstance(phase_tol_turns, float):
        raise TypeError("float phase tolerance forbidden; use Fraction, int, str, or Decimal")
    ctx_ceil = Context(prec=max(50, decimal_digits + 16), rounding=ROUND_CEILING)
    ctx_out = Context(prec=decimal_digits, rounding=ROUND_CEILING)

    if isinstance(phase_tol_turns, Decimal):
        tau = phase_tol_turns
    else:
        frac_tau = Fraction(phase_tol_turns)
        tau = ctx_ceil.divide(Decimal(frac_tau.numerator), Decimal(frac_tau.denominator))
    if tau < Decimal(0):
        raise ValueError("phase_tol_turns must be >= 0")

    root_eps = sqrt_upper(epsilon_upper, ctx_ceil)
    term1 = ctx_ceil.multiply(Decimal(2), root_eps)
    term2 = ctx_ceil.multiply(ctx_ceil.multiply(Decimal(2), PI_UPPER), tau)
    bound = ctx_ceil.add(term1, term2)
    if bound >= Decimal(1):
        return Decimal(1)
    return ctx_out.create_decimal(bound)


def z_pow_synthesis_cost(
    circuit: KickMixCircuit,
    total_rot_error_neg_log2: int = 30,
) -> dict[str, Any]:
    """Clifford+T synthesis cost breakdown for `Z_POW` and Toffoli-like operations.

    Uses standard Gidney/Jones `4 T` per Toffoli-like (`CCX`/`CCZ`/`CSWAP`) with clean ancilla,
    `1 T` per `Z_POW(±1/4)` (`T`/`T†`), and gridsynth upper bound `ceil(3 * log2(M / eps_rot))`
    per finer rotation (`level >= 4`).
    """
    stats = circuit.operation_statistics()
    m_rot = stats.num_rotations
    if m_rot == 0:
        t_per_rot = 0
    else:
        # ceil(3 * (total_rot_error_neg_log2 + ceil(log2(m_rot))))
        extra_bits = (m_rot - 1).bit_length()
        t_per_rot = 3 * (total_rot_error_neg_log2 + extra_bits)
    rot_t_total = m_rot * t_per_rot
    toff_t_total = 4 * stats.num_toffoli_like
    direct_t = stats.num_t_like
    return {
        "toffoli_like": stats.num_toffoli_like,
        "toffoli_t_cost": toff_t_total,
        "direct_t_gates": direct_t,
        "fine_rotations": m_rot,
        "t_per_fine_rotation": t_per_rot,
        "rotation_t_cost": rot_t_total,
        "total_t_upper": toff_t_total + direct_t + rot_t_total,
        "rotation_budget_neg_log2": total_rot_error_neg_log2,
    }


def dqi_repetition_factor(epsilon_upper: Decimal, decimal_digits: int = 18) -> Decimal | None:
    """Outward-rounded post-selection repetition overhead `1 / (1 - epsilon_upper)^2`
    when `epsilon_upper < 1` (Khattar et al. Section IV).
    """
    if epsilon_upper >= Decimal(1):
        return None
    ctx_floor = Context(prec=max(48, decimal_digits + 16), rounding=ROUND_FLOOR)
    ctx_ceil = Context(prec=decimal_digits, rounding=ROUND_CEILING)
    one_minus = ctx_floor.subtract(Decimal(1), epsilon_upper)
    if one_minus <= Decimal(0):
        return None
    denom = ctx_floor.multiply(one_minus, one_minus)
    if denom <= Decimal(0):
        return None
    return ctx_ceil.divide(Decimal(1), denom)


def _min_circular_half_turn_diff(a: Fraction, b: Fraction) -> Fraction:
    """Shortest distance modulo 2 half-turns in `[0, 1]` half-turns."""
    diff = (a - b) % Fraction(2, 1)
    if diff > Fraction(1, 1):
        diff = Fraction(2, 1) - diff
    return diff


@dataclass(frozen=True)
class TrajectoryResult:
    shots: int
    failures: int
    dirty_resets: int
    epsilon_upper: Decimal
    delta_neg_log2: int
    phase_tol_turns: Fraction
    coherent_trace_distance_upper: Decimal | None
    coherent_bound_valid: bool
    z_pow_budget: dict[str, Any]
    dqi_repetition_upper: Decimal | None
    fiat_shamir_seed_hex: str
    first_failure_witness: dict[str, Any] | None

    @property
    def passed(self) -> bool:
        return self.failures == 0 and self.dirty_resets == 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "shots": self.shots,
            "failures": self.failures,
            "dirty_resets": self.dirty_resets,
            "epsilon_upper": str(self.epsilon_upper),
            "delta_neg_log2": self.delta_neg_log2,
            "phase_tol_turns": f"{self.phase_tol_turns.numerator}/{self.phase_tol_turns.denominator}",
            "coherent_trace_distance_upper": (
                str(self.coherent_trace_distance_upper)
                if self.coherent_trace_distance_upper is not None
                else None
            ),
            "coherent_bound_valid": self.coherent_bound_valid,
            "z_pow_budget": self.z_pow_budget,
            "dqi_repetition_upper": (
                str(self.dqi_repetition_upper) if self.dqi_repetition_upper is not None else None
            ),
            "fiat_shamir_seed_hex": self.fiat_shamir_seed_hex,
            "first_failure_witness": self.first_failure_witness,
        }


OracleFn = Callable[[dict[str, int]], tuple[dict[str, int], Fraction | int] | dict[str, int]]


def simulate_trajectories(
    circuit: KickMixCircuit,
    *,
    shots: int = 1024,
    oracle: OracleFn | None = None,
    input_sampler: Callable[[FiatShamirPRNG, dict[str, int]], dict[str, int]] | None = None,
    salt: str = "ftqc_uncompute-fs-v1",
    delta_neg_log2: int = 128,
    phase_tol_turns: Fraction | int | str = 0,
    branch_injective: bool | None = None,
    allow_global_phase: bool = False,
    check_scratch_clean: bool = True,
    check_classical_clean: bool = False,
) -> TrajectoryResult:
    """Execute 64-lane bit-packed Fiat-Shamir trajectories on `circuit`.

    If `oracle` is omitted, checks the identity / compute-uncompute contract
    (all public registers returned to their sampled initial values with zero phase).
    """
    if shots < 1:
        raise ValueError("shots must be >= 1")
    if isinstance(phase_tol_turns, float):
        raise TypeError("float phase tolerance forbidden")
    tol_turns = Fraction(phase_tol_turns)
    tol_half_turns = Fraction(2, 1) * tol_turns

    prng = FiatShamirPRNG(fiat_shamir_seed(circuit, salt=salt))

    q_reg_specs = [
        (r.name if r.name is not None else f"r{r.reg_id}", tuple(w.id for w in r.wires if isinstance(w, Qubit)))
        for r in circuit.registers
        if r.kind == "q"
    ]
    if not circuit.registers and circuit.num_qubits > 0:
        q_reg_specs = [("q_all", tuple(range(circuit.num_qubits)))]

    reg_widths = {name: len(qids) for name, qids in q_reg_specs}
    scratch_qids = circuit.scratch_qubits() if circuit.registers else ()

    failures = 0
    dirty_resets = 0
    first_witness: dict[str, Any] | None = None
    global_phase_ref: Fraction | None = None

    nq = circuit.num_qubits
    nb = circuit.num_bits

    shot_offset = 0
    while shot_offset < shots:
        batch = min(LANE_WIDTH, shots - shot_offset)
        batch_mask = (1 << batch) - 1

        # Sample `batch` inputs deterministically
        sampled_inputs: list[dict[str, int]] = []
        q_state = [0] * nq
        c_state = [0] * nb
        for lane in range(batch):
            if input_sampler is not None:
                inp = input_sampler(prng, reg_widths)
            else:
                inp = {name: prng.randbelow(1 << width) for name, width in reg_widths.items()}
            sampled_inputs.append(inp)
            for name, qids in q_reg_specs:
                val = inp[name]
                for bit_idx, qid in enumerate(qids):
                    if (val >> bit_idx) & 1:
                        q_state[qid] |= 1 << lane

        sign_mask = 0
        fine_phases: list[Fraction] = [Fraction(0, 1) for _ in range(batch)]
        cond_stack: list[int] = []
        dirty_lane_mask = 0
        dirty_op_idx: int | None = None
        dirty_qubit: int | None = None

        for op_idx, op in enumerate(circuit.ops):
            kind = op.kind
            if kind in ANNOTATION_OPS:
                continue
            if kind == "PUSH_CONDITION":
                assert op.condition is not None
                top = cond_stack[-1] if cond_stack else batch_mask
                cond_stack.append(top & c_state[op.condition])
                continue
            if kind == "POP_CONDITION":
                if cond_stack:
                    cond_stack.pop()
                continue

            active = cond_stack[-1] if cond_stack else batch_mask
            if op.condition is not None:
                active &= c_state[op.condition]
            if active == 0:
                continue

            if kind == "X":
                q_state[op.qubits[0]] ^= active
            elif kind == "CX":
                c, t = op.qubits
                q_state[t] ^= active & q_state[c]
            elif kind == "CCX":
                c1, c2, t = op.qubits
                q_state[t] ^= active & q_state[c1] & q_state[c2]
            elif kind == "SWAP":
                q1, q2 = op.qubits
                diff = active & (q_state[q1] ^ q_state[q2])
                q_state[q1] ^= diff
                q_state[q2] ^= diff
            elif kind == "CSWAP":
                c, q1, q2 = op.qubits
                diff = active & q_state[c] & (q_state[q1] ^ q_state[q2])
                q_state[q1] ^= diff
                q_state[q2] ^= diff
            elif kind == "NEG":
                sign_mask ^= active
            elif kind == "Z":
                sign_mask ^= active & q_state[op.qubits[0]]
            elif kind == "CZ":
                c1, c2 = op.qubits
                sign_mask ^= active & q_state[c1] & q_state[c2]
            elif kind == "CCZ":
                c1, c2, c3 = op.qubits
                sign_mask ^= active & q_state[c1] & q_state[c2] & q_state[c3]
            elif kind == "Z_POW":
                assert op.angle is not None
                fired = active & q_state[op.qubits[0]]
                if op.angle == 1 or op.angle == -1:
                    sign_mask ^= fired
                elif op.angle != 0 and fired:
                    ang = op.angle
                    m = fired
                    while m:
                        lsb = m & -m
                        lane = lsb.bit_length() - 1
                        fine_phases[lane] = (fine_phases[lane] + ang) % Fraction(2, 1)
                        m ^= lsb
            elif kind == "R":
                q = op.qubits[0]
                dirty = active & q_state[q]
                if dirty:
                    dirty_resets += dirty.bit_count()
                    dirty_lane_mask |= dirty
                    if dirty_op_idx is None:
                        dirty_op_idx = op_idx
                        dirty_qubit = q
                    # KickMix dirty-reset trajectory semantics: random ±1 phase kick!
                    rand_flip = prng.next_u64() & dirty
                    sign_mask ^= rand_flip
                q_state[q] &= ~active
            elif kind == "HMR":
                q = op.qubits[0]
                b = op.bits[0]
                meas = prng.next_u64() & active
                sign_mask ^= q_state[q] & meas
                c_state[b] = (c_state[b] & ~active) | meas
                q_state[q] &= ~active
            elif kind == "BIT_STORE0":
                c_state[op.bits[0]] &= ~active
            elif kind == "BIT_STORE1":
                c_state[op.bits[0]] = (c_state[op.bits[0]] & ~active) | active
            elif kind == "BIT_INVERT":
                c_state[op.bits[0]] ^= active

        # Evaluate each lane in the batch
        for lane in range(batch):
            inp = sampled_inputs[lane]
            if oracle is None:
                exp_regs = inp
                exp_phase = Fraction(0, 1)
            else:
                raw_exp = oracle(inp)
                if isinstance(raw_exp, tuple):
                    exp_regs, raw_ph = raw_exp
                    exp_phase = Fraction(raw_ph) % Fraction(2, 1)
                else:
                    exp_regs = raw_exp
                    exp_phase = Fraction(0, 1)

            got_regs: dict[str, int] = {}
            for name, qids in q_reg_specs:
                val = 0
                for bit_idx, qid in enumerate(qids):
                    if (q_state[qid] >> lane) & 1:
                        val |= 1 << bit_idx
                got_regs[name] = val

            got_phase = (fine_phases[lane] + Fraction((sign_mask >> lane) & 1, 1)) % Fraction(2, 1)
            if allow_global_phase:
                if global_phase_ref is None:
                    global_phase_ref = (got_phase - exp_phase) % Fraction(2, 1)
                exp_eff = (exp_phase + global_phase_ref) % Fraction(2, 1)
            else:
                exp_eff = exp_phase

            phase_err = _min_circular_half_turn_diff(got_phase, exp_eff)
            scratch_dirty_q: int | None = None
            if check_scratch_clean:
                for sq in scratch_qids:
                    if (q_state[sq] >> lane) & 1:
                        scratch_dirty_q = sq
                        break
            classical_dirty_b: int | None = None
            if check_classical_clean:
                for cb_idx in range(nb):
                    if (c_state[cb_idx] >> lane) & 1:
                        classical_dirty_b = cb_idx
                        break

            lane_dirty_reset = bool((dirty_lane_mask >> lane) & 1)
            reg_mismatch = got_regs != exp_regs
            phase_mismatch = phase_err > tol_half_turns

            if lane_dirty_reset or reg_mismatch or phase_mismatch or scratch_dirty_q is not None or classical_dirty_b is not None:
                failures += 1
                if first_witness is None:
                    first_witness = {
                        "shot_index": shot_offset + lane,
                        "inputs": inp,
                        "expected_outputs": exp_regs,
                        "got_outputs": got_regs,
                        "expected_half_turns": f"{exp_eff.numerator}/{exp_eff.denominator}",
                        "got_half_turns": f"{got_phase.numerator}/{got_phase.denominator}",
                        "dirty_reset": lane_dirty_reset,
                        "dirty_reset_op_index": dirty_op_idx if lane_dirty_reset else None,
                        "dirty_reset_qubit": dirty_qubit if lane_dirty_reset else None,
                        "uncleaned_scratch_qubit": scratch_dirty_q,
                        "uncleaned_classical_bit": classical_dirty_b,
                    }

        shot_offset += batch

    eps_up = clopper_pearson_upper_bound(failures, shots, delta_neg_log2=delta_neg_log2)
    if branch_injective is None:
        from .hmr import verify_hmr_circuit
        branch_injective = verify_hmr_circuit(circuit, check_scratch_clean=check_scratch_clean).branch_injective
    coh_valid = bool(branch_injective and dirty_resets == 0)
    coh_up = coherent_trace_distance_bound(eps_up, phase_tol_turns=tol_turns) if coh_valid else None
    z_budget = z_pow_synthesis_cost(circuit)
    dqi_rep = dqi_repetition_factor(eps_up)

    return TrajectoryResult(
        shots=shots,
        failures=failures,
        dirty_resets=dirty_resets,
        epsilon_upper=eps_up,
        delta_neg_log2=delta_neg_log2,
        phase_tol_turns=tol_turns,
        coherent_trace_distance_upper=coh_up,
        coherent_bound_valid=coh_valid,
        z_pow_budget=z_budget,
        dqi_repetition_upper=dqi_rep,
        fiat_shamir_seed_hex=prng.seed_hex,
        first_failure_witness=first_witness,
    )
