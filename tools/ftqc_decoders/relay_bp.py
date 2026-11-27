"""Clean-room zero-float Disordered-Memory Relay Belief Propagation (`Relay-BP`).

Implements the message-passing equations of:
  * T. Muller et al. (IBM Quantum), "Improved belief propagation is sufficient for
    real-time decoding of quantum memory", arXiv:2506.01779 (2025).
  * Applied to Bivariate Bicycle (`BB72` [[72,12,6]], `BB144` [[144,12,12]]) gross
    codes of S. Bravyi et al., Nature 627, 778-782 (2024), arXiv:2308.07915.

Core mechanism (zero Gaussian elimination, zero floats):
  Standard BP stagnates or oscillates on symmetric quantum trapping sets `e_1 + e_2 = g`
  (where `g` is a stabilizer generator) because every variable in `supp(g)` sees identical
  local graph neighborhoods and prior weights. `Relay-BP` breaks these symmetries via:
    1. Leg 0 (`pre_iters` rounds): standard Normalized Min-Sum (`gamma = 0`).
    2. Relay legs `r = 1 .. num_relays` (`iters_per_relay` rounds each):
       - Draws deterministic per-variable rational memory strengths
         `gamma_v in [gamma_min, gamma_max]` (including negative `gamma_v < 0`, which
         repels stagnant marginals away from zero-crossings!).
       - Sets reference bias `L_ref` to the clamped posterior from the previous leg
         (or mixed with original channel priors `L0`), updating variable bias each round:
           `L_bias[v] = (1 - gamma_v) * L_ref[v] + gamma_v * L_posterior_prev[v]`.
    3. Collects up to `stop_n_solutions` valid syndrome-matching solutions across legs
       and returns the minimum-soft-weight recovery.

Clean-room MIT implementation in 100% Python stdlib (`int` & `Fraction`, zero floats).
"""

from __future__ import annotations

from fractions import Fraction
import hashlib

from .bp_osd import DecodeResult, run_bp_nms
from .decimal_bounds import RefusalError
from .dem import DecodingProblem, bit_positions


def _deterministic_memory_numerators(
    seed: int,
    leg: int,
    num_vars: int,
    num_min: int,
    num_max: int,
) -> list[int]:
    """Deterministic platform-independent integer PRNG over `[num_min, num_max]` via SHA-256."""
    span = num_max - num_min + 1
    out: list[int] = []
    counter = 0
    while len(out) < num_vars:
        digest = hashlib.sha256(f"relay_bp:{seed}:{leg}:{counter}".encode("ascii")).digest()
        for k in range(0, 32, 4):
            if len(out) >= num_vars:
                break
            word = int.from_bytes(digest[k : k + 4], "little")
            out.append(num_min + (word % span))
        counter += 1
    return out


def decode_relay_bp(
    problem: DecodingProblem,
    syndrome: int,
    *,
    pre_iters: int = 20,
    num_relays: int = 16,
    iters_per_relay: int = 20,
    gamma_min: Fraction = Fraction(-1, 4),
    gamma_max: Fraction = Fraction(3, 4),
    alpha: Fraction = Fraction(3, 4),
    stop_n_solutions: int = 1,
    seed: int = 20260601,
) -> DecodeResult:
    """Decode `syndrome` via clean-room zero-float `Relay-BP` (no `OSD`/`LSD` elimination)."""
    if syndrome < 0 or (syndrome >> problem.num_checks):
        raise RefusalError("SYNDROME-OUT-OF-RANGE", f"Syndrome exceeds {problem.num_checks} checks")
    if gamma_min > gamma_max:
        raise RefusalError("INVALID-GAMMA-WINDOW", f"gamma_min={gamma_min} > gamma_max={gamma_max}")
    if stop_n_solutions < 1:
        raise RefusalError("INVALID-STOP-SOLUTIONS", "stop_n_solutions must be >= 1")

    if syndrome == 0:
        return DecodeResult(
            decoder="relay_bp",
            error_mask=0,
            converged=True,
            observable_mask=0,
            hamming_weight=0,
            soft_weight=0,
            stats={"leg_converged": 0, "total_iterations": 0, "solutions_found": 1},
        )

    # Leg 0: initial standard NMS pre-run
    e0, conv0, post_prev, it0 = run_bp_nms(
        problem, syndrome, max_iters=pre_iters, schedule="parallel", alpha=alpha
    )
    solutions: list[tuple[int, int, int, int]] = []  # (soft_weight, hamming_weight, error_mask, leg_id)
    if conv0:
        solutions.append((problem.soft_weight_of(e0), e0.bit_count(), e0, 0))
        if len(solutions) >= stop_n_solutions:
            return DecodeResult(
                decoder="relay_bp",
                error_mask=e0,
                converged=True,
                observable_mask=problem.observable_of(e0),
                hamming_weight=e0.bit_count(),
                soft_weight=problem.soft_weight_of(e0),
                stats={"leg_converged": 0, "total_iterations": it0, "solutions_found": 1},
            )

    m = problem.num_checks
    n = problem.num_mechanisms
    check_vars = [bit_positions(r) for r in problem.check_rows]
    var_checks = [bit_positions(c) for c in problem.mech_detectors]
    pos_in_check: list[dict[int, int]] = [
        {v: idx for idx, v in enumerate(check_vars[c])} for c in range(m)
    ]
    syn_bits = [(syndrome >> c) & 1 for c in range(m)]
    l0 = list(problem.weights)
    cap_msg = max(l0, default=1024) * 32

    # Common integer denominator for gamma disorder grid (`G_DEN = 256` mixed with exact bounds)
    g_den = 256 * gamma_min.denominator * gamma_max.denominator
    g_num_min = int(gamma_min * g_den)
    g_num_max = int(gamma_max * g_den)
    a_num = alpha.numerator
    a_den = alpha.denominator

    total_iters = it0
    last_leg_posterior = list(post_prev)

    for leg in range(1, num_relays + 1):
        gammas = _deterministic_memory_numerators(seed, leg, n, g_num_min, g_num_max)
        # Relay reference bias blends channel prior `l0` with clamped previous-leg posterior
        l_ref = [
            max(-cap_msg // 2, min(cap_msg // 2, (l0[v] + last_leg_posterior[v]) // 2))
            for v in range(n)
        ]
        l_post = list(l_ref)
        msg_v2c: list[list[int]] = [[l_ref[v] for v in check_vars[c]] for c in range(m)]
        msg_c2v: list[list[int]] = [[0 for _ in check_vars[c]] for c in range(m)]

        leg_hit = False
        for _it in range(1, iters_per_relay + 1):
            total_iters += 1
            for c in range(m):
                vs = check_vars[c]
                deg = len(vs)
                if deg == 0:
                    continue
                v2c_row = msg_v2c[c]
                c2v_row = msg_c2v[c]
                total_neg = syn_bits[c]
                min1 = cap_msg + 1
                min2 = cap_msg + 1
                min1_idx = -1
                for k_i in range(deg):
                    val = v2c_row[k_i]
                    if val < 0:
                        total_neg ^= 1
                        mag = -val
                    else:
                        mag = val
                    if mag < min1:
                        min2 = min1
                        min1 = mag
                        min1_idx = k_i
                    elif mag < min2:
                        min2 = mag
                for k_i in range(deg):
                    val = v2c_row[k_i]
                    other_neg = total_neg ^ (1 if val < 0 else 0)
                    best_mag = min2 if k_i == min1_idx else min1
                    sc = (a_num * best_mag) // a_den
                    if sc > cap_msg:
                        sc = cap_msg
                    c2v_row[k_i] = -sc if other_neg else sc

            e_hat = 0
            for v in range(n):
                g_num = gammas[v]
                # Disordered memory prior: `(1 - gamma_v) * l_ref[v] + gamma_v * l_post_prev[v]`
                mem_prior = ((g_den - g_num) * l_ref[v] + g_num * l_post[v]) // g_den
                tot = mem_prior
                for c in var_checks[v]:
                    tot += msg_c2v[c][pos_in_check[c][v]]
                if tot > cap_msg:
                    tot = cap_msg
                elif tot < -cap_msg:
                    tot = -cap_msg
                l_post[v] = tot
                if tot < 0:
                    e_hat |= 1 << v

                for c in var_checks[v]:
                    idx_cv = pos_in_check[c][v]
                    raw_out = tot - msg_c2v[c][idx_cv]
                    if raw_out > cap_msg:
                        raw_out = cap_msg
                    elif raw_out < -cap_msg:
                        raw_out = -cap_msg
                    msg_v2c[c][idx_cv] = raw_out

            if problem.syndrome_of(e_hat) == syndrome:
                solutions.append((problem.soft_weight_of(e_hat), e_hat.bit_count(), e_hat, leg))
                leg_hit = True
                break

        last_leg_posterior = l_post
        if leg_hit and len(solutions) >= stop_n_solutions:
            break

    if solutions:
        solutions.sort()
        best_sw, best_hw, best_e, best_leg = solutions[0]
        return DecodeResult(
            decoder="relay_bp",
            error_mask=best_e,
            converged=True,
            observable_mask=problem.observable_of(best_e),
            hamming_weight=best_hw,
            soft_weight=best_sw,
            stats={
                "leg_converged": best_leg,
                "total_iterations": total_iters,
                "solutions_found": len(solutions),
            },
        )

    e_fail = 0
    for v in range(n):
        if last_leg_posterior[v] < 0:
            e_fail |= 1 << v
    return DecodeResult(
        decoder="relay_bp",
        error_mask=e_fail,
        converged=False,
        observable_mask=problem.observable_of(e_fail),
        hamming_weight=e_fail.bit_count(),
        soft_weight=problem.soft_weight_of(e_fail),
        stats={
            "leg_converged": None,
            "total_iterations": total_iters,
            "solutions_found": 0,
        },
    )
