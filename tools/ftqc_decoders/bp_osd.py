"""Zero-float Belief Propagation (`BP`), Ordered Statistics Decoding (`BP+OSD`),
and Localized Statistics Decoding (`BP+LSD`) for quantum LDPC & spacetime DEMs.

Literature & Algorithmic Lineage:
  1. Iterative quantum BP & symmetric degenerate trapping sets:
     D. Poulin and Y. Chung, "On the iterative decoding of sparse quantum codes",
     Quantum Inf. Comput. 8(10), 987-1000 (2008), arXiv:0801.1241.
  2. Normalized Min-Sum (`NMS`) message passing:
     J. Chen, A. Dholakia, E. Eleftheriou, M. P. C. Fossorier, and X.-Y. Hu,
     IEEE Trans. Commun. 53(8), 1288-1299 (2005).
  3. Ordered Statistics Decoding (`OSD-0`, `OSD-E`, `OSD-CS` combination sweep):
     M. P. C. Fossorier and S. Lin, IEEE Trans. Inf. Theory 41(5), 1379-1396 (1995);
     P. Panteleev and G. Kalachev, "Degenerate Quantum LDPC Codes With Good Finite
     Length Performance", Quantum 5, 585 (2021), arXiv:1904.02703;
     J. Roffe, D. R. White, S. Burton, and E. Campbell, Phys. Rev. Research 2,
     043423 (2020), arXiv:2005.07016.
  4. Localized Statistics Decoding (`BP+LSD`, cluster-local on-the-fly elimination):
     T. Hillmann, L. Berent, A. O. Quintavalle, J. Eisert, R. Wille, and J. Roffe,
     "Localized statistics decoding for quantum low-density parity-check codes",
     Nature Commun. 16, 8214 (2025), arXiv:2406.18655.

Purity invariant:
  ZERO `float` or `complex` literals/calls appear anywhere in this module.
  * `method='nms'` runs purely in scaled `int` arithmetic with rational damping `alpha`.
  * `method='sum_product'` runs purely in `fractions.Fraction` probability-bias
    coordinates `delta = 1 - 2p in (-1, 1)` (exact rational algebra without `tanh`/`ln`).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Any

from .decimal_bounds import RefusalError
from .dem import DecodingProblem, bit_positions


@dataclass(frozen=True)
class DecodeResult:
    """Standardized deterministic result across all `ftqc_decoders` algorithms."""

    decoder: str
    error_mask: int
    converged: bool
    observable_mask: int
    hamming_weight: int
    soft_weight: int
    stats: dict[str, Any]


def run_bp_nms(
    problem: DecodingProblem,
    syndrome: int,
    *,
    max_iters: int = 30,
    schedule: str = "serial",
    alpha: Fraction = Fraction(3, 4),
) -> tuple[int, bool, list[int], int]:
    """Scaled-integer Normalized Min-Sum (`NMS`) Belief Propagation.

    Returns `(hard_error_mask, converged, posterior_llr_list, iterations_used)`.
    """
    if schedule not in ("parallel", "serial"):
        raise RefusalError("INVALID-BP-SCHEDULE", f"schedule must be 'parallel' or 'serial', got {schedule!r}")
    if not (Fraction(0) < alpha <= Fraction(1)):
        raise RefusalError("INVALID-NMS-ALPHA", f"NMS damping alpha must lie in (0, 1], got {alpha}")
    if max_iters < 1:
        raise RefusalError("INVALID-BP-ITERS", f"max_iters must be >= 1, got {max_iters}")

    m = problem.num_checks
    n = problem.num_mechanisms
    if syndrome == 0:
        return 0, True, list(problem.weights), 0

    check_vars = [bit_positions(r) for r in problem.check_rows]
    var_checks = [bit_positions(c) for c in problem.mech_detectors]
    syn_bits = [(syndrome >> c) & 1 for c in range(m)]
    l0 = list(problem.weights)

    a_num = alpha.numerator
    a_den = alpha.denominator
    cap_msg = max(l0, default=1024) * 32

    # Edge messages: `msg_v2c[c][idx_in_c]` and `msg_c2v[c][idx_in_c]`
    pos_in_check: list[dict[int, int]] = [
        {v: idx for idx, v in enumerate(check_vars[c])} for c in range(m)
    ]
    msg_v2c: list[list[int]] = [[l0[v] for v in check_vars[c]] for c in range(m)]
    msg_c2v: list[list[int]] = [[0 for _ in check_vars[c]] for c in range(m)]
    posterior = list(l0)

    def update_check(c: int) -> None:
        vs = check_vars[c]
        deg = len(vs)
        if deg == 0:
            return
        v2c_row = msg_v2c[c]
        c2v_row = msg_c2v[c]
        total_neg_parity = syn_bits[c]
        min1 = cap_msg + 1
        min2 = cap_msg + 1
        min1_idx = -1
        for k_i in range(deg):
            val = v2c_row[k_i]
            if val < 0:
                total_neg_parity ^= 1
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
            other_parity = total_neg_parity ^ (1 if val < 0 else 0)
            best_mag = min2 if k_i == min1_idx else min1
            scaled = (a_num * best_mag) // a_den
            if scaled > cap_msg:
                scaled = cap_msg
            c2v_row[k_i] = -scaled if other_parity else scaled

    for it in range(1, max_iters + 1):
        # Roffe et al. (Phys. Rev. Research 2, 043423, 2020, Eq. 9) dyadic warm-up
        # `alpha_t = min(alpha, 1 - 2^-t)` prevents Iteration-1 parallel overshoot!
        it_cap = min(it, 10)
        dyn_alpha = min(alpha, Fraction((1 << it_cap) - 1, 1 << it_cap))
        a_num = dyn_alpha.numerator
        a_den = dyn_alpha.denominator
        if schedule == "parallel":
            for c in range(m):
                update_check(c)
            e_hat = 0
            for v in range(n):
                tot = l0[v]
                for c in var_checks[v]:
                    tot += msg_c2v[c][pos_in_check[c][v]]
                posterior[v] = tot
                if tot < 0:
                    e_hat |= 1 << v
            if problem.syndrome_of(e_hat) == syndrome:
                return e_hat, True, posterior, it
            for v in range(n):
                tot = posterior[v]
                for c in var_checks[v]:
                    idx_cv = pos_in_check[c][v]
                    raw_out = tot - msg_c2v[c][idx_cv]
                    if raw_out > cap_msg:
                        raw_out = cap_msg
                    elif raw_out < -cap_msg:
                        raw_out = -cap_msg
                    msg_v2c[c][idx_cv] = raw_out
        else:
            # Serial (variable-layered) schedule (Roffe et al., Phys. Rev. Research 2, 043423, 2020)
            for c in range(m):
                update_check(c)
            for v in range(n):
                vcs = var_checks[v]
                tot = l0[v]
                for c in vcs:
                    tot += msg_c2v[c][pos_in_check[c][v]]
                posterior[v] = tot
                for c in vcs:
                    idx_cv = pos_in_check[c][v]
                    raw_out = tot - msg_c2v[c][idx_cv]
                    if raw_out > cap_msg:
                        raw_out = cap_msg
                    elif raw_out < -cap_msg:
                        raw_out = -cap_msg
                    msg_v2c[c][idx_cv] = raw_out
                    update_check(c)

            e_hat = 0
            for v in range(n):
                tot = l0[v]
                for c in var_checks[v]:
                    tot += msg_c2v[c][pos_in_check[c][v]]
                posterior[v] = tot
                if tot < 0:
                    e_hat |= 1 << v
            if problem.syndrome_of(e_hat) == syndrome:
                return e_hat, True, posterior, it

    e_final = 0
    for v in range(n):
        if posterior[v] < 0:
            e_final |= 1 << v
    return e_final, (problem.syndrome_of(e_final) == syndrome), posterior, max_iters


def run_bp_sum_product_rational(
    problem: DecodingProblem,
    syndrome: int,
    *,
    max_iters: int = 20,
    max_denominator: int = 10**9,
) -> tuple[int, bool, list[int], int]:
    """Exact-rational Sum-Product BP in probability-bias coordinates `delta = 1 - 2p in Q`.

    Uses the algebraic identity `tanh(L/2) = 1 - 2p`, avoiding all transcendental
    functions (`tanh`, `atanh`, `ln`, `exp`) and all floats! Posterior biases are
    converted to scaled integer reliability scores for downstream OSD/LSD ranking.
    """
    m = problem.num_checks
    n = problem.num_mechanisms
    if syndrome == 0:
        return 0, True, list(problem.weights), 0

    check_vars = [bit_positions(r) for r in problem.check_rows]
    var_checks = [bit_positions(c) for c in problem.mech_detectors]
    pos_in_check: list[dict[int, int]] = [
        {v: idx for idx, v in enumerate(check_vars[c])} for c in range(m)
    ]
    syn_signs = [-1 if ((syndrome >> c) & 1) else 1 for c in range(m)]

    # Clamp helper keeping rational bias strictly inside (-1, 1) with bounded denominator
    eps = Fraction(1, max_denominator)
    hi_clamp = Fraction(1, 1) - eps
    lo_clamp = -hi_clamp

    def clamp_bias(d: Fraction) -> Fraction:
        d_lim = d.limit_denominator(max_denominator)
        if d_lim >= hi_clamp:
            return hi_clamp
        if d_lim <= lo_clamp:
            return lo_clamp
        return d_lim

    delta0 = [clamp_bias(Fraction(1, 1) - 2 * p) for p in problem.priors]
    ratio0 = [(Fraction(1, 1) + d) / (Fraction(1, 1) - d) for d in delta0]

    msg_v2c: list[list[Fraction]] = [[delta0[v] for v in check_vars[c]] for c in range(m)]
    msg_c2v: list[list[Fraction]] = [[Fraction(0, 1) for _ in check_vars[c]] for c in range(m)]
    post_bias = list(delta0)

    for it in range(1, max_iters + 1):
        # Check-to-variable bias update: prefix/suffix rational products (O(deg) without division by 0!)
        for c in range(m):
            vs = check_vars[c]
            deg = len(vs)
            if deg == 0:
                continue
            row_in = msg_v2c[c]
            pref = [Fraction(1, 1)] * (deg + 1)
            suff = [Fraction(1, 1)] * (deg + 1)
            for k_i in range(deg):
                pref[k_i + 1] = (pref[k_i] * row_in[k_i]).limit_denominator(max_denominator)
            for k_i in range(deg - 1, -1, -1):
                suff[k_i] = (suff[k_i + 1] * row_in[k_i]).limit_denominator(max_denominator)
            s_sign = syn_signs[c]
            for k_i in range(deg):
                prod_other = pref[k_i] * suff[k_i + 1]
                msg_c2v[c][k_i] = clamp_bias(prod_other * s_sign)

        # Variable posterior & variable-to-check likelihood-ratio update
        e_hat = 0
        for v in range(n):
            vcs = var_checks[v]
            deg_v = len(vcs)
            factors = []
            for c in vcs:
                d_cv = msg_c2v[c][pos_in_check[c][v]]
                factors.append(((Fraction(1, 1) + d_cv) / (Fraction(1, 1) - d_cv)).limit_denominator(max_denominator))
            pref_r = [Fraction(1, 1)] * (deg_v + 1)
            suff_r = [Fraction(1, 1)] * (deg_v + 1)
            for k_i in range(deg_v):
                pref_r[k_i + 1] = (pref_r[k_i] * factors[k_i]).limit_denominator(max_denominator)
            for k_i in range(deg_v - 1, -1, -1):
                suff_r[k_i] = (suff_r[k_i + 1] * factors[k_i]).limit_denominator(max_denominator)

            r_tot = (ratio0[v] * pref_r[deg_v]).limit_denominator(max_denominator)
            d_post = clamp_bias((r_tot - 1) / (r_tot + 1))
            post_bias[v] = d_post
            if d_post < 0:
                e_hat |= 1 << v

            for k_i, c in enumerate(vcs):
                r_ext = (ratio0[v] * pref_r[k_i] * suff_r[k_i + 1]).limit_denominator(max_denominator)
                msg_v2c[c][pos_in_check[c][v]] = clamp_bias((r_ext - 1) / (r_ext + 1))

        if problem.syndrome_of(e_hat) == syndrome:
            scores = [int(d * 1_000_000) for d in post_bias]
            return e_hat, True, scores, it

    e_final = 0
    for v, d_post in enumerate(post_bias):
        if d_post < 0:
            e_final |= 1 << v
    scores = [int(d * 1_000_000) for d in post_bias]
    return e_final, (problem.syndrome_of(e_final) == syndrome), scores, max_iters


def _osd_postprocess(
    problem: DecodingProblem,
    syndrome: int,
    posterior_scores: list[int],
    *,
    osd_order: int = 0,
    osd_search_depth: int = 24,
) -> tuple[int, bool, dict[str, Any]]:
    """Global Ordered Statistics Decoding (`OSD-0` & `OSD-CS(osd_order)`).

    Ranks mechanisms from least reliable (`min posterior_scores[j]`) to most reliable,
    pivots `H` into column-RREF over `GF(2)`, solves `OSD-0`, and sweeps combinations
    up to `osd_order` on the `osd_search_depth` least-reliable non-pivot columns.
    """
    n = problem.num_mechanisms
    # Order columns from MOST LIKELY IN ERROR (smallest / most negative score) to MOST RELIABLE
    ordered_cols = sorted(range(n), key=lambda j: (posterior_scores[j], j))

    # Incremental column-space RREF over check-space bitmasks (`0..2^m - 1`):
    # Each basis entry at `row_pivot_bit` stores `(col_vec_in_check_space, mechanism_xor_mask_inducing_it)`.
    basis_vec: dict[int, tuple[int, int]] = {}
    pivot_cols: list[int] = []
    free_cols: list[int] = []
    free_kernel_masks: list[int] = []  # Each entry lies in ker(H): `(1 << free_j) ^ pivot_mask`

    for col_j in ordered_cols:
        v = problem.mech_detectors[col_j]
        comb_mask = 1 << col_j
        for piv_bit in sorted(basis_vec.keys(), reverse=True):
            if (v >> piv_bit) & 1:
                bv, bm = basis_vec[piv_bit]
                v ^= bv
                comb_mask ^= bm
        if v == 0:
            free_cols.append(col_j)
            free_kernel_masks.append(comb_mask)
        else:
            new_piv = v.bit_length() - 1
            # Reduce existing basis vectors at `new_piv` for clean RREF
            for pb in list(basis_vec.keys()):
                bv, bm = basis_vec[pb]
                if (bv >> new_piv) & 1:
                    basis_vec[pb] = (bv ^ v, bm ^ comb_mask)
            basis_vec[new_piv] = (v, comb_mask)
            pivot_cols.append(col_j)

    # Reduce `syndrome` against `basis_vec` to obtain the unique `OSD-0` solution on `pivot_cols`
    rem_syn = syndrome
    e_osd0 = 0
    for piv_bit in sorted(basis_vec.keys(), reverse=True):
        if (rem_syn >> piv_bit) & 1:
            bv, bm = basis_vec[piv_bit]
            rem_syn ^= bv
            e_osd0 ^= bm

    if rem_syn != 0:
        # Syndrome lies outside colspace(H) (infeasible syndrome!)
        return 0, False, {
            "post_stage": "osd_infeasible",
            "eliminated_rank": len(pivot_cols),
            "candidates_evaluated": 0,
        }

    best_e = e_osd0
    best_key = (problem.soft_weight_of(best_e), best_e.bit_count(), best_e)
    candidates = 1

    if osd_order >= 1 and free_kernel_masks:
        # Combination Sweep (`OSD-CS`, Panteleev & Kalachev 2021):
        # 1) All single flips across the first `min(len(free), max(osd_search_depth * 2, 60))` free columns
        w1_limit = min(len(free_kernel_masks), max(osd_search_depth * 2, 60))
        for idx_f in range(w1_limit):
            cand = e_osd0 ^ free_kernel_masks[idx_f]
            candidates += 1
            ck = (problem.soft_weight_of(cand), cand.bit_count(), cand)
            if ck < best_key:
                best_key = ck
                best_e = cand

        # 2) Higher-order combinations `2 .. osd_order` on the top `osd_search_depth` least-reliable free columns
        depth = min(len(free_kernel_masks), osd_search_depth)
        top_kernels = free_kernel_masks[:depth]
        for w in range(2, osd_order + 1):
            if w > depth:
                break
            for comb_indices in combinations(range(depth), w):
                cand = e_osd0
                for idx_f in comb_indices:
                    cand ^= top_kernels[idx_f]
                candidates += 1
                ck = (problem.soft_weight_of(cand), cand.bit_count(), cand)
                if ck < best_key:
                    best_key = ck
                    best_e = cand

    return best_e, True, {
        "post_stage": "osd0" if osd_order == 0 else f"osd_cs_{osd_order}",
        "eliminated_rank": len(pivot_cols),
        "max_cluster_rank": len(pivot_cols),
        "candidates_evaluated": candidates,
    }


def _lsd_postprocess(
    problem: DecodingProblem,
    syndrome: int,
    posterior_scores: list[int],
    *,
    lsd_order: int = 0,
    lsd_search_depth: int = 16,
) -> tuple[int, bool, dict[str, Any]]:
    """Localized Statistics Decoding (`BP+LSD`, Hillmann et al., Nat. Commun. 16, 8214, 2025, arXiv:2406.18655).

    Grows syndrome clusters starting from violated checks `supp(syndrome)` by greedily
    absorbing the least-reliable boundary fault mechanisms, maintaining per-cluster
    on-the-fly `GF(2)` echelon factorizations, and solving independently per cluster!
    """
    m = problem.num_checks
    check_vars = [bit_positions(r) for r in problem.check_rows]

    # Disjoint-set union over active check clusters
    parent_c = [-1] * m  # -1 means check is not yet in any active cluster
    cluster_checks: dict[int, set[int]] = {}
    cluster_mechs: dict[int, list[int]] = {}
    cluster_basis: dict[int, dict[int, tuple[int, int]]] = {}
    cluster_free_kernels: dict[int, list[int]] = {}
    cluster_target_syn: dict[int, int] = {}
    used_mech = [False] * problem.num_mechanisms

    def find_c(c: int) -> int:
        root = c
        while parent_c[root] != root:
            root = parent_c[root]
        cur = c
        while cur != root:
            nxt = parent_c[cur]
            parent_c[cur] = root
            cur = nxt
        return root

    def add_mech_to_cluster(root: int, mech_j: int) -> None:
        cluster_mechs[root].append(mech_j)
        v = problem.mech_detectors[mech_j]
        comb_mask = 1 << mech_j
        bmap = cluster_basis[root]
        for piv_bit in sorted(bmap.keys(), reverse=True):
            if (v >> piv_bit) & 1:
                bv, bm = bmap[piv_bit]
                v ^= bv
                comb_mask ^= bm
        if v == 0:
            cluster_free_kernels[root].append(comb_mask)
        else:
            new_piv = v.bit_length() - 1
            for pb in list(bmap.keys()):
                bv, bm = bmap[pb]
                if (bv >> new_piv) & 1:
                    bmap[pb] = (bv ^ v, bm ^ comb_mask)
            bmap[new_piv] = (v, comb_mask)

    def merge_clusters(r1: int, r2: int) -> int:
        if r1 == r2:
            return r1
        if len(cluster_checks[r1]) < len(cluster_checks[r2]):
            r1, r2 = r2, r1
        parent_c[r2] = r1
        cluster_checks[r1].update(cluster_checks[r2])
        cluster_target_syn[r1] ^= cluster_target_syn[r2]
        combined_mechs = sorted(cluster_mechs[r1] + cluster_mechs[r2], key=lambda j: (posterior_scores[j], j))
        cluster_mechs[r1] = []
        cluster_basis[r1] = {}
        cluster_free_kernels[r1] = []
        for mj in combined_mechs:
            add_mech_to_cluster(r1, mj)
        del cluster_checks[r2]
        del cluster_mechs[r2]
        del cluster_basis[r2]
        del cluster_free_kernels[r2]
        del cluster_target_syn[r2]
        return r1

    def ensure_check_cluster(c: int) -> int:
        if parent_c[c] != -1:
            return find_c(c)
        parent_c[c] = c
        cluster_checks[c] = {c}
        cluster_mechs[c] = []
        cluster_basis[c] = {}
        cluster_free_kernels[c] = []
        cluster_target_syn[c] = (1 << c) if ((syndrome >> c) & 1) else 0
        return c

    def is_cluster_valid(root: int) -> bool:
        rem = cluster_target_syn[root]
        bmap = cluster_basis[root]
        for piv_bit in sorted(bmap.keys(), reverse=True):
            if (rem >> piv_bit) & 1:
                rem ^= bmap[piv_bit][0]
        return rem == 0

    for d_c in bit_positions(syndrome):
        ensure_check_cluster(d_c)

    # Grow invalid clusters by adding the single most-unreliable boundary mechanism
    while True:
        invalid_roots = [r for r in sorted(cluster_checks.keys()) if not is_cluster_valid(r)]
        if not invalid_roots:
            break

        progressed = False
        for r in list(invalid_roots):
            if r not in cluster_checks or is_cluster_valid(r):
                continue
            # Find best unused mechanism incident to any check in `cluster_checks[r]`
            best_j: int | None = None
            best_score: tuple[int, int] | None = None
            for c_in in cluster_checks[r]:
                for j_cand in check_vars[c_in]:
                    if not used_mech[j_cand]:
                        sc = (posterior_scores[j_cand], j_cand)
                        if best_score is None or sc < best_score:
                            best_score = sc
                            best_j = j_cand
            if best_j is None:
                continue
            used_mech[best_j] = True
            progressed = True
            # Ensure all checks touched by `best_j` belong to the same merged cluster!
            cur_root = find_c(r)
            for c_touch in bit_positions(problem.mech_detectors[best_j]):
                r_other = ensure_check_cluster(c_touch)
                cur_root = merge_clusters(cur_root, r_other)
            add_mech_to_cluster(cur_root, best_j)

        if not progressed:
            # Infeasible syndrome component!
            return 0, False, {
                "post_stage": "lsd_infeasible",
                "num_clusters": len(cluster_checks),
                "eliminated_rank": sum(len(b) for b in cluster_basis.values()),
                "max_cluster_rank": max((len(b) for b in cluster_basis.values()), default=0),
                "candidates_evaluated": 0,
            }

    # Higher-order LSD-CS (`lsd_order >= 1`, Hillmann et al. arXiv:2406.18655):
    # Absorb up to `lsd_search_depth` least-reliable unused mechanisms incident to each
    # active cluster. Merging valid/zero-syndrome check clusters preserves validity
    # unconditionally (`s_1 + s_2 in span(M_1 U M_2 U {j})`) while admitting cancelled checks!
    if lsd_order >= 1:
        for init_root in list(sorted(cluster_checks.keys())):
            if init_root not in cluster_checks:
                continue
            for _extra in range(lsd_search_depth):
                r_cur = find_c(init_root)
                best_j_inc: int | None = None
                best_sc_inc: tuple[int, int] | None = None
                for c_in in cluster_checks[r_cur]:
                    for j_cand in check_vars[c_in]:
                        if not used_mech[j_cand]:
                            sc = (posterior_scores[j_cand], j_cand)
                            if best_sc_inc is None or sc < best_sc_inc:
                                best_sc_inc = sc
                                best_j_inc = j_cand
                if best_j_inc is None:
                    break
                used_mech[best_j_inc] = True
                for c_touch in bit_positions(problem.mech_detectors[best_j_inc]):
                    r_other = ensure_check_cluster(c_touch)
                    r_cur = merge_clusters(r_cur, r_other)
                add_mech_to_cluster(r_cur, best_j_inc)

    # Solve independently inside each valid cluster!
    total_e = 0
    total_candidates = 0
    total_rank = 0
    max_c_rank = 0

    for root in sorted(cluster_checks.keys()):
        bmap = cluster_basis[root]
        c_rank = len(bmap)
        total_rank += c_rank
        if c_rank > max_c_rank:
            max_c_rank = c_rank
        rem = cluster_target_syn[root]
        e_c0 = 0
        for piv_bit in sorted(bmap.keys(), reverse=True):
            if (rem >> piv_bit) & 1:
                bv, bm = bmap[piv_bit]
                rem ^= bv
                e_c0 ^= bm
        assert rem == 0

        best_ec = e_c0
        best_ck = (problem.soft_weight_of(best_ec), best_ec.bit_count(), best_ec)
        total_candidates += 1

        fk = cluster_free_kernels[root]
        if lsd_order >= 1 and fk:
            for k_mask in fk:
                cand = e_c0 ^ k_mask
                total_candidates += 1
                ck = (problem.soft_weight_of(cand), cand.bit_count(), cand)
                if ck < best_ck:
                    best_ck = ck
                    best_ec = cand
            depth = min(len(fk), lsd_search_depth)
            for w in range(2, lsd_order + 1):
                if w > depth:
                    break
                for comb_idx in combinations(range(depth), w):
                    cand = e_c0
                    for idx_f in comb_idx:
                        cand ^= fk[idx_f]
                    total_candidates += 1
                    ck = (problem.soft_weight_of(cand), cand.bit_count(), cand)
                    if ck < best_ck:
                        best_ck = ck
                        best_ec = cand

        total_e ^= best_ec

    return total_e, True, {
        "post_stage": "lsd0" if lsd_order == 0 else f"lsd_cs_{lsd_order}",
        "num_clusters": len(cluster_checks),
        "eliminated_rank": total_rank,
        "max_cluster_rank": max_c_rank,
        "candidates_evaluated": total_candidates,
    }


def decode_bp_osd(
    problem: DecodingProblem,
    syndrome: int,
    *,
    bp_method: str = "nms",
    schedule: str = "serial",
    max_iters: int = 30,
    alpha: Fraction = Fraction(3, 4),
    post_method: str = "osd",
    osd_order: int = 0,
    osd_search_depth: int = 24,
) -> DecodeResult:
    """Decode `syndrome` using zero-float `BP` (`nms` or `sum_product`) optionally
    followed by `OSD` (`post_method='osd'`) or localized `LSD` (`post_method='lsd'`),
    or pure `BP` (`post_method='none'`).
    """
    if syndrome < 0 or (syndrome >> problem.num_checks):
        raise RefusalError("SYNDROME-OUT-OF-RANGE", f"Syndrome exceeds {problem.num_checks} checks")
    if post_method not in ("none", "osd", "lsd"):
        raise RefusalError("INVALID-POST-METHOD", f"post_method must be 'none', 'osd', or 'lsd', got {post_method!r}")

    if bp_method == "nms":
        e_bp, conv_bp, scores, iters = run_bp_nms(
            problem, syndrome, max_iters=max_iters, schedule=schedule, alpha=alpha
        )
    elif bp_method == "sum_product":
        e_bp, conv_bp, scores, iters = run_bp_sum_product_rational(
            problem, syndrome, max_iters=max_iters
        )
    else:
        raise RefusalError("INVALID-BP-METHOD", f"bp_method must be 'nms' or 'sum_product', got {bp_method!r}")

    label = (
        f"bp_only:{bp_method}"
        if post_method == "none"
        else (f"bp_osd{osd_order}:{bp_method}" if post_method == "osd" else f"bp_lsd{osd_order}:{bp_method}")
    )

    if post_method == "none" or (conv_bp and e_bp.bit_count() <= 1):
        return DecodeResult(
            decoder=label,
            error_mask=e_bp,
            converged=conv_bp,
            observable_mask=problem.observable_of(e_bp),
            hamming_weight=e_bp.bit_count(),
            soft_weight=problem.soft_weight_of(e_bp),
            stats={
                "bp_method": bp_method,
                "schedule": schedule if bp_method == "nms" else "parallel",
                "bp_iterations": iters,
                "bp_converged_directly": conv_bp,
                "post_stage": "none",
                "eliminated_rank": 0,
                "max_cluster_rank": 0,
                "candidates_evaluated": 1,
            },
        )

    if post_method == "osd":
        e_post, conv_post, p_stats = _osd_postprocess(
            problem, syndrome, scores, osd_order=osd_order, osd_search_depth=osd_search_depth
        )
    else:
        e_post, conv_post, p_stats = _lsd_postprocess(
            problem, syndrome, scores, lsd_order=osd_order, lsd_search_depth=osd_search_depth
        )

    # If BP also found a valid solution, keep whichever valid recovery has strictly smaller soft weight!
    if conv_bp and (
        (not conv_post)
        or (problem.soft_weight_of(e_bp), e_bp.bit_count(), e_bp)
        < (problem.soft_weight_of(e_post), e_post.bit_count(), e_post)
    ):
        e_post = e_bp
        conv_post = True

    return DecodeResult(
        decoder=label,
        error_mask=e_post,
        converged=conv_post and (problem.syndrome_of(e_post) == syndrome),
        observable_mask=problem.observable_of(e_post) if conv_post else 0,
        hamming_weight=e_post.bit_count() if conv_post else 0,
        soft_weight=problem.soft_weight_of(e_post) if conv_post else 0,
        stats={
            "bp_method": bp_method,
            "schedule": schedule if bp_method == "nms" else "parallel",
            "bp_iterations": iters,
            "bp_converged_directly": False,
            **p_stats,
        },
    )
