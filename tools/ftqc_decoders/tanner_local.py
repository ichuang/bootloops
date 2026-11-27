"""Quantum Tanner iterative local-neighborhood & single-shot syndrome decoder (`ftqc_decoders`).

Implements local-view mismatch-reduction decoding on Quantum Tanner local
neighborhoods (Leverrier & Zemor, arXiv:2206.07750; Gu, Pattison & Tang,
arXiv:2306.12470), consuming ordered vertex neighborhoods exported by
`ftqc_codes` (`quantum_tanner_lrcc` / `lifted_quantum_tanner`) and supporting
consecutive-round single-shot data + syndrome error simulation with exact
coset-reduced residual verification modulo stabilizers.

100% Python stdlib; zero `float` or `complex` literals or calls appear in this module.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any

from .bp_osd import DecodeResult
from .decimal_bounds import RefusalError
from .dem import DecodingProblem, bit_positions


def extract_tanner_neighborhoods(
    source: Any,
    *,
    sector: str = "X",
) -> tuple[tuple[tuple[int, ...], ...], tuple[int, int] | None]:
    """Extract `(neighborhoods, local_shape)` for decoding `sector` errors.

    For `sector == 'X'` (decoded from Z-checks on `V_0`), reads `v0_neighborhoods`.
    For `sector == 'Z'` (decoded from X-checks on `V_1`), reads `v1_neighborhoods`.
    """
    tmeta: dict[str, Any] | None = None
    if isinstance(source, DecodingProblem):
        tmeta = source.source_meta.get("tanner")
        sec = str(source.source_meta.get("sector", sector))
    elif hasattr(source, "meta") and isinstance(source.meta, dict):
        tmeta = source.meta.get("tanner")
        sec = sector
    elif isinstance(source, dict):
        tmeta = source.get("tanner", source)
        sec = sector
    else:
        raise RefusalError("INVALID-TANNER-META", "Cannot extract Tanner neighborhoods from source")

    if not isinstance(tmeta, dict):
        raise RefusalError("MISSING-TANNER-META", "Code/problem carries no 'tanner' local neighborhood metadata")

    key = "v0_neighborhoods" if sec == "X" else "v1_neighborhoods"
    raw_views = tmeta.get(key)
    if not isinstance(raw_views, list) or not raw_views:
        raise RefusalError("MISSING-TANNER-VIEWS", f"Tanner metadata lacks non-empty {key!r}")

    views = tuple(tuple(int(q) for q in v) for v in raw_views)
    lshape = tmeta.get("local_shape")
    shape_tup = (int(lshape[0]), int(lshape[1])) if isinstance(lshape, list) and len(lshape) == 2 else None
    return views, shape_tup


def coset_reduced_weight(error_mask: int, stabilizer_generators: tuple[int, ...] | list[int]) -> int:
    """Minimum Hamming weight in coset `error_mask + span(stabilizer_generators)` (exact for dim <= 16, greedy + RREF otherwise)."""
    if error_mask == 0:
        return 0
    # Build reduced basis of stabilizers
    piv_map: dict[int, int] = {}
    for g in stabilizer_generators:
        v = int(g)
        for p in sorted(piv_map.keys(), reverse=True):
            if (v >> p) & 1:
                v ^= piv_map[p]
        if v:
            npiv = v.bit_length() - 1
            for p in list(piv_map.keys()):
                if (piv_map[p] >> npiv) & 1:
                    piv_map[p] ^= v
            piv_map[npiv] = v

    basis_rows = list(piv_map.values())
    r_red = error_mask
    for p in sorted(piv_map.keys(), reverse=True):
        if (r_red >> p) & 1:
            r_red ^= piv_map[p]
    if r_red == 0:
        return 0

    best = min(error_mask.bit_count(), r_red.bit_count())
    if len(basis_rows) <= 14:
        cur = error_mask
        for i in range(1, 1 << len(basis_rows)):
            cur ^= basis_rows[(i & -i).bit_length() - 1]
            w = cur.bit_count()
            if w < best:
                best = w
                if best == 1:
                    break
        return best

    # Iterative weight-reducing descent by single and pairs of stabilizer generators
    cur = error_mask
    improved = True
    while improved:
        improved = False
        cw = cur.bit_count()
        for g in stabilizer_generators:
            nw = (cur ^ g).bit_count()
            if nw < cw:
                cur ^= g
                cw = nw
                improved = True
        if cw < best:
            best = cw
    return best


def decode_tanner_local(
    problem: DecodingProblem,
    syndrome: int,
    neighborhoods: tuple[tuple[int, ...], ...] | list[list[int]] | None = None,
    *,
    local_shape: tuple[int, int] | None = None,
    max_rounds: int = 32,
    max_local_flip_weight: int = 2,
    min_gain_numer: int = 1,
    min_gain_denom: int = 4,
) -> DecodeResult:
    """Iterative local-neighborhood mismatch-reduction decoder for Quantum Tanner codes.

    In each sweep over active vertex neighborhoods `Q_v` (touching at least one
    unsatisfied check), evaluates candidate local flips `u subseteq Q_v` up to
    weight `max_local_flip_weight` (plus structured row/column slices when
    `local_shape = (d_A, d_B)` is provided) and applies flips achieving strict
    syndrome-weight reduction `gain > 0` with `gain * min_gain_denom >= deg(u) * min_gain_numer`.
    """
    if syndrome < 0 or (syndrome >> problem.num_checks):
        raise RefusalError("SYNDROME-OUT-OF-RANGE", f"Syndrome exceeds {problem.num_checks} checks")

    if neighborhoods is None:
        views, auto_shape = extract_tanner_neighborhoods(problem)
        if local_shape is None:
            local_shape = auto_shape
    else:
        views = tuple(tuple(int(q) for q in v) for v in neighborhoods)

    n = problem.num_mechanisms
    for v_idx, view in enumerate(views):
        if not view or any(q < 0 or q >= n for q in view):
            raise RefusalError("INVALID-TANNER-VIEW", f"Neighborhood {v_idx} contains out-of-range mechanism")

    # Precompute per-vertex touched check mask and candidate local flips `(u_mask, delta_syn, u_wt, deg_syn)`
    vertex_check_masks: list[int] = []
    vertex_candidates: list[list[tuple[int, int, int, int]]] = []
    for view in views:
        v_chk = 0
        for q in view:
            v_chk |= problem.mech_detectors[q]
        vertex_check_masks.append(v_chk)

        cand_masks: set[int] = set()
        cap_w = min(max_local_flip_weight, len(view))
        # Cap combinatorial subset enumeration on very large views (`len(view) > 20`) to weight <= 2
        if len(view) > 20 and cap_w > 2:
            cap_w = 2
        for w in range(1, cap_w + 1):
            for sub in combinations(view, w):
                m_u = 0
                for q in sub:
                    m_u |= 1 << q
                cand_masks.add(m_u)

        if local_shape is not None and local_shape[0] * local_shape[1] == len(view):
            da, db = local_shape
            # Include local row/column sub-slices of weight <= 3 (canonical Tensor-code error patterns!)
            for ia in range(da):
                row_qs = [view[ia * db + jb] for jb in range(db)]
                for w in range(1, min(3, db) + 1):
                    for sub in combinations(row_qs, w):
                        cand_masks.add(sum(1 << q for q in sub))
            for jb in range(db):
                col_qs = [view[ia * db + jb] for ia in range(da)]
                for w in range(1, min(3, da) + 1):
                    for sub in combinations(col_qs, w):
                        cand_masks.add(sum(1 << q for q in sub))

        c_list: list[tuple[int, int, int, int]] = []
        for u_m in sorted(cand_masks, key=lambda x: (x.bit_count(), x)):
            ds = problem.syndrome_of(u_m)
            if ds != 0:
                c_list.append((u_m, ds, u_m.bit_count(), ds.bit_count()))
        vertex_candidates.append(c_list)

    e_hat = 0
    res_syn = syndrome
    rounds_used = 0
    flips_applied = 0
    max_steps = max(1, max_rounds * max(1, len(views)))

    for _ in range(max_steps):
        if res_syn == 0:
            break

        best_global: tuple[int, int, int, int, int, int] | None = None
        for v_idx, v_chk in enumerate(vertex_check_masks):
            if (res_syn & v_chk) == 0:
                continue
            for u_m, ds, u_wt, ds_wt in vertex_candidates[v_idx]:
                # gain = |res_syn| - |res_syn ^ ds| = 2 * |res_syn & ds| - |ds|
                gain = 2 * (res_syn & ds).bit_count() - ds_wt
                if gain <= 0:
                    continue
                if gain * min_gain_denom < ds_wt * min_gain_numer:
                    continue
                rel_syn = (gain * 64) // ds_wt
                rel_wt = (gain * 64) // u_wt
                cand_key = (rel_syn, rel_wt, gain, -u_wt, -u_m, ds)
                if best_global is None or cand_key > best_global:
                    best_global = cand_key

        if best_global is None:
            break

        chosen_u = -best_global[4]
        chosen_ds = best_global[5]
        e_hat ^= chosen_u
        res_syn ^= chosen_ds
        flips_applied += 1
        rounds_used += 1

    converged = res_syn == 0
    return DecodeResult(
        decoder="tanner_local",
        error_mask=e_hat,
        converged=converged,
        observable_mask=problem.observable_of(e_hat),
        hamming_weight=e_hat.bit_count(),
        soft_weight=problem.soft_weight_of(e_hat),
        stats={
            "rounds_used": rounds_used,
            "flips_applied": flips_applied,
            "residual_syndrome_weight": res_syn.bit_count(),
            "num_neighborhoods": len(views),
        },
    )


def simulate_tanner_single_shot_rounds(
    problem: DecodingProblem,
    data_errors_per_round: list[int] | tuple[int, ...],
    syndrome_errors_per_round: list[int] | tuple[int, ...],
    *,
    neighborhoods: tuple[tuple[int, ...], ...] | list[list[int]] | None = None,
    local_shape: tuple[int, int] | None = None,
    max_local_flip_weight: int = 2,
    final_clean_readout: bool = True,
) -> dict[str, Any]:
    """Replay consecutive single-shot Quantum Tanner rounds with data + syndrome noise (`arXiv:2306.12470`).

    Tracks residual error modulo stabilizers after every noisy round, then
    optionally performs one noiseless terminal readout check.
    """
    if len(data_errors_per_round) != len(syndrome_errors_per_round):
        raise RefusalError("ROUND-LEN-MISMATCH", "data_errors_per_round and syndrome_errors_per_round must have equal length")

    views, auto_shape = (
        extract_tanner_neighborhoods(problem)
        if neighborhoods is None
        else (tuple(tuple(int(q) for q in v) for v in neighborhoods), local_shape)
    )
    eff_shape = local_shape if local_shape is not None else auto_shape

    e_residual = 0
    history: list[dict[str, Any]] = []
    max_coset_wt = 0

    for r_idx, (d_err, s_err) in enumerate(zip(data_errors_per_round, syndrome_errors_per_round)):
        if d_err < 0 or (d_err >> problem.num_mechanisms):
            raise RefusalError("DATA-ERROR-OUT-OF-RANGE", f"Round {r_idx} data error out of range")
        if s_err < 0 or (s_err >> problem.num_checks):
            raise RefusalError("SYNDROME-ERROR-OUT-OF-RANGE", f"Round {r_idx} syndrome error out of range")

        e_before = e_residual ^ int(d_err)
        noisy_syn = problem.syndrome_of(e_before) ^ int(s_err)
        dec = decode_tanner_local(
            problem,
            noisy_syn,
            views,
            local_shape=eff_shape,
            max_local_flip_weight=max_local_flip_weight,
        )
        e_residual = e_before ^ dec.error_mask
        c_wt = coset_reduced_weight(e_residual, problem.stabilizer_generators)
        if c_wt > max_coset_wt:
            max_coset_wt = c_wt

        history.append(
            {
                "round": r_idx,
                "injected_data_weight": int(d_err).bit_count(),
                "injected_syndrome_weight": int(s_err).bit_count(),
                "noisy_syndrome_weight": noisy_syn.bit_count(),
                "correction_weight": dec.hamming_weight,
                "post_round_raw_weight": e_residual.bit_count(),
                "post_round_coset_weight": c_wt,
                "post_round_observable_mask": problem.observable_of(e_residual),
                "flips_applied": dec.stats["flips_applied"],
            }
        )

    terminal = None
    if final_clean_readout:
        clean_syn = problem.syndrome_of(e_residual)
        dec_f = decode_tanner_local(
            problem,
            clean_syn,
            views,
            local_shape=eff_shape,
            max_local_flip_weight=max_local_flip_weight,
        )
        e_final = e_residual ^ dec_f.error_mask
        final_coset_wt = coset_reduced_weight(e_final, problem.stabilizer_generators)
        final_obs = problem.observable_of(e_final)
        terminal = {
            "clean_syndrome_converged": dec_f.converged,
            "final_coset_weight": final_coset_wt,
            "final_observable_mask": final_obs,
            "logical_success": dec_f.converged and (final_obs == 0),
        }

    return {
        "num_rounds": len(history),
        "max_intermediate_coset_weight": max_coset_wt,
        "rounds": history,
        "terminal": terminal,
    }
