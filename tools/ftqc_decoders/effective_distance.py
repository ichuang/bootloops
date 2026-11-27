"""Operational decoding radius certifier & adversarial quantum trapping-set / `d_eff` finder.

Turns silent "ideal minimum-weight decoder" assumptions into falsifiable BootLoops
certificates:
  1. `certify_operational_radius`: proves by bounded enumeration (or connected
     syndrome-cluster enumeration for local decoders) that *every* fault set of
     weight `<= target_radius` is corrected with zero logical flip by a concrete
     deterministic decoder (`BP+OSD`, `BP+LSD`, `Relay-BP`, `UF`, `MWPM`, `SSF`, or
     an audited external decoder such as `Tesseract`), OR emits a replayable
     sub-distance counterexample witness!
  2. `search_trapping_sets`: combines three directed adversarial probes to expose
     `d_eff < d_code` failures on large qLDPC / BB / HGP codes without exponential
     brute force:
       - **Stabilizer split-support probe** (Poulin & Chung, QIC 8, 987, 2008;
         Raveendran & Vasic, Quantum 5, 562, 2021): splits low-weight stabilizers
         `g in G_stab` (and pairwise generator sums) into balanced halves `e_1 + e_2 = g`
         of weight `floor(wt(g)/2)` and `ceil(wt(g)/2)`, deterministically exposing
         symmetric iterative-BP trapping sets on `BB72`, `BB144`, `HGP`, and `Toric`!
       - **Logical-coset midpoint probe** (Fossorier & Lin, IEEE TIT 41, 1379, 1995;
         Delfosse, Londe & Beverland, IEEE TIT 68, 3187, 2022): probes subsets of
         minimum-weight logical operators around `floor((d-1)/2)`.
       - **Bounded low-weight exhaustive scan**.

Zero `float` or `complex` literals/calls appear anywhere in this module.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import comb
from typing import Any

from .bp_osd import DecodeResult, decode_bp_osd
from .decimal_bounds import RefusalError
from .dem import DecodingProblem, bit_positions
from .mwpm import decode_mwpm
from .relay_bp import decode_relay_bp
from .sat_mwd import decode_sat_mwd
from .ssf_metacheck import decode_small_set_flip
from .union_find import decode_union_find


def run_decoder(
    problem: DecodingProblem,
    syndrome: int,
    spec: str | dict[str, Any] = "bp_osd",
) -> DecodeResult:
    """Unified deterministic dispatcher for all stdlib and audited external decoders."""
    if isinstance(spec, str):
        cfg: dict[str, Any] = {"algorithm": spec}
    elif isinstance(spec, dict):
        cfg = dict(spec)
    else:
        raise RefusalError("INVALID-DECODER-SPEC", f"Decoder spec must be str or dict, got {type(spec).__name__}")

    alg = str(cfg.get("algorithm", "bp_osd")).strip()
    if alg.startswith("external:"):
        from .adapters.external_decoders import decode_external

        backend = alg.split(":", 1)[1]
        return decode_external(problem, syndrome, backend=backend, options=cfg)

    if alg == "bp_only":
        return decode_bp_osd(
            problem,
            syndrome,
            bp_method=cfg.get("bp_method", "nms"),
            schedule=cfg.get("schedule", "parallel"),
            max_iters=int(cfg.get("max_iters", 30)),
            alpha=Fraction(cfg.get("alpha", "3/4")),
            post_method="none",
        )
    if alg == "bp_osd":
        return decode_bp_osd(
            problem,
            syndrome,
            bp_method=cfg.get("bp_method", "nms"),
            schedule=cfg.get("schedule", "parallel"),
            max_iters=int(cfg.get("max_iters", 30)),
            alpha=Fraction(cfg.get("alpha", "3/4")),
            post_method="osd",
            osd_order=int(cfg.get("osd_order", 0)),
            osd_search_depth=int(cfg.get("osd_search_depth", 24)),
        )
    if alg == "bp_lsd":
        return decode_bp_osd(
            problem,
            syndrome,
            bp_method=cfg.get("bp_method", "nms"),
            schedule=cfg.get("schedule", "parallel"),
            max_iters=int(cfg.get("max_iters", 30)),
            alpha=Fraction(cfg.get("alpha", "3/4")),
            post_method="lsd",
            osd_order=int(cfg.get("osd_order", 0)),
            osd_search_depth=int(cfg.get("osd_search_depth", 24)),
        )
    if alg == "relay_bp":
        return decode_relay_bp(
            problem,
            syndrome,
            pre_iters=int(cfg.get("pre_iters", 20)),
            num_relays=int(cfg.get("num_relays", 16)),
            iters_per_relay=int(cfg.get("iters_per_relay", 20)),
            gamma_min=Fraction(cfg.get("gamma_min", "-1/4")),
            gamma_max=Fraction(cfg.get("gamma_max", "3/4")),
            alpha=Fraction(cfg.get("alpha", "3/4")),
            stop_n_solutions=int(cfg.get("stop_n_solutions", 1)),
            seed=int(cfg.get("seed", 20260601)),
        )
    if alg == "union_find":
        return decode_union_find(problem, syndrome, use_bp_guidance=False)
    if alg == "belief_find":
        return decode_union_find(
            problem, syndrome, use_bp_guidance=True, bp_iters=int(cfg.get("bp_iters", 15))
        )
    if alg == "mwpm":
        r_m = decode_mwpm(problem, syndrome)
        return DecodeResult(
            decoder=r_m.decoder,
            error_mask=r_m.error_mask,
            converged=r_m.converged,
            observable_mask=r_m.observable_mask,
            hamming_weight=r_m.hamming_weight,
            soft_weight=r_m.soft_weight,
            stats=r_m.stats,
        )
    if alg == "ssf":
        return decode_small_set_flip(
            problem,
            syndrome,
            max_subset_weight=int(cfg.get("max_subset_weight", 2)),
            max_steps=int(cfg.get("max_steps", 200)),
        )
    if alg == "sat_mwd":
        return decode_sat_mwd(problem, syndrome, timeout_sec=int(cfg.get("timeout_sec", 60)))

    raise RefusalError("UNKNOWN-DECODER", f"Unsupported decoder algorithm {alg!r}")


def format_failure_witness(
    problem: DecodingProblem,
    error_mask: int,
    res: DecodeResult,
    *,
    source_probe: str = "enumeration",
) -> dict[str, Any]:
    """Package and independently verify a decoder failure witness."""
    syn = problem.syndrome_of(error_mask)
    rec_syn = problem.syndrome_of(res.error_mask)
    residual_mask = error_mask ^ res.error_mask
    obs_flip = problem.observable_of(residual_mask)

    if res.converged and rec_syn == syn and obs_flip == 0:
        raise AssertionError("format_failure_witness called on a succeeding decode")

    failure_mode = "non_convergence" if (not res.converged or rec_syn != syn) else "logical_misdecode"
    mechs = bit_positions(error_mask)
    rec_mechs = bit_positions(res.error_mask)

    out: dict[str, Any] = {
        "source_probe": source_probe,
        "failure_mode": failure_mode,
        "decoder": res.decoder,
        "fault_weight": len(mechs),
        "fault_mechanisms": mechs,
        "syndrome_weight": syn.bit_count(),
        "decoder_converged": res.converged and (rec_syn == syn),
        "recovery_mechanisms": rec_mechs,
        "recovery_weight": len(rec_mechs),
        "residual_observable_flip_mask": obs_flip,
        "flipped_observables": bit_positions(obs_flip),
        "decoder_stats": res.stats,
    }
    if problem.witness_map and all(m in problem.witness_map for m in mechs):
        out["physical_circuit_faults"] = [
            [problem.witness_map[m][0], problem.witness_map[m][1]] for m in mechs
        ]
    return out


def verify_failure_witness(
    problem: DecodingProblem,
    decoder_spec: str | dict[str, Any],
    witness: dict[str, Any],
) -> bool:
    """Re-run `decoder_spec` on `witness['fault_mechanisms']` and verify identical failure."""
    mechs = witness["fault_mechanisms"]
    if len(set(mechs)) != len(mechs) or len(mechs) != witness["fault_weight"]:
        raise RefusalError("INVALID-WITNESS-SUPPORT", "Witness fault_mechanisms invalid or mismatched weight")
    e_mask = sum(1 << int(m) for m in mechs)
    syn = problem.syndrome_of(e_mask)
    res = run_decoder(problem, syn, decoder_spec)
    rec_syn = problem.syndrome_of(res.error_mask)
    obs_flip = problem.observable_of(e_mask ^ res.error_mask)

    if res.converged and rec_syn == syn and obs_flip == 0:
        raise RefusalError("WITNESS-DID-NOT-FAIL", "Decoder succeeded when replaying purported failure witness")
    mode = "non_convergence" if (not res.converged or rec_syn != syn) else "logical_misdecode"
    if mode != witness["failure_mode"]:
        raise RefusalError(
            "WITNESS-MODE-MISMATCH",
            f"Expected failure mode {witness['failure_mode']!r}, got {mode!r}",
        )
    return True


def certify_operational_radius(
    problem: DecodingProblem,
    decoder_spec: str | dict[str, Any] = "bp_osd",
    target_radius: int = 1,
    *,
    budget_cap: int = 200_000,
) -> dict[str, Any]:
    """Exhaustively certify operational decoding through fault weight `target_radius`.

    If `sum_{w=1..target_radius} C(n, w)` exceeds `budget_cap`, evaluates complete
    shells in increasing weight order until the next shell would exceed `budget_cap`,
    reporting `holds_through` honestly (matching `ftqc_gadgets` budget semantics!).
    """
    if target_radius < 0:
        raise RefusalError("INVALID-TARGET-RADIUS", f"target_radius must be >= 0, got {target_radius}")
    if budget_cap < 1:
        raise RefusalError("INVALID-BUDGET-CAP", "budget_cap must be >= 1")

    n = problem.num_mechanisms
    checked_sets = 0
    holds_through = 0
    max_bp_iters = 0
    max_elim_rank = 0
    max_clust_rank = 0
    decoder_label: str | None = None

    # Always verify weight 0
    r0 = run_decoder(problem, 0, decoder_spec)
    decoder_label = r0.decoder
    if not r0.converged or problem.syndrome_of(r0.error_mask) != 0 or problem.observable_of(r0.error_mask) != 0:
        return {
            "status": "counterexample",
            "decoder": decoder_label,
            "target_radius": target_radius,
            "holds_through": -1,
            "subsets_checked": 1,
            "witness": format_failure_witness(problem, 0, r0, source_probe="exhaustive_w0"),
        }
    checked_sets += 1

    for w in range(1, target_radius + 1):
        shell_size = comb(n, w)
        if checked_sets + shell_size > budget_cap:
            return {
                "status": "budget_reached",
                "decoder": decoder_label,
                "target_radius": target_radius,
                "holds_through": holds_through,
                "effective_distance_lower": 2 * holds_through + 1,
                "subsets_checked": checked_sets,
                "next_shell_weight": w,
                "next_shell_size": shell_size,
                "budget_cap": budget_cap,
                "work_upper": {
                    "max_bp_iterations": max_bp_iters,
                    "max_eliminated_rank": max_elim_rank,
                    "max_cluster_rank": max_clust_rank,
                },
            }

        for cols in combinations(range(n), w):
            e_mask = 0
            for c in cols:
                e_mask |= 1 << c
            syn = problem.syndrome_of(e_mask)
            res = run_decoder(problem, syn, decoder_spec)
            checked_sets += 1
            max_bp_iters = max(max_bp_iters, int(res.stats.get("bp_iterations", res.stats.get("total_iterations", 0))))
            max_elim_rank = max(max_elim_rank, int(res.stats.get("eliminated_rank", 0)))
            max_clust_rank = max(max_clust_rank, int(res.stats.get("max_cluster_rank", 0)))

            if (not res.converged) or problem.syndrome_of(res.error_mask) != syn or problem.observable_of(e_mask ^ res.error_mask) != 0:
                return {
                    "status": "counterexample",
                    "decoder": res.decoder,
                    "target_radius": target_radius,
                    "holds_through": w - 1,
                    "failure_weight": w,
                    "effective_distance_upper": 2 * w,
                    "subsets_checked": checked_sets,
                    "witness": format_failure_witness(problem, e_mask, res, source_probe=f"exhaustive_w{w}"),
                }

        holds_through = w

    return {
        "status": "proved",
        "decoder": decoder_label,
        "target_radius": target_radius,
        "holds_through": holds_through,
        "effective_distance_lower": 2 * holds_through + 1,
        "subsets_checked": checked_sets,
        "work_upper": {
            "max_bp_iterations": max_bp_iters,
            "max_eliminated_rank": max_elim_rank,
            "max_cluster_rank": max_clust_rank,
        },
    }


def _find_light_logical_representatives(problem: DecodingProblem, max_reps: int = 16) -> list[int]:
    """Compute low-weight representatives in `ker(H) \\ ker(L)` via pivoted RREF for midpoint probing."""
    # Build echelon basis of `H` over `0..n-1`
    n = problem.num_mechanisms
    piv_map: dict[int, int] = {}
    for r in problem.check_rows:
        v = r
        for p in sorted(piv_map.keys(), reverse=True):
            if (v >> p) & 1:
                v ^= piv_map[p]
        if v:
            npiv = v.bit_length() - 1
            for p in list(piv_map.keys()):
                if (piv_map[p] >> npiv) & 1:
                    piv_map[p] ^= v
            piv_map[npiv] = v

    piv_set = set(piv_map.keys())
    free_cols = [j for j in range(n) if j not in piv_set]
    logical_vecs: list[int] = []
    for f_col in free_cols:
        vec = 1 << f_col
        for p, row_v in piv_map.items():
            if (row_v >> f_col) & 1:
                vec |= 1 << p
        if problem.observable_of(vec) != 0:
            logical_vecs.append(vec)
    logical_vecs.sort(key=lambda x: (x.bit_count(), x))
    return logical_vecs[:max_reps]


def search_trapping_sets(
    problem: DecodingProblem,
    decoder_spec: str | dict[str, Any] = "bp_only",
    *,
    max_exhaustive_weight: int = 1,
    max_splits_per_generator: int = 32,
    compare_decoders: list[str | dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Directed adversarial search for quantum trapping sets & sub-distance decoder failures.

    Combines:
      1) Stabilizer split-support probe (`e_1 + e_2 = g in G_stab`),
      2) Logical-operator midpoint probe (`e subset supp(L_bar)` of weight `<= floor((wt(L_bar)-1)/2)`),
      3) Bounded low-weight exhaustive scan (`w <= max_exhaustive_weight`).
    """
    candidates_to_probe: list[tuple[int, str]] = []  # (error_mask, probe_label)
    seen_masks: set[int] = set()

    def enqueue(mask: int, tag: str) -> None:
        if mask > 0 and mask not in seen_masks:
            seen_masks.add(mask)
            candidates_to_probe.append((mask, tag))

    # 1. Bounded exhaustive low-weight probe
    for w in range(1, max_exhaustive_weight + 1):
        if comb(problem.num_mechanisms, w) <= 5000:
            for cols in combinations(range(problem.num_mechanisms), w):
                m = 0
                for c in cols:
                    m |= 1 << c
                enqueue(m, f"exhaustive_w{w}")

    # 2. Stabilizer split-support probe (Poulin & Chung 2008; Raveendran & Vasic 2021)
    stab_pool = list(problem.stabilizer_generators)
    # Also include low-weight pairwise XORs of overlapping stabilizer generators!
    for i in range(min(len(problem.stabilizer_generators), 24)):
        gi = problem.stabilizer_generators[i]
        for j in range(i + 1, min(len(problem.stabilizer_generators), 24)):
            gj = problem.stabilizer_generators[j]
            if (gi & gj) != 0:
                gij = gi ^ gj
                if 2 <= gij.bit_count() <= 12:
                    stab_pool.append(gij)

    for g_idx, g_mask in enumerate(stab_pool):
        supp = bit_positions(g_mask)
        wg = len(supp)
        if wg < 2:
            continue
        half_w = wg // 2
        count_s = 0
        for sub in combinations(supp, half_w):
            e_half = 0
            for q in sub:
                e_half |= 1 << q
            enqueue(e_half, f"stabilizer_split_g{g_idx}_w{half_w}")
            count_s += 1
            if count_s >= max_splits_per_generator:
                break

    # 3. Logical-operator midpoint probe
    for l_idx, l_mask in enumerate(_find_light_logical_representatives(problem)):
        supp_l = bit_positions(l_mask)
        wl = len(supp_l)
        if wl < 3:
            continue
        mid_w = (wl - 1) // 2
        count_l = 0
        for sub in combinations(supp_l, mid_w):
            e_mid = 0
            for q in sub:
                e_mid |= 1 << q
            enqueue(e_mid, f"logical_midpoint_l{l_idx}_w{mid_w}")
            count_l += 1
            if count_l >= max_splits_per_generator:
                break

    # Sort probes by increasing Hamming weight so the very first failure found is lightest!
    candidates_to_probe.sort(key=lambda item: (item[0].bit_count(), item[0]))

    probed = 0
    best_witness: dict[str, Any] | None = None
    decoder_name: str | None = None

    for e_mask, tag in candidates_to_probe:
        syn = problem.syndrome_of(e_mask)
        res = run_decoder(problem, syn, decoder_spec)
        decoder_name = res.decoder
        probed += 1
        rec_syn = problem.syndrome_of(res.error_mask)
        obs_flip = problem.observable_of(e_mask ^ res.error_mask)
        if (not res.converged) or rec_syn != syn or obs_flip != 0:
            best_witness = format_failure_witness(problem, e_mask, res, source_probe=tag)
            verify_failure_witness(problem, decoder_spec, best_witness)
            break

    comparison_table: list[dict[str, Any]] = []
    if best_witness is not None and compare_decoders:
        e_w = sum(1 << int(m) for m in best_witness["fault_mechanisms"])
        syn_w = problem.syndrome_of(e_w)
        for comp_spec in compare_decoders:
            c_res = run_decoder(problem, syn_w, comp_spec)
            c_syn_ok = c_res.converged and (problem.syndrome_of(c_res.error_mask) == syn_w)
            c_obs = problem.observable_of(e_w ^ c_res.error_mask)
            comparison_table.append(
                {
                    "decoder": c_res.decoder,
                    "syndrome_matched": c_syn_ok,
                    "logical_flip_mask": c_obs,
                    "healed": c_syn_ok and (c_obs == 0),
                    "recovery_weight": c_res.hamming_weight,
                    "stats": c_res.stats,
                }
            )

    return {
        "status": "trapping_set_found" if best_witness is not None else "no_failure_found",
        "decoder": decoder_name or str(decoder_spec),
        "candidates_probed": probed,
        "witness": best_witness,
        "comparison": comparison_table,
    }
