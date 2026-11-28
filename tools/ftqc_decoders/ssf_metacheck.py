"""Small-Set-Flip (`SSF`) single-shot HGP decoder & two-stage Metacheck syndrome repair.

Literature & Algorithmic Lineage:
  1. Small-Set-Flip (`SSF`) linear-time decoder for quantum expander / HGP codes:
     A. Leverrier, J.-P. Tillich, and G. Zemor, "Quantum Expander Codes",
     IEEE 56th FOCS, pp. 810-824 (2015), arXiv:1504.00822.
  2. Single-shot HGP fault-tolerance & residual weight proportionality:
     O. Fawzi, A. Grospellier, and A. Leverrier, "Constant overhead quantum
     fault-tolerance with quantum expander codes", IEEE 59th FOCS (2018) /
     Commun. ACM 64(1), 106-114 (2020), arXiv:1808.03821.
  3. Metacheck single-shot syndrome repair + data recovery:
     E. T. Campbell, "A theory of single-shot error correction for adversarial
     noise", Quantum Sci. Technol. 4, 025006 (2019), arXiv:1805.09271;
     A. Cowtan, Z. He, D. J. Williamson, and T. J. Yoder, "Parallel Logical
     Measurements via Quantum Code Surgery", arXiv:2510.14895
     (matching `ftqc_extractors/surgery_sequence.py` mapping-cone metachecks!).

Zero `float` or `complex` literals/calls appear anywhere in this module.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from typing import Any

from .bp_osd import DecodeResult, decode_bp_osd
from .decimal_bounds import RefusalError
from .dem import DecodingProblem, bit_positions
from .relay_bp import decode_relay_bp
from .union_find import decode_union_find


def decode_small_set_flip(
    problem: DecodingProblem,
    syndrome: int,
    *,
    max_subset_weight: int = 2,
    max_steps: int = 200,
) -> DecodeResult:
    """Decode `syndrome` via Leverrier-Tillich-Zemor Small-Set-Flip (`SSF`).

    Searches singletons `{j}` (Sipser-Spielman bit-flip) together with nonempty subsets
    `F subseteq supp(g)` of size `<= max_subset_weight` drawn from `problem.stabilizer_generators`
    (the dual CSS generators `H_X` when decoding `X` errors from `H_Z` syndromes),
    greedily applying the flip `F` maximizing exact rational relative syndrome drop
    `(|s| - |s ^ H F|) / |F| > 0`.
    """
    if syndrome < 0 or (syndrome >> problem.num_checks):
        raise RefusalError("SYNDROME-OUT-OF-RANGE", f"Syndrome exceeds {problem.num_checks} checks")
    if max_subset_weight < 1:
        raise RefusalError("INVALID-SSF-WEIGHT", "max_subset_weight must be >= 1")

    if syndrome == 0:
        return DecodeResult(
            decoder="ssf",
            error_mask=0,
            converged=True,
            observable_mask=0,
            hamming_weight=0,
            soft_weight=0,
            stats={"flip_steps": 0, "residual_syndrome_weight": 0, "candidate_small_sets": 0},
        )

    # Precompute candidate `(flip_mask, flip_size, induced_syndrome)` small sets!
    seen_masks: set[int] = set()
    candidates: list[tuple[int, int, int]] = []

    # 1. All singletons (size 1)
    for j in range(problem.num_mechanisms):
        m_mask = 1 << j
        ds = problem.mech_detectors[j]
        if ds != 0:
            seen_masks.add(m_mask)
            candidates.append((m_mask, 1, ds))

    # 2. Subsets of size `2 .. max_subset_weight` inside each dual generator support `g`
    if max_subset_weight >= 2 and problem.stabilizer_generators:
        for g in problem.stabilizer_generators:
            supp_g = bit_positions(g)
            max_w = min(max_subset_weight, len(supp_g))
            for w in range(2, max_w + 1):
                for sub in combinations(supp_g, w):
                    f_mask = 0
                    f_syn = 0
                    for q in sub:
                        f_mask |= 1 << q
                        f_syn ^= problem.mech_detectors[q]
                    if f_syn != 0 and f_mask not in seen_masks:
                        seen_masks.add(f_mask)
                        candidates.append((f_mask, w, f_syn))

    cur_syn = syndrome
    cur_wt = cur_syn.bit_count()
    e_hat = 0
    steps = 0

    for _step in range(max_steps):
        if cur_wt == 0:
            break
        best_f = 0
        best_fsyn = 0
        best_gain = Fraction(0, 1)
        best_new_wt = cur_wt

        for f_mask, f_sz, f_syn in candidates:
            # Skip sets disjoint from `cur_syn` (`(cur_syn & f_syn) == 0` cannot reduce weight!)
            if (cur_syn & f_syn) == 0:
                continue
            new_wt = (cur_syn ^ f_syn).bit_count()
            drop = cur_wt - new_wt
            if drop > 0:
                gain = Fraction(drop, f_sz)
                if gain > best_gain:
                    best_gain = gain
                    best_f = f_mask
                    best_fsyn = f_syn
                    best_new_wt = new_wt

        if best_f == 0:
            break
        e_hat ^= best_f
        cur_syn ^= best_fsyn
        cur_wt = best_new_wt
        steps += 1

    return DecodeResult(
        decoder="ssf",
        error_mask=e_hat,
        converged=(cur_syn == 0),
        observable_mask=problem.observable_of(e_hat),
        hamming_weight=e_hat.bit_count(),
        soft_weight=problem.soft_weight_of(e_hat),
        stats={
            "flip_steps": steps,
            "residual_syndrome_weight": cur_wt,
            "candidate_small_sets": len(candidates),
        },
    )


def decode_with_metachecks(
    problem: DecodingProblem,
    noisy_syndrome: int,
    *,
    syndrome_prior: str | Fraction = "1/1000",
    meta_decoder: str = "bp_osd",
    data_decoder: str = "bp_osd",
) -> dict[str, Any]:
    """Two-stage single-shot / surgery-cone metacheck syndrome repair + data recovery
    (Campbell 2019, arXiv:1805.09271; Cowtan, He, Williamson & Yoder, arXiv:2510.14895).

    Stage 1: Given `noisy_syndrome = (H * e_data) ^ e_meas`, computes metacheck defect
             `m_syn = M * noisy_syndrome = M * e_meas` (since `M * H == 0 mod 2`), and
             decodes `e_meas_hat` on `H_meta = M` so `repaired_syndrome = noisy_syndrome ^ e_meas_hat`
             satisfies `M * repaired_syndrome == 0`.
    Stage 2: Decodes `e_data_hat` on `(H, L)` against `repaired_syndrome` (or if `repaired_syndrome`
             remains outside `colspace(H)`, reports residual syndrome weight cleanly).
    """
    if not problem.metacheck_rows:
        raise RefusalError("MISSING-METACHECKS", f"Problem {problem.name!r} declares no metacheck_rows")
    if noisy_syndrome < 0 or (noisy_syndrome >> problem.num_checks):
        raise RefusalError("SYNDROME-OUT-OF-RANGE", f"Syndrome exceeds {problem.num_checks} checks")

    meta_prob = DecodingProblem.from_matrices(
        problem.metacheck_rows,
        (),
        problem.num_checks,
        name=f"{problem.name}:metachecks",
        default_prior=syndrome_prior,
    )
    m_syn = problem.metacheck_syndrome_of(noisy_syndrome)

    def _run_named(prob_inst: DecodingProblem, syn_val: int, kind: str) -> DecodeResult:
        if kind == "bp_osd":
            return decode_bp_osd(prob_inst, syn_val, post_method="osd", osd_order=1)
        if kind == "bp_lsd":
            return decode_bp_osd(prob_inst, syn_val, post_method="lsd", osd_order=1)
        if kind == "relay_bp":
            return decode_relay_bp(prob_inst, syn_val)
        if kind == "union_find":
            return decode_union_find(prob_inst, syn_val, use_bp_guidance=True)
        if kind == "ssf":
            return decode_small_set_flip(prob_inst, syn_val)
        raise RefusalError("UNKNOWN-DECODER", f"Unsupported metacheck stage decoder {kind!r}")

    meta_res = _run_named(meta_prob, m_syn, meta_decoder)
    repaired_syndrome = noisy_syndrome ^ meta_res.error_mask
    repaired_meta_defect = problem.metacheck_syndrome_of(repaired_syndrome)

    data_res = _run_named(problem, repaired_syndrome, data_decoder)
    residual_syndrome = repaired_syndrome ^ problem.syndrome_of(data_res.error_mask)

    return {
        "meta_decoder": meta_res.decoder,
        "data_decoder": data_res.decoder,
        "metacheck_defect_weight": m_syn.bit_count(),
        "meta_converged": meta_res.converged and (repaired_meta_defect == 0),
        "syndrome_repair_mask": meta_res.error_mask,
        "syndrome_repair_weight": meta_res.hamming_weight,
        "repaired_syndrome": repaired_syndrome,
        "data_converged": data_res.converged,
        "data_error_mask": data_res.error_mask,
        "data_hamming_weight": data_res.hamming_weight,
        "observable_mask": data_res.observable_mask,
        "residual_syndrome_weight": residual_syndrome.bit_count(),
    }
