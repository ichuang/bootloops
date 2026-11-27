"""Zero-float outward-rounded Logical Error Rate (`LER`) enclosures for `ftqc_decoders`.

Solves the global union-bound vacuity barrier (`C(L, t+1) * p^(t+1) >= 1` on large
qLDPC / surface-code schedules) via three mathematically distinct, clearly labeled
modes in `ftqc_decoders-receipt-1`:

  1. `exact_shells` (`evidence_kind = 'deterministic_exact_shells'`):
     Exhaustively evaluates decoder failures `F_w` across weight shells
     `w = verified_radius + 1 .. max_exact_weight` and bounds the remaining binomial
     tail `1 - sum_{w=0..W} C(n, w) p^w (1 - p)^(n - w)` in exact `fractions.Fraction`,
     emitting outward-rounded `Decimal` strings `[ler_lo, ler_hi]` (`ROUND_FLOOR` /
     `ROUND_CEILING`) with runtime containment assertions!

  2. `peierls_cluster` (`evidence_kind = 'analytical_peierls_cluster_bound'`):
     Applies the connected-syndrome-cluster percolation bound of:
       * A. A. Kovalev and L. P. Pryadko, Phys. Rev. A 87, 020304(R) (2013), arXiv:1208.2317;
       * D. Gottesman, Quantum Inf. Comput. 14, 1338-1372 (2014), arXiv:1310.2984.
     When disjoint clusters in the fault-interaction graph (max degree `D`) decode
     independently and every connected component of size `<= t` is corrected, a failure
     requires at least one connected component of size `s >= t + 1`. Using the strict
     rational upper bound `e < 27183 / 10000`, the number of size-`s` connected animals
     is at most `n * (e_hi * D)^(s - 1)`, yielding for `e_hi * D * p < 1`:
         `LER <= n * (e_hi * D)^t * p^(t + 1) / (1 - e_hi * D * p)`,
     giving non-vacuous (`<< 1`) analytical upper bounds at `n > 10^4` where `C(n, t+1) p^(t+1) = 1`!

  3. `stratified_shells` (`evidence_kind = 'empirical_stratified_hoeffding'`):
     Rare-event weight-shell importance splitting (S. Bravyi and A. Vargo, Phys. Rev. A
     88, 062308, 2013, arXiv:1308.6270) with deterministic SHA-256 shell sampling and
     one-sided Hoeffding confidence margins (`ROUND_CEILING` via `Decimal.ln` / `sqrt`),
     never conflated with deterministic proofs.

Zero `float` or `complex` literals/calls appear anywhere in this module.
"""

from __future__ import annotations

from decimal import Context, Decimal, ROUND_CEILING
from fractions import Fraction
import hashlib
from itertools import combinations
from math import comb
from typing import Any

from .decimal_bounds import (
    RefusalError,
    decimal_str_upper,
    fraction_to_outward_strs,
    parse_exact_probability,
)
from .dem import DecodingProblem
from .effective_distance import run_decoder

# Strict rational upper bound on Euler's constant: `e = 2.718281828... < 27183 / 10000`
EULER_E_UPPER = Fraction(27183, 10000)


def _binomial_shell_mass(n: int, w: int, p: Fraction) -> Fraction:
    return comb(n, w) * (p**w) * ((Fraction(1, 1) - p) ** (n - w))


def _tail_probability_exact_or_upper(n: int, max_w: int, p: Fraction) -> tuple[Fraction, str]:
    """Return `(tail_upper_fraction, method)` bounding `Pr[Binomial(n, p) > max_w]`."""
    if max_w >= n:
        return Fraction(0, 1), "zero_above_n"
    if n <= 1500:
        kept = Fraction(0, 1)
        q = Fraction(1, 1) - p
        # Compute `(1 - p)^n` once and step via exact ratio `p / (1 - p)` when `p < 1`
        if q == 0:
            return Fraction(1, 1), "exact_complement"
        term = q**n
        kept += term
        ratio = p / q
        for w in range(1, max_w + 1):
            term = term * Fraction(n - w + 1, w) * ratio
            kept += term
        tail = Fraction(1, 1) - kept
        if tail < 0:
            tail = Fraction(0, 1)
        return tail, "exact_binomial_complement"

    # Large-`n` fast rational union bound `Pr[X >= max_w + 1] <= min(1, C(n, max_w + 1) * p^(max_w + 1))`
    ub = min(Fraction(1, 1), comb(n, max_w + 1) * (p ** (max_w + 1)))
    return ub, "markov_subset_union_upper"


def enclose_exact_shells(
    problem: DecodingProblem,
    decoder_spec: str | dict[str, Any],
    physical_p: str | int | Fraction,
    *,
    verified_radius: int = 0,
    max_exact_weight: int = 2,
    budget_cap: int = 200_000,
    digits: int = 30,
) -> dict[str, Any]:
    """Deterministic `[ler_lo, ler_hi]` enclosure by exact shell enumeration + exact tail."""
    p = parse_exact_probability(physical_p, allow_zero=True, allow_one=False, max_half=True)
    if verified_radius < 0 or max_exact_weight < verified_radius:
        raise RefusalError("INVALID-SHELL-RANGE", "Require 0 <= verified_radius <= max_exact_weight")

    n = problem.num_mechanisms
    max_w = min(max_exact_weight, n)
    total_subsets = sum(comb(n, w) for w in range(verified_radius + 1, max_w + 1))
    if total_subsets > budget_cap:
        raise RefusalError(
            "SHELL-BUDGET-EXCEEDED",
            f"Exact shells {verified_radius + 1}..{max_w} on n={n} require {total_subsets} > budget_cap={budget_cap}",
        )

    shell_reports: list[dict[str, Any]] = []
    exact_fail_prob = Fraction(0, 1)

    for w in range(verified_radius + 1, max_w + 1):
        fails = 0
        total_w = comb(n, w)
        for cols in combinations(range(n), w):
            e_mask = 0
            for c in cols:
                e_mask |= 1 << c
            syn = problem.syndrome_of(e_mask)
            res = run_decoder(problem, syn, decoder_spec)
            if (not res.converged) or problem.syndrome_of(res.error_mask) != syn or problem.observable_of(e_mask ^ res.error_mask) != 0:
                fails += 1
        mass_w = _binomial_shell_mass(n, w, p)
        contrib = Fraction(fails, total_w) * mass_w if total_w > 0 else Fraction(0, 1)
        exact_fail_prob += contrib
        shell_reports.append(
            {
                "weight": w,
                "subsets": total_w,
                "failures": fails,
                "conditional_failure_fraction": f"{fails}/{total_w}",
            }
        )

    tail_hi, tail_method = _tail_probability_exact_or_upper(n, max_w, p)
    ler_lo_frac = exact_fail_prob
    ler_hi_frac = min(Fraction(1, 1), exact_fail_prob + tail_hi)

    lo_str, _ = fraction_to_outward_strs(ler_lo_frac, digits=digits)
    _, hi_str = fraction_to_outward_strs(ler_hi_frac, digits=digits)
    if Fraction(lo_str) > Fraction(hi_str):
        raise AssertionError("Outward interval inversion in enclose_exact_shells")

    return {
        "mode": "exact_shells",
        "evidence_kind": "deterministic_exact_shells",
        "physical_p": f"{p.numerator}/{p.denominator}",
        "verified_radius": verified_radius,
        "max_exact_weight": max_w,
        "tail_method": tail_method,
        "shells": shell_reports,
        "ler_lo": lo_str,
        "ler_hi": hi_str,
        "naive_global_union_hi": fraction_to_outward_strs(
            min(Fraction(1, 1), comb(n, verified_radius + 1) * (p ** (verified_radius + 1))),
            digits=digits,
        )[1],
    }


def enclose_peierls_cluster(
    problem: DecodingProblem,
    physical_p: str | int | Fraction,
    *,
    verified_component_radius: int,
    cluster_degree_override: int | None = None,
    digits: int = 30,
) -> dict[str, Any]:
    """Analytical connected-cluster Peierls upper bound (Kovalev & Pryadko 2013; Gottesman 2014).

    Bounds the probability that an independent Bernoulli(`p`) fault configuration on
    `problem` contains at least one connected component of size `>= t + 1` in the
    syndrome-interaction graph (max degree `D`), compared directly against the naive
    global union bound `min(1, C(n, t+1) * p^(t+1))`.
    """
    p = parse_exact_probability(physical_p, allow_zero=True, allow_one=False, max_half=True)
    if verified_component_radius < 0:
        raise RefusalError("INVALID-COMPONENT-RADIUS", "verified_component_radius must be >= 0")

    t = verified_component_radius
    n = problem.num_mechanisms
    deg = cluster_degree_override if cluster_degree_override is not None else max(1, problem.max_interaction_degree())
    if deg < 1:
        raise RefusalError("INVALID-CLUSTER-DEGREE", f"Cluster graph degree must be >= 1, got {deg}")

    naive_union = min(Fraction(1, 1), comb(n, t + 1) * (p ** (t + 1)))
    _, naive_hi_str = fraction_to_outward_strs(naive_union, digits=digits)

    if p == 0:
        return {
            "mode": "peierls_cluster",
            "evidence_kind": "analytical_peierls_cluster_bound",
            "physical_p": "0/1",
            "verified_component_radius": t,
            "interaction_degree": deg,
            "subthreshold_condition_met": True,
            "ler_lo": "0",
            "ler_hi": "0",
            "naive_global_union_hi": "0",
        }

    branch_factor = EULER_E_UPPER * deg * p
    if branch_factor >= Fraction(1, 1):
        return {
            "mode": "peierls_cluster",
            "evidence_kind": "analytical_peierls_cluster_bound",
            "physical_p": f"{p.numerator}/{p.denominator}",
            "verified_component_radius": t,
            "interaction_degree": deg,
            "subthreshold_condition_met": False,
            "branch_factor": f"{branch_factor.numerator}/{branch_factor.denominator}",
            "ler_lo": "0",
            "ler_hi": "1",
            "naive_global_union_hi": naive_hi_str,
        }

    # Exact rational series sum: `n * (e_hi * D)^t * p^(t + 1) / (1 - e_hi * D * p)`
    num_series = Fraction(n, 1) * ((EULER_E_UPPER * deg) ** t) * (p ** (t + 1))
    peierls_hi = min(Fraction(1, 1), num_series / (Fraction(1, 1) - branch_factor))
    _, hi_str = fraction_to_outward_strs(peierls_hi, digits=digits)

    return {
        "mode": "peierls_cluster",
        "evidence_kind": "analytical_peierls_cluster_bound",
        "assumption": "disjoint interaction-graph fault clusters decode independently and every component of size <= verified_component_radius succeeds",
        "physical_p": f"{p.numerator}/{p.denominator}",
        "verified_component_radius": t,
        "interaction_degree": deg,
        "subthreshold_condition_met": True,
        "branch_factor": f"{branch_factor.numerator}/{branch_factor.denominator}",
        "ler_lo": "0",
        "ler_hi": hi_str,
        "naive_global_union_hi": naive_hi_str,
    }


def _deterministic_subset_sample(n: int, w: int, seed: int, shot_idx: int) -> int:
    """Deterministic SHA-256 sampling of a uniform size-`w` subset bitmask from `0..n-1`."""
    chosen: set[int] = set()
    counter = 0
    while len(chosen) < w:
        digest = hashlib.sha256(f"shell:{seed}:{w}:{shot_idx}:{counter}".encode("ascii")).digest()
        for k in range(0, 32, 4):
            if len(chosen) >= w:
                break
            idx = int.from_bytes(digest[k : k + 4], "little") % n
            chosen.add(idx)
        counter += 1
    mask = 0
    for c in chosen:
        mask |= 1 << c
    return mask


def enclose_stratified_shells(
    problem: DecodingProblem,
    decoder_spec: str | dict[str, Any],
    physical_p: str | int | Fraction,
    *,
    verified_radius: int = 0,
    max_shell_weight: int = 4,
    shots_per_shell: int = 200,
    failure_delta: str | Fraction = "1/10000",
    seed: int = 20261126,
    digits: int = 30,
) -> dict[str, Any]:
    """Empirical rare-event weight-shell splitting (Bravyi & Vargo 2013) with outward-rounded
    one-sided Hoeffding confidence bounds (`evidence_kind = 'empirical_stratified_hoeffding'`).
    """
    p = parse_exact_probability(physical_p, allow_zero=False, allow_one=False, max_half=True)
    delta = parse_exact_probability(failure_delta, allow_zero=False, allow_one=False)
    if shots_per_shell < 10:
        raise RefusalError("TOO-FEW-SHOTS", "shots_per_shell must be >= 10")
    if verified_radius < 0 or max_shell_weight <= verified_radius:
        raise RefusalError("INVALID-SHELL-RANGE", "Require 0 <= verified_radius < max_shell_weight")

    n = problem.num_mechanisms
    max_w = min(max_shell_weight, n)
    num_shells = max_w - verified_radius

    # Union-bound `delta / num_shells` across shells; Hoeffding margin `sqrt(ln(num_shells / delta) / (2 * N))`
    # evaluated strictly upward (`ROUND_CEILING`) at high working precision with zero floats!
    ctx_hi = Context(prec=max(80, digits + 30), rounding=ROUND_CEILING)
    inv_delta_per_shell = Decimal(num_shells * delta.denominator) / Decimal(delta.numerator)
    ln_term = ctx_hi.ln(inv_delta_per_shell)
    eps_hoeffding = ctx_hi.sqrt(ctx_hi.divide(ln_term, Decimal(2 * shots_per_shell)))

    total_emp_d = Decimal(0)
    total_upper_d = Decimal(0)
    shell_rows: list[dict[str, Any]] = []

    for w in range(verified_radius + 1, max_w + 1):
        fails = 0
        for s_idx in range(shots_per_shell):
            e_mask = _deterministic_subset_sample(n, w, seed, s_idx)
            syn = problem.syndrome_of(e_mask)
            res = run_decoder(problem, syn, decoder_spec)
            if (not res.converged) or problem.syndrome_of(res.error_mask) != syn or problem.observable_of(e_mask ^ res.error_mask) != 0:
                fails += 1

        mass_w = _binomial_shell_mass(n, w, p)
        _, mass_hi_s = fraction_to_outward_strs(mass_w, digits=digits + 10)
        mass_hi_d = Decimal(mass_hi_s)

        emp_rate_d = ctx_hi.divide(Decimal(fails), Decimal(shots_per_shell))
        rate_hi_d = min(Decimal(1), ctx_hi.add(emp_rate_d, eps_hoeffding))
        total_emp_d = ctx_hi.add(total_emp_d, ctx_hi.multiply(mass_hi_d, emp_rate_d))
        total_upper_d = ctx_hi.add(total_upper_d, ctx_hi.multiply(mass_hi_d, rate_hi_d))
        shell_rows.append(
            {
                "weight": w,
                "shots": shots_per_shell,
                "failures": fails,
                "empirical_shell_rate": f"{fails}/{shots_per_shell}",
                "hoeffding_rate_hi": decimal_str_upper(rate_hi_d, digits=digits),
            }
        )

    tail_frac, tail_method = _tail_probability_exact_or_upper(n, max_w, p)
    _, tail_hi_s = fraction_to_outward_strs(tail_frac, digits=digits + 10)
    total_upper_d = min(Decimal(1), ctx_hi.add(total_upper_d, Decimal(tail_hi_s)))

    return {
        "mode": "stratified_shells",
        "evidence_kind": "empirical_stratified_hoeffding",
        "confidence_level": f"1 - {delta.numerator}/{delta.denominator}",
        "seed": seed,
        "physical_p": f"{p.numerator}/{p.denominator}",
        "verified_radius": verified_radius,
        "max_shell_weight": max_w,
        "shots_per_shell": shots_per_shell,
        "tail_method": tail_method,
        "shells": shell_rows,
        "empirical_point_estimate_hi": decimal_str_upper(min(Decimal(1), total_emp_d), digits=digits),
        "ler_hi": decimal_str_upper(total_upper_d, digits=digits),
    }
