"""Syndrome-equivariant decoder ensemble on original, weighted mechanism columns."""

from functools import lru_cache

from ftqc_codes.interchange import transport_word
from .bp_osd import DecodeResult
from .decimal_bounds import RefusalError


@lru_cache(maxsize=128)
def row_transform(rows, inverse):
    """Express pulled-back checks in the supplied row basis, including dependent rows."""
    pivots = {}

    def reduce(row, combination):
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                break
            other, how = pivots[pivot]
            row ^= other
            combination ^= how
        return row, combination

    for i, row in enumerate(rows):
        remaining, combination = reduce(row, 1 << i)
        if remaining:
            pivots[remaining.bit_length() - 1] = remaining, combination
    transform = []
    for row in rows:
        remaining, combination = reduce(transport_word(row, inverse), 0)
        if remaining:
            raise RefusalError(
                "INVALID-AUTOMORPHISM",
                "permutation changes the check or observable row space",
            )
        transform.append(combination)
    return tuple(transform)


def action_permutations(problem, actions):
    """Pull sector-preserving single-site data permutations out of signed CSS receipts."""
    from ftqc_codes import automorphisms, certify, codes

    if problem.source_meta.get("kind") != "css_code" or not isinstance(actions, list):
        raise RefusalError(
            "INVALID-AUTOMORPHISM", "code actions require a static CSS sector"
        )
    x_sector = problem.source_meta["sector"] == "X"
    checks, stabs = list(problem.check_rows), list(problem.stabilizer_generators)
    hx, hz = (stabs, checks) if x_sector else (checks, stabs)
    n = problem.num_mechanisms
    code = codes.CSSCode(problem.name, n, hx, hz)
    if code.sha256() != problem.source_meta.get("code_sha256"):
        raise RefusalError(
            "INVALID-AUTOMORPHISM", "decoder CSS coordinates differ from the code"
        )
    perms = []
    for action in actions:
        try:
            automorphisms.verify(code, action)
        except (certify.Refusal, ValueError) as exc:
            raise RefusalError("INVALID-AUTOMORPHISM", str(exc)) from exc
        perm = []
        for q in range(n):
            pauli = (1 << q, 0, 0) if x_sector else (0, 1 << q, 0)
            x, z, _ = automorphisms.conjugate(pauli, action["steps"])
            moved, other = (x, z) if x_sector else (z, x)
            if other or moved.bit_count() != 1:
                raise RefusalError(
                    "INVALID-AUTOMORPHISM",
                    "action does not permute single-site sector errors",
                )
            perm.append(moved.bit_length() - 1)
        perms.append(perm)
    return perms


def decode(problem, syndrome, cfg):
    from .effective_distance import run_decoder

    base = cfg.get("base", "bp_osd")
    if (
        base.get("algorithm") if isinstance(base, dict) else base
    ) == "automorphism_ensemble":
        raise RefusalError(
            "RECURSIVE-ENSEMBLE", "base decoder must be a single decoder"
        )
    n = problem.num_mechanisms
    identity = list(range(n))
    if ("permutations" in cfg) == ("actions" in cfg):
        raise RefusalError(
            "INVALID-AUTOMORPHISM", "supply permutations or signed code actions"
        )
    perms = (
        action_permutations(problem, cfg["actions"])
        if "actions" in cfg
        else cfg["permutations"]
    )
    if not isinstance(perms, list):
        raise RefusalError("INVALID-AUTOMORPHISM", "expected mechanism permutations")
    best = None
    seen = set()
    successes = 0
    for perm in [identity, *perms]:
        if (
            not isinstance(perm, list)
            or len(perm) != n
            or any(type(q) is not int for q in perm)
            or set(perm) != set(identity)
        ):
            raise RefusalError("INVALID-AUTOMORPHISM", "expected a mechanism bijection")
        if tuple(perm) in seen:
            continue
        seen.add(tuple(perm))
        inverse = [0] * n
        for i, j in enumerate(perm):
            if (
                problem.priors[i] != problem.priors[j]
                or problem.weights[i] != problem.weights[j]
            ):
                raise RefusalError(
                    "INVALID-AUTOMORPHISM",
                    "permutation changes a mechanism prior or weight",
                )
            inverse[j] = i
        pullback = tuple(inverse)
        transform = row_transform(problem.check_rows, pullback)
        row_transform(problem.check_rows + problem.observable_rows, pullback)
        moved = sum(
            ((row & syndrome).bit_count() % 2) << i for i, row in enumerate(transform)
        )
        result = run_decoder(problem, moved, base)
        correction = transport_word(result.error_mask, pullback)
        valid = result.converged and problem.syndrome_of(correction) == syndrome
        successes += int(valid)
        weight = problem.soft_weight_of(correction)
        score = (not valid, weight, correction.bit_count(), correction)
        if best is None or score < best[0]:
            best = (score, len(seen) - 1)
    (failed, weight, hamming, correction), member = best
    return DecodeResult(
        "automorphism_ensemble",
        correction,
        not failed,
        problem.observable_of(correction),
        hamming,
        weight,
        {
            "members": len(seen),
            "converged_members": successes,
            "selected_member": member,
        },
    )
