"""Search the signed Pauli orbit under supplied code actions; retain a replayable word."""
from collections import deque

from ftqc_codes import action_domain, automorphisms as aut, basis
from .algebra import Limit


def orbit(actions, seed, budget):
    if type(budget) is not int or budget < 1:
        raise Limit("signed measurement orbit requires a positive state budget")
    words = {seed: ()}
    queue = deque([seed])
    while queue:
        current = queue.popleft()
        for i, action in enumerate(actions):
            dest = aut.logical_image(action, current)
            if dest not in words:
                if len(words) >= budget:
                    raise Limit("signed measurement orbit exceeds budget")
                words[dest] = (*words[current], i)
                queue.append(dest)
    return words


def analyze(spec, budget):
    code = action_domain.load(spec["code"])
    paired = spec["basis"]
    basis.verify_basis_receipt(code, paired)
    actions = spec["actions"]
    for action in actions:
        aut.verify(code, action)
        if action["basis_sha256"] != paired["basis_sha256"]:
            raise ValueError("generator uses another logical basis")
    k = paired["k"]
    seed = basis.parse_pauli(spec["seed"], k)
    if not (seed[0] or seed[1]):
        raise ValueError("measurement seed must be nontrivial")
    physical = basis.parse_pauli(basis.embed_logical_pauli(paired, *seed), code.n)
    words = orbit(actions, seed, budget)
    settings = []
    group = basis._SignedStabilizerGroup(code.n, basis._code_signed_checks(code))
    for requested in spec["targets"]:
        target = basis.parse_pauli(requested, k)
        if target not in words:
            settings.append(dict(target=requested, reachable=False))
            continue
        steps = [s for i in words[target] for s in actions[i]["steps"]]
        image = aut.conjugate(physical, steps)
        desired = basis.parse_pauli(basis.embed_logical_pauli(paired, *target), code.n)
        x, z, phase = basis.mul_pauli(desired, image)
        if phase % 2 or not group.contains_exact((x, z, phase // 2)):
            raise ValueError("signed physical orbit disagrees with logical target")
        settings.append(dict(target=requested, reachable=True, word=list(words[target]),
                             physical=basis.pauli_from_xz(code.n, *image)))
    return dict(status="proved", basis_sha256=paired["basis_sha256"],
                orbit_size=len(words), settings=settings, scope="signed measurement settings")

