"""Parallel, basis-compatible Pauli readouts via coherent repetition branches."""

from ftqc_codes.gf2 import supp
from .surgery import pauli


def branch(n, checks, targets, rounds=1):
    """Separate shared sites using a reversible fanout tree in the local Pauli basis.

    Private leaves are disentangled and checked before release. Physical faults
    on the tree belong to this instrument, including its inverse.
    """
    if type(rounds) != int or rounds < 1 or not targets:
        raise ValueError("expected targets and positive rounds")
    axes = ["I"] * n
    consumers = [[] for _ in range(n)]
    for j, target in enumerate(targets):
        p = pauli(target, n)
        for q in supp(p.x | p.z):
            a = target[q + 1]
            if axes[q] not in ("I", a):
                raise ValueError("incompatible single-site Pauli bases")
            axes[q] = a
            consumers[q].append(j)
    ticks = []

    def layer(ops):
        if ops:
            ticks.append(ops)

    phase = [q for q, a in enumerate(axes) if a == "Y"]
    rot = [q for q, a in enumerate(axes) if a in ("X", "Y")]
    layer([["SDG", q] for q in phase])
    layer([["H", q] for q in rot])
    leaves = []
    wire = n
    for q, users in enumerate(consumers):
        extra = max(0, len(users) - 1)
        leaves.append([q, *range(wire, wire + extra)])
        wire += extra
    layer([["R", q] for q in range(n, wire)])
    tree = []
    step = 1
    while any(len(group) > step for group in leaves):
        row = [
            ["CX", group[j], group[step + j]]
            for group in leaves
            for j in range(min(step, len(group) - step))
        ]
        tree.append(row)
        layer(row)
        step *= 2
    supports = [[] for _ in targets]
    for users, group in zip(consumers, leaves):
        for j, q in zip(users, group):
            supports[j].append(q)
    if any(not group for group in supports):
        raise ValueError("identity target")
    meters = list(range(wire, wire + len(targets)))
    detectors = []
    for r in range(rounds):
        layer([["R", m] for m in meters])
        for i in range(max(map(len, supports))):
            layer(
                [
                    ["CX", group[i], m]
                    for group, m in zip(supports, meters)
                    if i < len(group)
                ]
            )
        layer([["M", m] for m in meters])
        if r:
            detectors.extend([[j, j + r * len(meters)] for j in range(len(meters))])
    for row in reversed(tree):
        layer(row)
    layer([["M", q] for q in range(n, wire)])
    detectors.extend([[rounds * len(meters) + j] for j in range(wire - n)])
    layer([["H", q] for q in rot])
    layer([["S", q] for q in phase])
    return dict(
        n=n,
        checks=checks,
        targets=targets,
        ticks=ticks,
        detectors=detectors,
        readouts=[
            dict(records=[j], sign=int(w[0] == "-")) for j, w in enumerate(targets)
        ],
        terminal_syndrome=True,
    )
