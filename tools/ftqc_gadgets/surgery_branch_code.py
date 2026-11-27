"""Finite CSS branching cylinders and X-basis unbranching for parallel Z products.

The brute-force tree uses the same restricted Tanner graph on each branch
(Cowtan et al., arXiv:2503.05003, Sec. 2 and App. D).
"""

from ftqc_codes.gf2 import Echelon, dot, supp
from .surgery_merge import word


def cylinders(n, hx, hz, targets, rounds):
    from .surgery import pauli
    from .surgery_geometry import check_rows

    if type(rounds) != int or rounds < 1 or not targets:
        raise ValueError("expected positive branching rounds and targets")
    check_rows(hx, n)
    check_rows(hz, n)
    logicals = [pauli(t, n) for t in targets]
    if any(v.x or not v.z or any(dot(v.z, h) for h in hx) for v in logicals):
        raise ValueError("branching requires commuting logical Z supports")
    sites = sorted(set(q for v in logicals for q in supp(v.z)))
    incident = [r for r, h in enumerate(hx) if any((h >> q) & 1 for q in sites)]
    # For each tree edge create an odd layer (one wire per incident X check)
    # and an even layer (one fresh copy of each vertex in the support).
    wire = n
    xrows = list(hx)
    zrows = list(hz)
    paths = {}
    copy_sites = []
    source = {q: q for q in sites}
    pending = [(source, {r: r for r in incident}, list(range(len(targets))), [])]
    leaves = {}
    while pending:
        parent, xids, indices, ancestors = pending.pop()
        groups = (
            [indices]
            if len(indices) == 1
            else [indices[: len(indices) // 2], indices[len(indices) // 2 :]]
        )
        for subset in groups:
            edge = {r: wire + i for i, r in enumerate(incident)}
            wire += len(edge)
            child = {q: wire + i for i, q in enumerate(sites)}
            wire += len(child)
            own = []
            for q in sites:
                row = (
                    (1 << parent[q])
                    | (1 << child[q])
                    | sum(1 << edge[r] for r in incident if hx[r] >> q & 1)
                )
                zrows.append(row)
                own.append(len(zrows) - 1)
                copy_sites.append((child[q], q))
            nxt = {}
            for r in incident:
                xrows[xids[r]] |= 1 << edge[r]
                nxt[r] = len(xrows)
                xrows.append(
                    (1 << edge[r]) | sum(1 << child[q] for q in sites if hx[r] >> q & 1)
                )
            path = [*ancestors, (child, edge, own)]
            if len(subset) == 1:
                leaves[subset[0]] = child
                paths[subset[0]] = path
            else:
                pending.append((child, nxt, subset, path))
    if wire - Echelon(xrows).rank - Echelon(zrows).rank != n - Echelon(
        hx
    ).rank - Echelon(hz).rank or any(dot(x, z) for x in xrows for z in zrows):
        raise ValueError("branch code does not preserve the memory dimension")
    eq = Echelon(zrows)
    leaf_masks = [
        sum(1 << leaves[j][q] for q in supp(v.z)) for j, v in enumerate(logicals)
    ]
    for v, leaf in zip(logicals, leaf_masks):
        if (v.z ^ leaf) not in eq:
            raise ValueError(
                "leaf representative is not equivalent to the requested parity"
            )
    ticks = [[["RX", q] for q in range(n, wire)]]
    nr = 0
    detectors = []

    def check(row, axis, meter):
        nonlocal nr
        x = axis == "X"
        ticks.append([["RX" if x else "R", meter]])
        for q in supp(row):
            ticks.append([["CX", meter, q] if x else ["CX", q, meter]])
        ticks.append([["MX" if x else "M", meter]])
        nr += 1
        return nr - 1

    for t in range(rounds):
        zs = [check(r, "Z", wire) for r in zrows]
        xs = [check(r, "X", wire) for r in xrows]
        if t == 0:
            first = zs
            detectors.extend([r] for r in zs[: len(hz)] + xs)
        else:
            detectors.extend([[a, b] for a, b in zip(prior, zs + xs)])
        prior = zs + xs
    # Each leaf reports its representative in the measured branching frame.
    # The products of interface checks on a root-to-leaf path relate it to the root.
    reads = [dict(records=[], sign=v.s) for v in logicals]
    meters = list(range(wire, wire + len(logicals)))
    for t in range(rounds):
        ticks.append([["R", q] for q in meters])
        supports = [supp(mask) for mask in leaf_masks]
        for k in range(max(map(len, supports))):
            ticks.append(
                [["CX", qs[k], m] for qs, m in zip(supports, meters) if k < len(qs)]
            )
        ticks.append([["M", m] for m in meters])
        for j, v in enumerate(logicals):
            if t:
                detectors.append([public[j], nr + j])
            else:
                reads[j]["records"] = [nr + j] + [
                    first[rows[i]]
                    for _, _, rows in paths[j]
                    for i, q in enumerate(sites)
                    if v.z >> q & 1
                ]
        if t == 0:
            public = list(range(nr, nr + len(meters)))
        nr += len(meters)
    ticks.append([["MX", q] for q in range(n, wire)])

    def split(q):
        return nr + q - n

    # X-only checks on the new patch give deterministic split detectors.
    for row in xrows[len(hx) :]:
        detectors.append([split(q) for q in supp(row)])
    for q, root in copy_sites:
        ticks.append([["Z_IF", root, split(q)]])
    spec = dict(
        n=n,
        checks=[word(n, r, "X") for r in hx] + [word(n, r, "Z") for r in hz],
        targets=targets,
        ticks=ticks,
        readouts=reads,
        detectors=detectors,
        terminal_syndrome=True,
    )
    return dict(
        spec=spec,
        geometry=dict(
            method="branch-code",
            n=wire,
            hx=xrows,
            hz=zrows,
            leaf_representatives=leaf_masks,
        ),
    )
