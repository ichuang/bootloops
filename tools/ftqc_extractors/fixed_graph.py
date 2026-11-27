"""Check finite fixed-coupler extractor hypotheses and synthesize a selected signed measurement."""

from itertools import combinations, product
from math import comb

from ftqc_audit.stab import Gen
from ftqc_codes.gf2 import Echelon
from .algebra import bits, kernel, parity, transpose, Limit
from .block_reading import combine, rows
from .extractor import incidence, matching


def pauli(word, n):
    if (
        not isinstance(word, str)
        or len(word) != n + 1
        or word[0] not in "+-"
        or set(word[1:]) - set("IXYZ")
    ):
        raise ValueError("expected signed Pauli, wire zero first")
    return Gen(
        sum(1 << i for i, p in enumerate(word[1:]) if p in "XY"),
        sum(1 << i for i, p in enumerate(word[1:]) if p in "ZY"),
        int(word[0] == "-"),
    )


def signed_basis(gens, n):
    packed = [g.x | (g.z << n) for g in gens]
    if any(
        parity(g.x, h.z) != parity(g.z, h.x)
        for i, g in enumerate(gens)
        for h in gens[:i]
    ):
        raise ValueError("anticommuting checks")
    for rel in kernel(transpose(packed, 2 * n), len(gens)):
        acc = Gen(0, 0)
        for i in bits(rel):
            acc.mul(gens[i])
        if acc.s:
            raise ValueError("inconsistent stabilizer signs")
    return Echelon(packed)


def base_witness(n, gens, basis, d, budget):
    for w in range(1, min(n + 1, d)):
        budget -= comb(n, w) * 3**w
        if budget < 0:
            raise Limit("signed stabilizer distance search")
        for support in combinations(range(n), w):
            for labels in product((1, 2, 3), repeat=w):
                x = sum((p & 1) << q for q, p in zip(support, labels))
                z = sum((p >> 1) << q for q, p in zip(support, labels))
                if x | (z << n) not in basis and all(
                    parity(x, g.z) == parity(z, g.x) for g in gens
                ):
                    return [x, z]


def analyze(s, budget):
    n, nv, edges = s["n"], s["vertices"], s["edges"]
    if type(n) is not int or n < 1 or type(nv) is not int or nv < 1:
        raise ValueError("invalid dimensions")
    inc = incidence(nv, edges)
    ne = len(edges)
    ports = s["ports"]
    if len(ports) != n or any(type(v) is not int or not -1 <= v < nv for v in ports):
        raise ValueError("invalid port mapping")
    active = [q for q, v in enumerate(ports) if v >= 0]
    if not active or len({ports[q] for q in active}) != len(active):
        raise ValueError("ports must be injective")
    cycles = rows(s["cycles"], ne)
    if (
        len(cycles) != ne - nv + 1
        or Echelon(cycles).rank != len(cycles)
        or any(parity(c, a) for c in cycles for a in inc)
    ):
        raise ValueError("cycle list is not a basis")
    gens = [pauli(w, n) for w in s["checks"]]
    basis = signed_basis(gens, n)
    pools = rows(s["routes"], ne)
    if len(pools) != len(gens):
        raise ValueError("one physical edge pool per check required")
    # Every even subset is a sum of pairs from a fixed root to the remaining terminals.
    paths = []
    for g, pool in zip(gens, pools):
        sites = [q for q in active if ((g.x | g.z) >> q) & 1]
        choices = list(bits(pool))
        routes = {}
        for q in sites[1:]:
            try:
                path = matching(
                    nv,
                    [edges[i] for i in choices],
                    (1 << ports[q]) ^ (1 << ports[sites[0]]),
                )
            except KeyError:
                return dict(
                    status="counterexample",
                    obligation="even-subset routing",
                    pair=[sites[0], q],
                )
            routes[q] = sum(1 << choices[i] for i in bits(path))
        paths.append(routes)
    # Check-side row masks specify every physical coupler of the fixed Tanner graph.
    wiring = [(g.x | g.z) | (pool << n) for g, pool in zip(gens, pools)]
    wiring += [
        (a << n) | sum(1 << q for q in active if ports[q] == v)
        for v, a in enumerate(inc)
    ]
    wiring += [c << n for c in cycles]
    if "couplers" in s and rows(s["couplers"], n + ne) != wiring:
        raise ValueError("coupler inventory differs from Definition 32")
    degree = max(v.bit_count() for v in [*wiring, *transpose(wiring, n + ne)])
    if type(s["max_degree"]) is not int or degree > s["max_degree"]:
        return dict(status="counterexample", obligation="coupler degree", degree=degree)
    d = s["target"]
    if nv - 1 >= budget.bit_length() or (1 << (nv - 1)) > budget:
        raise Limit("relative expansion cut budget")
    endpoints = transpose(inc, ne)
    portmask = sum(1 << ports[q] for q in active)
    for shore in range(1 << (nv - 1)):
        count = (shore & portmask).bit_count()
        cut = sum(parity(shore, e) for e in endpoints)
        if cut < min(d, count, len(active) - count):
            return dict(
                status="counterexample",
                obligation="relative expansion",
                shore=shore,
                cut=cut,
            )
    bad = base_witness(n, gens, basis, d, budget)
    if bad:
        return dict(status="counterexample", obligation="memory distance", logical=bad)
    L = pauli(s["logical"], n)
    if (L.x | (L.z << n)) in basis or any(
        parity(L.x, g.z) != parity(L.z, g.x) for g in gens
    ):
        raise ValueError("selected Pauli must be a nontrivial memory logical")
    if any(ports[q] < 0 for q in bits(L.x | L.z)):
        raise ValueError("logical representative extends beyond partial support")
    merged = []
    for g, r in zip(gens, paths):
        anti = (g.x & L.z) ^ (g.z & L.x)
        path = combine(
            list(r.values()), sum(1 << i for i, q in enumerate(r) if anti >> q & 1)
        )
        merged.append(Gen(g.x | (path << n), g.z, g.s))
    for v, a in enumerate(inc):
        m = sum(1 << q for q in active if ports[q] == v)
        merged.append(Gen(L.x & m, (L.z & m) | (a << n)))
    merged.extend(Gen(c << n, 0) for c in cycles)
    if any((g.x | g.z) & ~allowed for g, allowed in zip(merged, wiring, strict=True)):
        raise ValueError("selected setting requires an uninstalled coupler")
    if signed_basis(merged, n + ne).rank != basis.rank + ne + 1:
        raise ValueError("incorrect merged logical dimension")
    return dict(
        status="proved",
        model="fixed-coupler graph; phenomenological Theorem 11",
        lower=d,
        rounds=[d, d, d],
        support=active,
        full=len(active) == n,
        physical_check_couplers=wiring,
        max_degree=degree,
        merged_checks=[[g.x, g.z, g.s] for g in merged],
        readout_sign=L.s,
        split_X_paths=[
            matching(nv, edges, 1 ^ (1 << ports[q])) for q in bits(L.x | L.z)
        ],
        split_sites=list(bits(L.x | L.z)),
    )
