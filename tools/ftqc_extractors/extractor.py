"""Graph surgery and the canonical cyclic construction of arXiv:2606.03507v2."""

from functools import reduce
from operator import xor
from .algebra import (
    bits,
    cyclic,
    css_dimension,
    echelon,
    hgp,
    hgp_parameters,
    kernel,
    parity,
    reduce_row,
    seed,
    span,
    Limit,
)


def incidence(n, edges):
    rows = [0] * n
    seen = set()
    for e, (u, v) in enumerate(edges):
        if not 0 <= u < n or not 0 <= v < n or u == v or tuple(sorted((u, v))) in seen:
            raise ValueError("expected simple graph on the declared vertices")
        seen.add(tuple(sorted((u, v))))
        rows[u] |= 1 << e
        rows[v] |= 1 << e
    if len(echelon(rows)) != n - 1:
        raise ValueError("extractor graph disconnected")
    return rows


def matching(vertices, edges, boundary):
    """Return a tree T-join by expressing the desired boundary in incidence columns."""
    basis = {}
    for i, (u, v) in enumerate(edges):
        s, chain = (1 << u) ^ (1 << v), 1 << i
        while s and s.bit_length() in basis:
            b, c = basis[s.bit_length()]
            s, chain = s ^ b, chain ^ c
        if s:
            basis[s.bit_length()] = s, chain
    chain = 0
    while boundary:
        b, c = basis[boundary.bit_length()]
        boundary, chain = boundary ^ b, chain ^ c
    return chain


def column(spec):
    a, na = seed(spec["a"])
    b, nb = seed(spec["b"])
    j, c = spec["column"], spec["word"]
    if not 0 <= j < nb or not reduce_row(1 << j, echelon(b)):
        raise ValueError("column is trivial in the cokernel")
    if not 0 < c < 2**na or any(parity(c, row) for row in a):
        raise ValueError("word must be a nonzero seed codeword")
    graph = spec["edges"]
    nv = spec["vertices"]
    inc = incidence(nv, graph)
    ports = spec["ports"]
    if (
        len(ports) != na
        or sorted(ports[i] for i in bits(c)) != sorted(set(ports[i] for i in bits(c)))
        or any(not 0 <= p < nv for p in ports)
    ):
        raise ValueError("invalid single-column ports")
    chains = []
    for i, row in enumerate(a):
        ends = sum(1 << ports[u] for u in bits(row & c))
        if "paths" in spec:
            if len(spec["paths"]) != len(a):
                raise ValueError("path count")
            path = spec["paths"][i]
        else:
            path = matching(nv, graph, ends)
        if (
            not 0 <= path < 2 ** len(graph)
            or sum(parity(path, r) << v for v, r in enumerate(inc)) != ends
            or (not ends and path)
        ):
            raise ValueError("deformation boundary mismatch")
        chains.append(path)
    n, x, z = hgp(a, na, b, nb)
    vertex = [r << n for r in inc]
    for u in bits(c):
        vertex[ports[u]] |= 1 << (u * nb + j)
    cycles = kernel(inc, len(graph))
    merged = (
        n + len(graph),
        [s ^ ((chains[i // nb] << n) if i % nb == j else 0) for i, s in enumerate(x)]
        + [s << n for s in cycles],
        z + vertex,
    )
    if css_dimension(merged) != css_dimension((n, x, z)) - 1 or reduce(
        xor, vertex, 0
    ) != sum(1 << (u * nb + j) for u in bits(c)):
        raise ValueError("measurement lost or added a logical degree of freedom")
    sub_vertex = [r << na for r in inc]
    for u in bits(c):
        sub_vertex[ports[u]] |= 1 << u
    reduced = (
        na + len(graph),
        [r ^ (p << na) for r, p in zip(a, chains)] + [s << na for s in cycles],
        sub_vertex,
    )
    css_dimension(reduced)
    return merged, reduced


def full_graph(n, k, edges):
    def vertex(axis, copy, site):
        return (axis * k + copy) * n + site

    out = []
    for axis in range(2):
        for i in range(k):
            out.extend([vertex(axis, i, u), vertex(axis, i, v)] for u, v in edges)
            if i:
                out.extend(
                    [vertex(axis, i - 1, u), vertex(axis, i, u)] for u in range(n)
                )
    out.extend([vertex(0, j, i), vertex(1, i, j)] for j in range(k) for i in range(k))
    return out


def canonical(spec):
    """Theorem application to an implicitly specified canonical family, App. E.1/F."""
    n, p, edges = spec["n"], spec["polynomial"], spec["edges"]
    h = cyclic(n, p)
    words = span(kernel(h, n))
    k = len(words).bit_length() - 1
    if (
        not k
        or p.bit_length() - 1 != k
        or len({w % (1 << k) for w in words}) != len(words)
    ):
        raise ValueError("invalid cyclic check polynomial or information set")
    d = min(w.bit_count() for w in words if w)
    if n > 20 or 2**n * len(words) ** 2 > 200000000:
        raise Limit("canonical graph cut enumeration")
    inc = incidence(n, edges)
    # All vertex products modulo the all-vertex product. Complements exchange
    # c' and c'+c, both included below. Equation (26) handles the other columns.
    connectivity = len(edges)
    for shore in range(1 << (n - 1)):
        cut = reduce(xor, (inc[u] for u in bits(shore)), 0).bit_count()
        if shore:
            connectivity = min(connectivity, cut)
        if cut >= d:
            continue
        for c in words[1:]:
            for other in words[1:]:
                if c != other and (other ^ (shore & c)).bit_count() + cut < d:
                    return dict(
                        status="counterexample",
                        claim="single-column preservation",
                        word=c,
                        other=other,
                        shore=shore,
                        weight=(other ^ (shore & c)).bit_count() + cut,
                    )
    graph = full_graph(n, k, edges)
    if "assembly" in spec and spec["assembly"] != graph:
        raise ValueError("assembly differs from the canonical labeled graph")
    incidence(2 * n * k, graph)
    return dict(
        status="proved",
        model="canonical Appendix E.1/F phenomenological protocol",
        theorem="2606.03507v2 Theorem 22",
        lower=min(d, k * k, k * connectivity),
        base=hgp_parameters(h, n, h, n),
        seed_dimension=k,
        single_columns=len(words) - 1,
        connectivity=connectivity,
        vertices=2 * n * k,
        edges=len(graph),
        rounds_per_stage=d,
    )


def check_column(spec, budget, replay=False):
    """Use the capped Z reduction on a genuine single-column construction."""
    from itertools import combinations
    from math import comb
    from ftqc_codes import gf2
    from .algebra import distance_below

    merged, reduced = column(spec)
    a, na = seed(spec["a"])
    b, nb = seed(spec["b"])
    base = hgp_parameters(a, na, b, nb)
    threshold = spec["target"]
    result = dict(claim="merged CSS distance", base=base)
    if threshold > base["d"]:
        witness = distance_below(merged, threshold, budget)
        result["method"] = "direct merged-code enumeration"
    else:
        # An X logical projects to a surviving base X logical: an original
        # X stabilizer can be stripped using its path check and the remaining
        # ancilla cycle is a stabilizer. For Z, use Lemma 15 / Theorem 17.
        width, x, z = reduced
        result["method"] = (
            "2606.03507v2 Appendix F.1: capped Z subcode and X projection"
        )
        result["reduced_qubits"] = width
        witness = None
        zb = gf2.Echelon(z)
        dual = gf2.nullspace(z, width) if replay else []
        for w in range(1, min(threshold, width + 1)):
            budget -= comb(width, w)
            if budget < 0:
                raise Limit("reduced Z distance enumeration")
            for pos in combinations(range(width), w):
                v = sum(1 << i for i in pos)
                outside = any(gf2.dot(v, r) for r in dual) if replay else v not in zb
                if outside and not any(parity(v, r) for r in x):
                    physical = [
                        i * nb + spec["column"] if i < na else base["n"] + i - na
                        for i in pos
                    ]
                    lifted = sum(1 << i for i in physical)
                    if lifted in gf2.Echelon(merged[2]) or any(
                        parity(lifted, r) for r in merged[1]
                    ):
                        raise ValueError("failed to lift reduced witness")
                    witness = dict(axis="Z", qubits=physical)
                    break
            if witness:
                break
    if witness:
        return dict(
            result,
            status="counterexample",
            witness=witness,
            upper=len(witness["qubits"]),
        )
    return dict(result, status="proved", lower=threshold)
