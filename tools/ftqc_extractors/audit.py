"""Cross-check the cyclic certificate on the actual Definition 16 CSS matrices."""

from ftqc_codes import gf2
from .extractor import column


def all_sums(rows):
    values = [0]
    for r in rows:
        values.extend([v ^ r for v in values])
    return values


def audit_canonical(s, result):
    n, p, edges = s["n"], s["polynomial"], s["edges"]
    strings = [
        "".join(str((p >> ((j - i) % n)) & 1) for j in range(n)) for i in range(n)
    ]
    checks = [int(row[::-1], 2) for row in strings]
    words = all_sums(gf2.nullspace(checks, n))
    k = (len(words) - 1).bit_length()
    d = min(w.bit_count() for w in words if w)
    if n > 20 or 2**n * len(words) ** 2 > 200000000:
        raise ValueError("replay exceeds canonical budget")
    mincut = min(
        sum(((u >> i) ^ (u >> j)) & 1 for i, j in edges) for u in range(1, 1 << (n - 1))
    )
    for c in words[1:]:
        spec = dict(
            a=strings,
            b=strings,
            vertices=n,
            ports=list(range(n)),
            edges=edges,
            column=0,
            word=c,
        )
        _, (width, x, z) = column(spec)
        # Distinct surviving pure-data cosets must exhaust the Z logical quotient.
        zb = gf2.Echelon(z)
        if any(gf2.dot(w, r) for w in words for r in x):
            raise ValueError("invalid reduced logical representative")
        count = width - gf2.rank(x) - zb.rank
        reps = {zb.reduce(w): w for w in words}
        if len(reps) != 2**count or c not in zb:
            raise ValueError("incomplete reduced quotient")
        reps.pop(0)
        # Include every Z stabilizer representative. Each vertex product acts
        # on edge qubits as a cut, bounding all members of that coset.
        for op in all_sums([z[i] for i in gf2.independent_rows(z)]):
            if (op >> n).bit_count() >= d:
                continue
            if any((op ^ w).bit_count() < d for w in reps.values()):
                raise ValueError("reduced subcode violates bound")
    if result["lower"] != min(d, k * k, k * mincut) or result["connectivity"] != mincut:
        raise ValueError("incorrect full-extractor bound")


def audit_hgp(spec, result):
    from .algebra import hgp, seed, transpose

    a, na = seed(spec["a"])
    b, nb = seed(spec["b"])
    n, hx, hz = hgp(a, na, b, nb)
    logicals = []
    for checks, stabilizers in [(hz, hx), (hx, hz)]:
        space = gf2.Echelon(stabilizers)
        logicals.append([w for w in gf2.nullspace(checks, n) if space.add(w)])
    xs, zs = logicals
    if len(xs) != result["k"] or len(xs) != len(zs) or n != result["n"]:
        raise ValueError("logical dimension mismatch")
    for i, op in enumerate(xs):
        pivot = next((j for j in range(i, len(zs)) if gf2.dot(op, zs[j])), None)
        if pivot is None:
            raise ValueError("singular logical pairing")
        zs[i], zs[pivot] = zs[pivot], zs[i]
        for j in range(len(zs)):
            if j != i and gf2.dot(op, zs[j]):
                zs[j] ^= zs[i]
    if any(
        gf2.dot(x, z) != (i == j) for i, x in enumerate(xs) for j, z in enumerate(zs)
    ):
        raise ValueError("logical symplectic basis failed")
    minima = []
    for h1, n1, h2, n2 in [
        (a, na, b, nb),
        (transpose(a, na), len(a), transpose(b, nb), len(b)),
    ]:
        bases = [gf2.nullspace(h1, n1), gf2.nullspace(h2, n2)]
        if all(bases):
            minima += [min(v.bit_count() for v in all_sums(K) if v) for K in bases]
    if min(minima, default=None) != result["d"]:
        raise ValueError("seed distance mismatch")
