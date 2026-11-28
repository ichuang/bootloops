"""One-location responses of a parallel verified-cat Pauli check.

The cat is supplied in an ideal GHZ state. Each leg controls P on one
data qubit and is then read in X. Each tuple gives a data Pauli after the
coupling and the flip of this shot's parity relative to the incoming data.
"""

from itertools import product


def cat_faults(width, check, axis):
    for q in range(width):
        bit = 1 << q
        if not check & bit:
            for x, z in ((1, 0), (0, 1), (1, 1)):
                yield ("idle", q, x, z), x * bit, z * bit, 0
            continue
        for mx, mz in ((1, 0), (0, 1), (1, 1)):
            yield (
                ("cat-leg", q, mx, mz),
                bit * mx * (axis == "X"),
                bit * mx * (axis == "Z"),
                mz,
            )
        for mx, mz, x, z in product(range(2), repeat=4):
            if mx or mz or x or z:
                yield ("controlled-" + axis, q, mx, mz, x, z), bit * x, bit * z, mz
        yield ("read-X", q), 0, 0, 1


def conjugate_ccz(x, linear, edges, sites):
    """Diagonal difference f(v+x)+f(v) over GF(2), up to a global sign.

    Residual is X^x D(v) with D a signed quadratic phase. Keep every CZ.
    """
    a, b, c = sites
    for q, r, s in ((a, b, c), (b, a, c), (c, a, b)):
        if x >> q & 1:
            edge = tuple(sorted((r, s)))
            edges.symmetric_difference_update([edge])
        if x >> r & 1 and x >> s & 1:
            linear ^= 1 << q
    return linear


def z_expansion(linear, edges):
    """Signed integer Pauli coefficients, common denominator 2^len(edges)."""
    terms = {linear: 1}
    for a, b in sorted(edges):
        merged = {}
        for word, weight in terms.items():
            for mask, sign in (
                (0, 1),
                (1 << a, 1),
                (1 << b, 1),
                ((1 << a) | (1 << b), -1),
            ):
                v = word ^ mask
                merged[v] = merged.get(v, 0) + sign * weight
        terms = {v: k for v, k in merged.items() if k}
    return terms
