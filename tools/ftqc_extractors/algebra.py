"""Binary row spaces, CSS codes and hypergraph products. Bit i is qubit i."""

from itertools import combinations
from math import comb


class Limit(Exception):
    """A calculation exceeded its explicit work budget."""


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def parity(a, b):
    return (a & b).bit_count() % 2


def reduce_row(x, basis):
    for pivot in sorted(basis, reverse=True):
        if x >> pivot & 1:
            x ^= basis[pivot]
    return x


def echelon(rows):
    basis = {}
    for x in rows:
        x = reduce_row(x, basis)
        if x:
            basis[x.bit_length() - 1] = x
    return basis


def kernel(rows, n):
    basis = echelon(rows)
    result = []
    for j in range(n):
        if j in basis:
            continue
        x = 1 << j
        for p in sorted(basis):
            x |= parity(x, basis[p]) << p
        result.append(x)
    return result


def span(basis, limit=1 << 20):
    if 2 ** len(basis) > limit:
        raise Limit("kernel enumeration budget")
    words = [0]
    for row in basis:
        words += [v ^ row for v in words]
    return words


def transpose(rows, width):
    cols = [0] * width
    for i, row in enumerate(rows):
        for j in bits(row):
            cols[j] |= 1 << i
    return cols


def seed(rows):
    if (
        not isinstance(rows, list)
        or not rows
        or not all(isinstance(r, str) and r and set(r) <= set("01") for r in rows)
    ):
        raise ValueError("seed requires nonempty binary row strings")
    n = len(rows[0])
    if any(len(r) != n for r in rows):
        raise ValueError("ragged seed")
    return [int(r[::-1], 2) for r in rows], n


def hgp(a, na, b, nb):
    ma, mb = len(a), len(b)
    offset = na * nb
    at, bt = transpose(a, na), transpose(b, nb)
    x = [
        (sum(1 << (u * nb + j) for u in bits(a[i])) | (bt[j] << (offset + i * mb)))
        for i in range(ma)
        for j in range(nb)
    ]
    z = [
        (b[j] << (i * nb)) | sum(1 << (offset + u * mb + j) for u in bits(at[i]))
        for i in range(na)
        for j in range(mb)
    ]
    return offset + ma * mb, x, z


def css_dimension(code):
    n, x, z = code
    if any(parity(a, b) for a in x for b in z):
        raise ValueError("anticommuting checks")
    return n - len(echelon(x)) - len(echelon(z))


def hgp_parameters(a, na, b, nb):
    distances = []
    for left, nl, right, nr in [
        (a, na, b, nb),
        (transpose(a, na), len(a), transpose(b, nb), len(b)),
    ]:
        ka, kb = kernel(left, nl), kernel(right, nr)
        if ka and kb:
            for basis in (ka, kb):
                distances.append(min(v.bit_count() for v in span(basis) if v))
    code = hgp(a, na, b, nb)
    return dict(n=code[0], k=css_dimension(code), d=min(distances, default=None))


def distance_below(code, threshold, budget=2000000):
    n, x, z = code
    xb, zb = echelon(x), echelon(z)
    for w in range(1, min(threshold, n + 1)):
        budget -= comb(n, w)
        if budget < 0:
            raise Limit("direct CSS distance enumeration")
        for positions in combinations(range(n), w):
            op = sum(1 << i for i in positions)
            for axis, checks, stabilizers in [("X", z, xb), ("Z", x, zb)]:
                if not any(parity(op, c) for c in checks) and reduce_row(
                    op, stabilizers
                ):
                    return dict(axis=axis, qubits=list(positions))
    return None


def cyclic(n, polynomial):
    if not 1 < n <= 128 or not 0 < polynomial < 2**n or polynomial % 2 == 0:
        raise ValueError("invalid cyclic check polynomial")
    mask = 2**n - 1
    return [((polynomial << i) | (polynomial >> (n - i))) & mask for i in range(n)]
