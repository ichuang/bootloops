"""Exact GF(2), GF(2^m), and Algebraic Normal Form (ANF) kernel.

100% Python standard library, MIT-licensed, zero float literals.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


# Standard minimal-weight / Conway-compatible primitive polynomials over GF(2)
# for degrees 1..16 (bit i = coefficient of x^i; primitive element = x = 2 for m >= 2).
_DEFAULT_PRIMITIVE_POLYS: dict[int, int] = {
    1: 0b11,                 # x + 1
    2: 0b111,                # x^2 + x + 1
    3: 0b1011,               # x^3 + x + 1
    4: 0b10011,              # x^4 + x + 1
    5: 0b100101,             # x^5 + x^2 + 1
    6: 0b1000011,            # x^6 + x + 1
    7: 0b10000011,           # x^7 + x + 1
    8: 0b100011101,          # x^8 + x^4 + x^3 + x^2 + 1
    9: 0b1000010001,         # x^9 + x^4 + 1
    10: 0b10000001001,       # x^10 + x^3 + 1
    11: 0b100000000101,      # x^11 + x^2 + 1
    12: 0b1000001010011,     # x^12 + x^6 + x^4 + x + 1
    13: 0b10000000011011,    # x^13 + x^4 + x^3 + x + 1
    14: 0b100000000101011,   # x^14 + x^5 + x^3 + x + 1
    15: 0b1000000000000011,  # x^15 + x + 1
    16: 0b10001000000001011, # x^16 + x^12 + x^3 + x + 1
}


def poly_deg(p: int) -> int:
    """Degree of a polynomial over GF(2) encoded as a non-negative integer (-1 for 0)."""
    if p < 0:
        raise ValueError("GF(2) polynomial bitmask must be non-negative")
    return p.bit_length() - 1


def poly_mul_gf2(a: int, b: int) -> int:
    """Carryless product of two polynomials over GF(2)."""
    res = 0
    x, y = int(a), int(b)
    while y:
        if y & 1:
            res ^= x
        x <<= 1
        y >>= 1
    return res


def poly_divmod_gf2(a: int, b: int) -> tuple[int, int]:
    """Euclidean division (q, r) of GF(2) polynomials with a = q * b ^ r."""
    if b <= 0:
        raise ValueError("divisor polynomial must be nonzero")
    deg_b = poly_deg(b)
    q = 0
    r = int(a)
    while r and poly_deg(r) >= deg_b:
        shift = poly_deg(r) - deg_b
        q ^= 1 << shift
        r ^= b << shift
    return q, r


def poly_mod_gf2(a: int, mod: int) -> int:
    """Remainder of GF(2) polynomial `a` modulo `mod`."""
    return poly_divmod_gf2(a, mod)[1]


def poly_gcd_gf2(a: int, b: int) -> int:
    """Greatest common divisor of two GF(2) polynomials."""
    x, y = int(a), int(b)
    while y:
        x, y = y, poly_mod_gf2(x, y)
    return x


def poly_powmod_gf2(base: int, exp: int, mod: int) -> int:
    """Compute `base ** exp (mod mod)` in GF(2)[x]."""
    if exp < 0:
        raise ValueError("exponent must be non-negative")
    res = 1
    cur = poly_mod_gf2(base, mod)
    e = int(exp)
    while e:
        if e & 1:
            res = poly_mod_gf2(poly_mul_gf2(res, cur), mod)
        cur = poly_mod_gf2(poly_mul_gf2(cur, cur), mod)
        e >>= 1
    return res


def prime_factors(n: int) -> list[int]:
    """Distinct prime factors of positive integer `n` in ascending order."""
    if n < 1:
        raise ValueError("n must be >= 1")
    out: list[int] = []
    m = int(n)
    d = 2
    while d * d <= m:
        if m % d == 0:
            out.append(d)
            while m % d == 0:
                m //= d
        d = 3 if d == 2 else d + 2
    if m > 1:
        out.append(m)
    return out


def is_irreducible_gf2(poly: int) -> bool:
    """Rabin irreducibility test for a polynomial over GF(2)."""
    m = poly_deg(poly)
    if m < 1:
        return False
    if m == 1:
        return True
    if (poly & 1) == 0:
        return False
    x_poly = poly_mod_gf2(0b10, poly)
    for q in prime_factors(m):
        # Compute x^(2^(m/q)) mod poly
        h = x_poly
        for _ in range(m // q):
            h = poly_mod_gf2(poly_mul_gf2(h, h), poly)
        if poly_gcd_gf2(poly, h ^ x_poly) != 1:
            return False
    h = x_poly
    for _ in range(m):
        h = poly_mod_gf2(poly_mul_gf2(h, h), poly)
    return (h ^ x_poly) == 0


def _find_primitive_poly(degree: int) -> int:
    if degree in _DEFAULT_PRIMITIVE_POLYS:
        return _DEFAULT_PRIMITIVE_POLYS[degree]
    if not 1 <= degree <= 32:
        raise ValueError(f"GF2Field degree must lie in [1, 32], got {degree}")
    order = (1 << degree) - 1
    pf = prime_factors(order)
    base = (1 << degree) | 1
    for mid in range(0, 1 << (degree - 1)):
        cand = base | (mid << 1)
        if not is_irreducible_gf2(cand):
            continue
        if all(poly_powmod_gf2(0b10, order // q, cand) != 1 for q in pf):
            return cand
    raise RuntimeError(f"failed to locate primitive polynomial of degree {degree}")


@dataclass(frozen=True)
class GF2Field:
    """Finite field GF(2^m) = GF(2)[x] / (modulus)."""

    degree: int
    modulus: int
    primitive_element: int

    def __init__(
        self,
        degree: int,
        modulus: int | None = None,
        primitive_element: int | None = None,
    ) -> None:
        if not 1 <= int(degree) <= 32:
            raise ValueError(f"degree must be in 1..32, got {degree}")
        deg = int(degree)
        mod = _find_primitive_poly(deg) if modulus is None else int(modulus)
        if poly_deg(mod) != deg or not is_irreducible_gf2(mod):
            raise ValueError(f"modulus {bin(mod)} is not an irreducible degree-{deg} polynomial")
        order = (1 << deg) - 1
        pf = prime_factors(order)
        if primitive_element is None:
            prim = 1 if deg == 1 else 0b10
            if deg > 1 and not all(poly_powmod_gf2(prim, order // q, mod) != 1 for q in pf):
                for c in range(2, 1 << deg):
                    if all(poly_powmod_gf2(c, order // q, mod) != 1 for q in pf):
                        prim = c
                        break
        else:
            prim = int(primitive_element)
            if not (1 <= prim < (1 << deg)):
                raise ValueError("primitive_element out of field range")
            if deg > 1 and not all(poly_powmod_gf2(prim, order // q, mod) != 1 for q in pf):
                raise ValueError(f"{prim} is not primitive in GF(2^{deg})")
        object.__setattr__(self, "degree", deg)
        object.__setattr__(self, "modulus", mod)
        object.__setattr__(self, "primitive_element", prim)

    @property
    def order(self) -> int:
        return 1 << self.degree

    def is_element(self, a: int) -> bool:
        return 0 <= int(a) < (1 << self.degree)

    def add(self, a: int, b: int) -> int:
        return int(a) ^ int(b)

    def mul(self, a: int, b: int) -> int:
        if not (self.is_element(a) and self.is_element(b)):
            raise ValueError("field operand out of range")
        return poly_mod_gf2(poly_mul_gf2(a, b), self.modulus)

    def square(self, a: int) -> int:
        return self.mul(a, a)

    def pow(self, a: int, exponent: int) -> int:
        if not self.is_element(a):
            raise ValueError("field operand out of range")
        if exponent < 0:
            return self.pow(self.invert(a), -exponent)
        return poly_powmod_gf2(a, exponent, self.modulus)

    def invert(self, a: int) -> int:
        if not self.is_element(a) or a == 0:
            raise ZeroDivisionError("cannot invert 0 in GF(2^m)")
        if self.degree == 1:
            return 1
        return poly_powmod_gf2(a, (1 << self.degree) - 2, self.modulus)

    def div(self, a: int, b: int) -> int:
        return self.mul(a, self.invert(b))

    def mul_const_matrix(self, factor: int) -> list[int]:
        """Return m rows (as m-bit ints): row k has bit i set iff (a * factor)_k depends on a_i."""
        if not self.is_element(factor):
            raise ValueError("factor out of range")
        cols = [self.mul(1 << i, factor) for i in range(self.degree)]
        rows: list[int] = []
        for k in range(self.degree):
            row = 0
            for i in range(self.degree):
                if (cols[i] >> k) & 1:
                    row |= 1 << i
            rows.append(row)
        return rows

    def mul_bilinear_tensor(self) -> tuple[tuple[int, ...], ...]:
        """Return `tensor[k][i]` as an m-bit mask in `j` such that:

        `(a * b)_k = parity_i( a_i * parity_j( tensor[k][i]_j * b_j ) )`.
        """
        m = self.degree
        basis_prod = [[self.mul(1 << i, 1 << j) for j in range(m)] for i in range(m)]
        tensor: list[tuple[int, ...]] = []
        for k in range(m):
            mat_k: list[int] = []
            for i in range(m):
                mask_j = 0
                for j in range(m):
                    if (basis_prod[i][j] >> k) & 1:
                        mask_j |= 1 << j
                mat_k.append(mask_j)
            tensor.append(tuple(mat_k))
        return tuple(tensor)


# ---------------------------------------------------------------------------
# Binary linear algebra over GF(2) (rows packed as Python ints)
# ---------------------------------------------------------------------------


def gf2_rref(rows: Iterable[int], ncols: int) -> tuple[list[int], list[int]]:
    """Reduced row echelon form over GF(2); returns `(rref_rows, pivot_cols)`."""
    mat = [int(r) & ((1 << ncols) - 1) for r in rows]
    pivots: list[int] = []
    r_idx = 0
    for col in range(ncols):
        bit = 1 << col
        pivot = None
        for i in range(r_idx, len(mat)):
            if mat[i] & bit:
                pivot = i
                break
        if pivot is None:
            continue
        mat[r_idx], mat[pivot] = mat[pivot], mat[r_idx]
        piv_row = mat[r_idx]
        for i in range(len(mat)):
            if i != r_idx and (mat[i] & bit):
                mat[i] ^= piv_row
        pivots.append(col)
        r_idx += 1
        if r_idx == len(mat):
            break
    return [r for r in mat if r != 0], pivots


def gf2_rank(rows: Iterable[int], ncols: int) -> int:
    return len(gf2_rref(rows, ncols)[1])


def gf2_invert_matrix(rows: list[int], n: int) -> list[int]:
    """Invert an n x n binary matrix given as `n` n-bit row masks (bit j = col j)."""
    if len(rows) != n:
        raise ValueError("matrix must be square")
    aug = [((rows[i] & ((1 << n) - 1)) | (1 << (n + i))) for i in range(n)]
    for col in range(n):
        bit = 1 << col
        piv = next((r for r in range(col, n) if aug[r] & bit), None)
        if piv is None:
            raise ValueError("GF(2) matrix is singular")
        aug[col], aug[piv] = aug[piv], aug[col]
        prow = aug[col]
        for r in range(n):
            if r != col and (aug[r] & bit):
                aug[r] ^= prow
    return [(aug[i] >> n) & ((1 << n) - 1) for i in range(n)]


def gf2_lup_decompose(rows: list[int], n: int) -> list[tuple[str, int, int]]:
    """Decompose an invertible n x n GF(2) linear map `y_i = parity(rows[i] & x)`
    into reversible elementary row operations `[("SWAP", i, j) | ("CX", ctrl, dst), ...]`
    such that executing the returned gates in order maps `|x> -> |A x>`.
    """
    inv_rows = gf2_invert_matrix(rows, n)
    # Row reduction on A to I: A = E_1^-1 ... E_k^-1 => executing E_k ... E_1 on |x> gives |A x>.
    work = [r & ((1 << n) - 1) for r in rows]
    _ = inv_rows  # certified invertible above
    ops_rev: list[tuple[str, int, int]] = []
    for col in range(n):
        bit = 1 << col
        piv = next(r for r in range(col, n) if work[r] & bit)
        if piv != col:
            work[col], work[piv] = work[piv], work[col]
            ops_rev.append(("SWAP", col, piv))
        for r in range(n):
            if r != col and (work[r] & bit):
                work[r] ^= work[col]
                ops_rev.append(("CX", col, r))
    ops_rev.reverse()
    return ops_rev


# ---------------------------------------------------------------------------
# Algebraic Normal Form (ANF) over GF(2)[v_0, ..., v_{k-1}]
# Represented as `frozenset[int]` of nonzero monomial variable-bitmasks.
# `0` (empty bitmask) represents the constant 1!
# ---------------------------------------------------------------------------

ANF = frozenset[int]

ANF_ZERO: ANF = frozenset()
ANF_ONE: ANF = frozenset({0})


def anf_var(var_index: int) -> ANF:
    if var_index < 0:
        raise ValueError("variable index must be >= 0")
    return frozenset({1 << var_index})


def anf_const(val: int | bool) -> ANF:
    return ANF_ONE if bool(val) else ANF_ZERO


def anf_xor(a: ANF, b: ANF) -> ANF:
    return a ^ b


def anf_not(a: ANF) -> ANF:
    return a ^ ANF_ONE


def anf_and(a: ANF, b: ANF, max_degree: int | None = None) -> ANF:
    """Multiply two multilinear GF(2) polynomials (`x_i^2 = x_i`)."""
    acc: set[int] = set()
    for ma in a:
        for mb in b:
            m = ma | mb
            if max_degree is not None and m.bit_count() > max_degree:
                raise ValueError(f"ANF degree {m.bit_count()} exceeds cap {max_degree}")
            if m in acc:
                acc.remove(m)
            else:
                acc.add(m)
    return frozenset(acc)


def anf_degree(a: ANF) -> int:
    if not a:
        return -1
    return max(m.bit_count() for m in a)


def anf_vars_mask(a: ANF) -> int:
    mask = 0
    for m in a:
        mask |= m
    return mask


def anf_eval(a: ANF, assignment_mask: int) -> int:
    acc = 0
    for m in a:
        if (assignment_mask & m) == m:
            acc ^= 1
    return acc
