"""Exact zero-float cyclotomic, dyadic, and outward-rounded rational arithmetic.

Implements:
- ZSqrt2: exact ring Z[sqrt(2)] with integer-squaring sign/order comparisons.
- ZOmega8: exact cyclotomic ring Z[omega_8] (omega = exp(i*pi/4), omega^4 = -1),
  conjugation, norm |z|^2 in Z[sqrt(2)], and real-subring projection.
- QSqrt2: exact field Q(sqrt(2)) with exact division, order comparisons, and
  certifiably outward-rounded Decimal upper/lower enclosures via integer isqrt.
- PolyQ: exact univariate polynomials over Q (fractions.Fraction) for weight
  enumerators, acceptance probabilities, and conditional error polynomials.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
from math import comb, isqrt
from typing import Iterable, Sequence

DEFAULT_SIG = 50


def parse_exact_fraction(value: int | str | Fraction) -> Fraction:
    """Parse an integer, rational string ('a/b' or decimal '1e-3'), or Fraction without float."""
    if isinstance(value, Fraction):
        return value
    if type(value) is int:
        return Fraction(value)
    if isinstance(value, str):
        s = value.strip()
        if not s:
            raise ValueError("Empty rational string")
        if "/" in s:
            num_s, den_s = s.split("/", 1)
            return Fraction(int(num_s.strip()), int(den_s.strip()))
        return Fraction(Decimal(s))
    raise TypeError(f"Float or unsupported numeric type rejected by zero-float gate: {type(value)!r}")


def _fraction_to_decimal_exact_bound(x: Fraction, sig: int, rounding: str) -> Decimal:
    if sig < 1:
        raise ValueError("sig must be >= 1")
    if x == 0:
        return Decimal(0)
    ctx = Context(prec=sig, rounding=rounding)
    # Divide with extra guard precision first, then round once at target sig in the same direction.
    guard = Context(prec=sig + 24, rounding=rounding)
    approx = guard.divide(Decimal(x.numerator), Decimal(x.denominator))
    return ctx.plus(approx)


def round_up_decimal(x: Fraction | int, sig: int = DEFAULT_SIG) -> Decimal:
    """Smallest Decimal with `sig` significant digits that is >= x."""
    return _fraction_to_decimal_exact_bound(Fraction(x), sig, ROUND_CEILING)


def round_down_decimal(x: Fraction | int, sig: int = DEFAULT_SIG) -> Decimal:
    """Largest Decimal with `sig` significant digits that is <= x."""
    return _fraction_to_decimal_exact_bound(Fraction(x), sig, ROUND_FLOOR)


def round_up_str(x: Fraction | int, sig: int = DEFAULT_SIG) -> str:
    return format(round_up_decimal(x, sig), "f" if Fraction(x) == 0 else "g")


def round_down_str(x: Fraction | int, sig: int = DEFAULT_SIG) -> str:
    return format(round_down_decimal(x, sig), "f" if Fraction(x) == 0 else "g")


def _sqrt2_rational_bounds(digits: int) -> tuple[Fraction, Fraction]:
    """Return exact (lo, hi) in Q with lo <= sqrt(2) <= hi and hi - lo == 10**(-digits)."""
    scale = 10**digits
    target = 2 * scale * scale
    lo_num = isqrt(target)
    if lo_num * lo_num == target:
        return Fraction(lo_num, scale), Fraction(lo_num, scale)
    return Fraction(lo_num, scale), Fraction(lo_num + 1, scale)


@dataclass(frozen=True)
class ZSqrt2:
    """Exact element a + b*sqrt(2) in Z[sqrt(2)]."""

    a: int
    b: int = 0

    def __post_init__(self) -> None:
        if type(self.a) is not int or type(self.b) is not int:
            raise TypeError("ZSqrt2 coefficients must be exact Python ints")

    def __add__(self, other: ZSqrt2 | int) -> ZSqrt2:
        if type(other) is int:
            return ZSqrt2(self.a + other, self.b)
        if isinstance(other, ZSqrt2):
            return ZSqrt2(self.a + other.a, self.b + other.b)
        return NotImplemented

    def __radd__(self, other: int) -> ZSqrt2:
        return self.__add__(other)

    def __sub__(self, other: ZSqrt2 | int) -> ZSqrt2:
        if type(other) is int:
            return ZSqrt2(self.a - other, self.b)
        if isinstance(other, ZSqrt2):
            return ZSqrt2(self.a - other.a, self.b - other.b)
        return NotImplemented

    def __neg__(self) -> ZSqrt2:
        return ZSqrt2(-self.a, -self.b)

    def __mul__(self, other: ZSqrt2 | int) -> ZSqrt2:
        if type(other) is int:
            return ZSqrt2(self.a * other, self.b * other)
        if isinstance(other, ZSqrt2):
            return ZSqrt2(
                self.a * other.a + 2 * self.b * other.b,
                self.a * other.b + self.b * other.a,
            )
        return NotImplemented

    def __rmul__(self, other: int) -> ZSqrt2:
        return self.__mul__(other)

    def galois_conj(self) -> ZSqrt2:
        """Send sqrt(2) -> -sqrt(2)."""
        return ZSqrt2(self.a, -self.b)

    def norm_z(self) -> int:
        """Field norm a^2 - 2*b^2 in Z."""
        return self.a * self.a - 2 * self.b * self.b

    def sign(self) -> int:
        """Exact sign in {-1, 0, +1} via integer squaring (zero floats)."""
        a, b = self.a, self.b
        if a == 0 and b == 0:
            return 0
        if a >= 0 and b >= 0:
            return 1
        if a <= 0 and b <= 0:
            return -1
        if a > 0 and b < 0:
            return 1 if a * a > 2 * b * b else -1
        # a < 0 and b > 0
        return 1 if 2 * b * b > a * a else -1

    def is_zero(self) -> bool:
        return self.a == 0 and self.b == 0

    def is_positive(self) -> bool:
        return self.sign() > 0

    def is_nonnegative(self) -> bool:
        return self.sign() >= 0

    def __lt__(self, other: ZSqrt2 | int) -> bool:
        return (self - other).sign() < 0

    def __le__(self, other: ZSqrt2 | int) -> bool:
        return (self - other).sign() <= 0

    def __gt__(self, other: ZSqrt2 | int) -> bool:
        return (self - other).sign() > 0

    def __ge__(self, other: ZSqrt2 | int) -> bool:
        return (self - other).sign() >= 0

    def to_qsqrt2(self) -> QSqrt2:
        return QSqrt2(Fraction(self.a), Fraction(self.b))

    def to_dict(self) -> dict[str, int]:
        return {"a": self.a, "b": self.b}


@dataclass(frozen=True)
class QSqrt2:
    """Exact element a + b*sqrt(2) in Q(sqrt(2))."""

    a: Fraction
    b: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", parse_exact_fraction(self.a))
        object.__setattr__(self, "b", parse_exact_fraction(self.b))

    @classmethod
    def from_int(cls, n: int) -> QSqrt2:
        return cls(Fraction(n), Fraction(0))

    def __add__(self, other: QSqrt2 | Fraction | int) -> QSqrt2:
        if isinstance(other, QSqrt2):
            return QSqrt2(self.a + other.a, self.b + other.b)
        f = parse_exact_fraction(other)
        return QSqrt2(self.a + f, self.b)

    def __radd__(self, other: Fraction | int) -> QSqrt2:
        return self.__add__(other)

    def __sub__(self, other: QSqrt2 | Fraction | int) -> QSqrt2:
        if isinstance(other, QSqrt2):
            return QSqrt2(self.a - other.a, self.b - other.b)
        f = parse_exact_fraction(other)
        return QSqrt2(self.a - f, self.b)

    def __rsub__(self, other: Fraction | int) -> QSqrt2:
        return (-self).__add__(other)

    def __neg__(self) -> QSqrt2:
        return QSqrt2(-self.a, -self.b)

    def __mul__(self, other: QSqrt2 | Fraction | int) -> QSqrt2:
        if isinstance(other, QSqrt2):
            return QSqrt2(
                self.a * other.a + 2 * self.b * other.b,
                self.a * other.b + self.b * other.a,
            )
        f = parse_exact_fraction(other)
        return QSqrt2(self.a * f, self.b * f)

    def __rmul__(self, other: Fraction | int) -> QSqrt2:
        return self.__mul__(other)

    def inv(self) -> QSqrt2:
        if self.a == 0 and self.b == 0:
            raise ZeroDivisionError("Division by zero in Q(sqrt(2))")
        denom = self.a * self.a - 2 * self.b * self.b
        return QSqrt2(self.a / denom, -self.b / denom)

    def __truediv__(self, other: QSqrt2 | Fraction | int) -> QSqrt2:
        if isinstance(other, QSqrt2):
            return self * other.inv()
        f = parse_exact_fraction(other)
        if f == 0:
            raise ZeroDivisionError("Division by zero in Q(sqrt(2))")
        return QSqrt2(self.a / f, self.b / f)

    def sign(self) -> int:
        if self.a == 0 and self.b == 0:
            return 0
        common_den = self.a.denominator * self.b.denominator
        ia = self.a.numerator * self.b.denominator
        ib = self.b.numerator * self.a.denominator
        assert common_den > 0
        return ZSqrt2(ia, ib).sign()

    def __lt__(self, other: QSqrt2 | Fraction | int) -> bool:
        return (self - other).sign() < 0

    def __le__(self, other: QSqrt2 | Fraction | int) -> bool:
        return (self - other).sign() <= 0

    def __gt__(self, other: QSqrt2 | Fraction | int) -> bool:
        return (self - other).sign() > 0

    def __ge__(self, other: QSqrt2 | Fraction | int) -> bool:
        return (self - other).sign() >= 0

    def enclosure_fractions(self, digits: int = DEFAULT_SIG + 12) -> tuple[Fraction, Fraction]:
        """Exact rational [lo, hi] enclosing a + b*sqrt(2)."""
        if self.b == 0:
            return self.a, self.a
        s_lo, s_hi = _sqrt2_rational_bounds(digits)
        if self.b > 0:
            return self.a + self.b * s_lo, self.a + self.b * s_hi
        return self.a + self.b * s_hi, self.a + self.b * s_lo

    def upper_decimal(self, sig: int = DEFAULT_SIG) -> Decimal:
        _, hi = self.enclosure_fractions(sig + 16)
        return round_up_decimal(hi, sig)

    def lower_decimal(self, sig: int = DEFAULT_SIG) -> Decimal:
        lo, _ = self.enclosure_fractions(sig + 16)
        return round_down_decimal(lo, sig)

    def to_dict(self) -> dict[str, str]:
        return {"a": str(self.a), "b": str(self.b)}


@dataclass(frozen=True)
class ZOmega8:
    """Exact element c0 + c1*w + c2*w^2 + c3*w^3 in Z[omega_8], w = exp(i*pi/4), w^4 = -1."""

    c0: int = 0
    c1: int = 0
    c2: int = 0
    c3: int = 0

    def __post_init__(self) -> None:
        for c in (self.c0, self.c1, self.c2, self.c3):
            if type(c) is not int:
                raise TypeError("ZOmega8 coefficients must be exact Python ints")

    @classmethod
    def zero(cls) -> ZOmega8:
        return cls(0, 0, 0, 0)

    @classmethod
    def one(cls) -> ZOmega8:
        return cls(1, 0, 0, 0)

    @classmethod
    def omega_pow(cls, k: int) -> ZOmega8:
        r = k % 8
        sgn = 1 if r < 4 else -1
        idx = r % 4
        arr = [0, 0, 0, 0]
        arr[idx] = sgn
        return cls(arr[0], arr[1], arr[2], arr[3])

    @classmethod
    def phase_sum(cls, exponents: Iterable[int]) -> ZOmega8:
        acc = [0, 0, 0, 0]
        for e in exponents:
            r = e % 8
            if r < 4:
                acc[r] += 1
            else:
                acc[r - 4] -= 1
        return cls(acc[0], acc[1], acc[2], acc[3])

    def is_zero(self) -> bool:
        return (self.c0 | self.c1 | self.c2 | self.c3) == 0

    def __add__(self, other: ZOmega8 | int) -> ZOmega8:
        if type(other) is int:
            return ZOmega8(self.c0 + other, self.c1, self.c2, self.c3)
        if isinstance(other, ZOmega8):
            return ZOmega8(
                self.c0 + other.c0,
                self.c1 + other.c1,
                self.c2 + other.c2,
                self.c3 + other.c3,
            )
        return NotImplemented

    def __radd__(self, other: int) -> ZOmega8:
        return self.__add__(other)

    def __sub__(self, other: ZOmega8 | int) -> ZOmega8:
        if type(other) is int:
            return ZOmega8(self.c0 - other, self.c1, self.c2, self.c3)
        if isinstance(other, ZOmega8):
            return ZOmega8(
                self.c0 - other.c0,
                self.c1 - other.c1,
                self.c2 - other.c2,
                self.c3 - other.c3,
            )
        return NotImplemented

    def __neg__(self) -> ZOmega8:
        return ZOmega8(-self.c0, -self.c1, -self.c2, -self.c3)

    def __mul__(self, other: ZOmega8 | int) -> ZOmega8:
        if type(other) is int:
            return ZOmega8(self.c0 * other, self.c1 * other, self.c2 * other, self.c3 * other)
        if isinstance(other, ZOmega8):
            a0, a1, a2, a3 = self.c0, self.c1, self.c2, self.c3
            b0, b1, b2, b3 = other.c0, other.c1, other.c2, other.c3
            return ZOmega8(
                a0 * b0 - a1 * b3 - a2 * b2 - a3 * b1,
                a0 * b1 + a1 * b0 - a2 * b3 - a3 * b2,
                a0 * b2 + a1 * b1 + a2 * b0 - a3 * b3,
                a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0,
            )
        return NotImplemented

    def __rmul__(self, other: int) -> ZOmega8:
        return self.__mul__(other)

    def conj(self) -> ZOmega8:
        """Complex conjugate: w* = -w^3, (w^2)* = -w^2, (w^3)* = -w."""
        return ZOmega8(self.c0, -self.c3, -self.c2, -self.c1)

    def is_real(self) -> bool:
        return self.c2 == 0 and self.c1 + self.c3 == 0

    def to_zsqrt2(self) -> ZSqrt2:
        """Project a real Z[omega_8] element (c2=0, c3=-c1) to c0 + c1*sqrt(2) since w - w^3 = sqrt(2)."""
        if not self.is_real():
            raise ValueError(f"ZOmega8 element {self} is not real")
        return ZSqrt2(self.c0, self.c1)

    def abs_sq(self) -> ZSqrt2:
        """Exact squared modulus z * conj(z) in Z[sqrt(2)]."""
        return (self * self.conj()).to_zsqrt2()

    def coeffs(self) -> tuple[int, int, int, int]:
        return (self.c0, self.c1, self.c2, self.c3)


@dataclass(frozen=True)
class PolyQ:
    """Exact univariate polynomial sum_i coeffs[i] * p^i over Fraction."""

    coeffs: tuple[Fraction, ...]

    def __post_init__(self) -> None:
        cleaned = [parse_exact_fraction(c) for c in self.coeffs]
        while len(cleaned) > 1 and cleaned[-1] == 0:
            cleaned.pop()
        if not cleaned:
            cleaned = [Fraction(0)]
        object.__setattr__(self, "coeffs", tuple(cleaned))

    @classmethod
    def from_weight_enumerator(cls, n: int, weights: Sequence[int]) -> PolyQ:
        """Convert weight counts W_w (0 <= w <= m) on n independent sites into exact polynomial in p:
        sum_{w=0}^m W_w * p^w * (1 - p)^(n - w).
        """
        out = [Fraction(0)] * (n + 1)
        for w, count in enumerate(weights):
            if count == 0:
                continue
            rem = n - w
            for j in range(rem + 1):
                deg = w + j
                term = count * comb(rem, j) * (-1 if (j & 1) else 1)
                out[deg] += Fraction(term)
        return cls(tuple(out))

    def eval_exact(self, p: Fraction | int | str) -> Fraction:
        pf = parse_exact_fraction(p)
        acc = Fraction(0)
        for c in reversed(self.coeffs):
            acc = acc * pf + c
        return acc

    def leading_nonzero(self, min_degree: int = 0) -> tuple[int, Fraction] | None:
        for deg in range(min_degree, len(self.coeffs)):
            if self.coeffs[deg] != 0:
                return deg, self.coeffs[deg]
        return None


def eval_binomial_tail_exact(n: int, t: int, q: Fraction) -> Fraction:
    """Exact probability Pr[Bin(n, q) > t] in Fraction."""
    if not (0 <= q <= 1):
        raise ValueError(f"Probability q={q} outside [0, 1]")
    if t < 0:
        return Fraction(1)
    if t >= n or q == 0:
        return Fraction(0)
    if q == 1:
        return Fraction(1)
    one_minus_q = Fraction(1) - q
    # Sum smaller side for speed and exactness.
    if n - (t + 1) <= t + 1:
        acc = Fraction(0)
        for k in range(t + 1, n + 1):
            acc += comb(n, k) * (q**k) * (one_minus_q ** (n - k))
        return acc
    acc = Fraction(0)
    for k in range(0, t + 1):
        acc += comb(n, k) * (q**k) * (one_minus_q ** (n - k))
    return Fraction(1) - acc
