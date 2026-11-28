"""Directed square root for nonnegative Decimal bounds."""

from decimal import Context, Decimal
from fractions import Fraction


def sqrt_upper(value: Decimal, context: Context) -> Decimal:
    # Decimal.sqrt always rounds to nearest. Check the squared result exactly.
    root = context.sqrt(value)
    return context.next_plus(root) if Fraction(root) ** 2 < Fraction(value) else root
