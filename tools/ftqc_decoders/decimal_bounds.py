"""Zero-float exact rational and outward-rounded Decimal bounds for `ftqc_decoders`.

All probability, log-odds, and logical-error-rate interval calculations in
`ftqc_decoders` route through this module:
  * Exact probabilities live in `fractions.Fraction` (`0 <= p <= 1`), rejecting
    `float` and `bool` inputs (`RefusalError("FLOAT-INPUT-REFUSED", ...)`).
  * Integer log-odds weights `floor(scale * ln((1 - p) / p))` are evaluated via
    arbitrary-precision `decimal.Context.ln` with zero `float` literals or calls.
  * Every exported `Decimal` string pair `(lo_str, hi_str)` uses explicit
    `ROUND_FLOOR` / `ROUND_CEILING` contexts (`ctx.plus`, `.copy_negate()`,
    `.copy_abs()`) and checks runtime containment `Fraction(lo) <= exact <= Fraction(hi)`.
"""

from __future__ import annotations

from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction


class RefusalError(ValueError):
    """Structured refusal carrying a machine-checkable BootLoops error code."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"[{code}] {message}")
        self.code = code
        self.detail = message


def parse_exact_probability(
    value: str | int | Fraction,
    *,
    allow_zero: bool = True,
    allow_one: bool = True,
    max_half: bool = False,
) -> Fraction:
    """Parse an exact probability in `[0, 1]` without `float` conversion."""
    if isinstance(value, bool) or not isinstance(value, (str, int, Fraction)):
        raise RefusalError(
            "FLOAT-INPUT-REFUSED",
            f"Probability must be str, int, or Fraction (never float/bool); got {type(value).__name__}",
        )
    if isinstance(value, Fraction):
        frac = value
    elif isinstance(value, int):
        frac = Fraction(value, 1)
    else:
        text = value.strip()
        if not text:
            raise RefusalError("INVALID-PROBABILITY", "Empty probability string")
        try:
            if "/" in text:
                num_s, den_s = text.split("/", 1)
                frac = Fraction(int(num_s.strip()), int(den_s.strip()))
            else:
                frac = Fraction(text)
        except (ValueError, ZeroDivisionError) as exc:
            raise RefusalError("INVALID-PROBABILITY", f"Cannot parse probability {value!r}: {exc}") from exc

    if frac < 0 or frac > 1:
        raise RefusalError("PROBABILITY-OUT-OF-RANGE", f"Probability {frac} lies outside [0, 1]")
    if not allow_zero and frac == 0:
        raise RefusalError("ZERO-PROBABILITY-REFUSED", "Strictly positive probability required")
    if not allow_one and frac == 1:
        raise RefusalError("UNIT-PROBABILITY-REFUSED", "Probability strictly below 1 required")
    if max_half and frac > Fraction(1, 2):
        raise RefusalError("PROBABILITY-ABOVE-HALF", f"Error prior {frac} exceeds 1/2")
    return frac


def xor_prob(p1: Fraction, p2: Fraction) -> Fraction:
    """Exact Bernoulli parity-sum probability `p1 * (1 - p2) + p2 * (1 - p1)`."""
    return p1 + p2 - 2 * p1 * p2


def scaled_log_odds_int(p: Fraction, scale: int = 1024) -> int:
    """Deterministic integer log-odds `floor(scale * ln((1 - p) / p))` via `Decimal.ln`.

    Matches fixed-point FPGA/ASIC message-passing quantization (Chen et al., IEEE
    Trans. Commun. 53, 1288, 2005; Roffe et al., arXiv:2005.07016) with zero floats.
    """
    if type(scale) is not int or scale <= 0:
        raise RefusalError("INVALID-SCALE", f"Log-odds scale must be a positive int, got {scale!r}")
    if p <= 0:
        return scale * 64
    if p >= Fraction(1, 2):
        if p == Fraction(1, 2):
            return 0
        if p >= 1:
            return -(scale * 64)
    ctx = Context(prec=50, rounding=ROUND_FLOOR)
    one_minus_p = Fraction(1, 1) - p
    # Ratio `(1 - p) / p = (one_minus_p.numerator * p.denominator) / (one_minus_p.denominator * p.numerator)`
    num_d = Decimal(one_minus_p.numerator * p.denominator)
    den_d = Decimal(one_minus_p.denominator * p.numerator)
    ratio = ctx.divide(num_d, den_d)
    ln_val = ctx.ln(ratio)
    scaled = ctx.multiply(Decimal(scale), ln_val)
    w = int(scaled.to_integral_value(rounding=ROUND_FLOOR))
    if p < Fraction(1, 2) and w < 1:
        return 1
    return w


def fraction_to_outward_strs(val: Fraction, digits: int = 30) -> tuple[str, str]:
    """Convert exact nonnegative `Fraction` to `(lo_str, hi_str)` outward `Decimal` bounds."""
    if val < 0:
        raise RefusalError("NEGATIVE-BOUND", f"Expected nonnegative Fraction, got {val}")
    if type(digits) is not int or digits < 6:
        raise RefusalError("INVALID-PRECISION", f"Precision digits must be int >= 6, got {digits!r}")

    work_prec = max(100, digits + 40)
    ctx_work_lo = Context(prec=work_prec, rounding=ROUND_FLOOR)
    ctx_work_hi = Context(prec=work_prec, rounding=ROUND_CEILING)
    ctx_out_lo = Context(prec=digits, rounding=ROUND_FLOOR)
    ctx_out_hi = Context(prec=digits, rounding=ROUND_CEILING)

    num_d = Decimal(val.numerator)
    den_d = Decimal(val.denominator)
    lo = ctx_out_lo.plus(ctx_work_lo.divide(num_d, den_d))
    hi = ctx_out_hi.plus(ctx_work_hi.divide(num_d, den_d))

    lo_s = "0" if lo == 0 else str(lo)
    hi_s = "0" if hi == 0 else str(hi)
    if not (Fraction(lo_s) <= val <= Fraction(hi_s)):
        raise AssertionError(f"Outward interval containment violated for {val}: [{lo_s}, {hi_s}]")
    return lo_s, hi_s


def decimal_str_upper(val: Decimal, digits: int = 30) -> str:
    """Round a nonnegative `Decimal` strictly upward (`ROUND_CEILING`) at `digits` precision."""
    if not val.is_finite() or val < 0:
        raise RefusalError("INVALID-DECIMAL-UPPER", f"Expected finite nonnegative Decimal, got {val!r}")
    ctx_hi = Context(prec=digits, rounding=ROUND_CEILING)
    out = ctx_hi.plus(val)
    return "0" if out == 0 else str(out)
