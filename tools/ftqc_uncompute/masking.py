"""Superposition masking and streaming traveling-eigenstate ladder verifiers.

Implements formal validation and outward-rounded `Decimal` error bounds from:
  - Van Kirk & Gidney, "Streaming Regev's Quantum Factoring Algorithm"
    (arXiv:2610.06733v1, Sections 3-5 & Appendix B)
  - Gidney, "Approximate encoded permutations and piecewise quantum adders" (2019) [runways]
  - Gidney, "Halving the cost of quantum addition" (2018) [phase-gradient catalysis]

Refusal codes:
  - `MASK-NON-BIJECTIVE-GARBAGE`    — branch update is not a bijection on the mask group
  - `MASK-SUPPORT-OVERFLOW-HAZARD`  — average-case truncated subroutine overflows on full masked support
  - `LADDER-NON-UNIMODULAR-STEP`    — streaming conjugation step has |det(M)| != 1
  - `LADDER-RESIDUAL-ENTANGLEMENT`  — streamed control bit fails to decouple from work registers

100% Python standard library, MIT-licensed, zero float literals or float() calls.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
import math
from typing import Any, Sequence


VALID_MASK_KINDS = frozenset({
    "mult_mod_n",      # |M_N> = (1/sqrt(phi(N))) sum_{gcd(r,N)=1} |r>
    "uniform_mod_n",   # |R_N> = (1/sqrt(N)) sum_{r=0}^{N-1} |r>  (approxs |M_N> without knowing phi(N))
    "add_mod_2w",      # |+>^{otimes w}
    "phase_grad",      # |grad_w> = 2^{-w/2} sum_{y=0}^{2^w - 1} exp(-2*pi*i*y / 2^w) |y>
    "carry_runway",    # 2^{-w_pad/2} sum_{r=0}^{2^{w_pad}-1} |r> absorbing bounded carries <= delta_max
})


def euler_totient_from_factors(n: int, prime_factors_list: Sequence[int] | None = None) -> int:
    """Exact Euler totient `phi(n)` from `prime_factors_list` (or trial division up to `10^7`)."""
    if n < 2:
        raise ValueError("modulus N must be >= 2")
    if prime_factors_list is not None:
        phi = n
        for p in sorted(set(int(x) for x in prime_factors_list)):
            if p < 2 or n % p != 0:
                raise ValueError(f"{p} is not a prime factor of {n}")
            phi = (phi // p) * (p - 1)
        return phi

    rem = int(n)
    phi = int(n)
    d = 2
    limit = 10_000_000
    while d * d <= rem and d <= limit:
        if rem % d == 0:
            phi = (phi // d) * (d - 1)
            while rem % d == 0:
                rem //= d
        d = 3 if d == 2 else d + 2
    if d * d > rem:
        if rem > 1:
            phi = (phi // rem) * (rem - 1)
        return phi
    raise ValueError(
        f"Modulus N={n} has large unfactored composite part; supply `prime_factors_list` or `min_prime_factor` bound"
    )


def uniform_residue_mask_bound(
    modulus: int,
    *,
    runs: int = 1,
    prime_factors_list: Sequence[int] | None = None,
    min_prime_factor: int | None = None,
    decimal_digits: int = 24,
) -> dict[str, Any]:
    """Outward-rounded (`ROUND_CEILING`) trace / total-variation distance bound when
    approximating `|M_N>` by the public uniform residue superposition `|R_N>`
    (Van Kirk & Gidney 2026, Section 3.2):

        `D_tr(|R_N>, |M_N>) = sqrt(1 - phi(N) / N)`,
        `d_TV_total <= min(1, runs * sqrt(1 - phi(N) / N))`.

    Supports both exact `phi(N)` (for test moduli or known factorizations) and factor-free
    cryptographic bounds when all prime divisors of `N` satisfy `p >= min_prime_factor`.
    """
    if modulus < 2:
        raise ValueError("modulus must be >= 2")
    if runs < 1:
        raise ValueError("runs must be >= 1")

    ctx_floor = Context(prec=max(56, decimal_digits + 20), rounding=ROUND_FLOOR)
    ctx_ceil = Context(prec=max(56, decimal_digits + 20), rounding=ROUND_CEILING)
    ctx_out = Context(prec=decimal_digits, rounding=ROUND_CEILING)

    if min_prime_factor is not None and prime_factors_list is None:
        b_min = int(min_prime_factor)
        if b_min < 2:
            raise ValueError("min_prime_factor must be >= 2")
        # Number of prime factors k <= floor(log_{b_min}(modulus))
        max_k = 1
        prod = b_min
        while prod * b_min <= modulus:
            prod *= b_min
            max_k += 1
        # By Bernoulli's inequality, 1 - (1 - 1/b_min)^max_k <= max_k / b_min without
        # floating/decimal subtraction cancellation at 1024+ bit primes!
        bernoulli_hi = ctx_ceil.divide(Decimal(max_k), Decimal(b_min))
        ratio_lo = ctx_floor.power(
            ctx_floor.divide(Decimal(b_min - 1), Decimal(b_min)),
            Decimal(max_k),
        )
        direct_hi = ctx_ceil.subtract(Decimal(1), ratio_lo)
        one_minus_hi = bernoulli_hi if bernoulli_hi < direct_hi else direct_hi
        method = f"factor_free_min_prime>={b_min}"
        phi_exact: int | None = None
    else:
        phi_exact = euler_totient_from_factors(modulus, prime_factors_list)
        # 1 - phi(N)/N = (N - phi(N)) / N
        one_minus_hi = ctx_ceil.divide(Decimal(modulus - phi_exact), Decimal(modulus))
        method = "exact_totient"

    single_dtr = ctx_ceil.sqrt(one_minus_hi)
    total_dtv = ctx_ceil.multiply(Decimal(runs), single_dtr)
    if total_dtv > Decimal(1):
        total_dtv = Decimal(1)

    return {
        "modulus": modulus,
        "phi_n": phi_exact,
        "method": method,
        "single_mask_trace_distance_upper": str(+ctx_out.create_decimal(single_dtr)),
        "runs": runs,
        "total_tv_distance_upper": str(+ctx_out.create_decimal(total_dtv)),
    }


def carry_runway_mask_bound(
    runway_bits: int,
    max_added_delta: int,
    *,
    num_additions: int = 1,
    decimal_digits: int = 24,
) -> dict[str, Any]:
    """Outward-rounded trace distance bound for a `runway_bits`-qubit carry runway
    absorbing `num_additions` additions each of magnitude `<= max_added_delta`:

        `D_tr_single = sqrt(1 - (1 - delta / 2^w_pad)^2)`,
        `D_tr_total <= min(1, num_additions * D_tr_single)`.
    """
    if runway_bits < 1:
        raise ValueError("runway_bits must be >= 1")
    if max_added_delta < 0:
        raise ValueError("max_added_delta must be >= 0")
    if num_additions < 1:
        raise ValueError("num_additions must be >= 1")

    two_w = 1 << runway_bits
    if max_added_delta >= two_w:
        return {
            "runway_bits": runway_bits,
            "max_added_delta": max_added_delta,
            "num_additions": num_additions,
            "single_step_trace_distance_upper": "1",
            "total_trace_distance_upper": "1",
        }
    if max_added_delta == 0:
        return {
            "runway_bits": runway_bits,
            "max_added_delta": 0,
            "num_additions": num_additions,
            "single_step_trace_distance_upper": "0",
            "total_trace_distance_upper": "0",
        }

    ctx_floor = Context(prec=max(56, decimal_digits + 20), rounding=ROUND_FLOOR)
    ctx_ceil = Context(prec=max(56, decimal_digits + 20), rounding=ROUND_CEILING)
    ctx_out = Context(prec=decimal_digits, rounding=ROUND_CEILING)

    overlap_lo = ctx_floor.divide(Decimal(two_w - max_added_delta), Decimal(two_w))
    ov_sq_lo = ctx_floor.multiply(overlap_lo, overlap_lo)
    one_minus_hi = ctx_ceil.subtract(Decimal(1), ov_sq_lo)
    single_hi = ctx_ceil.sqrt(one_minus_hi)
    tot_hi = ctx_ceil.multiply(Decimal(num_additions), single_hi)
    if tot_hi > Decimal(1):
        tot_hi = Decimal(1)

    return {
        "runway_bits": runway_bits,
        "max_added_delta": max_added_delta,
        "num_additions": num_additions,
        "single_step_trace_distance_upper": str(+ctx_out.create_decimal(single_hi)),
        "total_trace_distance_upper": str(+ctx_out.create_decimal(tot_hi)),
    }


@dataclass(frozen=True)
class MaskFinding:
    code: str
    message: str
    witness: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "message": self.message,
            "witness": self.witness,
        }


@dataclass(frozen=True)
class MaskVerificationResult:
    valid: bool
    mask_kind: str
    bijective_on_group: bool
    worst_case_support_safe: bool
    trace_distance_upper: Decimal
    details: dict[str, Any]
    findings: tuple[MaskFinding, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "valid": self.valid,
            "mask_kind": self.mask_kind,
            "bijective_on_group": self.bijective_on_group,
            "worst_case_support_safe": self.worst_case_support_safe,
            "trace_distance_upper": str(self.trace_distance_upper),
            "details": self.details,
            "findings": [f.to_dict() for f in self.findings],
        }


def verify_mask_spec(
    mask_kind: str,
    *,
    modulus: int | None = None,
    width: int | None = None,
    multipliers: Sequence[int] = (),
    runs: int = 1,
    prime_factors_list: Sequence[int] | None = None,
    min_prime_factor: int | None = None,
    runway_bits: int | None = None,
    max_added_delta: int = 1,
    num_additions: int = 1,
    subroutine_kind: str | None = None,
    subroutine_steps: int | None = None,
) -> MaskVerificationResult:
    """Validate a superposition mask declaration against group bijection and full-support safety."""
    if mask_kind not in VALID_MASK_KINDS:
        raise ValueError(f"Unknown mask_kind {mask_kind!r}; expected one of {sorted(VALID_MASK_KINDS)}")

    findings: list[MaskFinding] = []
    bijective = True
    support_safe = True
    details: dict[str, Any] = {}
    d_tr_up = Decimal(0)

    if mask_kind in ("mult_mod_n", "uniform_mod_n"):
        if modulus is None or modulus < 2:
            raise ValueError(f"{mask_kind} requires modulus >= 2")
        n_bits = (modulus - 1).bit_length()
        details["modulus"] = modulus
        details["bit_length"] = n_bits

        for idx, g_val in enumerate(multipliers):
            g_mod = int(g_val) % modulus
            g_gcd = math.gcd(g_mod, modulus) if g_mod != 0 else modulus
            if g_mod == 0 or g_gcd != 1:
                bijective = False
                findings.append(
                    MaskFinding(
                        code="MASK-NON-BIJECTIVE-GARBAGE",
                        message=(
                            f"Multiplier g[{idx}]={g_val} has gcd({g_val}, {modulus})={g_gcd} != 1; "
                            f"r -> g*r mod N collapses cosets and entangles mask with data"
                        ),
                        witness={"multiplier_index": idx, "multiplier": g_val, "modulus": modulus, "gcd": g_gcd},
                    )
                )

        # Van Kirk & Gidney 2026 Section 5.4.1 & Appendix B:
        # Masked / Fibonacci-conjugated states populate FULL worst-case Euclidean pairs,
        # requiring `>= 2 * n - 1` iterations for Dialog/Kaliski division/multiplication
        # (average-case `1.413 * n` truncation fails!).
        if subroutine_kind in ("euclidean_dialog", "kaliski_mod_inv", "truncated_gcd"):
            req_steps = 2 * n_bits - 1
            details["required_worst_case_steps"] = req_steps
            if subroutine_steps is not None and subroutine_steps < req_steps:
                support_safe = False
                findings.append(
                    MaskFinding(
                        code="MASK-SUPPORT-OVERFLOW-HAZARD",
                        message=(
                            f"Subroutine {subroutine_kind!r} truncated to {subroutine_steps} steps < "
                            f"worst-case {req_steps} (2*n - 1); superposition-masked registers have full "
                            f"support over worst-case Euclidean chains (Van Kirk & Gidney 2026 Sec. 5.4.1)"
                        ),
                        witness={
                            "subroutine_kind": subroutine_kind,
                            "configured_steps": subroutine_steps,
                            "required_worst_case_steps": req_steps,
                            "bit_length": n_bits,
                        },
                    )
                )

        if mask_kind == "uniform_mod_n":
            ub = uniform_residue_mask_bound(
                modulus,
                runs=runs,
                prime_factors_list=prime_factors_list,
                min_prime_factor=min_prime_factor,
            )
            details.update(ub)
            d_tr_up = Decimal(ub["total_tv_distance_upper"])
        else:
            details["exact_invariant_state"] = "|M_N>"
            d_tr_up = Decimal(0)

    elif mask_kind in ("add_mod_2w", "phase_grad"):
        w = width if width is not None else 1
        if w < 1:
            raise ValueError(f"{mask_kind} requires width >= 1")
        details["width"] = w
        details["group"] = f"Z_{1 << w}"
        d_tr_up = Decimal(0)

    elif mask_kind == "carry_runway":
        rw = runway_bits if runway_bits is not None else (width if width is not None else 0)
        if rw < 1:
            raise ValueError("carry_runway requires runway_bits >= 1")
        rb = carry_runway_mask_bound(rw, max_added_delta, num_additions=num_additions)
        details.update(rb)
        d_tr_up = Decimal(rb["total_trace_distance_upper"])
        if max_added_delta >= (1 << rw):
            bijective = False
            findings.append(
                MaskFinding(
                    code="MASK-SUPPORT-OVERFLOW-HAZARD",
                    message=f"Carry delta {max_added_delta} overflows {rw}-bit runway (2^{rw}={1 << rw})",
                    witness={"runway_bits": rw, "max_added_delta": max_added_delta},
                )
            )

    return MaskVerificationResult(
        valid=(len(findings) == 0),
        mask_kind=mask_kind,
        bijective_on_group=bijective,
        worst_case_support_safe=support_safe,
        trace_distance_upper=d_tr_up,
        details=details,
        findings=tuple(findings),
    )


# ---------------------------------------------------------------------------
# Symbolic GL_2(Z) Streaming Traveling-Eigenstate / Fibonacci Ladder Verifier
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LadderStep:
    """One step in a streamed masked exponentiation / IQPE ladder.

    Supported `op` kinds:
      - `"fib_fwd"`       : `(x, y) -> (y, x * y mod N)`   (`M <- [[0, 1], [1, 1]] @ M`)
      - `"fib_inv"`       : `(x, y) -> (y * x^-1, x)`      (`M <- [[-1, 1], [1, 0]] @ M`)
      - `"swap"`          : `(x, y) -> (y, x)`             (`M <- [[0, 1], [1, 0]] @ M`)
      - `"ctrl_mul_reg0"` : controlled on streamed bit `control_id`, `x <- x * base^power mod N`
      - `"ctrl_mul_reg1"` : controlled on streamed bit `control_id`, `y <- y * base^power mod N`
      - `"bad_non_unimodular"` : e.g. out-of-place/non-invertible `[[1, 1], [1, 1]]` for negative tests
    """

    op: str
    control_id: str | None = None
    base_val: int | None = None
    power: int = 1
    matrix: tuple[tuple[int, int], tuple[int, int]] | None = None

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "LadderStep":
        mat_raw = d.get("matrix")
        mat = (
            ((int(mat_raw[0][0]), int(mat_raw[0][1])), (int(mat_raw[1][0]), int(mat_raw[1][1])))
            if mat_raw is not None
            else None
        )
        return cls(
            op=str(d["op"]),
            control_id=str(d["control_id"]) if d.get("control_id") is not None else None,
            base_val=int(d["base_val"]) if d.get("base_val") is not None else None,
            power=int(d.get("power", 1)),
            matrix=mat,
        )

    def to_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {"op": self.op, "power": self.power}
        if self.control_id is not None:
            out["control_id"] = self.control_id
        if self.base_val is not None:
            out["base_val"] = self.base_val
        if self.matrix is not None:
            out["matrix"] = [list(self.matrix[0]), list(self.matrix[1])]
        return out


@dataclass(frozen=True)
class LadderVerificationResult:
    valid: bool
    final_matrix: tuple[tuple[int, int], tuple[int, int]]
    determinant: int
    control_phase_weights: dict[str, int]
    reverse_uncompute_eliminated: bool
    trace_distance_upper: Decimal
    findings: tuple[MaskFinding, ...]
    mask_result: MaskVerificationResult

    def to_dict(self) -> dict[str, Any]:
        return {
            "valid": self.valid,
            "final_matrix": [list(self.final_matrix[0]), list(self.final_matrix[1])],
            "determinant": self.determinant,
            "control_phase_weights": dict(sorted(self.control_phase_weights.items())),
            "reverse_uncompute_eliminated": self.reverse_uncompute_eliminated,
            "trace_distance_upper": str(self.trace_distance_upper),
            "findings": [f.to_dict() for f in self.findings],
            "mask_result": self.mask_result.to_dict(),
        }


def _mat2_mul(
    a: tuple[tuple[int, int], tuple[int, int]],
    b: tuple[tuple[int, int], tuple[int, int]],
) -> tuple[tuple[int, int], tuple[int, int]]:
    return (
        (a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]),
        (a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]),
    )


def _mat2_det(m: tuple[tuple[int, int], tuple[int, int]]) -> int:
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]


def verify_streaming_fibonacci_ladder(
    modulus: int,
    eigen_base: int,
    steps: Sequence[LadderStep],
    *,
    mask_kind: str = "mult_mod_n",
    expected_control_weights: dict[str, int] | None = None,
    runs: int = 1,
    prime_factors_list: Sequence[int] | None = None,
    min_prime_factor: int | None = None,
    subroutine_kind: str | None = None,
    subroutine_steps: int | None = None,
) -> LadderVerificationResult:
    """Symbolically verify a streamed Fibonacci / traveling-eigenstate ladder over `Z_N^* x Z_N^*`.

    Tracks `(r, g^s)` exponent state `(r^{M_00} g^{M_01 s + E_0(c)}, r^{M_10} g^{M_11 s + E_1(c)})`
    and proves via unimodular change of variables `M in GL_2(Z)` (`det M == +-1`) that:
      1. Work registers land in `U_M(|M_N> |u_theta>)`, completely unentangled from every streamed
         control bit `c_k` without needing a reverse uncomputation pass!
      2. Each streamed control bit `c_k` acquires exact phase kickback `exp(2*pi*i * theta * w_k)`,
         where `w_k = det(M) * (-M_10 * E_0[c_k] + M_00 * E_1[c_k])`.
    """
    multipliers = [eigen_base] + [st.base_val for st in steps if st.base_val is not None]
    mask_res = verify_mask_spec(
        mask_kind,
        modulus=modulus,
        multipliers=multipliers,
        runs=runs,
        prime_factors_list=prime_factors_list,
        min_prime_factor=min_prime_factor,
        subroutine_kind=subroutine_kind,
        subroutine_steps=subroutine_steps,
    )

    findings: list[MaskFinding] = list(mask_res.findings)
    m_cur: tuple[tuple[int, int], tuple[int, int]] = ((1, 0), (0, 1))
    e0: dict[str, int] = {}
    e1: dict[str, int] = {}

    step_mats = {
        "fib_fwd": ((0, 1), (1, 1)),
        "fib_inv": ((-1, 1), (1, 0)),
        "swap": ((0, 1), (1, 0)),
    }

    for idx, st in enumerate(steps):
        if st.op in step_mats or st.matrix is not None:
            smat = st.matrix if st.matrix is not None else step_mats[st.op]
            det_s = _mat2_det(smat)
            if abs(det_s) != 1:
                findings.append(
                    MaskFinding(
                        code="LADDER-NON-UNIMODULAR-STEP",
                        message=f"Step {idx} ({st.op}) has determinant {det_s} != +-1; destroys group bijection",
                        witness={"step_index": idx, "matrix": [list(smat[0]), list(smat[1])], "det": det_s},
                    )
                )
            m_cur = _mat2_mul(smat, m_cur)
            all_keys = set(e0) | set(e1)
            new_e0: dict[str, int] = {}
            new_e1: dict[str, int] = {}
            for k in all_keys:
                v0 = e0.get(k, 0)
                v1 = e1.get(k, 0)
                nv0 = smat[0][0] * v0 + smat[0][1] * v1
                nv1 = smat[1][0] * v0 + smat[1][1] * v1
                if nv0 != 0:
                    new_e0[k] = nv0
                if nv1 != 0:
                    new_e1[k] = nv1
            e0, e1 = new_e0, new_e1
        elif st.op in ("ctrl_mul_reg0", "ctrl_mul_reg1"):
            cid = st.control_id or f"c{idx}"
            base = eigen_base if st.base_val is None else st.base_val
            if base % modulus != eigen_base % modulus:
                # Check if base is an explicit power of eigen_base or general coprime
                if math.gcd(base, modulus) != 1:
                    pass  # already caught in mask_res
            if st.op == "ctrl_mul_reg0":
                e0[cid] = e0.get(cid, 0) + st.power
            else:
                e1[cid] = e1.get(cid, 0) + st.power
        else:
            raise ValueError(f"Unknown LadderStep op {st.op!r}")

    det_final = _mat2_det(m_cur)
    weights: dict[str, int] = {}
    if abs(det_final) == 1:
        for cid in sorted(set(e0) | set(e1)):
            w_c = det_final * (-m_cur[1][0] * e0.get(cid, 0) + m_cur[0][0] * e1.get(cid, 0))
            weights[cid] = w_c
    else:
        findings.append(
            MaskFinding(
                code="LADDER-RESIDUAL-ENTANGLEMENT",
                message=f"Final ladder matrix has det={det_final} != +-1; streamed controls remain entangled with mask",
                witness={"final_matrix": [list(m_cur[0]), list(m_cur[1])], "det": det_final},
            )
        )

    if expected_control_weights is not None and abs(det_final) == 1:
        for cid, exp_w in expected_control_weights.items():
            got_w = weights.get(cid, 0)
            if got_w != exp_w:
                findings.append(
                    MaskFinding(
                        code="LADDER-PHASE-WEIGHT-MISMATCH",
                        message=f"Control {cid!r} kicked back exponent weight {got_w}, expected {exp_w}",
                        witness={"control_id": cid, "got_weight": got_w, "expected_weight": exp_w},
                    )
                )

    is_valid = len(findings) == 0
    return LadderVerificationResult(
        valid=is_valid,
        final_matrix=m_cur,
        determinant=det_final,
        control_phase_weights=weights,
        reverse_uncompute_eliminated=(abs(det_final) == 1 and mask_res.bijective_on_group),
        trace_distance_upper=mask_res.trace_distance_upper,
        findings=tuple(findings),
        mask_result=mask_res,
    )
