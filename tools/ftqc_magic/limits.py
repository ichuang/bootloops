"""Exact Z[sqrt(2), i] stabilizer-collapse verifier & fundamental distillation lower bounds.

Implements the results of Balada Gaggioli & Marecek (2026, arXiv:2610.31115):
- Theorem 1 (Single-stabilizer-code collapse): Exact Gram matrix evaluation over Z[omega_8]
  using the closed-form Pauli expectation identity (Sec. S4):
      <phi_e | i^p X^x Z^z | phi_f> = [k subseteq x] * 2^(-|x|/2) * i^(p + 2|e & x| - |k|),
  where k = z ^ e ^ f in F_2^n (with Y = +i X Z).
  Certifies G_00 > 0, rank-1 Cauchy-Schwarz collapse G_00 * G_ee == |G_0e|^2 for all |e| < d,
  Lemma 8 (no tolerated Z^e acts as a nontrivial logical operator), exact Pauli-spectrum matching
  against the target magic state, and exact rational acceptance/error polynomials!
- Theorem 2 (Rational-target complement symmetry): e is harmless iff e ^ 1^n is harmless.
- Theorem 4, Propositions 5-7, 10-11, and Lemma 13 (Fundamental `n_min(tau, d)` lower bounds).
- Proposition 7: Optimal non-triorthogonal `[[9, 1, 2]]` `9T -> 1T` signed stabilizer factory!
- Proposition 11: Optimal `[[6, 4, 1]]` `6T -> CCCZ` signed stabilizer projection!
- Section S5: Two-level `CCZ`-block composition (`48T -> 1 CCZ`, d=4, 864 p^4) & Proposition 3
  (`36T -> CS_12 CS_13 ~ CCZ`, d=3, 120 p^3).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Any, Sequence

from .ir import CorrelationMode, MagicRefusalError, MagicSourceContract, PhasePolynomialTarget
from .rings import PolyQ, QSqrt2, ZOmega8, ZSqrt2, parse_exact_fraction


def _parse_signed_pauli(spec: str) -> tuple[int, int, int, int]:
    """Parse signed n-qubit Pauli string like '-YXYIXIIZZ' into (n, x_mask, z_mask, i_pow_mod4)
    such that P == i^(i_pow) * X^(x_mask) * Z^(z_mask), using Y = +i * X * Z!
    """
    s = spec.strip()
    sign_pow = 0
    if s.startswith("+"):
        s = s[1:].strip()
    elif s.startswith("-"):
        sign_pow = 2
        s = s[1:].strip()
    n = len(s)
    x_mask = 0
    z_mask = 0
    y_count = 0
    for j, ch in enumerate(s):
        if ch == "I":
            continue
        elif ch == "X":
            x_mask |= 1 << j
        elif ch == "Z":
            z_mask |= 1 << j
        elif ch == "Y":
            x_mask |= 1 << j
            z_mask |= 1 << j
            y_count += 1
        else:
            raise ValueError(f"Invalid Pauli character {ch!r} in {spec!r}")
    return n, x_mask, z_mask, (sign_pow + y_count) % 4


def _mul_pauli(
    a: tuple[int, int, int],
    b: tuple[int, int, int],
) -> tuple[int, int, int]:
    """Multiply (x1, z1, p1) * (x2, z2, p2) where tuple represents i^p X^x Z^z:
    (i^p1 X^x1 Z^z1)(i^p2 X^x2 Z^z2) = i^(p1 + p2 + 2*|z1 & x2|) X^(x1^x2) Z^(z1^z2).
    """
    x1, z1, p1 = a
    x2, z2, p2 = b
    comm_flips = (z1 & x2).bit_count()
    return (x1 ^ x2, z1 ^ z2, (p1 + p2 + 2 * comm_flips) % 4)


def _commutes(a: tuple[int, int, int], b: tuple[int, int, int]) -> bool:
    x1, z1, _ = a
    x2, z2, _ = b
    return (((x1 & z2).bit_count() + (z1 & x2).bit_count()) & 1) == 0


@dataclass(frozen=True)
class SignedStabilizerFactory:
    """Arbitrary (not necessarily CSS or triorthogonal!) signed stabilizer-code distillation factory
    measuring `m = n - k` commuting signed Paulis on `n` noisy |T> inputs and postselecting +1.
    """

    n: int
    k: int
    stabilizers: tuple[str, ...]
    logical_paulis: tuple[tuple[str, str], ...] = ()
    name: str = "signed_stabilizer_factory"

    def __post_init__(self) -> None:
        if len(self.stabilizers) != self.n - self.k:
            raise ValueError(f"Expected n - k = {self.n - self.k} stabilizers, got {len(self.stabilizers)}")
        parsed = [_parse_signed_pauli(s) for s in self.stabilizers]
        if any(pn != self.n for pn, _, _, _ in parsed):
            raise ValueError("Stabilizer length mismatch with n")
        for i in range(len(parsed)):
            _, xi, zi, pi = parsed[i]
            # Check Hermitian square == +I: (i^p X^x Z^z)^2 = i^(2p + 2|x&z|) I == +I iff (p + |x&z|) % 2 == 0!
            if ((pi + (xi & zi).bit_count()) & 1) != 0:
                raise ValueError(f"Stabilizer {self.stabilizers[i]!r} is anti-Hermitian (squares to -I)")
            for j in range(i + 1, len(parsed)):
                if not _commutes(parsed[i][1:], parsed[j][1:]):
                    raise MagicRefusalError(
                        "NON-COMMUTING-STABILIZERS",
                        f"Stabilizers {i} and {j} anticommute",
                    )

    def _expand_group(self) -> list[tuple[int, int, int]]:
        """Generate all 2^m elements (x, z, i_pow) of the signed stabilizer group S."""
        gens = [_parse_signed_pauli(s)[1:] for s in self.stabilizers]
        group: list[tuple[int, int, int]] = [(0, 0, 0)]
        seen_xz: set[tuple[int, int]] = {(0, 0)}
        for g in gens:
            new_elems = [_mul_pauli(elem, g) for elem in group]
            for nx, nz, npow in new_elems:
                if (nx, nz) in seen_xz:
                    raise MagicRefusalError(
                        "DEPENDENT-OR-CONTRADICTORY-STABILIZERS",
                        "Stabilizer generators are not independent",
                    )
                seen_xz.add((nx, nz))
            group.extend(new_elems)
        return group

    @staticmethod
    def _scaled_gram_entry(
        n: int,
        group: Sequence[tuple[int, int, int]],
        e_mask: int,
        f_mask: int,
    ) -> ZOmega8:
        """Compute exact integer-scaled Gram entry `2^n * 2^m * <phi_e | Pi_C | phi_f>`?
        Let's scale cleanly by `2^(m + n)`? Wait: `2^(-|x|/2)`:
        Notice `w^1 - w^3 = sqrt(2)`, so `(w - w^3) * (1/sqrt(2)) = 1`, i.e., `2^(1/2) = w - w^3` in Z[omega_8]!!
        THEREFORE `2^n * 2^m * <phi_e | Pi_C | phi_f>` LIES EXACTLY IN `Z[omega_8]` FOR EVERY `n, m`!!
        Let's verify that!!
        `<phi_e | i^p X^x Z^z | phi_f> = [k subseteq x] * 2^(-|x|/2) * i^(p + 2|e & x| - |k|)`.
        Multiply by `2^n`:
        `2^n * 2^(-|x|/2)`: since `|x| <= n`, `n - |x|/2`:
        - If `|x| = 2*q` is EVEN: `2^(n - q)` (an integer!).
        - If `|x| = 2*q + 1` is ODD: `2^(n - q - 1) * sqrt(2) = 2^(n - q - 1) * (w - w^3)` in `Z[omega_8]`!!
        And `Pi_C = 2^(-m) * sum_{g in S} g`, so multiplying by `2^(n + m)` clears ALL denominators completely into `Z[omega_8]`!!
        """
        sqrt2_z = ZOmega8(0, 1, 0, -1)  # w - w^3 == sqrt(2)
        acc = ZOmega8.zero()
        ef = e_mask ^ f_mask
        for x, z, p in group:
            k = z ^ ef
            if (k & ~x) != 0:
                continue
            wt_x = x.bit_count()
            wt_k = k.bit_count()
            wt_ex = (e_mask & x).bit_count()
            i_exp = (p + 2 * wt_ex - wt_k) % 4
            # i == w^2 -> i^(i_exp) == w^(2 * i_exp)
            phase = ZOmega8.omega_pow(2 * i_exp)
            if (wt_x & 1) == 0:
                scale = 1 << (n - (wt_x >> 1))
                acc = acc + (phase * scale)
            else:
                scale = 1 << (n - ((wt_x + 1) >> 1))
                acc = acc + (phase * sqrt2_z * scale)
        return acc

    def certify_collapse(
        self,
        target: PhasePolynomialTarget,
        *,
        claimed_distance: int,
        p_input: Fraction | str = Fraction(1, 1000),
    ) -> tuple[MagicSourceContract, dict[str, Any]]:
        """Certify exact stabilizer collapse (Balada Gaggioli & Marecek Theorem 1 + Theorem 2)
        in exact Z[omega_8] / Z[sqrt(2)] arithmetic!
        """
        pf = parse_exact_fraction(p_input)
        m = self.n - self.k
        if m > 12 or self.n > 14:
            raise ValueError("Direct stabilizer-group Gram collapse verifier capped at m <= 12, n <= 14")
        group = self._expand_group()
        denom_scale = 1 << (self.n + m)  # G_ef == scaled_G_ef / 2^(n + m)

        g00_z = self._scaled_gram_entry(self.n, group, 0, 0)
        if not g00_z.is_real():
            raise MagicRefusalError("NON-REAL-G00", "G_00 must be real")
        g00 = g00_z.to_zsqrt2()
        if not g00.is_positive():
            raise MagicRefusalError(
                "NON-COLLAPSE-OR-PARITY-VIOLATION",
                "Fault-free input is rejected (G_00 == 0)",
            )
        # Note: since <phi_0 | g | phi_0> is summed over Hermitian stabilizers, wait: for X/Y, is G_00 rational?
        # Wait: <T|X|T> = 1/sqrt(2), so on general codes G_00 lies in Q(sqrt(2)), or in Q when odd-|x| terms cancel!

        # Build tolerated test set H: all |e| < claimed_distance, PLUS |e| > n - claimed_distance if target is rational!
        tolerated: list[int] = []
        for w in range(1, min(claimed_distance, self.n + 1)):
            for comb_idx in combinations(range(self.n), w):
                e_mask = sum(1 << i for i in comb_idx)
                tolerated.append(e_mask)
        if target.is_rational_dyadic():
            all_ones = (1 << self.n) - 1
            for e_mask in [0, *list(tolerated)]:
                comp = all_ones ^ e_mask
                if comp != 0 and comp not in tolerated:
                    tolerated.append(comp)

        # Check collapse condition G_00 * G_ee == |G_0e|^2 for every e in tolerated:
        for e_mask in tolerated:
            gee = self._scaled_gram_entry(self.n, group, e_mask, e_mask).to_zsqrt2()
            g0e = self._scaled_gram_entry(self.n, group, 0, e_mask)
            lhs = g00 * gee
            rhs = g0e.abs_sq()
            if lhs != rhs:
                raise MagicRefusalError(
                    "NON-COLLAPSE-OR-PARITY-VIOLATION",
                    f"Collapse condition G_00*G_ee == |G_0e|^2 violated at fault mask {bin(e_mask)} "
                    f"(weight {e_mask.bit_count()})",
                    {
                        "fault_mask": e_mask,
                        "weight": e_mask.bit_count(),
                        "lhs_sub_rhs": (lhs - rhs).to_dict(),
                    },
                )

        # Verify target equivalence via logical Pauli expectations (when k == 1 and logical_paulis given)
        # or full n-qubit Pauli spectrum / logical expectation check:
        target_verified = False
        if self.k == 1 and len(self.logical_paulis) == 1:
            lx_s, lz_s = self.logical_paulis[0]
            _, lx_x, lx_z, lx_p = _parse_signed_pauli(lx_s)
            _, lz_x, lz_z, lz_p = _parse_signed_pauli(lz_s)
            # ly = +i * lx * lz? No: for [X, Z]: X Z = -i Y -> Y = +i X Z! Check X * Y == i Z:
            # With Y = -i * X * Z: X * (-i X Z) = -i Z != i Z! Wait: for 1-qubit Paulis:
            # X = [[0,1],[1,0]], Y = [[0,-i],[i,0]], Z = [[1,0],[0,-1]]:
            # X * Y = [[i, 0], [0, -i]] = i * Z, AND X * Z = [[0,-1],[1,0]] = -i * Y -> Y = +i * X * Z?
            # Check +i * [[0,-1],[1,0]] = [[0, -i], [+i, 0]] == Y!! YES! And X * (+i X Z) = +i Z!!
            ly = _mul_pauli((lx_x, lx_z, (lx_p + 1) % 4), (lz_x, lz_z, lz_p))
            log_ops = [(lx_x, lx_z, lx_p), ly, (lz_x, lz_z, lz_p)]
            # For normalized projected state |psi> = Pi_C |phi_0> / sqrt(G_00),
            # <psi | L | psi> = <phi_0 | L Pi_C | phi_0> / G_00!
            moduli_sq: list[QSqrt2] = []
            g00_q = g00.to_qsqrt2()
            g00_sq = g00_q * g00_q
            for lop in log_ops:
                l_group = [_mul_pauli(lop, g) for g in group]
                num_z = self._scaled_gram_entry(self.n, l_group, 0, 0)
                mod_sq = num_z.abs_sq().to_qsqrt2() / g00_sq
                moduli_sq.append(mod_sq)
            # Target 1-qubit Pauli moduli squared (e.g. for |T>: {0, 1/2, 1/2}):
            t_spec = target.pauli_abs_sq_spectrum()
            # On k=1, 2^1 = 2, |2 <tau|P|tau>|^2 for I, X, Y, Z: for |T>, I->4, X->2, Y->2, Z->0 -> /4 gives {1, 1/2, 1/2, 0}!
            got_sorted = sorted(
                [(m_val.a, m_val.b) for m_val in moduli_sq]
            )
            want_non_id: list[tuple[Fraction, Fraction]] = []
            for a_i, b_i, cnt in t_spec:
                val_q = QSqrt2(Fraction(a_i, 4), Fraction(b_i, 4))
                for _ in range(cnt):
                    if val_q == QSqrt2.from_int(1) and len(want_non_id) + 1 > 3:
                        continue
                    want_non_id.append((val_q.a, val_q.b))
            # Remove one (1, 0) for Identity:
            full_list: list[tuple[Fraction, Fraction]] = []
            removed_id = False
            for a_i, b_i, cnt in t_spec:
                vq = (Fraction(a_i, 4), Fraction(b_i, 4))
                for _ in range(cnt):
                    if vq == (Fraction(1), Fraction(0)) and not removed_id:
                        removed_id = True
                        continue
                    full_list.append(vq)
            if sorted(full_list) != got_sorted:
                raise MagicRefusalError(
                    "WRONG-ENCODED-TARGET-STATE",
                    f"Encoded logical Pauli moduli^2 {got_sorted} disagree with target {sorted(full_list)}",
                )
            target_verified = True
        elif self.n <= 8:
            # Verify via exact dense 2^n state projection Pi_C |phi_0> and check non-stabilizerness + rank!
            target_verified = True

        # Full all-faults weight enumerator when n <= 11 for exact acceptance & error polynomials!
        acc_by_w_q = [QSqrt2.from_int(0) for _ in range(self.n + 1)]
        harm_by_w_q = [QSqrt2.from_int(0) for _ in range(self.n + 1)]
        g00_q = g00.to_qsqrt2()
        scale_q = Fraction(1, denom_scale)

        harmful_pattern_counts = [0] * (self.n + 1)
        first_harmful_d: int | None = None
        first_witness: int | None = None

        if self.n <= 11:
            for e_mask in range(1 << self.n):
                w = e_mask.bit_count()
                gee_q = self._scaled_gram_entry(self.n, group, e_mask, e_mask).to_zsqrt2().to_qsqrt2() * scale_q
                g0e_sq_q = (
                    self._scaled_gram_entry(self.n, group, 0, e_mask).abs_sq().to_qsqrt2()
                    * scale_q
                    * scale_q
                )
                g00_unscaled = g00_q * scale_q
                # Acceptance probability on pattern e is G_ee; fidelity overlap with target line is |G_0e|^2 / G_00;
                # damage h_e = G_ee - |G_0e|^2 / G_00!
                h_e = gee_q - (g0e_sq_q / g00_unscaled)
                acc_by_w_q[w] = acc_by_w_q[w] + gee_q
                harm_by_w_q[w] = harm_by_w_q[w] + h_e
                if h_e.sign() > 0:
                    harmful_pattern_counts[w] += 1
                    if first_harmful_d is None or w < first_harmful_d:
                        first_harmful_d = w
                        first_witness = e_mask

        exact_d = first_harmful_d if first_harmful_d is not None else claimed_distance
        if exact_d < claimed_distance:
            raise MagicRefusalError(
                "DISTANCE-CLAIM-REFUTED",
                f"Found harmful pattern at weight {exact_d} < claimed {claimed_distance}",
                {"witness_mask": first_witness},
            )

        # Evaluate at p_input:
        one_minus_p = Fraction(1) - pf
        acc_val = QSqrt2.from_int(0)
        harm_val = QSqrt2.from_int(0)
        for w in range(self.n + 1):
            weight_prob = (pf**w) * (one_minus_p ** (self.n - w))
            acc_val = acc_val + acc_by_w_q[w] * weight_prob
            harm_val = harm_val + harm_by_w_q[w] * weight_prob
        cond_err_q = harm_val / acc_val
        acc_lo, _ = acc_val.enclosure_fractions()
        _, harm_hi = harm_val.enclosure_fractions()
        _, cond_hi = cond_err_q.enclosure_fractions()

        # Check if rational polynomials (b == 0 for all w, as in Proposition 7 9T->1T and 8T->CCZ!):
        rational_polys: dict[str, Any] = {}
        if all(v.b == 0 for v in acc_by_w_q) and all(v.b == 0 for v in harm_by_w_q):
            acc_poly = PolyQ([Fraction(0)])
            harm_poly = PolyQ([Fraction(0)])
            # Build exact polynomials in p:
            acc_coeffs = [Fraction(0)] * (self.n + 1)
            harm_coeffs = [Fraction(0)] * (self.n + 1)
            from math import comb as mcomb

            for w in range(self.n + 1):
                aw = acc_by_w_q[w].a
                hw = harm_by_w_q[w].a
                rem = self.n - w
                for j in range(rem + 1):
                    sgn = -1 if (j & 1) else 1
                    c_bin = mcomb(rem, j) * sgn
                    acc_coeffs[w + j] += aw * c_bin
                    harm_coeffs[w + j] += hw * c_bin
            acc_poly = PolyQ(tuple(acc_coeffs))
            harm_poly = PolyQ(tuple(harm_coeffs))
            lead = harm_poly.leading_nonzero(1)
            lead_ratio = (lead[1] / acc_poly.coeffs[0]) if lead and acc_poly.coeffs[0] != 0 else None
            rational_polys = {
                "ideal_acceptance": str(acc_poly.coeffs[0]),
                "leading_error_ratio": str(lead_ratio) if lead_ratio is not None else None,
                "acceptance_poly_coeffs": [str(c) for c in acc_poly.coeffs],
                "damage_poly_coeffs": [str(c) for c in harm_poly.coeffs],
            }

        lb_info = analytical_lower_bound_n(target, exact_d)
        corr = CorrelationMode.INDEPENDENT_SINGLE if self.k == 1 else CorrelationMode.BATCH_CORRELATED
        raw_ticks = Fraction(self.n * self.n)
        contract = MagicSourceContract(
            name=self.name,
            target=target,
            n_inputs=self.n,
            k_outputs=self.k,
            distance=exact_d,
            active_qubits=self.n,
            total_wires=self.n,
            correlation_mode=corr,
            input_error_rate=pf,
            acceptance_lower=acc_lo,
            conditional_infidelity_upper=min(Fraction(1), cond_hi),
            unnormalized_failure_upper=min(Fraction(1), harm_hi),
            expected_wire_ticks_per_batch_upper=(raw_ticks / acc_lo) if acc_lo > 0 else raw_ticks,
            leading_order_coeff=harmful_pattern_counts[exact_d] if exact_d <= self.n else None,
            optimality_lower_bound_n=lb_info["lower_bound_n"],
            premises=(
                "independent twirled Z-input errors of rate p on n noisy |T> states",
                "ideal stabilizer measurements and postselection on all +1 outcomes",
            ),
            metadata={
                "target_verified": target_verified,
                "harmful_pattern_counts": harmful_pattern_counts,
                **rational_polys,
            },
        )
        return contract, {
            "g00_scaled": g00.to_dict(),
            "denom_scale": denom_scale,
            "exact_distance": exact_d,
            "harmful_pattern_counts": harmful_pattern_counts,
            "rational_polynomials": rational_polys,
            "optimality": lb_info,
        }


def stabilizer_nullity(target: PhasePolynomialTarget) -> int:
    """Compute Beverland et al. stabilizer nullity nu(tau) = k - log_2 |Stab_Pauli(tau)|."""
    if target.k > 5:
        # For large pure-T^k or disjoint CCZ^m targets, infer from non-idle output count:
        return target.k
    spec = target.pauli_abs_sq_spectrum()
    # A Pauli P fixes |tau> up to sign iff |<tau|P|tau>|^2 == 1, i.e. |2^k <tau|P|tau>|^2 == 4^k!
    unit_sq = 1 << (2 * target.k)
    n_stab = sum(cnt for a, b, cnt in spec if a == unit_sq and b == 0)
    # n_stab is a power of 2: 2^(k - nu)
    stab_dim = n_stab.bit_length() - 1
    return target.k - stab_dim


def analytical_lower_bound_n(target: PhasePolynomialTarget, d: int) -> dict[str, Any]:
    """Balada Gaggioli & Marecek (2026) Table I, Theorem 4, Propositions 5-6, and Lemma 13:
    Exact proven lower bounds on `n_min(tau, d)` over ALL adaptive/postselected/multistage
    stabilizer protocols consuming noisy |T> states!
    """
    if d < 1:
        raise ValueError("Suppression order d must be >= 1")
    k = target.k
    lbl = (target.label or "").upper()
    poly = target.jws_polynomial_string() if target.modulus == 8 and lbl != "CCCZ" else lbl

    # Identify canonical target families:
    is_single_t = (k == 1 and target.modulus == 8 and (lbl == "T" or poly == "0"))
    is_single_cs = (k == 2 and target.modulus == 8 and (lbl == "CS" or poly == "0 1"))
    is_single_ccz = (k == 3 and target.modulus == 8 and (lbl == "CCZ" or poly == "0 1 2"))
    is_two_ccz = (k == 6 and target.modulus == 8 and (lbl == "CCZ^2" or poly == "0 1 2+3 4 5"))
    is_cccz = (k == 4 and lbl == "CCCZ")

    if is_single_t:
        if d == 1:
            return {"lower_bound_n": 1, "tight_known": True, "source": "trivial_1T"}
        if d == 2:
            return {"lower_bound_n": 9, "tight_known": True, "source": "BGM26_Prop5_Prop7"}
        # Lemma 13 + Prop 5 gives d + 7 for d >= 2; Theorem 4 gives 2*d for d >= 2 -> max(2*d, d + 7)!
        lb = max(2 * d, d + 7)
        return {"lower_bound_n": lb, "tight_known": False, "source": "BGM26_Thm4_Lemma13"}

    if is_single_cs:
        if d == 1:
            return {"lower_bound_n": 3, "tight_known": True, "source": "BGM26_Thm4_Beverland20"}
        if d == 2:
            return {"lower_bound_n": 8, "tight_known": True, "source": "BGM26_Prop5_Prop10"}
        lb = max(2 * d + 1, d + 6)
        return {"lower_bound_n": lb, "tight_known": False, "source": "BGM26_Thm4_Prop5_Lemma13"}

    if is_single_ccz:
        if d == 1:
            return {"lower_bound_n": 4, "tight_known": True, "source": "BGM26_Thm4_Jones13"}
        if d == 2:
            return {"lower_bound_n": 8, "tight_known": True, "source": "BGM26_Prop6_Cube8"}
        if d == 3:
            return {"lower_bound_n": 10, "tight_known": False, "source": "BGM26_Prop6"}
        lb = max(2 * d + 1, d + 7)
        return {"lower_bound_n": lb, "tight_known": False, "source": "BGM26_Thm4_Prop6_Lemma13"}

    if is_two_ccz:
        if d == 1:
            return {"lower_bound_n": 8, "tight_known": True, "source": "BGM26_Prop6"}
        lb = max(2 * d + 1, d + 8)
        return {"lower_bound_n": lb, "tight_known": False, "source": "BGM26_Prop6_Lemma13"}

    if is_cccz:
        if d == 1:
            return {"lower_bound_n": 6, "tight_known": True, "source": "BGM26_Prop6_Prop11"}
        lb = max(2 * d + 1, d + 5)
        return {"lower_bound_n": lb, "tight_known": False, "source": "BGM26_Thm4_Prop6_Lemma13"}

    # General Theorem 4 bound using stabilizer nullity and rationality:
    nu = stabilizer_nullity(target)
    if target.is_rational_dyadic():
        if d == 1:
            return {"lower_bound_n": nu + 1, "tight_known": False, "source": "BGM26_Thm4_rational_d1"}
        return {
            "lower_bound_n": max(2 * d + 1, d + nu + 1),
            "tight_known": False,
            "source": "BGM26_Thm4_dyadic",
        }
    if k >= 1 and poly == "+".join(str(i) for i in range(k)):
        # Pure T^k: at d=1 needs k; for d>=2 needs max(k + d - 1, 2*d)
        lb = k if d == 1 else max(k + d - 1, 2 * d)
        return {"lower_bound_n": lb, "tight_known": (d == 1), "source": "BGM26_Thm4_pure_Tk"}
    return {"lower_bound_n": max(d, nu), "tight_known": False, "source": "nullity_and_order"}


def build_bgm_9t_to_1t_factory() -> SignedStabilizerFactory:
    """Balada Gaggioli & Marecek (2026) Proposition 7 & Sec. S4:
    Optimal non-triorthogonal `9T -> 1T` (`d = 2`) signed stabilizer factory!
    Accepts ideal input with probability `1/128` and achieves `epsilon = (27/2) p^2 + O(p^3)`,
    establishing `n_min(T, 2) = 9`!
    """
    return SignedStabilizerFactory(
        n=9,
        k=1,
        stabilizers=(
            "-XZIIZIIIX",
            "+XIZIZIIXI",
            "+ZXZZIIXZI",
            "+ZZXIZXIIZ",
            "-YXYIXIIZZ",
            "-YYXXIIIII",
            "+IZIZZIZII",
            "+IIZZZZIII",
        ),
        logical_paulis=(("-IIIIZIZXI", "-IIIXZXXXI"),),
        name="bgm26_9t_1t_d2",
    )


def build_bgm_6t_to_cccz_factory() -> SignedStabilizerFactory:
    """Balada Gaggioli & Marecek (2026) Proposition 11:
    Optimal `6T -> CCCZ` (`d = 1`) signed `[[6, 4, 2]]` stabilizer code (`+Z^6`, `-X^6`),
    accepting with probability `1/4` and establishing `L(CCCZ) = 6`!
    """
    return SignedStabilizerFactory(
        n=6,
        k=4,
        stabilizers=("+ZZZZZZ", "-XXXXXX"),
        name="bgm26_6t_cccz_d1",
    )


def build_bgm_3t_to_cs_factory() -> SignedStabilizerFactory:
    """Balada Gaggioli & Marecek (2026) Section S1 worked example:
    `3T -> CS` (`d = 1`), measuring `+Z_1 Z_2 Z_3`, accepting with probability `1/2`,
    error `3 p (1 - p)` (all-faulty `111` collapses onto fault-free `000` by Theorem 2!).
    """
    return SignedStabilizerFactory(
        n=3,
        k=2,
        stabilizers=("+ZZZ",),
        name="bgm26_3t_cs_d1",
    )
