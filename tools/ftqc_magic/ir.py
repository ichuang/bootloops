"""Intermediate representations for magic-state targets and typed factory contracts.

Supports:
- PhasePolynomialTarget: canonical representation of diagonal gates in the 3rd
  (and 4th, e.g. CCCZ / sqrt(T)) level of the Clifford hierarchy, including:
  * Jain-Wills-Singh polynomial strings ('0+1', '0 1', '0 1 2', '0+1 2'),
  * named targets ('T', 'T^k', 'CS', 'CCZ', 'CCCZ', 'DET3', 'SQRT_T'),
  * exact truth-table evaluation mod 2^ell,
  * diagonal C_2 Clifford equivalence (up to Z, S, CZ, and global phase),
  * GL(k, 2) + diagonal Clifford equivalence class invariants, and
  * exact Pauli-expectation spectrum in Z[sqrt(2)] for general Clifford-orbit
    equivalence (e.g., CS_01 * CS_02 ~ CCZ_012).
- MagicSourceContract: typed source guarantee carrying correlation/poisoning tags
  (`INDEPENDENT_SINGLE`, `BATCH_CORRELATED`, `CATALYST_POISONED`) for safe composition.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from fractions import Fraction
from itertools import combinations, permutations, product
from typing import Any, Mapping, Sequence

from .rings import DEFAULT_SIG, PolyQ, QSqrt2, ZOmega8, ZSqrt2, parse_exact_fraction, round_down_str, round_up_str


class MagicRefusalError(ValueError):
    """Fail-closed refusal with machine-readable hazard code."""

    def __init__(self, code: str, message: str, witness: dict[str, Any] | None = None) -> None:
        super().__init__(f"[{code}] {message}")
        self.code = code
        self.detail = message
        self.witness = witness or {}

    def to_dict(self) -> dict[str, Any]:
        return {"error": "MagicRefusalError", "code": self.code, "detail": self.detail, "witness": self.witness}


class CorrelationMode(str, Enum):
    """Output error correlation discipline across delivered magic states."""

    INDEPENDENT_SINGLE = "INDEPENDENT_SINGLE"
    BATCH_CORRELATED = "BATCH_CORRELATED"
    CATALYST_POISONED = "CATALYST_POISONED"


@dataclass(frozen=True)
class PhasePolynomialTarget:
    """Diagonal magic gate U|x> = exp(2*pi*i * phases[x] / modulus) |x> on k qubits.

    Default `modulus = 8` covers all diagonal C_3 gates (T: +1, CS: +2, CCZ: +4)
    as well as C^r Z for any r >= 2 (phases in {0, 4} mod 8). `modulus = 16` covers sqrt(T).
    """

    k: int
    phases: tuple[int, ...]
    modulus: int = 8
    label: str = ""

    def __post_init__(self) -> None:
        if type(self.k) is not int or self.k < 1:
            raise ValueError("Target qubit count k must be a positive int")
        if type(self.modulus) is not int or self.modulus < 2 or (self.modulus & (self.modulus - 1)) != 0:
            raise ValueError("modulus must be a power of two >= 2")
        if len(self.phases) != (1 << self.k):
            raise ValueError(f"Expected {1 << self.k} phase entries for k={self.k}, got {len(self.phases)}")
        p0 = self.phases[0] % self.modulus
        normed = tuple((int(v) - p0) % self.modulus for v in self.phases)
        object.__setattr__(self, "phases", normed)

    @classmethod
    def from_jws_polynomial(cls, k: int, poly: str, *, label: str = "") -> PhasePolynomialTarget:
        """Parse Jain-Wills-Singh degree-<=3 phase polynomial string (e.g. '0+1', '0 1', '0 1 2').

        Each term has wire indices separated by spaces; '+' joins terms:
        - 1 wire  `a`     -> T_a       (+1 * x_a mod 8)
        - 2 wires `a b`   -> CS_{a,b}  (+2 * x_a * x_b mod 8)
        - 3 wires `a b c` -> CCZ_{abc} (+4 * x_a * x_b * x_c mod 8)
        """
        monomials: dict[tuple[int, ...], int] = {}
        cleaned = poly.strip()
        if cleaned:
            for raw_term in cleaned.split("+"):
                tok = raw_term.strip()
                if not tok:
                    continue
                wires = tuple(sorted(int(w) for w in tok.split()))
                if len(set(wires)) != len(wires):
                    raise ValueError(f"Duplicate wire in monomial {tok!r}")
                if not (1 <= len(wires) <= 3):
                    raise ValueError(f"JWS monomial degree must be 1..3, got {wires}")
                if any(w < 0 or w >= k for w in wires):
                    raise ValueError(f"Wire index out of range 0..{k - 1} in {tok!r}")
                monomials[wires] = (monomials.get(wires, 0) + 1) % 2
        table = [0] * (1 << k)
        for x in range(1 << k):
            acc = 0
            for wires, bit in monomials.items():
                if not bit:
                    continue
                if all((x >> w) & 1 for w in wires):
                    weight = 1 if len(wires) == 1 else (2 if len(wires) == 2 else 4)
                    acc = (acc + weight) % 8
            table[x] = acc
        return cls(k=k, phases=tuple(table), modulus=8, label=label or cleaned or "I")

    @classmethod
    def named(cls, spec: str) -> PhasePolynomialTarget:
        """Construct standard targets: 'T', 'T^k', 'CS', 'CCZ', 'CCZ^m', 'CCCZ', 'DET3', 'SQRT_T'."""
        s = spec.strip()
        u = s.upper()
        if u == "T":
            return cls.from_jws_polynomial(1, "0", label="T")
        if u.startswith("T^") or u.startswith("T**"):
            sep = "^" if "^" in u else "**"
            k = int(u.split(sep, 1)[1])
            return cls.from_jws_polynomial(k, "+".join(str(i) for i in range(k)), label=f"T^{k}")
        if u == "CS":
            return cls.from_jws_polynomial(2, "0 1", label="CS")
        if u.startswith("CS^"):
            m = int(u.split("^", 1)[1])
            terms = "+".join(f"{2 * i} {2 * i + 1}" for i in range(m))
            return cls.from_jws_polynomial(2 * m, terms, label=f"CS^{m}")
        if u == "CCZ":
            return cls.from_jws_polynomial(3, "0 1 2", label="CCZ")
        if u.startswith("CCZ^"):
            m = int(u.split("^", 1)[1])
            terms = "+".join(f"{3 * i} {3 * i + 1} {3 * i + 2}" for i in range(m))
            return cls.from_jws_polynomial(3 * m, terms, label=f"CCZ^{m}")
        if u == "CCCZ":
            phases = tuple(4 if x == 15 else 0 for x in range(16))
            return cls(k=4, phases=phases, modulus=8, label="CCCZ")
        if u in ("DET3", "TOF#", "DET"):
            # 3x3 determinant polynomial over F_2 on qubits 0..8 (row-major): 6 permutation triples.
            triples = [
                (0, 4, 8),
                (0, 5, 7),
                (1, 3, 8),
                (1, 5, 6),
                (2, 3, 7),
                (2, 4, 6),
            ]
            poly = "+".join(f"{a} {b} {c}" for a, b, c in triples)
            return cls.from_jws_polynomial(9, poly, label="DET3")
        if u in ("SQRT_T", "T_1_2", "PI_16"):
            return cls(k=1, phases=(0, 1), modulus=16, label="SQRT_T")
        raise ValueError(f"Unknown named magic target {spec!r}")

    def non_clifford_moments_mod2(self) -> dict[tuple[int, ...], int]:
        """For a modulus-8 C_3 diagonal gate, extract its canonical degree-<=3 JWS moments mod 2.

        Raises ValueError if the gate has degree >= 4 non-Clifford content (e.g. CCCZ) or modulus != 8.
        """
        if self.modulus != 8:
            raise ValueError("C_3 moment extraction requires modulus == 8")
        # Mobius / integer ANF over {0,1}^k modulo 8:
        # phases[x] = sum_{S subseteq [k]} c_S * prod_{i in S} x_i (mod 8)
        anf = self.integer_anf_mod_modulus()
        moments: dict[tuple[int, ...], int] = {}
        for mask, coeff in enumerate(anf):
            if coeff == 0:
                continue
            deg = mask.bit_count()
            wires = tuple(i for i in range(self.k) if (mask >> i) & 1)
            if deg == 0:
                continue
            elif deg == 1:
                # Odd coeff -> T; even coeff -> S/Z Clifford
                if coeff & 1:
                    moments[wires] = 1
            elif deg == 2:
                # Coeff must be even for C_3; coeff % 4 == 2 -> CS; coeff % 4 == 0 -> CZ Clifford
                if coeff & 1:
                    raise ValueError(f"Degree-2 term {wires} has non-C_3 coefficient {coeff} mod 8")
                if (coeff % 4) == 2:
                    moments[wires] = 1
            elif deg == 3:
                # Coeff must be 0 or 4 mod 8 for C_3; 4 -> CCZ
                if (coeff % 4) != 0:
                    raise ValueError(f"Degree-3 term {wires} has non-C_3 coefficient {coeff} mod 8")
                if coeff == 4:
                    moments[wires] = 1
            else:
                raise ValueError(f"Degree-{deg} term {wires} lies above C_3 (coeff={coeff} mod 8)")
        return moments

    def integer_anf_mod_modulus(self) -> list[int]:
        """Exact multilinear expansion phases[x] = sum_S c_S prod_{i in S} x_i (mod modulus)."""
        m = self.modulus
        n_states = 1 << self.k
        # Inclusion-exclusion over Z_m: c_S = sum_{T subseteq S} (-1)^(|S|-|T|) phases[T] mod m
        coeffs = list(self.phases)
        for i in range(self.k):
            step = 1 << i
            for mask in range(n_states):
                if mask & step:
                    coeffs[mask] = (coeffs[mask] - coeffs[mask ^ step]) % m
        return coeffs

    def jws_polynomial_string(self) -> str:
        moments = self.non_clifford_moments_mod2()
        ordered = sorted(moments.keys(), key=lambda t: (len(t), t))
        return "+".join(" ".join(str(w) for w in wires) for wires in ordered)

    def is_rational_dyadic(self) -> bool:
        """True iff U|+>^k amplitudes up to global phase lie in Q(i), i.e. all phases are even in Z_8
        (CS, CCZ, CCCZ, DET3, etc. - subject to Balada Gaggioli complement symmetry!).
        """
        scale = 8 if self.modulus == 8 else self.modulus
        step = scale // 4
        return all((p % step) == 0 for p in self.phases)

    def diagonal_clifford_difference(self, other: PhasePolynomialTarget) -> dict[str, Any] | None:
        """If `self` and `other` differ only by a diagonal C_2 Clifford, return the correction
        `{"global_phase_mod8": c0, "s_powers": [...], "cz_pairs": [...]}` taking `self` -> `other`;
        otherwise return None.
        """
        if self.k != other.k or self.modulus != other.modulus:
            return None
        if self.modulus != 8:
            # Lift or compare exact diff in Z_modulus: C_2 diagonal terms have unit pi/2 (modulus//4)
            # on deg 1 and pi (modulus//2) on deg 2.
            q_unit = self.modulus // 4
            z_unit = self.modulus // 2
        else:
            q_unit = 2
            z_unit = 4
        diff = PhasePolynomialTarget(
            k=self.k,
            phases=tuple((other.phases[x] - self.phases[x]) % self.modulus for x in range(1 << self.k)),
            modulus=self.modulus,
        )
        anf = diff.integer_anf_mod_modulus()
        s_powers = [0] * self.k
        cz_pairs: list[tuple[int, int]] = []
        for mask, coeff in enumerate(anf):
            if coeff == 0:
                continue
            deg = mask.bit_count()
            wires = [i for i in range(self.k) if (mask >> i) & 1]
            if deg == 1:
                if coeff % q_unit != 0:
                    return None
                s_powers[wires[0]] = (coeff // q_unit) % 4
            elif deg == 2:
                if coeff % z_unit != 0:
                    return None
                if (coeff // z_unit) % 2 == 1:
                    cz_pairs.append((wires[0], wires[1]))
            else:
                return None
        c0 = (other.phases[0] - self.phases[0]) % self.modulus
        return {
            "global_phase": c0,
            "s_powers": s_powers,
            "cz_pairs": cz_pairs,
        }

    def equivalent_up_to_diagonal_clifford(self, other: PhasePolynomialTarget) -> bool:
        return self.diagonal_clifford_difference(other) is not None

    def state_vector_zomega8(self) -> tuple[ZOmega8, ...]:
        """Exact unnormalized magic state 2^(k/2) U|+>^k in Z[omega_8]^(2^k) (requires modulus | 8)."""
        if 8 % self.modulus != 0:
            raise ValueError("ZOmega8 state vector requires modulus dividing 8")
        scale = 8 // self.modulus
        return tuple(ZOmega8.omega_pow(p * scale) for p in self.phases)

    def pauli_abs_sq_spectrum(self) -> tuple[tuple[int, int, int], ...]:
        """Exact Clifford-orbit invariant for k <= 5: sorted histogram of (|2^k <psi|P|psi>|^2 in Z[sqrt(2)])
        across all 4^k n-qubit Paulis P. Invariant under the ENTIRE n-qubit Clifford group C_2!
        """
        if self.k > 5:
            raise ValueError("Full 4^k Pauli spectrum restricted to k <= 5")
        psi = self.state_vector_zomega8()
        dim = 1 << self.k
        counts: dict[tuple[int, int], int] = {}
        for x_mask in range(dim):
            for z_mask in range(dim):
                # <psi | X(x_mask) Z(z_mask) | psi> (ignoring i^(x.z) phase since we take |.|^2!)
                acc = ZOmega8.zero()
                for ket in range(dim):
                    bra = ket ^ x_mask
                    z_sign = -1 if ((ket & z_mask).bit_count() & 1) else 1
                    term = psi[bra].conj() * psi[ket]
                    acc = acc + (term if z_sign > 0 else -term)
                norm_sq: ZSqrt2 = acc.abs_sq()
                key = (norm_sq.a, norm_sq.b)
                counts[key] = counts.get(key, 0) + 1
        return tuple((a, b, cnt) for (a, b), cnt in sorted(counts.items()))

    def cnot_diagonal_clifford_equivalent(self, other: PhasePolynomialTarget) -> bool:
        """Check equivalence under GL(k, 2) + X-bitflips + diagonal C_2 Cliffords (exact for k <= 3,
        and permutation + diagonal C_2 for k == 4). Covers e.g. CS_01 * CS_02 ~ CS * CCZ and JWS classes!
        """
        if self.k != other.k or self.modulus != other.modulus:
            return False
        if self.equivalent_up_to_diagonal_clifford(other):
            return True
        dim = 1 << self.k
        if self.k <= 3:
            # Enumerate all invertible k x k binary matrices A and bit shifts b in F_2^k.
            rows_pool = [r for r in range(1, dim)]
            for rows in permutations(rows_pool, self.k):
                # Check rank == k
                basis: list[int] = []
                ok = True
                for r in rows:
                    v = r
                    for b in basis:
                        v = min(v, v ^ b)
                    if v == 0:
                        ok = False
                        break
                    basis.append(v)
                    basis.sort(reverse=True)
                if not ok:
                    continue
                for shift in range(dim):
                    permuted = [0] * dim
                    for x in range(dim):
                        y = shift
                        for row_idx, r in enumerate(rows):
                            if ((r & x).bit_count() & 1):
                                y |= 1 << row_idx
                        permuted[x] = self.phases[y]
                    cand = PhasePolynomialTarget(k=self.k, phases=tuple(permuted), modulus=self.modulus)
                    if cand.equivalent_up_to_diagonal_clifford(other):
                        return True
            return False
        if self.k == 4:
            for perm in permutations(range(4)):
                for shift in range(dim):
                    permuted = [0] * dim
                    for x in range(dim):
                        y = shift
                        for i, pi in enumerate(perm):
                            if (x >> i) & 1:
                                y |= 1 << pi
                        permuted[x] = self.phases[y]
                    cand = PhasePolynomialTarget(k=4, phases=tuple(permuted), modulus=self.modulus)
                    if cand.equivalent_up_to_diagonal_clifford(other):
                        return True
        return False


@dataclass(frozen=True)
class MagicSourceContract:
    """Certified resource summary produced by a single factory, cultivation protocol, or pipeline."""

    name: str
    target: PhasePolynomialTarget
    n_inputs: int
    k_outputs: int
    distance: int
    active_qubits: int
    total_wires: int
    correlation_mode: CorrelationMode
    # Exact polynomials or certified rational bounds at physical/input error rate p:
    input_error_rate: Fraction
    acceptance_lower: Fraction
    conditional_infidelity_upper: Fraction
    unnormalized_failure_upper: Fraction
    expected_wire_ticks_per_batch_upper: Fraction
    leading_order_coeff: int | None = None
    joint_z_error_upper: Mapping[int, Fraction] = field(default_factory=dict)
    optimality_lower_bound_n: int | None = None
    premises: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def expected_wire_ticks_per_output_upper(self) -> Fraction:
        return self.expected_wire_ticks_per_batch_upper / self.k_outputs

    def to_receipt_summary(self, sig: int = DEFAULT_SIG) -> dict[str, Any]:
        def _safe_frac_str(q: Fraction, *, upper: bool) -> str:
            if q.numerator.bit_length() <= 1024 and q.denominator.bit_length() <= 1024:
                return f"{q.numerator}/{q.denominator}" if q.denominator != 1 else str(q.numerator)
            return round_up_str(q, sig) if upper else round_down_str(q, sig)

        gap = (
            self.n_inputs - self.optimality_lower_bound_n
            if self.optimality_lower_bound_n is not None
            else None
        )
        return {
            "name": self.name,
            "target_label": self.target.label or self.target.jws_polynomial_string(),
            "k_outputs": self.k_outputs,
            "n_inputs": self.n_inputs,
            "distance": self.distance,
            "total_wires": self.total_wires,
            "active_qubits": self.active_qubits,
            "correlation_mode": self.correlation_mode.value,
            "input_error_rate": _safe_frac_str(self.input_error_rate, upper=True),
            "acceptance_lower": _safe_frac_str(self.acceptance_lower, upper=False),
            "acceptance_lower_dec": round_down_str(self.acceptance_lower, sig),
            "conditional_infidelity_upper": _safe_frac_str(self.conditional_infidelity_upper, upper=True),
            "conditional_infidelity_upper_dec": round_up_str(self.conditional_infidelity_upper, sig),
            "unnormalized_failure_upper": _safe_frac_str(self.unnormalized_failure_upper, upper=True),
            "unnormalized_failure_upper_dec": round_up_str(self.unnormalized_failure_upper, sig),
            "expected_wire_ticks_per_output_upper": _safe_frac_str(
                self.expected_wire_ticks_per_output_upper(), upper=True
            ),
            "expected_wire_ticks_per_output_upper_dec": round_up_str(
                self.expected_wire_ticks_per_output_upper(), sig
            ),
            "leading_order_coeff": self.leading_order_coeff,
            "optimality_lower_bound_n": self.optimality_lower_bound_n,
            "optimality_gap_inputs": gap,
            "premises": list(self.premises),
        }

    def to_dict(self, sig: int = DEFAULT_SIG) -> dict[str, Any]:
        return self.to_receipt_summary(sig=sig)
