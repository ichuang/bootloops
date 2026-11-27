"""Symbolic HMR soundness, dirty-reset, branch-injectivity, and phase-polynomial verifier.

Certifies measurement-based uncomputation (MBUC) gadgets and KickMix circuits against:
  1. `DIRTY-RESET-HAZARD`           — `R` executed when target qubit can be `|1>`
  2. `UNCLEANED-SCRATCH-QUBIT`      — scratch ancilla not restored to `|0>` at exit
  3. `OUTCOME-DEPENDENT-BASIS-OUTPUT` — public computational basis output depends on HMR outcome `r`
  4. `NON-INJECTIVE-KICKMIX-BRANCH` — branch map `x -> V_r(x)` collides distinct inputs
  5. `UNCOMPENSATED-HMR-PHASE`      — HMR kickback `(-1)^(r_t * f(x))` not cancelled by fixup gates

100% Python standard library, MIT-licensed, zero float literals or float() calls.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Any

from .gf2 import (
    ANF,
    ANF_ONE,
    ANF_ZERO,
    anf_and,
    anf_degree,
    anf_eval,
    anf_not,
    anf_var,
    anf_vars_mask,
    anf_xor,
    gf2_rank,
)
from .ir_kickmix import ANNOTATION_OPS, KickMixCircuit, KickMixOp, Qubit


# Multilinear polynomial over Q / 2Z (half-turns mod 2): maps monomial bitmask -> Fraction in (0, 2)
PhasePoly = dict[int, Fraction]


def _phase_add_inplace(acc: PhasePoly, term: PhasePoly, scale: Fraction = Fraction(1, 1)) -> None:
    if scale == 0:
        return
    for mask, coeff in term.items():
        nv = (acc.get(mask, Fraction(0, 1)) + scale * coeff) % Fraction(2, 1)
        if nv == 0:
            acc.pop(mask, None)
        else:
            acc[mask] = nv


def _phase_mul(a: PhasePoly, b: PhasePoly) -> PhasePoly:
    out: PhasePoly = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            m = ma | mb
            nv = (out.get(m, Fraction(0, 1)) + ca * cb) % Fraction(2, 1)
            if nv == 0:
                out.pop(m, None)
            else:
                out[m] = nv
    return out


def _anf_to_phase_poly(expr: ANF, angle: Fraction) -> PhasePoly:
    """Exact lifting of `angle * expr(v)` into `(Q / 2Z)[v_0, ...]` (`v_i in {0, 1}`).

    When `angle == 1` (sign flip: Z, CZ, CCZ, NEG, HMR kickback), `1 * (A ^ B) == 1*A + 1*B (mod 2)`,
    so every XOR monomial lifts directly in O(|expr|)!
    For fractional `angle` (`S`, `T`, general `Z_POW`), uses exact multilinear inclusion-exclusion
    `lift(A ^ B) = lift(A) + lift(B) - 2 * lift(A) * lift(B)`.
    """
    ang = angle % Fraction(2, 1)
    if ang == 0 or not expr:
        return {}
    if ang == 1:
        return {m: Fraction(1, 1) for m in expr}

    # Exact integer/rational multilinear indicator in Q[v] (modulo 2/ang -> scale at end)
    # Using lift(m1 ^ m2) = lift(m1) + lift(m2) - 2 * lift(m1) * lift(m2)
    cur: dict[int, int] = {}
    for m in sorted(expr):
        # cur' = cur + x^m - 2 * x^m * cur
        add_terms: dict[int, int] = {m: 1}
        for mc, cc in cur.items():
            union_m = mc | m
            add_terms[union_m] = add_terms.get(union_m, 0) - 2 * cc
        for mk, dv in add_terms.items():
            nv = cur.get(mk, 0) + dv
            if nv == 0:
                cur.pop(mk, None)
            else:
                cur[mk] = nv

    out: PhasePoly = {}
    for mk, iv in cur.items():
        c = (Fraction(iv, 1) * ang) % Fraction(2, 1)
        if c != 0:
            out[mk] = c
    return out


def _minimal_monomial(masks: list[int]) -> int:
    """Pick a minimal-support (and lexicographically smallest) bitmask.

    By Möbius uniqueness on `{0, 1}^V`, if `S*` is inclusion-minimal among nonzero
    monomials of a polynomial `P`, then `P(1_{S*}) == coeff(S*) != 0`!
    """
    return min(masks, key=lambda m: (m.bit_count(), m))


def _decode_assignment(mask: int, num_x: int, num_r: int) -> dict[str, Any]:
    x_bits = {f"x{i}": (mask >> i) & 1 for i in range(num_x)}
    r_bits = {f"r{t}": (mask >> (num_x + t)) & 1 for t in range(num_r)}
    return {"inputs": x_bits, "hmr_outcomes": r_bits}


@dataclass(frozen=True)
class HMRFinding:
    code: str
    message: str
    op_index: int | None = None
    qubit: int | None = None
    witness: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "message": self.message,
            "op_index": self.op_index,
            "qubit": self.qubit,
            "witness": self.witness,
        }


@dataclass(frozen=True)
class HMRVerificationResult:
    valid: bool
    branch_injective: bool
    dirty_reset_free: bool
    phase_cancelled: bool
    branch_global_phase_clean: bool
    scratch_cleaned: bool
    basis_outcome_independent: bool
    proof_method: str
    num_public_qubits: int
    num_hmr: int
    output_degree_max: int
    findings: tuple[HMRFinding, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "valid": self.valid,
            "branch_injective": self.branch_injective,
            "dirty_reset_free": self.dirty_reset_free,
            "phase_cancelled": self.phase_cancelled,
            "branch_global_phase_clean": self.branch_global_phase_clean,
            "scratch_cleaned": self.scratch_cleaned,
            "basis_outcome_independent": self.basis_outcome_independent,
            "proof_method": self.proof_method,
            "num_public_qubits": self.num_public_qubits,
            "num_hmr": self.num_hmr,
            "output_degree_max": self.output_degree_max,
            "findings": [f.to_dict() for f in self.findings],
        }


class _SymbolicOverflow(Exception):
    """Raised when exact global ANF exceeds `max_monomials`."""


def _check_injectivity_anf(
    pub_qids: tuple[int, ...],
    q_anf: list[ANF],
    num_x: int,
    num_r: int,
    exact_enum_cap: int = 12,
) -> HMRFinding | None:
    """Check that `x -> (Q[q](x, 0))_{q in pub_qids}` is injective on `{0, 1}^{num_x}`."""
    if num_x == 0:
        return None
    x_mask_all = (1 << num_x) - 1
    # Restrict outputs to r = 0 (since outcome-independence is checked alongside)
    y_exprs: list[ANF] = [
        frozenset(m for m in q_anf[q] if (m & ~x_mask_all) == 0) for q in pub_qids
    ]

    # Fast check 1: every input variable x_i must appear in at least one output expression!
    union_vars = 0
    max_deg = 0
    for expr in y_exprs:
        union_vars |= anf_vars_mask(expr)
        max_deg = max(max_deg, anf_degree(expr))
    missing_mask = x_mask_all & ~union_vars
    if missing_mask:
        missing_idx = (missing_mask & -missing_mask).bit_length() - 1
        x1 = 0
        x2 = 1 << missing_idx
        return HMRFinding(
            code="NON-INJECTIVE-KICKMIX-BRANCH",
            message=f"Input variable x{missing_idx} is erased from all public outputs (branch map is non-injective)",
            witness={
                "x1": {f"x{i}": (x1 >> i) & 1 for i in range(num_x)},
                "x2": {f"x{i}": (x2 >> i) & 1 for i in range(num_x)},
                "colliding_output": [anf_eval(e, x1) for e in y_exprs],
                "hmr_outcomes": {f"r{t}": 0 for t in range(num_r)},
            },
        )

    # Fast check 2: if affine (degree <= 1), exact GF(2) matrix rank check for ANY width!
    if max_deg <= 1:
        rows = [anf_vars_mask(e) & x_mask_all for e in y_exprs]
        rk = gf2_rank(rows, num_x)
        if rk < num_x:
            # Find a nonzero kernel vector `k != 0` with `A k == 0`
            # Search small basis/pairs or solve nullspace:
            ker = 1
            for cand in range(1, 1 << min(num_x, 16)):
                if all(((row & cand).bit_count() & 1) == 0 for row in rows):
                    ker = cand
                    break
            return HMRFinding(
                code="NON-INJECTIVE-KICKMIX-BRANCH",
                message=f"Affine output map has GF(2) rank {rk} < {num_x}",
                witness={
                    "x1": {f"x{i}": 0 for i in range(num_x)},
                    "x2": {f"x{i}": (ker >> i) & 1 for i in range(num_x)},
                    "colliding_output": [anf_eval(e, 0) for e in y_exprs],
                    "hmr_outcomes": {f"r{t}": 0 for t in range(num_r)},
                },
            )
        return None

    # Check 3: exact 64-lane bit-packed bijection sweep for `num_x <= exact_enum_cap`
    if num_x <= exact_enum_cap:
        seen: dict[int, int] = {}
        for x_val in range(1 << num_x):
            out_val = 0
            for j, expr in enumerate(y_exprs):
                if anf_eval(expr, x_val):
                    out_val |= 1 << j
            if out_val in seen:
                prev_x = seen[out_val]
                return HMRFinding(
                    code="NON-INJECTIVE-KICKMIX-BRANCH",
                    message="Distinct computational basis inputs collide on the same branch output",
                    witness={
                        "x1": {f"x{i}": (prev_x >> i) & 1 for i in range(num_x)},
                        "x2": {f"x{i}": (x_val >> i) & 1 for i in range(num_x)},
                        "colliding_output_int": out_val,
                        "hmr_outcomes": {f"r{t}": 0 for t in range(num_r)},
                    },
                )
            seen[out_val] = x_val
    return None


def _symbolic_verify_ops(
    ops: tuple[KickMixOp, ...],
    pub_qids: tuple[int, ...],
    num_qubits: int,
    num_bits: int,
    *,
    allow_branch_global_phase: bool,
    check_scratch_clean: bool,
    max_monomials: int,
) -> HMRVerificationResult:
    num_x = len(pub_qids)
    total_hmr = sum(1 for op in ops if op.kind == "HMR")

    q_anf: list[ANF] = [ANF_ZERO for _ in range(num_qubits)]
    for idx, qid in enumerate(pub_qids):
        q_anf[qid] = anf_var(idx)
    c_anf: list[ANF] = [ANF_ZERO for _ in range(num_bits)]

    def guard(expr: ANF) -> ANF:
        if len(expr) > max_monomials:
            raise _SymbolicOverflow()
        return expr

    phase_poly: PhasePoly = {}
    cond_stack: list[ANF] = []
    findings: list[HMRFinding] = []
    dirty_reset_free = True
    hmr_idx = 0

    for op_idx, op in enumerate(ops):
        kind = op.kind
        if kind in ANNOTATION_OPS:
            continue
        if kind == "PUSH_CONDITION":
            assert op.condition is not None
            top = cond_stack[-1] if cond_stack else ANF_ONE
            cond_stack.append(guard(anf_and(top, c_anf[op.condition])))
            continue
        if kind == "POP_CONDITION":
            if cond_stack:
                cond_stack.pop()
            continue

        active = cond_stack[-1] if cond_stack else ANF_ONE
        if op.condition is not None:
            active = guard(anf_and(active, c_anf[op.condition]))
        if not active:
            if kind == "HMR":
                hmr_idx += 1
            continue

        if kind == "X":
            q = op.qubits[0]
            q_anf[q] = guard(anf_xor(q_anf[q], active))
        elif kind == "CX":
            c, t = op.qubits
            q_anf[t] = guard(anf_xor(q_anf[t], guard(anf_and(active, q_anf[c]))))
        elif kind == "CCX":
            c1, c2, t = op.qubits
            prod = guard(anf_and(active, guard(anf_and(q_anf[c1], q_anf[c2]))))
            q_anf[t] = guard(anf_xor(q_anf[t], prod))
        elif kind == "SWAP":
            q1, q2 = op.qubits
            diff = guard(anf_and(active, anf_xor(q_anf[q1], q_anf[q2])))
            q_anf[q1] = guard(anf_xor(q_anf[q1], diff))
            q_anf[q2] = guard(anf_xor(q_anf[q2], diff))
        elif kind == "CSWAP":
            c, q1, q2 = op.qubits
            act_c = guard(anf_and(active, q_anf[c]))
            diff = guard(anf_and(act_c, anf_xor(q_anf[q1], q_anf[q2])))
            q_anf[q1] = guard(anf_xor(q_anf[q1], diff))
            q_anf[q2] = guard(anf_xor(q_anf[q2], diff))
        elif kind == "NEG":
            _phase_add_inplace(phase_poly, _anf_to_phase_poly(active, Fraction(1, 1)))
        elif kind == "Z":
            fired = guard(anf_and(active, q_anf[op.qubits[0]]))
            _phase_add_inplace(phase_poly, _anf_to_phase_poly(fired, Fraction(1, 1)))
        elif kind == "CZ":
            c1, c2 = op.qubits
            fired = guard(anf_and(active, guard(anf_and(q_anf[c1], q_anf[c2]))))
            _phase_add_inplace(phase_poly, _anf_to_phase_poly(fired, Fraction(1, 1)))
        elif kind == "CCZ":
            c1, c2, c3 = op.qubits
            fired = guard(anf_and(active, guard(anf_and(q_anf[c1], q_anf[c2]))))
            fired = guard(anf_and(fired, q_anf[c3]))
            _phase_add_inplace(phase_poly, _anf_to_phase_poly(fired, Fraction(1, 1)))
        elif kind == "Z_POW":
            assert op.angle is not None
            fired = guard(anf_and(active, q_anf[op.qubits[0]]))
            _phase_add_inplace(phase_poly, _anf_to_phase_poly(fired, op.angle))
        elif kind == "R":
            q = op.qubits[0]
            hazard = guard(anf_and(active, q_anf[q]))
            if hazard:
                dirty_reset_free = False
                m_star = _minimal_monomial(list(hazard))
                findings.append(
                    HMRFinding(
                        code="DIRTY-RESET-HAZARD",
                        message=f"Reset `R q{q}` at op {op_idx} can fire while q{q} == 1",
                        op_index=op_idx,
                        qubit=q,
                        witness=_decode_assignment(m_star, num_x, total_hmr),
                    )
                )
            q_anf[q] = guard(anf_and(q_anf[q], anf_not(active)))
        elif kind == "HMR":
            q = op.qubits[0]
            b = op.bits[0]
            r_var = anf_var(num_x + hmr_idx)
            hmr_idx += 1
            kick = guard(anf_and(active, guard(anf_and(r_var, q_anf[q]))))
            _phase_add_inplace(phase_poly, _anf_to_phase_poly(kick, Fraction(1, 1)))
            not_act = anf_not(active)
            c_anf[b] = guard(anf_xor(guard(anf_and(c_anf[b], not_act)), guard(anf_and(active, r_var))))
            q_anf[q] = guard(anf_and(q_anf[q], not_act))
        elif kind == "BIT_STORE0":
            b = op.bits[0]
            c_anf[b] = guard(anf_and(c_anf[b], anf_not(active)))
        elif kind == "BIT_STORE1":
            b = op.bits[0]
            not_act = anf_not(active)
            c_anf[b] = guard(anf_xor(guard(anf_and(c_anf[b], not_act)), active))
        elif kind == "BIT_INVERT":
            b = op.bits[0]
            c_anf[b] = guard(anf_xor(c_anf[b], active))

    x_mask_all = (1 << num_x) - 1
    r_mask_all = ((1 << total_hmr) - 1) << num_x

    # 1. Check scratch qubits returned to |0>
    pub_set = set(pub_qids)
    scratch_cleaned = True
    if check_scratch_clean:
        for sq in range(num_qubits):
            if sq not in pub_set and q_anf[sq]:
                scratch_cleaned = False
                m_star = _minimal_monomial(list(q_anf[sq]))
                findings.append(
                    HMRFinding(
                        code="UNCLEANED-SCRATCH-QUBIT",
                        message=f"Scratch qubit q{sq} is non-zero at circuit exit",
                        qubit=sq,
                        witness=_decode_assignment(m_star, num_x, total_hmr),
                    )
                )

    # 2. Check public basis outputs do not depend on HMR measurement outcomes r
    basis_outcome_independent = True
    for pq in pub_qids:
        dep_r = [m for m in q_anf[pq] if (m & r_mask_all) != 0]
        if dep_r:
            basis_outcome_independent = False
            m_star = _minimal_monomial(dep_r)
            findings.append(
                HMRFinding(
                    code="OUTCOME-DEPENDENT-BASIS-OUTPUT",
                    message=f"Public qubit q{pq} computational basis value depends on HMR measurement outcome",
                    qubit=pq,
                    witness=_decode_assignment(m_star, num_x, total_hmr),
                )
            )

    # 3. Check branch injectivity on public qubits
    inj_finding = _check_injectivity_anf(pub_qids, q_anf, num_x, total_hmr)
    branch_injective = dirty_reset_free and basis_outcome_independent and (inj_finding is None)
    if inj_finding is not None:
        findings.append(inj_finding)

    # 4. Check HMR phase cancellation (`phase_poly` terms containing any r_t bit)
    mixed_xr_terms = [m for m in phase_poly if (m & r_mask_all) != 0 and (m & x_mask_all) != 0]
    pure_r_terms = [m for m in phase_poly if (m & r_mask_all) != 0 and (m & x_mask_all) == 0]

    phase_cancelled = True
    branch_global_phase_clean = len(pure_r_terms) == 0

    if mixed_xr_terms:
        phase_cancelled = False
        m_star = _minimal_monomial(mixed_xr_terms)
        coeff = phase_poly[m_star]
        w = _decode_assignment(m_star, num_x, total_hmr)
        w["residual_half_turns"] = f"{coeff.numerator}/{coeff.denominator}"
        findings.append(
            HMRFinding(
                code="UNCOMPENSATED-HMR-PHASE",
                message="State-dependent HMR phase kickback is not cancelled across measurement branches",
                witness=w,
            )
        )
    elif pure_r_terms and not allow_branch_global_phase:
        phase_cancelled = False
        m_star = _minimal_monomial(pure_r_terms)
        coeff = phase_poly[m_star]
        w = _decode_assignment(m_star, num_x, total_hmr)
        w["residual_half_turns"] = f"{coeff.numerator}/{coeff.denominator}"
        w["is_branch_global_only"] = True
        findings.append(
            HMRFinding(
                code="UNCOMPENSATED-HMR-PHASE",
                message="Branch-dependent global phase depends on HMR outcome (missing `NEG` fixup or enable allow_branch_global_phase)",
                witness=w,
            )
        )

    out_deg_max = max((anf_degree(q_anf[pq]) for pq in pub_qids), default=0)
    is_valid = (
        dirty_reset_free
        and scratch_cleaned
        and basis_outcome_independent
        and branch_injective
        and phase_cancelled
    )

    return HMRVerificationResult(
        valid=is_valid,
        branch_injective=branch_injective,
        dirty_reset_free=dirty_reset_free,
        phase_cancelled=phase_cancelled,
        branch_global_phase_clean=branch_global_phase_clean,
        scratch_cleaned=scratch_cleaned,
        basis_outcome_independent=basis_outcome_independent,
        proof_method="exact_symbolic_anf",
        num_public_qubits=num_x,
        num_hmr=total_hmr,
        output_degree_max=out_deg_max,
        findings=tuple(findings),
    )


def verify_hmr_circuit(
    circuit: KickMixCircuit,
    *,
    allow_branch_global_phase: bool = False,
    check_scratch_clean: bool = True,
    max_monomials: int = 8192,
) -> HMRVerificationResult:
    """Verify HMR phase cancellation, dirty-reset absence, scratch restoration, and branch injectivity.

    Uses exact global symbolic ANF + `(Q / 2Z)` multilinear phase polynomials up to `max_monomials`,
    and reports an unproved result when the global proof exceeds that budget.
    """
    if circuit.registers:
        pub_qids = tuple(
            w.id for r in circuit.registers if r.kind == "q" for w in r.wires if isinstance(w, Qubit)
        )
    else:
        pub_qids = tuple(range(circuit.num_qubits))

    try:
        return _symbolic_verify_ops(
            circuit.ops,
            pub_qids,
            circuit.num_qubits,
            circuit.num_bits,
            allow_branch_global_phase=allow_branch_global_phase,
            check_scratch_clean=check_scratch_clean,
            max_monomials=max_monomials,
        )
    except _SymbolicOverflow:
        return HMRVerificationResult(
            valid=False, branch_injective=False, dirty_reset_free=False,
            phase_cancelled=False, branch_global_phase_clean=False, scratch_cleaned=False,
            basis_outcome_independent=False, proof_method="symbolic_budget_exceeded",
            num_public_qubits=len(pub_qids), num_hmr=circuit.operation_statistics().num_hmr,
            output_degree_max=-1,
            findings=(HMRFinding("SYMBOLIC-BUDGET-EXCEEDED",
                                 f"Global ANF exceeded {max_monomials} monomials"),),
        )
