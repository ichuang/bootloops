"""Compositional HMR verification with cuts, and exact integer relations between wires.

`hmr.verify_hmr_circuit` writes every wire as one polynomial in the circuit's inputs.  That
is exact, but a ripple carry doubles at every bit, so wide arithmetic exceeds any budget.
This module proves the same obligations locally.

cut        The value of a wire is replaced by a fresh variable y, and its polynomial (in the
           inputs and earlier cut variables) is kept as the definition of y.  Everything
           after that is computed with y free.  An identity that holds for every value of
           the cut variables holds for their actual values, so every proof below is sound
           wherever the cuts are placed; cuts only decide whether a proof is found.  Cuts
           are placed automatically (a wire of degree two or more that is about to control
           a non-linear gate, or a wire whose polynomial has grown long) or by tags
           `INFO_TAG "@cut:q7"`.  A wire whose value depends on a measurement outcome is
           never cut, so definitions do not depend on outcomes.
unfolding  When an obligation is not an identity in the free variables, the most recent cut
           variable in it is replaced by its definition, and so on.  Definitions only refer
           to earlier variables, so this ends; when no cut variable is left, a polynomial
           that is not zero has an assignment of inputs and outcomes at which it is not zero.
checks     dirty resets (the reset wire is zero wherever the reset fires), scratch wires
           zero at exit, outputs independent of measurement outcomes, and cancellation of
           every outcome-dependent phase.  A branch is then injective: an outcome's
           kickback (-1)^(r q) can only be cancelled by diagonal fix-ups on the wires that
           remain, so q is a function of them and measuring it loses nothing.
relation   `prove_relation` shows that an integer combination of products of wire values
           (at entry, at exit, or at a tagged probe `INFO_TAG "@probe:name"`) is zero for
           every input: each value is lifted to an exact integer multilinear polynomial and
           the cut variables are eliminated latest first.  For carry chains each
           elimination is the identity a + b + k = (a xor b xor k) + 2 k', so the
           polynomial stays small.  A multilinear polynomial in the inputs that is zero on
           all 0/1 points has all coefficients zero, so the test at the end is exact.

100% Python standard library, exact integer and rational arithmetic only.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any, Iterable

from .gf2 import ANF, ANF_ONE, ANF_ZERO, anf_and, anf_degree, anf_not, anf_var, anf_xor
from .hmr import HMRFinding, PhasePoly, _SymbolicOverflow, _anf_to_phase_poly, _phase_add_inplace
from .ir_kickmix import ANNOTATION_OPS, KickMixCircuit, Qubit

IntPoly = dict[int, int]


@dataclass
class CutResult:
    valid: bool
    proof_method: str
    dirty_reset_free: bool
    scratch_cleaned: bool
    basis_outcome_independent: bool
    phase_cancelled: bool
    branch_injective: bool
    phase_free: bool
    num_inputs: int
    num_hmr: int
    num_cuts: int
    unfoldings: int
    largest_polynomial: int
    findings: tuple[HMRFinding, ...]
    outputs: dict[int, ANF] = field(default_factory=dict)
    probes: dict[str, dict[int, ANF]] = field(default_factory=dict)
    definitions: dict[int, ANF] = field(default_factory=dict)
    kinds: dict[int, str] = field(default_factory=dict)
    input_var: dict[int, int] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "valid": self.valid,
            "proof_method": self.proof_method,
            "dirty_reset_free": self.dirty_reset_free,
            "scratch_cleaned": self.scratch_cleaned,
            "basis_outcome_independent": self.basis_outcome_independent,
            "phase_cancelled": self.phase_cancelled,
            "branch_injective": self.branch_injective,
            "phase_free": self.phase_free,
            "num_inputs": self.num_inputs,
            "num_hmr": self.num_hmr,
            "num_cuts": self.num_cuts,
            "unfoldings": self.unfoldings,
            "largest_polynomial": self.largest_polynomial,
            "findings": [f.to_dict() for f in self.findings],
        }


def _substitute(f: ANF, var: int, definition: ANF) -> ANF:
    bit = 1 << var
    keep: set[int] = set()
    rest: set[int] = set()
    for m in f:
        if m & bit:
            rest.add(m ^ bit)
        else:
            keep.add(m)
    return anf_xor(frozenset(keep), anf_and(frozenset(rest), definition))


def int_lift(expr: ANF, budget: int) -> IntPoly:
    """The 0/1 value of an ANF as an integer multilinear polynomial."""
    cur: IntPoly = {}
    for m in sorted(expr):
        add: IntPoly = {m: 1}
        for mc, cc in cur.items():
            u = mc | m
            add[u] = add.get(u, 0) - 2 * cc
        for mk, dv in add.items():
            nv = cur.get(mk, 0) + dv
            if nv:
                cur[mk] = nv
            else:
                cur.pop(mk, None)
        if len(cur) > budget:
            raise _SymbolicOverflow()
    return cur


def int_mul(a: IntPoly, b: IntPoly, budget: int) -> IntPoly:
    out: IntPoly = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            m = ma | mb
            nv = out.get(m, 0) + ca * cb
            if nv:
                out[m] = nv
            else:
                out.pop(m, None)
        if len(out) > budget:
            raise _SymbolicOverflow()
    return out


def int_add(acc: IntPoly, term: IntPoly, scale: int = 1) -> None:
    for m, c in term.items():
        nv = acc.get(m, 0) + scale * c
        if nv:
            acc[m] = nv
        else:
            acc.pop(m, None)


class _Engine:
    def __init__(self, cut_degree: int, cut_terms: int, max_monomials: int) -> None:
        self.cut_degree, self.cut_terms, self.budget = cut_degree, cut_terms, max_monomials
        self.kinds: dict[int, str] = {}
        self.defs: dict[int, ANF] = {}
        self.unfoldings = 0
        self.largest = 0
        self.r_mask = 0
        self.y_mask = 0

    def new_var(self, kind: str) -> int:
        v = len(self.kinds)
        self.kinds[v] = kind
        if kind == "r":
            self.r_mask |= 1 << v
        elif kind == "y":
            self.y_mask |= 1 << v
        return v

    def guard(self, f: ANF) -> ANF:
        if len(f) > self.budget:
            raise _SymbolicOverflow()
        if len(f) > self.largest:
            self.largest = len(f)
        return f

    def latest_cut(self, masks: Iterable[int]) -> int | None:
        union = 0
        for m in masks:
            union |= m
        union &= self.y_mask
        return union.bit_length() - 1 if union else None

    def unfold(self, f: ANF) -> ANF:
        """f with every cut variable replaced by its definition, latest first, stopping
        early when it becomes zero."""
        while f:
            y = self.latest_cut(f)
            if y is None:
                break
            f = self.guard(_substitute(f, y, self.defs[y]))
            self.unfoldings += 1
        return f

    def witness(self, f: ANF) -> dict[str, Any]:
        """An assignment at which a cut-free non-zero polynomial is one (a minimal monomial)."""
        m = min(f, key=lambda x: (x.bit_count(), x))
        names: dict[str, int] = {}
        count = {"x": 0, "r": 0, "y": 0}
        out: dict[str, dict[str, int]] = {"inputs": {}, "hmr_outcomes": {}}
        for v in sorted(self.kinds):
            k = self.kinds[v]
            name = f"{k}{count[k]}"
            count[k] += 1
            names[name] = v
            if k == "x":
                out["inputs"][name] = (m >> v) & 1
            elif k == "r":
                out["hmr_outcomes"][name] = (m >> v) & 1
        return out


def verify_with_cuts(
    circuit: KickMixCircuit,
    *,
    inputs: Iterable[int] | None = None,
    zero_public: Iterable[int] = (),
    outputs: Iterable[int] | None = None,
    allow_branch_global_phase: bool = False,
    auto_cut: bool = True,
    cut_degree: int = 2,
    cut_terms: int = 8,
    max_monomials: int = 8192,
) -> CutResult:
    """Verify dirty resets, scratch cleanliness, outcome independence and phase cancellation of
    a KickMix circuit with cuts.

    The input qubits (`inputs`, or the register wires, or all qubits when there are no
    registers) enter with a free input variable each, except those in `zero_public`, which
    enter as |0> like every other qubit.  `outputs` are the qubits that may be non-zero at
    exit (default: the input qubits).
    """
    if inputs is not None:
        pub = tuple(inputs)
    elif circuit.registers:
        pub = tuple(w.id for r in circuit.registers if r.kind == "q" for w in r.wires if isinstance(w, Qubit))
    else:
        pub = tuple(range(circuit.num_qubits))
    zero_set = set(zero_public)
    out_set = set(pub) if outputs is None else set(outputs)
    eng = _Engine(cut_degree, cut_terms, max_monomials)
    q_anf: list[ANF] = [ANF_ZERO for _ in range(circuit.num_qubits)]
    input_var: dict[int, int] = {}
    for qid in pub:
        if qid not in zero_set:
            v = eng.new_var("x")
            input_var[qid] = v
            q_anf[qid] = anf_var(v)
    c_anf: list[ANF] = [ANF_ZERO for _ in range(circuit.num_bits)]
    phase: PhasePoly = {}
    findings: list[HMRFinding] = []
    probes: dict[str, dict[int, ANF]] = {}
    cond_stack: list[ANF] = []
    dirty_reset_free = True
    num_hmr = sum(1 for op in circuit.ops if op.kind == "HMR")
    g = eng.guard

    def cut(q: int) -> None:
        f = q_anf[q]
        if not f or f == ANF_ONE:
            return
        if len(f) == 1 and next(iter(f)).bit_count() == 1:
            return
        if any(m & eng.r_mask for m in f):
            return
        y = eng.new_var("y")
        eng.defs[y] = f
        q_anf[q] = anf_var(y)

    def tidy(q: int) -> None:
        if auto_cut and len(q_anf[q]) > eng.cut_terms:
            cut(q)

    def control(*qs: int) -> None:
        if auto_cut:
            for q in qs:
                if anf_degree(q_anf[q]) >= eng.cut_degree:
                    cut(q)

    def add_phase(fired: ANF, angle: Fraction) -> None:
        _phase_add_inplace(phase, _anf_to_phase_poly(fired, angle, eng.budget))
        if len(phase) > eng.budget:
            raise _SymbolicOverflow()

    try:
        for op_idx, op in enumerate(circuit.ops):
            kind = op.kind
            if kind == "INFO_TAG":
                tag = op.info_tag or ""
                if tag.startswith("@cut:q"):
                    cut(int(tag[6:]))
                elif tag.startswith("@probe:"):
                    probes[tag[7:]] = {q: f for q, f in enumerate(q_anf) if f}
                continue
            if kind in ANNOTATION_OPS:
                continue
            if kind == "PUSH_CONDITION":
                top = cond_stack[-1] if cond_stack else ANF_ONE
                cond_stack.append(g(anf_and(top, c_anf[op.condition])))
                continue
            if kind == "POP_CONDITION":
                if cond_stack:
                    cond_stack.pop()
                continue
            active = cond_stack[-1] if cond_stack else ANF_ONE
            if op.condition is not None:
                active = g(anf_and(active, c_anf[op.condition]))
            if not active:
                continue
            if kind == "X":
                q = op.qubits[0]
                q_anf[q] = g(anf_xor(q_anf[q], active))
            elif kind == "CX":
                c, t = op.qubits
                q_anf[t] = g(anf_xor(q_anf[t], g(anf_and(active, q_anf[c]))))
                tidy(t)
            elif kind == "CCX":
                c1, c2, t = op.qubits
                control(c1, c2)
                q_anf[t] = g(anf_xor(q_anf[t], g(anf_and(active, g(anf_and(q_anf[c1], q_anf[c2]))))))
                tidy(t)
            elif kind == "SWAP":
                q1, q2 = op.qubits
                diff = g(anf_and(active, anf_xor(q_anf[q1], q_anf[q2])))
                q_anf[q1] = g(anf_xor(q_anf[q1], diff))
                q_anf[q2] = g(anf_xor(q_anf[q2], diff))
            elif kind == "CSWAP":
                c, q1, q2 = op.qubits
                control(c)
                diff = g(anf_and(g(anf_and(active, q_anf[c])), anf_xor(q_anf[q1], q_anf[q2])))
                q_anf[q1] = g(anf_xor(q_anf[q1], diff))
                q_anf[q2] = g(anf_xor(q_anf[q2], diff))
                tidy(q1)
                tidy(q2)
            elif kind == "NEG":
                add_phase(active, Fraction(1, 1))
            elif kind == "Z":
                add_phase(g(anf_and(active, q_anf[op.qubits[0]])), Fraction(1, 1))
            elif kind == "CZ":
                c1, c2 = op.qubits
                add_phase(g(anf_and(active, g(anf_and(q_anf[c1], q_anf[c2])))), Fraction(1, 1))
            elif kind == "CCZ":
                c1, c2, c3 = op.qubits
                control(c1, c2, c3)
                fired = g(anf_and(g(anf_and(active, g(anf_and(q_anf[c1], q_anf[c2])))), q_anf[c3]))
                add_phase(fired, Fraction(1, 1))
            elif kind == "Z_POW":
                add_phase(g(anf_and(active, q_anf[op.qubits[0]])), op.angle)
            elif kind == "R":
                q = op.qubits[0]
                hazard = eng.unfold(g(anf_and(active, q_anf[q])))
                if hazard:
                    dirty_reset_free = False
                    findings.append(HMRFinding(
                        code="DIRTY-RESET-HAZARD",
                        message=f"Reset `R q{q}` at op {op_idx} can fire while q{q} == 1",
                        op_index=op_idx, qubit=q, witness=eng.witness(hazard)))
                q_anf[q] = g(anf_and(q_anf[q], anf_not(active)))
            elif kind == "HMR":
                q, b = op.qubits[0], op.bits[0]
                r = anf_var(eng.new_var("r"))
                add_phase(g(anf_and(active, g(anf_and(r, q_anf[q])))), Fraction(1, 1))
                not_act = anf_not(active)
                c_anf[b] = g(anf_xor(g(anf_and(c_anf[b], not_act)), g(anf_and(active, r))))
                q_anf[q] = g(anf_and(q_anf[q], not_act))
            elif kind == "BIT_STORE0":
                b = op.bits[0]
                c_anf[b] = g(anf_and(c_anf[b], anf_not(active)))
            elif kind == "BIT_STORE1":
                b = op.bits[0]
                not_act = anf_not(active)
                c_anf[b] = g(anf_xor(g(anf_and(c_anf[b], not_act)), active))
            elif kind == "BIT_INVERT":
                b = op.bits[0]
                c_anf[b] = g(anf_xor(c_anf[b], active))

        # scratch wires are zero at exit
        scratch_cleaned = True
        for q in range(circuit.num_qubits):
            if q in out_set:
                continue
            left = eng.unfold(q_anf[q])
            if left:
                scratch_cleaned = False
                findings.append(HMRFinding(
                    code="UNCLEANED-SCRATCH-QUBIT", message=f"Scratch qubit q{q} is non-zero at circuit exit",
                    qubit=q, witness=eng.witness(left)))
        # outputs do not depend on measurement outcomes
        basis_outcome_independent = True
        for q in sorted(out_set):
            if not any(m & eng.r_mask for m in q_anf[q]):
                continue
            dep = frozenset(m for m in eng.unfold(q_anf[q]) if m & eng.r_mask)
            if dep:
                basis_outcome_independent = False
                findings.append(HMRFinding(
                    code="OUTCOME-DEPENDENT-BASIS-OUTPUT",
                    message=f"Qubit q{q} computational basis value depends on HMR measurement outcome",
                    qubit=q, witness=eng.witness(dep)))
        # every outcome-dependent phase cancels; the rest of the phase is reported
        def unfold_phase(only_r: bool) -> None:
            while True:
                pick = [m for m in phase if m & eng.y_mask and (m & eng.r_mask or not only_r)]
                y = eng.latest_cut(pick)
                if y is None:
                    return
                bit = 1 << y
                lifted: dict[Fraction, PhasePoly] = {}
                for m in [m for m in pick if m & bit]:
                    c = phase.pop(m)
                    if c not in lifted:
                        lifted[c] = _anf_to_phase_poly(eng.defs[y], c, eng.budget)
                    rest = m ^ bit
                    _phase_add_inplace(phase, {mm | rest: cc for mm, cc in lifted[c].items()})
                eng.unfoldings += 1
                if len(phase) > eng.budget:
                    raise _SymbolicOverflow()

        unfold_phase(True)
        x_like = ~eng.r_mask
        mixed = [m for m in phase if m & eng.r_mask and m & x_like]
        pure_r = [m for m in phase if m & eng.r_mask and not m & x_like]
        phase_cancelled = True
        if mixed or (pure_r and not allow_branch_global_phase):
            phase_cancelled = False
            m_star = min(mixed or pure_r, key=lambda x: (x.bit_count(), x))
            w = eng.witness(frozenset({m_star}))
            c = phase[m_star]
            w["residual_half_turns"] = f"{c.numerator}/{c.denominator}"
            if not mixed:
                w["is_branch_global_only"] = True
            findings.append(HMRFinding(
                code="UNCOMPENSATED-HMR-PHASE",
                message="HMR phase kickback is not cancelled across measurement branches", witness=w))
        unfold_phase(False)
        phase_free = not any(not m & eng.r_mask for m in phase)
    except _SymbolicOverflow:
        return CutResult(
            valid=False, proof_method="symbolic_budget_exceeded", dirty_reset_free=False, scratch_cleaned=False,
            basis_outcome_independent=False, phase_cancelled=False, branch_injective=False, phase_free=False,
            num_inputs=len(input_var), num_hmr=num_hmr, num_cuts=len(eng.defs), unfoldings=eng.unfoldings,
            largest_polynomial=eng.largest,
            findings=(HMRFinding("SYMBOLIC-BUDGET-EXCEEDED",
                                 f"A polynomial exceeded {max_monomials} monomials even with cuts"),))
    branch_injective = dirty_reset_free and scratch_cleaned and basis_outcome_independent and phase_cancelled
    return CutResult(
        valid=branch_injective, proof_method="symbolic_anf_with_cuts", dirty_reset_free=dirty_reset_free,
        scratch_cleaned=scratch_cleaned, basis_outcome_independent=basis_outcome_independent,
        phase_cancelled=phase_cancelled, branch_injective=branch_injective, phase_free=phase_free,
        num_inputs=len(input_var), num_hmr=num_hmr, num_cuts=len(eng.defs), unfoldings=eng.unfoldings,
        largest_polynomial=eng.largest, findings=tuple(findings),
        outputs={q: q_anf[q] for q in sorted(out_set)}, probes=probes, definitions=dict(eng.defs),
        kinds=dict(eng.kinds), input_var=input_var)


def unfolded_outputs(result: CutResult, max_monomials: int = 8192) -> dict[int, ANF]:
    """Outputs as polynomials in the inputs alone (for comparison with the global method on
    small circuits)."""
    y_mask = sum(1 << v for v, k in result.kinds.items() if k == "y")
    out = {}
    for q, f in result.outputs.items():
        while True:
            union = 0
            for m in f:
                union |= m
            union &= y_mask
            if not union:
                break
            f = _substitute(f, union.bit_length() - 1, result.definitions[union.bit_length() - 1])
            if len(f) > max_monomials:
                raise _SymbolicOverflow()
        out[q] = f
    return out


@dataclass
class RelationResult:
    proved: bool
    reason: str
    terms: int
    eliminations: int
    largest_polynomial: int
    residual_monomials: int

    def to_dict(self) -> dict[str, Any]:
        return {"proved": self.proved, "reason": self.reason, "terms": self.terms,
                "eliminations": self.eliminations, "largest_polynomial": self.largest_polynomial,
                "residual_monomials": self.residual_monomials}


def prove_relation(
    result: CutResult,
    terms: Iterable[tuple[int, tuple[ANF, ...]]],
    *,
    max_monomials: int = 8192,
) -> RelationResult:
    """Show that sum_j coef_j * prod_k value(anf_jk) is zero for every input.  The ANFs are in
    the variables of `result` (inputs, cut variables); a term with no factors is a constant."""
    kinds, defs = result.kinds, result.definitions
    y_mask = sum(1 << v for v, k in kinds.items() if k == "y")
    r_mask = sum(1 << v for v, k in kinds.items() if k == "r")
    total: IntPoly = {}
    count = largest = eliminations = 0
    try:
        for coef, factors in terms:
            count += 1
            p: IntPoly = {0: 1}
            for f in factors:
                p = int_mul(p, int_lift(f, max_monomials), max_monomials)
            int_add(total, p, coef)
        largest = len(total)
        while True:
            union = 0
            for m in total:
                union |= m
            union &= y_mask
            if not union:
                break
            y = union.bit_length() - 1
            bit = 1 << y
            lifted = int_lift(defs[y], max_monomials)
            moved: IntPoly = {}
            for m in [m for m in total if m & bit]:
                moved[m ^ bit] = total.pop(m)
            int_add(total, int_mul(moved, lifted, max_monomials))
            eliminations += 1
            largest = max(largest, len(total))
            if len(total) > max_monomials:
                raise _SymbolicOverflow()
    except _SymbolicOverflow:
        return RelationResult(False, "symbolic budget exceeded", count, eliminations, largest, len(total))
    if any(m & r_mask for m in total):
        return RelationResult(False, "the relation depends on measurement outcomes", count, eliminations,
                              largest, len(total))
    if total:
        return RelationResult(False, "the relation is not an identity", count, eliminations, largest, len(total))
    return RelationResult(True, "identity after eliminating every cut variable", count, eliminations, largest, 0)


def prove_equal(result: CutResult, a: ANF, b: ANF, max_monomials: int = 8192) -> bool:
    """Whether two wire values are equal for every input and outcome (by unfolding)."""
    y_mask = sum(1 << v for v, k in result.kinds.items() if k == "y")
    f = anf_xor(a, b)
    while f:
        union = 0
        for m in f:
            union |= m
        union &= y_mask
        if not union:
            return False
        f = _substitute(f, union.bit_length() - 1, result.definitions[union.bit_length() - 1])
        if len(f) > max_monomials:
            raise _SymbolicOverflow()
    return True
