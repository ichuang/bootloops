"""Disjoint-lineage exact Kraus / dense-statevector channel verifier (`n <= 12`).

Evolves quantum superposition states across explicit `2^m` measurement-reset Kraus
branches using 100% Python standard-library `decimal.Decimal` linear algebra (zero
`float` or `complex` literals, zero shared code with `hmr.py` or `kickmix_sim.py`).

Computes exact mixed-state fidelity and exact trace distance `D_tr(rho, sigma) = 1/2 ||rho - sigma||_1`
via subspace Gram-Schmidt + algebraic Jacobi diagonalization.
"""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Context, Decimal, ROUND_CEILING, ROUND_HALF_EVEN
from fractions import Fraction
from typing import Any, Callable

from .ir_kickmix import ANNOTATION_OPS, KickMixCircuit, Qubit

Cplx = tuple[Decimal, Decimal]
ZERO_D = Decimal(0)
ONE_D = Decimal(1)
HALF_D = Decimal("0.5")
EPS_PRUNE = Decimal("1e-36")
CTX = Context(prec=44, rounding=ROUND_HALF_EVEN)
CTX_CEIL = Context(prec=24, rounding=ROUND_CEILING)

# High-precision pi for Z_POW Taylor series (50 digits)
_PI_DEC = Decimal("3.14159265358979323846264338327950288419716939937510")


def _c_add(a: Cplx, b: Cplx) -> Cplx:
    return (CTX.add(a[0], b[0]), CTX.add(a[1], b[1]))


def _c_sub(a: Cplx, b: Cplx) -> Cplx:
    return (CTX.subtract(a[0], b[0]), CTX.subtract(a[1], b[1]))


def _c_neg(a: Cplx) -> Cplx:
    return (a[0].copy_negate(), a[1].copy_negate())


def _c_mul(a: Cplx, b: Cplx) -> Cplx:
    re = CTX.subtract(CTX.multiply(a[0], b[0]), CTX.multiply(a[1], b[1]))
    im = CTX.add(CTX.multiply(a[0], b[1]), CTX.multiply(a[1], b[0]))
    return (re, im)


def _c_conj(a: Cplx) -> Cplx:
    return (a[0], a[1].copy_negate())


def _c_scale(s: Decimal, a: Cplx) -> Cplx:
    return (CTX.multiply(s, a[0]), CTX.multiply(s, a[1]))


def _c_abs2(a: Cplx) -> Decimal:
    return CTX.add(CTX.multiply(a[0], a[0]), CTX.multiply(a[1], a[1]))


def _exp_i_pi_turns(half_turns: Fraction) -> Cplx:
    """Exact/44-digit `(cos(pi * half_turns), sin(pi * half_turns))` in pure `Decimal`."""
    red = half_turns % Fraction(2, 1)
    if red == 0:
        return (ONE_D, ZERO_D)
    if red == 1:
        return (-ONE_D, ZERO_D)
    if red == Fraction(1, 2):
        return (ZERO_D, ONE_D)
    if red == Fraction(3, 2):
        return (ZERO_D, -ONE_D)

    # Reduce into [-1, 1] half-turns so |theta| <= pi
    if red > Fraction(1, 1):
        red -= Fraction(2, 1)
    theta = CTX.divide(
        CTX.multiply(_PI_DEC, Decimal(red.numerator)),
        Decimal(red.denominator),
    )
    theta2 = CTX.multiply(theta, theta)

    # Taylor series for cos(theta) and sin(theta)
    cos_sum = ONE_D
    cos_term = ONE_D
    sin_sum = theta
    sin_term = theta
    for k in range(1, 32):
        cos_term = CTX.divide(CTX.multiply(cos_term, theta2), Decimal((2 * k - 1) * (2 * k))).copy_negate()
        sin_term = CTX.divide(CTX.multiply(sin_term, theta2), Decimal((2 * k) * (2 * k + 1))).copy_negate()
        cos_sum = CTX.add(cos_sum, cos_term)
        sin_sum = CTX.add(sin_sum, sin_term)
        if cos_term.copy_abs() < Decimal("1e-42") and sin_term.copy_abs() < Decimal("1e-42"):
            break
    return (cos_sum, sin_sum)


def _vec_norm2(vec: list[Cplx]) -> Decimal:
    acc = ZERO_D
    for c in vec:
        if c[0] != ZERO_D or c[1] != ZERO_D:
            acc = CTX.add(acc, _c_abs2(c))
    return acc


def _vec_inner(bra: list[Cplx], ket: list[Cplx]) -> Cplx:
    """Compute `<bra | ket>`."""
    acc = (ZERO_D, ZERO_D)
    for a, b in zip(bra, ket):
        if (a[0] != ZERO_D or a[1] != ZERO_D) and (b[0] != ZERO_D or b[1] != ZERO_D):
            acc = _c_add(acc, _c_mul(_c_conj(a), b))
    return acc


@dataclass
class _Branch:
    inv_pow2: int  # density weight multiplier = 2 ** (-inv_pow2)
    c_bits: list[int]
    cond_stack: list[int]
    vec: list[Cplx]


def evolve_kraus_ensemble(
    circuit: KickMixCircuit,
    init_vec: list[Cplx],
    *,
    max_qubits: int = 12,
    max_branches: int = 1024,
) -> list[tuple[Decimal, list[Cplx]]]:
    """Evolve pure initial state `init_vec` (`||init_vec|| == 1`) through `circuit`.

    Returns `[(weight_b, unit_or_unnormalized_vec_b), ...]` such that the output
    quantum density matrix (after tracing out classical bits) is:
        `rho_out = sum_b weight_b * |vec_b><vec_b|`.
    """
    nq = circuit.num_qubits
    nb = circuit.num_bits
    if nq > max_qubits:
        raise ValueError(f"kraus_oracle supports up to {max_qubits} qubits, got {nq}")
    dim = 1 << nq
    if len(init_vec) != dim:
        raise ValueError(f"init_vec length {len(init_vec)} != 2^{nq} ({dim})")

    branches: list[_Branch] = [
        _Branch(inv_pow2=0, c_bits=[0] * nb, cond_stack=[], vec=list(init_vec))
    ]

    for op in circuit.ops:
        kind = op.kind
        if kind in ANNOTATION_OPS:
            continue

        next_branches: list[_Branch] = []
        for br in branches:
            if kind == "PUSH_CONDITION":
                assert op.condition is not None
                top = br.cond_stack[-1] if br.cond_stack else 1
                br.cond_stack.append(top & br.c_bits[op.condition])
                next_branches.append(br)
                continue
            if kind == "POP_CONDITION":
                if br.cond_stack:
                    br.cond_stack.pop()
                next_branches.append(br)
                continue

            active = br.cond_stack[-1] if br.cond_stack else 1
            if op.condition is not None:
                active &= br.c_bits[op.condition]
            if not active:
                next_branches.append(br)
                continue

            v = br.vec
            if kind == "X":
                bit = 1 << op.qubits[0]
                for idx in range(dim):
                    if (idx & bit) == 0:
                        j = idx | bit
                        v[idx], v[j] = v[j], v[idx]
                next_branches.append(br)
            elif kind == "CX":
                c_bit = 1 << op.qubits[0]
                t_bit = 1 << op.qubits[1]
                for idx in range(dim):
                    if (idx & c_bit) and (idx & t_bit) == 0:
                        j = idx | t_bit
                        v[idx], v[j] = v[j], v[idx]
                next_branches.append(br)
            elif kind == "CCX":
                c1_bit = 1 << op.qubits[0]
                c2_bit = 1 << op.qubits[1]
                t_bit = 1 << op.qubits[2]
                cmask = c1_bit | c2_bit
                for idx in range(dim):
                    if (idx & cmask) == cmask and (idx & t_bit) == 0:
                        j = idx | t_bit
                        v[idx], v[j] = v[j], v[idx]
                next_branches.append(br)
            elif kind == "SWAP":
                b1 = 1 << op.qubits[0]
                b2 = 1 << op.qubits[1]
                for idx in range(dim):
                    if (idx & b1) and not (idx & b2):
                        j = (idx ^ b1) | b2
                        v[idx], v[j] = v[j], v[idx]
                next_branches.append(br)
            elif kind == "CSWAP":
                cb = 1 << op.qubits[0]
                b1 = 1 << op.qubits[1]
                b2 = 1 << op.qubits[2]
                for idx in range(dim):
                    if (idx & cb) and (idx & b1) and not (idx & b2):
                        j = (idx ^ b1) | b2
                        v[idx], v[j] = v[j], v[idx]
                next_branches.append(br)
            elif kind == "NEG":
                for idx in range(dim):
                    v[idx] = _c_neg(v[idx])
                next_branches.append(br)
            elif kind == "Z":
                bit = 1 << op.qubits[0]
                for idx in range(dim):
                    if idx & bit:
                        v[idx] = _c_neg(v[idx])
                next_branches.append(br)
            elif kind == "CZ":
                cmask = (1 << op.qubits[0]) | (1 << op.qubits[1])
                for idx in range(dim):
                    if (idx & cmask) == cmask:
                        v[idx] = _c_neg(v[idx])
                next_branches.append(br)
            elif kind == "CCZ":
                cmask = (1 << op.qubits[0]) | (1 << op.qubits[1]) | (1 << op.qubits[2])
                for idx in range(dim):
                    if (idx & cmask) == cmask:
                        v[idx] = _c_neg(v[idx])
                next_branches.append(br)
            elif kind == "Z_POW":
                assert op.angle is not None
                bit = 1 << op.qubits[0]
                ph = _exp_i_pi_turns(op.angle)
                for idx in range(dim):
                    if idx & bit:
                        v[idx] = _c_mul(v[idx], ph)
                next_branches.append(br)
            elif kind == "BIT_STORE0":
                br.c_bits[op.bits[0]] = 0
                next_branches.append(br)
            elif kind == "BIT_STORE1":
                br.c_bits[op.bits[0]] = 1
                next_branches.append(br)
            elif kind == "BIT_INVERT":
                br.c_bits[op.bits[0]] ^= 1
                next_branches.append(br)
            elif kind == "R":
                # Kraus operators R_0 = |0><0| and R_1 = |0><1|
                qbit = 1 << op.qubits[0]
                v0: list[Cplx] = [(ZERO_D, ZERO_D)] * dim
                v1: list[Cplx] = [(ZERO_D, ZERO_D)] * dim
                nz0 = False
                nz1 = False
                for idx in range(dim):
                    if (idx & qbit) == 0:
                        a0 = v[idx]
                        a1 = v[idx | qbit]
                        v0[idx] = a0
                        v1[idx] = a1
                        if a0[0] != ZERO_D or a0[1] != ZERO_D:
                            nz0 = True
                        if a1[0] != ZERO_D or a1[1] != ZERO_D:
                            nz1 = True
                if nz0:
                    next_branches.append(
                        _Branch(br.inv_pow2, list(br.c_bits), list(br.cond_stack), v0)
                    )
                if nz1:
                    next_branches.append(
                        _Branch(br.inv_pow2, list(br.c_bits), list(br.cond_stack), v1)
                    )
            elif kind == "HMR":
                # Kraus operators K_0 = (|0><0| + |0><1|)/sqrt(2), K_1 = (|0><0| - |0><1|)/sqrt(2)
                qbit = 1 << op.qubits[0]
                b_idx = op.bits[0]
                v0 = [(ZERO_D, ZERO_D)] * dim
                v1 = [(ZERO_D, ZERO_D)] * dim
                nz0 = False
                nz1 = False
                for idx in range(dim):
                    if (idx & qbit) == 0:
                        a0 = v[idx]
                        a1 = v[idx | qbit]
                        s0 = _c_add(a0, a1)
                        s1 = _c_sub(a0, a1)
                        v0[idx] = s0
                        v1[idx] = s1
                        if abs(s0[0]) > EPS_PRUNE or abs(s0[1]) > EPS_PRUNE:
                            nz0 = True
                        if abs(s1[0]) > EPS_PRUNE or abs(s1[1]) > EPS_PRUNE:
                            nz1 = True
                if nz0:
                    cb0 = list(br.c_bits)
                    cb0[b_idx] = 0
                    next_branches.append(
                        _Branch(br.inv_pow2 + 1, cb0, list(br.cond_stack), v0)
                    )
                if nz1:
                    cb1 = list(br.c_bits)
                    cb1[b_idx] = 1
                    next_branches.append(
                        _Branch(br.inv_pow2 + 1, cb1, list(br.cond_stack), v1)
                    )
            else:
                next_branches.append(br)

        if len(next_branches) > max_branches:
            raise ValueError(f"Kraus branch count {len(next_branches)} exceeded cap {max_branches}")
        branches = next_branches

    ensemble: list[tuple[Decimal, list[Cplx]]] = []
    for br in branches:
        wt = CTX.power(HALF_D, Decimal(br.inv_pow2)) if br.inv_pow2 > 0 else ONE_D
        ensemble.append((wt, br.vec))
    return ensemble


def _jacobi_symmetric_eigenvalues(mat: list[list[Decimal]], sweeps: int = 40) -> list[Decimal]:
    """Algebraic Jacobi diagonalization of a small real symmetric `Decimal` matrix."""
    m = len(mat)
    if m == 0:
        return []
    a = [list(row) for row in mat]
    tol = Decimal("1e-36")
    for _ in range(sweeps):
        # Find largest off-diagonal entry
        max_val = ZERO_D
        p_idx, q_idx = 0, 1
        for i in range(m):
            for j in range(i + 1, m):
                v = abs(a[i][j])
                if v > max_val:
                    max_val = v
                    p_idx, q_idx = i, j
        if max_val <= tol:
            break
        app = a[p_idx][p_idx]
        aqq = a[q_idx][q_idx]
        apq = a[p_idx][q_idx]
        diff = CTX.subtract(aqq, app)
        if abs(diff) <= Decimal("1e-40"):
            t = ONE_D if apq >= ZERO_D else -ONE_D
        else:
            tau = CTX.divide(diff, CTX.multiply(Decimal(2), apq))
            root = CTX.sqrt(CTX.add(ONE_D, CTX.multiply(tau, tau)))
            denom = CTX.add(abs(tau), root)
            t = CTX.divide(ONE_D if tau >= ZERO_D else -ONE_D, denom)
        c = CTX.divide(ONE_D, CTX.sqrt(CTX.add(ONE_D, CTX.multiply(t, t))))
        s = CTX.multiply(t, c)

        # Update matrix entries
        a[p_idx][p_idx] = CTX.subtract(app, CTX.multiply(t, apq))
        a[q_idx][q_idx] = CTX.add(aqq, CTX.multiply(t, apq))
        a[p_idx][q_idx] = ZERO_D
        a[q_idx][p_idx] = ZERO_D
        for r in range(m):
            if r != p_idx and r != q_idx:
                arp = a[r][p_idx]
                arq = a[r][q_idx]
                new_rp = CTX.subtract(CTX.multiply(c, arp), CTX.multiply(s, arq))
                new_rq = CTX.add(CTX.multiply(s, arp), CTX.multiply(c, arq))
                a[r][p_idx] = new_rp
                a[p_idx][r] = new_rp
                a[r][q_idx] = new_rq
                a[q_idx][r] = new_rq
    return [a[i][i] for i in range(m)]


def exact_ensemble_trace_distance(
    ens_a: list[tuple[Decimal, list[Cplx]]],
    ens_b: list[tuple[Decimal, list[Cplx]]],
) -> Decimal:
    """Exact trace distance `1/2 ||rho_A - rho_B||_1` in zero-float `Decimal`."""
    # Orthonormalize active span via Modified Gram-Schmidt
    basis: list[list[Cplx]] = []
    for wt, vec in ens_a + ens_b:
        if wt <= ZERO_D:
            continue
        cur = list(vec)
        for bvec in basis:
            proj = _vec_inner(bvec, cur)
            if abs(proj[0]) > EPS_PRUNE or abs(proj[1]) > EPS_PRUNE:
                for idx in range(len(cur)):
                    cur[idx] = _c_sub(cur[idx], _c_mul(proj, bvec[idx]))
        n2 = _vec_norm2(cur)
        if n2 > Decimal("1e-32"):
            inv_n = CTX.divide(ONE_D, CTX.sqrt(n2))
            basis.append([_c_scale(inv_n, z) for z in cur])

    d = len(basis)
    if d == 0:
        return ZERO_D

    # Build Hermitian difference Delta = rho_A - rho_B in basis: d x d complex
    re_h = [[ZERO_D] * d for _ in range(d)]
    im_h = [[ZERO_D] * d for _ in range(d)]

    def accumulate(ens: list[tuple[Decimal, list[Cplx]]], sign: Decimal) -> None:
        for wt, vec in ens:
            if wt <= ZERO_D:
                continue
            s_wt = CTX.multiply(sign, wt)
            coeffs = [_vec_inner(basis[p], vec) for p in range(d)]
            for p in range(d):
                for q in range(d):
                    prod = _c_mul(coeffs[p], _c_conj(coeffs[q]))
                    re_h[p][q] = CTX.add(re_h[p][q], CTX.multiply(s_wt, prod[0]))
                    im_h[p][q] = CTX.add(im_h[p][q], CTX.multiply(s_wt, prod[1]))

    accumulate(ens_a, ONE_D)
    accumulate(ens_b, -ONE_D)

    # Real symmetric embedding S = [[R, -M], [M, R]] of size 2d x 2d
    sym = [[ZERO_D] * (2 * d) for _ in range(2 * d)]
    for p in range(d):
        for q in range(d):
            r_pq = CTX.multiply(HALF_D, CTX.add(re_h[p][q], re_h[q][p]))
            m_pq = CTX.multiply(HALF_D, CTX.subtract(im_h[p][q], im_h[q][p]))
            sym[p][q] = r_pq
            sym[p + d][q + d] = r_pq
            sym[p][q + d] = m_pq.copy_negate()
            sym[p + d][q] = m_pq

    eigs = _jacobi_symmetric_eigenvalues(sym)
    sum_abs = ZERO_D
    for ev in eigs:
        sum_abs = CTX.add(sum_abs, ev.copy_abs())
    # D_tr = 1/4 * sum |lambda(S)|
    d_tr = CTX.multiply(Decimal("0.25"), sum_abs)
    if d_tr < Decimal("1e-30"):
        return ZERO_D
    return +CTX_CEIL.create_decimal(d_tr)


@dataclass(frozen=True)
class KrausVerificationResult:
    exact_match: bool
    trace_distance: Decimal
    min_fidelity: Decimal
    trace_preservation_error: Decimal
    branch_count_max: int
    num_qubits: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "exact_match": self.exact_match,
            "trace_distance": str(self.trace_distance),
            "min_fidelity": str(self.min_fidelity),
            "trace_preservation_error": str(self.trace_preservation_error),
            "branch_count_max": self.branch_count_max,
            "num_qubits": self.num_qubits,
        }


OracleSpec = Callable[[dict[str, int]], tuple[dict[str, int], Fraction | int] | dict[str, int]]


def verify_channel_kraus(
    circuit: KickMixCircuit,
    *,
    oracle: OracleSpec | None = None,
    max_qubits: int = 12,
    tolerance: Decimal = Decimal("1e-24"),
) -> KrausVerificationResult:
    """Verify `circuit` as a quantum channel on full coherent superpositions (`n <= 12`).

    Tests both the uniform superposition `|+>^{n_pub} |0>^{n_scratch}` and a
    basis-distinguishing complex superposition `sum_x ((x + 1) + i) |x> |0>`,
    computing exact mixed-state fidelity and trace distance against `oracle`.
    """
    nq = circuit.num_qubits
    dim = 1 << nq
    q_reg_specs = [
        (r.name if r.name is not None else f"r{r.reg_id}", tuple(w.id for w in r.wires if isinstance(w, Qubit)))
        for r in circuit.registers
        if r.kind == "q"
    ]
    if not circuit.registers and nq > 0:
        q_reg_specs = [("q_all", tuple(range(nq)))]

    pub_qids = tuple(qid for _, qids in q_reg_specs for qid in qids)
    n_pub = len(pub_qids)
    d_pub = 1 << n_pub

    def pack_reg_dict(x_pub: int) -> tuple[int, dict[str, int]]:
        full_idx = 0
        reg_vals: dict[str, int] = {}
        bit_cursor = 0
        for name, qids in q_reg_specs:
            val = 0
            for j, qid in enumerate(qids):
                bval = (x_pub >> (bit_cursor + j)) & 1
                if bval:
                    val |= 1 << j
                    full_idx |= 1 << qid
            reg_vals[name] = val
            bit_cursor += len(qids)
        return full_idx, reg_vals

    def pack_out_idx(out_regs: dict[str, int]) -> int:
        full_idx = 0
        for name, qids in q_reg_specs:
            val = out_regs[name]
            for j, qid in enumerate(qids):
                if (val >> j) & 1:
                    full_idx |= 1 << qid
        return full_idx

    # Build two probe pairs (|psi_in>, |phi_target>): uniform + basis-distinguishing complex
    probes: list[tuple[list[Cplx], list[Cplx]]] = []
    for mode in ("uniform", "distinguishing"):
        psi_in: list[Cplx] = [(ZERO_D, ZERO_D)] * dim
        phi_tgt: list[Cplx] = [(ZERO_D, ZERO_D)] * dim
        norm2 = ZERO_D
        for x_pub in range(d_pub):
            amp: Cplx = (ONE_D, ZERO_D) if mode == "uniform" else (Decimal(x_pub + 1), ONE_D)
            norm2 = CTX.add(norm2, _c_abs2(amp))
            in_idx, in_regs = pack_reg_dict(x_pub)
            psi_in[in_idx] = amp

            if oracle is None:
                out_regs = in_regs
                ph_turns = Fraction(0, 1)
            else:
                res = oracle(in_regs)
                if isinstance(res, tuple):
                    out_regs, raw_ph = res
                    ph_turns = Fraction(raw_ph) % Fraction(2, 1)
                else:
                    out_regs = res
                    ph_turns = Fraction(0, 1)
            out_idx = pack_out_idx(out_regs)
            ph_cplx = _exp_i_pi_turns(ph_turns)
            phi_tgt[out_idx] = _c_add(phi_tgt[out_idx], _c_mul(amp, ph_cplx))

        inv_norm = CTX.divide(ONE_D, CTX.sqrt(norm2))
        psi_in = [_c_scale(inv_norm, z) for z in psi_in]
        phi_tgt = [_c_scale(inv_norm, z) for z in phi_tgt]
        probes.append((psi_in, phi_tgt))

    max_d_tr = ZERO_D
    min_fid = ONE_D
    max_tp_err = ZERO_D
    max_branches = 0

    for psi_in, phi_tgt in probes:
        ens = evolve_kraus_ensemble(circuit, psi_in, max_qubits=max_qubits)
        max_branches = max(max_branches, len(ens))

        # Check trace preservation & fidelity <phi_tgt | rho | phi_tgt>
        tr_sum = ZERO_D
        fid = ZERO_D
        for wt, bvec in ens:
            tr_sum = CTX.add(tr_sum, CTX.multiply(wt, _vec_norm2(bvec)))
            ov = _vec_inner(phi_tgt, bvec)
            fid = CTX.add(fid, CTX.multiply(wt, _c_abs2(ov)))

        tp_err = CTX.subtract(ONE_D, tr_sum).copy_abs()
        if tp_err > max_tp_err:
            max_tp_err = tp_err
        if fid < min_fid:
            min_fid = fid

        # Fast-path when fidelity == 1 within machine precision (1e-36)
        if CTX.subtract(ONE_D, fid).copy_abs() <= Decimal("1e-36") and tp_err <= Decimal("1e-36"):
            d_tr = ZERO_D
        else:
            d_tr = exact_ensemble_trace_distance(ens, [(ONE_D, phi_tgt)])
        if d_tr > max_d_tr:
            max_d_tr = d_tr

    clean_fid = ONE_D if CTX.subtract(ONE_D, min_fid).copy_abs() <= Decimal("1e-36") else +CTX.create_decimal(min_fid)
    clean_tp = ZERO_D if max_tp_err <= Decimal("1e-36") else +CTX_CEIL.create_decimal(max_tp_err)
    exact = max_d_tr <= tolerance and clean_tp <= tolerance

    return KrausVerificationResult(
        exact_match=exact,
        trace_distance=max_d_tr,
        min_fidelity=clean_fid,
        trace_preservation_error=clean_tp,
        branch_count_max=max_branches,
        num_qubits=nq,
    )
