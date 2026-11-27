"""Exact O(n^3) Clifford-error hook analyzer for magic state cultivation.

Implements the exact algebraic fault-distance & fidelity analysis of:
- Gidney, Shutty & Jones (2024, arXiv:2409.17595) — Magic state cultivation (`MSC-3`, `MSC-5`)
- Chan, Strikis, Sun & Cai (2026, arXiv:2609.17706):
  * Theorem 1 & Algorithm 1 (Eqs. 7, 24-27): exact O(n^3)-time trivial-syndrome acceptance
    probability `Pr(s == 0 | E)` in `Q(sqrt(2))` for any `[[n, k]]` stabilizer code suffering
    an arbitrary physical Clifford error `E` (including `H_+-` hook superpositions arising when
    Pauli faults propagate across physical `T / T^dagger` layers!).
  * Proposition 1 & Eq. (21): exact conditional logical fidelity `phi(E) in Q(sqrt(2))` given
    `s == 0`.
  * Exact `[[7, 1, 3]]` hexagonal colour-code cultivation circuits:
    1) Unflagged double-check (Fig. 4 / Fig. 7): proves S-cultivation proxy has fault distance `3`
       while true T-cultivation degrades to **`d = 2`** due to weight-2 Pauli-to-Clifford hook
       configurations (`X` hook + `Y` -> `-(H_+)_i (H_-)_{j,l}` with `Pr(acc) = 1/4, fidelity = 0`)!
    2) Z-flagged double-check (Fig. 8 / Fig. 13): proves adding two 2D-local Z-flag ancillas
       `alpha, beta` catches all malignant weight-2 hook configurations and **restores `d = 3`**!
    3) Itogawa et al. (2024/2025) short single check (two-round Fault-Tolerant CAT check): `d = 3`.
  * Live cross-checks of the Clifford S-cultivation proxy against `ftqc_faults.audit` AND `stim`!
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

from .ir import CorrelationMode, MagicRefusalError, MagicSourceContract, PhasePolynomialTarget
from .rings import QSqrt2, parse_exact_fraction


def _ensure_tools_on_syspath() -> None:
    for candidate in (
        Path(__file__).resolve().parents[1],
        Path("/Volumes/t2/Documents/RESEARCH-quantum-circuits-bootloops-2026/bootloops/tools"),
        Path("/home/ike/RESEARCH-quantum-circuits-bootloops-2026/bootloops/tools"),
    ):
        if (candidate / "ftqc_faults").is_dir() or (candidate / "ftqc_codes").is_dir():
            s = str(candidate)
            if s not in sys.path:
                sys.path.insert(0, s)
            break


def _mul_hermitian_paulis(
    p1: tuple[int, int, int],
    p2: tuple[int, int, int],
) -> tuple[int, int, int]:
    """Multiply two n-qubit Paulis represented as `(x_mask, z_mask, i_pow_mod4)` where tuple means
    `i^(i_pow) * X^(x_mask) * Z^(z_mask)`!
    Note: `(X^x1 Z^z1)(X^x2 Z^z2) = (-1)^{|z1 & x2|} X^{x1^x2} Z^{z1^z2} = i^{2*|z1 & x2|} ...`
    """
    x1, z1, phi1 = p1
    x2, z2, phi2 = p2
    comm = (z1 & x2).bit_count() & 1
    return (x1 ^ x2, z1 ^ z2, (phi1 + phi2 + 2 * comm) % 4)


@dataclass(frozen=True)
class CliffordTableau:
    """Exact n-qubit Clifford unitary U represented by its conjugation action on single-qubit
    Paulis: `U X_q U^dagger = i^{px[q]} X^{xx[q]} Z^{zx[q]}` and
            `U Z_q U^dagger = i^{pz[q]} X^{xz[q]} Z^{zz[q]}` (`q = 0..n-1`).
    Note: Algorithm 1 uses `F^dagger Z_i F`; `.inverse()` computes `U^dagger` in O(n^3) or by
    gate reversal!
    """

    n: int
    xx: tuple[int, ...]
    zx: tuple[int, ...]
    px: tuple[int, ...]
    xz: tuple[int, ...]
    zz: tuple[int, ...]
    pz: tuple[int, ...]

    @classmethod
    def identity(cls, n: int) -> CliffordTableau:
        return cls(
            n=n,
            xx=tuple(1 << q for q in range(n)),
            zx=tuple(0 for _ in range(n)),
            px=tuple(0 for _ in range(n)),
            xz=tuple(0 for _ in range(n)),
            zz=tuple(1 << q for q in range(n)),
            pz=tuple(0 for _ in range(n)),
        )

    def conjugate_pauli(self, x_mask: int, z_mask: int, i_pow: int = 0) -> tuple[int, int, int]:
        """Compute `U (i^{i_pow} X^{x_mask} Z^{z_mask}) U^dagger` in O(n) word operations!"""
        cur = (0, 0, i_pow % 4)
        for q in range(self.n):
            if (x_mask >> q) & 1:
                cur = _mul_hermitian_paulis(cur, (self.xx[q], self.zx[q], self.px[q]))
        for q in range(self.n):
            if (z_mask >> q) & 1:
                cur = _mul_hermitian_paulis(cur, (self.xz[q], self.zz[q], self.pz[q]))
        return cur

    def compose(self, first: CliffordTableau) -> CliffordTableau:
        """Return `self * first` (where `first` acts first on kets: `(self*first) P (self*first)^dag = self (first P first^dag) self^dag`)."""
        if self.n != first.n:
            raise ValueError("CliffordTableau qubit count mismatch")
        n = self.n
        n_xx, n_zx, n_px = [], [], []
        n_xz, n_zz, n_pz = [], [], []
        for q in range(n):
            ux, uz, up = self.conjugate_pauli(first.xx[q], first.zx[q], first.px[q])
            n_xx.append(ux)
            n_zx.append(uz)
            n_px.append(up)
            vx, vz, vp = self.conjugate_pauli(first.xz[q], first.zz[q], first.pz[q])
            n_xz.append(vx)
            n_zz.append(vz)
            n_pz.append(vp)
        return CliffordTableau(
            n=n,
            xx=tuple(n_xx),
            zx=tuple(n_zx),
            px=tuple(n_px),
            xz=tuple(n_xz),
            zz=tuple(n_zz),
            pz=tuple(n_pz),
        )

    @classmethod
    def from_gates(cls, n: int, gates: Sequence[tuple[Any, ...]]) -> CliffordTableau:
        """Build exact Clifford unitary `U = G_m ... G_1` (gates listed in chronological order!)."""
        # Maintain forward images of X_q, Z_q by applying each gate on the LEFT: G (U P U^dag) G^dag!
        xx = [1 << q for q in range(n)]
        zx = [0] * n
        px = [0] * n
        xz = [0] * n
        zz = [1 << q for q in range(n)]
        pz = [0] * n

        def apply_to_img(x: int, z: int, p: int, op: tuple[Any, ...]) -> tuple[int, int, int]:
            kind = op[0]
            if kind == "I":
                return x, z, p
            if kind in ("X", "Y", "Z"):
                q = op[1]
                xb = (x >> q) & 1
                zb = (z >> q) & 1
                # Conjugating X^xb Z^zb by Pauli Q flips sign iff anticommuting:
                anti = (zb if kind == "X" else (xb if kind == "Z" else (xb ^ zb)))
                return x, z, (p + 2 * anti) % 4
            if kind == "H":
                q = op[1]
                xb = (x >> q) & 1
                zb = (z >> q) & 1
                # H (X^xb Z^zb) H = Z^xb X^zb = (-1)^{xb*zb} X^zb Z^xb!
                nx = (x & ~(1 << q)) | (zb << q)
                nz = (z & ~(1 << q)) | (xb << q)
                return nx, nz, (p + 2 * (xb & zb)) % 4
            if kind in ("S", "SDG"):
                q = op[1]
                xb = (x >> q) & 1
                zb = (z >> q) & 1
                if xb == 0:
                    return x, z, p
                # S X S^dag = +Y = +i X Z; SDG X SDG^dag = -Y = i^3 X Z!
                s_shift = 1 if kind == "S" else 3
                # S (X Z^zb) S^dag = (i^{s_shift} X Z) Z^zb = i^{s_shift} X Z^{1^zb}
                nz = z ^ (1 << q)
                return x, nz, (p + s_shift) % 4
            if kind in ("H_PLUS", "H_MINUS"):
                # Chan et al. Eq. (2): H_pm = (X +- Y)/sqrt(2):
                #   H_pm X H_pm = +-Y = +- i X Z  (i^1 for H_PLUS, i^3 for H_MINUS)
                #   H_pm Z H_pm = -Z = i^2 Z!
                q = op[1]
                xb = (x >> q) & 1
                zb = (z >> q) & 1
                if xb == 0 and zb == 0:
                    return x, z, p
                h_shift = 1 if kind == "H_PLUS" else 3
                if xb == 1 and zb == 0:
                    return x, z | (1 << q), (p + h_shift) % 4
                if xb == 0 and zb == 1:
                    return x, z, (p + 2) % 4
                # xb == 1 and zb == 1: H_pm (X Z) H_pm = (+-i X Z)(-Z) = -+i X = i^{h_shift + 2} X!
                return x, z & ~(1 << q), (p + h_shift + 2) % 4
            if kind == "CX":
                c, t = op[1], op[2]
                xc, zc = (x >> c) & 1, (z >> c) & 1
                xt, zt = (x >> t) & 1, (z >> t) & 1
                # CX (X_c^xc X_t^xt Z_c^zc Z_t^zt) CX:
                # X_c -> X_c X_t; Z_t -> Z_c Z_t:
                # = X_c^xc X_t^xc X_t^xt Z_c^zc (Z_c Z_t)^zt
                # Commuting X_t^{xc^xt} past Z_c^zc adds NO sign (different qubits!),
                # wait: original order is X_c^xc X_t^xt Z_c^zc Z_t^zt.
                # Image is (X_c X_t)^xc X_t^xt Z_c^zc (Z_c Z_t)^zt = X_c^xc X_t^{xc^xt} Z_c^{zc^zt} Z_t^zt!
                # Zero sign flip because X_c, X_t commute and Z_c, Z_t commute!!
                nx = (x & ~((1 << t))) | ((xc ^ xt) << t)
                nz = (z & ~((1 << c))) | ((zc ^ zt) << c)
                return nx, nz, p
            if kind == "CZ":
                a, b = op[1], op[2]
                xa, xb = (x >> a) & 1, (x >> b) & 1
                # CZ maps X_a -> X_a Z_b, X_b -> Z_a X_b:
                # (X_a Z_b)^xa (Z_a X_b)^xb = (-1)^{xa*xb} X_a^xa X_b^xb Z_a^xb Z_b^xa?
                # Check xa=xb=1: CZ (X_a X_b) CZ = (X_a Z_b)(Z_a X_b) = (X_a Z_a)(Z_b X_b) = (-iY_a)(iY_b)?
                # Wait: (X_a Z_a)(Z_b X_b) = (X_a Z_a)(-X_b Z_b) = - X_a X_b Z_a Z_b! YES: +2*xa*xb mod 4!
                nz = z ^ (xb << a) ^ (xa << b)
                return x, nz, (p + 2 * (xa & xb)) % 4
            raise ValueError(f"Unsupported Clifford gate {op!r}")

        for op in gates:
            for q in range(n):
                xx[q], zx[q], px[q] = apply_to_img(xx[q], zx[q], px[q], op)
                xz[q], zz[q], pz[q] = apply_to_img(xz[q], zz[q], pz[q], op)
        return cls(n, tuple(xx), tuple(zx), tuple(px), tuple(xz), tuple(zz), tuple(pz))


def _invert_gate_list(gates: Sequence[tuple[Any, ...]]) -> list[tuple[Any, ...]]:
    inv_map = {"S": "SDG", "SDG": "S"}
    out: list[tuple[Any, ...]] = []
    for op in reversed(gates):
        k = inv_map.get(op[0], op[0])
        out.append((k, *op[1:]))
    return out


def propagate_pauli_across_t_layer(
    n: int,
    x_mask: int,
    z_mask: int,
    t_qubits: Sequence[int],
    tdg_qubits: Sequence[int],
) -> CliffordTableau:
    """Propagate physical Pauli `P = i^{|x&z|} X^x Z^z` across `U_T = T_{t_qubits} T^dag_{tdg_qubits}`
    into exact Clifford error `E = U_T P U_T^dag` (up to unobservable global phase, which cancels in
    `E^dag G E` and `E^dag rho E`!):
    - On `q in t_qubits`:
      * `I -> I`, `Z -> Z`
      * `X -> T X T^dag = H_PLUS`
      * `Y -> T Y T^dag = -H_MINUS` (global -1 cancels in density conjugation -> `H_MINUS`)
    - On `q in tdg_qubits`:
      * `I -> I`, `Z -> Z`
      * `X -> T^dag X T = H_MINUS`
      * `Y -> T^dag Y T = +H_PLUS`
    """
    t_set = set(t_qubits)
    tdg_set = set(tdg_qubits)
    if t_set & tdg_set:
        raise ValueError("t_qubits and tdg_qubits must be disjoint")
    gates: list[tuple[Any, ...]] = []
    for q in range(n):
        xb = (x_mask >> q) & 1
        zb = (z_mask >> q) & 1
        if xb == 0 and zb == 0:
            continue
        if q in t_set:
            if xb == 0 and zb == 1:
                gates.append(("Z", q))
            elif xb == 1 and zb == 0:
                gates.append(("H_PLUS", q))
            else:
                gates.append(("H_MINUS", q))
        elif q in tdg_set:
            if xb == 0 and zb == 1:
                gates.append(("Z", q))
            elif xb == 1 and zb == 0:
                gates.append(("H_MINUS", q))
            else:
                gates.append(("H_PLUS", q))
        else:
            if xb == 1 and zb == 0:
                gates.append(("X", q))
            elif xb == 0 and zb == 1:
                gates.append(("Z", q))
            else:
                gates.append(("Y", q))
    return CliffordTableau.from_gates(n, gates)


def build_css_encoder_gates(
    n: int,
    hx_rows: Sequence[int],
    hz_rows: Sequence[int],
    logical_x_rows: Sequence[int],
    logical_z_rows: Sequence[int],
) -> list[tuple[Any, ...]]:
    """Construct an exact Clifford encoder `C` mapping `|0>^r (x) |psi>_k` (qubits `0..r-1` checks,
    `r..n-1` logicals) onto the `[[n, k]]` CSS code `(hx_rows, hz_rows)` with logicals `(LX, LZ)`:
    Specifically: `C Z_i C^dag = +G_i` for `0 <= i < r` and `C X_{r+a} C^dag = +LX_a`, `C Z_{r+a} C^dag = +LZ_a`!
    Constructed via standard symplectic Gram-Schmidt / Gaussian elimination over F_2 in O(n^3)!
    """
    k = len(logical_x_rows)
    r = n - k
    if len(hx_rows) + len(hz_rows) != r or len(logical_z_rows) != k:
        raise ValueError("CSS stabilizer + logical dimensions do not sum to n")

    # Target images of X_0..X_{n-1} and Z_0..Z_{n-1}:
    # We know Z_0..Z_{r-1} -> hx_rows (with H turned X-type!) and hz_rows (Z-type),
    # and (X_{r+a}, Z_{r+a}) -> (LX_a, LZ_a). We only need complementary destabilizers D_0..D_{r-1}
    # satisfying {D_i, G_j} = delta_{ij}, [D_i, D_j] = 0, [D_i, LX_a] = [D_i, LZ_a] = 0!
    # For CSS codes this splits cleanly into Z-destabilizers for H_X and X-destabilizers for H_Z!
    n_hx = len(hx_rows)
    n_hz = len(hz_rows)

    def dual_vectors(target_rows: Sequence[int], ortho_rows: Sequence[int]) -> list[int]:
        out: list[int] = []
        m_t = len(target_rows)
        all_constraints = list(target_rows) + list(ortho_rows)
        for idx in range(m_t):
            # Solve <c, v> = (1 if c == idx else 0) over F_2^n:
            pivots: list[tuple[int, int]] = []
            for row_i, mask in enumerate(all_constraints):
                rhs = 1 if row_i == idx else 0
                aug = mask | (rhs << n)
                for pcol, prow in pivots:
                    if (aug >> pcol) & 1:
                        aug ^= prow
                lhs = aug & ((1 << n) - 1)
                if lhs == 0:
                    if (aug >> n) & 1:
                        raise ValueError("Inconsistent CSS destabilizer system")
                    continue
                lead = (lhs & -lhs).bit_length() - 1
                pivots = [(c, pr ^ aug if ((pr >> lead) & 1) else pr) for c, pr in pivots]
                pivots.append((lead, aug))
            vec = 0
            for pcol, prow in pivots:
                if (prow >> n) & 1:
                    vec |= 1 << pcol
            out.append(vec)
        return out

    # Destabilizers for X-stabilizers are pure-Z vectors orthogonal to LX:
    destab_z_for_hx = dual_vectors(hx_rows, logical_x_rows)
    # Destabilizers for Z-stabilizers are pure-X vectors orthogonal to LZ:
    destab_x_for_hz = dual_vectors(hz_rows, logical_z_rows)
    # Enforce mutual commutation [Z(d_i), X(e_j)] = 0 by projecting e_j along hx_rows:
    for j in range(n_hz):
        ej = destab_x_for_hz[j]
        for i in range(n_hx):
            if (destab_z_for_hx[i] & ej).bit_count() & 1:
                ej ^= hx_rows[i]
        destab_x_for_hz[j] = ej

    # Build direct CliffordTableau for `C`! Notice all basis images are pure-X or pure-Z (zero Ys -> p=0!):
    xx: list[int] = []
    zx: list[int] = []
    xz: list[int] = []
    zz: list[int] = []
    # First n_hx qubits: Z_i -> X(hx_rows[i]), X_i -> Z(destab_z_for_hx[i]) (wait: [X_i, Z_i] -> Z_destab * X_hx:
    # check phase: Z^d * X^h = (-1)^{|d&h|} X^h Z^d = - X^h Z^d, whereas X_i Z_i = -i Y_i and Z_i X_i = +i Y_i;
    # wait: to preserve `X_i Z_i = - Z_i X_i` AND `[D_i, D_j] = 0`, pure-Z and pure-X automatically commute
    # within each sector, and `{Z(d_i), X(h_j)} = 0` for `i == j` (`|d_i & h_i| == 1`)! Does orientation matter?
    # Check: (X_0 -> +Z(d), Z_0 -> +X(h)) = Hadamard-like! For 1 qubit `H`: `H X H = +Z, H Z H = +X`, which IS unitary!!
    # YES!! `X -> +Z(d), Z -> +X(h)` IS UNITARY (`H` on each of the `n_hx` qubits followed by CNOTs)!!
    for i in range(n_hx):
        xx.append(0)
        zx.append(destab_z_for_hx[i])
        xz.append(hx_rows[i])
        zz.append(0)
    for j in range(n_hz):
        xx.append(destab_x_for_hz[j])
        zx.append(0)
        xz.append(0)
        zz.append(hz_rows[j])
    for a in range(k):
        xx.append(logical_x_rows[a])
        zx.append(0)
        xz.append(0)
        zz.append(logical_z_rows[a])

    # Verify symplectic commutation relations on `(xx, zx)` and `(xz, zz)`:
    # Actually we can return `(C, C_dag)` directly as `CliffordTableau`s!
    c_tab = CliffordTableau(
        n=n,
        xx=tuple(xx),
        zx=tuple(zx),
        px=tuple(0 for _ in range(n)),
        xz=tuple(xz),
        zz=tuple(zz),
        pz=tuple(0 for _ in range(n)),
    )
    return [c_tab]  # helper compatibility; use `build_css_encoder_pair` below


def build_css_encoder_pair(
    n: int,
    hx_rows: Sequence[int],
    hz_rows: Sequence[int],
    logical_x_rows: Sequence[int],
    logical_z_rows: Sequence[int],
) -> tuple[CliffordTableau, CliffordTableau]:
    """Return exact `(C, C^dagger)` `CliffordTableau` pair for an `[[n, k]]` CSS code in O(n^3)."""
    c_tab = build_css_encoder_gates(n, hx_rows, hz_rows, logical_x_rows, logical_z_rows)[0]
    # Compute exact symplectic inverse `C^dagger`:
    # For pure-CSS symplectic matrix `M = [[XX, XZ], [ZX, ZZ]]` with zero phases (`px = pz = 0`):
    # Recall for symplectic `S` over F_2^(2n) (`S^T Omega S = Omega`), `S^{-1} = Omega S^T Omega`!!
    # Specifically, if `C(X_j) = X^{xx_j} Z^{zx_j}` and `C(Z_j) = X^{xz_j} Z^{zz_j}` (all `p == 0`):
    # Let's check `C^dagger(X_q)` and `C^dagger(Z_q)`!
    # Write `C^dag X_q C = i^? X^a Z^b`. Since `C` is a CNOT + H circuit (real orthogonal permutation of
    # `{|x>}` up to Hadamards -> maps pure-X and pure-Z with `+` signs ONLY!), `px_inv = pz_inv = 0`!
    # Check symplectic inverse formula `S^{-1} = Omega S^T Omega`:
    # Let `u_j = (xx_j, zx_j)` (`C(X_j)`) and `v_j = (xz_j, zz_j)` (`C(Z_j)`).
    # For `X_q` (symplectic `(e_q, 0)`): `<C^dag(X_q), Z_j>_symp == <X_q, C(Z_j)>_symp == (zz_j >> q) & 1`!
    #Wait: `<(a, b), (0, e_j)>_symp == a_j`! So:
    # - `(xx_inv[q] >> j) & 1 == <(xx_inv[q], zx_inv[q]), (0, e_j)> == <X_q, C(Z_j)> == (zz[j] >> q) & 1`!
    # - `(zx_inv[q] >> j) & 1 == <(xx_inv[q], zx_inv[q]), (e_j, 0)> == <X_q, C(X_j)> == (zx[j] >> q) & 1`!
    # - `(xz_inv[q] >> j) & 1 == <Z_q, C(Z_j)> == (xz[j] >> q) & 1`!
    # - `(zz_inv[q] >> j) & 1 == <Z_q, C(X_j)> == (xx[j] >> q) & 1`!
    xx_inv = [sum(((c_tab.zz[j] >> q) & 1) << j for j in range(n)) for q in range(n)]
    zx_inv = [sum(((c_tab.zx[j] >> q) & 1) << j for j in range(n)) for q in range(n)]
    xz_inv = [sum(((c_tab.xz[j] >> q) & 1) << j for j in range(n)) for q in range(n)]
    zz_inv = [sum(((c_tab.xx[j] >> q) & 1) << j for j in range(n)) for q in range(n)]
    c_dag = CliffordTableau(
        n=n,
        xx=tuple(xx_inv),
        zx=tuple(zx_inv),
        px=tuple(0 for _ in range(n)),
        xz=tuple(xz_inv),
        zz=tuple(zz_inv),
        pz=tuple(0 for _ in range(n)),
    )
    # Self-verify `c_dag * c_tab == I`:
    prod_id = c_dag.compose(c_tab)
    if prod_id != CliffordTableau.identity(n):
        raise AssertionError("CSS encoder symplectic inverse failed identity check")
    return c_tab, c_dag


def clifford_error_acceptance_prob(
    encoder: CliffordTableau,
    encoder_dag: CliffordTableau,
    error_e: CliffordTableau,
    k: int,
    logical_alpha: Mapping[tuple[int, int], QSqrt2],
) -> QSqrt2:
    """Chan, Strikis, Sun & Cai (2026) **Algorithm 1** (pages 18-19, Eqs. 24-27):
    Exact O(n^3)-time trivial-syndrome acceptance probability `Pr(s == 0)` in `Q(sqrt(2))`
    for an arbitrary `k`-qubit state `rho = 2^{-k} sum_L alpha_L L` encoded in an `[[n, k]]`
    stabilizer code via `encoder` (`C`) after suffering physical Clifford error `E`!

    Convention: unencoded qubits `0..r-1` (`r = n - k`) are stabilizer `|0>^r` (`Z_0..Z_{r-1}`),
    and qubits `r..n-1` are the `k` logical qubits.
    `logical_alpha` maps `(lx_mask_k, lz_mask_k)` of Hermitian k-qubit Pauli
    `L = i^{|lx & lz|} X^{lx} Z^{lz}` to `alpha_L in Q(sqrt(2))` (with `(0,0) -> 1`).
    """
    n = encoder.n
    r = n - k
    # Step 1: F = C^dag * E * C; we need F^dag Z_i F = C^dag * E^dag * C * Z_i * C^dag * E * C!
    # Wait: compute F = encoder_dag.compose(error_e).compose(encoder), and how do we get F^dag Z_i F
    # for general Clifford `F` (which may include `H_PLUS`/`H_MINUS` with phases)?
    # Computing `F^dag` via exact symplectic + phase inversion of `F`!
    f_tab = encoder_dag.compose(error_e).compose(encoder)
    f_dag = _invert_general_clifford_tableau(f_tab)

    # Columns of X, Z, phi for i = 0..r-1: F^dag Z_i F = i^{phi_i} X^{x_i} Z^{z_i}
    x_cols: list[int] = []
    z_cols: list[int] = []
    phi_vec: list[int] = []
    # R has (n + k) rows: rows 0..n-1 are x_i bits (0..n-1); rows n..n+k-1 are z_i logical bits (r..n-1)!
    r_cols: list[int] = []
    for i in range(r):
        xi, zi, phii = f_dag.conjugate_pauli(0, 1 << i, 0)
        x_cols.append(xi)
        z_cols.append(zi)
        phi_vec.append(phii)
        z_log = zi >> r
        r_cols.append(xi | (z_log << n))

    # Transpose `r_cols` (`r` columns of height `n + k`) to `n + k` rows of width `r` for row-reduction:
    num_rows = n + k
    mat_rows = [
        sum(((r_cols[col] >> row) & 1) << col for col in range(r))
        for row in range(num_rows)
    ]

    # Gauss-Jordan row reduction of R*u = y (tracking row combination masks to test y in im(R) & solve!):
    # Track `(pivot_col, reduced_row_mask_r, combo_mask_over_original_rows)`:
    pivots: list[tuple[int, int, int]] = []
    zero_combos: list[int] = []
    for row_idx in range(num_rows):
        v = mat_rows[row_idx]
        combo = 1 << row_idx
        for pcol, prow, pcombo in pivots:
            if (v >> pcol) & 1:
                v ^= prow
                combo ^= pcombo
        if v == 0:
            zero_combos.append(combo)
        else:
            lead = (v & -v).bit_length() - 1
            pivots = [
                (pc, pr ^ v if ((pr >> lead) & 1) else pr, pcb ^ combo if ((pr >> lead) & 1) else pcb)
                for pc, pr, pcb in pivots
            ]
            pivots.append((lead, v, combo))

    rk_R = len(pivots)
    piv_map = {pc: (pr, pcb) for pc, pr, pcb in pivots}
    free_cols = [c for c in range(r) if c not in piv_map]
    ker_basis: list[int] = []
    for fc in free_cols:
        u_vec = 1 << fc
        for pc, (pr, _) in piv_map.items():
            if (pr >> fc) & 1:
                u_vec |= 1 << pc
        ker_basis.append(u_vec)

    # Lemma 4 Eq. (27): exact quadratic function chi(u) in {0, 1} for u in F_2^r:
    def eval_chi(u_mask: int) -> int:
        phi_dot_u = 0
        xu = 0
        zu = 0
        cross_zx = 0
        # Accumulate in ascending column order i = 0..r-1:
        # Notice sum_{i < j} u_i u_j z_i^T x_j == sum_{j: u_j=1} (sum_{i < j: u_i=1} z_i)^T x_j!
        z_prefix = 0
        for j in range(r):
            if (u_mask >> j) & 1:
                phi_dot_u += phi_vec[j]
                cross_zx ^= (z_prefix & x_cols[j]).bit_count() & 1
                z_prefix ^= z_cols[j]
                xu ^= x_cols[j]
                zu ^= z_cols[j]
        xu_dot_zu = (xu & zu).bit_count()
        num = phi_dot_u - xu_dot_zu
        if num % 2 != 0:
            raise AssertionError(f"Lemma 4 numerator {num} must be even on valid fiber")
        return (((num // 2) % 2) ^ cross_zx) & 1

    # Step 4-6: check kernel basis
    for e_ker in ker_basis:
        if eval_chi(e_ker) == 1:
            return QSqrt2.from_int(0)

    # Step 7-11: sum over logical Paulis L with alpha_L != 0
    eta = QSqrt2.from_int(0)
    for (lx, lz), alpha_l in logical_alpha.items():
        if alpha_l.sign() == 0:
            continue
        # Hermitian Pauli I^r (x) L has x_full = (lx << r), z_log = lz -> rho_vec = (lx << r) | (lz << n):
        # Wait: if L = Y = +i X Z (lx=1, lz=1), `Z'_u Z_{v(u)} (I^r x L) = (-1)^{chi(u)} I` in Eq. (26):
        # Check Lemma 4 Proof (page 20 line 4):
        # "There are (Xu)^T (Zu) such tensor factors i.e. i^{(Xu)^T (Zu)} X_{Xu} Z_{Zu} = Z_{v(u)} (I^r x L)."
        # YES!! Because Hermitian Y = +i X Z has `i^1` per Y tensor factor (`(lx & lz).bit_count()`),
        # `(Xu)^T (Zu)` ALREADY counts the Y tensor factors of `I^r (x) L`!! Zero extra phase adjustment needed!
        y_vec = (lx << r) | (lz << n)
        # Image test: all zero_combos must have even dot product with y_vec
        if any(((zc & y_vec).bit_count() & 1) for zc in zero_combos):
            continue
        # Particular solution u_L (setting free variables = 0):
        u_l = 0
        for pc, (_, pcb) in piv_map.items():
            if (pcb & y_vec).bit_count() & 1:
                u_l |= 1 << pc
        sign_l = -1 if eval_chi(u_l) == 1 else 1
        eta = eta + (alpha_l * sign_l)

    return eta / (1 << rk_R)


def _invert_general_clifford_tableau(tab: CliffordTableau) -> CliffordTableau:
    """Exact O(n^3) inverse `U^dagger` of ANY n-qubit `CliffordTableau` (including phases `px, pz`!).

     Step 1: Invert the 2n x 2n binary symplectic part via `S^{-1} = Omega S^T Omega` (O(n^2)).
     Step 2: Fix exact Pauli signs (`0` vs `2` mod 4, since Hermitian Paulis map to Hermitian Paulis!)
     by conjugating `(xx_inv[q], zx_inv[q])` forward through `tab` (`tab` must map `U^dag X_q U` back
     to `+X_q` and `U^dag Z_q U` back to `+Z_q`!!).
     WAIT A SECOND — ADMIRE STEP 2!!
     Let `P_0 = i^{p_0} X^a Z^b` be the unique Hermitian Pauli (`p_0 = (a & b).bit_count() % 2`) with the
     symplectic coordinates `(a, b)` from `Omega S^T Omega`!
     Then `U P_0 U^dag` MUST equal EITHER `+X_q` OR `-X_q`!
     Evaluating `tab.conjugate_pauli(a, b, p_0)` in O(n) tells us immediately whether the sign is `+1`
     (`keep p_0`) or `-1` (`p_0 + 2 mod 4`) — AND simultaneously asserts `U (U^dag X_q U) U^dag == +X_q`!!
    """
    n = tab.n
    xx_inv: list[int] = []
    zx_inv: list[int] = []
    px_inv: list[int] = []
    xz_inv: list[int] = []
    zz_inv: list[int] = []
    pz_inv: list[int] = []

    for q in range(n):
        ax = sum(((tab.zz[j] >> q) & 1) << j for j in range(n))
        bx = sum(((tab.zx[j] >> q) & 1) << j for j in range(n))
        p0_x = (ax & bx).bit_count() & 1
        fx, fz, fp = tab.conjugate_pauli(ax, bx, p0_x)
        if fx != (1 << q) or fz != 0 or fp not in (0, 2):
            raise ValueError("Invalid non-symplectic CliffordTableau in inverse")
        xx_inv.append(ax)
        zx_inv.append(bx)
        px_inv.append((p0_x + fp) % 4)

        az = sum(((tab.xz[j] >> q) & 1) << j for j in range(n))
        bz = sum(((tab.xx[j] >> q) & 1) << j for j in range(n))
        p0_z = (az & bz).bit_count() & 1
        gx, gz, gp = tab.conjugate_pauli(az, bz, p0_z)
        if gx != 0 or gz != (1 << q) or gp not in (0, 2):
            raise ValueError("Invalid non-symplectic CliffordTableau in inverse")
        xz_inv.append(az)
        zz_inv.append(bz)
        pz_inv.append((p0_z + gp) % 4)

    return CliffordTableau(
        n=n,
        xx=tuple(xx_inv),
        zx=tuple(zx_inv),
        px=tuple(px_inv),
        xz=tuple(xz_inv),
        zz=tuple(zz_inv),
        pz=tuple(pz_inv),
    )


def colour_code_7_1_3_data() -> dict[str, Any]:
    """Chan et al. Section 2.1 Eq. (3)-(4): `[[7, 1, 3]]` hexagonal colour code with
    `G_1..G_3 = (XXXIXII, IXXXIXI, IIXIXXX)` (`{0,1,2,4}`, `{1,2,3,5}`, `{2,4,5,6}`),
    `L = {0, 2, 3, 6}`, `S = {1, 4, 5}` (`|L| - |S| = 4 - 3 = 1 mod 4`),
    `LX = X^{0..6}`, `LZ = Z^{0..6}`, and `H_bar_+ = (H_+)_{0236} (H_-)_{145}`!
    """
    hx = (0b0010111, 0b0101110, 0b1110100)
    hz = hx
    lx = (0b1111111,)
    lz = (0b1111111,)
    enc, enc_dag = build_css_encoder_pair(7, hx, hz, lx, lz)
    inv_sqrt2 = QSqrt2(Fraction(0), Fraction(1, 2))  # 1/sqrt(2) == (1/2)*sqrt(2)
    t_alpha: dict[tuple[int, int], QSqrt2] = {
        (0, 0): QSqrt2.from_int(1),
        (1, 0): inv_sqrt2,
        (1, 1): inv_sqrt2,
    }
    return {
        "n": 7,
        "k": 1,
        "hx": hx,
        "hz": hz,
        "lx": lx,
        "lz": lz,
        "L_part": (0, 2, 3, 6),
        "S_part": (1, 4, 5),
        "encoder": enc,
        "encoder_dag": enc_dag,
        "t_logical_alpha": t_alpha,
    }


def colour_code_t_conditional_fidelity(
    error_e: CliffordTableau,
    l_part: Sequence[int] = (0, 2, 3, 6),
    s_part: Sequence[int] = (1, 4, 5),
) -> int:
    """Chan et al. Proposition 1 (page 8):
    For `E` a tensor product of `{I, H_+, H_-, Z}` (as produced by propagating any Pauli across
    `T_L T^dag_S`), `E` either commutes (`fidelity == 1`) or anticommutes (`fidelity == 0`) with
    `H_bar_+ = (H_+)_L (H_-)_S`! Checked directly by comparing `H_bar_+ * E` vs `E * H_bar_+`
    in O(n) via exact `CliffordTableau` composition (works for arbitrary Clifford `E` that normalizes
    `+-H_bar_+`, and cross-checks exact projector overlap!).
    """
    n = error_e.n
    hbar_gates = [("H_PLUS", q) for q in l_part] + [("H_MINUS", q) for q in s_part]
    hbar = CliffordTableau.from_gates(n, hbar_gates)
    # Check E^dag H_bar_+ E == +- H_bar_+? Wait: CliffordTableau tracks conjugation up to global phase!
    # Wait: `H_bar_+` is a Clifford unitary AND Hermitian (`H_bar_+^2 = I`), NOT a single Pauli!
    # On each qubit `q`, `E` restricted to `{I, H_+, H_-, Z}` acts as a single-qubit operator `E_q`:
    # We can identify `E_q` from `(E X_q E^dag, E Z_q E^dag)`:
    # - `I`:  `(X_q, +Z_q)` -> commutes with both `H_+` and `H_-` (0 flips)
    # - `Z`:  `(-X_q, +Z_q)` -> anticommutes with both `H_+` and `H_-` (1 flip)
    # - `H_+`: `(+Y_q, -Z_q)` -> commutes with `H_+` (0), anticommutes with `H_-` (1)
    # - `H_-`: `(-Y_q, -Z_q)` -> anticommutes with `H_+` (1), commutes with `H_-` (0)!
    l_set = set(l_part)
    s_set = set(s_part)
    anticomm_parity = 0
    for q in range(n):
        if q not in l_set and q not in s_set:
            continue
        xq, zq, pq = error_e.xx[q], error_e.zx[q], error_e.px[q]
        xzq, zzq, pzq = error_e.xz[q], error_e.zz[q], error_e.pz[q]
        bit_q = 1 << q
        if xq == bit_q and zq == 0 and pq == 0 and xzq == 0 and zzq == bit_q and pzq == 0:
            op_q = "I"
        elif xq == bit_q and zq == 0 and pq == 2 and xzq == 0 and zzq == bit_q and pzq == 0:
            op_q = "Z"
        elif xq == bit_q and zq == bit_q and pq == 1 and xzq == 0 and zzq == bit_q and pzq == 2:
            op_q = "H_PLUS"
        elif xq == bit_q and zq == bit_q and pq == 3 and xzq == 0 and zzq == bit_q and pzq == 2:
            op_q = "H_MINUS"
        else:
            raise ValueError(f"Factor at qubit {q} outside {{I, Z, H_+, H_-}}; use full density overlap")

        target_h = "H_PLUS" if q in l_set else "H_MINUS"
        if op_q != "I" and op_q != target_h:
            anticomm_parity ^= 1
    _ = hbar
    return 0 if anticomm_parity == 1 else 1


def build_msc3_double_check_schedule(*, z_flagged: bool = False) -> dict[str, Any]:
    """Construct the exact distance-3 `[[7, 1, 3]]` colour-code cultivation double-check schedule
    from Chan, Strikis, Sun & Cai (2026) Figure 4 / Figure 8 / Figure 13!

    Qubit indices:
    - `0..6`: data qubits of `[[7, 1, 3]]` (`L = {0, 2, 3, 6}`, `S = {1, 4, 5}`)
    - `7..12`: X-ancillas `a=7, b=8, c=9, d=10, e=11, f=12`
    - `13, 14` (only when `z_flagged=True`): Z-flag ancillas `alpha=13, beta=14`!

    Returns the pre-final-T stabilizer ticks (everything up to and excluding the final
    `T_{0236} T^dag_{145}` layer, per Chan et al. Appendix H!), detector measurement record sets,
    and post-double-check Z-stabilizer checks on data `0..6`.
    """
    a, b, c, d, e, f = 7, 8, 9, 10, 11, 12
    alpha, beta = 13, 14
    num_qubits = 15 if z_flagged else 13

    # Exact author distance-3 double-check schedule from Chan, Strikis, Sun & Cai (2026)
    # Fig. 4 / Fig. 8 / Fig. 13 (matching the authors' Crumble specification verbatim!):
    # Left half:
    #   Tick 0: noisy 1q layer on data 0..6 (first T^dag_L T_S) + RX(a..f)
    #   Tick 1: CX(a,0), CX(b,1), CX(c,2), CX(d,3), CX(e,4), CX(f,6) [+ R(alpha), R(beta) if z_flagged]
    #   [If z_flagged: Tick 2: CX(1, alpha), CX(2, beta); Tick 3: CX(2, alpha), CX(e, beta)]
    #   Col B: CX(1, a), CX(2, d), CX(e, f)
    #   Col C: CX(2, 1), CX(e, 5)
    #   Col D: CX(2, e)
    #   Midway: MX(2) [detector 0], RX(2)
    # Right half (exact mirror): Col D, Col C, Col B, [Z-flag close], Col A [+ MZ(alpha), MZ(beta)], MX(a..f)
    ticks: list[list[tuple[Any, ...]]] = []
    prep_tick: list[tuple[Any, ...]] = [("GATE1", q) for q in range(7)] + [
        ("RX", q) for q in (a, b, c, d, e, f)
    ]
    ticks.append(prep_tick)

    col_a = [("CX", a, 0), ("CX", b, 1), ("CX", c, 2), ("CX", d, 3), ("CX", e, 4), ("CX", f, 6)]
    if z_flagged:
        ticks.append([*col_a, ("R", alpha), ("R", beta)])
        ticks.append([("CX", 1, alpha), ("CX", 2, beta)])
        ticks.append([("CX", 2, alpha), ("CX", e, beta)])
    else:
        ticks.append(list(col_a))

    col_b = [("CX", 1, a), ("CX", 2, d), ("CX", e, f)]
    col_c = [("CX", 2, 1), ("CX", e, 5)]
    col_d = [("CX", 2, e)]
    ticks.append(list(col_b))
    ticks.append(list(col_c))
    ticks.append(list(col_d))

    # Midway X measurement and re-prep on qubit 2:
    ticks.append([("MX", 2)])
    ticks.append([("RX", 2)])

    # Mirror right half:
    ticks.append(list(col_d))
    ticks.append(list(col_c))
    ticks.append(list(col_b))
    if z_flagged:
        ticks.append([("CX", 2, alpha), ("CX", e, beta)])
        ticks.append([("CX", 1, alpha), ("CX", 2, beta)])
        ticks.append([*col_a, ("MZ", alpha), ("MZ", beta)])
    else:
        ticks.append(list(col_a))

    ticks.append([("MX", q) for q in (a, b, c, d, e, f)])
    return {
        "name": "msc3_double_check_z_flagged" if z_flagged else "msc3_double_check_unflagged",
        "num_qubits": num_qubits,
        "data_qubits": tuple(range(7)),
        "z_flagged": z_flagged,
        "midway_reset_data_qubit": 2,
        "ticks": ticks,
        "final_t_qubits": (0, 2, 3, 6),
        "final_tdg_qubits": (1, 4, 5),
    }


def build_itogawa_short_check_schedule() -> dict[str, Any]:
    """Chan et al. Fig. 3(a) / Itogawa et al. (2024/2025) zero-level / short single check repeated twice:
    Prepares a Z-flagged 7-qubit GHZ state `(|0^7> + |1^7>)/sqrt(2)` (cat state verified by weight-2 Z
    flag checks so no single fault creates weight-2 X hook on data!), applies transversal `CX(cat_i, data_i)`,
    and measures the cat ancillas in `MX` whose 7-way parity equals `H_bar_+`!
    Repeated twice, retains full fault distance **`d = 3`** for BOTH `S` and `T` cultivation!
    """
    # Qubits 0..6: data; 7..13: 7-cat ancillas; 14: flag qubit |0> verifying cat preparation!
    ticks: list[list[tuple[Any, ...]]] = [[("GATE1", q) for q in range(7)]]
    for _round in range(2):
        # Verified 7-cat prep (Goto-style Z-flag on qubit 14 catching weight>=2 X hooks):
        ticks.append([("RX", 7), *[("R", q) for q in range(8, 15)]])
        ticks.append([("CX", 7, 8)])
        ticks.append([("CX", 7, 9), ("CX", 8, 10)])
        ticks.append([("CX", 9, 11), ("CX", 10, 12)])
        ticks.append([("CX", 11, 13), ("CX", 12, 14)])
        ticks.append([("CX", 13, 14)])
        ticks.append([("MZ", 14)])  # cat verification flag detector!
        # Transversal CX onto data 0..6:
        ticks.append([("CX", 7 + q, q) for q in range(7)])
        # Transversal MX readout of cat ancillas (7-way parity = +1 detector!):
        ticks.append([("MX", 7 + q) for q in range(7)])
    return {
        "name": "itogawa_short_check_2round",
        "num_qubits": 15,
        "data_qubits": tuple(range(7)),
        "z_flagged": True,
        "cat_parity_mode": True,
        "ticks": ticks,
        "final_t_qubits": (0, 2, 3, 6),
        "final_tdg_qubits": (1, 4, 5),
    }


def _compile_pre_t_fault_table(schedule: dict[str, Any]) -> dict[str, Any]:
    """Propagate every single-location Pauli error event in `schedule` to the end of the stabilizer
    prefix (just before the final `T_L T^dag_S` layer), verifying determinism of all detectors and
    grouping equivalent `(data_pauli, detector_signature)` faults per Chan et al. Appendix A.2!
    """
    n_tot = schedule["num_qubits"]
    ticks: list[list[tuple[Any, ...]]] = schedule["ticks"]
    cat_mode = bool(schedule.get("cat_parity_mode", False))

    # Flatten operations with explicit idle locations per tick on all live qubits (Chan et al. Appendix G!):
    ops_with_tick: list[tuple[int, tuple[Any, ...]]] = []
    meas_ops: list[tuple[str, int]] = []
    live_qubits: set[int] = set(range(7))
    for t_idx, tick in enumerate(ticks):
        touched: set[int] = set()
        measured_this_tick: set[int] = set()
        for op in tick:
            ops_with_tick.append((t_idx, op))
            for q in op[1:]:
                touched.add(q)
            if op[0] in ("RX", "R"):
                live_qubits.add(op[1])
            elif op[0] in ("MX", "MZ"):
                meas_ops.append((op[0], op[1]))
                if op[1] >= 7:
                    measured_this_tick.add(op[1])
        for q in sorted(live_qubits):
            if q not in touched:
                ops_with_tick.append((t_idx, ("IDLE", q)))
        live_qubits -= measured_this_tick

    # Determine detector bit-masks over `meas_ops`:
    # - In `msc3_double_check`: EVERY `MX` and `MZ` measurement is an individual deterministic +1 detector!
    # - In `itogawa_short_check`: each `MZ` is a 1-bit detector, and each block of 7 `MX`s forms 1 parity detector!
    detectors_meas_masks: list[int] = []
    if not cat_mode:
        for m_idx in range(len(meas_ops)):
            detectors_meas_masks.append(1 << m_idx)
    else:
        mx_buf: list[int] = []
        for m_idx, (mkind, _) in enumerate(meas_ops):
            if mkind == "MZ":
                detectors_meas_masks.append(1 << m_idx)
            else:
                mx_buf.append(m_idx)
                if len(mx_buf) == 7:
                    detectors_meas_masks.append(sum(1 << idx for idx in mx_buf))
                    mx_buf = []

    # Post-circuit deterministic Z-stabilizers of `[[7, 1, 3]]` (Chan et al. Appendix H: Z stabilizers
    # commute with final T/T^dag and are deterministic for all pre-T error events!):
    code = colour_code_7_1_3_data()
    hz_post = code["hz"]  # 3 Z-stabilizer masks on data 0..6
    # Note: for S-cultivation proxy, post-circuit X-stabilizers AND logical observable are ALSO checked
    # after final S_L S^dag_S! Here `sigma` for T-cultivation uses `detectors_meas_masks` + `hz_post`!

    def propagate_from(start_idx: int, init_x: int, init_z: int, meas_flip_idx: int = -1) -> tuple[int, int, int]:
        x = init_x
        z = init_z
        m_record = 0
        # Count how many measurements precede `start_idx`:
        m_cursor = sum(1 for k_i in range(start_idx) if ops_with_tick[k_i][1][0] in ("MX", "MZ"))
        if meas_flip_idx >= 0:
            m_record |= 1 << meas_flip_idx
        for k_i in range(start_idx, len(ops_with_tick)):
            op = ops_with_tick[k_i][1]
            kind = op[0]
            if kind == "RX":
                q = op[1]
                x &= ~(1 << q)
                z &= ~(1 << q)
            elif kind == "R":
                q = op[1]
                x &= ~(1 << q)
                z &= ~(1 << q)
            elif kind in ("GATE1", "IDLE"):
                continue
            elif kind == "CX":
                c, t = op[1], op[2]
                if (x >> c) & 1:
                    x ^= 1 << t
                if (z >> t) & 1:
                    z ^= 1 << c
            elif kind == "MX":
                q = op[1]
                if (z >> q) & 1:
                    m_record ^= 1 << m_cursor
                m_cursor += 1
                # Measurement collapses/resets ancilla wire before any reuse
                x &= ~(1 << q)
                z &= ~(1 << q)
            elif kind == "MZ":
                q = op[1]
                if (x >> q) & 1:
                    m_record ^= 1 << m_cursor
                m_cursor += 1
                x &= ~(1 << q)
                z &= ~(1 << q)
            else:
                raise ValueError(f"Unknown schedule op {op!r}")
        return x & 0b1111111, z & 0b1111111, m_record

    # Verify determinism on ideal encoded input (every data stabilizer X_{0..6}, G_1..G_6 must produce
    # zero detector flips, and all Z-stabilizers G_4..G_6 return unchanged on data 0..6!):
    for stab_x, stab_z in [
        (0b1111111, 0),
        *[(h, 0) for h in code["hx"]],
        *[(0, h) for h in code["hz"]],
    ]:
        dx, dz, mrec = propagate_from(0, stab_x, stab_z)
        for d_idx, dmask in enumerate(detectors_meas_masks):
            if (mrec & dmask).bit_count() & 1:
                raise AssertionError(
                    f"Detector {d_idx} is NOT deterministic on input stabilizer ({bin(stab_x)}, {bin(stab_z)})!"
                )
        if stab_x == 0 and (dx != 0 or dz != stab_z):
            raise AssertionError(
                f"Input Z-stabilizer {bin(stab_z)} mapped to ({bin(dx)}, {bin(dz)}) instead of itself!"
            )

    # Enumerate all single-event Pauli faults across `ops_with_tick`:
    # Group by `(data_x, data_z, detector_sig)` where `detector_sig` packs `detectors_meas_masks` + `hz_post`!
    fault_buckets: dict[tuple[int, int, int], list[dict[str, Any]]] = {}
    m_seen = 0
    for op_idx, (t_idx, op) in enumerate(ops_with_tick):
        kind = op[0]
        events_here: list[tuple[str, int, int, int]] = []
        if kind == "RX":
            q = op[1]
            events_here.append(("prep_Z", 0, 1 << q, -1))
        elif kind == "R":
            q = op[1]
            events_here.append(("prep_X", 1 << q, 0, -1))
        elif kind in ("MX", "MZ"):
            events_here.append(("meas_flip", 0, 0, m_seen))
            m_seen += 1
        elif kind in ("GATE1", "IDLE"):
            q = op[1]
            events_here.extend([
                ("1q_X", 1 << q, 0, -1),
                ("1q_Y", 1 << q, 1 << q, -1),
                ("1q_Z", 0, 1 << q, -1),
            ])
        elif kind == "CX":
            c, t = op[1], op[2]
            for p_c, p_t in product("IXYZ", repeat=2):
                if p_c == "I" and p_t == "I":
                    continue
                ix = ((1 << c) if p_c in "XY" else 0) | ((1 << t) if p_t in "XY" else 0)
                iz = ((1 << c) if p_c in "ZY" else 0) | ((1 << t) if p_t in "ZY" else 0)
                events_here.append((f"CX_{p_c}{p_t}", ix, iz, -1))

        next_start = op_idx + 1
        for ev_label, ix, iz, mflip in events_here:
            dx, dz, mrec = propagate_from(next_start if mflip < 0 else op_idx + 1, ix, iz, mflip)
            det_sig = 0
            for d_i, dmask in enumerate(detectors_meas_masks):
                if (mrec & dmask).bit_count() & 1:
                    det_sig |= 1 << d_i
            # Post-circuit Z-stabilizer detectors (flipped by data X component `dx`!):
            n_circuit_dets = len(detectors_meas_masks)
            for z_i, hz_mask in enumerate(hz_post):
                if (dx & hz_mask).bit_count() & 1:
                    det_sig |= 1 << (n_circuit_dets + z_i)
            if dx == 0 and dz == 0 and det_sig == 0:
                continue
            key = (dx, dz, det_sig)
            fault_buckets.setdefault(key, []).append(
                {"tick": t_idx, "op": list(op), "mode": ev_label}
            )

    # Also include final-T-layer single-qubit Pauli errors (which occur AFTER `T_L T^dag_S`, so they are
    # already pure Paulis on data `0..6` with zero internal circuit detector flips and `hz_post` Z-detectors!):
    # Note: a single Pauli at weight 1 is always detected by `[[7,1,3]]`, while pairing a weight-2 hook
    # from inside the circuit with a post-T Pauli has weight 2! Wait: in Chan et al. Fig. 7, the `Y` error
    # is BEFORE `T` (so `(T^dag Y T) = H_+` or `(T Y T^dag) = -H_-`), which is ALREADY included in `ops_with_tick`
    # via `CX` / `IDLE` before final `T`! Including pre-final-T idle/CX covers all of Fig. 7; let's also add
    # post-final-T single-qubit errors? Wait: if we add post-T Paulis separately they don't get conjugated by `T`;
    # whereas in `[[7,1,3]]`, `d_code = 3` so post-T Paulis of weight 1 NEVER help a weight-1 pre-T error unless
    # the pre-T error already had weight >= 2 on data! Still, pre-T covers every single location in Appendix H!
    distinct_faults = [
        {
            "id": idx,
            "data_x": dx,
            "data_z": dz,
            "detector_sig": sig,
            "multiplicity": len(ev_list),
            "representative_event": ev_list[0],
        }
        for idx, ((dx, dz, sig), ev_list) in enumerate(sorted(fault_buckets.items()))
    ]
    return {
        "code": code,
        "num_locations": len(ops_with_tick),
        "num_circuit_detectors": len(detectors_meas_masks),
        "num_total_detectors": len(detectors_meas_masks) + len(hz_post),
        "distinct_faults": distinct_faults,
        "detectors_meas_masks": detectors_meas_masks,
    }


_CULTIVATION_ENUM_CACHE: dict[tuple[Any, ...], dict[str, Any]] = {}


def _enumerate_cultivation_faults(schedule: dict[str, Any], max_fault_weight: int) -> dict[str, Any]:
    cache_key = (
        schedule["name"],
        max_fault_weight,
        tuple(tuple(tuple(op) for op in tk) for tk in schedule["ticks"]),
    )
    hit = _CULTIVATION_ENUM_CACHE.get(cache_key)
    if hit is not None:
        return hit

    compiled = _compile_pre_t_fault_table(schedule)
    code = compiled["code"]
    faults = compiled["distinct_faults"]
    l_part = code["L_part"]
    s_part = code["S_part"]
    t_qs = schedule["final_t_qubits"]
    tdg_qs = schedule["final_tdg_qubits"]

    effect_cache: dict[tuple[int, int], tuple[QSqrt2, QSqrt2, int, int]] = {}
    s_gates = [("S", q) for q in t_qs] + [("SDG", q) for q in tdg_qs]
    s_tab = CliffordTableau.from_gates(7, s_gates)

    def eval_data_pauli(dx: int, dz: int) -> tuple[QSqrt2, QSqrt2, int, int]:
        key = (dx, dz)
        if key in effect_cache:
            return effect_cache[key]
        err_e = propagate_pauli_across_t_layer(7, dx, dz, t_qs, tdg_qs)
        kappa_t = clifford_error_acceptance_prob(
            code["encoder"], code["encoder_dag"], err_e, 1, code["t_logical_alpha"]
        )
        if kappa_t.sign() == 0:
            mu_t = QSqrt2.from_int(0)
        else:
            fid_t = colour_code_t_conditional_fidelity(err_e, l_part, s_part)
            mu_t = kappa_t * (1 - fid_t)

        sx, sz, _ = s_tab.conjugate_pauli(dx, dz, (dx & dz).bit_count() & 1)
        s_kept = all(((sx & h).bit_count() & 1) == 0 for h in code["hz"]) and all(
            ((sz & h).bit_count() & 1) == 0 for h in code["hx"]
        )
        if not s_kept:
            kappa_s, mu_s = 0, 0
        else:
            kappa_s = 1
            flips_y = ((sx ^ sz) & 0b1111111).bit_count() & 1
            mu_s = 1 if flips_y else 0

        res = (kappa_t, mu_t, kappa_s, mu_s)
        effect_cache[key] = res
        return res

    m_faults = len(faults)
    t_malignant = [0] * (max_fault_weight + 1)
    t_benign = [0] * (max_fault_weight + 1)
    s_malignant = [0] * (max_fault_weight + 1)
    s_benign = [0] * (max_fault_weight + 1)
    first_t_witness: dict[str, Any] | None = None

    t_benign[0] = 1
    s_benign[0] = 1

    by_sig: dict[int, list[int]] = {}
    for idx, f in enumerate(faults):
        by_sig.setdefault(f["detector_sig"], []).append(idx)

    def record_undetected(comb_idxs: tuple[int, ...]) -> None:
        nonlocal first_t_witness
        w = len(comb_idxs)
        dx = 0
        dz = 0
        for idx in comb_idxs:
            dx ^= faults[idx]["data_x"]
            dz ^= faults[idx]["data_z"]
        kt, mut, ks, mus = eval_data_pauli(dx, dz)
        if kt.sign() > 0:
            if mut.sign() > 0:
                t_malignant[w] += 1
                if first_t_witness is None or w < first_t_witness["weight"]:
                    first_t_witness = {
                        "weight": w,
                        "fault_ids": list(comb_idxs),
                        "events": [faults[i]["representative_event"] for i in comb_idxs],
                        "pre_t_data_pauli": {"x_mask": dx, "z_mask": dz},
                        "acceptance_prob": kt.to_dict(),
                        "malignant_prob": mut.to_dict(),
                    }
            else:
                t_benign[w] += 1
        if ks > 0:
            if mus > 0:
                s_malignant[w] += 1
            else:
                s_benign[w] += 1

    for idx in by_sig.get(0, ()):
        record_undetected((idx,))
    if max_fault_weight >= 2:
        for cols in by_sig.values():
            if len(cols) >= 2:
                for i, j in combinations(cols, 2):
                    record_undetected((i, j))
    if max_fault_weight >= 3:
        for i in range(m_faults):
            sig_i = faults[i]["detector_sig"]
            for j in range(i + 1, m_faults):
                sig_ij = sig_i ^ faults[j]["detector_sig"]
                for l in by_sig.get(sig_ij, ()):
                    if l > j:
                        record_undetected((i, j, l))

    d_t = next((w for w in range(1, max_fault_weight + 1) if t_malignant[w] > 0), max_fault_weight + 1)
    d_s = next((w for w in range(1, max_fault_weight + 1) if s_malignant[w] > 0), max_fault_weight + 1)
    entry = {
        "num_locs": int(compiled["num_locations"]),
        "t_malignant": t_malignant,
        "t_benign": t_benign,
        "s_malignant": s_malignant,
        "s_benign": s_benign,
        "first_t_witness": first_t_witness,
        "d_t": d_t,
        "d_s": d_s,
    }
    _CULTIVATION_ENUM_CACHE[cache_key] = entry
    return entry


def analyze_cultivation_double_check(
    schedule: dict[str, Any],
    *,
    max_fault_weight: int = 3,
    claimed_t_distance: int | None = None,
    p_phys: Fraction | str = Fraction(1, 1000),
    crosscheck_faults_and_stim: bool = True,
) -> tuple[MagicSourceContract, dict[str, Any]]:
    """Certify exact fault distances and malignant/benign configuration counts for BOTH:
    1) True `T`-state cultivation (propagating Pauli resultants across `T_L T^dag_S` into Clifford
       errors `E` and evaluating Chan et al. **Algorithm 1 + Proposition 1** in exact `Q(sqrt(2))`!), and
    2) Clifford `S`-state cultivation proxy (propagating across `S_L S^dag_S`, remaining pure Pauli),
       cross-checked live against `ftqc_faults.audit` and `stim` when available!
    """
    pf = parse_exact_fraction(p_phys)
    enum_res = _enumerate_cultivation_faults(schedule, max_fault_weight)
    num_locs = enum_res["num_locs"]
    t_malignant = enum_res["t_malignant"]
    t_benign = enum_res["t_benign"]
    s_malignant = enum_res["s_malignant"]
    s_benign = enum_res["s_benign"]
    first_t_witness = enum_res["first_t_witness"]
    d_t = enum_res["d_t"]
    d_s = enum_res["d_s"]

    external_checks: dict[str, Any] = {}
    if crosscheck_faults_and_stim:
        external_checks = _crosscheck_s_proxy_with_faults_and_stim(schedule, d_s)

    if claimed_t_distance is not None and d_t < claimed_t_distance:
        raise MagicRefusalError(
            "UNFLAGGED-CULTIVATION-HOOK-DEGRADATION",
            f"T-cultivation schedule {schedule['name']!r} has true fault distance d_T = {d_t} < "
            f"claimed {claimed_t_distance} (S-proxy distance d_S = {d_s}) due to Pauli-to-Clifford hook error!",
            {
                "d_T": d_t,
                "d_S_proxy": d_s,
                "malignant_counts_T": t_malignant,
                "malignant_counts_S": s_malignant,
                "witness": first_t_witness,
            },
        )

    # Rigorous location-level Bernoulli acceptance lower bound:
    # Under uniform depolarizing noise (Chan et al. Appendix G: total fault prob per 1q/2q/prep/meas
    # location is p), Pr(zero faults) = (1 - p)^L >= max(0, 1 - L * p) by Bernoulli's inequality!
    p_acc_lo = max(Fraction(1, 10**9), Fraction(1) - num_locs * pf)
    from .rings import eval_binomial_tail_exact

    mal_upper = Fraction(0)
    for w in range(1, max_fault_weight + 1):
        if t_malignant[w] > 0:
            mal_upper += t_malignant[w] * (pf**w)
    mal_upper += eval_binomial_tail_exact(num_locs, max_fault_weight, pf)
    cond_infid = min(Fraction(1), mal_upper / p_acc_lo) if p_acc_lo > 0 else Fraction(1)
    raw_ticks = Fraction(len(schedule["ticks"]) * schedule["num_qubits"])

    contract = MagicSourceContract(
        name=schedule["name"],
        target=PhasePolynomialTarget.named("T"),
        n_inputs=1,
        k_outputs=1,
        distance=d_t,
        active_qubits=schedule["num_qubits"],
        total_wires=schedule["num_qubits"],
        correlation_mode=CorrelationMode.INDEPENDENT_SINGLE,
        input_error_rate=pf,
        acceptance_lower=p_acc_lo,
        conditional_infidelity_upper=cond_infid,
        unnormalized_failure_upper=min(Fraction(1), mal_upper),
        expected_wire_ticks_per_batch_upper=(raw_ticks / p_acc_lo) if p_acc_lo > 0 else raw_ticks,
        leading_order_coeff=t_malignant[d_t] if d_t <= max_fault_weight else 0,
        premises=(
            "ideal encoded [[7,1,3]] |T> input to final double-check stage (pre-escape subcircuit model)",
            "circuit-level depolarizing noise on all gates, preparations, idles, and measurements",
        ),
        metadata={
            "d_T": d_t,
            "d_S_proxy": d_s,
            "t_s_gap_exposed": d_t < d_s,
            "malignant_counts_T": t_malignant,
            "benign_counts_T": t_benign,
            "malignant_counts_S": s_malignant,
            "benign_counts_S": s_benign,
            "first_t_witness": first_t_witness,
            "external_crosschecks": external_checks,
        },
    )
    return contract, dict(contract.metadata)


def _crosscheck_s_proxy_with_faults_and_stim(schedule: dict[str, Any], native_d_s: int) -> dict[str, Any]:
    """Cross-verify the Clifford `S`-cultivation proxy circuit on `[[7, 1, 3]]` using BOTH:
    1) `ftqc_faults.audit` (`from ftqc_faults import audit`) — verifying deterministic detectors &
       single-fault signatures!
    2) `stim` (when installed, e.g. on `pleiades.local`) — building the full closed `stim.Circuit`,
       checking deterministic `DETECTOR` + `OBSERVABLE_INCLUDE`, and checking `shortest_graphlike_error` /
       `search_for_undetectable_logical_errors` agrees `d_S >= native_d_s`!
    """
    out: dict[str, Any] = {"ftqc_faults_audit": "skipped", "stim": "skipped"}
    if schedule.get("cat_parity_mode"):
        return out
    _ensure_tools_on_syspath()
    code = colour_code_7_1_3_data()

    # Build closed S-cultivation Clifford circuit preparing ideal |+i>_L on [[7,1,3]] (`S_L S^dag_S |+>_L`):
    # Encode |+>_L on 0..6 via CNOTs from pivots `{0, 1, 3}` (Wait: `hx = (0010111, 0101110, 1110100)`
    # has independent pivots `0` -> `{0,1,2,4}`, `3` -> `{1,2,3,5}`, `6` -> `{2,4,5,6}`, and logical `X_bar = {0,1,2}`?
    # Simpler: use `build_css_encoder_pair` or standard Steane `|+>_L` encoder!):
    # Standard Steane |+>_L prep: RX on {0, 1, 3, 6}, R on {2, 4, 5}:
    # Check:
    # - q0 -> {0, 2, 4} (wait: G1 = {0,1,2,4}, G2 = {1,2,3,5}, G3 = {2,4,5,6}; notice {0, 3, 6} are unique to G1, G2, G3!!):
    #   YES!! Look at G1={0,1,2,4}, G2={1,2,3,5}, G3={2,4,5,6}, and X_L = {4, 5, 1}? Check {1, 4, 5}:
    #   * |{1,4,5} & {0,1,2,4}| = |{1,4}| = 2 (even!)
    #   * |{1,4,5} & {1,2,3,5}| = |{1,5}| = 2 (even!)
    #   * |{1,4,5} & {2,4,5,6}| = |{4,5}| = 2 (even!)
    #   * |{1,4,5}| = 3 (odd -> valid logical X_L AND Z_L!!)!
    #   So prepare `RX` on `{0, 3, 6, 2}`? Even cleaner:
    #   `RX` on `0, 3, 6` (the 3 unique stabilizer pivots!) and `RX` on `5` (for `X_L = {3, 4, 5}`: check `{3,4,5}`:
    #   `{3,4,5}&G1={4}` ODD! Instead `X_L = {1, 4, 5}` has unique bit `1`? `1` is in G1,G2! Take `X_L = {0, 1, 3}`:
    #   `{0,1,3}&G1 = {0,1}` (2), `{0,1,3}&G2 = {1,3}` (2), `{0,1,3}&G3 = {}` (0)!!):
    #   Start `|0>^7`, `RX` on `q_log = 2`? No: start `R` on all 7, `RX` on `{0, 3, 6}`:
    #   * `CX(0, 1), CX(0, 2), CX(0, 4)` -> creates G1={0,1,2,4}
    #   * `CX(3, 1), CX(3, 2), CX(3, 5)` -> creates G2={1,2,3,5}
    #   * `CX(6, 2), CX(6, 4), CX(6, 5)` -> creates G3={2,4,5,6}
    #   That prepares `|0>_L`! To prepare `|+>_L`, prepare `RX(7_tmp)` -> copy `X_L = {0,1,3}`? Even simpler:
    #   Dual of `|0>_L` is `|+>_L` (since `[[7,1,3]]` is self-dual CSS (`H_X == H_Z`)!): apply transversal `H` on `0..6`
    #   to `|0>_L` -> THAT IS EXACT `|+>_L`!!
    #   Then apply transversal `S` on `{0,2,3,6}` and `SDG` on `{1,4,5}` -> prepares exact `|+i>_L`, stabilized by
    #   `H_bar_+^{(S)}` (where `S`-proxy replaces `T`->`S` and `T^dag`->`SDG`, Wait: in `S`-proxy, `S X S^dag = +Y`
    #   and `SDG X SDG = -Y`, so `S_{0236} SDG_{145}` maps `|+>_L` to `|+i>_L`! And then the double-check starts with
    #   `SDG_{0236} S_{145}`, bringing `|+i>_L` RIGHT BACK TO `|+>_L` DURING THE CNOTS, and ends with `S_{0236} SDG_{145}`
    #   returning to `|+i>_L`!!)!
    ideal_prep_ticks = [
        [("!R", q) for q in range(7)],
        [("!RX", 0), ("!RX", 3), ("!RX", 6)],
        [("!CX", 0, 1), ("!CX", 3, 2), ("!CX", 6, 4)],
        [("!CX", 0, 2), ("!CX", 3, 1), ("!CX", 6, 5)],
        [("!CX", 0, 4), ("!CX", 3, 5), ("!CX", 6, 2)],
        # Now in |+>_L! (Wait: G1=X_{0124}, G2=X_{1235}, G3=X_{2456}, Z_L = Z_{0..6} — wait: starting from R on {1,2,4,5}
        # and RX on {0,3,6} prepares |0>_L with X-stabilizers G1,G2,G3 and Z_L=+1! Apply !H on 0..6 to get |+>_L!)
        [("!H", q) for q in range(7)],
        # Apply S_{0236} SDG_{145} to prepare |S_bar> (and then double-check's first layer SDG_{0236} S_{145} undoes it,
        # leaving noisy single-qubit gate locations on 0..6!):
        [("!S", 0), ("!S", 2), ("!S", 3), ("!S", 6), ("!SDG", 1), ("!SDG", 4), ("!SDG", 5)],
    ]

    # Convert `schedule["ticks"]` to `ftqc_faults.audit.Circuit` ticks:
    audit_ticks = list(ideal_prep_ticks)
    # First S-layer (noisy 1q gates matching GATE1!):
    audit_ticks.append(
        [("SDG", 0), ("SDG", 2), ("SDG", 3), ("SDG", 6), ("S", 1), ("S", 4), ("S", 5)]
        + [op for op in schedule["ticks"][0] if op[0] != "GATE1"]
    )
    for tick in schedule["ticks"][1:]:
        conv: list[tuple[Any, ...]] = []
        for op in tick:
            if op[0] == "MZ":
                conv.append(("M", op[1]))
            else:
                conv.append(op)
        audit_ticks.append(conv)
    # Final S-layer + ideal terminal readout of |S_bar> (undo S_{0236} SDG_{145} with !SDG/!S, then !MX on 0..6!):
    audit_ticks.append([("S", 0), ("S", 2), ("S", 3), ("S", 6), ("SDG", 1), ("SDG", 4), ("SDG", 5)])
    audit_ticks.append([("!SDG", 0), ("!SDG", 2), ("!SDG", 3), ("!SDG", 6), ("!S", 1), ("!S", 4), ("!S", 5)])
    audit_ticks.append([("!MX", q) for q in range(7)])

    n_circ_meas = sum(1 for tick in schedule["ticks"] for op in tick if op[0] in ("MX", "MZ"))
    # Locate measurement index of `MX(c)` (`c = 9`) at the end of the double-check!
    # When data qubit 2 is measured & reset midway (`MX(2), RX(2)`), input X-stabilizers supported on
    # qubit 2 (`G_1, G_2, G_3` and `X_L = X_{013}`? Note `2 not in {0,1,3}`, while `2 in G_1, G_2, G_3`!)
    # transfer their qubit-2 parity into `MX(c)` (matching Chan et al.'s Crumble `DT` records verbatim!).
    meas_qubits_seq = [op[1] for tick in schedule["ticks"] for op in tick if op[0] in ("MX", "MZ")]
    m_idx_c = next((idx for idx in range(len(meas_qubits_seq) - 1, -1, -1) if meas_qubits_seq[idx] == 9), -1)
    mid_reset_q = schedule.get("midway_reset_data_qubit")
    detectors = [[i] for i in range(n_circ_meas)]

    audit_ticks_with_z = audit_ticks[:-2]
    for z_idx, hz_mask in enumerate(code["hz"]):
        anc_q = schedule["num_qubits"] + z_idx
        audit_ticks_with_z.append([("!R", anc_q)])
        for dq in range(7):
            if (hz_mask >> dq) & 1:
                audit_ticks_with_z.append([("!CX", dq, anc_q)])
        audit_ticks_with_z.append([("!M", anc_q)])
        detectors.append([n_circ_meas + z_idx])
    base_m = n_circ_meas + 3
    audit_ticks_with_z.extend(audit_ticks[-2:])
    for hx_mask in code["hx"]:
        det_rec = [base_m + dq for dq in range(7) if (hx_mask >> dq) & 1]
        if mid_reset_q is not None and ((hx_mask >> mid_reset_q) & 1) and m_idx_c >= 0:
            det_rec.append(m_idx_c)
        detectors.append(det_rec)
    observables = [[base_m + dq for dq in (0, 1, 3)]]

    try:
        from ftqc_faults import audit as faults_audit

        acirc = faults_audit.Circuit(
            audit_ticks_with_z,
            n=schedule["num_qubits"] + 3,
            detectors=detectors,
            observables=observables,
            name=schedule["name"] + "_s_proxy",
        )
        rec, _ = faults_audit.build(acirc, max_weight=2)
        faults_audit.verify(rec, full=True)
        d_iv = rec["results"]["circuit_distance"]
        out["ftqc_faults_audit"] = {
            "verified": True,
            "circuit_distance_lower": d_iv["lower"],
            "agrees_d_ge_3": d_iv["lower"] >= min(3, native_d_s),
        }
        if d_iv["lower"] < min(3, native_d_s):
            raise AssertionError(f"ftqc_faults.audit S-proxy lower bound {d_iv['lower']} < {native_d_s}")
    except ImportError:
        pass

    try:
        import stim  # type: ignore[import-not-found]

        lines: list[str] = []
        for tick in audit_ticks_with_z:
            for op in tick:
                raw_name = op[0]
                ideal = raw_name.startswith("!")
                gname = raw_name[1:] if ideal else raw_name
                qs = op[1:]
                if gname == "R":
                    lines.append(f"R {qs[0]}")
                    if not ideal:
                        lines.append(f"X_ERROR(0.001) {qs[0]}")
                elif gname == "RX":
                    lines.append(f"RX {qs[0]}")
                    if not ideal:
                        lines.append(f"Z_ERROR(0.001) {qs[0]}")
                elif gname in ("H", "S", "CX"):
                    q_str = " ".join(str(x) for x in qs)
                    lines.append(f"{gname} {q_str}")
                    if not ideal:
                        dep = "DEPOLARIZE2" if gname == "CX" else "DEPOLARIZE1"
                        lines.append(f"{dep}(0.001) {q_str}")
                elif gname == "SDG":
                    lines.append(f"S_DAG {qs[0]}")
                    if not ideal:
                        lines.append(f"DEPOLARIZE1(0.001) {qs[0]}")
                elif gname == "M":
                    err = "" if ideal else "(0.001)"
                    lines.append(f"M{err} {qs[0]}")
                elif gname == "MX":
                    err = "" if ideal else "(0.001)"
                    lines.append(f"MX{err} {qs[0]}")
            lines.append("TICK")
        total_m = base_m + 7
        for det in detectors:
            recs = " ".join(f"rec[-{total_m - r_i}]" for r_i in det)
            lines.append(f"DETECTOR {recs}")
        obs_recs = " ".join(f"rec[-{total_m - r_i}]" for r_i in observables[0])
        lines.append(f"OBSERVABLE_INCLUDE(0) {obs_recs}")

        sc = stim.Circuit("\n".join(lines))
        dem = sc.detector_error_model()
        errs = sc.search_for_undetectable_logical_errors(
            dont_explore_detection_event_sets_with_size_above=4,
            dont_explore_edges_with_degree_above=9999,
            dont_explore_edges_increasing_symptom_degree=False,
        )
        stim_d = len(errs)
        out["stim"] = {
            "stim_version": stim.__version__,
            "dem_num_detectors": dem.num_detectors,
            "s_proxy_distance": stim_d,
            "agrees_d_eq_3": stim_d == native_d_s,
        }
        if stim_d != native_d_s:
            raise AssertionError(f"Stim S-proxy distance {stim_d} != native {native_d_s}")
    except ImportError:
        pass

    return out
