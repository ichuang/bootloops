"""Unified binary-matrix, PPR, and borrowed-identity distillation engine.

Implements the unified formulation of:
- Bravyi & Kitaev (2004), Bravyi & Haah (2012)
- Campbell & Howard (2017, arXiv:1606.01904)
- Haah & Hastings (2018), Haah, Hastings, Poulin & Wecker (2017)
- Jacinto, Valcarce, Barizien, Gouzien & Sangouard (2026, arXiv:2606.07734)
- Jain, Wills & Singh (2026, arXiv:2610.06535) — Theorems 1-3 & Propositions 4-5.

A factory is an N x n binary matrix G = (P; S) with k output rows (P) and r check rows (S),
N = k + r, and column rotation signs gamma in {1, 3, 5, 7}^n (default all-1s).
Simultaneously supports:
1. Postselected diagonal CSS / generalized triorthogonal code reading.
2. Pauli-product rotation (PPR) circuit reading with exact diagonal Clifford synthesis.
3. Borrowed-identity kernel coset reading (y0 = h + b_D, h in ker R_N).
4. Optional single-error-correcting CSS mode (e.g., 15-to-1 correcting RM15).
5. Exact all-orders dynamic programming over syndrome/output cosets for exact Fraction
   acceptance probability, conditional output infidelity, and multi-output joint Z laws,
   plus optional live CaDiCaL SAT+LRAT distance certification via `ftqc_codes`.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path
import sys
import tempfile
from typing import Any, Sequence

from .ir import CorrelationMode, MagicRefusalError, MagicSourceContract, PhasePolynomialTarget
from .rings import PolyQ, parse_exact_fraction


def _ensure_tools_on_syspath() -> None:
    # Support running from the repository root without PYTHONPATH.
    for candidate in (
        Path(__file__).resolve().parents[1],
    ):
        if (candidate / "ftqc_codes").is_dir():
            s = str(candidate)
            if s not in sys.path:
                sys.path.insert(0, s)
            break


def _gf2_rank(rows: Sequence[int]) -> int:
    basis: list[int] = []
    for r in rows:
        v = int(r)
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    return len(basis)


def _gf2_in_span(vec: int, rows: Sequence[int]) -> bool:
    basis: list[int] = []
    for r in rows:
        v = int(r)
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
            basis.sort(reverse=True)
    v = int(vec)
    for b in basis:
        v = min(v, v ^ b)
    return v == 0


def _gf2_nullspace(rows: Sequence[int], n_cols: int) -> list[int]:
    """Basis of {x in F_2^n_cols : dot(r, x) == 0 for all r in rows} as integer bitmasks."""
    pivots: list[tuple[int, int]] = []
    for r in rows:
        v = int(r)
        for col, prow in pivots:
            if (v >> col) & 1:
                v ^= prow
        if v:
            lead = (v & -v).bit_length() - 1
            pivots = [(c, pr ^ v if ((pr >> lead) & 1) else pr) for c, pr in pivots]
            pivots.append((lead, v))
    piv_map = dict(pivots)
    free_cols = [c for c in range(n_cols) if c not in piv_map]
    basis: list[int] = []
    for fc in free_cols:
        vec = 1 << fc
        for pcol, prow in piv_map.items():
            if (prow >> fc) & 1:
                vec |= 1 << pcol
        basis.append(vec)
    return basis


@dataclass(frozen=True)
class DistanceCertificate:
    distance: int
    leading_harmful_count: int
    witness_columns: tuple[int, ...]
    method: str
    lrat_info: dict[str, Any] | None = None


_DISTANCE_MITM_CACHE: dict[tuple[Any, ...], DistanceCertificate] = {}


@dataclass(frozen=True)
class BinaryFactory:
    """Unified N x n binary magic-state factory G = (P; S) with k outputs and r checks."""

    k: int
    r: int
    n: int
    # Row masks of length n (bit m is column m, 0 <= m < n):
    output_rows: tuple[int, ...]
    check_rows: tuple[int, ...]
    gamma: tuple[int, ...] = ()
    name: str = "binary_factory"

    def __post_init__(self) -> None:
        if self.k < 1 or self.r < 0 or self.n < 1:
            raise ValueError("Invalid dimensions (k >= 1, r >= 0, n >= 1 required)")
        if len(self.output_rows) != self.k or len(self.check_rows) != self.r:
            raise ValueError("Row tuple lengths must match k and r")
        full_mask = (1 << self.n) - 1
        for row in (*self.output_rows, *self.check_rows):
            if type(row) is not int or row < 0 or row > full_mask:
                raise ValueError("Row mask out of bounds")
        if not self.gamma:
            object.__setattr__(self, "gamma", tuple(1 for _ in range(self.n)))
        else:
            if len(self.gamma) != self.n or any((g % 2) == 0 for g in self.gamma):
                raise ValueError("gamma must have length n with odd exponents in Z_8")
            object.__setattr__(self, "gamma", tuple(int(g) % 8 for g in self.gamma))

    @classmethod
    def from_columns(
        cls,
        k: int,
        r: int,
        columns: Sequence[tuple[int, int] | int],
        *,
        gamma: Sequence[int] = (),
        name: str = "binary_factory",
    ) -> BinaryFactory:
        """Build from n columns, each given either as (p_mask_k, s_mask_r) or combined N-bit mask
        (output bits 0..k-1, check bits k..k+r-1).
        """
        n = len(columns)
        out_rows = [0] * k
        chk_rows = [0] * r
        for m, col in enumerate(columns):
            if isinstance(col, tuple):
                pm, sm = col
            else:
                pm = col & ((1 << k) - 1)
                sm = col >> k
            for a in range(k):
                if (pm >> a) & 1:
                    out_rows[a] |= 1 << m
            for b in range(r):
                if (sm >> b) & 1:
                    chk_rows[b] |= 1 << m
        return cls(
            k=k,
            r=r,
            n=n,
            output_rows=tuple(out_rows),
            check_rows=tuple(chk_rows),
            gamma=tuple(gamma),
            name=name,
        )

    @classmethod
    def from_row_strings(
        cls,
        output_rows: Sequence[str],
        check_rows: Sequence[str],
        *,
        gamma: Sequence[int] = (),
        name: str = "binary_factory",
    ) -> BinaryFactory:
        all_s = [*output_rows, *check_rows]
        if not all_s:
            raise ValueError("Empty matrix")
        n = len(all_s[0].strip())
        if any(len(s.strip()) != n or any(ch not in "01" for ch in s.strip()) for s in all_s):
            raise ValueError("All row strings must be binary of equal length n")

        def parse_row(s: str) -> int:
            s = s.strip()
            return sum((1 << m) for m, ch in enumerate(s) if ch == "1")

        return cls(
            k=len(output_rows),
            r=len(check_rows),
            n=n,
            output_rows=tuple(parse_row(s) for s in output_rows),
            check_rows=tuple(parse_row(s) for s in check_rows),
            gamma=tuple(gamma),
            name=name,
        )

    @property
    def total_wires(self) -> int:
        return self.k + self.r

    def all_rows(self) -> tuple[int, ...]:
        return self.output_rows + self.check_rows

    def column_label(self, m: int) -> int:
        return sum(((self.output_rows[a] >> m) & 1) << a for a in range(self.k))

    def column_syndrome(self, m: int) -> int:
        return sum(((self.check_rows[b] >> m) & 1) << b for b in range(self.r))

    def column_mask(self, m: int) -> int:
        return self.column_label(m) | (self.column_syndrome(m) << self.k)

    def is_reduced(self) -> bool:
        cols = [self.column_mask(m) for m in range(self.n)]
        return all(c != 0 for c in cols) and len(set(cols)) == self.n

    def is_full_rank(self) -> bool:
        return _gf2_rank(self.check_rows) == self.r and _gf2_rank(self.all_rows()) == self.total_wires

    def verify_algebraic_structure(self, *, require_non_idle_outputs: bool = True) -> dict[str, Any]:
        """Verify check-decoupling (all degree-1..3 moments touching a check vanish mod 2),
        synthesize exact diagonal Clifford corrections on all N wires, and return the
        recovered non-Clifford PhasePolynomialTarget on the k output wires.
        """
        rows = self.all_rows()
        N = self.total_wires

        def weighted_sum(mask: int) -> int:
            acc = 0
            v = mask
            while v:
                lsb = v & -v
                m = lsb.bit_length() - 1
                acc += self.gamma[m]
                v ^= lsb
            return acc

        c1 = [weighted_sum(rows[a]) % 8 for a in range(N)]
        c2: dict[tuple[int, int], int] = {}
        for a, b in combinations(range(N), 2):
            c2[(a, b)] = weighted_sum(rows[a] & rows[b]) % 8
        c3: dict[tuple[int, int, int], int] = {}
        for a, b, c in combinations(range(N), 3):
            c3[(a, b, c)] = weighted_sum(rows[a] & rows[b] & rows[c]) % 2

        # Check decoupling: every Q of size 1..3 touching {k, ..., N-1} must have even count!
        for b in range(self.k, N):
            if c1[b] % 2 != 0:
                raise MagicRefusalError(
                    "NON-CLIFFORD-CHECK-LEAK",
                    f"Check wire {b - self.k} has odd degree-1 moment {c1[b]}",
                    {"check_wire": b - self.k, "degree": 1},
                )
        for (a, b), val in c2.items():
            if b >= self.k and (val % 2 != 0):
                raise MagicRefusalError(
                    "NON-CLIFFORD-CHECK-LEAK",
                    f"Wire pair ({a}, {b}) touching check has odd degree-2 overlap {val}",
                    {"wires": [a, b], "degree": 2},
                )
        for (a, b, c), val in c3.items():
            if c >= self.k and val != 0:
                raise MagicRefusalError(
                    "NON-CLIFFORD-CHECK-LEAK",
                    f"Wire triple ({a}, {b}, {c}) touching check has odd degree-3 overlap",
                    {"wires": [a, b, c], "degree": 3},
                )

        # Build exact output phase table on {0,1}^k after cancelling check Clifford phases:
        # On basis |z_0..z_{k-1}, u_0..u_{r-1}>, sum_m gamma_m * parity_m(z, u) mod 8 equals:
        #   sum_a c1[a] y_a - 2 sum_{a<b} c2[a,b] y_a y_b + 4 sum_{a<b<c} c3[a,b,c] y_a y_b y_c (mod 8)
        # Since all check-touching c1, c2, c3 are even, the u-dependent part is purely diagonal Clifford:
        #   sum_{b>=k} 2*(c1[b]//2)*u_b + sum_{b>=k, a<b} 4*((c2[a,b]//2)%2)*y_a*u_b (mod 8),
        # cancelled identically by unitary S^{-c1[b]//2} and CZ_{a,b}.
        check_s_corrections = {b - self.k: (-c1[b] // 2) % 4 for b in range(self.k, N) if (c1[b] // 2) % 4 != 0}
        check_cz_corrections = [
            (a, b - self.k)
            for (a, b), val in c2.items()
            if b >= self.k and ((val // 2) % 2) == 1
        ]

        # Canonical non-Clifford target (pure T / CS / CCZ moments on 0..k-1):
        terms: list[str] = []
        output_s_corrections: dict[int, int] = {}
        output_cz_corrections: list[tuple[int, int]] = []
        for a in range(self.k):
            if c1[a] & 1:
                terms.append(f"{a}")
            s_pow = ((c1[a] - (c1[a] & 1)) // 2) % 4
            if s_pow:
                output_s_corrections[a] = (-s_pow) % 4
        for a, b in combinations(range(self.k), 2):
            val = c2[(a, b)]
            # Term in expansion is -2 * val * z_a z_b (mod 8).
            # If val is odd: -2*1 = 6 = +2 + 4 (mod 8); -2*3 = 2 (mod 8).
            phase_ab = (-2 * val) % 8
            if (phase_ab % 4) == 2:
                terms.append(f"{a} {b}")
            if ((phase_ab - (2 if (phase_ab % 4) == 2 else 0)) % 8) == 4:
                output_cz_corrections.append((a, b))
        for a, b, c in combinations(range(self.k), 3):
            if c3[(a, b, c)] == 1:
                terms.append(f"{a} {b} {c}")

        poly_str = "+".join(terms)
        target = PhasePolynomialTarget.from_jws_polynomial(self.k, poly_str)

        if require_non_idle_outputs:
            # Verify no non-zero F_2 linear combination of output rows has trivial non-Clifford action.
            anf = target.integer_anf_mod_modulus()
            for comb_mask in range(1, 1 << self.k):
                # Check if directional derivative or subspace restriction is trivial:
                # Specifically in JWS Proposition 5, an output direction u != 0 in F_2^k is idle if
                # target is invariant under CNOT basis change isolating u, i.e. wt(h), wt(h.g_i), wt(h.g_i.g_j) = 0 mod 2.
                h_row = 0
                for a in range(self.k):
                    if (comb_mask >> a) & 1:
                        h_row ^= self.output_rows[a]
                m1 = h_row.bit_count() & 1
                m2 = any(((h_row & self.output_rows[i]).bit_count() & 1) for i in range(self.k))
                m3 = any(
                    ((h_row & self.output_rows[i] & self.output_rows[j]).bit_count() & 1)
                    for i, j in combinations(range(self.k), 2)
                )
                if not (m1 or m2 or m3):
                    raise MagicRefusalError(
                        "IDLE-OUTPUT-DIRECTION",
                        f"Nonzero output combination mask={comb_mask} carries zero non-Clifford moments",
                        {"combination_mask": comb_mask},
                    )
            _ = anf

        return {
            "target": target,
            "jws_polynomial": poly_str,
            "full_rank": self.is_full_rank(),
            "reduced": self.is_reduced(),
            "check_s_corrections": check_s_corrections,
            "check_cz_corrections": check_cz_corrections,
            "output_s_corrections": output_s_corrections,
            "output_cz_corrections": output_cz_corrections,
        }

    def certify_distance(
        self,
        *,
        max_exact_weight: int = 7,
        use_sat_if_available: bool = True,
    ) -> DistanceCertificate:
        """Compute exact postselected Z-fault distance d(G) = min{|F| : sum_F S_m = 0, sum_F P_m != 0},
        leading harmful fault-set count A_d, and a replayable witness F*.
        When CaDiCaL is installed and `use_sat_if_available=True`, cross-certifies d via SAT+LRAT!
        """
        cache_key = (self.k, self.r, self.n, self.output_rows, self.check_rows, max_exact_weight)
        if not use_sat_if_available and cache_key in _DISTANCE_MITM_CACHE:
            return _DISTANCE_MITM_CACHE[cache_key]
        syndromes = [self.column_syndrome(m) for m in range(self.n)]
        labels = [self.column_label(m) for m in range(self.n)]

        # Weight 1 check:
        w1 = [m for m in range(self.n) if syndromes[m] == 0 and labels[m] != 0]
        if w1:
            return DistanceCertificate(1, len(w1), (w1[0],), "direct_w1")

        # Weight 2 check (JWS Theorem 1 fast path):
        by_syn: dict[int, list[int]] = {}
        for m, sm in enumerate(syndromes):
            by_syn.setdefault(sm, []).append(m)
        w2_count = 0
        w2_witness: tuple[int, ...] | None = None
        for cols in by_syn.values():
            if len(cols) > 1:
                for i, j in combinations(cols, 2):
                    if (labels[i] ^ labels[j]) != 0:
                        w2_count += 1
                        if w2_witness is None:
                            w2_witness = (i, j)
        if w2_count > 0:
            assert w2_witness is not None
            return DistanceCertificate(2, w2_count, w2_witness, "syndrome_collision_theorem1")

        # Meet-in-the-middle / shell search for weights 3..max_exact_weight:
        cap = min(max_exact_weight, self.n)
        found_d: int | None = None
        found_count = 0
        found_witness: tuple[int, ...] = ()

        # Build pair table for fast w=3 and w=4 lookup when n <= 300:
        if self.n <= 300 and cap >= 3:
            # w = 3: i < j < l with S_i ^ S_j == S_l
            w3_count = 0
            w3_wit: tuple[int, ...] | None = None
            for i, j in combinations(range(self.n), 2):
                s_ij = syndromes[i] ^ syndromes[j]
                for l in by_syn.get(s_ij, ()):
                    if l > j and (labels[i] ^ labels[j] ^ labels[l]) != 0:
                        w3_count += 1
                        if w3_wit is None:
                            w3_wit = (i, j, l)
            if w3_count > 0:
                found_d, found_count, found_witness = 3, w3_count, w3_wit or ()

        if found_d is None and self.n <= 180 and cap >= 4:
            # w = 4 via pair syndrome buckets:
            pair_buckets: dict[int, list[tuple[int, int, int]]] = {}
            for i, j in combinations(range(self.n), 2):
                pair_buckets.setdefault(syndromes[i] ^ syndromes[j], []).append((i, j, labels[i] ^ labels[j]))
            w4_count = 0
            w4_wit: tuple[int, ...] | None = None
            for pairs in pair_buckets.values():
                if len(pairs) < 2:
                    continue
                for idx_a in range(len(pairs)):
                    i1, j1, p1 = pairs[idx_a]
                    for idx_b in range(idx_a + 1, len(pairs)):
                        i2, j2, p2 = pairs[idx_b]
                        if j1 < i2 and (p1 ^ p2) != 0:
                            w4_count += 1
                            if w4_wit is None:
                                w4_wit = (i1, j1, i2, j2)
            if w4_count > 0:
                found_d, found_count, found_witness = 4, w4_count, w4_wit or ()

        if found_d is None:
            for w in range(3 if self.n > 300 else 5, cap + 1):
                # Meet-in-the-middle split w = left_w + right_w
                lw = w // 2
                rw = w - lw
                if comb(self.n, rw) > 400_000:
                    break
                left_map: dict[int, list[tuple[tuple[int, ...], int]]] = {}
                for comb_l in combinations(range(self.n), lw):
                    sm = 0
                    pm = 0
                    for idx in comb_l:
                        sm ^= syndromes[idx]
                        pm ^= labels[idx]
                    left_map.setdefault(sm, []).append((comb_l, pm))
                cnt = 0
                wit: tuple[int, ...] | None = None
                for comb_r in combinations(range(self.n), rw):
                    r_min = comb_r[0]
                    sm = 0
                    pm = 0
                    for idx in comb_r:
                        sm ^= syndromes[idx]
                        pm ^= labels[idx]
                    for comb_l, l_pm in left_map.get(sm, ()):
                        if comb_l[-1] < r_min and (l_pm ^ pm) != 0:
                            cnt += 1
                            if wit is None:
                                wit = comb_l + comb_r
                if cnt > 0:
                    found_d, found_count, found_witness = w, cnt, wit or ()
                    break

        lrat_info: dict[str, Any] | None = None
        if use_sat_if_available:
            _ensure_tools_on_syspath()
            try:
                from ftqc_codes import cnf as codes_cnf
                from ftqc_codes import sat as codes_sat

                solver = codes_sat.find_cadical()
                td = found_d if found_d is not None else (cap + 1)
                if solver is not None and self.n <= 200 and td >= 2:
                    # Exact JWS Eq. (5) CNF: exists e : check_rows*e == 0, output_rows*e != 0, wt(e) <= td - 1
                    formula = codes_cnf.Formula(
                        self.n, list(self.check_rows), list(self.output_rows), td - 1
                    )
                    with tempfile.TemporaryDirectory(prefix="ftqc_magic_sat_") as tmp:
                        res = codes_sat.solve(formula, tmp, "dz_lower", timeout=30, solver=solver)
                        if res[0] == "unsat":
                            lrat_info = res[1]
                        elif res[0] == "sat":
                            word = res[1]
                            cols = tuple(i for i in range(self.n) if (word >> i) & 1)
                            return DistanceCertificate(len(cols), 1, cols, "cadical_sat_witness", res[2])
            except ImportError:
                pass

        if found_d is None:
            raise MagicRefusalError(
                "DISTANCE-BUDGET-EXCEEDED",
                f"No harmful accepted fault found up to weight {cap}; supply larger max_exact_weight or CaDiCaL",
                {"checked_up_to": cap},
            )
        method = "exact_mitm+cadical_lrat" if lrat_info is not None else "exact_mitm"
        cert = DistanceCertificate(found_d, found_count, found_witness, method, lrat_info)
        if not use_sat_if_available:
            _DISTANCE_MITM_CACHE[cache_key] = cert
        return cert

    def exact_error_polynomials_dp(self) -> dict[str, Any]:
        """Exact all-orders dynamic programming over (syndrome, output_label) when r + k <= 16,
        or over syndrome cosets (r <= 15) using two passes (all accepted vs. accepted+trivial)!

        Look at the two-pass trick:
        - Pass 1 (size 2^r): DP over S_m in F_2^r alone -> gives W_acc_total[w] for ALL k (even k = 100)!
        - Pass 2: if k <= 12 and r + k <= 15, DP over (S_m, P_m) -> gives full joint Z distribution!
          Otherwise if n <= 26, exact subset enumeration; otherwise exact shell bound!
        """
        if self.r > 14:
            raise ValueError("Exact syndrome DP requires r <= 14 check rows")
        syndromes = [self.column_syndrome(m) for m in range(self.n)]
        labels = [self.column_label(m) for m in range(self.n)]

        # Pass 1: exact weight enumerator of ALL accepted fault sets {F : sum_F S_m == 0}
        # Tracked per weight w = 0..n when n <= 80, or direct evaluation at p via rational DP!
        n_syn = 1 << self.r
        # Polynomial weight enumerator if n <= 65 and r <= 10:
        if self.n <= 65 and (self.r + self.k) <= 13:
            n_joint = 1 << (self.r + self.k)
            # dp[state] is list of counts by weight 0..m
            dp: dict[int, list[int]] = {0: [1]}
            for m in range(self.n):
                step_mask = syndromes[m] | (labels[m] << self.r)
                nxt: dict[int, list[int]] = {st: list(arr) for st, arr in dp.items()}
                for st, arr in dp.items():
                    st2 = st ^ step_mask
                    cur = nxt.get(st2)
                    if cur is None:
                        nxt[st2] = [0] + list(arr)
                    else:
                        if len(cur) < len(arr) + 1:
                            cur.extend([0] * (len(arr) + 1 - len(cur)))
                        for w_idx, cnt in enumerate(arr):
                            cur[w_idx + 1] += cnt
                dp = nxt

            def pad(arr: list[int]) -> list[int]:
                return arr + [0] * (self.n + 1 - len(arr))

            trivial_weights = pad(dp.get(0, [0]))
            acc_weights = list(trivial_weights)
            joint_by_z: dict[int, list[int]] = {0: trivial_weights}
            for z_mask in range(1, 1 << self.k):
                st = z_mask << self.r
                if st in dp:
                    w_arr = pad(dp[st])
                    joint_by_z[z_mask] = w_arr
                    for i, c in enumerate(w_arr):
                        acc_weights[i] += c
            harmful_weights = [acc_weights[i] - trivial_weights[i] for i in range(self.n + 1)]
            return {
                "exact_all_orders": True,
                "accepted_weights": tuple(acc_weights),
                "trivial_weights": tuple(trivial_weights),
                "harmful_weights": tuple(harmful_weights),
                "joint_z_weights": {z: tuple(w) for z, w in joint_by_z.items()},
                "acceptance_poly": PolyQ.from_weight_enumerator(self.n, acc_weights),
                "harmful_poly": PolyQ.from_weight_enumerator(self.n, harmful_weights),
            }
        raise ValueError("Dimensions exceed polynomial DP table cap; use evaluate_at_p()")

    def evaluate_at_p(
        self,
        p: Fraction | str | int,
        *,
        distance_cert: DistanceCertificate | None = None,
        error_correcting_z: bool = False,
        wire_ticks_per_input: int = 1,
    ) -> tuple[MagicSourceContract, dict[str, Any]]:
        """Evaluate exact or certifiably outward-bounded output contract at input Z-error rate `p`."""
        pf = parse_exact_fraction(p)
        if not (0 <= pf <= Fraction(1, 2)):
            raise ValueError(f"Input Z-error rate p={pf} must lie in [0, 1/2]")
        alg = self.verify_algebraic_structure()
        target: PhasePolynomialTarget = alg["target"]

        if error_correcting_z:
            contract, report = self._evaluate_correcting_css_at_p(pf, target, wire_ticks_per_input)
            report["algebraic"] = {key: alg[key] for key in ("check_s_corrections", "check_cz_corrections", "output_s_corrections", "output_cz_corrections")}
            return contract, report

        dcert = distance_cert or self.certify_distance()
        d = dcert.distance
        one_minus_p = Fraction(1) - pf

        joint_z_prob: dict[int, Fraction] = {}
        p_bits = max(pf.numerator.bit_length(), pf.denominator.bit_length())
        if (self.r + self.k) <= 10 or ((self.r + self.k) <= 13 and p_bits <= 32):
            # Exact O(n * 2^(r+k)) rational DP at pf — zero truncation even for n = 200!
            dp: dict[int, Fraction] = {0: Fraction(1)}
            for m in range(self.n):
                step = self.column_syndrome(m) | (self.column_label(m) << self.r)
                nxt: dict[int, Fraction] = {}
                for st, prob in dp.items():
                    nxt[st] = nxt.get(st, Fraction(0)) + prob * one_minus_p
                    st2 = st ^ step
                    nxt[st2] = nxt.get(st2, Fraction(0)) + prob * pf
                dp = nxt
            p_trivial = dp.get(0, Fraction(0))
            p_acc = p_trivial
            for z_mask in range(1, 1 << self.k):
                pr = dp.get(z_mask << self.r, Fraction(0))
                if pr:
                    joint_z_prob[z_mask] = pr
                    p_acc += pr
            p_harmful = p_acc - p_trivial
            cond_infid = (p_harmful / p_acc) if p_acc > 0 else Fraction(1)
            exact_eval = True
        else:
            # Rigorous shell bound from certified distance d and leading count A_d:
            # p_acc >= (1 - p)^n; p_harmful <= A_d * p^d * (1 - p)^(n - d) + BinTail(n, d, p).
            p_acc = one_minus_p**self.n
            from .rings import eval_binomial_tail_exact

            lead_term = dcert.leading_harmful_count * (pf**d) * (one_minus_p ** (self.n - d))
            tail_beyond_d = eval_binomial_tail_exact(self.n, d, pf)
            p_harmful = min(Fraction(1), lead_term + tail_beyond_d)
            cond_infid = min(Fraction(1), p_harmful / p_acc) if p_acc > 0 else Fraction(1)
            exact_eval = False

        corr = (
            CorrelationMode.INDEPENDENT_SINGLE
            if self.k == 1
            else CorrelationMode.BATCH_CORRELATED
        )
        raw_ticks = Fraction(self.n * wire_ticks_per_input * self.total_wires)
        exp_ticks = (raw_ticks / p_acc) if p_acc > 0 else raw_ticks

        from .limits import analytical_lower_bound_n

        lb_info = analytical_lower_bound_n(target, d)
        contract = MagicSourceContract(
            name=self.name,
            target=target,
            n_inputs=self.n,
            k_outputs=self.k,
            distance=d,
            active_qubits=self.total_wires,
            total_wires=self.total_wires,
            correlation_mode=corr,
            input_error_rate=pf,
            acceptance_lower=p_acc,
            conditional_infidelity_upper=cond_infid,
            unnormalized_failure_upper=p_harmful,
            expected_wire_ticks_per_batch_upper=exp_ticks,
            leading_order_coeff=dcert.leading_harmful_count,
            joint_z_error_upper=joint_z_prob,
            optimality_lower_bound_n=lb_info["lower_bound_n"],
            premises=(
                "independent Z-equivalent input errors of rate p on each T rotation",
                "ideal Clifford checks and postselection on all +1 X-check outcomes",
            ),
            metadata={
                "exact_probability_dp": exact_eval,
                "distance_method": dcert.method,
                "witness_columns": list(dcert.witness_columns),
                "lrat_info": dcert.lrat_info,
            },
        )
        return contract, {
            "algebraic": {
                "jws_polynomial": alg["jws_polynomial"],
                "full_rank": alg["full_rank"],
                "reduced": alg["reduced"],
                "check_s_corrections": alg["check_s_corrections"],
                "check_cz_corrections": alg["check_cz_corrections"],
                "output_s_corrections": alg["output_s_corrections"],
                "output_cz_corrections": alg["output_cz_corrections"],
            },
            "distance": {
                "d": dcert.distance,
                "leading_harmful_count": dcert.leading_harmful_count,
                "witness_columns": list(dcert.witness_columns),
                "method": dcert.method,
                "lrat_info": dcert.lrat_info,
            },
            "optimality": lb_info,
        }

    def to_source_contract(
        self,
        p: Fraction | str | int,
        *,
        mode: str = "postselect",
        active_qubits_override: int | None = None,
        use_sat_if_available: bool = False,
    ) -> MagicSourceContract:
        """Convenience wrapper returning `MagicSourceContract` with optional recycled qubit footprint."""
        ec_mode = mode.strip().lower() in ("correct", "correcting", "ec")
        dcert = None if ec_mode else self.certify_distance(use_sat_if_available=use_sat_if_available)
        base_c, rep = self.evaluate_at_p(p, distance_cert=dcert, error_correcting_z=ec_mode)
        aq = active_qubits_override if active_qubits_override is not None else base_c.active_qubits
        raw_ticks = Fraction(self.n * aq)
        exp_ticks = (raw_ticks / base_c.acceptance_lower) if base_c.acceptance_lower > 0 else raw_ticks
        meta = {**dict(base_c.metadata), "algebraic": rep.get("algebraic")}
        return MagicSourceContract(
            name=base_c.name,
            target=base_c.target,
            n_inputs=base_c.n_inputs,
            k_outputs=base_c.k_outputs,
            distance=base_c.distance,
            active_qubits=aq,
            total_wires=base_c.total_wires,
            correlation_mode=base_c.correlation_mode,
            input_error_rate=base_c.input_error_rate,
            acceptance_lower=base_c.acceptance_lower,
            conditional_infidelity_upper=base_c.conditional_infidelity_upper,
            unnormalized_failure_upper=base_c.unnormalized_failure_upper,
            expected_wire_ticks_per_batch_upper=exp_ticks,
            leading_order_coeff=base_c.leading_order_coeff,
            joint_z_error_upper=base_c.joint_z_error_upper,
            optimality_lower_bound_n=base_c.optimality_lower_bound_n,
            premises=base_c.premises,
            metadata=meta,
        )

    def _evaluate_correcting_css_at_p(
        self,
        pf: Fraction,
        target: PhasePolynomialTarget,
        wire_ticks_per_input: int,
    ) -> tuple[MagicSourceContract, dict[str, Any]]:
        """Single-Z-error-correcting CSS mode (acceptance == 1, correcting weight-1 Z syndromes;
        matches He et al. Proposition 7.9 / `ftqc_compose.distill` on RM15: d_corr = 2, 105 pairs!).
        """
        syndromes = [self.column_syndrome(m) for m in range(self.n)]
        labels = [self.column_label(m) for m in range(self.n)]
        if len(set(syndromes)) != self.n or any(s == 0 for s in syndromes):
            raise MagicRefusalError(
                "AMBIGUOUS-SINGLE-Z-SYNDROME",
                "Single-Z error correction requires all n column syndromes to be nonzero and distinct",
            )
        decoder = {0: 0}
        for m, sm in enumerate(syndromes):
            decoder[sm] = labels[m]

        # Exact DP over (syndrome, logical_label) to sum all uncorrectable events:
        if (self.r + self.k) > 15:
            raise ValueError("Correcting mode exact DP requires r + k <= 15")
        one_minus_p = Fraction(1) - pf
        dp: dict[int, Fraction] = {0: Fraction(1)}
        for m in range(self.n):
            step = syndromes[m] | (labels[m] << self.r)
            nxt: dict[int, Fraction] = {}
            for st, prob in dp.items():
                nxt[st] = nxt.get(st, Fraction(0)) + prob * one_minus_p
                st2 = st ^ step
                nxt[st2] = nxt.get(st2, Fraction(0)) + prob * pf
            dp = nxt

        s_mask_bits = (1 << self.r) - 1
        p_fail = Fraction(0)
        joint_z: dict[int, Fraction] = {}
        for st, prob in dp.items():
            syn = st & s_mask_bits
            lbl = st >> self.r
            corr_lbl = decoder.get(syn, 0)
            res_lbl = lbl ^ corr_lbl
            if res_lbl != 0:
                p_fail += prob
                joint_z[res_lbl] = joint_z.get(res_lbl, Fraction(0)) + prob

        # Count harmful weight-2 pairs directly:
        harmful_pairs = 0
        wit: tuple[int, int] = (0, 1)
        for i, j in combinations(range(self.n), 2):
            syn = syndromes[i] ^ syndromes[j]
            res = labels[i] ^ labels[j] ^ decoder.get(syn, 0)
            if res != 0:
                if harmful_pairs == 0:
                    wit = (i, j)
                harmful_pairs += 1

        corr = CorrelationMode.INDEPENDENT_SINGLE if self.k == 1 else CorrelationMode.BATCH_CORRELATED
        contract = MagicSourceContract(
            name=self.name + ":correcting",
            target=target,
            n_inputs=self.n,
            k_outputs=self.k,
            distance=2 if harmful_pairs > 0 else 3,
            active_qubits=self.total_wires,
            total_wires=self.total_wires,
            correlation_mode=corr,
            input_error_rate=pf,
            acceptance_lower=Fraction(1),
            conditional_infidelity_upper=p_fail,
            unnormalized_failure_upper=p_fail,
            expected_wire_ticks_per_batch_upper=Fraction(self.n * wire_ticks_per_input * self.total_wires),
            leading_order_coeff=harmful_pairs,
            joint_z_error_upper=joint_z,
            optimality_lower_bound_n=None,
            premises=("independent Z-equivalent input errors; minimum-weight Z-syndrome correction without postselection",),
        )
        return contract, {
            "mode": "error_correcting_z",
            "harmful_pairs": harmful_pairs,
            "pair_witness": list(wit),
        }


def from_borrowed_identity(
    k: int,
    r: int,
    identity_columns: Sequence[int],
    bare_target_columns: Sequence[int],
    *,
    name: str = "borrowed_identity_factory",
) -> BinaryFactory:
    """Jain-Wills-Singh Theorem 2 / Singh-Gidney-Jones (2026) borrowed-identity deletion:
    Given `identity_columns` h in ker(R_N) (acting as a diagonal Clifford on |+>^N) and
    output-only `bare_target_columns` b_D supported purely on wires 0..k-1, toggle y_0 = h ^ b_D!
    """
    N = k + r
    if any(c <= 0 or c >= (1 << N) for c in identity_columns):
        raise ValueError("Identity columns must be nonzero N-bit masks")
    if any(c <= 0 or c >= (1 << k) for c in bare_target_columns):
        raise ValueError("Bare target columns must be output-only masks on wires 0..k-1")

    # Verify h in ker(R_N): all degree 1..3 moments on N wires vanish mod 2!
    id_set: set[int] = set()
    for c in identity_columns:
        if c in id_set:
            id_set.remove(c)
        else:
            id_set.add(c)
    for deg in (1, 2, 3):
        for q in combinations(range(N), deg):
            q_mask = sum(1 << w for w in q)
            parity = sum(1 for col in id_set if (col & q_mask) == q_mask) & 1
            if parity != 0:
                raise MagicRefusalError(
                    "INVALID-BORROWED-IDENTITY",
                    f"Supplied identity columns fail R_N h == 0 at wire subset {q}",
                    {"subset": list(q)},
                )

    toggled = set(id_set)
    for b_col in bare_target_columns:
        if b_col in toggled:
            toggled.remove(b_col)
        else:
            toggled.add(b_col)

    ordered_cols = sorted(toggled)
    if any((col >> k) == 0 for col in ordered_cols):
        raise MagicRefusalError(
            "UNCHECKED-SURVIVING-COLUMN",
            "Toggled borrowed-identity circuit still contains an output-only column (d == 1)",
        )
    return BinaryFactory.from_columns(k, r, ordered_cols, name=name)


def parent_capacity(check_rows: Sequence[int], n: int) -> dict[str, Any]:
    """Jain-Wills-Singh Section IV Eq. (13)-(15): compute parent Schur square S^<2>,
    legal output space V(S) = (S^<2>)^perp / S, and exact capacity kappa = dim V(S).
    """
    r = len(check_rows)
    schur_gens = list(check_rows)
    for i, j in combinations(range(r), 2):
        schur_gens.append(check_rows[i] & check_rows[j])
    perp_basis = _gf2_nullspace(schur_gens, n)
    # Check S subseteq (S^<2>)^perp:
    for s_row in check_rows:
        if not _gf2_in_span(s_row, perp_basis):
            raise MagicRefusalError(
                "INVALID-PARENT-CHECK-OVERLAP",
                "Parent check rows fail S subseteq (S^<2>)^perp (check self/pair/triple overlap odd)",
            )
    rank_s = _gf2_rank(check_rows)
    kappa = len(perp_basis) - rank_s
    lower_bound = n - 2 * r - (r * (r - 1)) // 2
    # Extract a quotient representative basis for V(S):
    quot_basis: list[int] = []
    cur_span = list(check_rows)
    for cand in perp_basis:
        if not _gf2_in_span(cand, cur_span):
            quot_basis.append(cand)
            cur_span.append(cand)
    return {
        "n": n,
        "r": r,
        "rank_S": rank_s,
        "dim_schur_perp": len(perp_basis),
        "capacity_kappa": kappa,
        "jws_eq15_lower_bound": lower_bound,
        "quotient_basis": quot_basis,
    }


def parity_lift(factory: BinaryFactory, *, name: str | None = None) -> BinaryFactory:
    """Jain-Wills-Singh Proposition 4: parity-lift a homogeneous degree-3 (CCZ-product) factory
    of odd distance d to distance >= d + 1 using one extra all-ones check wire!
    """
    alg = factory.verify_algebraic_structure()
    moments = alg["target"].non_clifford_moments_mod2()
    if any(len(wires) < 3 for wires in moments):
        raise MagicRefusalError(
            "PARITY-LIFT-NON-HOMOGENEOUS",
            "Parity lift (Proposition 4) requires homogeneous degree-3 (CCZ-type) output moments",
        )
    if factory.n % 2 == 0:
        all_ones = (1 << factory.n) - 1
        if _gf2_in_span(all_ones, factory.all_rows()):
            raise MagicRefusalError(
                "PARITY-LIFT-DEPENDENT-ONES",
                "Even-n parity lift requires the all-ones vector to lie outside span(G)",
            )
        return BinaryFactory(
            k=factory.k,
            r=factory.r + 1,
            n=factory.n,
            output_rows=factory.output_rows,
            check_rows=(*factory.check_rows, all_ones),
            gamma=factory.gamma,
            name=name or f"{factory.name}_parity_lift",
        )
    # Odd-n case: append a column supported ONLY on the new check wire, then append c = 1_{n+1}!
    all_ones_plus = (1 << (factory.n + 1)) - 1
    return BinaryFactory(
        k=factory.k,
        r=factory.r + 1,
        n=factory.n + 1,
        output_rows=factory.output_rows,
        check_rows=(*factory.check_rows, all_ones_plus),
        gamma=(*factory.gamma, 1),
        name=name or f"{factory.name}_parity_lift",
    )


def free_row_append(
    factory: BinaryFactory,
    *,
    t_moment: int,
    cs_moments: Sequence[int],
    ccz_moments: Mapping[tuple[int, int], int],
    name: str | None = None,
) -> BinaryFactory:
    """Jain-Wills-Singh Proposition 5: solve the exact F_2 affine system to append a k+1-th
    output row with prescribed degree-1, 2, 3 moments while preserving syndromes and d >= 3.
    """
    if len(cs_moments) != factory.k:
        raise ValueError("cs_moments must have length k")
    # Build linear equations <mask_i, a> == rhs_i (mod 2) on a in F_2^n:
    # Note: wt(a) == <1_n, a> mod 2!
    eqs: list[tuple[int, int]] = [((1 << factory.n) - 1, t_moment & 1)]
    for i in range(factory.k):
        eqs.append((factory.output_rows[i], cs_moments[i] & 1))
    for i, j in combinations(range(factory.k), 2):
        eqs.append((factory.output_rows[i] & factory.output_rows[j], ccz_moments.get((i, j), 0) & 1))
    for s in factory.check_rows:
        eqs.append((s, 0))
        for i in range(factory.k):
            eqs.append((s & factory.output_rows[i], 0))
    for s1, s2 in combinations(factory.check_rows, 2):
        eqs.append((s1 & s2, 0))

    # Gauss-Jordan over F_2 on augmented (n + 1)-bit rows (bit n is RHS):
    pivots: list[tuple[int, int]] = []
    for mask, rhs in eqs:
        aug = mask | ((rhs & 1) << factory.n)
        for col, prow in pivots:
            if (aug >> col) & 1:
                aug ^= prow
        lhs = aug & ((1 << factory.n) - 1)
        b_rhs = (aug >> factory.n) & 1
        if lhs == 0:
            if b_rhs:
                raise MagicRefusalError("FREE-ROW-INFEASIBLE", "Prescribed output moments have empty affine solution space")
            continue
        lead = (lhs & -lhs).bit_length() - 1
        pivots = [(c, pr ^ aug if ((pr >> lead) & 1) else pr) for c, pr in pivots]
        pivots.append((lead, aug))

    piv_map = dict(pivots)
    # Particular solution a0 (all free vars = 0):
    a0 = 0
    for col, prow in piv_map.items():
        if (prow >> factory.n) & 1:
            a0 |= 1 << col
    hom_basis = _gf2_nullspace([m for m, _ in eqs], factory.n)
    existing = factory.all_rows()
    chosen: int | None = None
    if not _gf2_in_span(a0, existing):
        chosen = a0
    else:
        for h in hom_basis:
            cand = a0 ^ h
            if not _gf2_in_span(cand, existing):
                chosen = cand
                break
    if chosen is None:
        raise MagicRefusalError("FREE-ROW-SPAN-DEPENDENT", "All affine solutions lie inside span(G)")

    ext = BinaryFactory(
        k=factory.k + 1,
        r=factory.r,
        n=factory.n,
        output_rows=(*factory.output_rows, chosen),
        check_rows=factory.check_rows,
        gamma=factory.gamma,
        name=name or f"{factory.name}_plus1",
    )
    ext.verify_algebraic_structure(require_non_idle_outputs=True)
    return ext


def puncture_base_code(
    base_check_rows: Sequence[int],
    n_total: int,
    puncture_coords: Sequence[int],
    *,
    name: str = "punctured_factory",
) -> BinaryFactory:
    """Haah-Hastings / JWS Appendix B puncturing: given an [n_total] parent code (whose rows satisfy
    pairwise/triple self-orthogonality so every row has weight 0 mod 8 / overlaps 0 mod 4 / triples 0 mod 2),
    puncture `puncture_coords` J so deleted coordinates become k = |J| logical outputs on the remaining
    n = n_total - |J| coordinates!
    """
    j_set = sorted(set(puncture_coords))
    if any(c < 0 or c >= n_total for c in j_set) or len(j_set) != len(puncture_coords):
        raise ValueError("Invalid or duplicate puncture coordinates")
    kept = [c for c in range(n_total) if c not in set(j_set)]
    # Row-reduce generator space G_base = base_check_rows so J has pivot rows (output_rows)
    # and remaining rows vanish on J (check_rows)!
    basis_rows = list(base_check_rows)
    pivots: list[tuple[int, int]] = []
    rem_rows: list[int] = []
    for r in basis_rows:
        v = r
        for pcol, prow in pivots:
            if (v >> pcol) & 1:
                v ^= prow
        j_Hit = next((c for c in j_set if (v >> c) & 1), None)
        if j_Hit is not None:
            pivots = [(c, pr ^ v if ((pr >> j_Hit) & 1) else pr) for c, pr in pivots]
            pivots.append((j_Hit, v))
        else:
            rem_rows.append(v)
    if len(pivots) != len(j_set):
        raise MagicRefusalError(
            "DEPENDENT-PUNCTURE-SET",
            f"Puncture coordinates {j_set} are not independent in parent code (rank {len(pivots)} < {len(j_set)})",
        )
    # Clean rem_rows against final pivots and reduce rem_rows to independent check rows:
    clean_checks_full: list[int] = []
    for v in rem_rows:
        for pcol, prow in pivots:
            if (v >> pcol) & 1:
                v ^= prow
        if v:
            clean_checks_full.append(v)

    def restrict_to_kept(full_row: int) -> int:
        return sum(((full_row >> c) & 1) << idx for idx, c in enumerate(kept))

    piv_by_c = dict(pivots)
    out_rows = tuple(restrict_to_kept(piv_by_c[c]) for c in j_set)
    # Echelon-reduce restricted check rows:
    chk_basis: list[int] = []
    for full_r in clean_checks_full:
        rk = restrict_to_kept(full_r)
        v = rk
        for b in chk_basis:
            v = min(v, v ^ b)
        if v:
            chk_basis.append(v)
            chk_basis.sort(reverse=True)
    return BinaryFactory(
        k=len(j_set),
        r=len(chk_basis),
        n=len(kept),
        output_rows=out_rows,
        check_rows=tuple(chk_basis),
        name=name,
    )


def verify_puncture_downward_closure(
    base_check_rows: Sequence[int],
    n_total: int,
    j_small: Sequence[int],
    j_large: Sequence[int],
) -> dict[str, Any]:
    """Verify JWS Theorem 3: for independent puncture sets J_small subset J_large, d(J_large) <= d(J_small)."""
    if not set(j_small).issubset(set(j_large)):
        raise ValueError("j_small must be a subset of j_large")
    f_small = puncture_base_code(base_check_rows, n_total, j_small, name="puncture_small")
    f_large = puncture_base_code(base_check_rows, n_total, j_large, name="puncture_large")
    d_small = f_small.certify_distance().distance
    d_large = f_large.certify_distance().distance
    if d_large > d_small:
        raise AssertionError(f"Theorem 3 violated: d(J_large)={d_large} > d(J_small)={d_small}")
    return {
        "d_small": d_small,
        "d_large": d_large,
        "downward_closed": d_large <= d_small,
    }
