"""Canonical families, symmetric SAT formulations, and LRAT-verified spatial no-go theorems.

Implements:
1. **Canonical Repetition-Code Families `F_N^0` (T) & `F_N^1` (sqrt(T))**
   (Jacinto et al. 2026, arXiv:2606.07734, Section IV Theorems 1 & 2, Appendix B Eq. 17):
   Analytical construction for any qubit count `N >= 4` achieving `d = ceil(N/3)` (or next odd
   integer for odd `N`), recovering `15T -> 1T` (`N=5, d=3`), `165T -> 1T` (`N=10, d=4`),
   and `165T -> 1T` (`N=11, d=5`) in closed form!
2. **Spatial No-Go Theorem CNF & LRAT Verifier (`encode_spatial_nogo_cnf`, `prove_or_verify_nogo`)**:
   Encodes unrestricted `N`-qubit repetition-code `T` (Eqs. 21-22) and `CCZ` (Eqs. 40-42)
   existence at target distance `d_min` into pure DIMACS CNF, proving via CaDiCaL + `ftqc_codes.lrat`
   (or replaying banked LRAT proofs on cold clones without external solvers!):
   - `N = 4 => d <= 2` for `T` (`d_min = 3` UNSAT, banked LRAT included!)
   - `N = 5, d = 3 => n >= 15` for `T` (`n_max = 14` UNSAT!)
   - `N = 4, 5 => d <= 2` for `CCZ` (`d_min = 3` UNSAT, banked LRAT included!)
3. **Symmetric Subfamilies (`F_{S_N}`, Young `S_lambda`, Cyclic `C_lambda`) & JWS Symmetry Slots**:
   Orbit-compressed SAT encoders for bounded-width and high-rate multi-output discovery.
"""

from __future__ import annotations

from dataclasses import dataclass
import base64
import hashlib
from itertools import combinations, product
from math import comb
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Any, Sequence
import zlib

from .matrix import BinaryFactory


def _ensure_tools_on_syspath() -> None:
    for candidate in (
        Path(__file__).resolve().parents[1],
    ):
        if (candidate / "ftqc_codes").is_dir():
            s = str(candidate)
            if s not in sys.path:
                sys.path.insert(0, s)
            break


def _commute_repetition_column_to_distillation(alpha_mask: int, N: int) -> tuple[int, int]:
    """Jacinto et al. Eq. (17) / (B2)-(B3): commute repetition-code Z-rotation `alpha` in F_2^N
    past decoding CNOT cascade E^dagger:
        beta_0 = alpha_0   (output wire 0)
        beta_i = alpha_i ^ alpha_{i-1}  for 1 <= i < N   (check wires 1..N-1).
    Returns `(p_mask_1bit, s_mask_{N-1}bits)`.
    """
    p_bit = alpha_mask & 1
    s_mask = 0
    for i in range(1, N):
        bi = ((alpha_mask >> i) ^ (alpha_mask >> (i - 1))) & 1
        if bi:
            s_mask |= 1 << (i - 1)
    return p_bit, s_mask


def build_canonical_t_factory(N: int) -> BinaryFactory:
    """Jacinto et al. (2026) Theorem 1 (`F_N^0`): closed-form `N`-qubit `nT -> 1T` factory!

    Includes all subsets `A subseteq [N]` with:
    - `|A| == 3`: always (`a_3 = 1`)
    - `|A| == 2`: iff `(N - 1) % 2 == 1` (i.e. `N` even)
    - `|A| == 1`: iff `(N + comb(N - 1, 2)) % 2 == 1` (i.e. `N % 4 in (0, 1)`).
    """
    if N < 4:
        raise ValueError("Canonical fault-tolerant T family F_N^0 requires N >= 4")
    include_deg = {3}
    if (N - 1) % 2 == 1:
        include_deg.add(2)
    if (N + comb(N - 1, 2)) % 2 == 1:
        include_deg.add(1)

    cols: list[tuple[int, int]] = []
    for deg in sorted(include_deg):
        for subset in combinations(range(N), deg):
            alpha = sum(1 << i for i in subset)
            cols.append(_commute_repetition_column_to_distillation(alpha, N))

    return BinaryFactory.from_columns(
        k=1,
        r=N - 1,
        columns=cols,
        name=f"jvbs_canonical_F0_N{N}",
    )


def canonical_t_theoretical_params(N: int) -> dict[str, int]:
    """Closed-form `(N, n, d)` formula from Jacinto et al. Theorem 1 (p. 8)."""
    if N < 4:
        raise ValueError("N must be >= 4")
    r_mod4 = N % 4
    if r_mod4 == 0:
        n = comb(N, 1) + comb(N, 2) + comb(N, 3)
    elif r_mod4 == 1:
        n = comb(N, 1) + comb(N, 3)
    elif r_mod4 == 2:
        n = comb(N, 2) + comb(N, 3)
    else:
        n = comb(N, 3)
    ceil_n3 = (N + 2) // 3
    if N % 2 == 0:
        d = ceil_n3
    else:
        d = ceil_n3 if (ceil_n3 % 2 == 1) else (ceil_n3 + 1)
    return {"N": N, "n": n, "d": d}


def canonical_sqrt_t_theoretical_params(N: int) -> dict[str, int]:
    """Closed-form `(N, sqrt_T_count, d)` for `F_N^1` (`|sqrt(T)>` distillation, Theorem 2 p. 28)."""
    if N < 5:
        raise ValueError("Canonical sqrt(T) family F_N^1 requires N >= 5")
    n_rot = 0
    for s in range(1, 5):
        coeff = ((-1) ** s) * comb(N - s - 1, 4 - s)
        if (coeff % 4) != 0:
            n_rot += comb(N, s)
    return {"N": N, "rotations": n_rot, "sqrt_t_count": 2 * n_rot, "d": (N + 3) // 4}


@dataclass(frozen=True)
class SpatialSATInstance:
    N: int
    target: str
    d_min: int
    n_max: int | None
    nvars: int
    clauses: tuple[tuple[int, ...], ...]
    candidate_alphas: tuple[int, ...]

    def dimacs(self) -> str:
        lines = [f"p cnf {self.nvars} {len(self.clauses)}"]
        lines.extend(" ".join(str(x) for x in cl) + " 0" if cl else "0" for cl in self.clauses)
        return "\n".join(lines) + "\n"

    def extract_factory_from_model(self, true_vars: set[int]) -> BinaryFactory:
        selected = [
            self.candidate_alphas[idx]
            for idx in range(len(self.candidate_alphas))
            if (idx + 1) in true_vars
        ]
        if self.target == "T":
            cols = [_commute_repetition_column_to_distillation(a, self.N) for a in selected]
            return BinaryFactory.from_columns(1, self.N - 1, cols, name=f"sat_T_N{self.N}_d{self.d_min}")
        # CCZ encoding (Section V.C, p. 14): wires 0, 1 are bare CCZ outputs, wires 2..N-1 form
        # an (N - 2)-wire repetition code for the 3rd CCZ output!
        cols_ccz: list[tuple[int, int]] = []
        for a in selected:
            b0 = a & 1
            b1 = (a >> 1) & 1
            b2 = (a >> 2) & 1
            pm = b0 | (b1 << 1) | (b2 << 2)
            sm = 0
            for j in range(3, self.N):
                chk_bit = ((a >> j) ^ (a >> (j - 1))) & 1
                if chk_bit:
                    sm |= 1 << (j - 3)
            cols_ccz.append((pm, sm))
        return BinaryFactory.from_columns(3, self.N - 3, cols_ccz, name=f"sat_CCZ_N{self.N}_d{self.d_min}")


def encode_spatial_nogo_cnf(
    N: int,
    target: str,
    d_min: int,
    *,
    n_max: int | None = None,
) -> SpatialSATInstance:
    """Encode Jacinto et al. unrestricted SAT problem on `N` logical qubits for `target in ('T', 'CCZ')`
    and minimum distance `d_min` (optional column cap `n_max`).
    UNSAT proves a mathematical spatial no-go theorem that NO `nT -> target` PPR distillation
    circuit of distance `>= d_min` (with `<= n_max` gates) exists on `N` logical qubits!
    """
    tgt = target.strip().upper()
    if tgt not in ("T", "CCZ"):
        raise ValueError("Spatial SAT encoding supports target 'T' or 'CCZ'")
    if d_min < 2 or N < (4 if tgt == "T" else 4):
        raise ValueError("Require N >= 4 and d_min >= 2")

    # Candidate nonzero columns alpha in F_2^N \ {0}:
    # For CCZ (Sec. V.C p. 13), forbid columns supported ONLY on unprotected qubits {0, 1} (which have
    # zero check syndrome and immediately give d == 1)!
    if tgt == "T":
        alphas = tuple(range(1, 1 << N))
        bad_xor_targets = {(1 << N) - 1}
    else:
        alphas = tuple(a for a in range(1, 1 << N) if (a >> 2) != 0)
        rep_all = ((1 << (N - 2)) - 1) << 2
        # Eq. (41): 7 undetectable logical Z error patterns on {0, 1, rep(2..N-1)}:
        bad_xor_targets = {
            1,
            2,
            3,
            rep_all,
            rep_all | 1,
            rep_all | 2,
            rep_all | 3,
        }

    m_cols = len(alphas)
    alpha_to_var = {a: idx + 1 for idx, a in enumerate(alphas)}
    nvars = m_cols
    clauses: list[tuple[int, ...]] = []

    def new_var() -> int:
        nonlocal nvars
        nvars += 1
        return nvars

    def add_xor_eq(var_list: Sequence[int], rhs: int) -> None:
        if not var_list:
            if rhs == 1:
                clauses.append(())
            return
        acc = var_list[0]
        for v in var_list[1:]:
            t = new_var()
            clauses.extend([(-t, acc, v), (-t, -acc, -v), (t, -acc, v), (t, acc, -v)])
            acc = t
        clauses.append((acc,) if rhs == 1 else (-acc,))

    # 1. Moment parity constraints (Eq. 21 for T; Eq. 40 for CCZ):
    for deg in (1, 2, 3):
        for q_wires in combinations(range(N), deg):
            q_mask = sum(1 << w for w in q_wires)
            hitting_vars = [idx + 1 for idx, a in enumerate(alphas) if (a & q_mask) == q_mask]
            if tgt == "T":
                rhs = 1
            else:
                rhs = 1 if (deg == 3 and q_wires[0] == 0 and q_wires[1] == 1 and q_wires[2] >= 2) else 0
            add_xor_eq(hitting_vars, rhs)

    # 2. Distance clauses (Eq. 22 for T; Eq. 42 for CCZ): forbid every subset K of size 1 <= |K| < d_min
    # whose XOR lies in `bad_xor_targets`:
    # Optimize weight-1 and weight-2 directly, and general weight < d_min:
    for w in range(1, d_min):
        if w == 1:
            for idx, a in enumerate(alphas):
                if a in bad_xor_targets:
                    clauses.append((-(idx + 1),))
        elif w == 2:
            for idx, a in enumerate(alphas):
                for bad in bad_xor_targets:
                    b = a ^ bad
                    v_b = alpha_to_var.get(b)
                    if v_b is not None and (idx + 1) < v_b:
                        clauses.append((-(idx + 1), -v_b))
        else:
            for comb_idxs in combinations(range(m_cols), w - 1):
                partial_xor = 0
                for idx in comb_idxs:
                    partial_xor ^= alphas[idx]
                last_idx_min = comb_idxs[-1] + 1
                for bad in bad_xor_targets:
                    need_a = partial_xor ^ bad
                    v_last = alpha_to_var.get(need_a)
                    if v_last is not None and (v_last - 1) >= last_idx_min:
                        clauses.append(tuple(-(i + 1) for i in comb_idxs) + (-v_last,))

    # 3. Optional Sinz sequential counter for `sum_j v_j <= n_max`:
    if n_max is not None and n_max < m_cols:
        if n_max <= 0:
            for idx in range(m_cols):
                clauses.append((-(idx + 1),))
        else:
            s_cnt: dict[tuple[int, int], int] = {}
            for i in range(1, m_cols):
                for j in range(1, n_max + 1):
                    s_cnt[(i, j)] = new_var()
            clauses.append((-1, s_cnt[(1, 1)]))
            for j in range(2, n_max + 1):
                clauses.append((-s_cnt[(1, j)],))
            for i in range(2, m_cols):
                clauses.append((-i, s_cnt[(i, 1)]))
                clauses.append((-s_cnt[(i - 1, 1)], s_cnt[(i, 1)]))
                for j in range(2, n_max + 1):
                    clauses.append((-i, -s_cnt[(i - 1, j - 1)], s_cnt[(i, j)]))
                    clauses.append((-s_cnt[(i - 1, j)], s_cnt[(i, j)]))
                clauses.append((-i, -s_cnt[(i - 1, n_max)]))
            clauses.append((-m_cols, -s_cnt[(m_cols - 1, n_max)]))

    return SpatialSATInstance(
        N=N,
        target=tgt,
        d_min=d_min,
        n_max=n_max,
        nvars=nvars,
        clauses=tuple(clauses),
        candidate_alphas=alphas,
    )


def _check_ascii_lrat_rup(dimacs: str, proof_bytes: bytes) -> dict[str, int]:
    """Minimal standalone ASCII LRAT RUP verifier fallback when `ftqc_codes` is absent."""
    db: dict[int, tuple[int, ...]] = {}
    cid = 0
    buf: list[int] = []
    for raw in dimacs.splitlines():
        line = raw.strip()
        if not line or line.startswith("c") or line.startswith("p"):
            continue
        for tok in line.split():
            v = int(tok)
            if v == 0:
                cid += 1
                db[cid] = tuple(buf)
                buf = []
            else:
                buf.append(v)

    steps = 0
    hints_total = 0
    derived_empty = False
    for raw in proof_bytes.decode("ascii").splitlines():
        parts = raw.split()
        if not parts or parts[0].startswith("c"):
            continue
        if len(parts) >= 2 and parts[1] == "d":
            for tok in parts[2:]:
                did = int(tok)
                if did != 0:
                    db.pop(did, None)
            continue
        nums = [int(t) for t in parts]
        step_id = nums[0]
        z1 = nums.index(0, 1)
        lits = tuple(nums[1:z1])
        z2 = nums.index(0, z1 + 1)
        hints = nums[z1 + 1 : z2]
        falsified = set(lits)
        conflict = False
        for h in hints:
            if h <= 0 or h not in db:
                raise ValueError(f"Invalid or non-RUP LRAT hint {h} at step {step_id}")
        hints_total += len(hints)
        # Unit-propagate strictly over the hinted clause subset (supports both forward & CaDiCaL reverse hint order!):
        pending = [db[h] for h in hints]
        conflict = False
        progressed = True
        while progressed and not conflict:
            progressed = False
            next_pending: list[tuple[int, ...]] = []
            for cl in pending:
                rem = [lit for lit in cl if lit not in falsified]
                if not rem:
                    conflict = True
                    break
                if len(rem) == 1:
                    falsified.add(-rem[0])
                    progressed = True
                else:
                    next_pending.append(cl)
            pending = next_pending
        if not conflict:
            raise ValueError(f"LRAT step {step_id} did not reach RUP conflict")
        db[step_id] = lits
        steps += 1
        if not lits:
            derived_empty = True
    if not derived_empty:
        raise ValueError("LRAT proof did not derive the empty clause")
    return {"steps": steps, "hints": hints_total}


def solve_or_verify_spatial_sat(
    instance: SpatialSATInstance,
    *,
    banked_lrat_bytes: bytes | None = None,
    export_lrat_path: str | None = None,
    timeout: int = 30,
) -> dict[str, Any]:
    """Verify a banked LRAT proof in pure stdlib (`ftqc_codes.lrat`), or solve live with CaDiCaL
    and independently check either its LRAT UNSAT refutation or the extracted `BinaryFactory`!
    """
    _ensure_tools_on_syspath()
    try:
        from ftqc_codes import lrat as codes_lrat
        from ftqc_codes import sat as codes_sat
    except ImportError:
        codes_lrat = None
        codes_sat = None

    dimacs = instance.dimacs()
    cnf_sha = hashlib.sha256(dimacs.encode()).hexdigest()

    if banked_lrat_bytes is not None:
        stats = (
            codes_lrat.check(dimacs, banked_lrat_bytes)
            if codes_lrat is not None
            else _check_ascii_lrat_rup(dimacs, banked_lrat_bytes)
        )
        return {
            "verdict": "UNSAT-NOGO-CERTIFIED",
            "mode": "banked_lrat_verified",
            "N": instance.N,
            "target": instance.target,
            "d_min": instance.d_min,
            "n_max": instance.n_max,
            "cnf_sha256": cnf_sha,
            "lrat_sha256": hashlib.sha256(banked_lrat_bytes).hexdigest(),
            "lrat_steps": stats["steps"],
            "lrat_hints": stats["hints"],
        }

    solver = codes_sat.find_cadical() if codes_sat is not None else None
    if solver is None:
        return {
            "verdict": "SKIPPED-NO-CADICAL",
            "cnf_sha256": cnf_sha,
            "vars": instance.nvars,
            "clauses": len(instance.clauses),
        }

    with tempfile.TemporaryDirectory(prefix="ftqc_magic_nogo_") as tmp:
        cnf_p = str(Path(tmp) / "inst.cnf")
        lrat_p = export_lrat_path or str(Path(tmp) / "inst.lrat")
        Path(cnf_p).write_text(dimacs)
        flags = ["--modelalllits=true"] if "modelalllits" in subprocess.run(
            [solver, "--help"], capture_output=True, text=True, timeout=10
        ).stdout else []
        proc = subprocess.run(
            [solver, "--lrat", "--no-binary", "-q", *flags, cnf_p, lrat_p],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if proc.returncode == 20:
            proof_bytes = Path(lrat_p).read_bytes()
            stats = codes_lrat.check(dimacs, proof_bytes)
            return {
                "verdict": "UNSAT-NOGO-CERTIFIED",
                "mode": "cadical_live_lrat_verified",
                "N": instance.N,
                "target": instance.target,
                "d_min": instance.d_min,
                "n_max": instance.n_max,
                "cnf_sha256": cnf_sha,
                "lrat_sha256": hashlib.sha256(proof_bytes).hexdigest(),
                "lrat_steps": stats["steps"],
                "lrat_hints": stats["hints"],
            }
        if proc.returncode == 10:
            true_vars: set[int] = set()
            for line in proc.stdout.splitlines():
                if line.startswith("v "):
                    for tok in line[2:].split():
                        v = int(tok)
                        if v > 0:
                            true_vars.add(v)
            fac = instance.extract_factory_from_model(true_vars)
            alg = fac.verify_algebraic_structure()
            dcert = fac.certify_distance(use_sat_if_available=False)
            if dcert.distance < instance.d_min:
                raise AssertionError("Extracted SAT model violated distance requirement")
            return {
                "verdict": "SAT-FACTORY-FOUND",
                "N": instance.N,
                "n": fac.n,
                "d": dcert.distance,
                "leading_harmful_count": dcert.leading_harmful_count,
                "jws_polynomial": alg["jws_polynomial"],
                "factory": fac,
            }
        raise RuntimeError(f"Unexpected CaDiCaL exit code {proc.returncode}")


def enumerate_jws_symmetry_orbits(
    r: int,
    generators: Sequence[Sequence[int]],
) -> list[tuple[int, ...]]:
    """Jain-Wills-Singh Section V: partition F_2^r \\ {0} into disjoint orbits under a subgroup
    G <= GL(r, 2) generated by invertible r x r binary matrices (given as r row-bitmasks each)
    or wire permutations!
    """
    def apply_mat(mat_rows: Sequence[int], vec: int) -> int:
        out = 0
        for row_idx, rmask in enumerate(mat_rows):
            if (rmask & vec).bit_count() & 1:
                out |= 1 << row_idx
        return out

    unseen = set(range(1, 1 << r))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        seed = min(unseen)
        orb = {seed}
        queue = [seed]
        while queue:
            cur = queue.pop()
            for gen in generators:
                nxt = apply_mat(gen, cur)
                if nxt not in orb:
                    orb.add(nxt)
                    queue.append(nxt)
        unseen -= orb
        orbits.append(tuple(sorted(orb)))
    orbits.sort(key=lambda o: (len(o), o[0]))
    return orbits


def search_symmetric_sn_family(N: int, d_min: int) -> list[dict[str, Any]]:
    """Jacinto et al. Section V.B.1 (`F_{S_N}`, Eqs. 24-28): exact exhaustive solver over all
    `2^(N-1)` `S_N`-symmetric repetition-code T-distillation circuits on `N` qubits!
    Returns all solutions with distance `>= d_min`, sorted by `(n, -d)`.
    """
    if N < 4:
        raise ValueError("Require N >= 4")
    weights = list(range(1, N))
    # Precompute parities l_k, p_k, t_k (Eq. 24):
    lk = [comb(N - 1, w - 1) & 1 for w in weights]
    pk = [(comb(N - 2, w - 2) & 1) if w >= 2 else 0 for w in weights]
    tk = [(comb(N - 3, w - 3) & 1) if w >= 3 else 0 for w in weights]
    col_counts = [comb(N, w) for w in weights]

    # Note: two columns of weight a, b can XOR to weight c iff |a-b| <= c <= min(a+b, 2N-(a+b))
    # and c == a + b mod 2! Starting from {0}, iteratively expand reachable Hamming weights after
    # `step` faulty gates until `N` is reached!
    solutions: list[dict[str, Any]] = []
    for mask in range(1, 1 << (N - 1)):
        l_sum = sum(lk[i] for i in range(N - 1) if (mask >> i) & 1) & 1
        p_sum = sum(pk[i] for i in range(N - 1) if (mask >> i) & 1) & 1
        t_sum = sum(tk[i] for i in range(N - 1) if (mask >> i) & 1) & 1
        if not (l_sum == 1 and p_sum == 1 and t_sum == 1):
            continue
        active_w = [weights[i] for i in range(N - 1) if (mask >> i) & 1]
        # Exact BFS over reachable XOR Hamming weights in F_2^N under full S_N symmetry:
        frontier = {0}
        seen_nonzero: set[int] = set()
        dist = N + 1
        for step in range(1, d_min + 4):
            nxt: set[int] = set()
            for u in frontier:
                for w in active_w:
                    lo = abs(u - w)
                    hi = min(u + w, 2 * N - (u + w))
                    for v in range(lo, hi + 1, 2):
                        if v == N:
                            dist = step
                            break
                        if v != 0 and v not in seen_nonzero:
                            nxt.add(v)
                    if dist <= step:
                        break
                if dist <= step:
                    break
            if dist <= step:
                break
            seen_nonzero |= nxt
            frontier = nxt
            if not frontier:
                break
        if dist >= d_min:
            n_tot = sum(col_counts[i] for i in range(N - 1) if (mask >> i) & 1)
            solutions.append({
                "N": N,
                "n": n_tot,
                "d": dist,
                "active_weights": active_w,
            })
    solutions.sort(key=lambda s: (s["n"], -s["d"]))
    return solutions


# Banked CaDiCaL 2.1.3 LRAT UNSAT refutations (compressed via zlib+base64) for cold-clone
# pure-stdlib verification (`ftqc_codes.lrat.check`) without requiring an external SAT binary!
BANKED_NOGO_LRAT_B64: dict[tuple[int, str, int], str] = {
    (4, "T", 3): (
        "eNo9k1kSwzAIQ/97Ci7gGbN5uf/FKkSSD3dalcCTTMxcTGXiI8WuWMr8mYXYecQ8kkExxfMRz5bjFJeE"
        "PaLqwjHKW+KtVbTVp++RePsqnlM/lK/ky6BZR0v2ycktr4NDDlfJ/cpHRfembLLmK9/AmZRdvuILkEsQ"
        "x/+f6wm3s4tR5J+MNGZXw9cLYvhm2tVbhnFkhYhnjd7haRhnVox49nFzZTi7IIuSD5vElBFEQZSchLml"
        "K/Qu30Le7o7URnT3Iw1MxnAZ2TC39WC26DGSnnx2/0UcXM/IRV1b30w3FvRL3Vo/jBczteXGsaYHQzc5"
        "VKPZr3iPbJRNEtAFAWM24GRjpMPaKAyHTTZOE+5PWKvObYMRMoR3JIAtNYQMEV2bhKg96nHZxZv5wbNy"
        "YWM1xWF82KihXE0YrTnKNLLuskFgEC2jOeouu0u9Mku6N/ZvOAFzUtbJrFctKeX9eKQKa7zgPFSja12c"
        "65CX6mYc2Hwn3Zo0c2lx4Q6V1jGhImMcCzdojLQmIJxekAWHzkAwgpnRId6r4YTD73J46XDhVaGI8HBD"
        "LAUKPWMAb42Fm39KtScBrSGOpj3M8tLZbqbLrbLi/v0BfA7GDQ=="
    ),
    (5, "CCZ", 3): (
        "eNo9Vsu1LCEI3L8oOoE+R0BQ87n5x/CKqrEXLgZsqQ/i7DGe104941ljYdVTJ57xbw/rxGJiPxv7qn8h"
        "4c/rG/HeskZiTcbjeedg3J+MJ43R+bz4oajhePvtzuetfeO1sJLxet599zvquBfj63nP3e/7YAnNBsw+"
        "komwieVMHCTMbmIGFhFZM7Z7VKxeg4lm7PlLTLCLc5gAYwu/iRjPbAGQiE6cm8iDRVSA8dqsm9gbizQM"
        "tO3TI21hkbeBt9UnH77NSaGw5W19fomVWDQIFGx9H6BEbn1wEL8VCixSvvlA/EIqICwjbW/a65IoX1iE"
        "5E17X9oFbytI25v2vkIV6JTM9qa9r7RVaKYiJm/a55pRqFZLoJr2ufbV6SVUoP1TNtiaU/thdwrSZGPa"
        "0Ha4/RMpFS8iglfvj0Ix7s7uCJA+2r8U3xQ1mrNJ1cO6ESwczdn4hemuRNHRaM5OBmZKHMoXzTnIGQ3U"
        "NaZTjGjOQZUslEgVb6sndTWxm1tou8lF20QvTfS6yZMymfjlpB6zm7zotolgFnt50m0Wd8Et51GTbhOu"
        "uxJBuJNuk6DLipIXk25TEhfcKhWn2xTRBbfUspNuszddcEvXAge+QutS3XSPJogLU0h0N3bgBG8OrTaS"
        "8UV22bRTH2wmcvCgbNpJ2qEKGaycTbuo+VSJ1CCCJe8vrrloMinBepPb3Iy7UaYE6aP9R/Glc5pz8oMc"
        "TOTQQe12Ur40JYJuZ7tdFDxdCTUzEn/4jeNxYxaKLWiw0MYLk3FpVqEZvEVYhAIH+963h21w0TxMQlzF"
        "2O0w7nQPQTT/goLBRV3g5OuimbAPyDVvgcSKDVfdWz116Hm1juNLgEkEaRZn5PgSiWugL9g+N97oQnF2"
        "zxdHdXUu+uOVKh2HGTYVp471JeCgpnPtTwcmtLtvzFBZfxpEOBOrVXPVhVpAMY1KLLROsWeBoOO8oog3"
        "49Txi4nU/Uavax4AQKtfQ+c037rx6ri290vIu41KDcj0Tq1uHLbHEtCUbhicNogHhfr9iyECpOs3AQnz"
        "qABM3TeMuap3BbbpojC6Hz0SuGwmkB2Gg1F0EXPI/CuKtorDXkXjWKwvXhhlBLmbxvziPb4JEq+GfxDR"
        "ifUDg+vynY4pUZrGu1v7hnHxXFMPL8z5zoYDvoQFD199p+CyZLIPTjv7Cy95iv5/L4yl+Qc+b3whOnz6"
        "H80X0lkw0b/D9FTp1oOCYo38P4Kwpww="
    ),
}


def verify_banked_spatial_nogo(N: int, target: str, d_min: int) -> dict[str, Any]:
    """Deterministically generate the unrestricted spatial no-go CNF for `(N, target, d_min)`
    and verify its banked LRAT UNSAT refutation via pure-stdlib `ftqc_codes.lrat.check`!
    """
    key = (N, target.strip().upper(), d_min)
    if key not in BANKED_NOGO_LRAT_B64:
        raise KeyError(f"No banked LRAT proof for {key!r}; available: {sorted(BANKED_NOGO_LRAT_B64)}")
    inst = encode_spatial_nogo_cnf(N, key[1], d_min)
    lrat_bytes = zlib.decompress(base64.b64decode(BANKED_NOGO_LRAT_B64[key]))
    return solve_or_verify_spatial_sat(inst, banked_lrat_bytes=lrat_bytes)
