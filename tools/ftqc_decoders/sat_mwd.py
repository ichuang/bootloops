"""Proof-producing inhomogeneous `SAT+LRAT` Exact Minimum-Weight Decoder (`MWD`).

Extends the cardinality-bounded CNF + `CaDiCaL` + `ftqc_codes.lrat` architecture
to inhomogeneous syndrome equations on arbitrary hypergraph / qLDPC codes:

    exists e in GF(2)^n :   H * e == s (mod 2)   and   wt(e) <= w.

Literature & Algorithmic Lineage:
  1. SAT / MaxSAT reductions for non-graphlike QEC decoding:
     L. Berent, L. Burgholzer, P.-J. H. S. Derks, J. Eisert, and R. Wille,
     "Decoding quantum color codes with MaxSAT", Quantum 8, 1506 (2024),
     arXiv:2303.14237 (`cda-tum/mqt-qecc`).
  2. Sequential-counter cardinality encoding:
     C. Sinz, "Towards an Optimal CNF Encoding of Boolean Cardinality Constraints",
     CP 2005, LNCS 3709, pp. 827-831.
  3. Certified LRAT UNSAT verification (solver is untrusted):
     L. Cruz-Filipe, M. Heule, W. Hunt, M. Kaufmann, and P. Schneider-Kamp,
     CADE 2017, verified in-tree by `ftqc_codes.lrat.check`; solver: `CaDiCaL`
     (A. Biere et al.).

Zero `float` or `complex` literals/calls appear anywhere in this module (wall-clock
diagnostics are reported in integer milliseconds).
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any

from .bp_osd import DecodeResult, decode_bp_osd
from .decimal_bounds import RefusalError
from .dem import DecodingProblem, bit_positions


def _ensure_tools_on_path() -> None:
    tools_dir = str(Path(__file__).resolve().parents[1])
    if tools_dir not in sys.path:
        sys.path.insert(0, tools_dir)


def find_cadical() -> str | None:
    """Locate `cadical` via `FTQC_CODES_CADICAL`, `QECCERT_CADICAL`, or `PATH`."""
    cand = os.environ.get("FTQC_CODES_CADICAL") or os.environ.get("QECCERT_CADICAL")
    if cand:
        return cand if os.path.isfile(cand) and os.access(cand, os.X_OK) else None
    return shutil.which("cadical")


class SyndromeFormula:
    """DIMACS CNF formula satisfiable iff `exists e in GF(2)^n: H e == s and wt(e) <= w`.

    Proof-carrying direction (`extend`): every `e` with `holds(e)` induces a full
    Boolean assignment `extend(e)` satisfying every clause of `self.clauses`. Thus an
    LRAT-verified UNSAT proof of `dimacs()` certifies that NO solution of weight `<= w`
    exists!
    """

    def __init__(self, problem: DecodingProblem, syndrome: int, max_weight: int) -> None:
        if max_weight < 0:
            raise RefusalError("INVALID-WEIGHT-BOUND", f"max_weight must be >= 0, got {max_weight}")
        if syndrome < 0 or (syndrome >> problem.num_checks):
            raise RefusalError("SYNDROME-OUT-OF-RANGE", f"Syndrome exceeds {problem.num_checks} checks")

        self.problem = problem
        self.n = problem.num_mechanisms
        self.syndrome = syndrome
        self.w = max_weight
        self.nvars = self.n
        self.clauses: list[list[int]] = []
        self._defs: list[tuple[int, int, int]] = []  # (t, a, b): x_t <-> (x_a xor x_b)
        self._cnt: dict[tuple[int, int], int] = {}
        self._build()

    def _new_var(self) -> int:
        self.nvars += 1
        return self.nvars

    def _xor_chain(self, row_mask: int) -> int | None:
        supp = [j + 1 for j in bit_positions(row_mask)]
        if not supp:
            return None
        acc = supp[0]
        for x_lit in supp[1:]:
            t = self._new_var()
            self._defs.append((t, acc, x_lit))
            self.clauses.extend(
                [[-t, acc, x_lit], [-t, -acc, -x_lit], [t, -acc, x_lit], [t, acc, -x_lit]]
            )
            acc = t
        return acc

    def _build(self) -> None:
        n = self.n
        w = self.w
        # 1. Inhomogeneous parity constraints `H_c . e == s_c (mod 2)`
        for c_idx, row_mask in enumerate(self.problem.check_rows):
            target = (self.syndrome >> c_idx) & 1
            lit = self._xor_chain(row_mask)
            if lit is None:
                if target == 1:
                    self.clauses.append([])  # Contradictory empty parity row = 1
            else:
                self.clauses.append([lit if target == 1 else -lit])

        # 2. Sinz (2005) sequential cardinality counter `wt(e) <= w`
        if w == 0:
            for i in range(1, n + 1):
                self.clauses.append([-i])
            return
        if w >= n or n <= 1:
            return

        s = self._cnt
        for i in range(1, n):
            for j in range(1, w + 1):
                s[i, j] = self._new_var()
        self.clauses.append([-1, s[1, 1]])
        for j in range(2, w + 1):
            self.clauses.append([-s[1, j]])
        for i in range(2, n):
            self.clauses.append([-i, s[i, 1]])
            self.clauses.append([-s[i - 1, 1], s[i, 1]])
            for j in range(2, w + 1):
                self.clauses.append([-i, -s[i - 1, j - 1], s[i, j]])
                self.clauses.append([-s[i - 1, j], s[i, j]])
            self.clauses.append([-i, -s[i - 1, w]])
        self.clauses.append([-n, -s[n - 1, w]])

    def holds(self, error_mask: int) -> bool:
        return (
            0 <= error_mask < (1 << self.n)
            and self.problem.syndrome_of(error_mask) == self.syndrome
            and error_mask.bit_count() <= self.w
        )

    def extend(self, error_mask: int) -> list[bool]:
        """Extend `error_mask` into a complete assignment `val[1..nvars]`."""
        val = [False] * (self.nvars + 1)
        for j in range(self.n):
            val[j + 1] = bool((error_mask >> j) & 1)
        for t, a, b in self._defs:
            val[t] = val[a] != val[b]
        run = 0
        for i in range(1, self.n):
            if val[i]:
                run += 1
            for j in range(1, self.w + 1):
                if (i, j) in self._cnt:
                    val[self._cnt[i, j]] = run >= j
        return val

    def satisfied_by(self, val: list[bool]) -> bool:
        return all(
            any(val[lit] if lit > 0 else (not val[-lit]) for lit in cl)
            for cl in self.clauses
        )

    def dimacs(self) -> str:
        lines = [f"p cnf {self.nvars} {len(self.clauses)}"]
        lines.extend(" ".join(str(x) for x in cl) + " 0" if cl else "0" for cl in self.clauses)
        return "\n".join(lines) + "\n"


def build_weight0_lrat_proof(formula: SyndromeFormula) -> bytes:
    """Construct an exact ASCII LRAT UNSAT proof in strict unit-propagation order
    when `formula.w == 0` and `syndrome != 0`!

    Traces unit propagation across `formula.clauses` to record only the active unit
    clause IDs leading to the conflict clause, producing a minimal LRAT proof verified
    by `ftqc_codes.lrat.check` on every cold clone without requiring `cadical`.
    """
    if formula.w != 0 or formula.syndrome == 0:
        raise RefusalError("WEIGHT0-LRAT-PRECONDITION", "Analytic RUP proof helper requires w == 0 and syndrome != 0")

    assign: dict[int, bool] = {}  # var -> True/False
    used_clause_ids: list[int] = []
    conflict = False

    while not conflict:
        progressed = False
        for cid_0, cl in enumerate(formula.clauses):
            cid = cid_0 + 1
            if cid in used_clause_ids:
                continue
            satisfied = False
            unassigned: list[int] = []
            for lit in cl:
                v = abs(lit)
                if v in assign:
                    if assign[v] == (lit > 0):
                        satisfied = True
                        break
                else:
                    unassigned.append(lit)
            if satisfied:
                continue
            if len(unassigned) == 0:
                used_clause_ids.append(cid)
                conflict = True
                progressed = True
                break
            if len(unassigned) == 1:
                u_lit = unassigned[0]
                assign[abs(u_lit)] = (u_lit > 0)
                used_clause_ids.append(cid)
                progressed = True
        if not progressed and not conflict:
            raise RefusalError("RUP-NOT-UNSAT", "Formula did not conflict by unit propagation")

    hints = " ".join(str(cid) for cid in used_clause_ids)
    line = f"{len(formula.clauses) + 1} 0 {hints} 0\n"
    return line.encode("ascii")


def verify_lrat_optimality(
    problem: DecodingProblem,
    syndrome: int,
    optimal_error_mask: int,
    lrat_proof_bytes: bytes,
) -> dict[str, Any]:
    """Independently certify via `ftqc_codes.lrat.check` that `optimal_error_mask` is a
    strictly minimum-Hamming-weight solution to `H * e == syndrome`.
    """
    _ensure_tools_on_path()
    from ftqc_codes import lrat

    if problem.syndrome_of(optimal_error_mask) != syndrome:
        raise RefusalError("INVALID-MWD-WITNESS", "Candidate recovery does not match target syndrome")
    w_star = optimal_error_mask.bit_count()
    if w_star == 0:
        return {"optimal_weight": 0, "lrat_verified": True, "lrat_steps": 0}

    refute_formula = SyndromeFormula(problem, syndrome, w_star - 1)
    dimacs_txt = refute_formula.dimacs()
    try:
        st = lrat.check(dimacs_txt, lrat_proof_bytes)
    except Exception as exc:
        raise RefusalError("LRAT-PROOF-REJECTED", f"LRAT optimality proof rejected: {exc}") from exc

    return {
        "optimal_weight": w_star,
        "refuted_weight_upper": w_star - 1,
        "cnf_sha256": hashlib.sha256(dimacs_txt.encode("utf-8")).hexdigest(),
        "lrat_sha256": hashlib.sha256(lrat_proof_bytes).hexdigest(),
        "lrat_steps": st["steps"],
        "lrat_hints": st["hints"],
        "lrat_verified": True,
    }


def _solve_cnf_with_cadical(
    formula: SyndromeFormula,
    solver_bin: str,
    timeout_sec: int = 60,
) -> tuple[str, int | None, bytes | None, dict[str, Any]]:
    """Invoke `CaDiCaL` in a scratch directory and check SAT model / LRAT UNSAT proof."""
    _ensure_tools_on_path()
    from ftqc_codes import lrat

    dimacs_txt = formula.dimacs()
    help_out = subprocess.run([solver_bin, "--help"], capture_output=True, text=True, timeout=10).stdout
    extra_flags = ["--modelalllits=true"] if "modelalllits" in help_out else []

    with tempfile.TemporaryDirectory(prefix="ftqc_decoders_sat_") as tmp:
        cnf_p = os.path.join(tmp, "syn.cnf")
        lrat_p = os.path.join(tmp, "syn.lrat")
        Path(cnf_p).write_text(dimacs_txt)
        t0 = time.monotonic_ns()
        proc = subprocess.run(
            [solver_bin, "--lrat", "--no-binary", "-q", *extra_flags, cnf_p, lrat_p],
            capture_output=True,
            text=True,
            timeout=timeout_sec,
        )
        elapsed_ms = (time.monotonic_ns() - t0) // 1_000_000

        meta = {
            "cnf_sha256": hashlib.sha256(dimacs_txt.encode("utf-8")).hexdigest(),
            "vars": formula.nvars,
            "clauses": len(formula.clauses),
            "solve_ms": elapsed_ms,
        }
        if proc.returncode == 20:
            proof_bytes = Path(lrat_p).read_bytes()
            st = lrat.check(dimacs_txt, proof_bytes)
            meta.update(
                lrat_sha256=hashlib.sha256(proof_bytes).hexdigest(),
                lrat_steps=st["steps"],
                lrat_hints=st["hints"],
            )
            return "unsat", None, proof_bytes, meta

        if proc.returncode == 10:
            word = 0
            seen: set[int] = set()
            for line in proc.stdout.splitlines():
                if line.startswith("v "):
                    for tok in line[2:].split():
                        v = int(tok)
                        seen.add(abs(v))
                        if 0 < v <= formula.n:
                            word |= 1 << (v - 1)
            if not all(j in seen for j in range(1, formula.n + 1)):
                raise RefusalError("CADICAL-PARTIAL-MODEL", "CaDiCaL printed an incomplete variable assignment")
            if not formula.holds(word):
                raise RefusalError("CADICAL-INVALID-MODEL", "CaDiCaL model failed direct syndrome/weight verification")
            return "sat", word, None, meta

        raise RefusalError("CADICAL-SOLVER-ERROR", f"Unexpected CaDiCaL exit code {proc.returncode}")


def decode_sat_mwd(
    problem: DecodingProblem,
    syndrome: int,
    *,
    timeout_sec: int = 60,
    solver: str | None = None,
) -> DecodeResult:
    """Exact Minimum-Hamming-Weight Decoder certified by `ftqc_codes.lrat`.

    * If `syndrome == 0`, returns `0` immediately (`optimal_weight = 0`).
    * If a weight-1 column matches `syndrome`, emits & checks the analytic `w=0` LRAT
      RUP refutation in pure Python (no external solver required!).
    * Otherwise drives `CaDiCaL` upward from `w = 1` (bounded above by `BP+OSD` witness)
      and verifies the `w* - 1` LRAT UNSAT proof via `ftqc_codes.lrat.check`.
    """
    if syndrome < 0 or (syndrome >> problem.num_checks):
        raise RefusalError("SYNDROME-OUT-OF-RANGE", f"Syndrome exceeds {problem.num_checks} checks")

    if syndrome == 0:
        return DecodeResult(
            decoder="sat_mwd",
            error_mask=0,
            converged=True,
            observable_mask=0,
            hamming_weight=0,
            soft_weight=0,
            stats={"lrat_certified_optimal": True, "optimal_weight": 0, "solver": "analytic"},
        )

    # Pure-Python certified fast path when a weight-1 mechanism matches `syndrome`!
    best_w1: tuple[int, int] | None = None
    for j in range(problem.num_mechanisms):
        if problem.mech_detectors[j] == syndrome:
            cand = (problem.weights[j], j)
            if best_w1 is None or cand < best_w1:
                best_w1 = cand
    if best_w1 is not None:
        e1 = 1 << best_w1[1]
        f0 = SyndromeFormula(problem, syndrome, 0)
        proof0 = build_weight0_lrat_proof(f0)
        cert_info = verify_lrat_optimality(problem, syndrome, e1, proof0)
        return DecodeResult(
            decoder="sat_mwd",
            error_mask=e1,
            converged=True,
            observable_mask=problem.observable_of(e1),
            hamming_weight=1,
            soft_weight=problem.soft_weight_of(e1),
            stats={"lrat_certified_optimal": True, "solver": "analytic_rup_lrat", **cert_info},
        )

    solver_bin = solver or find_cadical()
    if solver_bin is None:
        raise RefusalError(
            "CADICAL-NOT-FOUND",
            "decode_sat_mwd for weight >= 2 requires CaDiCaL on PATH or FTQC_CODES_CADICAL",
        )

    # Seed an upper bound via BP+OSD so we know feasibility upfront
    seed_res = decode_bp_osd(problem, syndrome, post_method="osd", osd_order=2)
    upper_w = seed_res.hamming_weight if seed_res.converged else problem.num_mechanisms
    best_word = seed_res.error_mask if seed_res.converged else None
    last_unsat_meta: dict[str, Any] = {}

    for w in range(1, upper_w + 1):
        if w == upper_w and best_word is not None and best_word.bit_count() == w:
            return DecodeResult(
                decoder="sat_mwd",
                error_mask=best_word,
                converged=True,
                observable_mask=problem.observable_of(best_word),
                hamming_weight=w,
                soft_weight=problem.soft_weight_of(best_word),
                stats={
                    "lrat_certified_optimal": True,
                    "optimal_weight": w,
                    "solver": "cadical+lrat",
                    "last_unsat": last_unsat_meta,
                },
            )
        f_w = SyndromeFormula(problem, syndrome, w)
        status, word, _proof, meta = _solve_cnf_with_cadical(f_w, solver_bin, timeout_sec=timeout_sec)
        if status == "sat":
            assert word is not None
            return DecodeResult(
                decoder="sat_mwd",
                error_mask=word,
                converged=True,
                observable_mask=problem.observable_of(word),
                hamming_weight=word.bit_count(),
                soft_weight=problem.soft_weight_of(word),
                stats={
                    "lrat_certified_optimal": True,
                    "optimal_weight": word.bit_count(),
                    "solver": "cadical+lrat",
                    "sat_meta": meta,
                    "last_unsat": last_unsat_meta,
                },
            )
        last_unsat_meta = meta

    return DecodeResult(
        decoder="sat_mwd",
        error_mask=0,
        converged=False,
        observable_mask=0,
        hamming_weight=0,
        soft_weight=0,
        stats={"lrat_certified_optimal": False, "infeasible_syndrome": True, "last_unsat": last_unsat_meta},
    )
