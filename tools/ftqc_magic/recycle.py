"""Mid-circuit check-qubit recycling & active-patch footprint minimizer.

Implements Section VI and Appendix E of Jacinto, Valcarce, Barizien, Gouzien &
Sangouard (2026, arXiv:2606.07734), extended to general multi-output binary and
Galois-binarized factories (`k >= 1`, `gamma in {1,3,5,7}^n`):

1. **Generator Weight Minimization (`minimize_generator_weights`)**:
   - Output-row coset search `r_a ^ sum_{i in S} c_i` minimizing Hamming weight of each
     output row so output qubits are allocated as late as possible (since output qubits
     remain active through column `n - 1`).
   - Greedy pairwise check-row weight reduction `c_i <- c_i ^ c_{i'}` shortening active
     check intervals `[f_i, l_i]`.
2. **Branch-and-Bound Column Permutation Optimizer (`compress_active_qubits`)**:
   - Finds a column permutation `pi` minimizing peak active logical qubits
     `A(G) = max_t |{i : f_i <= t <= l_tilde_i}|` (`l_tilde_i = n - 1` for output rows,
     `l_i` for check rows).
   - Emits an explicit physical patch schedule (`INIT_PLUS`, `PPR`, `MX_POSTSELECT`) on
     `A(G)` physical patches and **independently replays** the schedule to certify that
     the recycled circuit realizes the exact target factory `[[n, k, d]]`.
3. **SAT + LRAT Column-Ordering Optimality Verifier (`certify_active_qubits_sat`)**:
   - Linearizes Eqs. (E5)-(E6) into CNF (bijection `x_{c,t}`, prefix `o_{i,t}`, suffix
     `u_{i,t}`, active `a_{i,t} = o_{i,t} & u_{i,t}`, and Sinz cardinality `sum_i a_{i,t} <= B`),
     proving via CaDiCaL + `ftqc_codes.lrat` when `B = A(G) - 1` is UNSAT!
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
import tempfile
from typing import Any, Sequence

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


@dataclass(frozen=True)
class RecycledSchedule:
    """Certified mid-circuit check-recycled implementation of a `BinaryFactory`."""

    original_wires: int
    active_qubits: int
    column_order: tuple[int, ...]
    reduced_factory: BinaryFactory
    active_trace: tuple[int, ...]
    # Physical patch events: ("INIT_PLUS", phys_q, logical_wire), ("PPR", phys_mask, gamma),
    # ("MX_CHECK", phys_q, check_idx)
    events: tuple[tuple[Any, ...], ...]
    output_physical_qubits: tuple[int, ...]
    bnb_exhaustive: bool
    lrat_lower_bound_info: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "original_wires": self.original_wires,
            "active_qubits": self.active_qubits,
            "qubits_saved": self.original_wires - self.active_qubits,
            "column_order": list(self.column_order),
            "active_trace": list(self.active_trace),
            "output_physical_qubits": list(self.output_physical_qubits),
            "event_count": len(self.events),
            "bnb_exhaustive": self.bnb_exhaustive,
            "lrat_lower_bound_info": self.lrat_lower_bound_info,
        }


def minimize_generator_weights(factory: BinaryFactory, *, reduce_checks: bool = True) -> BinaryFactory:
    """Jacinto et al. Appendix E.2: minimize Hamming weights of output rows modulo check span
    and (when `reduce_checks=True`) greedily reduce check-row weights pairwise while preserving
    the exact CSS code & target!
    """
    chk = list(factory.check_rows)
    r = len(chk)

    # 1. Optional greedy pairwise check-row weight reduction:
    improved = reduce_checks
    while improved:
        improved = False
        for i in range(r):
            wi = chk[i].bit_count()
            best_j = -1
            best_w = wi
            for j in range(r):
                if i == j:
                    continue
                cand_w = (chk[i] ^ chk[j]).bit_count()
                if cand_w < best_w:
                    best_w = cand_w
                    best_j = j
            if best_j >= 0:
                chk[i] ^= chk[best_j]
                improved = True

    # Also check 3-way check combinations when r <= 15 for tighter intervals:
    if reduce_checks and r <= 15:
        improved3 = True
        while improved3:
            improved3 = False
            for i in range(r):
                wi = chk[i].bit_count()
                for j in range(r):
                    if j == i:
                        continue
                    for l in range(j + 1, r):
                        if l == i:
                            continue
                        cand = chk[i] ^ chk[j] ^ chk[l]
                        if cand.bit_count() < wi:
                            chk[i] = cand
                            wi = cand.bit_count()
                            improved3 = True

    # 2. Exact coset search for each output row when r <= 15 (else greedy 1/2-flip descent):
    new_out: list[int] = []
    if r <= 15:
        # Precompute all 2^r check combinations via Gray code:
        span_vals = [0] * (1 << r)
        for mask in range(1, 1 << r):
            lsb = mask & -mask
            bit_idx = lsb.bit_length() - 1
            span_vals[mask] = span_vals[mask ^ lsb] ^ chk[bit_idx]
        for out_r in factory.output_rows:
            best_row = out_r
            best_w = out_r.bit_count()
            for s_comb in span_vals:
                cand = out_r ^ s_comb
                cw = cand.bit_count()
                if cw < best_w:
                    best_w = cw
                    best_row = cand
            new_out.append(best_row)
    else:
        for out_r in factory.output_rows:
            cur = out_r
            step_imp = True
            while step_imp:
                step_imp = False
                for c_row in chk:
                    if (cur ^ c_row).bit_count() < cur.bit_count():
                        cur ^= c_row
                        step_imp = True
            new_out.append(cur)

    return BinaryFactory(
        k=factory.k,
        r=factory.r,
        n=factory.n,
        output_rows=tuple(new_out),
        check_rows=tuple(chk),
        gamma=factory.gamma,
        name=factory.name,
    )


def evaluate_active_trace(factory: BinaryFactory, column_order: Sequence[int]) -> tuple[int, tuple[int, ...]]:
    """Compute exact active-qubit trace `a(t)` (`0 <= t < n`) and peak `A(G)` under `column_order`."""
    n = factory.n
    N = factory.total_wires
    if len(column_order) != n or sorted(column_order) != list(range(n)):
        raise ValueError("column_order must be a permutation of 0..n-1")
    rows = factory.all_rows()
    first_pos = [n] * N
    last_pos = [-1] * N
    for t, col in enumerate(column_order):
        for i in range(N):
            if (rows[i] >> col) & 1:
                if t < first_pos[i]:
                    first_pos[i] = t
                if t > last_pos[i]:
                    last_pos[i] = t
    # Output rows 0..k-1 stay active until t = n - 1 once opened!
    for a in range(factory.k):
        if last_pos[a] >= 0:
            last_pos[a] = n - 1

    trace: list[int] = []
    for t in range(n):
        cnt = sum(1 for i in range(N) if first_pos[i] <= t <= last_pos[i])
        trace.append(cnt)
    return (max(trace) if trace else 0), tuple(trace)


def _greedy_column_order(factory: BinaryFactory, start_seed: int = 0) -> list[int]:
    """Fast greedy heuristic: delay opening output wires and prioritize closing already-open check wires."""
    n = factory.n
    k = factory.k
    N = factory.total_wires
    rows = factory.all_rows()
    col_masks = [factory.column_mask(c) for c in range(n)]
    rem_ones = [rows[i].bit_count() for i in range(N)]
    opened_mask = 0
    unplaced = set(range(n))
    order: list[int] = []

    for step in range(n):
        best_c = -1
        best_key: tuple[int, int, int, int] | None = None
        for c in unplaced:
            cm = col_masks[c]
            new_open = cm & ~opened_mask
            new_opened_mask = opened_mask | cm
            # Active during step `c` is everyone in `new_opened_mask` not yet closed BEFORE `c`:
            active_now = new_opened_mask.bit_count()
            # Who closes AFTER `c`? Any check i >= k with rem_ones[i] == 1 and bit i set in cm!
            closes_after = 0
            opens_output = (new_open & ((1 << k) - 1)).bit_count()
            rem_score = 0
            v = cm >> k
            idx = k
            while v:
                if v & 1:
                    if rem_ones[idx] == 1:
                        closes_after += 1
                    rem_score += rem_ones[idx]
                v >>= 1
                idx += 1
            # Lexicographic priority:
            # 1) minimize active_now, 2) avoid opening output wires early, 3) maximize closes_after,
            # 4) prefer columns touching nearly-finished check wires!
            tie = (c + start_seed) % n
            key = (active_now, opens_output, -closes_after, rem_score * n + tie)
            if best_key is None or key < best_key:
                best_key = key
                best_c = c

        unplaced.remove(best_c)
        order.append(best_c)
        cm = col_masks[best_c]
        opened_mask |= cm
        # Decrement remaining 1s and release closed check wires from `opened_mask`:
        for i in range(k, N):
            if (cm >> i) & 1:
                rem_ones[i] -= 1
                if rem_ones[i] == 0:
                    opened_mask &= ~(1 << i)
    return order


def _build_and_verify_physical_schedule(
    reduced: BinaryFactory,
    column_order: Sequence[int],
) -> tuple[int, tuple[int, ...], tuple[tuple[Any, ...], ...], tuple[int, ...]]:
    """Allocate physical qubits `0..A-1` greedily from a free pool as wires open/close and
    independently replay the physical events to prove exact equivalence to `reduced`!
    """
    peak, trace = evaluate_active_trace(reduced, column_order)
    n = reduced.n
    k = reduced.k
    N = reduced.total_wires
    rows = reduced.all_rows()
    first_pos = [n] * N
    last_pos = [-1] * N
    for t, col in enumerate(column_order):
        for i in range(N):
            if (rows[i] >> col) & 1:
                if t < first_pos[i]:
                    first_pos[i] = t
                if t > last_pos[i]:
                    last_pos[i] = t
    for a in range(k):
        if last_pos[a] >= 0:
            last_pos[a] = n - 1

    free_pool = list(range(peak - 1, -1, -1))
    wire_to_phys: dict[int, int] = {}
    events: list[tuple[Any, ...]] = []

    for t, col in enumerate(column_order):
        # 1. Initialize any wire whose first appearance is `t`:
        for i in range(N):
            if first_pos[i] == t:
                if not free_pool:
                    raise AssertionError("Physical pool exhausted below certified peak")
                pq = free_pool.pop()
                wire_to_phys[i] = pq
                events.append(("INIT_PLUS", pq, i))
        # 2. Apply PPR for `col`:
        phys_mask = 0
        for i in range(N):
            if (rows[i] >> col) & 1:
                phys_mask |= 1 << wire_to_phys[i]
        events.append(("PPR", phys_mask, reduced.gamma[col], col))
        # 3. Measure & recycle any check wire `i >= k` whose last appearance is `t`:
        for i in range(k, N):
            if last_pos[i] == t:
                pq = wire_to_phys.pop(i)
                events.append(("MX_CHECK", pq, i - k))
                free_pool.append(pq)
                free_pool.sort(reverse=True)

    # Handle any degenerate output wire that was never touched (cannot happen on non-idle factories):
    out_phys: list[int] = []
    for a in range(k):
        if a not in wire_to_phys:
            pq = free_pool.pop()
            wire_to_phys[a] = pq
            events.append(("INIT_PLUS", pq, a))
        out_phys.append(wire_to_phys[a])

    # Independent replay check: map each physical qubit's current logical wire identity through `events`
    # and reconstruct the N x n matrix to verify exact equality with `reduced`!
    phys_owner: dict[int, int] = {}
    reconstructed_cols: dict[int, tuple[int, int]] = {}
    closed_checks: set[int] = set()
    for ev in events:
        kind = ev[0]
        if kind == "INIT_PLUS":
            _, pq, w_id = ev
            if pq in phys_owner:
                raise AssertionError(f"Physical qubit {pq} re-initialized while still live!")
            phys_owner[pq] = w_id
        elif kind == "PPR":
            _, pmask, g_val, orig_col = ev
            logical_mask = 0
            for pq in range(peak):
                if (pmask >> pq) & 1:
                    if pq not in phys_owner:
                        raise AssertionError(f"PPR touches uninitialized physical qubit {pq}")
                    logical_mask |= 1 << phys_owner[pq]
            reconstructed_cols[orig_col] = (logical_mask, g_val)
        elif kind == "MX_CHECK":
            _, pq, chk_id = ev
            owner = phys_owner.pop(pq)
            if owner != k + chk_id:
                raise AssertionError("Check measurement wire mismatch")
            closed_checks.add(chk_id)

    if len(closed_checks) != sum(1 for i in range(k, N) if rows[i] != 0):
        raise AssertionError("Not all active check wires were measured and postselected")
    for c in range(n):
        got_mask, got_g = reconstructed_cols[c]
        if got_mask != reduced.column_mask(c) or got_g != reduced.gamma[c]:
            raise AssertionError(f"Recycled replay mismatch at column {c}")

    return peak, trace, tuple(events), tuple(out_phys)


def _search_fixed_basis(
    reduced: BinaryFactory,
    *,
    max_bnb_nodes: int,
) -> tuple[int, list[int], bool]:
    n = reduced.n
    k = reduced.k
    N = reduced.total_wires
    rows = reduced.all_rows()
    col_masks = [reduced.column_mask(c) for c in range(n)]
    total_ones = [rows[i].bit_count() for i in range(N)]

    best_order = list(range(n))
    best_peak, _ = evaluate_active_trace(reduced, best_order)
    for seed in range(min(n, 20)):
        cand_order = _greedy_column_order(reduced, start_seed=seed)
        cand_peak, _ = evaluate_active_trace(reduced, cand_order)
        if cand_peak < best_peak:
            best_peak = cand_peak
            best_order = cand_order

    trivial_lb = k + (1 if reduced.r > 0 else 0)
    nodes_visited = 0
    budget_hit = False

    rem_counts = list(total_ones)
    partial: list[int] = []
    used = [False] * n
    out_bits_mask = (1 << k) - 1

    def explore(opened_mask: int, running_peak: int) -> None:
        nonlocal best_peak, best_order, nodes_visited, budget_hit
        if running_peak >= best_peak or best_peak <= trivial_lb:
            return
        if len(partial) == n:
            best_peak = running_peak
            best_order = list(partial)
            return
        nodes_visited += 1
        if nodes_visited > max_bnb_nodes:
            budget_hit = True
            return

        candidates: list[tuple[int, int, int, int, int]] = []
        for c in range(n):
            if used[c]:
                continue
            cm = col_masks[c]
            next_open = opened_mask | cm
            act_now = next_open.bit_count()
            if act_now >= best_peak:
                continue
            opens_out = ((cm & ~opened_mask) & out_bits_mask).bit_count()
            closes = 0
            rem_sum = 0
            v = cm >> k
            idx = k
            while v:
                if v & 1:
                    if rem_counts[idx] == 1:
                        closes += 1
                    rem_sum += rem_counts[idx]
                v >>= 1
                idx += 1
            candidates.append((act_now - closes, opens_out, act_now, rem_sum, c))

        candidates.sort()
        for _, _, act_now, _, c in candidates:
            if running_peak >= best_peak or budget_hit:
                break
            new_peak = act_now if act_now > running_peak else running_peak
            if new_peak >= best_peak:
                continue
            cm = col_masks[c]
            next_open = opened_mask | cm
            closed_mask = 0
            v = cm >> k
            idx = k
            while v:
                if v & 1:
                    rem_counts[idx] -= 1
                    if rem_counts[idx] == 0:
                        closed_mask |= 1 << idx
                v >>= 1
                idx += 1
            used[c] = True
            partial.append(c)
            explore(next_open & ~closed_mask, new_peak)
            partial.pop()
            used[c] = False
            v = cm >> k
            idx = k
            while v:
                if v & 1:
                    rem_counts[idx] += 1
                v >>= 1
                idx += 1

    if best_peak > trivial_lb:
        explore(0, 0)
    return best_peak, best_order, budget_hit


def compress_active_qubits(
    factory: BinaryFactory,
    *,
    max_bnb_nodes: int = 80_000,
    certify_lrat: bool = False,
) -> RecycledSchedule:
    """Minimize active logical qubits `A(G)` via coset row-weight reduction + multi-start greedy +
    exact branch-and-bound (`max_bnb_nodes` budget), plus optional CaDiCaL SAT+LRAT optimality proof!
    """
    orig_alg = factory.verify_algebraic_structure()
    N = factory.total_wires

    # Try both output-only coset reduction (preserving block-structured check rows as in Fig. 4 of
    # Jacinto et al.!) and check+output reduction, keeping whichever yields fewer active qubits!
    red_out_only = minimize_generator_weights(factory, reduce_checks=False)
    red_both = minimize_generator_weights(factory, reduce_checks=True)
    candidates_red = [red_out_only] if red_both == red_out_only else [red_out_only, red_both]

    half_budget = max(1, max_bnb_nodes // len(candidates_red))
    best_reduced = red_out_only
    best_peak = N + 1
    best_order: list[int] = list(range(factory.n))
    any_budget_hit = False

    for cand_red in candidates_red:
        red_alg = cand_red.verify_algebraic_structure()
        if not orig_alg["target"].cnot_diagonal_clifford_equivalent(red_alg["target"]):
            raise AssertionError("Row weight reduction altered target equivalence class")
        c_peak, c_ord, b_hit = _search_fixed_basis(cand_red, max_bnb_nodes=half_budget)
        if c_peak < best_peak:
            best_peak = c_peak
            best_order = c_ord
            best_reduced = cand_red
            any_budget_hit = b_hit

    reduced = best_reduced
    n = reduced.n
    budget_hit = any_budget_hit

    lrat_info: dict[str, Any] | None = None
    if certify_lrat and best_peak > 1 and n <= 20:
        lrat_info = certify_active_qubits_sat(reduced, best_peak - 1)

    peak, trace, events, out_phys = _build_and_verify_physical_schedule(reduced, best_order)
    assert peak == best_peak
    return RecycledSchedule(
        original_wires=N,
        active_qubits=peak,
        column_order=tuple(best_order),
        reduced_factory=reduced,
        active_trace=trace,
        events=events,
        output_physical_qubits=out_phys,
        bnb_exhaustive=not budget_hit,
        lrat_lower_bound_info=lrat_info,
    )


def certify_active_qubits_sat(
    factory: BinaryFactory,
    max_active_target: int,
    *,
    timeout: int = 30,
) -> dict[str, Any] | None:
    """Encode Jacinto et al. Appendix E.4 (Eqs. E5-E6) as pure CNF and run CaDiCaL + `ftqc_codes.lrat`
    to prove UNSAT (no permutation of `factory` achieves `<= max_active_target` active qubits) or
    extract a better permutation if SAT!
    """
    _ensure_tools_on_syspath()
    try:
        from ftqc_codes import lrat as codes_lrat
        from ftqc_codes import sat as codes_sat
    except ImportError:
        return None

    solver = codes_sat.find_cadical()
    if solver is None:
        return None

    n = factory.n
    k = factory.k
    N = factory.total_wires
    rows = factory.all_rows()
    clauses: list[list[int]] = []
    nvars = 0

    def new_var() -> int:
        nonlocal nvars
        nvars += 1
        return nvars

    def add_at_most_one(lits: Sequence[int]) -> None:
        # Sequential Sinz AMO for linear clause count:
        m = len(lits)
        if m <= 1:
            return
        if m <= 4:
            for i in range(m):
                for j in range(i + 1, m):
                    clauses.append([-lits[i], -lits[j]])
            return
        s = [new_var() for _ in range(m - 1)]
        clauses.append([-lits[0], s[0]])
        for idx in range(1, m - 1):
            clauses.append([-lits[idx], s[idx]])
            clauses.append([-s[idx - 1], s[idx]])
            clauses.append([-lits[idx], -s[idx - 1]])
        clauses.append([-lits[m - 1], -s[m - 2]])

    def add_at_most_w(lits: Sequence[int], w: int) -> None:
        m = len(lits)
        if w >= m:
            return
        if w <= 0:
            for lit in lits:
                clauses.append([-lit])
            return
        s: dict[tuple[int, int], int] = {}
        for i in range(m - 1):
            for j in range(1, w + 1):
                s[(i, j)] = new_var()
        clauses.append([-lits[0], s[(0, 1)]])
        for j in range(2, w + 1):
            clauses.append([-s[(0, j)]])
        for i in range(1, m - 1):
            clauses.append([-lits[i], s[(i, 1)]])
            clauses.append([-s[(i - 1, 1)], s[(i, 1)]])
            for j in range(2, w + 1):
                clauses.append([-lits[i], -s[(i - 1, j - 1)], s[(i, j)]])
                clauses.append([-s[(i - 1, j)], s[(i, j)]])
            clauses.append([-lits[i], -s[(i - 1, w)]])
        clauses.append([-lits[m - 1], -s[(m - 2, w)]])

    # Assignment variables x[c][t] (0 <= c, t < n):
    x = [[new_var() for _ in range(n)] for _ in range(n)]
    for c in range(n):
        clauses.append(list(x[c]))
        add_at_most_one(x[c])
    for t in range(n):
        col_t = [x[c][t] for c in range(n)]
        clauses.append(col_t)
        add_at_most_one(col_t)

    # Opened o[i][t], still-used u[i][t], active a[i][t]:
    o = [[new_var() for _ in range(n)] for _ in range(N)]
    u = [[new_var() for _ in range(n)] for _ in range(N)]
    a = [[new_var() for _ in range(n)] for _ in range(N)]

    for i in range(N):
        supp_cols = [c for c in range(n) if (rows[i] >> c) & 1]
        for t in range(n):
            # If column c in supp(i) is placed at t, then o[i][t] = True and u[i][t] = True:
            for c in supp_cols:
                clauses.append([-x[c][t], o[i][t]])
                clauses.append([-x[c][t], u[i][t]])
            if t > 0:
                # Monotone prefix: o[i][t-1] => o[i][t]
                clauses.append([-o[i][t - 1], o[i][t]])
            if t + 1 < n:
                # Monotone suffix: u[i][t+1] => u[i][t]
                clauses.append([-u[i][t + 1], u[i][t]])
            if i < k:
                # Output qubit stays used through end: o[i][t] => a[i][t]
                clauses.append([-o[i][t], a[i][t]])
            else:
                # Check qubit active when both opened and still used: o[i][t] & u[i][t] => a[i][t]
                clauses.append([-o[i][t], -u[i][t], a[i][t]])

    for t in range(n):
        add_at_most_w([a[i][t] for i in range(N)], max_active_target)

    dimacs = f"p cnf {nvars} {len(clauses)}\n" + "\n".join(
        " ".join(str(l) for l in cl) + " 0" for cl in clauses
    ) + "\n"

    import hashlib
    import subprocess

    with tempfile.TemporaryDirectory(prefix="ftqc_magic_recycle_sat_") as tmp:
        cnf_p = str(Path(tmp) / "recycle.cnf")
        lrat_p = str(Path(tmp) / "recycle.lrat")
        Path(cnf_p).write_text(dimacs)
        proc = subprocess.run(
            [solver, "--lrat", "--no-binary", "-q", cnf_p, lrat_p],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if proc.returncode == 20:
            proof_bytes = Path(lrat_p).read_bytes()
            stats = codes_lrat.check(dimacs, proof_bytes)
            return {
                "verdict": "UNSAT-OPTIMALITY-PROVED",
                "refuted_active_target": max_active_target,
                "certified_min_active_qubits": max_active_target + 1,
                "cnf_sha256": hashlib.sha256(dimacs.encode()).hexdigest(),
                "lrat_sha256": hashlib.sha256(proof_bytes).hexdigest(),
                "lrat_steps": stats["steps"],
                "lrat_hints": stats["hints"],
            }
        return {"verdict": f"CADICAL-EXIT-{proc.returncode}", "refuted_active_target": max_active_target}
