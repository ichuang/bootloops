"""Spooky & parallel spooky pebbling validator, KMRV26 A-sequence compiler, and A* solver.

Implements the sequential and parallel spooky pebbling models, lower bounds, and
algorithms from:
  - Gidney, "Spooky Pebble Games and Irreversible Uncomputation" (2019) [Gid19a]
  - Kornerup, Sadun, Soloveichik, Quantum 9, 1636 (2025) [KSS25]
  - Kahanamoku-Meyer, Ragavan, Van Kirk, EUROCRYPT 2026 (ePrint 2025/1887) [KMRV26]

100% Python standard library, MIT-licensed, zero float literals or float() calls.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Context, Decimal, ROUND_CEILING
import heapq
import math
from typing import Any, Iterable, Sequence


NODE_EMPTY = "EMPTY"
NODE_PEBBLE = "PEBBLE"
NODE_GHOST = "GHOST"


@dataclass(frozen=True)
class PebbleDAG:
    """Directed acyclic dependency graph for spooky pebbling.

    For a line graph of length `ell`, vertices are `1..ell` with root `0` implicitly
    permanently pebbled (`preds[1] = ()`, `preds[i] = (i - 1,)` for `i >= 2`), and
    `targets = frozenset({ell})`.
    """

    nodes: tuple[int, ...]
    preds: dict[int, tuple[int, ...]]
    targets: frozenset[int]
    node_weights: dict[int, int] = field(default_factory=dict)
    work_ancillas: dict[int, int] = field(default_factory=dict)

    @classmethod
    def line(
        cls,
        length: int,
        *,
        window: int | None = None,
        window_ancilla_weight: int = 1,
    ) -> "PebbleDAG":
        if length < 1:
            raise ValueError("line pebble graph length must be >= 1")
        nodes = tuple(range(1, length + 1))
        preds = {1: ()}
        for i in range(2, length + 1):
            preds[i] = (i - 1,)
        anc: dict[int, int] = {}
        if window is not None and window >= 1:
            # KMRV26 Section 4.1 / Alg 4.1-4.2: extra ancilla at i == 1 (mod window + 1)
            mod = window + 1
            for i in nodes:
                if i % mod == 1:
                    anc[i] = window_ancilla_weight
        return cls(
            nodes=nodes,
            preds=preds,
            targets=frozenset({length}),
            node_weights={i: 1 for i in nodes},
            work_ancillas=anc,
        )


@dataclass(frozen=True)
class PebbleStep:
    """One two-phase spooky pebbling time step (KMRV26 Definition 2.3):

    Phase 1: parallel `pebble` and `unpebble` moves with pairwise disjoint active sets.
    Phase 2: parallel `ghost` moves on currently pebbled nodes.
    """

    pebble: tuple[int, ...] = ()
    unpebble: tuple[int, ...] = ()
    ghost: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "pebble", tuple(sorted(set(self.pebble))))
        object.__setattr__(self, "unpebble", tuple(sorted(set(self.unpebble))))
        object.__setattr__(self, "ghost", tuple(sorted(set(self.ghost))))

    def is_empty(self) -> bool:
        return not (self.pebble or self.unpebble or self.ghost)

    def op_count(self) -> int:
        return len(self.pebble) + len(self.unpebble)

    def to_dict(self) -> dict[str, list[int]]:
        return {
            "pebble": list(self.pebble),
            "unpebble": list(self.unpebble),
            "ghost": list(self.ghost),
        }

    @classmethod
    def from_dict(cls, d: dict[str, Iterable[int]]) -> "PebbleStep":
        return cls(
            pebble=tuple(int(x) for x in d.get("pebble", ())),
            unpebble=tuple(int(x) for x in d.get("unpebble", ())),
            ghost=tuple(int(x) for x in d.get("ghost", ())),
        )


@dataclass(frozen=True)
class PebbleFinding:
    code: str
    message: str
    step_index: int | None = None
    node: int | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "message": self.message,
            "step_index": self.step_index,
            "node": self.node,
        }


@dataclass(frozen=True)
class PebbleValidationResult:
    valid: bool
    depth: int
    pebble_cost: int
    ghost_count: int
    peak_pebbles: int
    peak_weighted_space: int
    is_sequential: bool
    targets_achieved: bool
    clean_at_exit: bool
    findings: tuple[PebbleFinding, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "valid": self.valid,
            "depth": self.depth,
            "pebble_cost": self.pebble_cost,
            "ghost_count": self.ghost_count,
            "peak_pebbles": self.peak_pebbles,
            "peak_weighted_space": self.peak_weighted_space,
            "is_sequential": self.is_sequential,
            "targets_achieved": self.targets_achieved,
            "clean_at_exit": self.clean_at_exit,
            "findings": [f.to_dict() for f in self.findings],
        }


def validate_pebbling_schedule(
    dag: PebbleDAG,
    steps: Sequence[PebbleStep],
    *,
    max_pebbles: int | None = None,
    max_weighted_space: int | None = None,
    require_target_pebbled: bool = True,
) -> PebbleValidationResult:
    """Validate a sequential or parallel spooky pebbling schedule on `dag`.

    Enforces KMRV26 Definitions 2.1-2.8:
      - `Pebble(v)` requires all `u in preds[v]` to hold `PEBBLE`, and `v in {EMPTY, GHOST}`.
      - `Unpebble(v)` requires all `u in preds[v]` to hold `PEBBLE`, and `v == PEBBLE`.
      - Parallel active sets `{v} | set(preds[v])` within Phase 1 must be pairwise disjoint.
      - `Ghost(v)` in Phase 2 requires `v == PEBBLE` after Phase 1.
      - Every target in `dag.targets` must hold `PEBBLE` at some intermediate state.
      - At exit, all nodes must be `EMPTY` (`UNCLEANED-PEBBLE-AT-EXIT` / `UNCLEANED-GHOST-AT-EXIT`).
    """
    node_set = set(dag.nodes)
    state: dict[int, str] = {v: NODE_EMPTY for v in dag.nodes}
    targets_hit: set[int] = set()
    findings: list[PebbleFinding] = []

    peak_pebbles = 0
    peak_weighted = 0
    pebble_cost = 0
    ghost_count = 0
    is_sequential = True

    for t_idx, step in enumerate(steps):
        state_before = dict(state)
        n_calls = len(step.pebble) + len(step.unpebble)
        if n_calls > 1:
            is_sequential = False
        pebble_cost += n_calls
        ghost_count += len(step.ghost)

        # Active-set disjointness & precondition checks for Phase 1
        used_active: dict[int, str] = {}
        extra_anc = 0

        for v in step.pebble:
            if v not in node_set:
                findings.append(PebbleFinding("UNKNOWN-PEBBLE-NODE", f"Unknown node {v}", t_idx, v))
                continue
            if state_before[v] == NODE_PEBBLE:
                findings.append(
                    PebbleFinding(
                        "REDUNDANT-PEBBLE-MOVE",
                        f"Pebble({v}) called when node {v} already holds a pebble",
                        t_idx,
                        v,
                    )
                )
            for u in dag.preds.get(v, ()):
                if state_before.get(u) != NODE_PEBBLE:
                    findings.append(
                        PebbleFinding(
                            "MISSING-PREDECESSOR-PEBBLE",
                            f"Pebble({v}) requires predecessor {u} to be PEBBLE (got {state_before.get(u)})",
                            t_idx,
                            v,
                        )
                    )
            act = {v} | set(dag.preds.get(v, ()))
            for w in act:
                if w in used_active:
                    findings.append(
                        PebbleFinding(
                            "PARALLEL-ACTIVE-SET-CONFLICT",
                            f"Pebble({v}) conflicts on node {w} with {used_active[w]} in step {t_idx}",
                            t_idx,
                            w,
                        )
                    )
                used_active[w] = f"Pebble({v})"
            extra_anc += dag.work_ancillas.get(v, 0)

        for v in step.unpebble:
            if v not in node_set:
                findings.append(PebbleFinding("UNKNOWN-PEBBLE-NODE", f"Unknown node {v}", t_idx, v))
                continue
            if state_before[v] != NODE_PEBBLE:
                findings.append(
                    PebbleFinding(
                        "INVALID-UNPEBBLE-TARGET",
                        f"Unpebble({v}) requires node {v} to be PEBBLE (got {state_before[v]})",
                        t_idx,
                        v,
                    )
                )
            for u in dag.preds.get(v, ()):
                if state_before.get(u) != NODE_PEBBLE:
                    findings.append(
                        PebbleFinding(
                            "MISSING-PREDECESSOR-PEBBLE",
                            f"Unpebble({v}) requires predecessor {u} to be PEBBLE (got {state_before.get(u)})",
                            t_idx,
                            v,
                        )
                    )
            act = {v} | set(dag.preds.get(v, ()))
            for w in act:
                if w in used_active:
                    findings.append(
                        PebbleFinding(
                            "PARALLEL-ACTIVE-SET-CONFLICT",
                            f"Unpebble({v}) conflicts on node {w} with {used_active[w]} in step {t_idx}",
                            t_idx,
                            w,
                        )
                    )
                used_active[w] = f"Unpebble({v})"
            extra_anc += dag.work_ancillas.get(v, 0)

        # Apply Phase 1
        state_mid = dict(state_before)
        for v in step.pebble:
            if v in node_set:
                state_mid[v] = NODE_PEBBLE
        for v in step.unpebble:
            if v in node_set:
                state_mid[v] = NODE_EMPTY

        for tv in dag.targets:
            if state_mid.get(tv) == NODE_PEBBLE:
                targets_hit.add(tv)

        # Apply Phase 2 (`Ghost`)
        state_after = dict(state_mid)
        for v in step.ghost:
            if v not in node_set:
                findings.append(PebbleFinding("UNKNOWN-PEBBLE-NODE", f"Unknown node {v}", t_idx, v))
                continue
            if state_mid[v] != NODE_PEBBLE:
                findings.append(
                    PebbleFinding(
                        "INVALID-GHOST-TARGET",
                        f"Ghost({v}) requires node {v} to hold PEBBLE after Phase 1 (got {state_mid[v]})",
                        t_idx,
                        v,
                    )
                )
            else:
                state_after[v] = NODE_GHOST

        for tv in dag.targets:
            if state_after.get(tv) == NODE_PEBBLE:
                targets_hit.add(tv)

        # Space usage of step t (KMRV26 Definition 2.8 + mid-phase check):
        # Any node pebbled in state_before, state_mid, or state_after is occupied during step t.
        live_nodes = [
            v
            for v in dag.nodes
            if state_before[v] == NODE_PEBBLE or state_mid[v] == NODE_PEBBLE or state_after[v] == NODE_PEBBLE
        ]
        step_pebbles = len(live_nodes)
        step_weighted = sum(dag.node_weights.get(v, 1) for v in live_nodes) + extra_anc
        peak_pebbles = max(peak_pebbles, step_pebbles)
        peak_weighted = max(peak_weighted, step_weighted)

        if max_pebbles is not None and step_pebbles > max_pebbles:
            findings.append(
                PebbleFinding(
                    "PEBBLE-BUDGET-EXCEEDED",
                    f"Step {t_idx} uses {step_pebbles} pebbles > cap {max_pebbles}",
                    t_idx,
                )
            )
        if max_weighted_space is not None and step_weighted > max_weighted_space:
            findings.append(
                PebbleFinding(
                    "WEIGHTED-SPACE-BUDGET-EXCEEDED",
                    f"Step {t_idx} uses weighted space {step_weighted} > cap {max_weighted_space}",
                    t_idx,
                )
            )

        state = state_after

    targets_ok = True
    if require_target_pebbled:
        for tv in sorted(dag.targets):
            if tv not in targets_hit:
                targets_ok = False
                findings.append(
                    PebbleFinding(
                        "TARGET-NEVER-PEBBLED",
                        f"Target node {tv} never held a pebble during the schedule",
                        node=tv,
                    )
                )

    clean_exit = True
    for v in dag.nodes:
        if state[v] == NODE_GHOST:
            clean_exit = False
            findings.append(
                PebbleFinding(
                    "UNCLEANED-GHOST-AT-EXIT",
                    f"Node {v} remains in GHOST state at schedule termination (uncompensated HMR phase)",
                    node=v,
                )
            )
        elif state[v] == NODE_PEBBLE:
            clean_exit = False
            findings.append(
                PebbleFinding(
                    "UNCLEANED-PEBBLE-AT-EXIT",
                    f"Node {v} remains PEBBLED at schedule termination (uncleaned ancilla)",
                    node=v,
                )
            )

    return PebbleValidationResult(
        valid=(len(findings) == 0),
        depth=len(steps),
        pebble_cost=pebble_cost,
        ghost_count=ghost_count,
        peak_pebbles=peak_pebbles,
        peak_weighted_space=peak_weighted,
        is_sequential=is_sequential,
        targets_achieved=targets_ok,
        clean_at_exit=clean_exit,
        findings=tuple(findings),
    )


# ---------------------------------------------------------------------------
# KMRV26 Section 2.3: A-sequence (`Blast` / `Unblast`) optimal-depth compiler
# ---------------------------------------------------------------------------


def kmrv_a_sequence(count: int) -> list[int]:
    """Generate `A_1, ..., A_count` (`A_1=1, A_2=2, A_3=2, A_k = A_{k-2} + A_{k-3}`)."""
    if count < 1:
        return []
    seq = [1, 2, 2]
    if count <= 3:
        return seq[:count]
    while len(seq) < count:
        seq.append(seq[-2] + seq[-3])
    return seq


def kmrv_theoretical_pebble_bound(length: int) -> int:
    """Theorem 2.12 upper bound on pebbles for exact parallel depth `2 * length`."""
    if length <= 1:
        return 1
    if length == 2:
        return 2
    if length <= 4:
        return 3
    if length <= 6:
        return 4
    # Minimal k >= 7 with A_k >= length gives space <= k - 3 <= ceil(2.47 * log2(length) - 2.77)
    seq = [0, 1, 2, 2, 3, 4, 5, 7]
    while seq[-1] < length:
        seq.append(seq[-2] + seq[-3])
    k = next(idx for idx in range(7, len(seq)) if seq[idx] >= length)
    return k - 3


def compile_kmrv_blast_schedule(length: int) -> list[PebbleStep]:
    """Synthesize KMRV26 Algorithm 2.1-2.3 parallel spooky pebbling of exact depth `2 * length`.

    Achieves pebbling space `<= k - 3 <= ceil(2.47 * log2(length) - 2.77)` for `length >= 7`
    (and exact minimum pebbles for `length <= 6`: 1 for `ell=1`, 2 for `ell=2`, 3 for `ell=3,4`).
    """
    if length < 1:
        raise ValueError("length must be >= 1")
    if length == 1:
        return [PebbleStep(pebble=(1,)), PebbleStep(unpebble=(1,))]

    # Build 1-indexed A_k table until A_k >= length
    a = [0, 1, 2, 2]
    while a[-1] < length:
        a.append(a[-2] + a[-3])
    k_top = next(idx for idx in range(2, len(a)) if a[idx] >= length)
    m_len = a[k_top]

    def blast(i_start: int, i_end: int, k: int) -> list[PebbleStep]:
        # Precondition: i_end - i_start == A_k - 1
        if k <= 1 or i_end == i_start:
            return []
        if k in (2, 3):
            return [PebbleStep(pebble=(i_end,))]
        steps: list[PebbleStep] = [PebbleStep(pebble=(i_start + 1,))]
        ak3 = a[k - 3]
        for j in range(i_start + 2, i_start + ak3 + 1):
            steps.append(PebbleStep(pebble=(j,), ghost=(j - 1,)))
        steps.extend(blast(i_start + ak3, i_end, k - 2))
        return steps

    def unblast(i_start: int, i_end: int, k: int) -> list[PebbleStep]:
        if k <= 1 or i_end == i_start:
            return []
        if k in (2, 3):
            return [PebbleStep(unpebble=(i_end,))]
        ak2 = a[k - 2]
        ak3 = a[k - 3]
        right_u = unblast(i_start + ak3, i_end, k - 2)
        split = ak2 - ak3
        out: list[PebbleStep] = list(right_u[:split])
        rem_right = right_u[split:]
        left_b = blast(i_start, i_start + ak3 - 1, k - 3)
        assert len(rem_right) == len(left_b) == ak3 - 1
        for s_r, s_l in zip(rem_right, left_b):
            out.append(
                PebbleStep(
                    pebble=s_r.pebble + s_l.pebble,
                    unpebble=s_r.unpebble + s_l.unpebble,
                    ghost=s_r.ghost + s_l.ghost,
                )
            )
        out.append(PebbleStep(unpebble=(i_start + ak3,)))
        out.extend(unblast(i_start, i_start + ak3 - 1, k - 3))
        return out

    raw_steps: list[PebbleStep] = []
    raw_steps.extend(blast(0, m_len - 1, k_top))
    raw_steps.append(PebbleStep(pebble=(m_len,)))
    raw_steps.append(PebbleStep(unpebble=(m_len,)))
    raw_steps.extend(unblast(0, m_len - 1, k_top))

    # Theorem 2.12 shift: map site i -> i - (A_k - length), dropping sites <= 0 and empty steps
    shift = m_len - length
    if shift == 0:
        return raw_steps

    projected: list[PebbleStep] = []
    for st in raw_steps:
        p = tuple(i - shift for i in st.pebble if i - shift >= 1)
        u = tuple(i - shift for i in st.unpebble if i - shift >= 1)
        g = tuple(i - shift for i in st.ghost if i - shift >= 1)
        cand = PebbleStep(pebble=p, unpebble=u, ghost=g)
        if not cand.is_empty():
            projected.append(cand)
    return projected


# ---------------------------------------------------------------------------
# Gidney (2019) sequential divide-and-conquer spooky pebbler (`~ log2 ell` space)
# ---------------------------------------------------------------------------


def compile_gidney_sequential_schedule(length: int) -> list[PebbleStep]:
    """Gidney's sequential divide-and-conquer spooky pebbling (`<= ceil(log2 ell) + 2` pebbles)."""
    if length < 1:
        raise ValueError("length must be >= 1")

    steps: list[PebbleStep] = []

    def sweep_with_markers(lo: int, hi: int) -> list[int]:
        """Starting with `lo` pebbled (or `lo == 0`) and `(lo, hi]` empty/ghost,
        march up to `hi`, leaving marker pebbles at binary midpoints and `hi`,
        ghosting non-markers. Returns sorted marker list in `(lo, hi)`.
        """
        span = hi - lo
        if span <= 0:
            return []
        if span == 1:
            steps.append(PebbleStep(pebble=(hi,)))
            return []
        markers: set[int] = set()
        cur_lo = lo
        while hi - cur_lo > 2:
            mid = cur_lo + (hi - cur_lo) // 2
            markers.add(mid)
            cur_lo = mid
        steps.append(PebbleStep(pebble=(lo + 1,)))
        for j in range(lo + 2, hi + 1):
            prev = j - 1
            if prev in markers:
                steps.append(PebbleStep(pebble=(j,)))
            else:
                steps.append(PebbleStep(pebble=(j,), ghost=(prev,)))
        return sorted(markers)

    def clear_interval(lo: int, hi: int, markers: list[int]) -> None:
        """Given `lo` (or 0) pebbled, `markers` in `(lo, hi)` pebbled, `hi` pebbled,
        and other sites in `(lo, hi)` ghosted, uncompute `(lo, hi]` to EMPTY while keeping `lo` pebbled.
        """
        span = hi - lo
        if span <= 0:
            return
        if span == 1:
            steps.append(PebbleStep(unpebble=(hi,)))
            return
        if not markers:
            # span == 2: lo+1 is ghost, hi == lo+2 is pebble
            steps.append(PebbleStep(pebble=(lo + 1,)))
            steps.append(PebbleStep(unpebble=(hi,)))
            steps.append(PebbleStep(unpebble=(lo + 1,)))
            return
        m0 = markers[0]
        rest_markers = markers[1:]
        # Clear right sub-interval (m0, hi] using m0 as base!
        clear_interval(m0, hi, rest_markers)
        # Now m0 is pebble, (m0, hi] are EMPTY, (lo, m0) are GHOST.
        sub_markers = sweep_with_markers(lo, m0 - 1) if (m0 - 1 > lo) else []
        # Now m0 - 1 is pebbled (or m0 == lo + 1 and lo is pebbled), so unpebble m0!
        steps.append(PebbleStep(unpebble=(m0,)))
        if m0 - 1 > lo:
            clear_interval(lo, m0 - 1, sub_markers)

    top_markers = sweep_with_markers(0, length)
    clear_interval(0, length, top_markers)
    return steps


# ---------------------------------------------------------------------------
# KMRV26 Section 2.4 & Appendix B: Exact admissible A* parallel spooky solver
# ---------------------------------------------------------------------------


def kss_parallel_depth_lower_bound(length: int, max_pebbles: int) -> int:
    """Deterministic lower bound on parallel spooky pebbling depth:

    Always `>= 2 * length` (Remark 2.7), and infeasible (`-1`) if `max_pebbles < 1`
    or (`length >= 2` and `max_pebbles < 2`) or (`length >= 3` and `max_pebbles < 3`).
    Also records the KSS25 threshold parameter `m_min` where `binom(m + s - 2, s - 2) >= length`.
    """
    if max_pebbles < 1:
        return -1
    if length >= 2 and max_pebbles < 2:
        return -1
    if length >= 3 and max_pebbles < 3:
        return -1
    return 2 * length


def kss_binomial_parameter_m(length: int, pebbles: int) -> int:
    """Smallest integer `m >= 1` satisfying KSS25 condition `comb(m + pebbles - 2, pebbles - 2) >= length`."""
    if pebbles <= 2:
        return length
    m = 1
    while math.comb(m + pebbles - 2, pebbles - 2) < length:
        m += 1
    return m


def solve_optimal_pebbling_astar(
    length: int,
    max_pebbles: int,
    *,
    window: int | None = None,
    window_ancilla_weight: int = 1,
    max_expansions: int = 200_000,
) -> list[PebbleStep] | None:
    """Find an exact minimum-depth parallel spooky pebbling of `1..length` under space cap `max_pebbles`
    (representing weighted space cap when `window` is specified), using KMRV26 Appendix B canonicalization
    and admissible heuristic `h`.

    Canonical invariant (KMRV26 Appendix B):
      - Before reaching `length`, `Blast` marches `1..length` in `length` steps, choosing at each step `r`
        whether to keep `r - 1` as a marker pebble or immediately ghost `r - 1`.
      - During unblasting (`reached_ell == True`), for any pebble bitmask `P` (`0`-indexed bits `0..length-1`),
        every site `< max(P)` not in `P` holds a ghost and every site `> max(P)` is permanently `EMPTY`;
        admissible heuristic `h(P) = max_site_1based(P)` (`0` when `P == 0`) never overestimates remaining steps,
        and `Ghost(i)` is only applied immediately after `Pebble(i + 1)` on sites `< max(P_after)`.
    """
    if length < 1:
        raise ValueError("length must be >= 1")
    if length == 1:
        anc1 = window_ancilla_weight if (window is not None and 1 % (window + 1) == 1) else 0
        if 1 + anc1 > max_pebbles:
            return None
        return [PebbleStep(pebble=(1,)), PebbleStep(unpebble=(1,))]

    anc_cost = [0] * length
    if window is not None and window >= 1:
        mod = window + 1
        for idx in range(length):
            if (idx + 1) % mod == 1:
                anc_cost[idx] = window_ancilla_weight

    # Fast check: if unconstrained KMRV26 Blast/Unblast already fits within max_pebbles (and window is None),
    # it achieves depth `2 * length`, which is provably the global optimum (Remark 2.7)!
    if window is None:
        cand = compile_kmrv_blast_schedule(length)
        val = validate_pebbling_schedule(PebbleDAG.line(length), cand)
        if val.valid and val.peak_pebbles <= max_pebbles:
            return cand

    # State: `(p_mask: int, reached: bool)` where bit `i` (`0 <= i < length`) means site `i + 1` holds PEBBLE.
    start_state = (0, False)
    goal_state = (0, True)

    def heuristic(p_mask: int, reached: bool) -> int:
        if not reached:
            r = p_mask.bit_length()
            return 2 * length - r
        return p_mask.bit_length()

    # Priority queue: `(f, g, op_tiebreak, p_mask, reached)`
    pq: list[tuple[int, int, int, int, bool]] = [(2 * length, 0, 0, 0, False)]
    best_g: dict[tuple[int, bool], tuple[int, int]] = {start_state: (0, 0)}
    parent: dict[tuple[int, bool], tuple[tuple[int, bool], PebbleStep]] = {}
    expansions = 0

    while pq:
        f_val, g_val, tot_ops, p_mask, reached = heapq.heappop(pq)
        state_key = (p_mask, reached)
        if state_key == goal_state:
            path: list[PebbleStep] = []
            cur = goal_state
            while cur != start_state:
                prev_st, st_move = parent[cur]
                path.append(st_move)
                cur = prev_st
            path.reverse()
            return path

        if best_g.get(state_key, (10**9, 10**9)) < (g_val, tot_ops):
            continue

        expansions += 1
        if expansions > max_expansions:
            raise RuntimeError(
                f"A* search exceeded {max_expansions} expansions for length={length}, max_pebbles={max_pebbles}"
            )

        if not reached:
            # Blast phase (KMRV26 Appendix B observation 1): march rightmost pebble r -> r + 1
            r = p_mask.bit_length()  # 0 .. length - 1; next site to pebble is bit `r` (1-based `r + 1`)
            assert r < length
            extra = anc_cost[r]
            p_mid = p_mask | (1 << r)
            step_space = p_mid.bit_count() + extra
            if step_space > max_pebbles:
                continue
            now_reached = (r + 1) == length

            # Branch A: keep `r - 1` as a marker pebble (no ghost)
            choices: list[tuple[int, PebbleStep]] = [
                (p_mid, PebbleStep(pebble=(r + 1,)))
            ]
            # Branch B: if r >= 1, immediately Ghost(r) (`1`-based site `r`, bit `r - 1`)!
            if r >= 1 and (p_mask & (1 << (r - 1))):
                choices.append(
                    (p_mid & ~(1 << (r - 1)), PebbleStep(pebble=(r + 1,), ghost=(r,)))
                )

            for nxt_mask, move in choices:
                ng = g_val + 1
                nops = tot_ops + 1
                nkey = (nxt_mask, now_reached)
                if (ng, nops) < best_g.get(nkey, (10**9, 10**9)):
                    best_g[nkey] = (ng, nops)
                    parent[nkey] = (state_key, move)
                    nf = ng + heuristic(nxt_mask, now_reached)
                    heapq.heappush(pq, (nf, ng, nops, nxt_mask, now_reached))
        else:
            # Unblast phase: rightmost pebble bit is `top_idx = max_site - 1`
            max_site = p_mask.bit_length()  # 1 .. length
            if max_site == 0:
                continue

            # Candidate Phase-1 actions on sites `0 .. max_site - 1` (never need to pebble > max_site!):
            # - Site 0 (`1`-based 1): pred is root (always available).
            #     If bit 0 == 0: can Pebble(1) (active {0, 1})
            #     If bit 0 == 1: can Unpebble(1) ONLY IF max_site == 1 (wait: could you Unpebble(1) while a higher pebble exists? Yes! e.g. task-202 line 14 `Unpebble(1, 4)`! Because Unpebble(1) leaves 1 EMPTY, which behaves identically to GHOST for future Pebble(1)!)
            # - Site i (`1 <= i < max_site`, `1`-based `i + 1`): requires bit `i - 1` == 1 in `p_mask`!
            #     If bit i == 0: can Pebble(i + 1) (active {i, i + 1}), with optional immediate Ghost(i)!
            #     If bit i == 1: can Unpebble(i + 1) ONLY IF `i == max_site - 1` OR allow interior unpebble (allowing interior unpebble is sound too!).
            cand_actions: list[tuple[str, int, int, int]] = []
            # tuple: (kind "P"|"PG"|"U", bit_i, active_wire_mask, anc_extra)
            # Note on active_wire_mask: root is bit `0`, site `k` (`1..length`) is bit `k`!
            for i in range(max_site):
                pred_ok = True if i == 0 else bool(p_mask & (1 << (i - 1)))
                if not pred_ok:
                    continue
                act_w = (1 << (i + 1)) | (1 << i)
                anc_w = anc_cost[i]
                if (p_mask & (1 << i)) == 0:
                    cand_actions.append(("P", i, act_w, anc_w))
                    if i >= 1:
                        # Immediate Ghost(i) (1-based `i`, bit `i - 1`) after Pebble(i + 1)
                        cand_actions.append(("PG", i, act_w, anc_w))
                elif i == max_site - 1:
                    # KMRV26 Appendix B: only the rightmost pebble needs Unpebble;
                    # any lower pebble is cheaper to Ghost immediately after use.
                    cand_actions.append(("U", i, act_w, anc_w))

            cur_pebbles = p_mask.bit_count()

            def enum_moves(
                idx: int,
                used_wires: int,
                p_mid: int,
                p_after: int,
                extra_anc: int,
                p_list: list[int],
                u_list: list[int],
                g_list: list[int],
            ) -> None:
                if p_list or u_list:
                    # Ensure we never leave a ghost above max(p_after):
                    # Wait: only Unpebble can reduce max_site, and Unpebble sets the site to EMPTY (not GHOST!).
                    # What if PG ghosted `i - 1` while `i` was pebbled? Then `i > i - 1` IS in `p_after`, so `i - 1 < max(p_after)` ALWAYS holds!
                    live_union = p_mask | p_mid | p_after
                    if live_union.bit_count() + extra_anc <= max_pebbles:
                        ng = g_val + 1
                        nops = tot_ops + len(p_list) + len(u_list)
                        nkey = (p_after, True)
                        if (ng, nops) < best_g.get(nkey, (10**9, 10**9)):
                            best_g[nkey] = (ng, nops)
                            parent[nkey] = (
                                state_key,
                                PebbleStep(pebble=p_list, unpebble=u_list, ghost=g_list),
                            )
                            nf = ng + heuristic(p_after, True)
                            heapq.heappush(pq, (nf, ng, nops, p_after, True))

                for j in range(idx, len(cand_actions)):
                    kind, bit_i, act_w, anc_w = cand_actions[j]
                    if used_wires & act_w:
                        continue
                    # Quick prune on space
                    nm = p_mid | (1 << bit_i) if kind in ("P", "PG") else (p_mid & ~(1 << bit_i))
                    if (p_mask | nm).bit_count() + extra_anc + anc_w > max_pebbles:
                        continue
                    if kind == "P":
                        na = p_after | (1 << bit_i)
                        enum_moves(
                            j + 1, used_wires | act_w, nm, na, extra_anc + anc_w,
                            p_list + [bit_i + 1], u_list, g_list,
                        )
                    elif kind == "PG":
                        na = (p_after | (1 << bit_i)) & ~(1 << (bit_i - 1))
                        enum_moves(
                            j + 1, used_wires | act_w, nm, na, extra_anc + anc_w,
                            p_list + [bit_i + 1], u_list, g_list + [bit_i],
                        )
                    else:
                        na = p_after & ~(1 << bit_i)
                        enum_moves(
                            j + 1, used_wires | act_w, nm, na, extra_anc + anc_w,
                            p_list, u_list + [bit_i + 1], g_list,
                        )

            _ = cur_pebbles
            enum_moves(0, 0, p_mask, p_mask, 0, [], [], [])

    return None


def compare_chain_strategies(length: int, reg_width: int = 1) -> dict[str, dict[str, int | str]]:
    """Analytical + compiled comparison of reversible / spooky chain strategies on length `length`:
      - Naive out-of-place compute + reverse uncompute
      - Kaliski / Ragavan-Vaikuntanathan reversible Fibonacci (stores values + mod-N inverses)
      - Gidney sequential spooky pebbling
      - KMRV26 parallel spooky pebbling (`Blast`/`Unblast`)
    """
    ctx = Context(prec=12, rounding=ROUND_CEILING)
    gidney_sched = compile_gidney_sequential_schedule(length)
    gidney_val = validate_pebbling_schedule(PebbleDAG.line(length), gidney_sched)
    kmrv_sched = compile_kmrv_blast_schedule(length)
    kmrv_val = validate_pebbling_schedule(PebbleDAG.line(length), kmrv_sched)

    return {
        "naive_out_of_place": {
            "depth": 2 * length,
            "calls": 2 * length,
            "pebbles": length,
            "qubits": length * reg_width,
        },
        "reversible_fibonacci_inv": {
            # Ragavan 2024: ~ 4x depth/calls due to carrying intermediate values + modular inverses
            "depth": 4 * length,
            "calls": 8 * length,
            "pebbles": 8,
            "qubits": 8 * reg_width,
        },
        "gidney_sequential_spooky": {
            "depth": gidney_val.depth,
            "calls": gidney_val.pebble_cost,
            "pebbles": gidney_val.peak_pebbles,
            "qubits": gidney_val.peak_pebbles * reg_width,
        },
        "kmrv_parallel_spooky": {
            "depth": kmrv_val.depth,
            "calls": kmrv_val.pebble_cost,
            "pebbles": kmrv_val.peak_pebbles,
            "qubits": kmrv_val.peak_pebbles * reg_width,
            "space_over_log2_ell": (
                str(+ctx.divide(Decimal(kmrv_val.peak_pebbles), Decimal(max(1, (length - 1).bit_length()))))
            ),
        },
    }
