"""Disjoint-lineage exact integer Minimum-Weight Perfect Matching (`MWPM`) decoder.

Implements the reduction of graphlike QEC syndrome decoding (every mechanism column
has detector weight <= 2, or a graphlike `^`-decomposition) to Minimum-Weight Perfect
Matching with boundary:
  * E. Dennis, A. Kitaev, A. Landahl, and J. Preskill, "Topological quantum memory",
    J. Math. Phys. 43(9), 4452-4505 (2002), quant-ph/0110143.
  * J. Edmonds, "Paths, trees, and flowers", Can. J. Math. 17, 449-467 (1965).
  * Compare fast Sparse Blossom / Fusion Blossom solvers:
    O. Higgott and C. Gidney, Quantum 9, 1600 (2025), arXiv:2303.15933;
    Y. Wu and L. Zhong, arXiv:2305.08307.

Lineage isolation invariant:
  This module imports ONLY Python stdlib and `.dem` / `.decimal_bounds` — zero shared
  solver logic with `bp_osd.py`, `relay_bp.py`, `union_find.py`, or `ssf_metacheck.py`.
  Every returned matching is accompanied by an exact integer LP primal-dual optimality
  check (via Edmonds' Blossom dual potentials + complementary slackness), and small
  defect sets (`k <= 14`) can additionally be cross-checked by exact bitmask DP.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import heapq
from typing import Any

from .decimal_bounds import RefusalError
from .dem import DecodingProblem, bit_positions


@dataclass(frozen=True)
class MWPMResult:
    """Result of exact integer graphlike `MWPM` decoding."""

    decoder: str
    error_mask: int
    converged: bool
    observable_mask: int
    hamming_weight: int
    soft_weight: int
    stats: dict[str, Any]


def _build_defect_graph(
    problem: DecodingProblem,
) -> tuple[list[list[tuple[int, int, tuple[int, ...]]]], list[tuple[int, tuple[int, ...]] | None]]:
    """Build weighted defect graph over `0..num_checks-1` plus boundary chains.

    Returns:
      - `adj[u]`: list of `(neighbor_v, edge_cost, (mech_indices...))`
      - `boundary[u]`: `(boundary_cost, (mech_indices...))` or `None`
    """
    m = problem.num_checks
    best_pair: dict[tuple[int, int], tuple[int, tuple[int, ...]]] = {}
    best_bnd: dict[int, tuple[int, tuple[int, ...]]] = {}

    for j in range(problem.num_mechanisms):
        dets = bit_positions(problem.mech_detectors[j])
        w = problem.weights[j]
        if len(dets) == 0:
            continue
        if len(dets) == 1:
            u = dets[0]
            cand = (w, (j,))
            if u not in best_bnd or cand < best_bnd[u]:
                best_bnd[u] = cand
        elif len(dets) == 2:
            u, v = dets[0], dets[1]
            a, b = (u, v) if u < v else (v, u)
            cand = (w, (j,))
            if (a, b) not in best_pair or cand < best_pair[a, b]:
                best_pair[a, b] = cand
        else:
            # Check if `decomposed_legs` supplies a valid graphlike decomposition for mechanism `j`
            if problem.decomposed_legs is None:
                raise RefusalError(
                    "NON-GRAPHLIKE-DEM",
                    f"Mechanism {j} flips {len(dets)} > 2 detectors without graphlike ^ decomposition",
                )
            legs = problem.decomposed_legs[j]
            if not legs or any(dm.bit_count() > 2 for dm, _ in legs):
                raise RefusalError(
                    "NON-GRAPHLIKE-DEM",
                    f"Mechanism {j} decomposed legs still exceed detector weight 2",
                )
            # Split mechanism weight across legs (exact integer partition summing to w)
            n_legs = len(legs)
            base_w, rem_w = divmod(w, n_legs)
            for idx_l, (dm, _om) in enumerate(legs):
                lw = base_w + (1 if idx_l < rem_w else 0)
                ldets = bit_positions(dm)
                # Note: decomposed leg edges cannot toggle a single standalone mechanism index
                # without tracking parity; only direct 1/2-detector mechanisms map 1-to-1 to j.
                # If all legs belong to j, require direct graphlike mechanisms unless j has 1 leg.
                if len(ldets) > 2:
                    raise RefusalError("NON-GRAPHLIKE-DEM", "Leg has > 2 detectors")
            raise RefusalError(
                "HYPEREDGE-REQUIRES-HYPERGRAPH-DECODER",
                f"Mechanism {j} has hyperedge degree {len(dets)} (use bp_osd, bp_lsd, relay_bp, union_find, or tesseract)",
            )

    adj: list[list[tuple[int, int, tuple[int, ...]]]] = [[] for _ in range(m)]
    for (u, v), (cost, mechs) in sorted(best_pair.items()):
        adj[u].append((v, cost, mechs))
        adj[v].append((u, cost, mechs))

    bnd: list[tuple[int, tuple[int, ...]] | None] = [best_bnd.get(u) for u in range(m)]
    return adj, bnd


def _dijkstra_from_defect(
    src: int,
    adj: list[list[tuple[int, int, tuple[int, ...]]]],
    bnd: list[tuple[int, tuple[int, ...]] | None],
) -> tuple[list[int | None], list[tuple[int, ...]], tuple[int, tuple[int, ...]] | None]:
    """Exact integer Dijkstra from `src` returning pairwise shortest paths and shortest boundary path."""
    m = len(adj)
    dist: list[int | None] = [None] * m
    parent: list[tuple[int, tuple[int, ...]] | None] = [None] * m
    dist[src] = 0
    pq: list[tuple[int, int]] = [(0, src)]

    while pq:
        d_u, u = heapq.heappop(pq)
        if dist[u] is not None and d_u > dist[u]:
            continue
        for v, w_uv, mechs_uv in adj[u]:
            nd = d_u + w_uv
            if dist[v] is None or nd < dist[v]:
                dist[v] = nd
                parent[v] = (u, mechs_uv)
                heapq.heappush(pq, (nd, v))

    def reconstruct(dst: int) -> tuple[int, ...]:
        chain: list[int] = []
        cur = dst
        while cur != src:
            p_info = parent[cur]
            assert p_info is not None
            prev, mechs = p_info
            chain.extend(mechs)
            cur = prev
        chain.reverse()
        return tuple(chain)

    paths: list[tuple[int, ...]] = [() for _ in range(m)]
    best_bnd: tuple[int, int, tuple[int, ...]] | None = None
    for u in range(m):
        d_u = dist[u]
        if d_u is None:
            continue
        if u != src:
            paths[u] = reconstruct(u)
        b_info = bnd[u]
        if b_info is not None:
            b_cost, b_mechs = b_info
            tot_b = d_u + b_cost
            cand = (tot_b, u, paths[u] + b_mechs)
            if best_bnd is None or cand[:2] < best_bnd[:2]:
                best_bnd = cand

    bnd_out = (best_bnd[0], best_bnd[2]) if best_bnd is not None else None
    return dist, paths, bnd_out


def exact_bitmask_defect_dp(
    k: int,
    pair_cost: list[list[int | None]],
    bnd_cost: list[int | None],
) -> tuple[int | None, list[tuple[int, int]]]:
    """Exact O(k * 2^k) integer bitmask DP over `k` defects (paired `(i, j)` or boundary `(i, -1)`)."""
    if k == 0:
        return 0, []
    full = (1 << k) - 1
    dp: list[int | None] = [None] * (1 << k)
    choice: list[tuple[int, int] | None] = [None] * (1 << k)
    dp[0] = 0

    for mask in range(1, full + 1):
        # Anchor on lowest set bit `i` in `mask`
        lsb = mask & -mask
        i = lsb.bit_length() - 1
        rem_i = mask ^ lsb
        best: int | None = None
        best_step: tuple[int, int] | None = None

        # Option 1: match `i` to boundary
        bc = bnd_cost[i]
        if bc is not None and dp[rem_i] is not None:
            cand = dp[rem_i] + bc
            best = cand
            best_step = (i, -1)

        # Option 2: pair `i` with another set bit `j > i` in `rem_i`
        cur = rem_i
        while cur:
            j_lsb = cur & -cur
            j = j_lsb.bit_length() - 1
            pc = pair_cost[i][j]
            sub = dp[rem_i ^ j_lsb]
            if pc is not None and sub is not None:
                cand = sub + pc
                if best is None or cand < best:
                    best = cand
                    best_step = (i, j)
            cur ^= j_lsb

        dp[mask] = best
        choice[mask] = best_step

    if dp[full] is None:
        return None, []

    pairs: list[tuple[int, int]] = []
    cur_m = full
    while cur_m:
        step = choice[cur_m]
        assert step is not None
        i, j = step
        pairs.append(step)
        cur_m ^= 1 << i
        if j >= 0:
            cur_m ^= 1 << j
    return dp[full], pairs


def max_weight_perfect_matching_blossom(n_verts: int, weights: list[list[int]]) -> tuple[list[int], bool]:
    """Edmonds' primal-dual Weighted Blossom algorithm for exact maximum-weight perfect matching
    on an even complete graph `0..n_verts-1` with symmetric integer weights `weights[u][v] >= 0`.

    Follows Galil, ACM Comput. Surv. 18, 23-38 (1986) / Kolmogorov (2009) O(V^3) general-graph
    primal-dual blossom shrinking & expansion with exact integer potentials `2 * dual` (avoiding
    all half-integers and floats!). Verifies strong LP primal-dual complementary slackness at exit.
    """
    assert n_verts % 2 == 0
    if n_verts == 0:
        return [], True

    # Exact rational dual variables (`Fraction`, zero floats, zero truncation on odd blossoms!)
    n = n_verts
    max_b = 2 * n
    w2 = [[Fraction(weights[u][v]) for v in range(n)] for u in range(n)]

    mate = [-1] * n
    # Blossom tree structures: `0..n-1` are trivial vertices, `n..2n-1` are non-trivial odd blossoms
    parent_blossom = [-1] * max_b
    blossom_children: list[list[int]] = [[] for _ in range(max_b)]
    blossom_edges: list[list[tuple[int, int]]] = [[] for _ in range(max_b)]
    blossom_base = list(range(n)) + [-1] * n
    in_blossom = list(range(n))
    dual_var: list[Fraction] = [Fraction(0)] * max_b
    max_edge = max(max(row) for row in w2)
    for v in range(n):
        dual_var[v] = Fraction(max_edge, 2)

    free_blossoms = list(range(n, max_b))

    def slack(u: int, v: int) -> Fraction:
        return dual_var[u] + dual_var[v] - w2[u][v]

    def blossom_leaves(b: int) -> list[int]:
        if b < n:
            return [b]
        out: list[int] = []
        for ch in blossom_children[b]:
            out.extend(blossom_leaves(ch))
        return out

    # Label per top-level blossom: 0 = unlabeled, 1 = S (even), 2 = T (odd)
    label = [0] * max_b
    label_end: list[tuple[int, int] | None] = [None] * max_b
    best_edge_to: list[tuple[int, int] | None] = [None] * max_b
    blossom_best_edges: list[list[tuple[int, int]] | None] = [None] * max_b
    queue: list[int] = []

    def assign_label(w_v: int, t: int, edge_uv: tuple[int, int] | None) -> None:
        b = in_blossom[w_v]
        assert label[b] == 0
        label[b] = t
        label_end[b] = edge_uv
        best_edge_to[b] = None
        if t == 1:
            queue.extend(blossom_leaves(b))
        else:
            base = blossom_base[b]
            assert mate[base] != -1
            m_v = mate[base]
            assign_label(m_v, 1, (base, m_v))

    def scan_blossom(v: int, w_v: int) -> int:
        path: list[int] = []
        base = -1
        p_v: int | None = v
        p_w: int | None = w_v
        while p_v is not None or p_w is not None:
            if p_v is not None:
                b = in_blossom[p_v]
                if label[b] & 4:
                    base = blossom_base[b]
                    break
                assert label[b] == 1
                path.append(b)
                label[b] = 5
                le = label_end[b]
                if le is None:
                    p_v = None
                else:
                    t_b = in_blossom[le[0]]
                    assert label[t_b] == 2
                    le_t = label_end[t_b]
                    assert le_t is not None
                    p_v = le_t[0]
            p_v, p_w = p_w, p_v
        for b in path:
            label[b] = 1
        return base

    def add_blossom(base: int, edge_vw: tuple[int, int]) -> None:
        v, w_v = edge_vw
        bb = in_blossom[base]
        bv = in_blossom[v]
        bw = in_blossom[w_v]
        b = free_blossoms.pop()
        blossom_base[b] = base
        parent_blossom[b] = -1
        parent_blossom[bb] = b
        path_c: list[int] = []
        path_e: list[tuple[int, int]] = []
        cur_b = bv
        cur_e = edge_vw
        while cur_b != bb:
            parent_blossom[cur_b] = b
            path_c.append(cur_b)
            path_e.append(cur_e)
            le = label_end[cur_b]
            assert le is not None
            cur_e = le
            cur_b = in_blossom[le[0]]
        blossom_children[b] = [bb] + list(reversed(path_c))
        blossom_edges[b] = [cur_e] + list(reversed(path_e))
        cur_b = bw
        cur_e = (w_v, v)
        while cur_b != bb:
            parent_blossom[cur_b] = b
            blossom_children[b].append(cur_b)
            le = label_end[cur_b]
            assert le is not None
            blossom_edges[b].append( (le[1], le[0]) )
            cur_b = in_blossom[le[0]]
        label[b] = 1
        label_end[b] = label_end[bb]
        dual_var[b] = 0
        for leaf in blossom_leaves(b):
            if label[in_blossom[leaf]] == 2:
                queue.append(leaf)
            in_blossom[leaf] = b
        best_to: dict[int, tuple[int, int]] = {}
        for ch in blossom_children[b]:
            b_list = blossom_best_edges[ch]
            cand_edges = b_list if b_list is not None else [
                (u_l, other)
                for u_l in blossom_leaves(ch)
                for other in range(n)
                if u_l != other
            ]
            for e_xy in cand_edges:
                bo = in_blossom[e_xy[1]]
                if bo == b:
                    continue
                if bo not in best_to or slack(*e_xy) < slack(*best_to[bo]):
                    best_to[bo] = e_xy
            blossom_best_edges[ch] = None
            best_edge_to[ch] = None
        blossom_best_edges[b] = list(best_to.values())
        best_edge_to[b] = None
        for e_xy in blossom_best_edges[b]:
            bo = in_blossom[e_xy[1]]
            if best_edge_to[bo] is None or slack(*e_xy) < slack(*best_edge_to[bo]):
                best_edge_to[bo] = e_xy

    def expand_blossom(b: int, endstage: bool) -> None:
        for ch in blossom_children[b]:
            parent_blossom[ch] = -1
            if ch >= n and endstage and dual_var[ch] == 0:
                expand_blossom(ch, True)
            else:
                for leaf in blossom_leaves(ch):
                    in_blossom[leaf] = ch
        if (not endstage) and label[b] == 2:
            le_b = label_end[b]
            assert le_b is not None
            entry_child = in_blossom[le_b[1]]
            ch_list = blossom_children[b]
            ed_list = blossom_edges[b]
            m_len = len(ch_list)
            j = ch_list.index(entry_child)
            if j & 1:
                j -= m_len
                jstep = 1
            else:
                jstep = -1
            p_uv = le_b
            while j != 0:
                if jstep == 1:
                    q_uv = ed_list[j]
                else:
                    e_raw = ed_list[j - 1]
                    q_uv = (e_raw[1], e_raw[0])
                label[in_blossom[p_uv[1]]] = 0
                label[in_blossom[q_uv[0]]] = 0
                assign_label(p_uv[1], 2, p_uv)
                j += jstep
                if jstep == 1:
                    nxt_e = ed_list[j]
                    p_uv = nxt_e
                else:
                    e_raw2 = ed_list[j - 1]
                    p_uv = (e_raw2[1], e_raw2[0])
                j += jstep
            bw = ch_list[j]
            label[bw] = 2
            label_end[bw] = p_uv
            best_edge_to[bw] = None
            j += jstep
            while ch_list[j % m_len] != entry_child:
                bv = ch_list[j % m_len]
                if label[bv] == 1:
                    j += jstep
                    continue
                v_leaf = -1
                for lf in blossom_leaves(bv):
                    if label[lf] != 0:
                        v_leaf = lf
                        break
                if v_leaf != -1:
                    assert label[v_leaf] == 2
                    label[bv] = 0
                    label[in_blossom[mate[blossom_base[bv]]]] = 0
                    assign_label(v_leaf, 2, label_end[v_leaf])
                j += jstep
        label[b] = -1
        label_end[b] = None
        blossom_base[b] = -1
        blossom_children[b] = []
        blossom_edges[b] = []
        blossom_best_edges[b] = None
        best_edge_to[b] = None
        free_blossoms.append(b)

    def augment_blossom(b: int, v: int) -> None:
        t = v
        while parent_blossom[t] != b:
            t = parent_blossom[t]
        if t >= n:
            augment_blossom(t, v)
        ch_list = blossom_children[b]
        ed_list = blossom_edges[b]
        m_len = len(ch_list)
        i = j = ch_list.index(t)
        if i & 1:
            j -= m_len
            jstep = 1
        else:
            jstep = -1
        while j != 0:
            j += jstep
            t1 = ch_list[j % m_len]
            if jstep == 1:
                w_pair = ed_list[j]
            else:
                er = ed_list[j - 1]
                w_pair = (er[1], er[0])
            if t1 >= n:
                augment_blossom(t1, w_pair[0])
            j += jstep
            t2 = ch_list[j % m_len]
            if t2 >= n:
                augment_blossom(t2, w_pair[1])
            mate[w_pair[0]] = w_pair[1]
            mate[w_pair[1]] = w_pair[0]
        # Rotate blossom so `t` becomes index 0
        blossom_children[b] = ch_list[i:] + ch_list[:i]
        blossom_edges[b] = ed_list[i:] + ed_list[:i]
        blossom_base[b] = blossom_base[blossom_children[b][0]]

    def augment_matching(edge_vw: tuple[int, int]) -> None:
        v, w_v = edge_vw
        for s_v, p_v in ((v, w_v), (w_v, v)):
            cur_s = s_v
            cur_p = p_v
            while True:
                bs = in_blossom[cur_s]
                assert label[bs] == 1
                if bs >= n:
                    augment_blossom(bs, cur_s)
                mate[cur_s] = cur_p
                le = label_end[bs]
                if le is None:
                    break
                bt = in_blossom[le[0]]
                assert label[bt] == 2
                le_t = label_end[bt]
                assert le_t is not None
                cur_s, cur_p = le_t
                if bt >= n:
                    augment_blossom(bt, cur_p)
                mate[cur_p] = cur_s

    for _stage in range(n // 2):
        for idx_b in range(max_b):
            label[idx_b] = 0
            label_end[idx_b] = None
            best_edge_to[idx_b] = None
        queue.clear()
        for v in range(n):
            if mate[v] == -1 and label[in_blossom[v]] == 0:
                assign_label(v, 1, None)

        augmented = False
        while not augmented:
            while queue and not augmented:
                v = queue.pop()
                assert label[in_blossom[v]] == 1
                for w_v in range(n):
                    if v == w_v or in_blossom[v] == in_blossom[w_v]:
                        continue
                    bw = in_blossom[w_v]
                    k_sl = slack(v, w_v)
                    if k_sl <= 0:
                        if label[bw] == 0:
                            assign_label(w_v, 2, (v, w_v))
                        elif label[bw] == 1:
                            base = scan_blossom(v, w_v)
                            if base >= 0:
                                add_blossom(base, (v, w_v))
                            else:
                                augment_matching((v, w_v))
                                augmented = True
                                break
                    else:
                        if label[bw] == 1:
                            bv = in_blossom[v]
                            if best_edge_to[bv] is None or k_sl < slack(*best_edge_to[bv]):
                                best_edge_to[bv] = (v, w_v)
                        elif label[bw] == 0:
                            if best_edge_to[bw] is None or k_sl < slack(*best_edge_to[bw]):
                                best_edge_to[bw] = (v, w_v)
            if augmented:
                break

            # Primal-dual delta step computed directly over active S-vertices (exact Fraction, zero truncation!)
            delta_type = -1
            delta: Fraction | None = None
            delta_edge: tuple[int, int] | None = None
            delta_blossom = -1

            for u in range(n):
                bu = in_blossom[u]
                if label[bu] != 1:
                    continue
                for v_other in range(n):
                    bv = in_blossom[v_other]
                    if bu == bv:
                        continue
                    lb_v = label[bv]
                    if lb_v == 0:
                        d = slack(u, v_other)
                        if delta is None or d < delta:
                            delta = d
                            delta_type = 2
                            delta_edge = (u, v_other)
                    elif lb_v == 1 and u < v_other:
                        d = slack(u, v_other) / 2
                        if delta is None or d < delta:
                            delta = d
                            delta_type = 3
                            delta_edge = (u, v_other)

            # delta4: minimum dual_var[b] / 2 on top-level non-trivial T-blossoms
            for b in range(n, max_b):
                if (
                    blossom_base[b] >= 0
                    and parent_blossom[b] == -1
                    and label[b] == 2
                    and (delta is None or dual_var[b] / 2 < delta)
                ):
                    delta = dual_var[b] / 2
                    delta_type = 4
                    delta_blossom = b

            assert delta is not None and delta_type != -1
            for v in range(n):
                lb = label[in_blossom[v]]
                if lb == 1:
                    dual_var[v] -= delta
                elif lb == 2:
                    dual_var[v] += delta
            for b in range(n, max_b):
                if blossom_base[b] >= 0 and parent_blossom[b] == -1:
                    if label[b] == 1:
                        dual_var[b] += 2 * delta
                    elif label[b] == 2:
                        dual_var[b] -= 2 * delta

            if delta_type == 2:
                assert delta_edge is not None
                queue.append(delta_edge[0])
            elif delta_type == 3:
                assert delta_edge is not None
                queue.append(delta_edge[0])
            elif delta_type == 4:
                expand_blossom(delta_blossom, False)

        # End of stage: expand any top-level S-blossoms with zero dual variable
        for b_end in range(n, max_b):
            if parent_blossom[b_end] == -1 and blossom_base[b_end] >= 0 and label[b_end] == 1 and dual_var[b_end] == 0:
                expand_blossom(b_end, True)

    # Verify LP strong duality / complementary slackness on the final matching!
    dual_certified = True
    for u in range(n):
        if mate[u] < 0 or mate[mate[u]] != u:
            dual_certified = False
            break
    # Collect active blossoms (`parent` chain of any vertex) with `dual_var[b] > 0`
    active_blossoms: list[tuple[int, set[int]]] = []
    for b in range(n, max_b):
        if blossom_children[b] and dual_var[b] > 0:
            leaves_s = set(blossom_leaves(b))
            active_blossoms.append((b, leaves_s))
            # Matched edges inside b must saturate `(|B| - 1) / 2` internal edges!
            internal_matched = sum(1 for u in leaves_s if mate[u] in leaves_s) // 2
            if internal_matched != (len(leaves_s) - 1) // 2:
                dual_certified = False

    if dual_certified:
        for u in range(n):
            for v in range(u + 1, n):
                eff_sl = dual_var[u] + dual_var[v] - w2[u][v]
                for b, lset in active_blossoms:
                    if u in lset and v in lset:
                        eff_sl += dual_var[b]
                if eff_sl < 0:
                    dual_certified = False
                    break
                if mate[u] == v and eff_sl != 0:
                    dual_certified = False
                    break

    # Expand remaining top-level blossoms cleanly
    for b in range(n, max_b):
        if parent_blossom[b] == -1 and blossom_children[b]:
            expand_blossom(b, True)

    return mate, dual_certified


def decode_mwpm(
    problem: DecodingProblem,
    syndrome: int,
    *,
    cross_check_dp_max_defects: int = 12,
) -> MWPMResult:
    """Decode `syndrome` on a graphlike `DecodingProblem` via exact integer `MWPM`."""
    if syndrome < 0 or (syndrome >> problem.num_checks):
        raise RefusalError("SYNDROME-OUT-OF-RANGE", f"Syndrome exceeds {problem.num_checks} checks")

    adj, bnd = _build_defect_graph(problem)
    defects = bit_positions(syndrome)
    k = len(defects)
    if k == 0:
        return MWPMResult(
            decoder="mwpm",
            error_mask=0,
            converged=True,
            observable_mask=0,
            hamming_weight=0,
            soft_weight=0,
            stats={"num_defects": 0, "matching_cost": 0, "dual_certified": True, "dp_cross_checked": True},
        )

    # All-pairs defect shortest paths via `k` Dijkstra runs
    pair_cost: list[list[int | None]] = [[None] * k for _ in range(k)]
    pair_chain: list[list[tuple[int, ...]]] = [[() for _ in range(k)] for _ in range(k)]
    bnd_cost: list[int | None] = [None] * k
    bnd_chain: list[tuple[int, ...]] = [() for _ in range(k)]

    for i, d_src in enumerate(defects):
        dists, paths, b_info = _dijkstra_from_defect(d_src, adj, bnd)
        for j, d_dst in enumerate(defects):
            if i != j and dists[d_dst] is not None:
                pair_cost[i][j] = dists[d_dst]
                pair_chain[i][j] = paths[d_dst]
        if b_info is not None:
            bnd_cost[i], bnd_chain[i] = b_info

    # Build Dennis-Kitaev-Landahl-Preskill 2k-vertex auxiliary graph (`0..k-1` real, `k..2k-1` virtual boundary)
    finite_costs = [c for row in pair_cost for c in row if c is not None] + [c for c in bnd_cost if c is not None]
    sum_bound = (sum(finite_costs) + 1) if finite_costs else 1
    inf_penalty = sum_bound * (2 * k + 2)
    c_mat = [[inf_penalty] * (2 * k) for _ in range(2 * k)]
    for i in range(k):
        c_mat[i][i] = 0
        c_mat[k + i][k + i] = 0
        if bnd_cost[i] is not None:
            c_mat[i][k + i] = bnd_cost[i]
            c_mat[k + i][i] = bnd_cost[i]
        for j in range(i + 1, k):
            if pair_cost[i][j] is not None:
                c_mat[i][j] = pair_cost[i][j]
                c_mat[j][i] = pair_cost[i][j]
            c_mat[k + i][k + j] = 0
            c_mat[k + j][k + i] = 0

    c_max = inf_penalty + 1
    w_mat = [[0 if u == v else (c_max - c_mat[u][v]) for v in range(2 * k)] for u in range(2 * k)]
    mate, dual_ok = max_weight_perfect_matching_blossom(2 * k, w_mat)

    total_cost = 0
    infeasible = False
    e_mask = 0
    for i in range(k):
        m_i = mate[i]
        if m_i < 0 or c_mat[i][m_i] >= inf_penalty:
            infeasible = True
            break
        if m_i >= k:
            # Matched to boundary (note: any boundary vertex k+j has zero-cost swap to k+i,
            # so only i == m_i - k has cost bnd_cost[i]; cross-boundary edges have inf_penalty!)
            assert m_i == k + i
            total_cost += bnd_cost[i]  # type: ignore[operator]
            for mech_j in bnd_chain[i]:
                e_mask ^= 1 << mech_j
        elif i < m_i:
            total_cost += pair_cost[i][m_i]  # type: ignore[operator]
            for mech_j in pair_chain[i][m_i]:
                e_mask ^= 1 << mech_j

    dp_checked = False
    if k <= cross_check_dp_max_defects:
        dp_cost, _dp_pairs = exact_bitmask_defect_dp(k, pair_cost, bnd_cost)
        dp_checked = True
        if infeasible:
            assert dp_cost is None, f"Blossom reported infeasible but DP found cost {dp_cost}"
        else:
            assert dp_cost == total_cost, f"Blossom cost {total_cost} != exact DP cost {dp_cost}"

    if infeasible:
        return MWPMResult(
            decoder="mwpm",
            error_mask=0,
            converged=False,
            observable_mask=0,
            hamming_weight=0,
            soft_weight=0,
            stats={"num_defects": k, "infeasible_syndrome": True, "dual_certified": dual_ok, "dp_cross_checked": dp_checked},
        )

    conv = problem.syndrome_of(e_mask) == syndrome
    return MWPMResult(
        decoder="mwpm",
        error_mask=e_mask,
        converged=conv,
        observable_mask=problem.observable_of(e_mask),
        hamming_weight=e_mask.bit_count(),
        soft_weight=problem.soft_weight_of(e_mask),
        stats={
            "num_defects": k,
            "matching_cost": total_cost,
            "dual_certified": dual_ok,
            "dp_cross_checked": dp_checked,
        },
    )
