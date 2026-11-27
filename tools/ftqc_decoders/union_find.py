"""Peeling, Hypergraph-`GF(2)` Union-Find (`UF`), and `BeliefFind` syndrome decoders.

Literature & Algorithmic Lineage:
  1. Linear-time spanning-forest peeling on erasure/cluster graphs:
     N. Delfosse and G. Zemor, "Linear-time maximum likelihood decoding of surface
     codes over the quantum erasure channel", Phys. Rev. Research 2, 033042 (2020),
     arXiv:1703.01517.
  2. Almost-linear Union-Find topological decoder:
     N. Delfosse and N. H. Nickerson, "Almost-linear time decoding algorithm for
     topological codes", Quantum 5, 595 (2021), arXiv:1709.06218.
  3. Hypergraph / qLDPC Union-Find with local cluster Gaussian elimination:
     N. Delfosse, V. Londe, and M. E. Beverland, "Toward a Union-Find decoder for
     quantum LDPC codes", IEEE Trans. Inf. Theory 68(5), 3187-3199 (2022),
     arXiv:2103.08049.
  4. `BeliefFind` (BP-guided Union-Find cluster growth):
     J. Roffe et al., `quantumgizmos/ldpc` (`BeliefFindDecoder`), arXiv:2005.07016.

Zero `float` or `complex` literals/calls appear anywhere in this module.
"""

from __future__ import annotations

from collections import deque

from .bp_osd import DecodeResult, run_bp_nms
from .decimal_bounds import RefusalError
from .dem import DecodingProblem, bit_positions


def _peel_graphlike_cluster(
    cluster_checks: set[int],
    cluster_mechs: list[int],
    mech_detectors: tuple[int, ...],
    syndrome: int,
) -> tuple[int, bool]:
    """Linear-time spanning-tree peeling on a valid graphlike (`wt(col) <= 2`) cluster.

    Follows Delfosse & Zemor (Phys. Rev. Research 2, 033042, 2020):
      - Builds a spanning tree over `cluster_checks`, rooting at a check with a
        boundary edge (`wt == 1`) when one exists (so any residual root defect exits
        cleanly into the boundary edge!).
      - Peels leaves bottom-up in reverse BFS order.
    """
    if not cluster_checks:
        return 0, True

    adj: dict[int, list[tuple[int, int]]] = {c: [] for c in cluster_checks}
    bnd_edge: dict[int, int] = {}

    for mj in cluster_mechs:
        dets = [d for d in bit_positions(mech_detectors[mj]) if d in cluster_checks]
        if len(dets) == 1:
            u = dets[0]
            if u not in bnd_edge:
                bnd_edge[u] = mj
        elif len(dets) == 2:
            u, v = dets[0], dets[1]
            adj[u].append((v, mj))
            adj[v].append((u, mj))
        elif len(dets) > 2:
            return 0, False

    # Choose root: prefer a vertex with a boundary mechanism if available
    sorted_checks = sorted(cluster_checks)
    root = next((c for c in sorted_checks if c in bnd_edge), sorted_checks[0])

    # BFS spanning tree from `root`
    visited = {root}
    bfs_order = [root]
    parent_edge: dict[int, tuple[int, int]] = {}  # child -> (parent_check, mech_j)
    q: deque[int] = deque([root])
    while q:
        u = q.popleft()
        for v, mj in adj[u]:
            if v not in visited:
                visited.add(v)
                parent_edge[v] = (u, mj)
                bfs_order.append(v)
                q.append(v)

    if len(visited) != len(cluster_checks):
        return 0, False

    defect = {c: ((syndrome >> c) & 1) for c in cluster_checks}
    e_mask = 0
    for u in reversed(bfs_order[1:]):
        if defect[u]:
            p_u, mj = parent_edge[u]
            e_mask ^= 1 << mj
            defect[u] ^= 1
            defect[p_u] ^= 1

    if defect[root]:
        if root not in bnd_edge:
            return 0, False
        e_mask ^= 1 << bnd_edge[root]
        defect[root] = 0

    return e_mask, True


def decode_union_find(
    problem: DecodingProblem,
    syndrome: int,
    *,
    use_bp_guidance: bool = False,
    bp_iters: int = 15,
) -> DecodeResult:
    """Decode `syndrome` via `Union-Find` (`use_bp_guidance=False`) or `BeliefFind`
    (`use_bp_guidance=True`).

    * On graphlike clusters (`wt(col) <= 2`), checks cluster validity in `O(1)` via
      parity/boundary flags and recovers via linear-time spanning-tree peeling!
    * On hypergraph clusters (`wt(col) >= 3`, Delfosse-Londe-Beverland 2022), falls back
      cleanly to local per-cluster `GF(2)` column-echelon elimination.
    """
    if syndrome < 0 or (syndrome >> problem.num_checks):
        raise RefusalError("SYNDROME-OUT-OF-RANGE", f"Syndrome exceeds {problem.num_checks} checks")

    dec_name = "belief_find" if use_bp_guidance else "union_find"
    if syndrome == 0:
        return DecodeResult(
            decoder=dec_name,
            error_mask=0,
            converged=True,
            observable_mask=0,
            hamming_weight=0,
            soft_weight=0,
            stats={
                "bp_guided": use_bp_guidance,
                "num_clusters": 0,
                "peeled_clusters": 0,
                "gf2_clusters": 0,
                "max_cluster_checks": 0,
            },
        )

    prior_rank = list(problem.weights)
    bp_used = 0
    if use_bp_guidance:
        e_bp, conv_bp, scores, bp_used = run_bp_nms(problem, syndrome, max_iters=bp_iters, schedule="serial")
        if conv_bp:
            return DecodeResult(
                decoder=dec_name,
                error_mask=e_bp,
                converged=True,
                observable_mask=problem.observable_of(e_bp),
                hamming_weight=e_bp.bit_count(),
                soft_weight=problem.soft_weight_of(e_bp),
                stats={
                    "bp_guided": True,
                    "bp_iterations": bp_used,
                    "bp_converged_directly": True,
                    "num_clusters": 0,
                    "peeled_clusters": 0,
                    "gf2_clusters": 0,
                    "max_cluster_checks": 0,
                },
            )
        prior_rank = scores

    graphlike = problem.is_graphlike()
    m = problem.num_checks
    check_vars = [bit_positions(r) for r in problem.check_rows]

    parent_c = [-1] * m
    c_checks: dict[int, set[int]] = {}
    c_mechs: dict[int, list[int]] = {}
    c_parity: dict[int, int] = {}
    c_has_bnd: dict[int, bool] = {}
    c_syn_mask: dict[int, int] = {}
    c_basis: dict[int, dict[int, tuple[int, int]]] = {}
    used_mech = [False] * problem.num_mechanisms

    def find_c(c: int) -> int:
        r = c
        while parent_c[r] != r:
            r = parent_c[r]
        cur = c
        while cur != r:
            nxt = parent_c[cur]
            parent_c[cur] = r
            cur = nxt
        return r

    def add_basis_col(root: int, mj: int) -> None:
        v = problem.mech_detectors[mj]
        comb_m = 1 << mj
        bmap = c_basis[root]
        for pb in sorted(bmap.keys(), reverse=True):
            if (v >> pb) & 1:
                bv, bm = bmap[pb]
                v ^= bv
                comb_m ^= bm
        if v != 0:
            npiv = v.bit_length() - 1
            for pb in list(bmap.keys()):
                bv, bm = bmap[pb]
                if (bv >> npiv) & 1:
                    bmap[pb] = (bv ^ v, bm ^ comb_m)
            bmap[npiv] = (v, comb_m)

    def ensure_cluster(c: int) -> int:
        if parent_c[c] != -1:
            return find_c(c)
        parent_c[c] = c
        bit = (syndrome >> c) & 1
        c_checks[c] = {c}
        c_mechs[c] = []
        c_parity[c] = bit
        c_has_bnd[c] = False
        c_syn_mask[c] = (1 << c) if bit else 0
        c_basis[c] = {}
        return c

    def merge_c(r1: int, r2: int) -> int:
        if r1 == r2:
            return r1
        if len(c_checks[r1]) < len(c_checks[r2]):
            r1, r2 = r2, r1
        parent_c[r2] = r1
        c_checks[r1].update(c_checks[r2])
        c_parity[r1] ^= c_parity[r2]
        c_has_bnd[r1] = c_has_bnd[r1] or c_has_bnd[r2]
        c_syn_mask[r1] ^= c_syn_mask[r2]
        if graphlike:
            c_mechs[r1].extend(c_mechs[r2])
        else:
            for mj in sorted(c_mechs[r2], key=lambda j: (prior_rank[j], j)):
                c_mechs[r1].append(mj)
                add_basis_col(r1, mj)
        del c_checks[r2], c_mechs[r2], c_parity[r2], c_has_bnd[r2], c_syn_mask[r2], c_basis[r2]
        return r1

    def is_valid(root: int) -> bool:
        if graphlike:
            return c_has_bnd[root] or (c_parity[root] == 0)
        rem = c_syn_mask[root]
        bmap = c_basis[root]
        for pb in sorted(bmap.keys(), reverse=True):
            if (rem >> pb) & 1:
                rem ^= bmap[pb][0]
        return rem == 0

    for dc in bit_positions(syndrome):
        ensure_cluster(dc)

    # Half-edge growth progress for graphlike UF (Delfosse & Nickerson, Quantum 5, 595, 2021, §2.2):
    # Each edge requires 2 half-edge increments to be absorbed; an edge with TWO invalid
    # cluster endpoints receives +2 in one step (completing before single-endpoint boundary edges!).
    edge_half_progress = [0] * problem.num_mechanisms

    while True:
        invalid = [r for r in sorted(c_checks.keys()) if not is_valid(r)]
        if not invalid:
            break
        progressed = False

        if (not use_bp_guidance) and graphlike:
            invalid_set = set(invalid)
            touched_edges: set[int] = set()
            for r in invalid:
                for c_in in c_checks[r]:
                    for mj in check_vars[c_in]:
                        if not used_mech[mj]:
                            touched_edges.add(mj)
            if not touched_edges:
                return DecodeResult(
                    decoder=dec_name,
                    error_mask=0,
                    converged=False,
                    observable_mask=0,
                    hamming_weight=0,
                    soft_weight=0,
                    stats={"bp_guided": False, "infeasible_syndrome": True},
                )
            newly_grown: list[int] = []
            for mj in touched_edges:
                dets = bit_positions(problem.mech_detectors[mj])
                # Count how many endpoints of `mj` currently lie in an invalid cluster
                inc_invalid = sum(
                    1 for d in dets if parent_c[d] != -1 and find_c(d) in invalid_set
                )
                edge_half_progress[mj] += max(1, inc_invalid)
                progressed = True
                if edge_half_progress[mj] >= 2:
                    newly_grown.append(mj)

            for mj in sorted(newly_grown, key=lambda j: (prior_rank[j], j)):
                used_mech[mj] = True
                dets = bit_positions(problem.mech_detectors[mj])
                if not dets:
                    continue
                r0 = ensure_cluster(dets[0])
                for d_other in dets[1:]:
                    r0 = merge_c(r0, ensure_cluster(d_other))
                c_mechs[r0].append(mj)
                if len(dets) == 1:
                    c_has_bnd[r0] = True
        else:
            for r in invalid:
                if r not in c_checks or is_valid(r):
                    continue
                best_j: int | None = None
                best_sc: tuple[int, int] | None = None
                for c_in in c_checks[r]:
                    for mj in check_vars[c_in]:
                        if not used_mech[mj]:
                            sc = (prior_rank[mj], mj)
                            if best_sc is None or sc < best_sc:
                                best_sc = sc
                                best_j = mj
                if best_j is None:
                    continue
                used_mech[best_j] = True
                progressed = True
                dets = bit_positions(problem.mech_detectors[best_j])
                r_cur = find_c(r)
                for dc in dets:
                    r_cur = merge_c(r_cur, ensure_cluster(dc))
                c_mechs[r_cur].append(best_j)
                if len(dets) == 1:
                    c_has_bnd[r_cur] = True
                if not graphlike:
                    add_basis_col(r_cur, best_j)

        if not progressed:
            return DecodeResult(
                decoder=dec_name,
                error_mask=0,
                converged=False,
                observable_mask=0,
                hamming_weight=0,
                soft_weight=0,
                stats={"bp_guided": use_bp_guidance, "infeasible_syndrome": True},
            )

    # Recover error inside each cluster (Peeling first if graphlike, else local GF(2))
    e_total = 0
    peeled_cnt = 0
    gf2_cnt = 0
    max_cc = 0

    for r in sorted(c_checks.keys()):
        if len(c_checks[r]) > max_cc:
            max_cc = len(c_checks[r])
        if c_syn_mask[r] == 0:
            continue
        if graphlike:
            e_c, peeled_ok = _peel_graphlike_cluster(
                c_checks[r], c_mechs[r], problem.mech_detectors, syndrome
            )
            if peeled_ok:
                e_total ^= e_c
                peeled_cnt += 1
                continue
        # Local GF(2) fallback (Delfosse-Londe-Beverland 2022)
        if not c_basis[r]:
            for mj in sorted(c_mechs[r], key=lambda j: (prior_rank[j], j)):
                add_basis_col(r, mj)
        rem = c_syn_mask[r]
        e_c = 0
        for pb in sorted(c_basis[r].keys(), reverse=True):
            if (rem >> pb) & 1:
                bv, bm = c_basis[r][pb]
                rem ^= bv
                e_c ^= bm
        if rem != 0:
            return DecodeResult(
                decoder=dec_name,
                error_mask=0,
                converged=False,
                observable_mask=0,
                hamming_weight=0,
                soft_weight=0,
                stats={"bp_guided": use_bp_guidance, "cluster_solve_failed": True},
            )
        e_total ^= e_c
        gf2_cnt += 1

    conv = problem.syndrome_of(e_total) == syndrome
    return DecodeResult(
        decoder=dec_name,
        error_mask=e_total,
        converged=conv,
        observable_mask=problem.observable_of(e_total),
        hamming_weight=e_total.bit_count(),
        soft_weight=problem.soft_weight_of(e_total),
        stats={
            "bp_guided": use_bp_guidance,
            "bp_iterations": bp_used,
            "bp_converged_directly": False,
            "num_clusters": len(c_checks),
            "peeled_clusters": peeled_cnt,
            "gf2_clusters": gf2_cnt,
            "max_cluster_checks": max_cc,
        },
    )
