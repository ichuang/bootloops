"""Leg G: symmetry-reduced minimum-weight exclusion via verified code automorphisms.

Adapts Okada & Kasai (`arXiv:2607.14091` Algorithm 2 & Theorem 2) with exact
progressive orbit-coverage certificates that hold for arbitrary (prime or
composite, Abelian or non-Abelian) verified automorphism groups:

  1. Verify every supplied generator `g_i` preserves `rowspace(H_X)` and
     `rowspace(H_Z)`.
  2. Partition `0..n-1` into exact orbits `O_0, ..., O_{m-1}` under `<g_i>`,
     recording an explicit generator-word witness mapping each `q in O_t` to
     its canonical root `r_t = min(O_t)`.
  3. Because a minimum-weight nontrivial logical operator `E*` (`|E*| = d`) has
     no nonempty proper subset in `ker(H)` (any stabilizer subset `S_1 subset E*`
     would yield a strictly lighter nontrivial logical `E* \\ S_1`, and any
     nontrivial logical subset `S_1 subset E*` would already have `|S_1| < d`),
     branch-and-bound syndrome expansion from `S = {r_t}` with progressive
     forbidden mask `F = O_0 union ... union O_{t-1}` is 100% sound for excluding
     all nontrivial logicals of weight `<= W`!
"""
from . import certify, codes, gf2, interchange, symmetry


def _compose_perm(p2, p1):
    """Apply `p1` then `p2`: `out[j] = p2[p1[j]]`."""
    return [p2[p1[j]] for j in range(len(p1))]


def verify_automorphism_orbits(code, aut_generators):
    """Verify automorphism generators and compute replayable orbit-root witnesses."""
    if not isinstance(aut_generators, (list, tuple)) or not aut_generators:
        raise codes.SpecError("aut_generators must be a non-empty list of permutations")
    gens = []
    for g in aut_generators:
        symmetry.verify_code_automorphism(code, g)
        gens.append([int(x) for x in g])

    # Also include inverses of generators so undirected orbit traversal maps any q -> root
    inv_gens = []
    for g in gens:
        inv = [0] * code.n
        for j, gj in enumerate(g):
            inv[gj] = j
        inv_gens.append(inv)

    visited = [False] * code.n
    orbits = []
    for start in range(code.n):
        if visited[start]:
            continue
        # BFS to find full orbit, then BFS backwards from root = min(orbit) so we have
        # an explicit forward word in `gens` sending any x in orbit -> root!
        comp = {start}
        q = [start]
        visited[start] = True
        while q:
            u = q.pop()
            for g in gens + inv_gens:
                v = g[u]
                if not visited[v]:
                    visited[v] = True
                    comp.add(v)
                    q.append(v)
        members = sorted(comp)
        root = members[0]
        # In a finite permutation group, `inv_gens[k]` equals a power of `gens[k]`, and
        # backward step along `gens[k]` (`g[x] == cur` <=> `x == inv_g[cur]`) means applying
        # `gens[k]` to `x` moves one step closer to `root`!
        word_to_root = {root: []}
        bfs = [root]
        head = 0
        while head < len(bfs):
            cur = bfs[head]
            head += 1
            for g_idx, inv_g in enumerate(inv_gens):
                prev = inv_g[cur]
                if prev not in word_to_root:
                    # Applying gens[g_idx] to `prev` produces `cur`, followed by word_to_root[cur] -> `root`!
                    word_to_root[prev] = [g_idx] + word_to_root[cur]
                    bfs.append(prev)
        if len(word_to_root) != len(members):
            raise certify.Refusal("automorphism orbit failed strong connectivity check")
        orbits.append(
            {
                "root": root,
                "members": members,
                "words_to_root": {str(m): word_to_root[m] for m in members},
            }
        )

    return {"generators": gens, "orbits": orbits}


def recheck_orbit_coverage(code, generators, orbits):
    """Cheap deterministic verifier for claimed automorphism generators and orbit partition."""
    if not generators or not orbits:
        raise certify.Refusal("Leg G certificate lacks generators or orbits")
    gens = []
    for g in generators:
        symmetry.verify_code_automorphism(code, g)
        gens.append([int(x) for x in g])

    covered = set()
    for orb in orbits:
        root = orb.get("root")
        members = orb.get("members")
        words = orb.get("words_to_root") or {}
        if not isinstance(members, list) or not members or root not in members:
            raise certify.Refusal("invalid orbit entry in Leg G certificate")
        for m in members:
            if m in covered or not (0 <= m < code.n):
                raise certify.Refusal("orbit partition has duplicate or out-of-range qubit")
            covered.add(m)
            seq = words.get(str(m))
            if not isinstance(seq, list):
                raise certify.Refusal(f"missing word_to_root witness for qubit {m}")
            cur = m
            for g_idx in seq:
                if not (0 <= g_idx < len(gens)):
                    raise certify.Refusal("generator index out of range in word_to_root")
                cur = gens[g_idx][cur]
            if cur != root:
                raise certify.Refusal(f"word_to_root for qubit {m} landed on {cur} != root {root}")
    if len(covered) != code.n:
        raise certify.Refusal(f"orbits cover {len(covered)} qubits, expected n={code.n}")
    return True


def orbit_min_weight(
    code,
    quantity,
    aut_generators,
    *,
    target=None,
    max_nodes=2_000_000,
):
    """Symmetry-reduced branch-and-bound exclusion for `d_x` or `d_z` (`arXiv:2607.14091` Alg. 2).

    Searches iteratively by weight `w = 1, 2, ...` up to `target - 1` (or `n` when
    `target is None`) so stopping on `max_nodes` still returns the highest
    completed lower bound `w`.
    """
    if quantity not in ("d_x", "d_z"):
        raise codes.SpecError("orbit_min_weight supports 'd_x' and 'd_z'")
    if not code.commutes():
        raise certify.Refusal("code checks do not commute")

    orb_cert = verify_automorphism_orbits(code, aut_generators)
    orbits = orb_cert["orbits"]
    n = code.n
    checks = code.hz if quantity == "d_x" else code.hx
    same_stabs = code.hx if quantity == "d_x" else code.hz
    stab_ech = gf2.Echelon(same_stabs)
    k = n - gf2.rank(code.hx) - gf2.rank(code.hz)
    if k == 0:
        return {
            "empty": True,
            "lower": None,
            "upper": None,
            "witness": None,
            "exact": True,
            "nodes": 0,
            "orbits": orbits,
            "generators": orb_cert["generators"],
        }

    # Deduplicate identical check rows for faster syndrome tracking
    uniq_checks = list(dict.fromkeys(r for r in checks if r != 0))
    m_chk = len(uniq_checks)
    col_syndromes = [0] * n
    for c_i, r in enumerate(uniq_checks):
        bit = 1 << c_i
        for q in gf2.supp(r):
            col_syndromes[q] |= bit
    max_col_wt = max((s.bit_count() for s in col_syndromes), default=1)
    if max_col_wt == 0:
        max_col_wt = 1

    max_w_goal = (target - 1) if target is not None else n
    nodes = [0]
    closed_stabs = [0]
    budget_hit = [None]
    found_witness = [None]

    def search_at_cap(w_cap):
        forb_global = 0
        for orb in orbits:
            r0 = orb["root"]
            s0 = 1 << r0
            f0 = forb_global | s0
            syn0 = col_syndromes[r0]

            def dfs(s_mask, wt_s, f_mask, syn_mask):
                if found_witness[0] is not None or budget_hit[0] is not None:
                    return
                nodes[0] += 1
                if nodes[0] > max_nodes:
                    budget_hit[0] = f"stopped at weight {w_cap}: exceeded max_nodes={max_nodes}"
                    return
                if syn_mask == 0:
                    if s_mask not in stab_ech:
                        found_witness[0] = s_mask
                    else:
                        closed_stabs[0] += 1
                    return
                rem = w_cap - wt_s
                if rem <= 0:
                    return
                syn_wt = syn_mask.bit_count()
                if (syn_wt + max_col_wt - 1) // max_col_wt > rem:
                    return

                # Inspect unsatisfied checks: pick smallest |A_c| and compute greedy disjoint packing b2
                best_c_avail = None
                best_deg = n + 1
                unsat_avails = []
                cur_syn = syn_mask
                while cur_syn:
                    lsb = cur_syn & -cur_syn
                    c_idx = lsb.bit_length() - 1
                    cur_syn ^= lsb
                    av = uniq_checks[c_idx] & ~f_mask
                    if av == 0:
                        return
                    deg = av.bit_count()
                    if deg < best_deg:
                        best_deg = deg
                        best_c_avail = av
                    unsat_avails.append((deg, av))

                unsat_avails.sort(key=lambda t: t[0])
                disjoint_cnt = 0
                used_vars = 0
                for _, av in unsat_avails:
                    if (used_vars & av) == 0:
                        used_vars |= av
                        disjoint_cnt += 1
                        if disjoint_cnt > rem:
                            return

                cur_f = f_mask
                for v in gf2.supp(best_c_avail):
                    cur_f |= 1 << v
                    dfs(s_mask | (1 << v), wt_s + 1, cur_f, syn_mask ^ col_syndromes[v])
                    if found_witness[0] is not None or budget_hit[0] is not None:
                        return

            dfs(s0, 1, f0, syn0)
            if found_witness[0] is not None or budget_hit[0] is not None:
                return
            for q_orb in orb["members"]:
                forb_global |= 1 << q_orb

    proven_lower = 1
    for w_try in range(1, max_w_goal + 1):
        search_at_cap(w_try)
        if budget_hit[0] is not None:
            break
        if found_witness[0] is not None:
            w_val = gf2.wt(found_witness[0])
            return {
                "empty": False,
                "lower": w_val,
                "upper": w_val,
                "witness": found_witness[0],
                "exact": True,
                "nodes": nodes[0],
                "closed_stabilizer_partials": closed_stabs[0],
                "total_roots": len(orbits),
                "orbits": orbits,
                "generators": orb_cert["generators"],
                "budget": None,
            }
        proven_lower = w_try + 1

    return {
        "empty": False,
        "lower": proven_lower,
        "upper": None,
        "witness": None,
        "exact": False,
        "nodes": nodes[0],
        "closed_stabilizer_partials": closed_stabs[0],
        "total_roots": len(orbits),
        "orbits": orbits,
        "generators": orb_cert["generators"],
        "budget": budget_hit[0],
    }


def infer_default_automorphisms(code):
    """Infer cyclic shift automorphisms from `code.meta` when available."""
    cmeta = code.meta.get("construction") or {}
    if cmeta.get("kind") == "lifted_product":
        lift = cmeta["seed_inputs"]["lift"]
        coord = cmeta["generator_convention"]["coordinate_order"]
        if lift >= 2:
            return [symmetry.cyclic_cell_shift_perm(code.n, lift, 1, coordinate_order=coord)]
    if "cpm_pp" in code.meta:
        lift = code.meta["cpm_pp"]["lift"]
        if lift >= 2:
            return [symmetry.cyclic_cell_shift_perm(code.n, lift, 1, coordinate_order="cell_major")]
    if code.meta.get("family") == "bb":
        l, m = code.meta["l"], code.meta["m"]
        N = l * m
        gens = []
        if l >= 2:
            px = [((i // m + 1) % l) * m + (i % m) for i in range(N)]
            gens.append(px + [N + x for x in px])
        if m >= 2:
            py = [(i // m) * m + ((i % m + 1) % m) for i in range(N)]
            gens.append(py + [N + x for x in py])
        if gens:
            return gens
    return None


def run_leg(receipt, code, aut_generators=None, targets=None, max_nodes=2_000_000):
    """Run Leg G on `d_x` and `d_z` and attach certified bounds to `receipt`."""
    targets = targets or {}
    gens = aut_generators if aut_generators is not None else infer_default_automorphisms(code)
    for q in ("d_x", "d_z"):
        ev = {
            "leg": "G",
            "quantity": q,
            "method": "symmetry-reduced syndrome branch-and-bound (Okada-Kasai orbit coverage)",
        }
        if gens is None:
            ev.update(status="skipped", reason="no code automorphisms supplied or inferred")
            receipt["evidence"].append(ev)
            continue
        res = orbit_min_weight(code, q, gens, target=targets.get(q), max_nodes=max_nodes)
        if res.get("empty"):
            ev.update(status="empty", reason="k = 0")
            receipt["evidence"].append(ev)
            continue
        if res["witness"] is not None:
            certify.add_witness(receipt, code, q, res["witness"])
        ev.update(
            status="proved",
            claim=f"{q} >= {res['lower']}",
            nodes=res["nodes"],
            closed_stabilizer_partials=res["closed_stabilizer_partials"],
            total_roots=res["total_roots"],
            generators=res["generators"],
            orbits=res["orbits"],
        )
        if res.get("budget"):
            ev["budget"] = res["budget"]
        certify.add_lower_bound(receipt, q, res["lower"], ev)
    return receipt


def recheck(code, quantity, ev, *, full=False, max_nodes=2_000_000):
    """Re-verify a Leg G evidence entry (cheap mode checks automorphisms + orbits; full re-runs search)."""
    claimed = int(ev["claim"].split(">=")[1])
    gens = ev.get("generators")
    orbits = ev.get("orbits")
    recheck_orbit_coverage(code, gens, orbits)
    if full:
        res = orbit_min_weight(code, quantity, gens, target=claimed, max_nodes=max_nodes)
        if res["lower"] < claimed:
            raise certify.Refusal(
                f"{quantity}: leg G re-run proves only {res['lower']}, receipt claims {claimed}"
            )
    return claimed
