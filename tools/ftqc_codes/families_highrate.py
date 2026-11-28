"""Paper family specializations for high-rate quantum LDPC codes.

Covers four families:
  1. HGP trivial-lift equivalence (`arXiv:0903.0566`, `arXiv:2012.04068`) and
     lift-sensitive finite-LP distance witnesses (`arXiv:2503.07567`).
  2. Two-Block Group-Algebra (`2BGA`) CSS codes (`arXiv:2306.16400`) over
     Abelian and non-Abelian finite groups, reconciled with `codes.bb`.
  3. Balanced-product orbit quotients (`arXiv:2012.09271`) and Quantum Tanner
     codes on Left-Right Cayley Complexes & lifted group actions
     (`arXiv:2202.13641`, `QuantumExpanders.jl` at `71058f87`), with ordered
     local neighborhoods and explicit `"checked"` vs `"assumed"` hypothesis tags.
  4. Pair-Partition CPM & APM ultra-high-rate CSS codes (`arXiv:2607.14091`,
     `arXiv:2609.30069`, `arXiv:2610.06749`), reproducing `[[184, 50]]`,
     `[[256, 68]]`, `[[576, 148]]`, and `[[1152, 580]]` with unreplayed large
     distances recorded under `paper_claims` and fail-closed non-CSS refusal.
"""
import itertools

from . import certify, codes, constructions, gf2, interchange


# ---------------------------------------------------------------- Row 1: HGP & finite LP
def trivial_lift_hgp_equivalence(h1, n1, h2, n2, *, name=None):
    """Build HGP via `constructions.hgp` and via `lifted_product` at `L=1` (`B = H_2^T`)."""
    c_hgp = constructions.hgp(h1, n1, h2, n2, name=name or "HGP")
    h1_rows = constructions._parse_seed_rows(h1, n1)
    h2_rows = constructions._parse_seed_rows(h2, n2)
    h2_t = codes._transpose(h2_rows, n2)
    a_mat = [[(r >> j) & 1 for j in range(n1)] for r in h1_rows]
    b_mat = [[(r >> j) & 1 for j in range(len(h2_rows))] for r in h2_t]
    c_lp1 = constructions.lifted_product(a_mat, b_mat, 1, name=name or "HGP")
    if c_hgp.sha256() != c_lp1.sha256():
        raise certify.Refusal("trivial lift L=1 did not match Tillich-Zemor HGP")
    return c_hgp, c_lp1


def raveendran_lp_instance(lift, *, name=None):
    """Finite-length 2x2 cyclic LP family exhibiting lift-dependent distance (`arXiv:2503.07567`).

    Uses asymmetric base polynomials `A = [[1+x, 1+x^2], [1+x^2, 1+x^3]]` modulo `x^lift - 1`
    whose cycle structure and minimum distance depend non-monotonically on `lift`.
    """
    if lift < 4:
        raise codes.SpecError("raveendran_lp_instance requires lift >= 4")
    a = [[[0, 1], [0, 2]], [[0, 2], [0, 3]]]
    b = [[[0, 1], [0, 2]], [[0, 2], [0, 3]]]
    code = constructions.lifted_product(
        a, b, lift, coeff_encoding="exponents", name=name or f"RaveendranLP(L={lift})"
    )
    code.meta["construction"]["provenance"] = {
        "paper": "Panteleev-Kalachev (2021) / Raveendran, Declercq & Vasic (2025)",
        "arxiv": ["2012.04068", "2503.07567"],
        "locator": f"2x2 cyclic LP at finite lift L={lift}",
    }
    return code


# ---------------------------------------------------------------- Row 2: Two-Block Group Algebra (2BGA)
def two_block_ga(
    group_size,
    left_perms,
    right_perms,
    a_elems,
    b_elems,
    *,
    shift_convention="row_to_lg_rb",
    name=None,
):
    """Two-Block Group Algebra (`2BGA`) CSS code `[A | B]`, `[B^T | A^T]` (`arXiv:2306.16400`).

    For any finite group `G` (Abelian or non-Abelian), left regular action `L_a`
    commutes with right regular action `R_b` (`L_a(R_b(g)) = a g b = R_b(L_a(g))`),
    so `H_X = [L(A) | R(B)]` and `H_Z = [R(B)^T | L(A)^T]` satisfy `H_X H_Z^T = 0`!
    """
    L, R, ident_e, inv, abelian = constructions.validate_group_action(
        group_size, left_perms, right_perms
    )
    _, _, a_mat, a_json, _ = constructions._normalize_group_matrix([[a_elems]], group_size, "a_elems")
    _, _, b_mat, b_json, _ = constructions._normalize_group_matrix([[b_elems]], group_size, "b_elems")
    supp_a = a_mat[0][0]
    supp_b = b_mat[0][0]

    use_inv = shift_convention == "row_to_inv"
    if shift_convention not in ("row_to_lg_rb", "row_to_inv"):
        raise codes.SpecError(f"unknown shift_convention {shift_convention!r}")

    am = [0] * group_size
    bm = [0] * group_size
    for row_g in range(group_size):
        ra = 0
        for g in supp_a:
            gg = inv[g] if use_inv else g
            ra ^= 1 << L[gg][row_g]
        am[row_g] = ra
        rb = 0
        for g in supp_b:
            gg = inv[g] if use_inv else g
            rb ^= 1 << R[gg][row_g]
        bm[row_g] = rb

    bt = codes._transpose(bm, group_size)
    at = codes._transpose(am, group_size)
    hx = [am[i] | (bm[i] << group_size) for i in range(group_size)]
    hz = [bt[i] | (at[i] << group_size) for i in range(group_size)]

    code = codes.CSSCode(
        name or f"2BGA(|G|={group_size})",
        2 * group_size,
        hx,
        hz,
        {
            "family": "two_block_ga",
            "two_block_ga": {
                "group_size": group_size,
                "abelian": abelian,
                "identity": ident_e,
                "a_elements": a_json[0][0],
                "b_elements": b_json[0][0],
                "shift_convention": shift_convention,
            },
            "provenance": {
                "paper": "Lin & Pryadko, Quantum two-block group algebra codes",
                "arxiv": "2306.16400",
            },
        },
    )
    if not code.commutes():
        raise certify.Refusal("2BGA checks do not commute")
    return code


def reconcile_bb_as_2bga(l, m, A, B, *, name=None):
    """Build `BB(l, m, A, B)` via `two_block_ga` over `Z_l x Z_m` and verify exact match with `codes.bb`."""
    l, m = int(l), int(m)
    bb_code = codes.bb(l, m, A, B, name=name)
    elems = [(i, j) for i in range(l) for j in range(m)]
    idx = {e: k for k, e in enumerate(elems)}
    perms = [
        [idx[((g[0] + h[0]) % l, (g[1] + h[1]) % m)] for h in elems]
        for g in elems
    ]

    def parse_monomials(terms):
        out = []
        for t in terms:
            t = t.strip().lower()
            k = int(t[1:])
            if t[0] == "x":
                out.append(idx[(k % l, 0)])
            else:
                out.append(idx[(0, k % m)])
        return out

    ga_code = two_block_ga(
        l * m,
        perms,
        perms,
        parse_monomials(A),
        parse_monomials(B),
        shift_convention="row_to_lg_rb",
        name=bb_code.name,
    )
    if ga_code.sha256() != bb_code.sha256():
        raise certify.Refusal("2BGA specialization did not match codes.bb under canonical Z_l x Z_m order")
    # Also verify explicit coordinate map from inverted shift convention `row_to_inv`!
    ga_inv = two_block_ga(
        l * m,
        perms,
        perms,
        parse_monomials(A),
        parse_monomials(B),
        shift_convention="row_to_inv",
        name=bb_code.name,
    )
    inv_g = [idx[((-e[0]) % l, (-e[1]) % m)] for e in elems]
    perm_q = inv_g + [l * m + x for x in inv_g]
    hx_m, hz_m, cmap = interchange.apply_coordinate_map(
        2 * l * m,
        ga_inv.hx,
        ga_inv.hz,
        {"perm_q": perm_q, "perm_hx": inv_g, "perm_hz": inv_g},
    )
    mapped = codes.CSSCode(bb_code.name, 2 * l * m, hx_m, hz_m)
    if mapped.sha256() != bb_code.sha256():
        raise certify.Refusal("2BGA inverted-shift coordinate map failed to match codes.bb")
    return ga_code, cmap


# ---------------------------------------------------------------- Row 3: Balanced Products & Quantum Tanner
def balanced_product(
    group_size,
    left_perms,
    right_perms,
    hx_a,
    na,
    right_act_a_q,
    right_act_a_chk,
    hz_b,
    nb,
    left_act_b_q,
    left_act_b_chk,
    *,
    name=None,
):
    """Breuckmann-Eberhardt balanced-product CSS quotient `C(A) (x)_G C(B)` (`arXiv:2012.09271`).

    `A` carries a validated right `G`-action (`(x . g1) . g2 = x . (g1 g2)`) commuting
    with `hx_a`, and `B` carries a validated left `G`-action (`g1 . (g2 . y) = (g1 g2) . y`)
    commuting with `hz_b`.  Quotients `C_i(A) x C_j(B)` by `(x . g, y) ~ (x, g . y)`
    (equivalently `(x . g, g^{-1} . y)`).
    """
    L, R, ident_e, inv, abelian = constructions.validate_group_action(
        group_size, left_perms, right_perms
    )
    ma = len(hx_a)
    mb = len(hz_b)
    ra_q = [interchange.validate_permutation(p, na, "right_act_a_q") for p in right_act_a_q]
    ra_c = [interchange.validate_permutation(p, ma, "right_act_a_chk") for p in right_act_a_chk]
    lb_q = [interchange.validate_permutation(p, nb, "left_act_b_q") for p in left_act_b_q]
    lb_c = [interchange.validate_permutation(p, mb, "left_act_b_chk") for p in left_act_b_chk]

    # Verify right-action law `(x . g1) . g2 == x . (g1 * g2)` where `g1 * g2 = L[g1][g2]`
    # and left-action law `g1 . (g2 . y) == (g1 * g2) . y`!
    for g1 in range(group_size):
        for g2 in range(group_size):
            prod = L[g1][g2]
            if any(ra_q[g2][ra_q[g1][x]] != ra_q[prod][x] for x in range(na)):
                raise certify.Refusal("right_act_a_q violates right group action law (x.g1).g2 == x.(g1*g2)")
            if any(ra_c[g2][ra_c[g1][c]] != ra_c[prod][c] for c in range(ma)):
                raise certify.Refusal("right_act_a_chk violates right group action law")
            if any(lb_q[g1][lb_q[g2][y]] != lb_q[prod][y] for y in range(nb)):
                raise certify.Refusal("left_act_b_q violates left group action law g1.(g2.y) == (g1*g2).y")
            if any(lb_c[g1][lb_c[g2][d]] != lb_c[prod][d] for d in range(mb)):
                raise certify.Refusal("left_act_b_chk violates left group action law")

    # Verify G-equivariance of classical differentials H_A and H_B
    for g in range(group_size):
        for c in range(ma):
            lhs = interchange.transport_word(hx_a[c], ra_q[g])
            rhs = hx_a[ra_c[g][c]]
            if lhs != rhs:
                raise certify.Refusal("H_A is not equivariant under right G-action")
        for d in range(mb):
            lhs = interchange.transport_word(hz_b[d], lb_q[g])
            rhs = hz_b[lb_c[g][d]]
            if lhs != rhs:
                raise certify.Refusal("H_B is not equivariant under left G-action")

    def build_orbits(size_x, act_r_x, size_y, act_l_y):
        pair_to_orb = {}
        orbits = []
        for x in range(size_x):
            for y in range(size_y):
                if (x, y) in pair_to_orb:
                    continue
                orb = sorted({(act_r_x[g][x], act_l_y[inv[g]][y]) for g in range(group_size)})
                if len(orb) != group_size:
                    raise certify.Refusal("balanced_product requires free G-orbits on product cells")
                oid = len(orbits)
                for p in orb:
                    pair_to_orb[p] = oid
                orbits.append(orb)
        return orbits, pair_to_orb

    q_vv_orbs, q_vv_map = build_orbits(na, ra_q, nb, lb_q)
    q_cc_orbs, q_cc_map = build_orbits(ma, ra_c, mb, lb_c)
    cx_orbs, _ = build_orbits(ma, ra_c, nb, lb_q)
    cz_orbs, _ = build_orbits(na, ra_q, mb, lb_c)

    n_vv = len(q_vv_orbs)
    n_total = n_vv + len(q_cc_orbs)
    hb_t = codes._transpose(hz_b, nb)
    ha_t = codes._transpose(hx_a, na)

    hx = []
    for orb in cx_orbs:
        c_rep, y_rep = orb[0]
        row = 0
        # Boundary on V_A x V_B: sum_{x in H_A[c_rep]} orbit(x, y_rep)
        for x in gf2.supp(hx_a[c_rep]):
            row ^= 1 << q_vv_map[(x, y_rep)]
        # Boundary on C_A x C_B: sum_{d in H_B^T[y_rep]} orbit(c_rep, d)
        for d in gf2.supp(hb_t[y_rep]):
            row ^= 1 << (n_vv + q_cc_map[(c_rep, d)])
        hx.append(row)

    hz = []
    for orb in cz_orbs:
        x_rep, d_rep = orb[0]
        row = 0
        for y in gf2.supp(hz_b[d_rep]):
            row ^= 1 << q_vv_map[(x_rep, y)]
        for c in gf2.supp(ha_t[x_rep]):
            row ^= 1 << (n_vv + q_cc_map[(c, d_rep)])
        hz.append(row)

    code = codes.CSSCode(
        name or f"BalancedProduct(|G|={group_size})",
        n_total,
        hx,
        hz,
        {
            "family": "balanced_product",
            "balanced_product": {
                "group_size": group_size,
                "abelian": abelian,
                "vv_orbits": len(q_vv_orbs),
                "cc_orbits": len(q_cc_orbs),
                "x_check_orbits": len(cx_orbs),
                "z_check_orbits": len(cz_orbs),
            },
            "provenance": {
                "paper": "Breuckmann & Eberhardt, Balanced Product Quantum Codes",
                "arxiv": "2012.09271",
            },
        },
    )
    if not code.commutes():
        raise certify.Refusal("balanced product checks do not commute")
    return code


def _check_classical_dual_pair(H, G, width, label):
    if not H or not G or any(len(r) != width for r in H + G):
        raise codes.SpecError(f"{label}: invalid local matrix dimensions for width {width}")
    h_rows = gf2.rows_from_lists(H)
    g_rows = gf2.rows_from_lists(G)
    if any(gf2.dot(hr, gr) for hr in h_rows for gr in g_rows):
        raise certify.Refusal(f"{label}: local H * G^T != 0 over GF(2)")
    if gf2.rank(h_rows) + gf2.rank(g_rows) != width:
        raise certify.Refusal(f"{label}: rank(H) + rank(G) != {width}")
    return h_rows, g_rows


def quantum_tanner_lrcc(
    group_size,
    left_perms,
    right_perms,
    gen_a,
    gen_b,
    h_a,
    g_a,
    h_b,
    g_b,
    *,
    name=None,
):
    """Quantum Tanner code on the Left-Right Cayley Complex `Gamma(G, A, B)` (`arXiv:2202.13641`).

    Matches `QuantumExpanders.jl` (`QuantumTannerCode` at `71058f87`):
      - Checks symmetric generating sets `A = A^{-1}`, `B = B^{-1}`, `e not in A, B`,
        Total Non-Conjugacy (`a g != g b`), joint generation `<A, B> = G`, and local
        duality `H_A G_A^T = 0`, `H_B G_B^T = 0`.
      - Places qubits on unordered squares `{(g,0), (ag,1), (gb,1), (agb,0)}` in exact
        `QuantumExpanders.jl` enumeration order (`n = |G| |A| |B| // 2`).
      - Records ordered `|A| x |B|` local neighborhoods at `V_0` (Z-checks) and `V_1`
        (X-checks) plus individual `"checked"` / `"assumed"` hypothesis tags.
    """
    L, R, ident_e, inv, abelian = constructions.validate_group_action(
        group_size, left_perms, right_perms
    )
    da, db = len(gen_a), len(gen_b)
    if da < 1 or db < 1 or len(set(gen_a)) != da or len(set(gen_b)) != db:
        raise codes.SpecError("QuantumTannerCode requires distinct generator lists A and B")
    if ident_e in gen_a or ident_e in gen_b:
        raise certify.Refusal("LRCC generating sets A and B must not contain the identity")
    if {inv[a] for a in gen_a} != set(gen_a):
        raise certify.Refusal("generating set A is not closed under inverses (A != A^{-1})")
    if {inv[b] for b in gen_b} != set(gen_b):
        raise certify.Refusal("generating set B is not closed under inverses (B != B^{-1})")

    # Total Non-Conjugacy (TNC): for all a in A, b in B, g in G, a*g != g*b
    for a in gen_a:
        for b in gen_b:
            for g in range(group_size):
                if L[a][g] == R[b][g]:
                    raise certify.Refusal(
                        f"Total Non-Conjugacy violated: a={a}, b={b}, g={g} satisfy a*g == g*b"
                    )

    # Check <A union B> == G by BFS from identity
    closure = {ident_e}
    queue = [ident_e]
    gens_union = list(gen_a) + list(gen_b)
    while queue:
        cur = queue.pop()
        for s in gens_union:
            nxt = L[s][cur]
            if nxt not in closure:
                closure.add(nxt)
                queue.append(nxt)
    generates_g = len(closure) == group_size
    if not generates_g:
        raise certify.Refusal("A union B does not generate the full group G")

    ha_rows, ga_rows = _check_classical_dual_pair(h_a, g_a, da, "C_A")
    hb_rows, gb_rows = _check_classical_dual_pair(h_b, g_b, db, "C_B")

    a_pos = {a: i for i, a in enumerate(gen_a)}
    b_pos = {b: j for j, b in enumerate(gen_b)}

    # Exact QuantumExpanders.jl square enumeration
    squares = []
    square_key_to_id = {}
    v0_views = [[None] * (da * db) for _ in range(group_size)]
    v1_views = [[None] * (da * db) for _ in range(group_size)]

    for g in range(group_size):
        for ai, a in enumerate(gen_a):
            for bi, b in enumerate(gen_b):
                ag = L[a][g]
                gb = R[b][g]
                agb = R[b][ag]
                ai_inv = a_pos[inv[a]]
                bi_inv = b_pos[inv[b]]

                # Type 0 square: V0={g, agb}, V1={ag, gb}
                key0 = (frozenset((g, agb)), frozenset((ag, gb)))
                if key0 not in square_key_to_id:
                    qid = len(squares)
                    square_key_to_id[key0] = qid
                    squares.append({"v0": sorted((g, agb)), "v1": sorted((ag, gb))})
                    v0_views[g][ai * db + bi] = qid
                    v0_views[agb][ai_inv * db + bi_inv] = qid
                    v1_views[gb][ai * db + bi_inv] = qid
                    v1_views[ag][ai_inv * db + bi] = qid

                # Type 1 square: V1={g, agb}, V0={ag, gb}
                key1 = (frozenset((ag, gb)), frozenset((g, agb)))
                if key1 not in square_key_to_id:
                    qid = len(squares)
                    square_key_to_id[key1] = qid
                    squares.append({"v0": sorted((ag, gb)), "v1": sorted((g, agb))})
                    v1_views[g][ai * db + bi] = qid
                    v1_views[agb][ai_inv * db + bi_inv] = qid
                    v0_views[gb][ai * db + bi_inv] = qid
                    v0_views[ag][ai_inv * db + bi] = qid

    n = len(squares)
    if n != (group_size * da * db) // 2:
        raise certify.Refusal(f"unexpected square count {n} != {group_size * da * db // 2}")
    if any(q is None for view in v0_views + v1_views for q in view):
        raise certify.Refusal("incomplete local square bijection phi_v")

    # Tensor local bases: beta_0 = G_A (x) G_B for Z-checks on V_0; beta_1 = H_A (x) H_B for X-checks on V_1
    beta0 = codes._kron(ga_rows, da, gb_rows, db)
    beta1 = codes._kron(ha_rows, da, hb_rows, db)

    hz = []
    v0_check_ranges = []
    for v in range(group_size):
        start = len(hz)
        view = v0_views[v]
        for b_mask in beta0:
            hz.append(gf2.from_supp(view[loc] for loc in gf2.supp(b_mask)))
        v0_check_ranges.append([start, len(hz)])

    hx = []
    v1_check_ranges = []
    for v in range(group_size):
        start = len(hx)
        view = v1_views[v]
        for b_mask in beta1:
            hx.append(gf2.from_supp(view[loc] for loc in gf2.supp(b_mask)))
        v1_check_ranges.append([start, len(hx)])

    hypotheses = {
        "group_axioms": "checked",
        "symmetric_A": "checked",
        "symmetric_B": "checked",
        "total_non_conjugacy": "checked",
        "jointly_generates_group": "checked",
        "local_classical_duality": "checked",
        "spectral_expansion_ramanujan": "assumed",
        "robustness_coboundary_expansion": "assumed",
    }

    code = codes.CSSCode(
        name or f"QuantumTanner(|G|={group_size},{da}x{db})",
        n,
        hx,
        hz,
        {
            "family": "quantum_tanner",
            "tanner": {
                "variant": "lrcc_squares",
                "group_size": group_size,
                "abelian": abelian,
                "local_shape": [da, db],
                "v0_neighborhoods": v0_views,
                "v1_neighborhoods": v1_views,
                "v0_check_ranges": v0_check_ranges,
                "v1_check_ranges": v1_check_ranges,
                "hypotheses": hypotheses,
            },
            "provenance": {
                "paper": "Leverrier & Zemor, Quantum Tanner codes (FOCS 2022)",
                "arxiv": "2202.13641",
                "reference_impl": "QuantumExpanders.jl@71058f87 QuantumTannerCode",
            },
        },
    )
    if not code.commutes():
        raise certify.Refusal("Quantum Tanner checks do not commute")
    return code


def lifted_quantum_tanner(
    group_size,
    left_perms,
    right_perms,
    multiset_a,
    multiset_b,
    h0,
    g0,
    h1,
    g1,
    h0p,
    g0p,
    h1p,
    g1p,
    *,
    name=None,
):
    """Lifted left-right Quantum Tanner code (`QuantumExpanders.jl` `QuantumTannerViaLeftRightActions`).

    Implements Eq. (2) of Leverrier & Zemor (`lifted_quantum_tanner.jl` at `71058f87`):
      `H_X = [H_0 (x) G_0' (x) I_{|G|} ; (H_1 (x) G_1' (x) I_{|G|}) L_A R_B]`,
      `H_Z = [(G_0 (x) H_1' (x) I_{|G|}) R_B ; (G_1 (x) H_0' (x) I_{|G|}) L_A]`,
    with `n = n_A * n_B * |G|` physical qubits ordered as `((i * n_B + j) * |G| + g)`.
    """
    L, R, ident_e, inv, abelian = constructions.validate_group_action(
        group_size, left_perms, right_perms
    )
    na, nb = len(multiset_a), len(multiset_b)
    h0_r, g0_r = _check_classical_dual_pair(h0, g0, na, "C_0(A)")
    h1_r, g1_r = _check_classical_dual_pair(h1, g1, na, "C_1(A)")
    h0p_r, g0p_r = _check_classical_dual_pair(h0p, g0p, nb, "C_0'(B)")
    h1p_r, g1p_r = _check_classical_dual_pair(h1p, g1p, nb, "C_1'(B)")

    n = na * nb * group_size

    def q_id(i_a, j_b, g_el):
        return (i_a * nb + j_b) * group_size + g_el

    # In QuantumExpanders.jl: _lambda(a) has entry (g, a*g) = 1; _rho(b) has entry (g, g*b^{-1}) = 1!
    def expand_block(mat_a, mat_b, use_la, use_rb):
        kron_ab = codes._kron(mat_a, na, mat_b, nb)
        rows_out = []
        # Order by group element g so each vertex g has a contiguous check slice & local neighborhood!
        for r_ab in kron_ab:
            supp_cells = gf2.supp(r_ab)
            for g in range(group_size):
                mask = 0
                for cell in supp_cells:
                    ia, jb = divmod(cell, nb)
                    g_col = g
                    if use_la:
                        g_col = L[multiset_a[ia]][g_col]
                    if use_rb:
                        g_col = R[inv[multiset_b[jb]]][g_col]
                    mask ^= 1 << q_id(ia, jb, g_col)
                rows_out.append(mask)
        return rows_out

    hx = expand_block(h0_r, g0p_r, False, False) + expand_block(h1_r, g1p_r, True, True)
    hz = expand_block(g0_r, h1p_r, False, True) + expand_block(g1_r, h0p_r, True, False)

    # Local views at each g in G for Z-checks (R_B slice and L_A slice) and X-checks (identity slice and L_A R_B slice)
    v0_views = []
    v1_views = []
    for g in range(group_size):
        v0_views.append(
            [q_id(ia, jb, R[inv[multiset_b[jb]]][g]) for ia in range(na) for jb in range(nb)]
        )
        v0_views.append(
            [q_id(ia, jb, L[multiset_a[ia]][g]) for ia in range(na) for jb in range(nb)]
        )
        v1_views.append([q_id(ia, jb, g) for ia in range(na) for jb in range(nb)])
        v1_views.append(
            [q_id(ia, jb, R[inv[multiset_b[jb]]][L[multiset_a[ia]][g]]) for ia in range(na) for jb in range(nb)]
        )

    code = codes.CSSCode(
        name or f"LiftedQT(|G|={group_size},{na}x{nb})",
        n,
        hx,
        hz,
        {
            "family": "lifted_quantum_tanner",
            "tanner": {
                "variant": "lifted_left_right_actions",
                "group_size": group_size,
                "abelian": abelian,
                "local_shape": [na, nb],
                "v0_neighborhoods": v0_views,
                "v1_neighborhoods": v1_views,
                "hypotheses": {
                    "group_axioms": "checked",
                    "commuting_left_right_translations": "checked",
                    "local_classical_duality": "checked",
                    "spectral_expansion": "assumed",
                    "robustness": "assumed",
                },
            },
            "provenance": {
                "paper": "Leverrier & Zemor, Small quantum Tanner codes / QuantumExpanders.jl",
                "reference_impl": "QuantumExpanders.jl@71058f87 QuantumTannerViaLeftRightActions",
            },
        },
    )
    if not code.commutes():
        raise certify.Refusal("lifted Quantum Tanner checks do not commute")
    return code


def quantum_expanders_d8_36_8_3():
    """Documented `QuantumExpanders.jl` `[[36, 8, 3]]` nonabelian `D_8` Quantum Tanner code."""
    # D_8 = <s, r | s^2 = r^4 = (sr)^2 = 1>, elements (s_exp in 0..1, r_exp in 0..3)
    # Multiplication: (s1, r1) * (s2, r2): since r * s = s * r^3, if s2 == 1: (s1 ^ 1, (-r1 + r2) % 4) else (s1, (r1 + r2) % 4)
    elems = [(se, re) for se in range(2) for re in range(4)]
    idx = {e: i for i, e in enumerate(elems)}

    def mul(x, y):
        s1, r1 = x
        s2, r2 = y
        if s2 == 0:
            return (s1, (r1 + r2) % 4)
        return (s1 ^ 1, (-r1 + r2) % 4)

    left_p = [[idx[mul(g, h)] for h in elems] for g in elems]
    right_p = [[idx[mul(h, g)] for h in elems] for g in elems]
    s_el = idx[(1, 0)]
    r1_el = idx[(0, 1)]
    r2_el = idx[(0, 2)]
    r3_el = idx[(0, 3)]
    sr_el = idx[mul((1, 0), (0, 1))]
    sr3_el = idx[mul((1, 0), (0, 3))]

    gen_a = [s_el, r1_el, r3_el]
    gen_b = [sr_el, sr3_el, r2_el]
    h_a = [[1, 0, 1], [0, 1, 1]]
    g_a = [[1, 1, 1]]
    h_b = [[1, 1, 1]]
    g_b = [[1, 1, 0], [1, 0, 1]]
    return quantum_tanner_lrcc(
        8, left_p, right_p, gen_a, gen_b, h_a, g_a, h_b, g_b, name="QE_D8_36_8_3"
    )


def quantum_expanders_s3_288_8():
    """Documented `QuantumExpanders.jl` `[[288, 8, (<=15, <=15)]]` nonabelian `S_3` Quantum Tanner code."""
    # S_3 permutations on (0, 1, 2) (1-based [1,2,3] in Julia)
    s3_elems = list(itertools.permutations(range(3)))
    s3_idx = {p: i for i, p in enumerate(s3_elems)}

    def comp(p, q):
        # Julia PermGroup multiplies left-to-right or right-to-left; regular representation satisfies L_p(q)=p*q, R_p(q)=q*p
        return tuple(q[p[i]] for i in range(3))

    left_p = [[s3_idx[comp(g, h)] for h in s3_elems] for g in s3_elems]
    right_p = [[s3_idx[comp(h, g)] for h in s3_elems] for g in s3_elems]

    e = s3_idx[(0, 1, 2)]
    c23 = s3_idx[(0, 2, 1)]
    c123 = s3_idx[(1, 2, 0)]
    c12 = s3_idx[(1, 0, 2)]
    c132 = s3_idx[(2, 0, 1)]

    multiset_a = [e, e, c23, c23, c123, c12, c12, c132]
    multiset_b = [e, e, c23, c123, c12, c132]

    h844 = [
        [1, 0, 0, 0, 0, 1, 1, 1],
        [0, 1, 0, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 0, 1],
        [0, 0, 0, 1, 1, 1, 1, 0],
    ]
    g844 = [
        [0, 1, 1, 1, 1, 0, 0, 0],
        [1, 0, 1, 1, 0, 1, 0, 0],
        [1, 1, 0, 1, 0, 0, 1, 0],
        [1, 1, 1, 0, 0, 0, 0, 1],
    ]
    h633 = [
        [1, 0, 0, 0, 1, 1],
        [0, 1, 0, 1, 0, 1],
        [0, 0, 1, 1, 1, 0],
    ]
    g633 = [
        [0, 1, 1, 1, 0, 0],
        [1, 0, 1, 0, 1, 0],
        [1, 1, 0, 0, 0, 1],
    ]
    # p1 = 1:8 (identity), p2 = [1, 2, 4, 3, 6, 5] (0-based [0, 1, 3, 2, 5, 4])
    p2 = [0, 1, 3, 2, 5, 4]
    h1p = [[r[j] for j in p2] for r in h633]
    g1p = [[r[j] for j in p2] for r in g633]

    code = lifted_quantum_tanner(
        6, left_p, right_p, multiset_a, multiset_b,
        h844, g844, h844, g844,
        h633, g633, h1p, g1p,
        name="QE_S3_288_8",
    )
    code.meta["paper_claims"] = {"d_upper": 15, "dx_upper": 15, "dz_upper": 15}
    return code


# ---------------------------------------------------------------- Row 4: Pair-Partition CPM & APM CSS Codes
def _has_even_multiplicity_partition(diffs):
    counts = {}
    for d in diffs:
        counts[d] = counts.get(d, 0) + 1
    return all(v % 2 == 0 for v in counts.values())


def cpm_pair_partition(
    lift,
    ex_shifts,
    ez_shifts,
    *,
    code_type="css",
    verify_pair_partition=True,
    name=None,
    provenance=None,
    paper_claims=None,
    extra_meta=None,
):
    """Expand a pair-partition CPM seed (`E_X`, `E_Z` mod `P`) into a binary `CSSCode`.

    Convention (`arXiv:2610.06749` Appendix F.1): entry `(r, c)` of `C_P^s` is `1`
    precisely when `r - c == s (mod P)`, and binary qubit index is `b * P + c`.
    """
    if str(code_type).lower() != "css":
        raise certify.Refusal(
            "non-CSS / halved pair-partition code requires a separate symplectic-code extension; "
            "cpm_pair_partition only constructs CSS codes"
        )
    if not isinstance(lift, int) or isinstance(lift, bool) or lift < 1:
        raise codes.SpecError("lift P must be a positive integer")
    if not ex_shifts or not ez_shifts or not ex_shifts[0] or not ez_shifts[0]:
        raise codes.SpecError("ex_shifts and ez_shifts must be non-empty 2D shift arrays")
    jx, L = len(ex_shifts), len(ex_shifts[0])
    jz = len(ez_shifts)
    if L % 2 != 0 and verify_pair_partition:
        raise certify.Refusal(f"pair-partition CPM requires even row degree L, got L={L}")
    if any(len(r) != L for r in ex_shifts) or any(len(r) != L for r in ez_shifts):
        raise codes.SpecError("ragged CPM shift array")

    ex_mod = [[int(s) % lift for s in row] for row in ex_shifts]
    ez_mod = [[int(s) % lift for s in row] for row in ez_shifts]

    if verify_pair_partition:
        for i in range(jx):
            for j in range(jz):
                # Overlap between block row i of H_X and block row j of H_Z at relative row shift delta:
                # exists c with r_x - c = ex[i][b] and r_z - c = ez[j][b] <=> r_z - r_x = ez[j][b] - ex[i][b] (mod P)!
                diffs = [(ez_mod[j][b] - ex_mod[i][b]) % lift for b in range(L)]
                if not _has_even_multiplicity_partition(diffs):
                    raise certify.Refusal(
                        f"pair-partition condition violated at block pair ({i}, {j}): "
                        f"modular shift differences {diffs} do not pair up"
                    )

    n = L * lift
    hx = []
    for i in range(jx):
        for r in range(lift):
            rx = 0
            for b in range(L):
                c = (r - ex_mod[i][b]) % lift
                rx |= 1 << (b * lift + c)
            hx.append(rx)

    hz = []
    for j in range(jz):
        for r in range(lift):
            rz = 0
            for b in range(L):
                c = (r - ez_mod[j][b]) % lift
                rz |= 1 << (b * lift + c)
            hz.append(rz)

    meta = {
        "family": "cpm_pair_partition",
        "cpm_pp": {
            "lift": lift,
            "J_x": jx,
            "J_z": jz,
            "L": L,
            "ex_shifts": ex_mod,
            "ez_shifts": ez_mod,
            "pair_partition_verified": bool(verify_pair_partition),
        },
        "provenance": dict(provenance or {}),
        "paper_claims": dict(paper_claims or {}),
    }
    if extra_meta:
        meta.update(extra_meta)

    code = codes.CSSCode(name or f"CPM_PP(J={jx},L={L},P={lift})", n, hx, hz, meta)
    if not code.commutes():
        raise certify.Refusal("expanded CPM matrices fail H_X * H_Z^T == 0")
    return code


def pp184_50():
    """Lee et al. (`arXiv:2609.30069v1` Eq. D10-D11) `[[184, 50, 10]]` CPM-PP CSS code (`P=23`)."""
    ex = [[0] * 8, [0, 12, 8, 21, 6, 1, 19, 15], [0, 9, 18, 11, 7, 17, 10, 4]]
    ez = [[0, 15, 7, 22, 7, 0, 22, 15], [0, 1, 2, 4, 0, 1, 2, 4], [0, 5, 17, 6, 6, 17, 5, 0]]
    return cpm_pair_partition(
        23,
        ex,
        ez,
        name="PP184",
        provenance={
            "paper": "Lee, Okada, Maskara, Kasai & Zhou, Design Principles for Ultra-High-Rate Quantum Codes",
            "version": "arXiv:2609.30069v1",
            "locator": "Appendix D, Eq. (D10)-(D11)",
        },
        paper_claims={"d": 10, "d_x": 10, "d_z": 10},
    )


def pp256_68():
    """Maskara et al. (`arXiv:2610.06749v1` Eq. F1, Table I & V) `[[256, 68, 12]]` CPM-PP CSS code (`P=32`)."""
    ex = [
        [9, 0, 15, 11, 16, 14, 29, 6],
        [10, 18, 23, 14, 6, 18, 26, 23],
        [23, 11, 3, 12, 16, 22, 3, 10],
    ]
    ez = [
        [23, 30, 3, 6, 30, 9, 17, 4],
        [11, 25, 17, 3, 7, 7, 20, 30],
        [31, 29, 12, 19, 24, 29, 12, 28],
    ]
    return cpm_pair_partition(
        32,
        ex,
        ez,
        name="PP256",
        provenance={
            "paper": "Maskara et al., Efficient Logic with Ultra-High-Rate Quantum Codes",
            "version": "arXiv:2610.06749v1",
            "locator": "Table I, Table V, and Eq. (F1)",
        },
        paper_claims={"d": 12, "d_x": 12, "d_z": 12},
        extra_meta={
            "pivot_partition": {
                "sx_blocks": [0, 1, 2],
                "sz_blocks": [3, 4, 5],
                "info_blocks": [6, 7],
                "pair_program": [[0, 1, 13]],
                "spatial_layers": 2,
                "spatial_f_exponents": [0, 1, 3],
                "spatial_g_exponents": [0],
            }
        },
    )


def pp576_148():
    """Maskara et al. (`arXiv:2610.06749v1` Eq. F2, Table I & V) `[[576, 148, 18]]` CPM-PP CSS code (`P=72`)."""
    ex = [
        [50, 69, 70, 40, 5, 39, 5, 0],
        [32, 62, 69, 51, 61, 18, 6, 25],
        [33, 37, 30, 47, 46, 63, 23, 23],
    ]
    ez = [
        [32, 22, 29, 49, 59, 48, 36, 25],
        [40, 67, 60, 70, 69, 37, 69, 30],
        [53, 17, 18, 29, 66, 45, 11, 3],
    ]
    return cpm_pair_partition(
        72,
        ex,
        ez,
        name="PP576",
        provenance={
            "paper": "Maskara et al., Efficient Logic with Ultra-High-Rate Quantum Codes",
            "version": "arXiv:2610.06749v1",
            "locator": "Table I, Table V, and Eq. (F2)",
        },
        paper_claims={"d": 18, "d_x": 18, "d_z": 18},
    )


def apm1152_580(*, verify_affine_crt=True):
    """Maskara et al. (`arXiv:2610.06749v1` Eq. F3-F7, Table I & V) `[[1152, 580, 12]]` rate-1/2 APM code.

    Also cross-checks the native affine permutation formula `c = a r + b (mod 96)` (Eq. F3)
    against the Chinese-remainder `3 x 32` circulant decomposition (Eq. F4-F5)!
    """
    def perm_block_crt(kind3, exp3, exp32):
        cols = [0] * 96
        for ar in range(3):
            ac = (2 - ar) % 3 if kind3 == "J" else (ar - exp3) % 3
            for jr in range(32):
                jc = (jr - exp32) % 32
                cols[32 * ar + jr] = 32 * ac + jc
        return cols

    f_specs = [("I", 0, 3), ("I", 0, 8), ("J", 0, 10), ("C", 2, 18), ("I", 0, 5), ("I", 0, 26)]
    g_specs = [("J", 0, 9), ("C", 1, 21), ("I", 0, 6), ("I", 0, 0), ("I", 0, 24), ("I", 0, 9)]
    F = [perm_block_crt(*s) for s in f_specs]
    G = [perm_block_crt(*s) for s in g_specs]

    if verify_affine_crt:
        # Verify Eq. (F3)-(F4) CRT coordinate map on all 12 blocks!
        y_seq = [0]
        for _ in range(31):
            y_seq.append((13 * y_seq[-1] + 3) % 32)
        if len(set(y_seq)) != 32:
            raise AssertionError("CRT y_j sequence is not a bijection on Z_32")
        crt_to_native = [0] * 96
        native_to_crt = [0] * 96
        for z in range(96):
            a_c = z % 3
            j_c = y_seq.index(z % 32)
            idx_crt = 32 * a_c + j_c
            crt_to_native[idx_crt] = z
            native_to_crt[z] = idx_crt

        ab_params = [(61, 15), (1, 24), (89, 62), (25, 22), (85, 93), (25, 78)]
        cd_params = [(5, 41), (85, 77), (73, 66), (1, 0), (1, 72), (37, 9)]
        for idx_i, (a_p, b_p) in enumerate(ab_params):
            reord = [native_to_crt[(a_p * crt_to_native[rc] + b_p) % 96] for rc in range(96)]
            if reord != F[idx_i]:
                raise certify.Refusal(f"APM1152 block F_{idx_i} CRT decomposition disagrees with Eq. (F3)")
        for idx_i, (c_p, d_p) in enumerate(cd_params):
            reord = [native_to_crt[(c_p * crt_to_native[rc] + d_p) % 96] for rc in range(96)]
            if reord != G[idx_i]:
                raise certify.Refusal(f"APM1152 block G_{idx_i} CRT decomposition disagrees with Eq. (F3)")

    def inv_perm(p):
        out = [0] * len(p)
        for r, c in enumerate(p):
            out[c] = r
        return out

    FT = [inv_perm(p) for p in F]
    GT = [inv_perm(p) for p in G]

    hx = []
    hz = []
    for rb in range(3):
        blocks_x = [F[(j - rb) % 6] for j in range(6)] + [G[(j - rb) % 6] for j in range(6)]
        blocks_z = [GT[(rb - j) % 6] for j in range(6)] + [FT[(rb - j) % 6] for j in range(6)]
        for r in range(96):
            rx = 0
            rz = 0
            for cb in range(12):
                rx |= 1 << (cb * 96 + blocks_x[cb][r])
                rz |= 1 << (cb * 96 + blocks_z[cb][r])
            hx.append(rx)
            hz.append(rz)

    code = codes.CSSCode(
        "APM1152",
        1152,
        hx,
        hz,
        {
            "family": "apm_pair_partition",
            "apm": {
                "cyclic_lift": 32,
                "block_shape_cyclic": [9, 36],
                "affine_crt_verified": bool(verify_affine_crt),
            },
            "provenance": {
                "paper": "Maskara et al., Efficient Logic with Ultra-High-Rate Quantum Codes",
                "version": "arXiv:2610.06749v1",
                "locator": "Table I, Table V, and Eq. (F3)-(F7)",
            },
            "paper_claims": {"d": 12, "d_x": 12, "d_z": 12},
        },
    )
    if not code.commutes():
        raise certify.Refusal("APM1152 matrices do not commute")
    return code
