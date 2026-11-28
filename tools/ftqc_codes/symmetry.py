"""Structured logical basis producers and symmetry / lift chain-map verifiers.

Implements:
  1. Odd-order Abelian lift spectroscopy (Lee, arXiv:2607.05386): Frobenius
     character packets, Kunneth packet bases, reciprocal trace-dual pairing,
     and binary cross-verification; checks odd-order premise fail-closed.
  2. Canonical cyclic lifted-product symmetry & seed-orbit coverage
     (Zheng et al., arXiv:2607.28605): physical automorphism verification,
     induced symplectic logical action, and exact translated seed-orbit rank.
  3. Group-extension lift chain maps (Hirasaki & Lee, arXiv:2607.28621):
     chain/cochain boundary commutation and induced logical transformations.
"""
import itertools
from math import gcd

from . import basis, bch, certify, codes, constructions, gf2, interchange


def _lcm(a, b):
    return (a // gcd(a, b)) * b


# ---------------------------------------------------------------- 1. GF(2^m) linear algebra & Lee spectroscopy
class _GF2m:
    """Exact arithmetic over GF(2^m) = GF(2)[t]/(poly)."""

    def __init__(self, m, poly=None):
        self.m = m
        self.poly = poly if poly is not None else bch.find_irreducible(m)
        self.order = (1 << m) - 1

    def add(self, a, b):
        return a ^ b

    def mul(self, a, b):
        if a == 0 or b == 0:
            return 0
        return bch._pmulmod(a, b, self.poly)

    def pow(self, a, e):
        if e == 0:
            return 1
        if a == 0:
            return 0
        return bch._ppowmod(a, e % self.order, self.poly)

    def inv(self, a):
        if a == 0:
            raise ZeroDivisionError("0 has no inverse in GF(2^m)")
        if self.m == 1:
            return 1
        return bch._ppowmod(a, self.order - 1, self.poly)

    def primitive_element(self):
        if self.m == 1:
            return 1
        factors = bch._prime_factors(self.order)
        for cand in range(2, 1 << self.m):
            if all(self.pow(cand, self.order // q) != 1 for q in factors):
                return cand
        raise AssertionError("no primitive element found in GF(2^m)")


def _rref_gf2m(mat, ncols, fld):
    """Left-to-right Gauss-Jordan RREF over `_GF2m`. Returns `(rref_rows, pivots)`."""
    rows = [list(r) for r in mat]
    nrows = len(rows)
    pivots = []
    r_idx = 0
    for c in range(ncols):
        piv = next((i for i in range(r_idx, nrows) if rows[i][c] != 0), None)
        if piv is None:
            continue
        if piv != r_idx:
            rows[r_idx], rows[piv] = rows[piv], rows[r_idx]
        inv_p = fld.inv(rows[r_idx][c])
        rows[r_idx] = [fld.mul(x, inv_p) for x in rows[r_idx]]
        for i in range(nrows):
            if i != r_idx and rows[i][c] != 0:
                factor = rows[i][c]
                rows[i] = [rows[i][j] ^ fld.mul(factor, rows[r_idx][j]) for j in range(ncols)]
        pivots.append(c)
        r_idx += 1
        if r_idx == nrows:
            break
    return [rows[i] for i in range(len(pivots))], pivots


def _kernel_gf2m(mat, ncols, fld):
    rref, pivots = _rref_gf2m(mat, ncols, fld)
    piv_set = set(pivots)
    out = []
    for f in range(ncols):
        if f in piv_set:
            continue
        v = [0] * ncols
        v[f] = 1
        for row, p in zip(rref, pivots):
            v[p] = row[f]
        out.append(v)
    return out


def _coker_complement_gf2m(mat, nrows, ncols, fld):
    """Direct-sum complement basis of `im(mat)` inside `K^nrows`."""
    mt = [[mat[i][j] for i in range(nrows)] for j in range(ncols)]
    _, col_pivots = _rref_gf2m(mt, nrows, fld)
    piv_set = set(col_pivots)
    out = []
    for i in range(nrows):
        if i not in piv_set:
            e = [0] * nrows
            e[i] = 1
            out.append(e)
    return out


def _invert_gf2m(sq, fld):
    d = len(sq)
    if d == 0:
        return []
    aug = [list(sq[i]) + [1 if i == j else 0 for j in range(d)] for i in range(d)]
    rref, pivots = _rref_gf2m(aug, 2 * d, fld)
    if pivots != list(range(d)):
        raise certify.Refusal("singular packet pairing matrix over extension field")
    return [r[d:] for r in rref]


def _dualize_vectors_gf2m(primal_vecs, dual_vecs, fld):
    """Transform `dual_vecs` so `primal_vecs[a] . new_dual[b] = delta_{a,b}` over `fld`."""
    d = len(primal_vecs)
    if len(dual_vecs) != d:
        raise certify.Refusal("primal/dual Kunneth dimensions disagree")
    if d == 0:
        return []
    gram = [
        [
            sum_xor([fld.mul(primal_vecs[a][q], dual_vecs[b][q]) for q in range(len(primal_vecs[a]))])
            for b in range(d)
        ]
        for a in range(d)
    ]
    ginv = _invert_gf2m(gram, fld)
    n_dim = len(dual_vecs[0])
    out = []
    for b in range(d):
        vec = [0] * n_dim
        for c in range(d):
            coeff = ginv[c][b]
            if coeff:
                for q in range(n_dim):
                    vec[q] ^= fld.mul(coeff, dual_vecs[c][q])
        out.append(vec)
    return out


def sum_xor(items):
    acc = 0
    for x in items:
        acc ^= x
    return acc


def _invert_gf2_matrix(mat):
    d = len(mat)
    aug = [
        sum((mat[i][j] & 1) << j for j in range(d)) | (1 << (d + i))
        for i in range(d)
    ]
    for c in range(d):
        piv = next((r for r in range(c, d) if (aug[r] >> c) & 1), None)
        if piv is None:
            raise certify.Refusal("degenerate trace Gram matrix")
        if piv != c:
            aug[c], aug[piv] = aug[piv], aug[c]
        for r in range(d):
            if r != c and ((aug[r] >> c) & 1):
                aug[r] ^= aug[c]
    return [[(aug[i] >> (d + j)) & 1 for j in range(d)] for i in range(d)]


def spectral_abelian_basis(
    orders,
    a,
    b,
    *,
    coeff_encoding="exponents",
    coordinate_order="cell_major",
    slot_prefix="s",
    name=None,
):
    """Lee's logical spectroscopy for lifted-product codes over odd-order Abelian groups (`arXiv:2607.05386`).

    `orders`: integer lift `L` (for cyclic `Z_L`) or list of cyclic factor orders `[l_1, ..., l_d]`
              for `G = Z_{l_1} x ... x Z_{l_d}`.
    Checks odd-order premise `|G| % 2 == 1` fail-closed.
    """
    if isinstance(orders, int) and not isinstance(orders, bool):
        ord_list = [orders]
        is_cyclic_1d = True
    elif isinstance(orders, (list, tuple)) and orders and all(
        isinstance(x, int) and not isinstance(x, bool) for x in orders
    ):
        ord_list = [int(x) for x in orders]
        is_cyclic_1d = False
    else:
        raise codes.SpecError("orders must be a positive integer or non-empty list of integers")

    if any(x < 1 for x in ord_list):
        raise codes.SpecError("all Abelian factor orders must be >= 1")

    g_size = 1
    exp_lcm = 1
    for x in ord_list:
        g_size *= x
        exp_lcm = _lcm(exp_lcm, x)

    if g_size % 2 == 0:
        raise certify.Refusal(
            f"spectral basis requires odd-order Abelian group (got |G|={g_size}, even); "
            "GF(2)[G] is non-semisimple in characteristic 2 — use derive_paired_basis instead"
        )

    elems = list(itertools.product(*(range(x) for x in ord_list)))
    elem_idx = {t: i for i, t in enumerate(elems)}

    if is_cyclic_1d:
        L = ord_list[0]
        code = constructions.lifted_product(
            a, b, L, coeff_encoding=coeff_encoding, coordinate_order=coordinate_order, name=name
        )
        m, na, a_supp, _, _ = constructions._normalize_cyclic_matrix(a, L, coeff_encoding, "a")
        r, kb, b_supp, _, _ = constructions._normalize_cyclic_matrix(b, L, coeff_encoding, "b")
    else:
        perms = [
            [
                elem_idx[tuple((g[j] + h[j]) % ord_list[j] for j in range(len(ord_list)))]
                for h in elems
            ]
            for g in elems
        ]
        # Convert tuple-exponent or index entries in a, b
        def conv_mat(mat, lbl):
            out = []
            for row in mat:
                orow = []
                for entry in row:
                    ids = []
                    for item in entry:
                        if isinstance(item, (list, tuple)):
                            t = tuple(int(v) % ord_list[j] for j, v in enumerate(item))
                            ids.append(elem_idx[t])
                        else:
                            ids.append(int(item))
                    orow.append(ids)
                out.append(orow)
            return out

        a_idx = conv_mat(a, "a")
        b_idx = conv_mat(b, "b")
        code = constructions.group_lifted_product(
            g_size, perms, perms, a_idx, b_idx, coordinate_order=coordinate_order, name=name
        )
        m, na, a_supp, _, _ = constructions._normalize_group_matrix(a_idx, g_size, "a")
        r, kb, b_supp, _, _ = constructions._normalize_group_matrix(b_idx, g_size, "b")

    if exp_lcm == 1:
        fld = _GF2m(1, 0b11)
        alpha = 1
    else:
        poly, alpha = bch.root_of_unity(exp_lcm)
        fld = _GF2m(bch._deg(poly), poly)
    prim = fld.primitive_element()

    # Enumerate Frobenius orbits (packets) on dual characters u in Z_{l_1} x ... x Z_{l_d}
    visited = set()
    packets = []
    for u in elems:
        if u in visited:
            continue
        orbit = []
        cur = u
        while cur not in orbit:
            orbit.append(cur)
            visited.add(cur)
            cur = tuple((2 * cur[j]) % ord_list[j] for j in range(len(ord_list)))
        packets.append(orbit)

    def chi_val(u, g_id):
        g = elems[g_id]
        exp_sum = 0
        for j, lj in enumerate(ord_list):
            exp_sum = (exp_sum + (exp_lcm // lj) * u[j] * g[j]) % exp_lcm
        return fld.pow(alpha, exp_sum)

    def neg_elem(g_id):
        g = elems[g_id]
        return elem_idx[tuple((-g[j]) % ord_list[j] for j in range(len(ord_list)))]

    total_cells = na * r + m * kb

    def q_idx(cell, g_id):
        return cell * g_size + g_id if coordinate_order == "cell_major" else g_id * total_cells + cell

    slots = []
    packet_summaries = []

    for pkt_id, orbit in enumerate(packets):
        u = orbit[0]
        m_deg = len(orbit)
        neg_u = tuple((-u[j]) % ord_list[j] for j in range(len(ord_list)))

        # Primitive central idempotent e_Omega in GF(2)[G]
        idem_supp = []
        for g_id in range(g_size):
            ng = neg_elem(g_id)
            tr = sum_xor(chi_val(v, ng) for v in orbit)
            if tr not in (0, 1):
                raise AssertionError("Frobenius packet character sum landed outside GF(2)")
            if tr == 1:
                idem_supp.append(g_id)

        # Evaluate A_u (m x na) and B_u (r x kb) over K_Omega subseteq GF(2^d)
        A_u = [
            [sum_xor(chi_val(u, g_id) for g_id in a_supp[i][u_c]) for u_c in range(na)]
            for i in range(m)
        ]
        B_u = [
            [sum_xor(chi_val(u, g_id) for g_id in b_supp[j][v_c]) for v_c in range(kb)]
            for j in range(r)
        ]
        A_ut = [[A_u[i][u_c] for i in range(m)] for u_c in range(na)]
        B_ut = [[B_u[j][v_c] for j in range(r)] for v_c in range(kb)]

        ker_A = _kernel_gf2m(A_u, na, fld)
        coker_A = _coker_complement_gf2m(A_u, m, na, fld)
        ker_B = _kernel_gf2m(B_u, kb, fld)
        coker_B = _coker_complement_gf2m(B_u, r, kb, fld)

        coker_At = _dualize_vectors_gf2m(ker_A, _coker_complement_gf2m(A_ut, na, m, fld), fld)
        ker_At = _dualize_vectors_gf2m(coker_A, _kernel_gf2m(A_ut, m, fld), fld)
        coker_Bt = _dualize_vectors_gf2m(ker_B, _coker_complement_gf2m(B_ut, kb, r, fld), fld)
        ker_Bt = _dualize_vectors_gf2m(coker_B, _kernel_gf2m(B_ut, r, fld), fld)

        # Subfield K_Omega basis {1, gamma, ..., gamma^{m_deg-1}} and trace-dual basis
        if m_deg == 1:
            theta = [1]
            theta_dual = [1]
        else:
            step = fld.order // ((1 << m_deg) - 1)
            gamma = fld.pow(prim, step)
            theta = [fld.pow(gamma, s) for s in range(m_deg)]

            def trace_k(val):
                acc = 0
                cur = val
                for _ in range(m_deg):
                    acc ^= cur
                    cur = fld.mul(cur, cur)
                if acc not in (0, 1):
                    raise AssertionError("subfield trace landed outside GF(2)")
                return acc

            t_gram = [[trace_k(fld.mul(theta[a_i], theta[b_i])) for b_i in range(m_deg)] for a_i in range(m_deg)]
            t_inv = _invert_gf2_matrix(t_gram)
            theta_dual = [
                sum_xor(theta[c_i] for c_i in range(m_deg) if t_inv[c_i][b_i])
                for b_i in range(m_deg)
            ]

        def trace_ko(val):
            acc = 0
            cur = val
            for _ in range(m_deg):
                acc ^= cur
                cur = fld.mul(cur, cur)
            return acc

        def lift_z(cell_vec, mult_scalar):
            word = 0
            for cell, coeff in enumerate(cell_vec):
                if coeff == 0:
                    continue
                sc = fld.mul(mult_scalar, coeff)
                for g_id in range(g_size):
                    bit = trace_ko(fld.mul(sc, chi_val(u, neg_elem(g_id))))
                    if bit:
                        word |= 1 << q_idx(cell, g_id)
            return word

        def lift_x(cell_vec, mult_scalar):
            word = 0
            for cell, coeff in enumerate(cell_vec):
                if coeff == 0:
                    continue
                sc = fld.mul(mult_scalar, coeff)
                for g_id in range(g_size):
                    bit = trace_ko(fld.mul(sc, chi_val(u, g_id)))
                    if bit:
                        word |= 1 << q_idx(cell, g_id)
            return word

        k_pkt_ext = len(ker_A) * len(coker_B) + len(coker_A) * len(ker_B)
        pkt_slots = 0

        # Sector I (left cells: na * r)
        for a_i in range(len(ker_A)):
            for b_i in range(len(coker_B)):
                z_cells = [0] * total_cells
                x_cells = [0] * total_cells
                for u_c in range(na):
                    for j_c in range(r):
                        z_cells[u_c * r + j_c] = fld.mul(ker_A[a_i][u_c], coker_B[b_i][j_c])
                        x_cells[u_c * r + j_c] = fld.mul(coker_At[a_i][u_c], ker_Bt[b_i][j_c])
                for f_i in range(m_deg):
                    zw = lift_z(z_cells, theta[f_i])
                    xw = lift_x(x_cells, theta_dual[f_i])
                    nm = f"{slot_prefix}_p{pkt_id}_L{a_i}_{b_i}_f{f_i}"
                    slots.append(
                        {
                            "name": nm,
                            "x": basis.pauli_from_xz(code.n, xw, 0, 0),
                            "z": basis.pauli_from_xz(code.n, 0, zw, 0),
                        }
                    )
                    pkt_slots += 1

        # Sector II (right cells: m * kb)
        for a_i in range(len(coker_A)):
            for b_i in range(len(ker_B)):
                z_cells = [0] * total_cells
                x_cells = [0] * total_cells
                for i_c in range(m):
                    for v_c in range(kb):
                        cell = na * r + i_c * kb + v_c
                        z_cells[cell] = fld.mul(coker_A[a_i][i_c], ker_B[b_i][v_c])
                        x_cells[cell] = fld.mul(ker_At[a_i][i_c], coker_Bt[b_i][v_c])
                for f_i in range(m_deg):
                    zw = lift_z(z_cells, theta[f_i])
                    xw = lift_x(x_cells, theta_dual[f_i])
                    nm = f"{slot_prefix}_p{pkt_id}_R{a_i}_{b_i}_f{f_i}"
                    slots.append(
                        {
                            "name": nm,
                            "x": basis.pauli_from_xz(code.n, xw, 0, 0),
                            "z": basis.pauli_from_xz(code.n, 0, zw, 0),
                        }
                    )
                    pkt_slots += 1

        packet_summaries.append(
            {
                "packet_id": pkt_id,
                "representative": list(u),
                "reciprocal_representative": list(neg_u),
                "self_reciprocal": neg_u in orbit,
                "degree": m_deg,
                "idempotent_support": idem_supp,
                "kunneth_dim_over_packet_field": k_pkt_ext,
                "binary_logical_slots": pkt_slots,
            }
        )

    # Cross-check with independent binary elimination!
    bin_basis = basis.derive_paired_basis(code)
    if len(slots) != bin_basis["k"]:
        raise certify.Refusal(
            f"spectral packet count {len(slots)} disagrees with binary elimination k={bin_basis['k']}"
        )

    spec_basis = basis.build_basis_receipt(
        code,
        slots,
        source="lee_abelian_spectroscopy",
        coordinate_order=coordinate_order,
        provenance={
            "paper": "Lee, Logical Spectroscopy of Abelian Lifted-Product Quantum Codes",
            "arxiv": "2607.05386",
            "orders": ord_list,
            "binary_crosscheck_k": bin_basis["k"],
        },
    )
    report = {
        "orders": ord_list,
        "group_size": g_size,
        "odd_order_checked": True,
        "packets": packet_summaries,
        "binary_crosscheck_k": bin_basis["k"],
        "basis_sha256": spec_basis["basis_sha256"],
    }
    return code, spec_basis, report


# ---------------------------------------------------------------- 2. Canonical cyclic LP symmetry & seed orbits
def cyclic_cell_shift_perm(n, lift, amount=1, coordinate_order="cell_major"):
    """Physical qubit permutation shifting each cyclic fiber `t -> (t + amount) % lift`."""
    if lift < 1 or n % lift != 0:
        raise codes.SpecError(f"lift {lift} does not divide n={n}")
    cells = n // lift
    perm = [0] * n
    for c in range(cells):
        for t in range(lift):
            nt = (t + amount) % lift
            if coordinate_order == "cell_major":
                perm[c * lift + t] = c * lift + nt
            elif coordinate_order == "fiber_major":
                perm[t * cells + c] = nt * cells + c
            else:
                raise codes.SpecError(f"unknown coordinate_order {coordinate_order!r}")
    return perm


def verify_code_automorphism(code, perm_q):
    """Verify `perm_q` is a qubit permutation preserving `rowspace(H_X)` and `rowspace(H_Z)`."""
    perm = interchange.validate_permutation(perm_q, code.n, "automorphism")
    ex = gf2.Echelon(code.hx)
    ez = gf2.Echelon(code.hz)
    tx = [interchange.transport_word(r, perm) for r in code.hx]
    tz = [interchange.transport_word(r, perm) for r in code.hz]
    if any(r not in ex for r in tx) or gf2.rank(tx) != ex.rank:
        raise certify.Refusal("permutation does not preserve X stabilizer space")
    if any(r not in ez for r in tz) or gf2.rank(tz) != ez.rank:
        raise certify.Refusal("permutation does not preserve Z stabilizer space")
    return {"n": code.n, "perm_q": perm, "preserves_hx": True, "preserves_hz": True}


def _transport_signed_pauli(p_tuple, perm_q):
    x, z, s = p_tuple
    return interchange.transport_word(x, perm_q), interchange.transport_word(z, perm_q), s


def induced_logical_action(code, basis_receipt, perm_q):
    """Compute and verify the symplectic logical transformation induced by `perm_q` on `basis_receipt`."""
    verify_code_automorphism(code, perm_q)
    basis.verify_basis_receipt(code, basis_receipt)
    perm = list(perm_q)
    n = code.n
    k = basis_receipt["k"]
    slots = basis_receipt["slots"]
    parsed_x = [basis.parse_pauli(s["x"], n) for s in slots]
    parsed_z = [basis.parse_pauli(s["z"], n) for s in slots]
    checks = basis._code_signed_checks(code)
    sgrp = basis._SignedStabilizerGroup(n, checks)

    def decompose(op):
        # For op commuting with all checks:
        # X-component at slot j is anticommutation with Z_j; Z-component at slot j is anticommutation with X_j.
        xm = 0
        zm = 0
        for j in range(k):
            if not basis.commute_pauli(op, parsed_z[j]):
                xm |= 1 << j
            if not basis.commute_pauli(op, parsed_x[j]):
                zm |= 1 << j
        emb0 = basis.parse_pauli(basis.embed_logical_pauli(basis_receipt, xm, zm, sign=0), n)
        dx, dz, ph = basis.mul_pauli(emb0, op)
        if ph % 2 != 0:
            raise certify.Refusal("automorphism image failed symplectic commutation")
        rx, rz, rs = sgrp.reduce((dx, dz, ph // 2))
        if rx or rz:
            raise certify.Refusal("automorphism image did not reduce into stabilizer group")
        return {"x_mask": xm, "z_mask": zm, "sign": rs}

    lx_action = [decompose(_transport_signed_pauli(parsed_x[i], perm)) for i in range(k)]
    lz_action = [decompose(_transport_signed_pauli(parsed_z[i], perm)) for i in range(k)]

    # Verify symplectic preservation directly
    for i in range(k):
        for j in range(k):
            comm_xz = (
                gf2.dot(lx_action[i]["x_mask"], lz_action[j]["z_mask"])
                ^ gf2.dot(lx_action[i]["z_mask"], lz_action[j]["x_mask"])
            ) == 0
            if comm_xz != (i != j):
                raise certify.Refusal("induced logical action is not symplectic")

    return {
        "code_sha256": code.sha256(),
        "basis_sha256": basis_receipt["basis_sha256"],
        "logical_x": lx_action,
        "logical_z": lz_action,
    }


def verify_seed_orbit_coverage(
    code,
    basis_receipt,
    seeds,
    lift,
    *,
    sector="Z",
    coordinate_order="cell_major",
):
    """Verify translated cyclic orbits of `seeds` and compute exact covered logical rank & slots (`arXiv:2607.28605`)."""
    basis.verify_basis_receipt(code, basis_receipt)
    if sector not in ("X", "Z"):
        raise codes.SpecError("sector must be 'X' or 'Z'")
    shift1 = cyclic_cell_shift_perm(code.n, lift, 1, coordinate_order=coordinate_order)
    verify_code_automorphism(code, shift1)

    k = basis_receipt["k"]
    slots = basis_receipt["slots"]
    opp_checks = code.hx if sector == "Z" else code.hz
    dual_ops = [
        basis.parse_pauli(s["x" if sector == "Z" else "z"], code.n)
        for s in slots
    ]

    seed_words = []
    for s in seeds:
        if isinstance(s, int):
            w = int(s)
        elif isinstance(s, list):
            w = gf2.from_supp(s)
        elif isinstance(s, str):
            px, pz, _ = basis.parse_pauli(s, code.n)
            w = pz if sector == "Z" else px
        else:
            raise codes.SpecError("seed must be bitmask int, support list, or signed Pauli string")
        if w <= 0 or (w >> code.n):
            raise certify.Refusal("seed support out of bounds")
        if any(gf2.dot(w, chk) for chk in opp_checks):
            raise certify.Refusal(f"seed anticommutes with opposite {('X' if sector == 'Z' else 'Z')} check")
        seed_words.append(w)

    ech = gf2.Echelon()
    orbit_records = []
    for s_idx, sw in enumerate(seed_words):
        translates = []
        for t in range(lift):
            perm_t = cyclic_cell_shift_perm(code.n, lift, t, coordinate_order=coordinate_order)
            tw = interchange.transport_word(sw, perm_t)
            if any(gf2.dot(tw, chk) for chk in opp_checks):
                raise certify.Refusal("translated seed failed check commutation")
            op_tuple = (0, tw, 0) if sector == "Z" else (tw, 0, 0)
            coord_mask = 0
            for j in range(k):
                if not basis.commute_pauli(op_tuple, dual_ops[j]):
                    coord_mask |= 1 << j
            ech.add(coord_mask)
            translates.append({"shift": t, "weight": gf2.wt(tw), "logical_mask": coord_mask})
        orbit_records.append({"seed_index": s_idx, "seed_weight": gf2.wt(sw), "translates": translates})

    # Determine which individual unit slots e_j lie inside the spanned subspace
    covered_slots = [
        slots[j]["name"] for j in range(k) if (1 << j) in ech
    ]
    return {
        "code_sha256": code.sha256(),
        "basis_sha256": basis_receipt["basis_sha256"],
        "sector": sector,
        "lift": lift,
        "k": k,
        "covered_rank": ech.rank,
        "full_coverage": ech.rank == k,
        "covered_unit_slots": covered_slots,
        "logical_span_basis_masks": list(ech.basis),
        "orbits": orbit_records,
    }


# ---------------------------------------------------------------- 3. Group-extension lift chain maps
def _matmul_gf2(A, n_mid, B, n_cols):
    """Multiply binary matrices `A` (`m x n_mid` row masks) and `B` (`n_mid x n_cols` row masks)."""
    out = []
    for r in A:
        acc = 0
        for j in gf2.supp(r):
            acc ^= B[j]
        out.append(acc)
    return out


def verify_chain_map(
    source_code,
    target_code,
    f1_rows,
    f0_rows,
    f2_rows,
    *,
    source_basis=None,
    target_basis=None,
):
    """Verify CSS chain map `(f_2, f_1, f_0)` from `source_code` to `target_code` (`arXiv:2607.28621`).

    Conventions:
      `f1_rows`: `target.n` row bitmasks over `source.n` columns (`C_1^src -> C_1^tgt`),
      `f0_rows`: `len(target.hx)` row bitmasks over `len(source.hx)` columns (`C_0^src -> C_0^tgt`),
      `f2_rows`: `len(target.hz)` row bitmasks over `len(source.hz)` columns (`C_2^src -> C_2^tgt`).
    Verifies:
      1. `H_X^tgt * f_1 == f_0 * H_X^src  (mod 2)`
      2. `f_1 * (H_Z^src)^T == (H_Z^tgt)^T * f_2  (mod 2)`
    and computes induced logical maps `M_Z: Z_src -> Z_tgt` and `M_X: X_tgt -> X_src`
    checking `M_X == M_Z^T`.
    """
    if not source_code.commutes() or not target_code.commutes():
        raise certify.Refusal("source and target codes must both commute")
    if (
        len(f1_rows) != target_code.n
        or len(f0_rows) != len(target_code.hx)
        or len(f2_rows) != len(target_code.hz)
    ):
        raise codes.SpecError("chain map row counts do not match target code dimensions")
    if any(r < 0 or (r >> source_code.n) for r in f1_rows):
        raise codes.SpecError("f1_rows entry out of range for source.n")
    if any(r < 0 or (r >> len(source_code.hx)) for r in f0_rows):
        raise codes.SpecError("f0_rows entry out of range for len(source.hx)")
    if any(r < 0 or (r >> len(source_code.hz)) for r in f2_rows):
        raise codes.SpecError("f2_rows entry out of range for len(source.hz)")

    lhs_1 = _matmul_gf2(target_code.hx, target_code.n, f1_rows, source_code.n)
    rhs_1 = _matmul_gf2(f0_rows, len(source_code.hx), source_code.hx, source_code.n)
    if lhs_1 != rhs_1:
        raise certify.Refusal("chain map fails X-boundary commutation: H_X^tgt * f_1 != f_0 * H_X^src")

    hz_src_t = codes._transpose(source_code.hz, source_code.n)
    hz_tgt_t = codes._transpose(target_code.hz, target_code.n)
    lhs_2 = _matmul_gf2(f1_rows, source_code.n, hz_src_t, len(source_code.hz))
    rhs_2 = _matmul_gf2(hz_tgt_t, len(target_code.hz), f2_rows, len(source_code.hz))
    if lhs_2 != rhs_2:
        raise certify.Refusal("chain map fails Z-boundary commutation: f_1 * (H_Z^src)^T != (H_Z^tgt)^T * f_2")

    out = {
        "source_sha256": source_code.sha256(),
        "target_sha256": target_code.sha256(),
        "commutes_boundary_1": True,
        "commutes_boundary_2": True,
    }

    if source_basis is not None and target_basis is not None:
        basis.verify_basis_receipt(source_code, source_basis)
        basis.verify_basis_receipt(target_code, target_basis)
        k_src = source_basis["k"]
        k_tgt = target_basis["k"]
        f1_t = codes._transpose(f1_rows, source_code.n)

        src_z = [basis.parse_pauli(s["z"], source_code.n)[1] for s in source_basis["slots"]]
        src_x = [basis.parse_pauli(s["x"], source_code.n)[0] for s in source_basis["slots"]]
        tgt_z = [basis.parse_pauli(s["z"], target_code.n)[1] for s in target_basis["slots"]]
        tgt_x = [basis.parse_pauli(s["x"], target_code.n)[0] for s in target_basis["slots"]]

        # Forward Z map: w_tgt = f_1(z_src[j]); its coordinate on tgt slot i is <tgt_x[i], w_tgt>
        mz = [[0] * k_src for _ in range(k_tgt)]
        ez_tgt = gf2.Echelon(target_code.hz)
        for j in range(k_src):
            w_tgt = sum((gf2.dot(f1_rows[q], src_z[j])) << q for q in range(target_code.n))
            rem = w_tgt
            for i in range(k_tgt):
                c = gf2.dot(tgt_x[i], w_tgt)
                mz[i][j] = c
                if c:
                    rem ^= tgt_z[i]
            if rem not in ez_tgt:
                raise certify.Refusal("f_1(Z_src) did not decompose onto target Z basis + stabilizers")

        # Pullback X map: w_src = f_1^T(x_tgt[i]); its coordinate on src slot j is <w_src, src_z[j]>
        mx = [[0] * k_tgt for _ in range(k_src)]
        ex_src = gf2.Echelon(source_code.hx)
        for i in range(k_tgt):
            w_src = sum((gf2.dot(f1_t[q], tgt_x[i])) << q for q in range(source_code.n))
            rem = w_src
            for j in range(k_src):
                c = gf2.dot(w_src, src_z[j])
                mx[j][i] = c
                if c:
                    rem ^= src_x[j]
            if rem not in ex_src:
                raise certify.Refusal("f_1^T(X_tgt) did not decompose onto source X basis + stabilizers")

        for i in range(k_tgt):
            for j in range(k_src):
                if mx[j][i] != mz[i][j]:
                    raise certify.Refusal("chain map adjointness failure: M_X != M_Z^T")

        out["induced_mz"] = mz
        out["induced_mx"] = mx
        out["rank_mz"] = gf2.rank([sum(mz[i][j] << j for j in range(k_src)) for i in range(k_tgt)])
    return out


def build_cyclic_extension_chain_map(
    a,
    b,
    base_lift,
    ext_factor,
    *,
    mode="pullback",
    coeff_encoding="bitmask",
):
    """Build and verify the canonical chain map for `Z_{ext_factor * base_lift} -> Z_{base_lift}` (`arXiv:2607.28621`).

    `mode="pullback"`: fiber-sum chain map `C(Z_L) -> C(Z_{m L})`.
    `mode="pushforward"`: quotient-projection chain map `C(Z_{m L}) -> C(Z_L)`.
    """
    if ext_factor < 1 or base_lift < 1:
        raise codes.SpecError("base_lift and ext_factor must be >= 1")
    ext_lift = base_lift * ext_factor

    # Normalize reduced exponents on base_lift and lift the same integer exponents into ext_lift?
    # Note: for Z_{mL} -> Z_L quotient pi(t) = t mod L to intertwine with A and B, each seed term
    # in C(Z_{mL}) uses lift of exponent e in 0..mL-1 whose mod-L reduction matches C(Z_L)!
    _, _, a_exps, _, _ = constructions._normalize_cyclic_matrix(a, base_lift, coeff_encoding, "a")
    _, _, b_exps, _, _ = constructions._normalize_cyclic_matrix(b, base_lift, coeff_encoding, "b")

    base_code = constructions.lifted_product(
        a_exps, b_exps, base_lift, coeff_encoding="exponents", name=f"LP(L={base_lift})"
    )
    ext_code = constructions.lifted_product(
        a_exps, b_exps, ext_lift, coeff_encoding="exponents", name=f"LP(L={ext_lift})"
    )

    def fiber_pullback_rows(cells):
        # Rows indexed by (cell, t_ext in 0..mL-1), bit set at (cell, t_ext % L)
        rows = []
        for c in range(cells):
            for t_ext in range(ext_lift):
                rows.append(1 << (c * base_lift + (t_ext % base_lift)))
        return rows

    def quotient_pushforward_rows(cells):
        # Rows indexed by (cell, t_base in 0..L-1), bits set at all t_ext with t_ext % L == t_base
        rows = []
        for c in range(cells):
            for t_base in range(base_lift):
                mask = 0
                for s in range(ext_factor):
                    mask |= 1 << (c * ext_lift + s * base_lift + t_base)
                rows.append(mask)
        return rows

    q_cells = base_code.n // base_lift
    hx_cells = len(base_code.hx) // base_lift
    hz_cells = len(base_code.hz) // base_lift

    base_b = basis.derive_paired_basis(base_code)
    ext_b = basis.derive_paired_basis(ext_code)

    if mode == "pullback":
        f1 = fiber_pullback_rows(q_cells)
        f0 = fiber_pullback_rows(hx_cells)
        f2 = fiber_pullback_rows(hz_cells)
        cert = verify_chain_map(base_code, ext_code, f1, f0, f2, source_basis=base_b, target_basis=ext_b)
        return base_code, ext_code, base_b, ext_b, (f1, f0, f2), cert
    if mode == "pushforward":
        f1 = quotient_pushforward_rows(q_cells)
        f0 = quotient_pushforward_rows(hx_cells)
        f2 = quotient_pushforward_rows(hz_cells)
        cert = verify_chain_map(ext_code, base_code, f1, f0, f2, source_basis=ext_b, target_basis=base_b)
        return ext_code, base_code, ext_b, base_b, (f1, f0, f2), cert
    raise codes.SpecError("mode must be 'pullback' or 'pushforward'")
