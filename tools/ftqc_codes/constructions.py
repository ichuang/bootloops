"""Replayable algebraic CSS code constructors and construction receipts.

Provides pure-stdlib builders for hypergraph products (`hgp`), cyclic lifted
products (`lifted_product`), general finite-group lifted products
(`group_lifted_product`), and `radial90`, plus pre-allocation size estimation
and versioned `ftqc_codes-construction-1` receipts alongside `ftqc_codes-receipt-1`.
"""
import copy
import hashlib
import json
import os

from . import certify, codes, gf2, interchange

SPEC_SCHEMA = "ftqc_codes-construction-spec-1"
RECEIPT_SCHEMA = "ftqc_codes-construction-1"

DEFAULT_MAX_N = 200_000
DEFAULT_MAX_ROWS = 200_000
DEFAULT_MAX_NNZ = 10_000_000


def _is_strict_int(x):
    return isinstance(x, int) and not isinstance(x, bool)


def _assert_no_floats(obj, path="root"):
    if isinstance(obj, float) or isinstance(obj, complex):
        raise codes.SpecError(f"non-exact float/complex value at {path}")
    if isinstance(obj, dict):
        for k, v in obj.items():
            _assert_no_floats(v, f"{path}.{k}")
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            _assert_no_floats(v, f"{path}[{i}]")


def _canonical_sha256(payload):
    _assert_no_floats(payload)
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def _bits(v, n):
    return "".join("1" if (v >> j) & 1 else "0" for j in range(n))


def _parse_seed_rows(rows, width):
    if not _is_strict_int(width) or width <= 0:
        raise codes.SpecError("seed width must be a positive integer")
    if not isinstance(rows, (list, tuple)) or not rows:
        raise codes.SpecError("seed matrix must be non-empty")
    out = []
    for r in rows:
        if _is_strict_int(r):
            if r < 0 or r >> width:
                raise codes.SpecError(f"seed row bitmask out of range for width {width}")
            out.append(int(r))
        elif isinstance(r, str):
            if len(r) != width or set(r) - {"0", "1"}:
                raise codes.SpecError(f"seed row string must be 0/1 of length {width}")
            out.append(int(r[::-1], 2))
        else:
            raise codes.SpecError("seed row must be an integer bitmask or 0/1 string")
    return out


def _normalize_cyclic_matrix(mat, lift, coeff_encoding, label):
    if not isinstance(mat, (list, tuple)) or not mat or not isinstance(mat[0], (list, tuple)) or not mat[0]:
        raise codes.SpecError(f"{label}: expected non-empty 2D seed matrix")
    nrows, ncols = len(mat), len(mat[0])
    if any(not isinstance(row, (list, tuple)) or len(row) != ncols for row in mat):
        raise codes.SpecError(f"{label}: ragged protograph matrix")
    exps_mat = []
    raw_json = []
    cancelled = 0
    for r in range(nrows):
        erow = []
        rjson = []
        for c in range(ncols):
            entry = mat[r][c]
            if coeff_encoding == "bitmask":
                if not _is_strict_int(entry) or entry < 0 or entry >> lift:
                    raise codes.SpecError(
                        f"{label}[{r}][{c}]: bitmask {entry!r} out of range for lift {lift}"
                    )
                erow.append(gf2.supp(entry))
                rjson.append(int(entry))
            elif coeff_encoding == "exponents":
                if not isinstance(entry, (list, tuple)):
                    raise codes.SpecError(
                        f"{label}[{r}][{c}]: exponents mode requires a list of exponents"
                    )
                mask = 0
                clean_entry = []
                for e in entry:
                    if not _is_strict_int(e) or not (0 <= e < lift):
                        raise codes.SpecError(
                            f"{label}[{r}][{c}]: exponent {e!r} out of range 0..{lift - 1}"
                        )
                    mask ^= 1 << e
                    clean_entry.append(int(e))
                exps = gf2.supp(mask)
                cancelled += len(clean_entry) - len(exps)
                erow.append(exps)
                rjson.append(clean_entry)
            else:
                raise codes.SpecError(f"unknown coeff_encoding {coeff_encoding!r}")
        exps_mat.append(erow)
        raw_json.append(rjson)
    return nrows, ncols, exps_mat, raw_json, cancelled


def validate_group_action(group_size, left_perms, right_perms):
    """Validate explicit left and right regular permutations on `0..group_size-1`."""
    if not _is_strict_int(group_size) or group_size < 1:
        raise codes.SpecError("group_size must be a positive integer")
    if (
        not isinstance(left_perms, (list, tuple))
        or len(left_perms) != group_size
        or not isinstance(right_perms, (list, tuple))
        or len(right_perms) != group_size
    ):
        raise codes.SpecError("left_perms and right_perms must each provide group_size permutations")
    L = [interchange.validate_permutation(p, group_size, f"left_perms[{g}]") for g, p in enumerate(left_perms)]
    R = [interchange.validate_permutation(p, group_size, f"right_perms[{g}]") for g, p in enumerate(right_perms)]

    ident = list(range(group_size))
    id_candidates = [g for g in range(group_size) if L[g] == ident and R[g] == ident]
    if len(id_candidates) != 1:
        raise certify.Refusal("group action has no unique two-sided identity permutation")
    e = id_candidates[0]

    for g in range(group_size):
        if L[g][e] != g or R[g][e] != g:
            raise certify.Refusal(
                f"group element {g} does not satisfy L_g(e) = R_g(e) = g for identity e={e}"
            )

    for g1 in range(group_size):
        for g2 in range(group_size):
            prod_l = L[g1][g2]
            prod_r = R[g2][g1]
            if prod_l != prod_r:
                raise certify.Refusal(
                    f"inconsistent left/right multiplication at ({g1}, {g2}): L={prod_l} vs R={prod_r}"
                )
            for h in range(group_size):
                if L[g1][L[g2][h]] != L[prod_l][h]:
                    raise certify.Refusal("left_perms violates group associativity / homomorphism")
                if R[g1][R[g2][h]] != R[R[g1][g2]][h]:
                    raise certify.Refusal("right_perms violates right regular anti-homomorphism")
                if L[g1][R[g2][h]] != R[g2][L[g1][h]]:
                    raise certify.Refusal("left_perms and right_perms do not commute")

    inv = [0] * group_size
    for g in range(group_size):
        invs = [h for h in range(group_size) if L[g][h] == e and L[h][g] == e]
        if len(invs) != 1:
            raise certify.Refusal(f"group element {g} lacks a two-sided inverse")
        inv[g] = invs[0]

    abelian = all(L[g] == R[g] for g in range(group_size))
    return L, R, e, inv, abelian


def _normalize_group_matrix(mat, group_size, label):
    if not isinstance(mat, (list, tuple)) or not mat or not isinstance(mat[0], (list, tuple)) or not mat[0]:
        raise codes.SpecError(f"{label}: expected non-empty 2D group-ring matrix")
    nrows, ncols = len(mat), len(mat[0])
    if any(not isinstance(row, (list, tuple)) or len(row) != ncols for row in mat):
        raise codes.SpecError(f"{label}: ragged group-ring matrix")
    supp_mat = []
    raw_json = []
    cancelled = 0
    for r in range(nrows):
        srow = []
        rjson = []
        for c in range(ncols):
            entry = mat[r][c]
            if not isinstance(entry, (list, tuple)):
                raise codes.SpecError(f"{label}[{r}][{c}]: expected list of group element indices")
            mask = 0
            clean_entry = []
            for g in entry:
                if not _is_strict_int(g) or not (0 <= g < group_size):
                    raise codes.SpecError(
                        f"{label}[{r}][{c}]: element {g!r} out of range 0..{group_size - 1}"
                    )
                mask ^= 1 << g
                clean_entry.append(int(g))
            elems = gf2.supp(mask)
            cancelled += len(clean_entry) - len(elems)
            srow.append(elems)
            rjson.append(clean_entry)
        supp_mat.append(srow)
        raw_json.append(rjson)
    return nrows, ncols, supp_mat, raw_json, cancelled


def estimate_expansion(
    spec,
    *,
    max_n=DEFAULT_MAX_N,
    max_rows=DEFAULT_MAX_ROWS,
    max_nnz=DEFAULT_MAX_NNZ,
):
    """Estimate expanded `(n, rows_hx, rows_hz, nnz_hx_upper, nnz_hz_upper, bytes_est)` before allocating."""
    if not isinstance(spec, dict) or "kind" not in spec:
        raise codes.SpecError("construction spec requires 'kind'")
    kind = spec["kind"]
    if kind == "hgp":
        na, nb = spec.get("na"), spec.get("nb")
        a = _parse_seed_rows(spec.get("a"), na)
        b = _parse_seed_rows(spec.get("b"), nb)
        ma, mb = len(a), len(b)
        n = na * nb + ma * mb
        rows_hx = ma * nb
        rows_hz = na * mb
        nnz_a = sum(gf2.wt(r) for r in a)
        nnz_b = sum(gf2.wt(r) for r in b)
        nnz_hx = nnz_a * nb + ma * nnz_b
        nnz_hz = na * nnz_b + nnz_a * mb
    elif kind == "lifted_product":
        lift = spec.get("lift")
        if not _is_strict_int(lift) or lift < 1:
            raise codes.SpecError("lifted_product requires integer lift >= 1")
        enc = spec.get("coeff_encoding", "bitmask")
        m, na, a_exps, _, _ = _normalize_cyclic_matrix(spec.get("a"), lift, enc, "a")
        r, kb, b_exps, _, _ = _normalize_cyclic_matrix(spec.get("b"), lift, enc, "b")
        n = (na * r + m * kb) * lift
        rows_hx = m * r * lift
        rows_hz = na * kb * lift
        terms_a = sum(len(a_exps[i][u]) for i in range(m) for u in range(na))
        terms_b = sum(len(b_exps[j][v]) for j in range(r) for v in range(kb))
        nnz_hx = (terms_a * r + m * terms_b) * lift
        nnz_hz = (na * terms_b + terms_a * kb) * lift
    elif kind == "group_lifted_product":
        gsize = spec.get("group_size")
        if not _is_strict_int(gsize) or gsize < 1:
            raise codes.SpecError("group_lifted_product requires positive integer group_size")
        m, na, a_supp, _, _ = _normalize_group_matrix(spec.get("a"), gsize, "a")
        r, kb, b_supp, _, _ = _normalize_group_matrix(spec.get("b"), gsize, "b")
        n = (na * r + m * kb) * gsize
        rows_hx = m * r * gsize
        rows_hz = na * kb * gsize
        terms_a = sum(len(a_supp[i][u]) for i in range(m) for u in range(na))
        terms_b = sum(len(b_supp[j][v]) for j in range(r) for v in range(kb))
        nnz_hx = (terms_a * r + m * terms_b) * gsize
        nnz_hz = (na * terms_b + terms_a * kb) * gsize
    elif kind == "radial90":
        n, rows_hx, rows_hz, nnz_hx, nnz_hz = 90, 45, 45, 270, 270
    else:
        raise codes.SpecError(f"unknown construction kind {kind!r}")

    bytes_per_row = (n + 7) // 8 + 28
    est_bytes = (rows_hx + rows_hz) * bytes_per_row
    est = {
        "kind": kind,
        "n": n,
        "rows_hx": rows_hx,
        "rows_hz": rows_hz,
        "nnz_hx_upper": nnz_hx,
        "nnz_hz_upper": nnz_hz,
        "estimated_bitrow_bytes": est_bytes,
    }
    if n > max_n or (rows_hx + rows_hz) > max_rows or (nnz_hx + nnz_hz) > max_nnz:
        raise certify.Refusal(
            f"construction expansion {est} exceeds budget "
            f"(max_n={max_n}, max_rows={max_rows}, max_nnz={max_nnz})"
        )
    return est


def hgp(
    a,
    na,
    b,
    nb,
    *,
    name=None,
    max_n=DEFAULT_MAX_N,
    max_rows=DEFAULT_MAX_ROWS,
    max_nnz=DEFAULT_MAX_NNZ,
):
    """Tillich-Zemor hypergraph product matching `codes.hgp` and `ftqc_extractors.algebra.hgp`."""
    a_rows = _parse_seed_rows(a, na)
    b_rows = _parse_seed_rows(b, nb)
    spec = {
        "schema": SPEC_SCHEMA,
        "kind": "hgp",
        "name": name or "HGP",
        "na": int(na),
        "nb": int(nb),
        "a": [_bits(r, na) for r in a_rows],
        "b": [_bits(r, nb) for r in b_rows],
    }
    est = estimate_expansion(spec, max_n=max_n, max_rows=max_rows, max_nnz=max_nnz)
    ma, mb = len(a_rows), len(b_rows)
    offset = na * nb
    at = codes._transpose(a_rows, na)
    bt = codes._transpose(b_rows, nb)
    hx = [
        sum(1 << (u * nb + j) for u in gf2.supp(a_rows[i])) | (bt[j] << (offset + i * mb))
        for i in range(ma)
        for j in range(nb)
    ]
    hz = [
        (b_rows[j] << (i * nb)) | sum(1 << (offset + u * mb + j) for u in gf2.supp(at[i]))
        for i in range(na)
        for j in range(mb)
    ]
    code = codes.CSSCode(
        name or "HGP",
        offset + ma * mb,
        hx,
        hz,
        {
            "family": "hgp",
            "construction": {
                "kind": "hgp",
                "seed_inputs": {"na": int(na), "nb": int(nb), "a": spec["a"], "b": spec["b"]},
                "generator_convention": {
                    "sector_order": "VV_left_CC_right",
                    "row_order": "row_major_kron",
                },
                "estimate": est,
            },
        },
    )
    if not code.commutes():
        raise certify.Refusal("HGP checks do not commute")
    return code


def lifted_product(
    a,
    b,
    lift,
    *,
    shift_convention="hx_neg_hz_pos",
    coordinate_order="cell_major",
    coeff_encoding="bitmask",
    name=None,
    max_n=DEFAULT_MAX_N,
    max_rows=DEFAULT_MAX_ROWS,
    max_nnz=DEFAULT_MAX_NNZ,
):
    """Cyclic lifted product over `GF(2)[x]/(x^lift - 1)`.

    `shift_convention="hx_neg_hz_pos"` matches `ftqc_extractors.surgery_products.lifted_product`
    row-for-row when `coordinate_order="cell_major"` and `coeff_encoding="bitmask"`.
    """
    if not _is_strict_int(lift) or lift < 1:
        raise codes.SpecError("lift must be a positive integer")
    if shift_convention not in ("hx_neg_hz_pos", "hx_pos_hz_neg"):
        raise codes.SpecError(f"unknown shift_convention {shift_convention!r}")
    if coordinate_order not in ("cell_major", "fiber_major"):
        raise codes.SpecError(f"unknown coordinate_order {coordinate_order!r}")

    m, na, a_exps, a_json, canc_a = _normalize_cyclic_matrix(a, lift, coeff_encoding, "a")
    r, kb, b_exps, b_json, canc_b = _normalize_cyclic_matrix(b, lift, coeff_encoding, "b")
    spec = {
        "schema": SPEC_SCHEMA,
        "kind": "lifted_product",
        "name": name or f"LP(L={lift})",
        "lift": lift,
        "coeff_encoding": coeff_encoding,
        "shift_convention": shift_convention,
        "coordinate_order": coordinate_order,
        "a": a_json,
        "b": b_json,
    }
    est = estimate_expansion(spec, max_n=max_n, max_rows=max_rows, max_nnz=max_nnz)

    total_cells = na * r + m * kb
    width = total_cells * lift
    hx_rev = shift_convention == "hx_neg_hz_pos"
    hz_rev = not hx_rev

    def q_idx(cell, fiber):
        if coordinate_order == "cell_major":
            return cell * lift + fiber
        return fiber * total_cells + cell

    def orbit(cell, shift, exps, reverse):
        v = 0
        for e in exps:
            fib = (shift - e if reverse else shift + e) % lift
            v ^= 1 << q_idx(cell, fib)
        return v

    hx = []
    for i in range(m):
        for j in range(r):
            for t in range(lift):
                row = 0
                for u in range(na):
                    row ^= orbit(u * r + j, t, a_exps[i][u], hx_rev)
                for v in range(kb):
                    row ^= orbit(na * r + i * kb + v, t, b_exps[j][v], hx_rev)
                hx.append(row)

    hz = []
    for u in range(na):
        for v in range(kb):
            for t in range(lift):
                row = 0
                for j in range(r):
                    row ^= orbit(u * r + j, t, b_exps[j][v], hz_rev)
                for i in range(m):
                    row ^= orbit(na * r + i * kb + v, t, a_exps[i][u], hz_rev)
                hz.append(row)

    code = codes.CSSCode(
        name or f"LP(L={lift})",
        width,
        hx,
        hz,
        {
            "family": "lifted_product",
            "construction": {
                "kind": "lifted_product",
                "seed_inputs": {
                    "lift": lift,
                    "shape_a": [m, na],
                    "shape_b": [r, kb],
                    "coeff_encoding": coeff_encoding,
                    "a": a_json,
                    "b": b_json,
                    "a_reduced_exponents": a_exps,
                    "b_reduced_exponents": b_exps,
                    "cancelled_terms": canc_a + canc_b,
                },
                "generator_convention": {
                    "ring": f"GF(2)[x]/(x^{lift}-1)",
                    "shift_convention": shift_convention,
                    "coordinate_order": coordinate_order,
                    "row_order": "block_cell_then_fiber",
                },
                "estimate": est,
            },
        },
    )
    if not code.commutes():
        raise certify.Refusal("lifted product checks do not commute")
    return code


def lifted_product_circulant_ref(
    a,
    b,
    lift,
    *,
    shift_convention="hx_neg_hz_pos",
    coordinate_order="cell_major",
    coeff_encoding="bitmask",
    name=None,
):
    """Independent reference builder assembling dense `lift x lift` circulant blocks."""
    m, na, a_exps, _, _ = _normalize_cyclic_matrix(a, lift, coeff_encoding, "a")
    r, kb, b_exps, _, _ = _normalize_cyclic_matrix(b, lift, coeff_encoding, "b")
    hx_rev = shift_convention == "hx_neg_hz_pos"
    hz_rev = not hx_rev

    def circulant_block(exps, reverse):
        blk = [[0] * lift for _ in range(lift)]
        for e in exps:
            for row_t in range(lift):
                col_s = (row_t - e if reverse else row_t + e) % lift
                blk[row_t][col_s] ^= 1
        return blk

    total_cells = na * r + m * kb
    width = total_cells * lift

    def col_map(cell, fib):
        return cell * lift + fib if coordinate_order == "cell_major" else fib * total_cells + cell

    hx_grid = [[0] * width for _ in range(m * r * lift)]
    for i in range(m):
        for j in range(r):
            r_base = (i * r + j) * lift
            for u in range(na):
                blk = circulant_block(a_exps[i][u], hx_rev)
                c_cell = u * r + j
                for rt in range(lift):
                    for cs in range(lift):
                        if blk[rt][cs]:
                            hx_grid[r_base + rt][col_map(c_cell, cs)] ^= 1
            for v in range(kb):
                blk = circulant_block(b_exps[j][v], hx_rev)
                c_cell = na * r + i * kb + v
                for rt in range(lift):
                    for cs in range(lift):
                        if blk[rt][cs]:
                            hx_grid[r_base + rt][col_map(c_cell, cs)] ^= 1

    hz_grid = [[0] * width for _ in range(na * kb * lift)]
    for u in range(na):
        for v in range(kb):
            r_base = (u * kb + v) * lift
            for j in range(r):
                blk = circulant_block(b_exps[j][v], hz_rev)
                c_cell = u * r + j
                for rt in range(lift):
                    for cs in range(lift):
                        if blk[rt][cs]:
                            hz_grid[r_base + rt][col_map(c_cell, cs)] ^= 1
            for i in range(m):
                blk = circulant_block(a_exps[i][u], hz_rev)
                c_cell = na * r + i * kb + v
                for rt in range(lift):
                    for cs in range(lift):
                        if blk[rt][cs]:
                            hz_grid[r_base + rt][col_map(c_cell, cs)] ^= 1

    return codes.CSSCode(
        name or f"LP_ref(L={lift})",
        width,
        gf2.rows_from_lists(hx_grid),
        gf2.rows_from_lists(hz_grid),
    )


def group_lifted_product(
    group_size,
    left_perms,
    right_perms,
    a,
    b,
    *,
    mult_convention="hx_inv_left_hz_left",
    coordinate_order="cell_major",
    name=None,
    max_n=DEFAULT_MAX_N,
    max_rows=DEFAULT_MAX_ROWS,
    max_nnz=DEFAULT_MAX_NNZ,
):
    """Panteleev-Kalachev lifted product over a finite group ring `GF(2)[G]`.

    Validates that `left_perms` and `right_perms` form commuting left/right
    regular representations of a finite group `G` before expanding.
    """
    if mult_convention not in ("hx_inv_left_hz_left", "hx_left_hz_inv_left"):
        raise codes.SpecError(f"unknown mult_convention {mult_convention!r}")
    if coordinate_order not in ("cell_major", "fiber_major"):
        raise codes.SpecError(f"unknown coordinate_order {coordinate_order!r}")

    L, R, ident_e, inv, abelian = validate_group_action(group_size, left_perms, right_perms)
    m, na, a_supp, a_json, canc_a = _normalize_group_matrix(a, group_size, "a")
    r, kb, b_supp, b_json, canc_b = _normalize_group_matrix(b, group_size, "b")

    spec = {
        "schema": SPEC_SCHEMA,
        "kind": "group_lifted_product",
        "name": name or f"GroupLP(|G|={group_size})",
        "group_size": group_size,
        "left_perms": L,
        "right_perms": R,
        "mult_convention": mult_convention,
        "coordinate_order": coordinate_order,
        "a": a_json,
        "b": b_json,
    }
    est = estimate_expansion(spec, max_n=max_n, max_rows=max_rows, max_nnz=max_nnz)

    total_cells = na * r + m * kb
    width = total_cells * group_size
    inv_on_hx = mult_convention == "hx_inv_left_hz_left"

    def q_idx(cell, g_elem):
        if coordinate_order == "cell_major":
            return cell * group_size + g_elem
        return g_elem * total_cells + cell

    hx = []
    for i in range(m):
        for j in range(r):
            for t in range(group_size):
                row = 0
                for u in range(na):
                    cell = u * r + j
                    for g in a_supp[i][u]:
                        gg = inv[g] if inv_on_hx else g
                        row ^= 1 << q_idx(cell, L[gg][t])
                for v in range(kb):
                    cell = na * r + i * kb + v
                    for g in b_supp[j][v]:
                        gg = inv[g] if inv_on_hx else g
                        row ^= 1 << q_idx(cell, R[gg][t])
                hx.append(row)

    hz = []
    for u in range(na):
        for v in range(kb):
            for t in range(group_size):
                row = 0
                for j in range(r):
                    cell = u * r + j
                    for g in b_supp[j][v]:
                        gg = g if inv_on_hx else inv[g]
                        row ^= 1 << q_idx(cell, R[gg][t])
                for i in range(m):
                    cell = na * r + i * kb + v
                    for g in a_supp[i][u]:
                        gg = g if inv_on_hx else inv[g]
                        row ^= 1 << q_idx(cell, L[gg][t])
                hz.append(row)

    code = codes.CSSCode(
        name or f"GroupLP(|G|={group_size})",
        width,
        hx,
        hz,
        {
            "family": "group_lifted_product",
            "construction": {
                "kind": "group_lifted_product",
                "seed_inputs": {
                    "group_size": group_size,
                    "identity": ident_e,
                    "abelian": abelian,
                    "left_perms": L,
                    "right_perms": R,
                    "shape_a": [m, na],
                    "shape_b": [r, kb],
                    "a": a_json,
                    "b": b_json,
                    "a_reduced_elements": a_supp,
                    "b_reduced_elements": b_supp,
                    "cancelled_terms": canc_a + canc_b,
                },
                "generator_convention": {
                    "ring": f"GF(2)[G], |G|={group_size}",
                    "mult_convention": mult_convention,
                    "coordinate_order": coordinate_order,
                    "row_order": "block_cell_then_group_element",
                },
                "estimate": est,
            },
        },
    )
    if not code.commutes():
        raise certify.Refusal("group lifted product checks do not commute")
    return code


def radial90(name="Radial90"):
    """90-qubit cyclic lifted-product radial code (`n=90, k=8`, digest `b656926c5959`)."""
    a = [[4, 8, 16], [2, 16, 2], [16, 8, 4]]
    b = [[4, 4, 1], [16, 1, 16], [2, 8, 1]]
    code = lifted_product(a, b, 5, name=name)
    code.meta["family"] = "radial90"
    code.meta["construction"]["provenance"] = {
        "paper": "Scruby, Hillmann, Pesah (radial codes) / He et al. (lifted surgery)",
        "arxiv": ["2406.14445", "2609.11723"],
        "locator": "3x3 protographs over GF(2)[x]/(x^5-1)",
    }
    return code


def from_spec(
    spec_or_path,
    *,
    base_dir=None,
    max_n=DEFAULT_MAX_N,
    max_rows=DEFAULT_MAX_ROWS,
    max_nnz=DEFAULT_MAX_NNZ,
):
    """Build a `CSSCode` from a tagged construction JSON spec (`ftqc_codes-construction-spec-1`)."""
    if isinstance(spec_or_path, str):
        full_path = os.path.abspath(spec_or_path)
        if base_dir is None:
            base_dir = os.path.dirname(full_path)
        with open(full_path) as f:
            spec = json.load(f)
    elif isinstance(spec_or_path, dict):
        spec = spec_or_path
    else:
        raise codes.SpecError("construction spec must be a dict or path")

    _assert_no_floats(spec)
    if spec.get("schema") != SPEC_SCHEMA:
        raise codes.SpecError(f"expected schema {SPEC_SCHEMA!r}, got {spec.get('schema')!r}")
    kind = spec.get("kind")
    nm = spec.get("name")
    if kind == "hgp":
        code = hgp(
            spec.get("a"),
            spec.get("na"),
            spec.get("b"),
            spec.get("nb"),
            name=nm,
            max_n=max_n,
            max_rows=max_rows,
            max_nnz=max_nnz,
        )
    elif kind == "lifted_product":
        code = lifted_product(
            spec.get("a"),
            spec.get("b"),
            spec.get("lift"),
            shift_convention=spec.get("shift_convention", "hx_neg_hz_pos"),
            coordinate_order=spec.get("coordinate_order", "cell_major"),
            coeff_encoding=spec.get("coeff_encoding", "bitmask"),
            name=nm,
            max_n=max_n,
            max_rows=max_rows,
            max_nnz=max_nnz,
        )
    elif kind == "group_lifted_product":
        code = group_lifted_product(
            spec.get("group_size"),
            spec.get("left_perms"),
            spec.get("right_perms"),
            spec.get("a"),
            spec.get("b"),
            mult_convention=spec.get("mult_convention", "hx_inv_left_hz_left"),
            coordinate_order=spec.get("coordinate_order", "cell_major"),
            name=nm,
            max_n=max_n,
            max_rows=max_rows,
            max_nnz=max_nnz,
        )
    elif kind == "radial90":
        estimate_expansion(spec, max_n=max_n, max_rows=max_rows, max_nnz=max_nnz)
        code = radial90(name=nm or "Radial90")
    elif kind == "cpm_pair_partition":
        from . import families_highrate
        code = families_highrate.cpm_pair_partition(
            spec.get("lift"),
            spec.get("ex_shifts"),
            spec.get("ez_shifts"),
            code_type=spec.get("code_type", "css"),
            verify_pair_partition=spec.get("verify_pair_partition", True),
            name=nm,
            provenance=spec.get("provenance"),
            paper_claims=spec.get("paper_claims"),
        )
    elif kind == "named_family":
        fam_spec = spec.get("family_spec")
        if not isinstance(fam_spec, str) or not fam_spec:
            raise codes.SpecError("named_family requires string 'family_spec'")
        code = codes.load(fam_spec)
        if nm:
            code.name = nm
        code.meta["construction"] = {
            "kind": "named_family",
            "seed_inputs": {"family_spec": fam_spec},
            "generator_convention": {"loader": "ftqc_codes.codes.load"},
        }
    else:
        raise codes.SpecError(f"unknown construction kind {kind!r}")

    if "provenance" in spec and spec["provenance"] is not None:
        code.meta.setdefault("construction", {})["provenance"] = copy.deepcopy(spec["provenance"])
    if "paper_claims" in spec and spec["paper_claims"] is not None:
        code.meta.setdefault("construction", {})["paper_claims"] = copy.deepcopy(spec["paper_claims"])
    return code


def weight_distributions(rows, n):
    """Compute `(row_weight_hist, col_weight_hist, max_row_weight)` with string keys."""
    r_counts = {}
    c_weights = [0] * n
    max_rw = 0
    for r in rows:
        w = gf2.wt(r)
        if w > max_rw:
            max_rw = w
        r_counts[w] = r_counts.get(w, 0) + 1
        for j in gf2.supp(r):
            c_weights[j] += 1
    c_counts = {}
    for cw in c_weights:
        c_counts[cw] = c_counts.get(cw, 0) + 1
    row_hist = {str(k): r_counts[k] for k in sorted(r_counts)}
    col_hist = {str(k): c_counts[k] for k in sorted(c_counts)}
    return row_hist, col_hist, max_rw


def _reconstruct_spec_from_code(code):
    cmeta = code.meta.get("construction")
    if cmeta:
        kind = cmeta["kind"]
        seeds = cmeta["seed_inputs"]
        conv = cmeta["generator_convention"]
        if kind == "hgp":
            return {
                "schema": SPEC_SCHEMA,
                "kind": "hgp",
                "name": code.name,
                "na": seeds["na"],
                "nb": seeds["nb"],
                "a": seeds["a"],
                "b": seeds["b"],
            }
        if kind == "lifted_product":
            return {
                "schema": SPEC_SCHEMA,
                "kind": "lifted_product",
                "name": code.name,
                "lift": seeds["lift"],
                "coeff_encoding": seeds["coeff_encoding"],
                "shift_convention": conv["shift_convention"],
                "coordinate_order": conv["coordinate_order"],
                "a": seeds["a"],
                "b": seeds["b"],
            }
        if kind == "group_lifted_product":
            return {
                "schema": SPEC_SCHEMA,
                "kind": "group_lifted_product",
                "name": code.name,
                "group_size": seeds["group_size"],
                "left_perms": seeds["left_perms"],
                "right_perms": seeds["right_perms"],
                "mult_convention": conv["mult_convention"],
                "coordinate_order": conv["coordinate_order"],
                "a": seeds["a"],
                "b": seeds["b"],
            }
        if kind == "named_family":
            return {
                "schema": SPEC_SCHEMA,
                "kind": "named_family",
                "name": code.name,
                "family_spec": seeds["family_spec"],
            }
    if code.meta.get("family") == "bb":
        m = code.meta
        fam_s = f"bb:{m['l']},{m['m']},{'+'.join(m['A'])},{'+'.join(m['B'])}"
        return {
            "schema": SPEC_SCHEMA,
            "kind": "named_family",
            "name": code.name,
            "family_spec": fam_s,
        }
    if "cpm_pp" in code.meta:
        cpp = code.meta["cpm_pp"]
        return {
            "schema": SPEC_SCHEMA,
            "kind": "cpm_pair_partition",
            "name": code.name,
            "lift": cpp["lift"],
            "ex_shifts": cpp["ex_shifts"],
            "ez_shifts": cpp["ez_shifts"],
            "verify_pair_partition": cpp.get("pair_partition_verified", True),
        }
    named_alias = {
        "APM1152": "apm1152",
        "QE_D8_36_8_3": "qt_d8_36",
        "QE_S3_288_8": "qt_s3_288",
    }.get(code.name)
    if named_alias and named_alias in codes.NAMED:
        return {
            "schema": SPEC_SCHEMA,
            "kind": "named_family",
            "name": code.name,
            "family_spec": named_alias,
        }
    return None


def build_construction_receipt(
    code_or_spec,
    *,
    static_receipt=None,
    provenance=None,
    paper_claims=None,
    simulation_estimates=None,
    reproduction_spec=None,
    base_dir=None,
):
    """Build a versioned `ftqc_codes-construction-1` receipt alongside `ftqc_codes-receipt-1`."""
    if isinstance(code_or_spec, codes.CSSCode):
        code = code_or_spec
    else:
        code = from_spec(code_or_spec, base_dir=base_dir)

    if not code.commutes():
        raise certify.Refusal("H_X and H_Z do not commute")

    ex, ez = gf2.Echelon(code.hx), gf2.Echelon(code.hz)
    k = code.n - ex.rank - ez.rank
    rw_hx, cw_hx, max_hx = weight_distributions(code.hx, code.n)
    rw_hz, cw_hz, max_hz = weight_distributions(code.hz, code.n)

    imp = code.meta.get("import")
    replay_spec = reproduction_spec or _reconstruct_spec_from_code(code)

    if replay_spec is not None:
        rebuilt = from_spec(replay_spec, base_dir=base_dir)
        if rebuilt.sha256() != code.sha256():
            raise certify.Refusal(
                f"reproduction_spec produces sha256 {rebuilt.sha256()}, expected {code.sha256()}"
            )
        ancestry = "reproduced"
        cmeta = rebuilt.meta.get("construction", {})
        seed_inputs = copy.deepcopy(cmeta.get("seed_inputs", replay_spec))
        gen_conv = copy.deepcopy(cmeta.get("generator_convention", {}))
        source_hashes = copy.deepcopy(imp["source_hashes"]) if imp else None
    elif imp is not None:
        ancestry = "supplied"
        seed_inputs = {
            "hx_supports": [gf2.supp(r) for r in code.hx],
            "hz_supports": [gf2.supp(r) for r in code.hz],
        }
        gen_conv = {
            "producer": copy.deepcopy(imp["producer"]),
            "index_base": imp["index_base"],
            "coordinate_map": copy.deepcopy(imp["coordinate_map"]),
        }
        source_hashes = copy.deepcopy(imp["source_hashes"])
    else:
        ancestry = "supplied"
        seed_inputs = {
            "hx_supports": [gf2.supp(r) for r in code.hx],
            "hz_supports": [gf2.supp(r) for r in code.hz],
        }
        gen_conv = {"producer": {"name": "direct_css_matrices", "version": certify.VERSION}}
        source_hashes = {"code_sha256": code.sha256()}

    prov = copy.deepcopy(
        provenance
        or (code.meta.get("construction") or {}).get("provenance")
        or code.meta.get("provenance")
        or (imp or {}).get("provenance")
        or {"paper": "unspecified", "version": "unspecified", "locator": "unspecified"}
    )
    p_claims = copy.deepcopy(
        paper_claims
        if paper_claims is not None
        else (code.meta.get("construction") or {}).get("paper_claims")
        or code.meta.get("paper_claims")
        or (imp or {}).get("paper_claims")
        or {}
    )
    sim_est = copy.deepcopy(
        simulation_estimates
        if simulation_estimates is not None
        else (imp or {}).get("simulation_estimates")
        or []
    )
    _assert_no_floats(p_claims, "paper_claims")
    _assert_no_floats(sim_est, "simulation_estimates")

    witness_upper = {"d": None, "d_x": None, "d_z": None, "witnesses": {}}
    proved_lower = {"d": None, "d_x": None, "d_z": None, "legs": []}
    static_ref = None

    for q, winfo in (code.meta.get("witnesses") or {}).items():
        if q in ("d_x", "d_z"):
            witness_upper[q] = winfo["weight"]
            witness_upper["witnesses"][q] = {
                "weight": winfo["weight"],
                "support": list(winfo["support"]),
                "partner_support": list(winfo["partner_support"]),
            }

    if static_receipt is not None:
        certify.verify(static_receipt)
        if static_receipt["code"]["sha256"] != code.sha256():
            raise certify.Refusal("static_receipt code sha256 does not match constructed code")
        sf = static_receipt["facts"]
        for q in ("d_x", "d_z"):
            if sf[q]["upper"] is not None:
                w_int = certify._unbits(sf[q]["witness"])
                p_int = certify._unbits(sf[q]["partner"])
                if witness_upper[q] is None or sf[q]["upper"] <= witness_upper[q]:
                    witness_upper[q] = sf[q]["upper"]
                    witness_upper["witnesses"][q] = {
                        "weight": sf[q]["upper"],
                        "support": gf2.supp(w_int),
                        "partner_support": gf2.supp(p_int),
                    }
            if sf[q]["lower"] is not None:
                proved_lower[q] = sf[q]["lower"]
        proved_lower["d"] = sf["d"]["lower"]
        proved_lower["legs"] = sorted(
            {ev["leg"] for ev in static_receipt["evidence"] if ev.get("status") == "proved"}
        )
        static_ref = {
            "schema": static_receipt["schema"],
            "sha256": _canonical_sha256(static_receipt),
        }

    ups = [u for u in (witness_upper["d_x"], witness_upper["d_z"]) if u is not None]
    witness_upper["d"] = min(ups) if ups else None

    construction_payload = {
        "name": code.name,
        "ancestry": ancestry,
        "provenance": prov,
        "seed_inputs": seed_inputs,
        "generator_convention": gen_conv,
        "replay_spec": replay_spec,
        "source_hashes": source_hashes,
    }
    c_sha = _canonical_sha256(construction_payload)

    return {
        "schema": RECEIPT_SCHEMA,
        "construction_sha256": c_sha,
        "payload": construction_payload,
        "code_sha256": code.sha256(),
        "metrics": {
            "n": code.n,
            "rows_hx": len(code.hx),
            "rows_hz": len(code.hz),
            "rank_hx": ex.rank,
            "rank_hz": ez.rank,
            "k": k,
            "rate": [k, code.n],
            "max_check_weight": max(max_hx, max_hz),
            "row_weights_hx": rw_hx,
            "col_weights_hx": cw_hx,
            "row_weights_hz": rw_hz,
            "col_weights_hz": cw_hz,
        },
        "claims": {
            "paper_claims": p_claims,
            "witness_upper": witness_upper,
            "proved_lower": proved_lower,
            "simulation_estimates": sim_est,
        },
        "static_receipt_ref": static_ref,
    }


def verify_construction_receipt(crec, *, static_receipt=None, base_dir=None):
    """Re-verify a `ftqc_codes-construction-1` receipt by re-expanding and checking invariants."""
    if not isinstance(crec, dict) or crec.get("schema") != RECEIPT_SCHEMA:
        raise certify.Refusal("unknown construction receipt schema")
    _assert_no_floats(crec)

    payload = crec.get("payload")
    if not isinstance(payload, dict):
        raise certify.Refusal("construction receipt missing payload")
    want_csha = _canonical_sha256(payload)
    if crec.get("construction_sha256") != want_csha:
        raise certify.Refusal("construction_sha256 does not match canonical payload")

    lines = [f"construction_sha256 {want_csha[:16]} matches payload"]
    ancestry = payload.get("ancestry")
    if ancestry not in ("reproduced", "supplied"):
        raise certify.Refusal(f"invalid ancestry {ancestry!r}")

    replay_spec = payload.get("replay_spec")
    if ancestry == "reproduced":
        if not isinstance(replay_spec, dict):
            raise certify.Refusal("reproduced ancestry requires replay_spec")
        code = from_spec(replay_spec, base_dir=base_dir)
        lines.append(f"re-expanded generator {replay_spec.get('kind')}")
    else:
        if not payload.get("source_hashes"):
            raise certify.Refusal("supplied ancestry requires source_hashes")
        seeds = payload.get("seed_inputs") or {}
        m_meta = crec.get("metrics") or {}
        n_decl = m_meta.get("n")
        hx_supp = seeds.get("hx_supports")
        hz_supp = seeds.get("hz_supports")
        if not isinstance(hx_supp, list) or not isinstance(hz_supp, list) or not _is_strict_int(n_decl):
            raise certify.Refusal("supplied construction receipt lacks sparse matrix supports")
        hx = interchange.parse_sparse_supports(hx_supp, len(hx_supp), n_decl, label="hx")
        hz = interchange.parse_sparse_supports(hz_supp, len(hz_supp), n_decl, label="hz")
        code = codes.CSSCode(payload.get("name", "supplied"), n_decl, hx, hz)
        lines.append("reloaded supplied sparse matrices")

    if code.sha256() != crec.get("code_sha256"):
        raise certify.Refusal(
            f"re-expanded code sha256 {code.sha256()} != receipt {crec.get('code_sha256')}"
        )
    if not code.commutes():
        raise certify.Refusal("re-expanded checks do not commute")
    lines.append(f"code_sha256 {code.sha256()[:16]} and commutation verified")

    ex, ez = gf2.Echelon(code.hx), gf2.Echelon(code.hz)
    k = code.n - ex.rank - ez.rank
    rw_hx, cw_hx, max_hx = weight_distributions(code.hx, code.n)
    rw_hz, cw_hz, max_hz = weight_distributions(code.hz, code.n)

    expected_metrics = {
        "n": code.n,
        "rows_hx": len(code.hx),
        "rows_hz": len(code.hz),
        "rank_hx": ex.rank,
        "rank_hz": ez.rank,
        "k": k,
        "rate": [k, code.n],
        "max_check_weight": max(max_hx, max_hz),
        "row_weights_hx": rw_hx,
        "col_weights_hx": cw_hx,
        "row_weights_hz": rw_hz,
        "col_weights_hz": cw_hz,
    }
    m = crec.get("metrics")
    if not isinstance(m, dict):
        raise certify.Refusal("missing metrics block")
    rate = m.get("rate")
    if not isinstance(rate, list) or len(rate) != 2 or not all(_is_strict_int(x) for x in rate):
        raise certify.Refusal("metrics.rate must be an exact [k, n] integer pair")
    for key, val in expected_metrics.items():
        if m.get(key) != val:
            raise certify.Refusal(f"metric {key}: receipt has {m.get(key)!r}, recomputed {val!r}")
    lines.append(f"metrics verified: n={code.n} k={k} rate={k}/{code.n} max_wt={expected_metrics['max_check_weight']}")

    claims = crec.get("claims")
    if not isinstance(claims, dict) or set(claims) != {
        "paper_claims",
        "witness_upper",
        "proved_lower",
        "simulation_estimates",
    }:
        raise certify.Refusal("claims block must separate paper_claims, witness_upper, proved_lower, simulation_estimates")

    wu = claims["witness_upper"]
    wmap = wu.get("witnesses") or {}
    for q in ("d_x", "d_z"):
        u = wu.get(q)
        if u is not None:
            if q not in wmap:
                raise certify.Refusal(f"witness_upper.{q}={u} lacks witness support")
            winfo = interchange._parse_witness_entry(
                {
                    "index_base": 0,
                    "coordinates": "native",
                    "support": wmap[q]["support"],
                    "partner_support": wmap[q].get("partner_support"),
                },
                code,
                q,
                list(range(code.n)),
            )
            if winfo["weight"] != u:
                raise certify.Refusal(f"witness_upper.{q} weight mismatch")
            lines.append(f"witness {q} <= {u} verified")
    ups = [u for u in (wu.get("d_x"), wu.get("d_z")) if u is not None]
    if wu.get("d") != (min(ups) if ups else None):
        raise certify.Refusal("witness_upper.d inconsistent with d_x, d_z")

    pl = claims["proved_lower"]
    has_lower = any(pl.get(q) is not None for q in ("d", "d_x", "d_z"))
    sref = crec.get("static_receipt_ref")
    if has_lower and sref is None:
        raise certify.Refusal("proved_lower claims a bound without static_receipt_ref")
    for q in ("d", "d_x", "d_z"):
        if pl.get(q) is not None and wu.get(q) is not None and pl[q] > wu[q]:
            raise certify.Refusal(f"proved_lower.{q}={pl[q]} exceeds witness_upper.{q}={wu[q]}")

    if static_receipt is not None:
        certify.verify(static_receipt)
        if sref is None or _canonical_sha256(static_receipt) != sref.get("sha256"):
            raise certify.Refusal("supplied static_receipt does not match static_receipt_ref")
        if static_receipt["code"]["sha256"] != code.sha256():
            raise certify.Refusal("supplied static_receipt is for a different code_sha256")
        sf = static_receipt["facts"]
        for q in ("d_x", "d_z"):
            if pl.get(q) != sf[q]["lower"]:
                raise certify.Refusal(f"proved_lower.{q} disagrees with static_receipt")
        if pl.get("d") != sf["d"]["lower"]:
            raise certify.Refusal("proved_lower.d disagrees with static_receipt")
        lines.append("linked static_receipt verified")

    return lines
