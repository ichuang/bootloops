"""Finite tensor and cyclic lifts for fast surgery protocols."""

from itertools import combinations
import random

from ftqc_codes import certify as _fc_certify, constructions as _fc_constructions
from ftqc_codes.gf2 import Echelon, dot, nullspace, supp
from ftqc_gadgets.surgery_geometry import transpose, check_rows, restricted_cone
from .algebra import hgp


def product_sequence(c, nc, d, nd, words, target, padding, degree):
    """Chang et al. Eq. (42), tensor the graph-to-classical-code map with id_D."""
    check_rows(c, nc)
    check_rows(d, nd)
    md = len(d)
    n, hx, hz = hgp(c, nc, transpose(d, nd), md)
    ops = []
    for word in words:
        if type(word) != int or not 0 < word < 1 << nc or any(dot(r, word) for r in c):
            raise ValueError("logical row is outside the classical kernel")
        sites = supp(word)
        edges = list(combinations(range(len(sites)), 2))
        incidence = [(1 << a) | (1 << b) for a, b in edges]
        # Each original check restricts to an even subset. Pair its terminals.
        dressing = [0] * len(edges)
        index = {p: i for i, p in enumerate(edges)}
        for i, row in enumerate(c):
            hits = [j for j, q in enumerate(sites) if row >> q & 1]
            for pair in zip(hits[::2], hits[1::2]):
                dressing[index[pair]] |= 1 << i
        cycles = nullspace(transpose(incidence, len(sites)), len(edges))
        av, ax, az = hgp(incidence, len(sites), transpose(d, nd), md)
        ax += [
            sum(1 << (len(sites) * md + e * nd + j) for e in supp(cy))
            for cy in cycles
            for j in range(nd)
        ]
        f1 = [1 << (q * md + j) for q in sites for j in range(md)]
        f1 += [
            sum(1 << (nc * md + i * nd + j) for i in supp(row))
            for row in dressing
            for j in range(nd)
        ]
        f0 = [
            sum(1 << (i * md + j) for i in supp(row))
            for row in dressing
            for j in range(md)
        ]
        f0 += [0] * (len(cycles) * nd)
        f2 = [1 << (q * nd + j) for q in sites for j in range(nd)]
        ops.append(
            dict(
                aux=dict(n=av, hx=ax, hz=az),
                f0=f0,
                f1=f1,
                f2=f2,
                cycles=nullspace(transpose(ax, av), len(ax)),
                metachecks=[col | (a << len(hz)) for col, a in zip(f2, az)],
            )
        )
    return dict(
        schema=1,
        kind="surgery_sequence",
        memory=dict(n=n, hx=hx, hz=hz),
        operations=ops,
        target=target,
        padding=padding,
        max_degree=degree,
    )


def lifted_product(a, b, lift):
    """C1 tensor D0 + C0 tensor D1, routed through `ftqc_codes.constructions.lifted_product`."""
    if type(lift) != int or lift < 1 or not a or not b or not a[0] or not b[0]:
        raise ValueError("empty lifted product")
    n_cols, k_cols = len(a[0]), len(b[0])
    for matrix, cols in ((a, n_cols), (b, k_cols)):
        if any(len(row) != cols for row in matrix):
            raise ValueError("ragged protograph")
        for row in matrix:
            check_rows(row, lift)
    try:
        code = _fc_constructions.lifted_product(
            a,
            b,
            lift,
            shift_convention="hx_neg_hz_pos",
            coordinate_order="cell_major",
            coeff_encoding="bitmask",
            max_n=10_000_000,
            max_rows=10_000_000,
            max_nnz=100_000_000,
        )
    except _fc_certify.Refusal as exc:
        raise ValueError("lifted square fails") from exc
    return code.n, code.hx, code.hz


def radial90():
    """[[90, 8, 10]] radial lifted-product code routed through `ftqc_codes.constructions.radial90`."""
    code = _fc_constructions.radial90()
    return code.n, code.hx, code.hz


def shift_word(word, n, lift, amount=1):
    if n % lift or word < 0 or word >> n:
        raise ValueError("invalid cyclic orbit dimensions")
    return sum(1 << (q // lift * lift + (q % lift + amount) % lift) for q in supp(word))


def invariant_logicals(checks, opposite, n, lift):
    """Orbit norm over an odd lift; independent invariant homology classes."""
    basis = Echelon(checks)
    selected = []
    for v in nullspace(opposite, n):
        orbit = 0
        for j in range(lift):
            orbit ^= shift_word(v, n, lift, j)
        if basis.add(orbit):
            selected.append(orbit)
    return selected


def radial_design():
    """A rank-one projected chain map in the completely fine radial block.

    f1 = v w^T, f0=f2=0 with invariant Z and X representatives. This is an
    explicit dense representative of the homotopy class, before optimization.
    """
    n, hx, hz = radial90()
    zs = invariant_logicals(hz, hx, n, 5)
    xs = invariant_logicals(hx, hz, n, 5)
    v = min(zs, key=int.bit_count)
    w = min((u for u in xs if dot(v, u)), key=int.bit_count)
    ports = [v if w >> q & 1 else 0 for q in range(n)]
    c = restricted_cone(n, hx, hz, ports, hx)
    if c["logical_qubits"] != 7 or len(c["measured"]) != 1:
        raise ValueError("wrong radial logical rank")
    from ftqc_gadgets.surgery_merge import word

    return dict(
        n=n,
        hx=hx,
        hz=hz,
        method="linear",
        ports=ports,
        aux=hx,
        axis="Z",
        targets=[word(n, v, "Z")],
        rounds=1,
        distance_target=2,
    )


def radial_homotopy_design(steps=240, seed=20260910):
    """Z_5-equivariant chain-homotopy sparsification in the homotopy class of `radial_design()`.

    Follows Berent, Cohen & Quintavalle (`arXiv:2609.11723` Section IV.E): searches the
    affine space of Z_5-equivariant chain maps homotopic to `f1 = v w^T` (`aux = hx`,
    `rounds = 1`, `merged_qubits = 135`, `merged_logicals = 7`), reducing maximum check
    weight, qubit degree, and physical CNOT count while measuring an equivalent Z logical.
    """
    dense = radial_design()
    n, hx, hz = dense["n"], dense["hx"], dense["hz"]
    z_ech = Echelon(hz)
    x_ech = Echelon(hx)
    dense_cone = restricted_cone(n, hx, hz, dense["ports"], hx)
    want_coset = z_ech.reduce(dense_cone["measured"][0])

    # Z_5-invariant stabilizer subspaces for reducing invariant v and w
    inv_hz = []
    for r in hz:
        orb = 0
        for s in range(5):
            orb ^= shift_word(r, n, 5, s)
        if orb:
            inv_hz.append(orb)
    inv_hx = []
    for r in hx:
        orb = 0
        for s in range(5):
            orb ^= shift_word(r, n, 5, s)
        if orb:
            inv_hx.append(orb)

    zs = invariant_logicals(hz, hx, n, 5)
    xs = invariant_logicals(hx, hz, n, 5)
    v0 = min(zs, key=int.bit_count)
    w0 = min((u for u in xs if dot(v0, u)), key=int.bit_count)

    # Greedy coset weight reduction within Z_5-invariant stabilizers
    for stabs_inv, name_tag in ((inv_hz, "v"), (inv_hx, "w")):
        cur = v0 if name_tag == "v" else w0
        improved = True
        while improved:
            improved = False
            for g in stabs_inv:
                if (cur ^ g).bit_count() < cur.bit_count():
                    cur ^= g
                    improved = True
        if name_tag == "v":
            v0 = cur
        else:
            w0 = cur

    best_ports = [v0 if (w0 >> q) & 1 else 0 for q in range(n)]
    best_cone = restricted_cone(n, hx, hz, best_ports, hx)
    if best_cone["logical_qubits"] != 7 or z_ech.reduce(best_cone["measured"][0]) != want_coset:
        best_ports = list(dense["ports"])
        best_cone = dense_cone

    def score(c_obj, p_list):
        return (
            max(c_obj["max_check"], c_obj["max_qubit_degree"]),
            c_obj["max_check"] + c_obj["max_qubit_degree"],
            sum(p.bit_count() for p in p_list),
        )

    best_sc = score(best_cone, best_ports)
    z_cycles = list(hz) + nullspace(hx, n)
    rng = random.Random(seed)

    for _ in range(steps):
        cand = list(best_ports)
        x_cell = rng.randrange(9)
        u = rng.choice(z_cycles)
        for s in range(5):
            xr = hx[x_cell * 5 + s]
            us = shift_word(u, n, 5, s)
            for q in supp(xr):
                cand[q] ^= us
        for cell in range(18):
            bw = cand[cell * 5]
            for zr in hz:
                if (bw ^ zr).bit_count() < bw.bit_count():
                    bw ^= zr
            for s in range(5):
                cand[cell * 5 + s] = shift_word(bw, n, 5, s)
        cc = restricted_cone(n, hx, hz, cand, hx)
        if (
            cc["logical_qubits"] == 7
            and len(cc["measured"]) == 1
            and z_ech.reduce(cc["measured"][0]) == want_coset
        ):
            sc = score(cc, cand)
            if sc < best_sc:
                best_ports, best_cone, best_sc = cand, cc, sc

    from ftqc_gadgets.surgery_merge import word

    return dict(
        n=n,
        hx=hx,
        hz=hz,
        method="linear",
        ports=best_ports,
        aux=hx,
        axis="Z",
        targets=[word(n, best_cone["measured"][0], "Z")],
        rounds=1,
        distance_target=2,
    )
