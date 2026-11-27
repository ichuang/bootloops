"""Hierarchical decoding of a concatenated CSS code, and what it corrects.

Setting: an outer CSS code on N positions, an inner CSS code with k = 1 on m qubits,
and their concatenation with the row layout of `codes.concatenate` (the inner checks
of block 0, ..., block N - 1, then the outer checks lifted through an inner logical
representative).  A column is the m qubits of one outer position.

The decoder for one sector (X patterns from Z-check values, or the mirror image):
  1. in each column, the lightest pattern with that column's inner check values;
  2. the outer check values are corrected for what step 1 did to the lifted checks
     and decoded by the lightest outer pattern;
  3. every column in that outer pattern receives an inner logical representative.

Claim (`certificate`).  Let t1 and t2 be radii within which the lightest inner and
outer patterns are unique.  If at most t2 columns hold more than t1 errors, the
decoder returns the error up to a stabilizer.
Proof.  In a column with at most t1 errors step 1 returns the error itself.  In any
other column the error plus the step-1 pattern has zero inner check values, so it is
an inner stabilizer or an inner stabilizer times the inner logical (k = 1).  The set
of columns of the second kind is contained in the bad columns, so it has at most t2
elements; the corrected outer check values are exactly the outer check values of
that set, and step 2 returns it.  Step 3 then leaves an inner stabilizer in every
column.
"""
from itertools import combinations
from math import comb

from . import certify, concat, gf2, minweight


def leader_table(checks, npos, cap=5_000_000):
    """{check values: position mask of a lightest pattern}, lightest first, ties by
    position order, and the radius within which the lightest pattern is unique.
    Raises ValueError when the table would need more than `cap` patterns of one weight."""
    col = [sum(1 << i for i, r in enumerate(checks) if (r >> p) & 1) for p in range(npos)]
    need = 1 << gf2.rank(checks)
    tab, w, clash = {0: 0}, 0, None
    while len(tab) < need and w < npos:
        w += 1
        if comb(npos, w) > cap:
            raise ValueError(f"decoder table needs patterns of weight {w} on {npos} positions")
        for pos in combinations(range(npos), w):
            s = m = 0
            for p in pos:
                s ^= col[p]
                m |= 1 << p
            if s in tab:
                if clash is None:
                    clash = w
            else:
                tab[s] = m
    if clash is None and len(tab) >= need and w < npos:
        clash = w + 1
    return tab, (w if clash is None else clash - 1)


class Hier:
    """Hierarchical decoder of one sector: 'X' corrects X patterns (Z checks), 'Z' the
    mirror image.  `checks` are the concatenated code's checks of that sector in the
    layout the decoder reads."""

    def __init__(self, outer, inner, sector, cap=5_000_000):
        if certify.structure(inner)["k"] != 1:
            raise certify.Refusal("hierarchical decoding needs an inner code with k = 1")
        self.N, self.m, self.sector = outer.n, inner.n, sector
        lx = concat._logical_rep(inner.hz, inner.hx, inner.n)
        lz = concat._logical_rep(inner.hx, inner.hz, inner.n)
        hx, hz = concat.expected_stabilizers(outer, inner)
        if sector == "X":
            self.hi, self.ho, self.same, self.dual, self.checks = inner.hz, outer.hz, lx, lz, hz
        elif sector == "Z":
            self.hi, self.ho, self.same, self.dual, self.checks = inner.hx, outer.hx, lz, lx, hx
        else:
            raise ValueError("sector is X or Z")
        self.inner_table, self.t1 = leader_table(self.hi, self.m, cap)
        self.outer_table, self.t2 = leader_table(self.ho, self.N, cap)

    def check_values(self, pattern):
        return sum(1 << i for i, r in enumerate(self.checks) if gf2.dot(r, pattern))

    def decode(self, values):
        """Pattern for the given check values (int over `checks`); None if a table has
        no entry."""
        ri = len(self.hi)
        mask = (1 << ri) - 1
        pattern = flipped = 0
        for i in range(self.N):
            e = self.inner_table.get((values >> (i * ri)) & mask)
            if e is None:
                return None
            pattern |= e << (i * self.m)
            if gf2.dot(self.dual, e):
                flipped |= 1 << i
        osyn = values >> (self.N * ri)
        for j, row in enumerate(self.ho):
            if gf2.dot(row, flipped):
                osyn ^= 1 << j
        o = self.outer_table.get(osyn)
        if o is None:
            return None
        for i in gf2.supp(o):
            pattern ^= self.same << (i * self.m)
        return pattern

    def bad_columns(self, pattern, t1=None):
        t1 = self.t1 if t1 is None else t1
        mask = (1 << self.m) - 1
        return [i for i in range(self.N) if ((pattern >> (i * self.m)) & mask).bit_count() > t1]


def kernel_radius(checks, n, requested, cap=5_000_000):
    """Prove uniqueness through the requested weight using the actual decoder matrix.

    First by search over the kernel.  If the search budget stops short of the request and
    the kernel is a cyclic code, the BCH bound is tried on it: its zeros are computed
    exactly from this matrix and the theorem is cited (bch.py).  The evidence says which
    of the two the bound rests on."""
    kernel = gf2.nullspace(checks, n)
    search = minweight.min_weight(kernel, n, target=2*requested+1, max_combos=cap)
    bound = n if search["empty"] else (search["lower"] - 1) // 2
    evidence = dict(kernel_lower=search["lower"], empty_kernel=search["empty"])
    if bound < requested and not search["empty"]:
        from . import bch
        try:
            r = bch.code_bound(kernel, n)
        except (bch.NotCyclic, ValueError):
            r = None
        if r is not None and (r["bound"] - 1) // 2 > bound:
            bound = (r["bound"] - 1) // 2
            evidence.update(kernel_lower=r["bound"], search_lower=search["lower"],
                            bound_from="BCH bound on the check kernel: zeros computed exactly, theorem cited",
                            bch={k: v for k, v in r.items() if k != "bound"})
    return bound, evidence


def radius(receipt, sector):
    """Radius within which the lightest pattern of a sector is unique, from a code
    receipt: (m - 1) // 2 with m the smaller of the proven lower bounds on the distance
    and on the lightest stabilizer of that type."""
    f = receipt["facts"]
    d, w = ("d_x", "wmin_x") if sector == "X" else ("d_z", "wmin_z")
    if f[d]["lower"] is None or f[w]["lower"] is None:
        raise certify.Refusal(f"the receipt proves no lower bound on {d} or {w}")
    requested = (min(f[d]["lower"], f[w]["lower"]) - 1) // 2
    from .codes import CSSCode
    code = CSSCode.from_json(receipt["code"])
    bound, _ = kernel_radius(code.hz if sector == "X" else code.hx, code.n, requested)
    if not 0 <= requested <= bound:
        raise certify.Refusal(f"{sector}: requested radius {requested} exceeds the proved kernel radius {bound}")
    return requested


def certificate(outer_receipt, inner_receipt):
    """What hierarchical decoding corrects, from the two code receipts."""
    from .codes import CSSCode
    certify.verify(outer_receipt)
    certify.verify(inner_receipt)
    inner = CSSCode.from_json(inner_receipt["code"])
    if certify.structure(inner)["k"] != 1:
        raise certify.Refusal("hierarchical decoding needs an inner code with k = 1")
    out = dict(schema="ftqc_codes-hier-1", outer=outer_receipt["code"]["name"],
               inner=inner_receipt["code"]["name"],
               outer_sha256=outer_receipt["code"]["sha256"], inner_sha256=inner_receipt["code"]["sha256"],
               n=outer_receipt["facts"]["n"] * inner_receipt["facts"]["n"], sectors={})
    for s in "XZ":
        out["sectors"][s] = dict(t1=radius(inner_receipt, s), t2=radius(outer_receipt, s))
    out["claim"] = ("in each sector, every error with at most t2 columns holding more than t1 "
                    "errors is returned up to a stabilizer")
    out["proof"] = __doc__.split("Proof.")[1].strip()
    return out
