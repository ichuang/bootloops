"""Leg F: the product bound for a concatenated CSS code, with the structure checked.

Setting.  An outer CSS code on N positions and an inner CSS code with k = 1 on
m qubits.  The concatenated code on N*m qubits (position i occupies qubits
i*m .. i*m+m-1) has as stabilizers the inner checks on every block and the
outer checks with each 1 replaced by an inner logical representative.

What is checked exactly: that the given code's X and Z stabilizer spaces are
exactly those of this construction (equality of row spaces, both types).
What is then concluded: d_x(code) >= d_x(inner) * d_x(outer), and likewise d_z.

Proof of the bound.  Let x be an X-type logical of the big code.  On each
block, x commutes with the inner Z checks, so (k_inner = 1) its restriction
is an inner X stabilizer plus c_i times the inner logical representative.
Blocks with c_i = 1 carry a nontrivial inner logical, of weight at least
d_x(inner).  Removing the inner stabilizer parts (they are stabilizers of the
big code) leaves the lift of o = (c_i).  The lift commutes with the lifted
outer Z checks, and a lifted Z check r overlaps it in |r meet o| mod 2
positions, so o is in ker H_Z(outer).  If o were in rowspace H_X(outer), the
lift would be a stabilizer and so would x.  Hence o is a nontrivial outer
logical, wt(o) >= d_x(outer), and wt(x) >= d_x(inner) * d_x(outer).

The outer and inner bounds are taken from their own receipts; the evidence
entry records the hashes of the codes used, and carries their matrices so the
structural check can be repeated by `certify.verify`.
"""
from . import certify, gf2
from .codes import CSSCode


def _logical_rep(checks, stabs, n):
    """A vector of ker(checks) outside rowspace(stabs)."""
    S = gf2.Echelon(stabs)
    for v in gf2.nullspace(checks, n):
        if v not in S:
            return v
    return None


def expected_stabilizers(outer, inner):
    m, N = inner.n, outer.n
    lx = _logical_rep(inner.hz, inner.hx, m)     # inner X logical
    lz = _logical_rep(inner.hx, inner.hz, m)     # inner Z logical
    if lx is None or lz is None:
        raise certify.Refusal("inner code has no logical operator")

    def build(inner_rows, outer_rows, rep):
        rows = [r << (i * m) for i in range(N) for r in inner_rows]
        for o in outer_rows:
            v = 0
            for i in gf2.supp(o):
                v |= rep << (i * m)
            rows.append(v)
        return rows
    return build(inner.hx, outer.hx, lx), build(inner.hz, outer.hz, lz)


def verify_structure(code, outer, inner):
    """Raise Refusal unless `code` is the concatenation of outer with inner."""
    if code.n != outer.n * inner.n:
        raise certify.Refusal("length is not outer length times inner length")
    si = certify.structure(inner)
    if si["k"] != 1:
        raise certify.Refusal("product bound needs an inner code with k = 1")
    certify.structure(outer)          # refuses non-commuting outer checks
    ex_x, ex_z = expected_stabilizers(outer, inner)
    for tag, have, want in (("X", code.hx, ex_x), ("Z", code.hz, ex_z)):
        a, b = gf2.Echelon(have), gf2.Echelon(want)
        if a.rank != b.rank or any(r not in b for r in a.basis):
            raise certify.Refusal(f"{tag} stabilizer space differs from the concatenation")


def run_leg(receipt, code, outer_receipt, inner_receipt):
    outer = CSSCode.from_json(outer_receipt["code"])
    inner = CSSCode.from_json(inner_receipt["code"])
    base = dict(leg="F", method="product bound for concatenated CSS codes; structure checked "
                "exactly, factor bounds from the outer and inner receipts",
                outer=outer_receipt["code"], inner=inner_receipt["code"])
    try:
        verify_structure(code, outer, inner)
    except certify.Refusal as e:
        for q in ("d_x", "d_z"):
            receipt["evidence"].append(dict(leg="F", quantity=q, status="skipped", reason=str(e)))
        return receipt
    for q in ("d_x", "d_z"):
        do, di = outer_receipt["facts"][q]["lower"], inner_receipt["facts"][q]["lower"]
        if do is None or di is None:
            receipt["evidence"].append(dict(leg="F", quantity=q, status="skipped",
                                            reason="a factor receipt has no proven bound"))
            continue
        ev = dict(base, quantity=q, status="proved", claim=f"{q} >= {do * di}",
                  outer_bound=do, inner_bound=di)
        certify.add_lower_bound(receipt, q, do * di, ev)
    return receipt


def recheck(code, q, ev):
    """Used by certify.verify: repeat the structural check and the arithmetic."""
    outer, inner = CSSCode.from_json(ev["outer"]), CSSCode.from_json(ev["inner"])
    verify_structure(code, outer, inner)
    if ev["outer_bound"] * ev["inner_bound"] != int(ev["claim"].split(">=")[1]):
        raise certify.Refusal(f"{q}: product bound arithmetic does not match its claim")
    return ev["outer_bound"], ev["inner_bound"]
