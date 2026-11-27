"""Clifford SWAP between a raw wire and a CSS logical port in a live block."""

from ftqc_codes.gf2 import supp
from .surgery import pauli, validate_basis


def occupied_swap(block, slot):
    n, basis = block["n"], block["basis"]
    if type(slot) is not int or not 0 <= slot < len(basis):
        raise ValueError("invalid occupied block slot")
    names = [str(i) for i in range(len(basis))]
    checks = [pauli(w, n) for w in block["checks"]]
    validate_basis(
        dict(
            input=[
                dict(name=name, x=x, z=z)
                for name, (x, z) in zip(names, basis, strict=True)
            ],
            surviving=names,
        ),
        checks,
        [],
        n,
    )
    x, z = [pauli(w, n) for w in basis[slot]]
    if x.z or z.x or x.s or z.s:
        raise ValueError(
            "swap requires positive CSS representatives at the selected slot"
        )
    fan = [["CX", n, q] for q in supp(x.x)]
    ops = fan + [["CX", q, n] for q in supp(z.z)] + fan
    raw = 1 << n

    def image(a, b):
        for _, c, t in ops:
            a ^= ((a >> c) & 1) << t
            b ^= ((b >> t) & 1) << c
        return a, b

    # Real CNOTs map the ordered X^a Z^b monomial with positive sign.
    # Identity on the spectator Pauli algebra also fixes its reference.
    for g in checks:
        if image(g.x, g.z) != (g.x, g.z):
            raise ValueError("swap disturbed a stabilizer")
    for i, pair in enumerate(basis):
        for axis, w in enumerate(pair):
            g = pauli(w, n)
            expected = (
                ((raw, 0) if axis == 0 else (0, raw)) if i == slot else (g.x, g.z)
            )
            if image(g.x, g.z) != expected:
                raise ValueError("logical spectator or port disturbed by swap")
    if image(raw, 0) != (x.x, 0) or image(0, raw) != (0, z.z):
        raise ValueError("raw reference did not transfer")
    return dict(
        schema="ftqc-gadgets-occupied-swap-1",
        slot=slot,
        raw=n,
        discarded=[n],
        operations=ops,
        checks=block["checks"],
        basis=basis,
        unchanged_logicals=[i for i in range(len(basis)) if i != slot],
    )
