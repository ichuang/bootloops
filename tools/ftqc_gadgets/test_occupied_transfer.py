import unittest

from ftqc_codes.gf2 import Echelon, dot, nullspace
from ftqc_extractors.algebra import hgp
from ftqc_extractors.block_map import word
from ftqc_gadgets.occupied_transfer import occupied_swap


def hgp450_block():
    # Cyclic [15,4] seed with offsets 0,3,4, as in the HGP450 T fixture.
    rows = [sum(1 << ((q + d) % 15) for d in (0, 3, 4)) for q in range(15)]
    n, hx, hz = hgp(rows, 15, rows, 15)

    def quotient(checks, dual):
        span = Echelon(checks)
        return [v for v in nullspace(dual, n) if span.add(v)]

    xs, zs = quotient(hx, hz), quotient(hz, hx)
    basis = []
    for x in xs:
        j = next(j for j, z in enumerate(zs) if dot(x, z))
        z = zs.pop(j)
        zs = [u ^ (z if dot(x, u) else 0) for u in zs]
        for i, (old, other) in enumerate(basis):
            if dot(old, z):
                basis[i] = (old ^ x, other)
        basis.append((x, z))
    return dict(
        n=n,
        checks=[word(n, x, 0) for x in hx] + [word(n, 0, z) for z in hz],
        basis=[[word(n, x, 0), word(n, 0, z)] for x, z in basis],
    )


class ScaleTransfer(unittest.TestCase):
    def test_named_hgp450_spectators(self):
        block = hgp450_block()
        proof = occupied_swap(block, 0)
        self.assertEqual(block["n"], 450)
        self.assertEqual(len(proof["unchanged_logicals"]), 31)
