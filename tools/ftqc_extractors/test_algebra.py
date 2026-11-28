import unittest
from .algebra import *


class HGPTests(unittest.TestCase):
    def test_reference(self):
        for a in [["11"], ["11", "11"], ["101", "011"], ["0"], ["1"], ["001", "001"]]:
            for b in [["11"], ["0", "0"], ["1", "1"]]:
                rows, n = seed(a)
                cols, m = seed(b)
                code = hgp(rows, n, cols, m)
                info = hgp_parameters(rows, n, cols, m)
                witness = distance_below(code, code[0] + 1)
                self.assertEqual(info["d"], len(witness["qubits"]) if witness else None)

    def test_paper(self):
        h = cyclic(15, 25)
        self.assertEqual(hgp_parameters(h, 15, h, 15), dict(n=450, k=32, d=8))
