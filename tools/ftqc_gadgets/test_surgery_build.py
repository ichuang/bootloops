import unittest
from ftqc_gadgets.surgery import Instrument
from ftqc_gadgets.surgery_build import construct, graph


class ConstructorTests(unittest.TestCase):
    def test_cycle_and_mixed_parity(self):
        for mode in ("homomorphic", "cone"):
            s = graph(4, ["+XXXX"], "-XYZI", [(0, 1), (1, 2), (0, 2)], mode, 2)
            Instrument(s)
            self.assertEqual(s["targets"], ["-XYZI"])

    def test_leak(self):
        with self.assertRaises(ValueError):
            Instrument(construct(2, [], "ZZ", [1, 2], [], [(3, 0)], "cone"))

    def test_overlapping_joint(self):
        for mode in ("cone", "homomorphic"):
            s = construct(4, [], "ZZZZ", [3, 6, 12], [7], [(3, 0), (6, 1)], mode)
            # row 111 hides one parity, leaving two independent products.
            Instrument(s)
            self.assertEqual(s["targets"], ["+ZIZI", "-IZIZ"])
