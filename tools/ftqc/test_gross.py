import unittest
from functools import reduce
from operator import xor
from .algebra import bits, css_dimension, echelon, reduce_row
from .gross import cell, gross_z, logical_action, shifted_z


class GrossConstructionTests(unittest.TestCase):
    def test_tanner_coordinates_and_check_six(self):
        memory, (n, x, z), ports, incident = gross_z()
        for i in range(12):
            for j in range(6):
                row = sum(1 << cell(i + a, j + b) for a, b in [(3, 0), (0, 1), (0, 2)])
                row |= sum(
                    1 << (72 + cell(i + a, j + b)) for a, b in [(0, 3), (1, 0), (2, 0)]
                )
                self.assertEqual(memory.hx[cell(i, j)], row)
        # Figure 9 labels on the four edges of the purple square numbered 6.
        self.assertEqual({incident[q - 144] for q in bits(x[-1])}, {4, 10, 19, 20})
        self.assertEqual((n, len(ports)), (162, 12))
        self.assertEqual(reduce(xor, ports), shifted_z())
        self.assertEqual(css_dimension((n, x[:-1], z)), 12)
        self.assertEqual(css_dimension((n, x, z)), 11)

    def test_surviving_logicals(self):
        action = logical_action()
        self.assertEqual((action["memory_k"], action["merged_k"]), (12, 11))
        _, (_, _, z), _, _ = gross_z()
        self.assertTrue(reduce_row(shifted_z(0, 1), echelon(z)))
