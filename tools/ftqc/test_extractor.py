import unittest
from .algebra import *
from .extractor import *


class SurgeryTests(unittest.TestCase):
    def test_column(self):
        spec = dict(
            a=["110", "011"],
            b=["11", "11"],
            column=0,
            word=7,
            vertices=3,
            ports=[0, 1, 2],
            edges=[[0, 1], [1, 2]],
        )
        merged, small = column(spec)
        self.assertIsNone(distance_below(merged, 2))
        self.assertIsNone(distance_below(small, 2))
        with self.assertRaises(ValueError):
            column(dict(spec, paths=[0, 0]))

    def test_nontrivial_reduced_sector(self):
        spec = dict(
            a=["111"],
            b=["11"],
            column=0,
            word=3,
            vertices=3,
            ports=[0, 1, 2],
            edges=[[0, 1], [1, 2]],
        )
        for n, x, z in column(spec):
            distance = next(
                w
                for w in range(1, n + 1)
                if any(
                    not any(parity(v, s) for s in x) and reduce_row(v, echelon(z))
                    for v in range(1 << n)
                    if v.bit_count() == w
                )
            )
            self.assertEqual(distance, 2)

    def test_cyclic(self):
        s = dict(n=3, polynomial=7, edges=[[0, 1], [1, 2], [2, 0]])
        self.assertEqual(canonical(s)["lower"], 2)
        with self.assertRaises(ValueError):
            canonical(dict(s, assembly=[]))
