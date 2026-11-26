import unittest
from itertools import combinations, product

from .supports import Family, evaluate, family_result, predicate


def majority(sites):
    return dict(op="family", sites=sites, minimal=list(map(list, combinations(sites, 2))))


class Supports(unittest.TestCase):
    def test_majorities(self):
        outer = majority(list("abc"))
        self.assertEqual(family_result(dict(kind="family", expression=outer))["enumerator"], {"2": 3})
        node = dict(op="substitute", outer=outer,
                    at={v: majority([v + str(i) for i in range(3)]) for v in "abc"})
        result = family_result(dict(kind="family", expression=node))
        self.assertEqual((result["enumerator"], result["checked_subsets"]), ({"4": 27}, 512))
        for values in product((0, 1), repeat=9):
            bad = {s for s, bit in zip(result["sites"], values) if bit}
            self.assertEqual(predicate(node, bad), sum(sum(values[j:j+3]) >= 2 for j in (0,3,6)) >= 2)
        node["at"]["b"] = node["at"]["a"]
        with self.assertRaisesRegex(ValueError, "shared physical"):
            evaluate(node)

    def test_family_order_and_units(self):
        a = Family("xy", ["x", "xy"])
        b = Family("xy", ["xy"])
        self.assertTrue(a.covers(b))
        self.assertFalse(b.covers(a))
        self.assertEqual(a.union(a).describe(), a.describe())
        self.assertTrue(Family([], [[]]).bad(set()))
        self.assertFalse(Family([], []).bad(set()))
        with self.assertRaises(ValueError):
            a.product(b)
        with self.assertRaises(ValueError):
            Family("x", ["y"])
