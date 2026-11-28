import random
import unittest
from fractions import Fraction
from itertools import combinations

from .__main__ import evaluate
from .regions import lower_hull, pareto_minimal, regions_result

F = Fraction
SITES = ["d1", "d2", "d3", "m1", "m2", "m3"]
MINIMAL = [["d1", "d2"], ["d2", "d3"], ["d1", "m1"], ["m1", "m2", "m3"]]
KINDS = {s: ("data" if s[0] == "d" else "meas") for s in SITES}


def spec(**more):
    family = {"kind": "family", "expression": {"op": "family", "sites": SITES, "minimal": MINIMAL}}
    return dict(kind="regions", source=family, kinds=KINDS, **more)


class RegionTests(unittest.TestCase):
    def test_vectors_hull_and_ties(self):
        r = regions_result(spec())
        self.assertEqual(r["status"], "proved")
        self.assertEqual(r["kinds"], ["data", "meas"])
        self.assertEqual(r["count_vectors"], [{"counts": [0, 3], "supports": 1}, {"counts": [1, 1], "supports": 1},
                                              {"counts": [2, 0], "supports": 2}])
        self.assertEqual(r["lower_hull"], [[0, 3], [1, 1], [2, 0]])
        self.assertEqual(r["ties_at_exponent_ratio"]["values"], ["1/2", "1"])
        self.assertEqual(evaluate(spec()), r)

    def test_exponents_against_exact_probability(self):
        """With independent faults the bad-support probability is supports * eps^tau to leading order."""
        rows = [{"data": "1", "meas": m} for m in ("1/4", "1/2", "3/4", "1", "3")]
        r = regions_result(spec(exponents=rows))
        want = [("3/4", [[0, 3]], 1), ("3/2", [[0, 3], [1, 1]], 2), ("7/4", [[1, 1]], 1),
                ("2", [[1, 1], [2, 0]], 3), ("2", [[2, 0]], 2)]
        self.assertEqual([(x["tau"], x["dominant"], x["supports"]) for x in r["exponent_vectors"]], want)
        N = 10 ** 4                                           # eps = N^-4, so every rate is rational
        for row, x in zip(rows, r["exponent_vectors"]):
            rate = {"data": F(1, N ** 4), "meas": F(1, N ** int(4 * F(row["meas"])))}
            exact = F(0)
            for k in range(len(SITES) + 1):
                for S in combinations(SITES, k):
                    if any(set(t) <= set(S) for t in MINIMAL):
                        pr = F(1)
                        for s in SITES:
                            pr *= rate[KINDS[s]] if s in S else 1 - rate[KINDS[s]]
                        exact += pr
            lead = x["supports"] * F(1, N ** int(4 * F(x["tau"])))
            self.assertLess(abs(exact / lead - 1), F(1, 100))

    def test_hull_against_brute_force(self):
        rng = random.Random(5)
        for _ in range(200):
            pts = list({(rng.randrange(6), rng.randrange(6)) for _ in range(rng.randrange(1, 9))})
            hull, ties = lower_hull(pts)
            best = set()
            for a in range(1, 40):
                for b in range(1, 40):
                    vals = {p: a * p[0] + b * p[1] for p in pts}
                    m = min(vals.values())
                    winners = [p for p in pts if vals[p] == m]
                    if len(winners) == 1:
                        best.add(winners[0])
            self.assertTrue(best <= set(hull))
            self.assertTrue(set(hull) <= set(pareto_minimal(pts)))
            self.assertEqual(len(ties), len(hull) - 1)
            self.assertEqual(ties, sorted(ties))
            for (p, q), t in zip(zip(hull, hull[1:]), ties):
                self.assertEqual(p[0] + t * p[1], q[0] + t * q[1])

    def test_three_kinds_and_refusals(self):
        kinds = dict(KINDS, d3="idle")
        r = regions_result(dict(spec(), kinds=kinds, exponents=[{"data": "1", "meas": "1", "idle": "1/10"}]))
        self.assertNotIn("lower_hull", r)
        self.assertEqual(r["exponent_vectors"][0]["tau"], "11/10")
        self.assertEqual(r["exponent_vectors"][0]["dominant"], [[1, 1, 0]])
        self.assertIn([1, 1, 0], r["minimal_vectors"])
        with self.assertRaisesRegex(ValueError, "every non-ideal"):
            regions_result(dict(spec(), kinds={k: v for k, v in KINDS.items() if k != "m3"}))
        with self.assertRaisesRegex(ValueError, "positive"):
            regions_result(spec(exponents=[{"data": "0", "meas": "1"}]))
        with self.assertRaisesRegex(ValueError, "one exponent for every kind"):
            regions_result(spec(exponents=[{"data": "1"}]))
        with self.assertRaisesRegex(ValueError, "exact number"):
            regions_result(spec(exponents=[{"data": 0.5, "meas": "1"}]))


if __name__ == "__main__":
    unittest.main()
