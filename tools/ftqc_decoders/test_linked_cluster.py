"""Linked-cluster enclosures: additivity of localized decoders and containment of exact values."""

import unittest
from fractions import Fraction
from itertools import combinations

from ftqc_decoders.decimal_bounds import RefusalError
from ftqc_decoders.dem import repetition_problem, toric_x_problem
from ftqc_decoders.effective_distance import run_decoder
from ftqc_decoders.linked_cluster import (
    _mechanism_power_graph,
    _Signs,
    additivity_counterexample,
    connected_sets,
    decode_localized,
    enclose_linked_cluster,
    syndrome_groups,
)


def exact_flip_probability(problem, base, detector_range, p):
    """(1 - lambda) / 2 by summing over every fault set."""
    signs = _Signs(problem, base, detector_range, 1)
    M = problem.num_mechanisms
    lam = Fraction(0)
    for w in range(M + 1):
        shell = sum(signs.sign(F) for F in combinations(range(M), w))
        lam += shell * p ** w * (1 - p) ** (M - w)
    return (1 - lam) / 2


class LinkedCluster(unittest.TestCase):
    def test_connected_sets_are_enumerated_once(self):
        prob = toric_x_problem(2)
        adj = _mechanism_power_graph(prob, 1)
        got = list(connected_sets(adj, 3))
        self.assertEqual(len(got), len({frozenset(c) for c in got}))
        want = 0
        for w in (1, 2, 3):
            for F in combinations(range(prob.num_mechanisms), w):
                seen, todo = {F[0]}, [F[0]]
                while todo:
                    v = todo.pop()
                    for u in F:
                        if u not in seen and u in adj[v]:
                            seen.add(u)
                            todo.append(u)
                want += len(seen) == w
        self.assertEqual(len(got), want)

    def test_groups_partition_the_syndrome(self):
        prob = repetition_problem(9)
        s = (1 << 0) | (1 << 1) | (1 << 6)
        groups = syndrome_groups(prob, s, 1)
        self.assertEqual(sorted(groups), [0b11, 1 << 6])
        self.assertEqual(syndrome_groups(prob, s, 8), [s])
        self.assertEqual(syndrome_groups(prob, 0, 1), [])

    def test_localized_decoder_is_additive(self):
        for prob, base in ((repetition_problem(7), "union_find"), (toric_x_problem(3), "mwpm"),
                           (toric_x_problem(3), "union_find")):
            for r in (0, 1):
                signs = _Signs(prob, base, r, 1)
                self.assertIsNone(additivity_counterexample(prob, signs.sign, r + 1, 3))

    def test_a_nonlocal_rule_is_caught(self):
        prob = repetition_problem(9)
        signs = _Signs(prob, "union_find", 1, 1)

        def planted(F):                       # depends on the global number of fired detectors
            s = 0
            for i in F:
                s ^= prob.mech_detectors[i]
            return -signs.sign(F) if s.bit_count() == 4 else signs.sign(F)

        self.assertIsNotNone(additivity_counterexample(prob, planted, 2, 3))

    def test_dispatcher_runs_the_localized_decoder(self):
        prob = toric_x_problem(3)
        s = prob.mech_detectors[0]
        a = run_decoder(prob, s, {"algorithm": "localized", "base": "mwpm", "detector_range": 1})
        b = decode_localized(prob, s, "mwpm", 1)
        self.assertEqual((a.error_mask, a.observable_mask), (b.error_mask, b.observable_mask))
        self.assertEqual(prob.syndrome_of(a.error_mask), s)

    def test_enclosure_contains_the_exact_value(self):
        prob = toric_x_problem(3)
        for p in (Fraction(1, 2000), Fraction(1, 10 ** 7)):
            exact = exact_flip_probability(prob, "mwpm", 1, p)
            r = enclose_linked_cluster(prob, "mwpm", p, detector_range=1, max_cluster_size=2)
            self.assertTrue(r["subthreshold_condition_met"])
            self.assertLessEqual(Fraction(r["ler_lo"]), exact)
            self.assertLessEqual(exact, Fraction(r["ler_hi"]))
        self.assertTrue(r["lower_bound_nontrivial"])                 # at the smaller rate
        self.assertGreater(Fraction(r["ler_lo"]) * 3, Fraction(r["ler_hi"]))
        self.assertEqual(r["corrected_through_size"], 1)

    def test_enclosures_at_two_orders_agree(self):
        prob = repetition_problem(7)
        p = Fraction(1, 10 ** 5)
        r2 = enclose_linked_cluster(prob, "union_find", p, detector_range=1, max_cluster_size=2)
        r4 = enclose_linked_cluster(prob, "union_find", p, detector_range=1, max_cluster_size=4)
        self.assertEqual(r4["corrected_through_size"], r2["corrected_through_size"])
        # both contain the exact value, so the two intervals intersect
        self.assertLessEqual(Fraction(r4["ler_lo"]), Fraction(r2["ler_hi"]))
        self.assertLessEqual(Fraction(r2["ler_lo"]), Fraction(r4["ler_hi"]))
        exact = exact_flip_probability(prob, "union_find", 1, p)
        for r in (r2, r4):
            self.assertLessEqual(Fraction(r["ler_lo"]), exact)
            self.assertLessEqual(exact, Fraction(r["ler_hi"]))

    def test_refusals_and_the_radius(self):
        prob = toric_x_problem(3)
        r = enclose_linked_cluster(prob, "mwpm", "1/10", max_cluster_size=1)
        self.assertFalse(r["subthreshold_condition_met"])
        self.assertEqual((r["ler_lo"], r["ler_hi"]), ("0", "1"))
        with self.assertRaises(RefusalError):
            enclose_linked_cluster(prob, "mwpm", "1/1000", max_cluster_size=3, budget_cap=10)
        with self.assertRaises(RefusalError):
            enclose_linked_cluster(prob, "mwpm", "1/1000", kappa=0)


if __name__ == "__main__":
    unittest.main()
