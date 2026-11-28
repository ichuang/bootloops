"""Local-stochastic filter studies, moment contracts in the tail bounds, and the closure study."""
import copy
import json
import unittest
from fractions import Fraction
from math import comb
from pathlib import Path

from . import filters as fl
from . import filters_ls as ls
from .__main__ import receipt
from .closure import closure_result
from .tails import cone_bound_result, event_bounds, ls_tail, up
from .test_filters import RAW, encoded

F = Fraction
EXAMPLES = Path(__file__).parent / "examples"
INNER = {k: {s: {"junk": "1/1000000000", "modes": [["1", "1/20000"]]} for s in "XZ"} for k in ("plus", "zero")}
GOLAY = dict(kind="filter-contract", level="physical", noise="local-stochastic", p="1/10000",
             tree={"n": 23, "wmin": {"X": 8, "Z": 7}, "sequence": "XZXZ"},
             leaf={"locations": 151, "X": {"junk": 70, "local": 5}, "Z": {"junk": 71, "local": 5}})


class MomentContracts(unittest.TestCase):
    def test_ls_tail(self):
        self.assertEqual(ls_tail(7, 1, F(1, 100), []), 21 * F(1, 100) ** 2)
        # one process with a single mode of weight one adds its rate
        self.assertEqual(ls_tail(7, 1, F(1, 100), [[(F(1), F(1, 50))]]), 21 * F(3, 100) ** 2)
        # weight below one: the process is absent with the rest of the weight
        got = ls_tail(7, 1, F(1, 100), [[(F(1, 4), F(1, 50))]])
        self.assertEqual(got, 21 * (F(1, 4) * F(3, 100) ** 2 + F(3, 4) * F(1, 100) ** 2))
        # against the sum over the ways of splitting two dirty coordinates between sources
        b, w, q = F(1, 100), F(1, 4), F(1, 50)
        self.assertEqual(got, 21 * (b * b + 2 * b * w * q + w * q * q))
        self.assertEqual(ls_tail(3, 3, F(1, 2), []), 0)
        self.assertEqual(ls_tail(7, 1, F(1), []), 1)

    def test_equivalents_are_a_special_case(self):
        prog = {"name": "small", "inputs": ["A"], "outputs": ["A"], "layers": [[["ECa", "A"]], [["ECb", "A"]]]}
        base = dict(kind="cone-bound", noise="local-stochastic", p="1/1000", inner={"n": 7, "t": 1}, program=prog)
        eq = {k: {s: {"equivalents": 3} for s in "XZ"} for k in ("plus", "zero")}
        md = {k: {s: {"junk": "0", "modes": [["1", "3/1000"]]} for s in "XZ"} for k in ("plus", "zero")}
        self.assertEqual(event_bounds(dict(base, contracts=eq)), event_bounds(dict(base, contracts=md)))
        # junk of the two ancillas of a correction is charged to its decoding event
        jk = {k: {s: {"junk": "1/100000", "modes": [["1", "3/1000"]]} for s in "XZ"} for k in ("plus", "zero")}
        with_junk, _ = event_bounds(dict(base, contracts=jk))
        without, _ = event_bounds(dict(base, contracts=md))
        ec = [a - b for a, b in zip(with_junk, without) if a != b]
        self.assertTrue(ec and all(d == F(2, 100000) for d in ec))


class FilterStudies(unittest.TestCase):
    def test_physical_tree(self):
        r = fl.filter_contract_result(GOLAY)
        self.assertEqual(r["status"], "conditional")
        self.assertEqual(r["noise"], "local-stochastic")
        self.assertTrue(any("p^|S|" in x for x in r["premises"]))
        self.assertFalse(any("independent" in x for x in r["premises"]))
        p = F(1, 10000)
        nu = 70 * p + 23 * 5 * p
        self.assertEqual(r["accept_lower"][0], fl.decimal(fl.down(fl.down(1 - 2 * nu - 46 * p), 20)))
        # close to the calculus for independent faults, and never below it
        ind = fl.filter_contract_result({k: v for k, v in GOLAY.items() if k != "noise"})
        for s in "XZ":
            a, b = F(r["contract"][s]["junk"]), F(ind["contract"][s]["pi"])
            self.assertGreaterEqual(a, b)
            self.assertLess(a, b * F(11, 10))
        self.assertEqual(receipt(GOLAY), receipt(copy.deepcopy(GOLAY)))
        self.assertEqual(fl.filter_contract_result(dict(GOLAY, p="1/100"))["status"], "no bound")
        with self.assertRaisesRegex(ValueError, "local-stochastic"):
            fl.filter_contract_result(dict(GOLAY, noise="gaussian"))

    def test_first_check_written_out(self):
        """One check of two leaves, every weight from the rule."""
        p, n, wmin = F(1, 1000), 7, 4
        lf = {s: ls.leaf(p, 6, 2) for s in "XZ"}
        step = dict(fresh=dict(X=2 * p, Z=2 * p), partner=2 * p, clean=1 - 14 * p)
        c, alpha = ls.check(lf, lf, "X", n, wmin, step)
        pi, q, y = 6 * p, 2 * p, 2 * p + 2 * p
        rho = 2 * (pi + n * q) + 14 * p
        self.assertEqual(alpha, fl.down(1 - rho))
        self.assertEqual(c["X"]["pi"], up((pi * pi + comb(n, wmin) * (q + y) ** wmin) / alpha))
        for k in (1, 2, 3):
            one_attempt = (2 * p + q * y) ** k + pi * ((2 * p + q) ** k - (2 * p) ** k) + pi * ((2 * p + y) ** k - (2 * p) ** k)
            self.assertGreaterEqual(ls.moment(ls.total(c["X"]), k), one_attempt / alpha)
            # normalizing the weights to one moves the retry factor into the rates: at most alpha^-k
            self.assertLess(ls.moment(ls.total(c["X"]), k), one_attempt / alpha ** k * F(101, 100))
        self.assertEqual(c["Z"]["pi"], up(2 * pi / alpha))

    def test_encoded_tree(self):
        r = fl.filter_contract_result(encoded(noise="local-stochastic", contracts=INNER))
        self.assertEqual(r["status"], "conditional")
        self.assertTrue(r["carry_last"])
        ind = fl.filter_contract_result(encoded())
        for s in "XZ":
            a, b = F(r["contract"][s]["junk"]), F(ind["contract"][s]["pi"])
            self.assertGreaterEqual(a, b)
            self.assertLess(a, b * F(11, 10))
            self.assertLess(a, F(r["raw"]["junk"][s]))                # purification helps
        with self.assertRaises(ValueError):                           # contracts of the other noise model
            fl.filter_contract_result(encoded(noise="local-stochastic"))


class Closure(unittest.TestCase):
    def setUp(self):
        self.spec = json.loads((EXAMPLES / "closure-two-level-steane.json").read_text())

    def test_banked_example(self):
        r = closure_result(self.spec)
        self.assertEqual(r["status"], "conditional")
        self.assertEqual(r["upper_sci"], "5.121e-13")
        self.assertEqual([s["ancilla"] for s in r["stages"]], ["plus", "zero", "plus2", "zero2"])
        self.assertEqual(r["cones"]["outer_decoding_events"], 4)
        self.assertIn("p^|S|", r["noise_premise"])
        self.assertFalse(any("contract" in x and "declared" in x for x in r["deterministic_premises"]))
        self.assertEqual(receipt(self.spec)["result"], r)

    def test_the_chain_is_the_three_studies(self):
        spec = self.spec
        common = dict(kind="filter-contract", noise="local-stochastic", p=spec["p"])
        contracts = {k: fl.filter_contract_result(dict(common, level="physical", **v))["contract"]
                     for k, v in spec["inner_ancillas"].items()}
        inner = dict(contracts)
        for k, v in spec["outer_ancillas"].items():
            contracts[k] = fl.filter_contract_result(dict(common, level="encoded", inner=spec["inner"],
                                                          contracts=inner, **v))["contract"]
        direct = cone_bound_result(dict(kind="cone-bound", noise="local-stochastic", p=spec["p"], inner=spec["inner"],
                                        outer=spec["outer"], contracts=contracts, program=spec["program"]))
        self.assertEqual(closure_result(spec)["upper"], direct["upper"])

    def test_scaling_and_monotonicity(self):
        at = lambda p: F(closure_result(dict(self.spec, p=p))["upper"])
        a, b = at("1/1000000"), at("1/10000000")
        self.assertGreater(a, b)
        self.assertLess(abs(a / b / 10 ** 4 - 1), F(1, 100))          # two levels of t = 1: fourth order
        worse = copy.deepcopy(self.spec)
        worse["inner_ancillas"]["zero"]["leaf"]["X"]["junk"] = 60
        self.assertGreater(F(closure_result(worse)["upper"]), a)

    def test_refusals(self):
        with self.assertRaisesRegex(ValueError, "no bound|rejection"):
            closure_result(dict(self.spec, p="1/100"))
        missing = copy.deepcopy(self.spec)
        del missing["outer_ancillas"]["plus2"]
        with self.assertRaisesRegex(ValueError, "no contract"):
            closure_result(missing)
        both = copy.deepcopy(self.spec)
        both["outer_ancillas"]["zero"] = both["outer_ancillas"]["zero2"]
        with self.assertRaisesRegex(ValueError, "inner or outer"):
            closure_result(both)
        with self.assertRaises(ValueError):
            closure_result(dict(self.spec, extra=1))


if __name__ == "__main__":
    unittest.main()
