import unittest

from .stochastic import stochastic_result
from .test_level import example


class StochasticTests(unittest.TestCase):
    def test_correlated_law_attains_bound(self):
        expr = {"op": "family", "sites": ["a", "b"], "minimal": [["a", "b"]]}
        s = {
            "kind": "local-stochastic",
            "p": "1/4",
            "source": {"kind": "family", "expression": expr},
            "law": [
                {"support": ["a", "b"], "mass": "1/16"},
                {"support": [], "mass": "15/16"},
            ],
        }
        r = stochastic_result(s)
        self.assertEqual(
            (r["upper"], r["noise_evidence"]["bad_event_mass"]), ("1/16", "1/16")
        )
        self.assertIsNone(r["also_bounds"])
        s["law"][0]["mass"] = "1/4"
        s["law"][1]["mass"] = "3/4"
        with self.assertRaisesRegex(ValueError, "subset condition"):
            stochastic_result(s)

    def test_bound_scope_and_clipping(self):
        s = {"kind": "local-stochastic", "p": "0.01", "source": example()}
        self.assertEqual(stochastic_result(s)["upper"], "27/100000000")
        self.assertIsNotNone(stochastic_result(s)["also_bounds"])
        s["source"]["output_check"] = "none"
        self.assertIsNone(stochastic_result(s)["also_bounds"])
        r = stochastic_result(dict(s, p=1))
        self.assertEqual((r["raw_bound"], r["upper"]), ("27", "1"))
        for p in [-1, True, 0.1, "1/0", "1.01"]:
            with self.assertRaises(ValueError):
                stochastic_result(dict(s, p=p))

    def test_ideal_sites(self):
        with self.assertRaisesRegex(ValueError, "ideal physical site"):
            stochastic_result(
                {
                    "kind": "local-stochastic",
                    "p": 1,
                    "source": example(),
                    "law": [{"support": ["vote"], "mass": 1}],
                }
            )

    def test_zero_rate_and_empty_generator(self):
        s = {
            "kind": "local-stochastic",
            "p": 0,
            "source": {
                "kind": "family",
                "expression": {"op": "family", "sites": [], "minimal": [[]]},
            },
        }
        self.assertEqual(stochastic_result(s)["upper"], "1")
        s["source"]["expression"]["minimal"] = []
        self.assertEqual(stochastic_result(s)["upper"], "0")
