import copy
import unittest
from fractions import Fraction

from . import filters as fl
from . import volumes
from .__main__ import receipt
from .test_filters import INNER, RAW, encoded

F = Fraction
TWO = {"name": "two", "inputs": ["D"], "outputs": ["D", "Q"],
       "layers": [[["NEW", "zero2", "Q"], ["NEW", "and2", "a", "b", "t"]], [["CX", "D", "Q"], ["MZ", "a", "m"], ["MZ", "b", "n"]],
                  [["ECa", "D"], ["ECa", "Q"], ["MZ", "t", "k"]], [["ECb", "D"], ["ECb", "Q"]], [["EC2a", "Q"]], [["EC2b", "Q"]]]}


def spec(**change):
    s = dict(kind="volume", inner={"n": 23}, outer={"n": 47}, program=copy.deepcopy(TWO),
             costs={"plus": "100", "zero": "200", "plus2": "5000", "zero2": "7000.5"}, uncosted=["and2"])
    s.update(change)
    return s


class Volumes(unittest.TestCase):
    def test_every_part_by_hand(self):
        r = volumes.volume_result(spec())
        # live at the end of each layer: D Q a b t | D Q t | D Q | D Q | D Q ~p | D Q
        ticks = 5 + 3 + 2 + 2 + 3 + 2
        self.assertEqual(r["counts"], dict(layers=6, block_ticks=ticks, peak_blocks=5, inner_corrections=2, outer_corrections=1,
                                           resources={"and2": 1, "plus2": 1, "zero2": 2}))
        self.assertEqual(r["parts"]["online"], str(23 * 47 * ticks))
        self.assertEqual(r["parts"]["inner_corrections"], str(2 * 47 * 300))
        self.assertEqual(r["parts"]["resources"], {"plus2": "5000", "zero2": "14001"})
        self.assertEqual(F(r["total"]), 23 * 47 * ticks + 2 * 47 * 300 + 5000 + 14001)
        self.assertEqual(r["uncosted"], {"and2": 1})
        self.assertEqual(r["peak_physical_qubits"], 5 * 23 * 47)
        self.assertEqual(r["status"], "computed")
        rec = receipt(spec())
        self.assertEqual(rec["result"], r)
        import json
        from pathlib import Path

        banked = json.loads((Path(__file__).parent / "examples" / "volume-two-level.json").read_text())
        self.assertEqual(banked, spec())

    def test_single_level_and_inner_resources(self):
        prog = {"name": "one", "inputs": [], "outputs": ["A", "C"],
                "layers": [[["NEW", "zero", "A", "B"], ["NEW", "std", "C"]], [["CX", "A", "B"]], [["MZ", "B"], ["ECa", "A"]],
                           [["ECb", "A"]]]}
        r = volumes.volume_result(dict(kind="volume", inner={"n": 7}, program=prog, costs={"plus": "10", "zero": "20"}))
        self.assertEqual(F(r["total"]), 7 * (3 + 3 + 2 + 2) + 30 + 2 * 20)       # std blocks are free
        self.assertEqual(r["uncosted"], {})

    def test_refusals(self):
        for bad in (spec(costs={"plus": "100", "zero": "200", "plus2": "5000"}),          # zero2 has no cost
                    spec(uncosted=[]),                                                  # and2 neither priced nor set aside
                    spec(uncosted=["and2", "zero2"]),                                   # both priced and set aside
                    spec(costs={"plus": "-1", "zero": "200", "plus2": "1", "zero2": "1"}),
                    spec(costs={"plus": "x", "zero": "200", "plus2": "1", "zero2": "1"}),
                    {k: v for k, v in spec().items() if k != "outer"},                  # two levels need the outer n
                    spec(inner={"n": 0}), spec(extra=1)):
            with self.assertRaises(ValueError):
                volumes.volume_result(bad)

    def test_encoded_tree_priced_from_inner_costs(self):
        costs = {"plus": "1000", "zero": "1200"}
        r = fl.filter_contract_result(encoded(inner_costs=costs))
        self.assertEqual(r["cost_unit"], "physical locations")
        raw, _, raw_total, _ = volumes.program_volume(RAW, 7, None, {k: F(v) for k, v in costs.items()})
        # the raw preparation: 7 blocks over 8 layers, 14 corrections, 3 |+> and 4 |0> blocks delivered
        self.assertEqual(raw["block_ticks"], 7 * 8)
        self.assertEqual(raw_total, 7 * 56 + 14 * 2200 + 3 * 1000 + 4 * 1200)
        self.assertEqual(F(r["priced"]["raw"]["locations"]), raw_total)
        step = 7 * 4 + 2200                                                            # A B | A | A, one correction
        self.assertEqual(F(r["priced"]["X"]["locations_per_column"]), step)
        acc = [F(a) for a in r["accept_lower"]]
        lower = (2 * (2 * raw_total + 7 * step) / 1 + 7 * step) / 1                    # with acceptance 1
        self.assertGreater(F(r["expected_cost"]), lower)
        self.assertLess(F(r["expected_cost"]), lower / (acc[0] * acc[0] * acc[1]) * F(1001, 1000))
        same = fl.filter_contract_result(encoded(costs={"raw": str(raw_total), "X": str(7 * step), "Z": str(7 * step)}))
        self.assertEqual(F(same["expected_cost"]), F(r["expected_cost"]))
        self.assertEqual(same["contract"], r["contract"])
        with self.assertRaises(ValueError):
            fl.filter_contract_result(encoded(inner_costs=costs, costs={"raw": "1", "X": "1", "Z": "1"}))
        with self.assertRaises(ValueError):
            fl.filter_contract_result(encoded(inner_costs={"plus": "1"}))


if __name__ == "__main__":
    unittest.main()
