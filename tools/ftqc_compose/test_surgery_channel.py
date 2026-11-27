import unittest
from .__main__ import evaluate
from .contracts import derive
from ftqc_gadgets.surgery_parallel import branch
from ftqc_gadgets.surgery import receipt
from ftqc_gadgets.test_surgery import edge


class SurgeryChannelTests(unittest.TestCase):
    def test_conditioned_law(self):
        spec = branch(3, [], ["+ZZI", "+ZIZ"], 1)
        # A single faulty shared input gate can flip both returned parities.
        spec["ticks"] = [
            [["!" + op[0], *op[1:]] for op in row] for row in spec["ticks"]
        ]
        spec["ticks"].insert(0, [["I", 0]])
        cfg = dict(
            kind="surgery-channel",
            receipt=receipt(spec, 0),
            incoming_weight=1,
            p="3/10",
            cap=1000,
        )
        out = evaluate(cfg)
        self.assertEqual(out["status"], "proved")
        from fractions import Fraction

        joint = sum(
            Fraction(v["p"])
            for v in out["clean_input"]["conditional"]
            if v["readout"] == 3
        )
        self.assertEqual(joint, Fraction(1, 5))
        self.assertEqual(out["envelope"][-1]["accepted_damage_upper"], "1")
        self.assertEqual(evaluate(dict(cfg, cap=1))["status"], "inconclusive")

    def test_open_input_leaf(self):
        cfg = dict(
            op="surgery",
            instance="signed",
            receipt=receipt(edge(), 0),
            wires=list("abcdefgh"),
            incoming=0,
        )
        self.assertEqual(
            derive(cfg)["evidence"]["joint_input_internal"]["status"], "proved"
        )
