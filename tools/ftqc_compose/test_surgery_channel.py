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

    def test_decoded_bound(self):
        from fractions import Fraction
        from ftqc_gadgets.surgery_faults import columns
        from ftqc_gadgets.surgery import Instrument
        from ftqc_decoders.surgery_signatures import single_fault_table, decoded_damage
        from itertools import product

        s = branch(2, [], ["+ZZ"], 1)
        # Two noisy sites give an exhaustive independent check of the tail bound.
        for row in s["ticks"]:
            for op in row:
                op[0] = "!" + op[0]
        s["ticks"][:0] = [[["I", 0]], [["I", 1]]]
        cfg = dict(kind="decoded-surgery", receipt=receipt(s, 0), incoming_weight=0, p="1/10")
        result = evaluate(cfg)
        groups, _ = columns(Instrument(s))
        table = single_fault_table(groups)
        law = [[((0,0,0), Fraction(9,10))] +
               [(sig, Fraction(1,10*len(group))) for _,sig in group] for group in groups]
        fail = Fraction(0)
        for (a,p), (b,q) in product(*law):
            if any(decoded_damage(table, tuple(i^j for i,j in zip(a,b)))):
                fail += p*q
        self.assertLessEqual(fail, Fraction(result["metrics"]["joint"]["upper"]))
        with self.assertRaises(ValueError):
            evaluate(dict(cfg, incoming_weight=1))

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
