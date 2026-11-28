import copy
import unittest
from importlib.util import find_spec
from fractions import Fraction
from .__main__ import receipt, verify


class ExportTests(unittest.TestCase):
    def test_joint_and_replay(self):
        if find_spec("stim") is None:
            self.skipTest("optional Stim")
        m = dict(
            type="closed-Y",
            case="bb-uniform",
            locations=[[25, [38, 254]], [59, [144]], [19, [64]]],
            frame=[0] * 4,
        )
        s = dict(schema=1, kind="injection_export", attempt="a", block="BB", model=m)
        r = receipt(s)
        verify(s, r)
        law = r["result"]
        both = sum(
            Fraction(v) for i, v in law["accepted_mass"].items() if int(i) & 12 == 12
        )
        self.assertEqual(both / Fraction(law["acceptance"]), Fraction(1, 14988))
        for key in ("frame", "slots", "preparation_sha256"):
            bad = copy.deepcopy(r)
            bad["result"][key] = None
            with self.assertRaises(ValueError):
                verify(s, bad)
        self.assertEqual(receipt(s, budget=1)["result"]["status"], "inconclusive")
