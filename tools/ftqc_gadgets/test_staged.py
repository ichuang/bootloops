import unittest
from . import staged, staged_switch


class StagedTests(unittest.TestCase):
    def test_switches_and_phase(self):
        for op in ("raise", "rm15-phase", "lower"):
            receipt = staged.analyze(staged_switch.example(op))
            self.assertEqual(receipt["result"]["status"], "proved")
            staged.verify(receipt)

    def test_sign_and_record(self):
        for op, key, value in (("rm15-phase", "exponents", [1]*15),
                               ("lower", "corrections", [0]*8),
                               ("raise", "corrections", [0]*8)):
            spec = staged_switch.example(op)
            spec[key] = value
            self.assertEqual(staged.analyze(spec)["result"]["status"], "counterexample")

