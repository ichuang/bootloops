import unittest
from fractions import Fraction
from .injection_batch import result


class JointTests(unittest.TestCase):
    def test_correlated_consumers_and_refusals(self):
        try:
            import stim
        except ImportError:
            self.skipTest("optional Stim")
        spec = dict(
            kind="injection-batch",
            case="bb-uniform",
            locations=[[25, [38, 254]], [59, [144]], [19, [64]]],
            frame=[0] * 4,
            slots=[2, 3],
            state_budget=500,
        )
        r = result(spec)
        p = list(map(Fraction, r["consumer_Z_masks"]))
        self.assertEqual(sum(p), 1)
        self.assertGreater(p[3], (p[1] + p[3]) * (p[2] + p[3]))
        self.assertLess(Fraction(r["acceptance"]), 1)
        for changes in (
            dict(state_budget=1),
            dict(locations=[[10000, [1]]]),
            dict(slots=[2, 2]),
            dict(frame=[1 << 500] * 4),
        ):
            with self.assertRaises(ValueError):
                result(dict(spec, **changes))
