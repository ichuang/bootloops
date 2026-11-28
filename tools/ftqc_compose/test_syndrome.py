import unittest
from fractions import Fraction
from itertools import product

from ftqc_compose.__main__ import receipt
from ftqc_compose.syndrome import syndrome_envelope


class SyndromeEnvelopeTests(unittest.TestCase):
    def test_one_effect_per_location(self):
        spec = dict(kind="syndrome-envelope", bits=2, p="1/10",
                    classes=[dict(count=1, effects=[0, 1, 2, 2])], event=[3])
        self.assertEqual(receipt(spec)["result"]["syndrome_event_upper"], "0")
        spec["classes"][0]["count"] = 2
        spec["event"] = [3, 3]
        result = receipt(spec)["result"]
        self.assertEqual(result["status"], "conditional")
        self.assertEqual(result["syndrome_event_upper"], "1/50")

    def test_exhaustive_supports_with_adversarial_effect_choice(self):
        modes = [(1, 2), (2,), (1, 3)]
        classes = [dict(count=1, effects=list(group)) for group in modes]
        p = Fraction(1, 3)
        envelope = syndrome_envelope(2, classes, str(p))
        for target in range(4):
            worst = Fraction(0)
            for faulty in product((0, 1), repeat=len(modes)):
                options = [((0,) + group if f else (0,)) for group, f in zip(modes, faulty)]
                if any(a ^ b ^ c == target for a, b, c in product(*options)):
                    k = sum(faulty)
                    worst += p**k * (1-p)**(len(modes)-k)
            self.assertGreaterEqual(envelope[target], worst)
        only = syndrome_envelope(1, [dict(count=1, effects=[1])], "1/3")
        self.assertLessEqual(p, only[1])
        self.assertLess(only[1] - p, Fraction(1, 10**48))

    def test_empty_and_invalid(self):
        self.assertEqual(syndrome_envelope(0, [], "1"), {0: 1})
        for count, effects, p in [(-1, [1], "0"), (True, [1], "1"),
                                  (0, [2], "0"), (1, [False], "1"), (1, [], 0.1)]:
            with self.subTest(count=count, effects=effects, p=p), self.assertRaises(ValueError):
                syndrome_envelope(1, [dict(count=count, effects=effects)], p)
