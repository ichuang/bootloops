import json
import unittest
from pathlib import Path

from faultcert.circuit import Circuit
from faultcert.model import Model

from .faultcert_bridge import Bridge, cat_prefix

EXAMPLES = Path(__file__).parent / "examples"


def example(name):
    return json.loads((EXAMPLES / (name + ".json")).read_text())


class CrossCheckTests(unittest.TestCase):
    def test_shared_predicate(self):
        cases = [("edge-extractor", "gate", 2), ("edge-hook", "gate", 1)]
        cases += [
            (name, arm, distance)
            for name, distance in (("t-gate", 2), ("t-hook", 1))
            for arm in ("gate", "X", "Y")
        ]
        for name, arm, expected in cases:
            with self.subTest(name=name, arm=arm):
                bridge = Bridge(example(name), arm)
                result = bridge.compare(2, 100000)
                self.assertEqual(
                    result["checked_modes"], sum(map(len, bridge.ft.modes))
                )
                self.assertEqual(result["distance"]["lower"], expected)
                self.assertEqual(result["distance"]["upper"], expected)
                self.assertTrue(result["distance"]["witness"]["undetected_harm"])

    def test_search_limits_and_mismatch(self):
        bridge = Bridge(example("edge-extractor"))
        result = bridge.compare(2, 1)["distance"]
        self.assertFalse(result["exact"])
        self.assertIsNone(result["upper"])
        with self.assertRaisesRegex(ValueError, "once"):
            bridge.replay([[0, 1], [0, 1]])
        with self.assertRaisesRegex(ValueError, "location or mode"):
            bridge.replay([[-1, 1]])
        with self.assertRaisesRegex(ValueError, "location or mode"):
            bridge.replay([[0, 7]])
        bridge.tests = []
        with self.assertRaisesRegex(ValueError, "disagree"):
            bridge.compare(2, 100000)

    def test_later_cat_round_and_boundaries(self):
        spec = example("t-gate")
        bridge = Bridge(spec)
        cat = cat_prefix(spec, "gate", 1)
        native = Circuit.from_json(cat)
        origin = cat["ftqc_origin"]
        start = origin["first_location"]
        self.assertEqual(start, bridge.ft.parts[1][0])
        self.assertGreater(start, 0)
        self.assertEqual(
            len(native.locations()), bridge.ft.parts[1][-1][-1] - start + 1
        )
        model = Model(native)
        self.assertEqual(len(model.detectors), 3)
        detector, _, residual = model.replay([(0, "Z")])
        self.assertEqual(detector, 0)
        self.assertNotEqual(residual, 0)
        with self.assertRaises(ValueError):
            cat_prefix(spec, "gate", -1)
        with self.assertRaises(ValueError):
            Bridge(dict(spec, resource="raw-T"))
        with self.assertRaises(ValueError):
            Bridge(example("edge-extractor"), "Y")
