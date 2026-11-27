import json
import unittest
from pathlib import Path
from .__main__ import receipt, verify
from .algebra import parity

EXAMPLES = Path(__file__).parent / "examples"


def example(name):
    return json.loads((EXAMPLES / (name + ".json")).read_text())


class SourceTests(unittest.TestCase):
    def test_signed_general_and_syndrome_frames(self):
        s = example("general-injection")
        result = receipt(s)
        verify(s, result)
        r = result["result"]
        self.assertEqual(r["status"], "proved")
        self.assertTrue(r["carrier_gates"])
        for axis, rows in r["syndrome_frame"].items():
            checks = s["memory"]["hx" if axis == "Z" else "hz"]
            for idx, c in rows:
                self.assertEqual(parity(checks[idx], c), 1)
                self.assertFalse(
                    any(
                        parity(p[0 if axis == "Z" else 1], c)
                        for p in r["representatives"]
                    )
                )

    def test_overlapping_logicals(self):
        for name, status in (
            ("bell-injection", "proved"),
            ("split-bell-injection", "counterexample"),
        ):
            s = example(name)
            proof = receipt(s)
            verify(s, proof)
            self.assertEqual(proof["result"]["status"], status)
