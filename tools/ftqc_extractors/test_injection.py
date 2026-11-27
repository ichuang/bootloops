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

    def test_reference_sources_and_postselection(self):
        for name, k, f in [("bb144-injection-4", 4, 32)] + [
            (f"hgp225-injection-{i}", i, n)
            for i, n in enumerate([14, 19, 26, 25, 27, 28], 1)
        ]:
            s = example(name)
            r = receipt(s)["result"]
            self.assertEqual(r["status"], "proved")
            self.assertEqual(len(r["carriers"]), k)
            self.assertEqual(sum(map(len, r["postselect"].values())), f)
        s = example("bb144-injection-4")
        good = set(receipt(s)["result"]["fixed_checks"]["X"])
        s["postselect"]["X"].append(next(i for i in range(72) if i not in good))
        with self.assertRaisesRegex(ValueError, "random"):
            receipt(s)

    def test_search_both_objectives(self):
        s = example("bell-injection")
        for objective in ("protected", "fixed"):
            spec = dict(
                schema=1,
                kind="injection_search",
                memory=s["memory"],
                logicals=s["logicals"],
                count=2,
                objective=objective,
            )
            r = receipt(spec)
            verify(spec, r)
            self.assertEqual(
                receipt(r["result"]["source"])["result"]["status"], "proved"
            )
            self.assertEqual(
                receipt(spec, budget=1)["result"]["status"], "inconclusive"
            )
