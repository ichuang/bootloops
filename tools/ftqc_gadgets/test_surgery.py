import copy
import json
import unittest
from pathlib import Path

from ftqc_gadgets.surgery import Instrument, receipt, verify
from ftqc_extractors.surgery_export import export

EXAMPLES = Path(__file__).resolve().parents[1] / "ftqc_extractors" / "examples"


def edge(name="edge-extractor"):
    return export(json.loads((EXAMPLES / (name + ".json")).read_text()))


class SurgeryTests(unittest.TestCase):
    def test_edge_and_hook(self):
        for name, d in [("edge-extractor", 2), ("edge-hook", 1)]:
            r = receipt(edge(name), 2)
            self.assertEqual((r["faults"]["lower"], r["faults"]["upper"]), (d, d))
            verify(r)
            r["faults"]["upper"] += 1
            with self.assertRaises(ValueError):
                verify(r)

    def test_private_records_and_signed_readout(self):
        s = edge()
        for field in ["feedback", "record", "sign"]:
            t = copy.deepcopy(s)
            if field == "feedback":
                t["ticks"].pop()
            elif field == "record":
                t["readouts"][0]["records"] = [0]
            else:
                t["targets"][0] = "-" + t["targets"][0][1:]
            with self.assertRaises(ValueError):
                Instrument(t)
        t["readouts"][0]["sign"] = 1
        Instrument(t)

    def test_open_budget(self):
        self.assertEqual(receipt(edge(), 2, 1)["faults"]["status"], "inconclusive")

    def test_named_survivor(self):
        from .surgery_parallel import branch

        s = branch(2, [], ["+ZI"])
        s["basis"] = dict(
            input=[
                dict(name="first", x="+XI", z="+ZI"),
                dict(name="second", x="+IX", z="+IZ"),
            ],
            surviving=["second"],
        )
        Instrument(s)
        s["basis"]["surviving"] = ["first"]
        with self.assertRaisesRegex(ValueError, "survivor"):
            Instrument(s)
