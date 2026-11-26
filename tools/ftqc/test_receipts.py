import copy
import json
import unittest
from pathlib import Path
from qeccert import sat
from .__main__ import receipt, verify
from .circuit import Schedule
from .proof import formula

EXAMPLES = Path(__file__).parent / "examples"


def example(name):
    return json.loads((EXAMPLES / (name + ".json")).read_text())


class ReceiptTests(unittest.TestCase):
    def test_edge(self):
        s = example("edge-extractor")
        for spec, d in [(s, 2), (example("edge-hook"), 1)]:
            schedule = Schedule(spec)
            self.assertEqual(schedule.events(), schedule.events(reverse=True))
            wide = dict(spec, target=3)
            r = receipt(wide)
            self.assertEqual((r["result"]["lower"], r["result"]["upper"]), (d, d))
            verify(wide, r)
            r["result"]["faults"][0][1] = 0
            with self.assertRaises(ValueError):
                verify(wide, r)

    def test_circuit_controls(self):
        s = example("edge-extractor")
        self.assertEqual(receipt(s, budget=1)["result"]["status"], "inconclusive")
        with self.assertRaises(ValueError):
            receipt(dict(s, readout=[0]))
        with self.assertRaises(ValueError):
            receipt(dict(s, detectors=[[0]]))
        with self.assertRaises(ValueError):
            receipt(dict(s, ops=s["ops"][:-1]))
        r = receipt(s)
        with self.assertRaises(ValueError):
            verify(dict(s, target=3), r)
        with self.assertRaises(ValueError):
            receipt(dict(s, noise="unspecified"))
        with self.assertRaises(ValueError):
            verify(s, receipt(s, budget=0))

    def test_cnf_direction(self):
        s = Schedule(example("edge-hook"))
        p, labels = formula(s, 2)
        for j, (t, m) in enumerate(labels):
            v = s.signature({t: m})
            self.assertEqual(
                p.satisfied_by(p.extend(1 << j)), bool(v and not v & s.mask)
            )
        # Two different modes of one gate may never be charged as two sites.
        p, labels = formula(s, 3)
        j = next(i for i in range(len(labels) - 1) if labels[i][0] == labels[i + 1][0])
        self.assertFalse(p.satisfied_by(p.extend(3 << j)))

    def test_banked_lrat(self):
        s = example("repeat-z")
        r = example("repeat-z-lrat")
        verify(s, r)
        broken = copy.deepcopy(r)
        broken["result"]["lrat"] = ""
        with self.assertRaises(ValueError):
            verify(s, broken)

    @unittest.skipUnless(sat.find_cadical(), "optional CaDiCaL generation")
    def test_solver(self):
        for name in ["edge-extractor", "edge-hook"]:
            s = example(name)
            r = receipt(s, "sat")
            self.assertIn(r["result"]["status"], ["proved", "counterexample"])
            verify(s, r)

    def test_static(self):
        s = example("column")
        verify(s, receipt(s))
        self.assertEqual(receipt(s)["result"]["lower"], 2)
        for base in [
            dict(a=s["a"], b=s["b"]),
            dict(a=["010", "110"], b=["01", "01"]),
            dict(a=["1"], b=["1"]),
        ]:
            code = dict(schema=1, kind="hgp", **base)
            verify(code, receipt(code))

    def test_reduced(self):
        s = example("triangle")
        r = receipt(s)
        verify(s, r)
        self.assertEqual(r["result"]["lower"], 2)
        r["result"]["single_columns"] += 1
        with self.assertRaises(ValueError):
            verify(s, r)
