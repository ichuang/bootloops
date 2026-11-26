import json
import unittest
from pathlib import Path
from .magic import (
    branches,
    check_classifier,
    check_readout_classifier,
    frame_observable,
)
from .__main__ import receipt, verify
from .t_faults import FaultCircuit


def example():
    return json.loads((Path(__file__).parent / "examples/t-gate.json").read_text())


class MagicTests(unittest.TestCase):
    def test_phase(self):
        check_classifier()
        check_readout_classifier()
        self.assertFalse(all(b["implements_T"] for b in branches(correction=False)))
        self.assertEqual(frame_observable(1, 0), ("Y", 1))
        self.assertEqual(frame_observable(1, 1), ("Y", 0))

    def test_schedule(self):
        s = example()
        for arm in (None, 0, 1):
            t = FaultCircuit(s, arm)
            self.assertEqual(t.forward, t.backward)
        r = receipt(s)
        self.assertEqual(
            [
                r["result"]["cases"][a]["exact_distance"]
                for a in ("gate", "read_X", "read_Y")
            ],
            [2] * 3,
        )
        verify(s, r)
        r["result"]["cases"]["gate"]["weight_two"][0][1] = 0
        with self.assertRaises(ValueError):
            verify(s, r)

    def test_controls(self):
        s = example()
        broken = dict(s, trailing_identity=[["CX", 4, 5]] * 2)
        self.assertEqual(
            receipt(broken)["result"]["cases"]["gate"]["evidence"]["upper"], 1
        )
        with self.assertRaises(ValueError):
            receipt(dict(s, logical_x=5))
        with self.assertRaises(ValueError):
            receipt(dict(s, trailing_identity=[["CX", 4, 5]]))
        self.assertEqual(receipt(s, budget=1)["result"]["status"], "inconclusive")
        src = dict(s, kind="t_source", resource="in-situ-all-plus")
        r = receipt(src)
        verify(src, r)
        self.assertEqual(r["result"]["exact_distance"], 1)

    def test_structural_bound_against_fault_search(self):
        rows = ["110", "011", "101"]
        reps = [
            [7 << (3 * i), sum(1 << (3 * j + i) for j in range(3))] for i in range(3)
        ]
        s = dict(
            example(),
            a=rows,
            b=rows,
            logical_x=7,
            logical_z=73,
            rounds=3,
            representatives=reps,
            verification_rounds=2,
            target=3,
        )
        cab = dict(s, kind="t_bounds")
        del cab["target"]
        r = receipt(cab)
        verify(cab, r)
        self.assertEqual(r["result"]["cases"]["gate"]["lower"], 3)
        finite = receipt(s)
        self.assertEqual(finite["result"]["cases"]["gate"]["evidence"]["lower"], 3)
        bad = dict(cab, representatives=[[1, 73]] + reps[1:])
        with self.assertRaises(ValueError):
            receipt(bad)
