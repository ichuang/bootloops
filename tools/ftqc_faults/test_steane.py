import unittest
from unittest.mock import patch

from . import steane


class SteaneBridgeTests(unittest.TestCase):
    def test_two_controls_and_full_coverage(self):
        for result in steane.receipt()["axes"]:
            self.assertEqual(result["cnots"], 15)
            self.assertEqual(result["locations"], 95)
            self.assertEqual(result["single_fault_modes"], 455)
            self.assertEqual(len(result["detectors"]), 4)
            self.assertEqual(result["equivalence"]["verdict"], "FAULT-EQUIV-PASS")

    def test_changed_output_pauli_is_detected(self):
        original = steane.propagate_single_fault

        def flipped(*args):
            x, z, bits = original(*args)
            return x ^ 1, z, bits

        with patch.object(steane, "propagate_single_fault", flipped):
            with self.assertRaisesRegex(ValueError, "single-fault disagreement"):
                steane.compare_axis("X")

    def test_changed_record_is_detected(self):
        original = steane.propagate_single_fault

        def swapped(*args):
            x, z, bits = original(*args)
            return x, z, (bits & ~3) | ((bits & 1) << 1) | ((bits & 2) >> 1)

        with patch.object(steane, "propagate_single_fault", swapped):
            with self.assertRaisesRegex(ValueError, "single-fault disagreement"):
                steane.compare_axis("Z")
