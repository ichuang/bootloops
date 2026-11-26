import unittest
from unittest.mock import patch

from ftqc_audit.circuit import Circuit
from ftqc_audit.record_model import RecordModel
from ftqc_extractors.__main__ import receipt

from .gross_residual import check_kernel, gross_result


class GrossResidualTests(unittest.TestCase):
    def test_random_response_is_only_a_displacement(self):
        c = Circuit([[("RX", 0)], [("M", 0)]], detectors=[], outputs=[])
        model = RecordModel(c, response_records=[0])
        self.assertEqual(model.replay([[0, "X"]]), (0, 1))
        c.detectors.append([0])
        with self.assertRaisesRegex(ValueError, "random declared"):
            RecordModel(c, response_records=[0])

    def test_rank_and_orthogonality(self):
        check_kernel([1, 2], [4], 3)
        for rows in ([1], [1, 4]):
            with self.assertRaisesRegex(ValueError, "quotient"):
                check_kernel(rows, [4], 3)

    def test_frame_and_distance_binding(self):
        source = {"schema": 1, "kind": "gross_circuit", "rounds": 1, "target": 1}
        spec = {
            "kind": "gross-residual",
            "extractor": source,
            "distance_receipt": receipt(source),
        }
        r = gross_result(spec)
        self.assertEqual(r["closed_output_through"], 0)
        out = r["open_output"]
        self.assertGreater(out["accepted_with_nontrivial_residual"], 0)
        self.assertGreater(out["accepted_with_frame_update"], 0)
        from ftqc_extractors.gross import port_frame

        def missing_frame(*args):
            paths, cycles = port_frame(*args)
            return [0] * len(paths), cycles

        with (
            patch("ftqc_compose.gross_residual.port_frame", missing_frame),
            self.assertRaisesRegex(ValueError, "terminal probe"),
        ):
            gross_result(spec)
        spec["distance_receipt"]["result"]["lower"] += 1
        with self.assertRaisesRegex(ValueError, "replay mismatch"):
            gross_result(spec)
