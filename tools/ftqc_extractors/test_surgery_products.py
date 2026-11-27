import unittest
from importlib.util import find_spec

from ftqc_gadgets.surgery import Instrument
from ftqc_gadgets.surgery_design import build
from ftqc_gadgets.surgery_faults import columns, residual_profile
from ftqc_gadgets.test_surgery import edge
from ftqc_extractors.__main__ import evaluate
from ftqc_extractors.installed_surgery import analyze
from ftqc_extractors.surgery_bench import two_block
from ftqc_extractors.surgery_products import product_sequence, radial_design


class ProductTests(unittest.TestCase):
    def test_reverse_transport(self):
        instrument = Instrument(edge())
        self.assertEqual(columns(instrument)[0], list(instrument.locations()))
        result = residual_profile(instrument, 1, 1, 100000)
        self.assertEqual(result["status"], "proved")
        self.assertEqual(len(result["completed"]), 4)
        self.assertEqual(
            residual_profile(instrument, 1, 1, 1)["status"], "inconclusive"
        )

    def test_padded(self):
        s = product_sequence([7], 3, [3, 3], 2, [3, 5], 2, [2, 2], 12)
        result = evaluate(s, "exact", 1000000, 60)
        self.assertEqual(result["distance_lower"], 2)
        s["padding"] = [2, 1]
        with self.assertRaisesRegex(ValueError, "both boundaries"):
            evaluate(s, "exact", 1000000, 60)

    @unittest.skipUnless(
        find_spec("stim"), "larger Choi probes require Stim for fast replay"
    )
    def test_radial(self):
        s = radial_design()
        r = build(s)
        self.assertEqual(r["geometry"]["cone"]["logical_qubits"], 7)
        Instrument(r["spec"])

    def test_encoded_handoff(self):
        self.assertEqual(analyze(two_block(), 100000)["gadget"]["faults"]["lower"], 2)
