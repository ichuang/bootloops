import copy
import unittest
from fractions import Fraction
from math import comb
from unittest.mock import patch

from ftqc_extractors.__main__ import receipt
from ftqc_extractors.magic import branches

from .logical_t import css_encoder, encoding, instrument, logical_t_result
from .test_distill import study


def example():
    gate = {
        "schema": 1,
        "kind": "t_gate",
        "a": ["11", "11"],
        "b": ["11", "11"],
        "logical_x": 3,
        "logical_z": 5,
        "rounds": 2,
        "target": 2,
        "resource": "ideal-encoded-T",
    }
    return {
        "kind": "distilled-logical-t",
        "factory": {"kind": "local-stochastic", "source": study(), "p": "1/100"},
        "gate": gate,
        "gate_receipt": receipt(gate),
        "handoff": "fresh-product-css",
        "encoder_p": 0,
        "online_p": 0,
    }


class LogicalTTests(unittest.TestCase):
    def test_ideal_and_noisy_handoff(self):
        s = example()
        r = logical_t_result(s)
        self.assertEqual(r["cases"]["gate"]["accepted_failure_upper"], "21/2000")
        self.assertEqual(r["fault_free_logical_T_half_diamond_upper"], "21/2000")
        self.assertEqual(r["encoder"]["logical_Z_witness"]["after"], "ZI")
        self.assertEqual(r["instrument"]["completeness_matrix_entries"], 16)
        s.update(encoder_p="1/1000000", online_p="1/100000")
        r = logical_t_result(s)
        delta = Fraction(21, 2000) + Fraction(21, 1000000)
        q = Fraction(s["online_p"])
        gate = r["cases"]["gate"]
        self.assertEqual(
            Fraction(gate["accepted_failure_upper"]),
            delta + comb(gate["locations"], 2) * q**2,
        )
        tail = sum(
            Fraction(r["cases"][k]["online_fault_tail_upper"])
            for k in ("read_X", "read_Y")
        )
        self.assertEqual(
            Fraction(r["adaptive_readout_accepted_failure_upper"]), delta + tail
        )
        s["encoder_p"] = 1
        self.assertEqual(
            logical_t_result(s)["cases"]["gate"]["accepted_failure_upper"], "1"
        )

    def test_wrong_port_phase_and_erasure(self):
        s = example()
        rows = copy.deepcopy(branches())
        rows[1]["columns"][1][1] ^= 4
        with self.assertRaisesRegex(ValueError, "phase"):
            instrument(rows)

        def omit_last(*args):
            return css_encoder(*args)[:-1]

        with (
            patch("ftqc_compose.logical_t.css_encoder", omit_last),
            self.assertRaisesRegex(ValueError, "erase"),
        ):
            encoding(s["gate"])
        s["gate"]["logical_x"] ^= 1
        with self.assertRaisesRegex(ValueError, "digest"):
            logical_t_result(s)

    def test_wrong_factory_and_failed_consumer(self):
        s = example()
        s["factory"]["source"]["output_check"] = "none"
        with self.assertRaisesRegex(ValueError, "15-to-1"):
            logical_t_result(s)
        s = example()
        s["gate"]["target"] = 3
        s["gate_receipt"] = receipt(s["gate"])
        with self.assertRaisesRegex(ValueError, "proved distance"):
            logical_t_result(s)
