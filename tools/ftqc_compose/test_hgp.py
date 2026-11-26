import unittest

from .__main__ import receipt
from .contracts import SCOPE, contract_result


class HGPContracts(unittest.TestCase):
    def test_composed_measurements(self):
        result = receipt({"kind": "hgp-repeat", "variant": "edge-extractor"})["result"]
        self.assertEqual(result["composition"]["bad_faults"]["enumerator"], {"2": 342})
        flat = result["flat_check"]
        self.assertEqual((result["status"], flat["distance"]), ("proved", 2))
        self.assertEqual(
            flat["counts"],
            {
                "weight_0": 1,
                "weight_1": 266,
                "weight_2": 33559,
                "accepted_harmful": 408,
            },
        )
        self.assertEqual(flat["witness"]["first"]["readout_flip"], 1)
        self.assertEqual(len(result["composition"]["outputs"]), 3)
        self.assertEqual(result["composition"]["premises"], [])

    def test_known_hook_fails_on_spectator(self):
        result = receipt({"kind": "hgp-repeat", "variant": "edge-hook"})["result"]
        self.assertEqual(result["status"], "counterexample")
        self.assertEqual(result["witness"]["fault"], [19, 8])
        self.assertEqual(result["flat_check"]["detectors"], [])
        self.assertEqual(result["flat_check"]["second"]["survivor_X_flip"], 1)
        self.assertTrue(result["witness"]["audit"]["undetected_harm"])

    def test_unproved_quantum_filter(self):
        spec = {
            "kind": "contract",
            "scope": SCOPE,
            "require_friendly": True,
            "require_closed": True,
            "program": {
                "op": "hgp-z",
                "instance": "meter",
                "variant": "edge-extractor",
            },
        }
        with self.assertRaisesRegex(ValueError, "friendly"):
            contract_result(spec)
