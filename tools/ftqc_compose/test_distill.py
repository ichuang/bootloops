import copy
import unittest
from itertools import combinations
from unittest.mock import patch

from .distill import code_check
from .level import level_result
from .stochastic import stochastic_result


def study():
    labels = [f"T{i:02d}" for i in range(1, 16)]
    sources = [
        {
            "id": k,
            "inputs": [],
            "outputs": ["wire-" + k],
            "gadget": {"op": "magic-source", "site": "prep-" + k},
        }
        for k in labels
    ]
    sink = {
        "id": "distiller",
        "inputs": ["wire-" + k for k in labels],
        "outputs": ["out"],
        "gadget": {"op": "rm15-ideal", "site": "ideal-processing"},
    }
    return {
        "kind": "level",
        "vertices": [*sources, sink],
        "physical_sites": ["prep-" + k for k in labels] + ["ideal-processing"],
        "inputs": {},
        "outputs": ["out"],
        "require_closed": True,
        "output_check": "distill15",
        "outer": {
            "op": "family",
            "sites": [*labels, "distiller"],
            "minimal": [list(t) for t in combinations(labels, 2)] + [["distiller"]],
        },
    }


class DistillationTests(unittest.TestCase):
    def test_code_and_channel(self):
        r = code_check()
        self.assertEqual(r["code"], [15, 1, 3])
        self.assertEqual(r["single_site_pauli_pairs"], 240)
        self.assertEqual(r["decoded_injection_kraus_checks"], 240)
        self.assertEqual(r["harmful_input_Z_pairs"], 105)
        self.assertEqual(r["pair_witness"]["decoder_Z_label"], 3)

    def test_wrong_injection_phase(self):
        # Omit the measurement-dependent S correction.
        with (
            patch(
                "ftqc_compose.distill.teleport_residual",
                lambda b, c, x, z: ((b ^ c ^ x) * (1 + 4 * z) - b) % 8,
            ),
            self.assertRaisesRegex(ValueError, "injection disagrees"),
        ):
            code_check()

    def test_distill_bound(self):
        s = study()
        r = level_result(s)
        self.assertEqual(r["bad_faults"]["enumerator"], {"2": 105})
        self.assertEqual(r["ideal_physical_sites"], ["ideal-processing"])
        bound = stochastic_result(
            {"kind": "local-stochastic", "p": "0.01", "source": s}
        )
        self.assertEqual(bound["upper"], "21/2000")
        self.assertTrue(bound["also_bounds"])
        for term in (["T01", "T02"], ["distiller"]):
            broken = copy.deepcopy(s)
            broken["outer"]["minimal"].remove(term)
            with self.assertRaisesRegex(ValueError, "resource pair or distiller"):
                level_result(broken)
