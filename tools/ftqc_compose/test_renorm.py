"""Tests for verified-gadget Definition 4.8 discharge and multivariate renormalization."""

import copy
from fractions import Fraction
import json
import sys
import tempfile
from pathlib import Path
from subprocess import run
import unittest

from ftqc_gadgets import certify as gcert
from ftqc_gadgets import library as glib
from ftqc_uncompute.certify import build_receipt_for_kmx
from ftqc_uncompute.ir_kickmix import KickMixCircuit

from .level import level_result
from .renorm import renormalize_result


def _mini_roles():
    return {
        "prep_zero": {
            "locations": 20,
            "role_locations": {"prep_zero": 8, "cnot": 8, "measure_z": 4},
            "support_pairs_by_role": {"cnot,cnot": 28, "prep_zero,cnot": 16, "cnot,measure_z": 8},
            "depolarizing_pairs_by_role": {"cnot,cnot": "14/3", "prep_zero,cnot": "8/3", "cnot,measure_z": "2"},
            "qubits": 11,
        },
        "cnot": {
            "locations": 40,
            "role_locations": {"prep_zero": 16, "cnot": 16, "measure_z": 8},
            "support_pairs_by_role": {"cnot,cnot": 120, "prep_zero,cnot": 64, "cnot,measure_z": 32},
            "depolarizing_pairs_by_role": {"cnot,cnot": "20", "prep_zero,cnot": "32/3", "cnot,measure_z": "8"},
            "qubits": 30,
        },
        "measure_z": {
            "locations": 15,
            "role_locations": {"prep_zero": 4, "cnot": 4, "measure_z": 7},
            "support_pairs_by_role": {"cnot,measure_z": 16, "measure_z,measure_z": 21},
            "depolarizing_pairs_by_role": {"cnot,measure_z": "4", "measure_z,measure_z": "7"},
            "qubits": 7,
        },
    }


class VerifiedGadgetAndRenormTests(unittest.TestCase):
    def test_verified_gadget_discharges_def_4_8_with_zero_premises(self):
        g_prep = glib.load("prep:steane,zero,goto,idles")
        rec_prep = gcert.build(g_prep, sector=True, coherent=[("prep", 1, 1)])
        g_meas = glib.load("measure:steane,Z")
        rec_meas = gcert.build(g_meas, sector=True, coherent=[("ec", 1, None)])

        spec = {
            "kind": "level",
            "inputs": {},
            "outputs": [],
            "physical_sites": ["p0", "m0"],
            "vertices": [
                {
                    "id": "u_prep",
                    "inputs": [],
                    "outputs": ["w0"],
                    "gadget": {
                        "op": "verified-gadget",
                        "receipt": rec_prep,
                        "ideal": "prep-zero",
                        "faults": {"op": "family", "sites": ["p0"], "minimal": [["p0"]]},
                    },
                },
                {
                    "id": "u_meas",
                    "inputs": ["w0"],
                    "outputs": [],
                    "gadget": {
                        "op": "verified-gadget",
                        "receipt": rec_meas,
                        "ideal": "measure-z",
                        "faults": {"op": "family", "sites": ["m0"], "minimal": [["m0"]]},
                    },
                },
            ],
            "outer": {"op": "family", "sites": ["u_prep", "u_meas"], "minimal": [["u_prep"], ["u_meas"]]},
            "output_check": "none",
            "require_closed": True,
        }
        res = level_result(spec)
        self.assertEqual(res["status"], "proved")
        self.assertEqual(res["premises"], [])
        self.assertEqual(res["bad_faults"]["enumerator"], {"1": 2})

        # Planted negative control: postselected gadget receipt is rejected by verified-gadget.
        g_post = glib.load("prep:steane,zero,tree=XZ")
        rec_post = gcert.build(g_post, sector=True, coherent=[("prep", 1, 1)])
        bad_spec = copy.deepcopy(spec)
        bad_spec["vertices"][0]["gadget"]["receipt"] = rec_post
        with self.assertRaisesRegex(ValueError, "coherent"):
            level_result(bad_spec)

    def test_renormalize_rm15_and_uncompute_workload(self):
        kmx_src = "\n".join([
            "REGISTER r0 \"in\"",
            "APPEND_TO_REGISTER q0 r0",
            "APPEND_TO_REGISTER q1 r0",
            "INFO_TAG \"@enter:init_and\"",
            "CCX q0 q1 q2",
            "INFO_TAG \"@exit:init_and\"",
            "INFO_TAG \"@enter:del_and\"",
            "HMR q2 b0",
            "CZ q0 q1 if b0",
            "INFO_TAG \"@exit:del_and\"",
            "",
        ])
        kmx_circ = KickMixCircuit.from_kmx(kmx_src)
        urec = build_receipt_for_kmx(kmx_circ, source_name="and_del.kmx", shots=32, run_kraus=True)

        base_spec = {
            "kind": "renormalize",
            "noise": "local-stochastic",
            "p": "1/10000",
            "roles": _mini_roles(),
            "cone_weights": "auto",
            "non_clifford": {
                "scheme": "rm15-distill",
                "harmful_pairs": 105,
                "total_slots": 15,
                "injection_multiplier": "2",
                "clifford_counts": {"prep_zero": 4, "cnot": 12, "measure_z": 4},
            },
            "max_levels": 4,
            "workload": {
                "logical_qubits": 3,
                "logical_counts": {"prep_zero": 3, "cnot": 10, "measure_z": 3, "magic": 4},
                "target_eps": "1/1000000",
                "uncompute": {"kmx": kmx_src, "receipt": urec},
            },
        }
        res_ls = renormalize_result(base_spec)
        self.assertEqual(res_ls["status"], "proved")
        self.assertTrue(res_ls["below_certified_threshold"])
        self.assertGreater(Fraction(res_ls["threshold_lower"]), Fraction(1, 10000))
        self.assertTrue(res_ls["workload"]["achieved"])
        self.assertTrue(res_ls["workload"]["uncompute"]["branch_injective"])

        dep_spec = dict(base_spec, noise="depolarizing")
        res_dep = renormalize_result(dep_spec)
        self.assertEqual(res_dep["status"], "proved")
        self.assertGreater(Fraction(res_dep["threshold_lower"]), Fraction(res_ls["threshold_lower"]))

        # Planted negative control 1: tampered kmx source fails uncompute binding.
        bad_kmx = copy.deepcopy(base_spec)
        bad_kmx["workload"]["uncompute"]["kmx"] = kmx_src.replace("CCX q0 q1 q2", "CX q0 q2")
        with self.assertRaisesRegex(ValueError, "hash does not match"):
            renormalize_result(bad_kmx)

        # Planted negative control 2: distance-1 raw encoder with positive linear floor
        # prevents asymptotic contraction to zero and returns counterexample below floor.
        floor_spec = {
            "kind": "renormalize",
            "noise": "local-stochastic",
            "p": "1/1000",
            "roles": _mini_roles(),
            "cone_weights": "auto",
            "non_clifford": {
                "scheme": "staged-gadget",
                "locations": 10,
                "role_locations": {"cnot": 10},
                "support_pairs_by_role": {"cnot,cnot": 45},
                "linear_floor_coeff": "3",
            },
            "max_levels": 5,
            "workload": {
                "logical_qubits": 2,
                "logical_counts": {"cnot": 2, "magic": 1},
                "target_eps": "1/1000000",
            },
        }
        res_floor = renormalize_result(floor_spec)
        self.assertEqual(res_floor["status"], "counterexample")
        self.assertEqual(res_floor["threshold_lower"], "0")
        self.assertFalse(res_floor["below_certified_threshold"])

        # CLI check + verify + tamper refusal:
        with tempfile.TemporaryDirectory() as td:
            sp = Path(td) / "spec.json"
            rp = Path(td) / "rec.json"
            sp.write_text(json.dumps(dep_spec))
            cmd = [sys.executable, "-m", "ftqc_compose"]
            run([*cmd, "check", sp, "-o", rp], check=True, capture_output=True)
            v_ok = run([*cmd, "verify", sp, rp], check=False, capture_output=True, text=True)
            self.assertEqual(v_ok.returncode, 0, v_ok.stderr)
            tampered = json.loads(rp.read_text())
            tampered["result"]["threshold_lower"] = "1/2"
            rp.write_text(json.dumps(tampered))
            v_bad = run([*cmd, "verify", sp, rp], check=False, capture_output=True, text=True)
            self.assertEqual(v_bad.returncode, 2)


if __name__ == "__main__":
    unittest.main()
