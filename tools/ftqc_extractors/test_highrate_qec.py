"""Unit tests for Milestone 4 addressable high-rate QEC surgery, shims, and composed studies."""

import unittest

from ftqc_codes import constructions as fc_constructions
from ftqc_extractors import addressable_surgery, algebra, surgery_products


class HighRateQECTests(unittest.TestCase):
    def test_compatibility_shims_exact_parity(self):
        # 1. algebra.hgp == constructions.hgp
        for a_str, b_str in [(["110", "011"], ["1011", "0110"]), (["11"], ["11", "11"])]:
            a, na = algebra.seed(a_str)
            b, nb = algebra.seed(b_str)
            tup = algebra.hgp(a, na, b, nb)
            c_ref = fc_constructions.hgp(a, na, b, nb)
            self.assertEqual(tup, (c_ref.n, c_ref.hx, c_ref.hz))

        # 2. surgery_products.lifted_product & radial90 == constructions
        lp_tup = surgery_products.lifted_product([[3]], [[3]], 3)
        lp_ref = fc_constructions.lifted_product([[3]], [[3]], 3)
        self.assertEqual(lp_tup, (lp_ref.n, lp_ref.hx, lp_ref.hz))
        self.assertTrue(lp_ref.sha256().startswith("c6860dd462fe"))

        r_tup = surgery_products.radial90()
        r_ref = fc_constructions.radial90()
        self.assertEqual(r_tup, (r_ref.n, r_ref.hx, r_ref.hz))
        self.assertTrue(r_ref.sha256().startswith("b656926c5959"))

    def test_lp6_2_addressable_end_to_end_and_negative_controls(self):
        res = addressable_surgery.replay_lp6_2_addressable_experiment()
        self.assertTrue(res["code_sha256"].startswith("c6860dd462fe"))
        self.assertNotEqual(res["base_basis_sha256"], res["dressed_basis_sha256"])
        self.assertEqual(res["base_schedule"]["branch_verification"]["total_verified_branches"], 64)
        self.assertEqual(res["dressed_schedule"]["branch_verification"]["total_verified_branches"], 64)
        self.assertEqual(res["base_schedule"]["faults"]["status"], "proved")
        self.assertEqual(res["dressed_schedule"]["faults"]["status"], "proved")
        self.assertEqual(
            set(res["negative_controls"]),
            {
                "slot_swap_without_conversion",
                "stale_code_digest",
                "stale_basis_digest",
                "reversed_record_sign",
                "omitted_boundary_idle",
            },
        )

    def test_translated_seed_substitution_on_lp8_4(self):
        comp = addressable_surgery.compare_translated_seed_vs_ordinary_surgery()
        self.assertEqual(comp["seed_orbit_certificate"]["covered_rank"], 2)
        self.assertEqual(comp["translated_seed_gadget"]["faults"]["status"], "proved")
        self.assertEqual(comp["ordinary_basis_surgery"]["faults"]["status"], "proved")
        self.assertLessEqual(
            comp["translated_seed_gadget"]["installed_couplers"],
            comp["ordinary_basis_surgery"]["installed_couplers"],
        )

    def test_full_composed_chain_faults_decoders_compose_uncompute(self):
        chain = addressable_surgery.run_full_composed_chain()
        self.assertTrue(chain["audit"]["verified"])
        self.assertEqual(chain["audit"]["num_detectors"], chain["decoders"]["dem_checks"])
        for alg_res in chain["decoders"]["algorithms"].values():
            self.assertTrue(alg_res["static_converged"])
            self.assertTrue(alg_res["circuit_dem_converged"])
        self.assertEqual(len(chain["compose"]["leaf1_surviving_basis"]), 1)
        self.assertEqual(chain["compose"]["leaf1_surviving_basis"][0]["name"], "q1")
        self.assertEqual(chain["compose"]["channel_status"], "proved")
        self.assertTrue(chain["compose"]["boundary_wires_single_charged"])
        self.assertLess(
            chain["uncompute_pricing"]["multitier_coupler"]["total_wall_ticks_ceil"],
            chain["uncompute_pricing"]["reconfigurable_shuttle"]["total_wall_ticks_ceil"],
        )

    def test_scaleup_studies_bb144_radial90_pp256(self):
        bb = addressable_surgery.build_bb144_decoder_study(max_splits=3)
        self.assertTrue(bb["code_sha256"].startswith("0ae98ac8a6e5"))
        self.assertTrue(all(bb["trapping_set_probe"]["healed_by"].values()))

        rad = addressable_surgery.build_radial90_map_comparison_study()
        self.assertTrue(rad["same_logical_coset_verified"])
        self.assertGreater(rad["degree_reduction"]["max_check_delta"], 0)
        self.assertGreater(rad["degree_reduction"]["cx_gate_delta"], 0)

        pp = addressable_surgery.build_pp256_orbit_surgery_study()
        self.assertTrue(pp["same_34_generator_signed_x_subgroup_verified"])
        self.assertEqual(pp["space_accounting"]["bare_code_rate_k_over_n"], "68/256")
        self.assertEqual(pp["space_accounting"]["parallel_merged_qubits_per_surviving_logical"], "576/34")
        self.assertEqual(pp["space_accounting"]["parallel_physical_qubits_per_surviving_logical"], "1120/34")
        self.assertIsNone(pp["claims_by_model"]["static"]["parallel_spatial_m2"]["proved_lower"])
        for dec_entry in pp["claims_by_model"]["decoded"].values():
            self.assertTrue(dec_entry["static_w3_logical_ok"])
            self.assertTrue(dec_entry["phenom_3round_logical_ok"])


if __name__ == "__main__":
    unittest.main()
