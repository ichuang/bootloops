"""Unit tests for Milestone 5 addressable `LP(6,2)` logical-T handoff and negative controls."""

import copy
from fractions import Fraction
import unittest

from ftqc_codes import basis as fc_basis, certify as fc_certify, codes as fc_codes
from ftqc_magic.certify import build_receipt as build_magic_receipt
from ftqc_compose.addressable_magic import (
    compare_lp6_2_vs_hgp8_logical_t_control,
    consume_addressable_lp_t,
    reject_unjustified_correlated_batch_as_independent,
)
from ftqc_compose.magic_correlations import interleaved


def _make_lp6_t_spec(mode="correcting"):
    lp6 = fc_codes.load("lp6_2")
    lp6_rec = fc_certify.build(lp6)
    lp6_basis = fc_basis.derive_paired_basis(lp6)
    fac = {"kind": "factory", "catalog": "bk_15_1_3_t", "p_phys": "1/1000", "mode": mode}
    fac_rec = build_magic_receipt(fac)
    return {
        "kind": "addressable-magic-t",
        "code_receipt": lp6_rec,
        "basis_receipt": lp6_basis,
        "data_slot": "q0",
        "resource_slot": "q1",
        "factory": fac,
        "factory_receipt": fac_rec,
        "attempts": [
            {"id": "att-1", "start": 0, "end": 15, "wires": list(range(6, 11)), "raw": 6, "accepted": True}
        ],
        "delivery_request": {"attempt": "att-1", "slot": "q1", "arrival_tick": 17, "accepted": True},
        "frame": "I",
        "correction": "I",
        "latency_ticks": 2,
        "transfer_p": "1/1000000",
        "online_p": "1/100000",
        "surgery_rounds": 2,
    }


class AddressableMagicHandoffTests(unittest.TestCase):
    def test_lp6_2_single_output_t_delivery_branches_and_hgp8_reconciliation(self):
        comp = compare_lp6_2_vs_hgp8_logical_t_control()
        for mode in ("correcting", "postselect"):
            self.assertTrue(comp[mode]["source_matches_hgp8_control"])
            self.assertEqual(comp[mode]["lp6_2_addressable"]["shared_block_qubits"], 6)
            self.assertEqual(comp[mode]["hgp8_control"]["data_plus_resource_qubits"], 16)
        # Postselected mode includes rejected retry occupancy and strictly < 1 acceptance
        s_post = _make_lp6_t_spec("postselect")
        s_post["attempts"].insert(
            0, {"id": "att-0", "start": 0, "end": 15, "wires": list(range(6, 11)), "raw": 6, "accepted": False}
        )
        s_post["attempts"][1]["start"] = 15
        s_post["attempts"][1]["end"] = 30
        s_post["delivery_request"]["arrival_tick"] = 32
        res = consume_addressable_lp_t(s_post)
        self.assertEqual(res["status"], "proved")
        self.assertEqual(res["instrument"]["branches_verified"], 4)
        self.assertEqual(res["slots"]["unchanged_slots_during_swap"], ["q0"])
        self.assertEqual(res["spacetime_reconciliation"]["rejected_retries_in_transcript"], 1)
        self.assertLess(
            Fraction(res["error_propagation"]["unnormalized_attempt_joint_total_failure_upper"]),
            Fraction(res["error_propagation"]["conditional_total_failure_upper"]),
        )

        # Verified Z-frame with matching Z-correction also passes (7 swap CNOTs + 2 frame ops = 9)
        s_z = _make_lp6_t_spec("correcting")
        s_z["frame"] = "Z"
        s_z["correction"] = "Z"
        self.assertEqual(consume_addressable_lp_t(s_z)["transfer"]["swap_operations"], 9)

    def test_negative_controls_ccz_reversed_frame_and_correlated_batch(self):
        base = _make_lp6_t_spec("correcting")

        # 1. Reject CCZ-for-T substitution
        bad_ccz = copy.deepcopy(base)
        bad_ccz["factory"] = {"kind": "factory", "catalog": "cube_8_3_2_ccz", "p_phys": "1/1000"}
        bad_ccz["factory_receipt"] = build_magic_receipt(bad_ccz["factory"])
        with self.assertRaisesRegex(ValueError, "one verified T"):
            consume_addressable_lp_t(bad_ccz)

        # 2. Reject unproved reversed frame (`frame='Z', correction='I'`)
        bad_frame = copy.deepcopy(base)
        bad_frame["frame"] = "Z"
        bad_frame["correction"] = "I"
        with self.assertRaisesRegex(ValueError, "frame"):
            consume_addressable_lp_t(bad_frame)

        # 3. Reject treating two BATCH_CORRELATED outputs as independent!
        # (a) Multi-output correlated factory (`bh_14_2_2_t` produces 2 correlated T states in one batch!)
        corr_fac = {"kind": "factory", "catalog": "bh_14_2_2_t", "p_phys": "1/1000", "mode": "postselect"}
        corr_rec = build_magic_receipt(corr_fac)
        with self.assertRaisesRegex(ValueError, "BATCH_CORRELATED|one verified T"):
            bad_batch = copy.deepcopy(base)
            bad_batch["factory"] = corr_fac
            bad_batch["factory_receipt"] = corr_rec
            consume_addressable_lp_t(bad_batch)

        with self.assertRaisesRegex(ValueError, "BATCH_CORRELATED|cannot consume"):
            reject_unjustified_correlated_batch_as_independent(corr_fac, corr_rec)

        # (b) Also verify `interleaved` rejects reusing the same underlying batch attempt across 2 slots!
        corr_rec_shared = copy.deepcopy(corr_rec)
        corr_rec_shared.setdefault("delivery", {})
        with self.assertRaises(ValueError):
            # Reusing identical instance id or shared source_attempts is refused
            interleaved(
                {
                    "kind": "magic-interleaved",
                    "independent_noise": True,
                    "instances": [
                        {"id": "same-attempt", "batch": "b0", "member": 0, "factory": corr_fac, "receipt": corr_rec, "data": "q0"},
                        {"id": "same-attempt", "batch": "b0", "member": 1, "factory": corr_fac, "receipt": corr_rec, "data": "q1"},
                    ],
                }
            )


if __name__ == "__main__":
    unittest.main()
