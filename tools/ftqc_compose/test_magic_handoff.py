import copy
import unittest
from collections import Counter
from fractions import Fraction

from ftqc_magic.catalog import get_catalog_factory
from ftqc_magic.certify import build_receipt
from ftqc_magic.delivery import read_delivery
from .magic_handoff import consume_t, transfer_for
from .test_logical_t import example


def study(mode):
    s = example()
    s.pop("handoff")
    s.update(
        kind="magic-logical-t",
        factory=dict(kind="factory", catalog="bk_15_1_3_t", p_phys="1/1000", mode=mode),
        attempt="factory-0/attempt-1",
        port=dict(block="HGP8/resource", pair=[3, 5], frame="I", correction="I"),
    )
    s["factory_receipt"] = build_receipt(s["factory"])
    s["transfer_receipt"] = transfer_for(s["gate"], s["port"])
    return s


def enumerate_masses(name, correct):
    fac = get_catalog_factory(name)
    columns = [(fac.column_syndrome(i), fac.column_label(i)) for i in range(fac.n)]
    decoder = dict(columns) if correct else {}
    counts = Counter()
    witnesses = {}
    for faults in range(1 << fac.n):
        syn = label = 0
        for i, (check, logical) in enumerate(columns):
            if faults >> i & 1:
                syn ^= check
                label ^= logical
        if correct or not syn:
            mask = label ^ decoder.get(syn, 0)
            counts[mask] += 999 ** (fac.n - faults.bit_count())
            if mask and mask not in witnesses:
                witnesses[mask] = faults
    return {str(m): str(Fraction(v, 1000**fac.n)) for m, v in counts.items()}, witnesses


class HandoffTests(unittest.TestCase):
    def test_both_source_laws_and_instruments(self):
        for mode, correct, d in [("correcting", True, 2), ("postselect", False, 3)]:
            s = study(mode)
            s["encoder_p"] = "1/1000000"
            r = consume_t(s)
            source = r["delivery"]
            self.assertEqual(source["distance"], d)
            self.assertEqual(
                source["accepted_z_mass"], enumerate_masses("bk_15_1_3_t", correct)[0]
            )
            self.assertEqual(source["phases"], [0, 1])
            self.assertEqual(r["instrument"]["transcripts"], 4)
            self.assertEqual(r["encoder_any_fault_upper"], "21/1000000")
            self.assertEqual(r["encoder"]["raw_to_gate_accepted_fault_distance"], 1)
            if correct:
                self.assertEqual(source["leading_count"], 105)
                self.assertLessEqual(
                    Fraction(source["conditional_infidelity_upper"]),
                    105 * Fraction(1, 1000) ** 2,
                )
            else:
                self.assertLess(Fraction(source["acceptance_lower"]), 1)
            s["factory_receipt"]["details"] = {"target": "T-dagger", "law": "clean"}
            self.assertEqual(read_delivery(s["factory"], s["factory_receipt"]), source)

    def test_destination_and_producer_mutations(self):
        s = study("correcting")
        s["port"].update(frame="Z", correction="Z")
        s["transfer_receipt"] = transfer_for(s["gate"], s["port"])
        consume_t(s)
        for path, key, value in [
            ("port", "correction", "I"),
            ("port", "pair", [5, 3]),
            ("factory_receipt", "source_sha256", "0"),
            ("transfer_receipt", "operations", []),
        ]:
            bad = copy.deepcopy(s)
            bad[path][key] = value
            with self.assertRaises(ValueError):
                consume_t(bad)
        s["factory_receipt"]["contract"]["target_label"] = "T-dagger"
        with self.assertRaises(ValueError):
            consume_t(s)

    def test_forward_same_encoded_port(self):
        s = study("postselect")
        source = read_delivery(s["factory"], s["factory_receipt"])
        s["encoded_origin"] = dict(
            attempt=s["attempt"],
            factory=source["producer"]["spec_sha256"],
            port=s["port"],
            transfer=s["transfer_receipt"],
            wire_order=list(range(8)),
        )
        self.assertEqual(consume_t(s)["encoder"]["operation_count"], 21)
        s["encoded_origin"]["wire_order"].reverse()
        with self.assertRaisesRegex(ValueError, "origin"):
            consume_t(s)
