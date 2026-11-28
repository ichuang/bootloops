import json
import unittest
from pathlib import Path
from fractions import Fraction

from ftqc_compose.__main__ import receipt
from ftqc_extractors.__main__ import (
    receipt as install_receipt,
    verify as verify_install,
)
from ftqc_gadgets import staged
from ftqc_magic.staged import prepare
from ftqc_uncompute.adapters.staged import price


def load(name):
    return json.loads(
        (Path(__file__).parent / "examples" / (name + ".json")).read_text()
    )


class StagedIntegration(unittest.TestCase):
    def test_compose_and_replay(self):
        study = load("switched-t")
        original = receipt(study)["result"]
        self.assertEqual(original["inputs"], original["outputs"])
        self.assertEqual(original["action"], "T")
        self.assertGreater(Fraction(original["failure_upper"]), 0)
        stage = staged.analyze(study["stages"][0])
        study["stages"][0] = stage
        self.assertEqual(receipt(study)["result"], original)
        stage["result"]["outputs"][0]["signs"] = "+-"
        with self.assertRaisesRegex(ValueError, "replay"):
            receipt(study)

    def test_bad_join_and_duplicate_locations(self):
        study = load("switched-t")
        study["stages"][1]["blocks"] = ["elsewhere"]
        with self.assertRaisesRegex(ValueError, "interface"):
            receipt(study)
        study["instances"][1] = study["instances"][0]
        with self.assertRaisesRegex(ValueError, "duplicate"):
            receipt(study)

    def test_decoded_state_and_price(self):
        study = load("switched-t")
        source = prepare(study)
        self.assertEqual(source["acceptance"], "1")
        self.assertEqual(source["port"]["code"], "steane")
        self.assertEqual(
            source["conditional_infidelity_upper"],
            receipt(study)["result"]["failure_upper"],
        )
        first, second = (price(study, n, 10, 20) for n in (0, 1))
        self.assertEqual(
            second["reserved_wire_ticks_upper"] - first["reserved_wire_ticks_upper"], 30
        )
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            price(study, -1, 0, 0)
        with self.assertRaisesRegex(ValueError, "closed Steane"):
            prepare(load("pieceable-ccz"))

    def test_installed_couplers(self):
        spec = load("switched-t")["stages"][0]
        g = staged.analyze(spec)
        w = g["result"]["workload"]
        sites = list(range(w["data"] + w["peak_cat_legs"]))
        edges = {
            (q, 15 + j)
            for h in spec["checks"]
            for j, q in enumerate(q for q in range(15) if h >> q & 1)
        }
        s = dict(
            schema=1,
            kind="staged_install",
            receipt=g,
            sites=sites,
            couplers=list(map(list, sorted(edges))),
            ccz_triples=[],
        )
        installed = install_receipt(s)
        verify_install(s, installed)
        s["couplers"].pop()
        with self.assertRaisesRegex(ValueError, "missing installed"):
            install_receipt(s)
