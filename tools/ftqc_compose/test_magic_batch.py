import json
import unittest
from fractions import Fraction
from pathlib import Path

from ftqc_magic.certify import build_receipt
from ftqc_magic.delivery import read_delivery
from .test_magic_handoff import enumerate_masses
from .magic_batch import consume_batch


def batch(source, groups):
    return dict(
        kind="magic-batch",
        factory=source,
        factory_receipt=build_receipt(source),
        batch="factory-A/attempt-1",
        ports=[
            dict(member=i, wires=g, data=[f"data.{w}" for w in g])
            for i, g in enumerate(groups)
        ],
    )


def qldpc_spec():
    return json.loads(
        (
            Path(__file__).resolve().parents[1]
            / "ftqc_magic/examples/qldpc_tricycle_quadcycle_ccz.json"
        ).read_text()
    )


class BatchTests(unittest.TestCase):
    def test_shared_t_fault(self):
        spec = batch(
            dict(kind="factory", catalog="bh_14_2_2_t", p_phys="1/1000"), [[0], [1]]
        )
        r = consume_batch(spec)
        masses, witness = enumerate_masses("bh_14_2_2_t", False)
        self.assertEqual(
            {str(a["mask"]): a["accepted_mass"] for a in r["conditional_joint_atoms"]},
            masses,
        )
        self.assertTrue(witness[3])
        self.assertEqual(r["source"]["accepted_fault_witness"]["output_mask"], 3)
        self.assertEqual(len(r["source"]["accepted_fault_witness"]["input_columns"]), 2)
        both = next(a for a in r["conditional_joint_atoms"] if a["mask"] == 3)
        self.assertEqual(both["data_Z"], ["data.0", "data.1"])
        self.assertGreater(Fraction(both["conditional_mass"]), 0)
        spec["ports"][1]["data"] = ["data.0"]
        with self.assertRaises(ValueError):
            consume_batch(spec)

    def test_ccz_branches_and_catalyst(self):
        source = qldpc_spec()
        spec = batch(source, [[0, 1, 2], [3, 4, 5]])
        r = consume_batch(spec)
        self.assertEqual(len(r["instrument"]["branches"]), 8)
        row = r["instrument"]["branches"][1]
        self.assertEqual(row["cz"], [[1, 2]])
        self.assertIsNone(r["conditional_joint_atoms"])
        delta = Fraction(r["conditional_joint_trace_distance_upper"])
        self.assertGreaterEqual(
            delta**2, Fraction(r["source"]["conditional_infidelity_upper"])
        )
        spec["ports"][1]["wires"].reverse()
        with self.assertRaises(ValueError):
            consume_batch(spec)
        source = dict(
            kind="pipeline",
            stages=[source, dict(kind="catalyzed_c2t", catalyst_infidelity="0")],
        )
        spec = batch(source, [[i] for i in range(4)])
        spec["catalyst"] = dict(id="cat.0", history=[], retired_after=spec["batch"])
        self.assertIsNone(consume_batch(spec)["conditional_joint_atoms"])
        spec["catalyst"]["history"] = ["previous-attempt"]
        with self.assertRaises(ValueError):
            consume_batch(spec)

    def test_unverified_correlation_transfer(self):
        src = dict(
            kind="pipeline",
            stages=[
                dict(kind="distillation", factory="bh_14_2_2_t"),
                dict(
                    kind="distillation",
                    factory="bk_15_1_3_t",
                    correlation_propagation="joint_3twirl",
                ),
            ],
        )
        native = build_receipt(src)
        with self.assertRaisesRegex(ValueError, "joint transfer"):
            read_delivery(src, native)

    def test_full_twirl_and_distinct_instances(self):
        from .magic_correlations import ccz_twirl, interleaved

        law = ccz_twirl({"0": "1/2", "9": "1/2"}, 2)
        self.assertEqual(law["0"], "1/2")
        self.assertEqual(len(law), 50)
        self.assertEqual({v for m, v in law.items() if m != "0"}, {"1/98"})
        ccz = batch(dict(kind="factory", catalog="cube_8_3_2_ccz"), [[0, 1, 2]])
        ccz["twirl"] = "joint_3twirl"
        atoms = consume_batch(ccz)["conditional_joint_atoms"]
        self.assertEqual(len(atoms), 8)
        self.assertEqual(len({a["conditional_mass"] for a in atoms if a["mask"]}), 1)
        s = batch(dict(kind="factory", catalog="bh_14_2_2_t"), [[0], [1]])
        rows = [
            dict(
                id=f"lane{i}",
                batch="shot1",
                member=i,
                factory=s["factory"],
                receipt=s["factory_receipt"],
                data=f"d{i}",
            )
            for i in range(2)
        ]
        spec = dict(kind="magic-interleaved", independent_noise=True, instances=rows)
        out = interleaved(spec)["conditional_data_Z"]
        self.assertEqual(
            Fraction(out["3"]) * Fraction(out["0"]),
            Fraction(out["1"]) * Fraction(out["2"]),
        )
        rows[1]["id"] = rows[0]["id"]
        with self.assertRaises(ValueError):
            interleaved(spec)
