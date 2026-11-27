import copy
import unittest
from importlib.util import find_spec
from fractions import Fraction
from ftqc_extractors.__main__ import receipt
from ftqc_extractors.injection import analyze
from ftqc_magic.catalog import get_catalog_factory
from .__main__ import evaluate


def study(mode="shared", p="1/1000", reject="0"):
    spec = dict(
        schema=1,
        kind="injection",
        method="general",
        memory=dict(n=6, hx=[63], hz=[63]),
        logicals=[[1 | (1 << j), (1 << j) | 32, 0, 0] for j in range(1, 5)],
    )
    carriers = analyze(spec)["carriers"]
    inputs = []
    for start in range(0, 15, 4):
        k = min(4, 15 - start)
        events = []
        if mode == "shared" and start == 0:
            events.append(dict(id="shared-123", sites=carriers[:3], p=p, reject=False))
        if mode in ("iid", "three-iid"):
            for j in range(k):
                if mode == "iid" or start + j < 3:
                    events.append(
                        dict(id=f"fault-{j}", sites=[carriers[j]], p=p, reject=False)
                    )
        if start == 0:
            events.append(dict(id="herald", sites=[], p=reject, reject=True))
        s = dict(
            schema=1,
            kind="injection_export",
            attempt=f"source-{start}",
            block=f"block-{start}",
            model=dict(type="raw-T", injection=spec, events=events),
        )
        inputs.append(
            dict(spec=s, receipt=receipt(s), ports=[[i, start + i] for i in range(k)])
        )
    return dict(kind="insitu-rm15", inputs=inputs, independent_attempts=True)


class DistillTests(unittest.TestCase):
    def test_all_masks(self):
        result = evaluate(study("iid"))["factory"]
        counts = [0, 0]
        for mask in range(32768):
            syndrome = 0
            for i in range(15):
                if mask >> i & 1:
                    syndrome ^= i + 1
            if not syndrome:
                counts[mask.bit_count() % 2] += 999 ** (15 - mask.bit_count())
        self.assertEqual(
            list(map(Fraction, result["accepted_mass"].values())),
            [Fraction(n, 1000**15) for n in counts],
        )
        native = get_catalog_factory("bk_15_1_3_t").to_source_contract("1/1000")
        self.assertEqual(Fraction(result["final_accept"]), native.acceptance_lower)
        self.assertEqual(
            Fraction(result["conditional_output_Z"]),
            native.conditional_infidelity_upper,
        )

    def test_common_fault_and_controls(self):
        s = study(reject="1/4")
        out = evaluate(s)["factory"]
        self.assertEqual(
            (
                out["source_accept"],
                out["final_accept"],
                out["accepted_wrong"],
                out["conditional_output_Z"],
            ),
            ("3/4", "3/4", "3/4000", "1/1000"),
        )
        self.assertNotEqual(
            evaluate(study("three-iid"))["factory"]["conditional_output_Z"], "1/1000"
        )
        self.assertIsNone(
            evaluate(study(reject="1"))["factory"]["conditional_output_Z"]
        )
        for change in ("reuse", "alias", "receipt", "dependence"):
            bad = copy.deepcopy(s)
            if change == "reuse":
                bad["inputs"][0]["ports"][1][0] = 0
            if change == "alias":
                bad["inputs"][0]["ports"][1][1] = 0
            if change == "receipt":
                bad["inputs"][0]["receipt"]["result"]["target"] = "Y"
            if change == "dependence":
                bad["independent_attempts"] = False
            with self.assertRaises(ValueError):
                evaluate(bad)

    def test_retry_and_handoff(self):
        from .__main__ import receipt as compose_receipt
        from .test_magic_handoff import study as handoff

        s = study(reject="1/4")
        out = evaluate(
            dict(
                kind="insitu-retries",
                factory=s,
                lane_ticks=[2, 3, 4, 5],
                factory_ticks=10,
                factory_width=5,
                transcript=["source-reject", "delivered"],
            )
        )
        self.assertEqual(out["expected_wire_ticks_per_delivery"], "210")
        self.assertEqual(out["wire_ticks"], 290)
        gate = handoff("postselect")
        gate.update(factory=s, factory_receipt=compose_receipt(s))
        r = evaluate(gate)
        self.assertEqual(r["instrument"]["transcripts"], 4)
        self.assertEqual(r["resource_trace_distance_upper"], "1/1000")

    def test_y_evidence_is_not_a_t_input(self):
        if find_spec("stim") is None:
            self.skipTest("optional Stim")
        from ftqc_extractors.injection_circuit import run

        circuit = run(dict(case="bb-uniform"))
        sample = dict(
            circuit_sha256=circuit["circuit_sha256"],
            postselected_detectors=circuit["postselect_width"],
            outputs=4,
            shots=1,
            accepted=1,
            accepted_mask_counts=[1] + [0] * 15,
        )
        for model in (
            dict(type="closed-Y", case="bb-uniform", locations=[], frame=[0] * 4),
            dict(type="sampled-Y", case="bb-uniform", sample=sample),
        ):
            s = study()
            row = s["inputs"][0]
            row["spec"]["model"] = model
            row["receipt"] = receipt(row["spec"])
            with self.assertRaisesRegex(ValueError, "exact T"):
                evaluate(s)
