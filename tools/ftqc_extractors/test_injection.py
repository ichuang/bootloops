import json
import unittest
from pathlib import Path
from .__main__ import receipt, verify
from .algebra import parity

EXAMPLES = Path(__file__).parent / "examples"


def example(name):
    return json.loads((EXAMPLES / (name + ".json")).read_text())


class SourceTests(unittest.TestCase):
    def test_signed_general_and_syndrome_frames(self):
        s = example("general-injection")
        result = receipt(s)
        verify(s, result)
        r = result["result"]
        self.assertEqual(r["status"], "proved")
        self.assertTrue(r["carrier_gates"])
        for axis, rows in r["syndrome_frame"].items():
            checks = s["memory"]["hx" if axis == "Z" else "hz"]
            for idx, c in rows:
                self.assertEqual(parity(checks[idx], c), 1)
                self.assertFalse(
                    any(
                        parity(p[0 if axis == "Z" else 1], c)
                        for p in r["representatives"]
                    )
                )

    def test_overlapping_logicals(self):
        for name, status in (
            ("bell-injection", "proved"),
            ("split-bell-injection", "counterexample"),
        ):
            s = example(name)
            proof = receipt(s)
            verify(s, proof)
            self.assertEqual(proof["result"]["status"], status)

    def test_reference_sources_and_postselection(self):
        for name, k, f in [("bb144-injection-4", 4, 32)] + [
            (f"hgp225-injection-{i}", i, n)
            for i, n in enumerate([14, 19, 26, 25, 27, 28], 1)
        ]:
            s = example(name)
            r = receipt(s)["result"]
            self.assertEqual(r["status"], "proved")
            self.assertEqual(len(r["carriers"]), k)
            self.assertEqual(sum(map(len, r["postselect"].values())), f)
        s = example("bb144-injection-4")
        good = set(receipt(s)["result"]["fixed_checks"]["X"])
        s["postselect"]["X"].append(next(i for i in range(72) if i not in good))
        with self.assertRaisesRegex(ValueError, "random"):
            receipt(s)

    def test_search_both_objectives(self):
        s = example("bell-injection")
        for objective in ("protected", "fixed"):
            spec = dict(
                schema=1,
                kind="injection_search",
                memory=s["memory"],
                logicals=s["logicals"],
                count=2,
                objective=objective,
            )
            r = receipt(spec)
            verify(spec, r)
            self.assertEqual(
                receipt(r["result"]["source"])["result"]["status"], "proved"
            )
            self.assertEqual(
                receipt(spec, budget=1)["result"]["status"], "inconclusive"
            )


class CircuitTests(unittest.TestCase):
    def setUp(self):
        try:
            import stim
        except ImportError:
            self.skipTest("optional Stim circuit benchmark")
        self.stim = stim

    def test_reference_and_random_syndrome(self):
        from .injection_circuit import paths, run
        from ftqc_audit.injection_records import parse

        stim = self.stim
        for case, coefficient, raw in [
            ("bb-uniform", "94/15", 4),
            ("bb-asymmetric", "8/3", 4),
            ("hgp-release", "48", 0),
        ]:
            result = run(dict(case=case))
            self.assertEqual(result["first_order_total"], coefficient)
            self.assertEqual(len(result["explicit_raw_flip_witnesses"]), raw)
        _, path, _, _ = paths("bb-uniform")
        circuit = stim.Circuit.from_file(path).without_noise()
        rec, ops, det, obs = parse(circuit)
        measured = circuit.compile_sampler(seed=445).sample(64)
        fixed = {row[0] for row in det[:32]}
        self.assertFalse(measured[:, sorted(fixed)].any())
        self.assertTrue(measured[:, sorted(set(range(144)) - fixed)].any())
        # Forward replay a raw X after injecting the first carrier.
        modified = stim.Circuit()
        for inst in circuit:
            modified.append(inst)
            if inst.name == "RY" and 64 in [t.value for t in inst.targets_copy()]:
                modified.append("X", 64)
        actual = modified.compile_sampler(seed=3).sample(2)
        detectors, outputs = circuit.compile_m2d_converter().convert(
            measurements=actual, separate_observables=True
        )
        self.assertFalse(detectors.any())
        self.assertEqual(outputs.tolist(), [[True, False, False, False]] * 2)
