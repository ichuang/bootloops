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
        from ftqc_faults.audit.injection_records import parse

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

    def test_sample_masks_and_silent_signatures(self):
        from ftqc_faults.audit.injection_statistics import summarize
        from .injection_circuit import paths, run
        from fractions import Fraction

        s = dict(
            outputs=2, shots=200, accepted=100, accepted_mask_counts=[90, 0, 0, 10]
        )
        r = summarize(s, 2, 10, 2)
        self.assertEqual(r["correlated_error"]["estimate"], 0.1)
        self.assertEqual(r["marginals"][0]["estimate"], 0.1)
        self.assertEqual(r["batch_error_divided_by_outputs"], 0.05)
        with self.assertRaises(ValueError):
            summarize(dict(s, accepted=99), 2, 10, 2)
        stim = self.stim
        for case in ("bb-uniform", "bb-asymmetric", "hgp-release"):
            _, path, _, p = paths(case)
            circuit = stim.Circuit.from_file(path)
            result = run(dict(case=case))
            k = result["outputs"]
            filt = stim.DetectorErrorModel()
            for mask in range(1, 1 << k):
                filt.append(
                    "error",
                    1,
                    [
                        stim.target_logical_observable_id(i)
                        for i in range(k)
                        if mask >> i & 1
                    ],
                )
            for err in circuit.explain_detector_error_model_errors(dem_filter=filt):
                mask = sum(1 << x.dem_target.val for x in err.dem_error_terms)
                expected = Fraction()
                for loc in err.circuit_error_locations:
                    op = loc.instruction_targets
                    scale = {"DEPOLARIZE1": 3, "DEPOLARIZE2": 15}.get(op.gate, 1)
                    expected += Fraction(str(op.args[0])) / scale / p
                self.assertEqual(
                    Fraction(result["first_order_silent_by_mask"][str(mask)]), expected
                )


class ObjectiveTests(unittest.TestCase):
    def test_fixed_and_protected_closures(self):
        try:
            from scipy.optimize import milp
        except ImportError:
            self.skipTest("optional integer optimizer")
        from .injection_optimize import optimize

        # On two unconstrained spectators a fixed-check optimum is exactly two.
        s = dict(
            schema=1,
            kind="injection",
            method="selected",
            memory=dict(n=3, hx=[2], hz=[4]),
            logicals=[[1, 1, 0, 0]],
            preparation=dict(carriers=[0], plus=[2], zero=[1], bells=[]),
        )
        a = optimize(s, "fixed")
        self.assertEqual(a["postselect"], {"X": [0], "Z": [0]})
        self.assertEqual(
            receipt(optimize(s, "protected"))["result"]["status"], "proved"
        )


class IndependentSourceTests(unittest.TestCase):
    def test_general_full_rate(self):
        for name, k in [("bb144", 12), ("hgp225", 9)]:
            proof = receipt(example(name + "-injection-general"))
            self.assertEqual(proof["result"]["status"], "proved")
            self.assertEqual(len(proof["result"]["carriers"]), k)

    def test_signed_y_after_physical_projection(self):
        # Explicit unnormalized ten-qubit amplitudes. Project first onto X^10
        # and then Z^10. The omitted common normalization cancels.
        for name, expect in [("bell-injection", 1), ("split-bell-injection", 0)]:
            s = example(name)
            prep = s["preparation"]
            v = []
            for w in range(1024):
                if any(w >> q & 1 for q in prep["zero"]) or any(
                    (w >> a ^ w >> b) & 1 for a, b in prep["bells"]
                ):
                    v.append(0)
                else:
                    v.append(1j ** sum(w >> q & 1 for q in prep["carriers"]))
            projected = [
                (v[w] + v[w ^ 1023]) if w.bit_count() % 2 == 0 else 0
                for w in range(1024)
            ]
            norm = sum(abs(a) ** 2 for a in projected)
            for (x, z, _, _), target in zip(s["logicals"], (expect, 1)):
                value = sum(
                    a.conjugate()
                    * 1j
                    * (-1) ** ((w ^ x) & z).bit_count()
                    * projected[w ^ x]
                    for w, a in enumerate(projected)
                )
                self.assertAlmostEqual(value / norm, target)
