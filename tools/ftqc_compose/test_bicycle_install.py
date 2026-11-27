import copy
import unittest
from importlib.util import find_spec
from .test_insitu import study
from .bicycle_install import install
from ftqc_gadgets.bicycle_ports import generators


class InstalledTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spec = dict(kind="bb144-rm15", source=study(reject="1/4"), latency=1)
        cls.proof = install(cls.spec)

    def test_all_masks_and_transfer(self):
        p = self.proof
        self.assertEqual(p["ideal"]["input_masks"], 32768)
        self.assertEqual(p["ideal"]["accept_counts"], [1024, 1024])
        self.assertEqual(p["ideal"]["clean_output_phases"], [0, 1])
        self.assertEqual(p["factory"]["conditional_output_Z"], "1/1000")
        self.assertEqual(len(p["ideal"]["branches"]), 16)
        self.assertEqual(len(generators()), 12)
        for op in p["steps"]:
            if op["controls"]:
                self.assertTrue(
                    all(
                        any(
                            prior["record"] == k and prior["end"] < op["start"]
                            for prior in p["steps"]
                        )
                        for k in op["controls"]
                    )
                )

    def test_mutated_installation(self):
        for label, key, value in [
            ("r0.frame", "arms", []),
            ("r1.pivot-plus", "start", 0),
            ("r0.source-pivot-ZZ", "start", 0),
        ]:
            steps = copy.deepcopy(self.proof["steps"])
            next(op for op in steps if op["label"] == label)[key] = value
            with self.assertRaisesRegex(ValueError, "timing"):
                install(dict(self.spec, schedule=steps))
        with self.assertRaisesRegex(ValueError, "coupler"):
            install(dict(self.spec, couplers=self.proof["installed_couplers"][1:]))

    def test_physical_probe(self):
        if find_spec("stim") is None:
            self.skipTest("optional physical tableau probe")
        from .bicycle_probe import physical_probe

        for seed in (1, 17):
            self.assertTrue(physical_probe(self.spec, seed)["spectator_Bell"])

    def test_local_faults_and_cost(self):
        from fractions import Fraction
        from .bicycle_price import price

        inp = dict(kind="bb144-rm15-price", installation=self.spec, faults=[])
        self.assertEqual(price(inp)["conditional_output_Z"], "1/1000")
        inp["faults"] = [
            dict(loc, p="1/8") for loc in self.proof["fault_locations"][:3]
        ]
        r = price(inp)
        q = Fraction(1, 1000)
        a = Fraction(3, 4)
        # Direct enumeration: among eight subsets only empty and the full
        # triple survive the four terminal checks, independently of the carrier triple.
        import itertools

        accept = bad = 0
        for bits in itertools.product(range(2), repeat=3):
            syn = 0
            for i, b in enumerate(bits):
                if b:
                    syn ^= i + 1
            if syn:
                continue
            mass = a * Fraction(7 ** (3 - sum(bits)), 512)
            accept += mass
            bad += mass * (1 - q if sum(bits) % 2 else q)
        self.assertEqual(
            (Fraction(r["final_accept"]), Fraction(r["accepted_wrong"])), (accept, bad)
        )
        self.assertEqual([w["final_check_mask"] for w in r["witnesses"]], [1, 2, 3])
        cost = r["inclusive"]
        prep = cost["source_reject_ticks"]
        total = cost["full_attempt_ticks"]
        self.assertEqual(
            Fraction(cost["expected_qubit_cycles_per_output"]),
            cost["peak_reserved_qubits"] * (prep + a * (total - prep)) / accept,
        )
