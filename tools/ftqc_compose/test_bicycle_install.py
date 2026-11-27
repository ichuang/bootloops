import copy
import unittest
from .test_insitu import study
from .bicycle_install import install
from ftqc_gadgets.bicycle_ports import generators, route


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
