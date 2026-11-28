import unittest
from . import staged, staged_switch


class StagedTests(unittest.TestCase):
    def test_switches_and_phase(self):
        for op in ("raise", "rm15-phase", "lower"):
            receipt = staged.analyze(staged_switch.example(op))
            self.assertEqual(receipt["result"]["status"], "proved")
            staged.verify(receipt)

    def test_sign_and_record(self):
        for op, key, value in (("rm15-phase", "exponents", [1]*15),
                               ("lower", "corrections", [0]*8),
                               ("raise", "corrections", [0]*8)):
            spec = staged_switch.example(op)
            spec[key] = value
            self.assertEqual(staged.analyze(spec)["result"]["status"], "counterexample")


    def test_ccz(self):
        from .staged_pieceable import example
        spec = example()
        result = staged.analyze(spec)["result"]
        self.assertEqual(result["status"], "proved", result.get("witness"))
        self.assertEqual(result["ideal"]["piece_lengths"], [9,6,6,6])
        self.assertGreater(result["faults"]["maximum_residual_cz_terms"], 0)
        spec["pieces"][2].pop()
        self.assertEqual(staged.analyze(spec)["result"]["status"], "counterexample")

    def test_conjugation(self):
        from ftqc_faults.staged import conjugate_ccz, z_expansion
        for x in range(8):
            edges = set()
            z = conjugate_ccz(x,0,edges,[0,1,2])
            if x == 1:
                self.assertEqual(edges, {(1,2)})
            for v in range(8):
                value = sum(k*(-1)**((v&m).bit_count()%2)
                            for m,k in z_expansion(z,edges).items())
                sign = (v==7) ^ ((v^x)==7) ^ (x==7)
                self.assertEqual(value, (-1)**sign * 2**len(edges))
