import unittest
from . import staged, staged_switch


class StagedTests(unittest.TestCase):
    def test_switches_and_phase(self):
        for op in ("raise", "rm15-phase", "lower"):
            receipt = staged.analyze(staged_switch.example(op))
            self.assertEqual(receipt["result"]["status"], "proved")
            staged.verify(receipt)

    def test_sign_and_record(self):
        for op, key, value in (
            ("rm15-phase", "exponents", [1] * 15),
            ("lower", "corrections", [0] * 8),
            ("raise", "corrections", [0] * 8),
        ):
            spec = staged_switch.example(op)
            spec[key] = value
            self.assertEqual(staged.analyze(spec)["result"]["status"], "counterexample")

    def test_ccz(self):
        from .staged_pieceable import example

        spec = example()
        result = staged.analyze(spec)["result"]
        self.assertEqual(result["status"], "proved", result.get("witness"))
        self.assertEqual(result["ideal"]["piece_lengths"], [9, 6, 6, 6])
        self.assertGreater(result["faults"]["maximum_residual_cz_terms"], 0)
        spec["pieces"][2].pop()
        self.assertEqual(staged.analyze(spec)["result"]["status"], "counterexample")

    def test_conjugation(self):
        from ftqc_faults.staged import conjugate_ccz, z_expansion

        for x in range(8):
            edges = set()
            z = conjugate_ccz(x, 0, edges, [0, 1, 2])
            if x == 1:
                self.assertEqual(edges, {(1, 2)})
            for v in range(8):
                value = sum(
                    k * (-1) ** ((v & m).bit_count() % 2)
                    for m, k in z_expansion(z, edges).items()
                )
                sign = (v == 7) ^ ((v ^ x) == 7) ^ (x == 7)
                self.assertEqual(value, (-1) ** sign * 2 ** len(edges))

    def test_intermediate_record_witness(self):
        from .staged_pieceable import example

        spec = dict(example(), record_history=False)
        rec = staged.analyze(spec)
        self.assertEqual(rec["result"]["status"], "counterexample")
        self.assertIn("second", rec["result"]["witness"])
        staged.verify(rec)

    def test_physical_check_against_stim(self):
        try:
            import stim
        except ImportError:
            self.skipTest("Stim optional cat-leg cross-check")
        from ftqc_faults.staged import cat_faults

        for axis in ("X", "Z"):
            gate = "CX" if axis == "X" else "CZ"
            for mode, x, z, flip in cat_faults(3, 5, axis):
                name, q, *bits = mode
                error = stim.PauliString(4)

                def pauli(a, b):
                    return {(0, 0): 0, (1, 0): 1, (1, 1): 2, (0, 1): 3}[a, b]

                if name == "read-X":
                    error[3] = 3
                elif name == "idle":
                    error[q] = pauli(*bits)
                elif name == "cat-leg":
                    error[3] = pauli(*bits)
                    error = error.after(
                        stim.Tableau.from_named_gate(gate), targets=[3, q]
                    )
                else:
                    error[3] = pauli(*bits[:2])
                    error[q] = pauli(*bits[2:])
                self.assertEqual(x, sum((error[i] in (1, 2)) << i for i in range(3)))
                self.assertEqual(z, sum((error[i] in (2, 3)) << i for i in range(3)))
                self.assertEqual(flip, int(error[3] in (2, 3)))

    def test_supplied_recovery(self):
        spec = staged_switch.example("rm15-phase")
        rows = staged.analyze(spec)["result"]["faults"]["recovery"]
        spec["terminal_recovery"] = rows
        self.assertEqual(staged.analyze(spec)["result"]["status"], "proved")
        rows[0][-1] ^= 32767
        failed = staged.analyze(spec)
        self.assertEqual(failed["result"]["status"], "counterexample")
        self.assertIn("origin", failed["result"]["witness"])
        staged.verify(failed)
