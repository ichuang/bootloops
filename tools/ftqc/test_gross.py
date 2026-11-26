import unittest
from functools import reduce
from operator import xor
from .algebra import bits, css_dimension, echelon, reduce_row
from .gross import cell, gross_z, logical_action, shifted_z


class GrossConstructionTests(unittest.TestCase):
    def test_tanner_coordinates_and_check_six(self):
        memory, (n, x, z), ports, incident = gross_z()
        for i in range(12):
            for j in range(6):
                row = sum(1 << cell(i + a, j + b) for a, b in [(3, 0), (0, 1), (0, 2)])
                row |= sum(
                    1 << (72 + cell(i + a, j + b)) for a, b in [(0, 3), (1, 0), (2, 0)]
                )
                self.assertEqual(memory.hx[cell(i, j)], row)
        # Figure 9 labels on the four edges of the purple square numbered 6.
        self.assertEqual({incident[q - 144] for q in bits(x[-1])}, {4, 10, 19, 20})
        self.assertEqual((n, len(ports)), (162, 12))
        self.assertEqual(reduce(xor, ports), shifted_z())
        self.assertEqual(css_dimension((n, x[:-1], z)), 12)
        self.assertEqual(css_dimension((n, x, z)), 11)

    def test_surviving_logicals(self):
        action = logical_action()
        self.assertEqual((action["memory_k"], action["merged_k"]), (12, 11))
        _, (_, _, z), _, _ = gross_z()
        self.assertTrue(reduce_row(shifted_z(0, 1), echelon(z)))


class CoverTests(unittest.TestCase):
    def test_cover_against_all_supports(self):
        from random import Random
        from . import css_cover
        from .algebra import kernel, parity

        rng = Random(144)
        for n in range(2, 9):
            for _ in range(6):
                checks = [rng.randrange(1 << n) for _ in range(n // 2)]
                stabs = kernel(checks, n)[::2]
                span = echelon(stabs)
                logicals = [
                    q
                    for q in range(1, 1 << n)
                    if not any(parity(q, h) for h in checks) and reduce_row(q, span)
                ]
                for bound in range(n + 1):
                    w, nodes, proof = css_cover.cover(n, checks, stabs, bound, 100000)
                    self.assertEqual(
                        bool(w), any(q.bit_count() <= bound for q in logicals)
                    )
                    if w:
                        self.assertIn(w, logicals)
                        self.assertLessEqual(w.bit_count(), bound)
                    else:
                        css_cover.verify_cover(
                            n, checks, stabs, bound, nodes, proof, nodes
                        )

    def test_false_packing_and_missing_branch(self):
        from . import css_cover

        # At root zero both unsatisfied checks have the same remaining site.
        for tokens in ([130, 0, 1], [130, 0, 0], [1], []):
            with self.assertRaises(ValueError):
                css_cover.verify_cover(
                    2, [3, 3], [], 2, 3, css_cover.encode(bytes(tokens)), 3
                )

    def test_gross_receipt(self):
        import json
        from pathlib import Path
        from .__main__ import receipt, verify
        from . import css_cover

        s = json.loads((Path(__file__).parent / "examples/gross-z.json").read_text())
        self.assertEqual(receipt(s, budget=4)["result"]["status"], "inconclusive")
        r = receipt(s)
        for code in r["result"]["codes"].values():
            self.assertEqual(
                [code["axes"][p]["distance"] for p in ("X", "Z")], [12, 12]
            )
        verify(s, r)
        axis = r["result"]["codes"]["memory"]["axes"]["X"]
        original = axis["witness"]
        axis["witness"] = [0]
        with self.assertRaisesRegex(ValueError, "logical"):
            verify(s, r)
        axis["witness"] = original
        axis["cover"] = css_cover.encode(b"")
        with self.assertRaisesRegex(ValueError, "incomplete"):
            verify(s, r)


class PortFrameTests(unittest.TestCase):
    def test_every_ideal_x_record(self):
        from .algebra import parity, transpose
        from .gross import measurement, port_frame

        _, _, ports, labels = gross_z()
        paths, cycles = port_frame(ports, 144, len(labels))
        edges = transpose([row >> 144 for row in ports], len(labels))
        for potential in range(0, 1 << len(ports), 2):
            record = sum(parity(potential, edge) << j for j, edge in enumerate(edges))
            self.assertFalse(any(parity(record, c) for c in cycles))
            self.assertEqual(
                sum(parity(record, p) << j for j, p in enumerate(paths)), potential
            )
        export = measurement()
        self.assertEqual(export["logical_z_record_parity"], list(range(12)))
        self.assertEqual(export["merged_code"]["n"], 162)

    def test_omitting_gauge_six_admits_a_light_logical(self):
        from .gross_distance import witness_mask

        _, (n, x, z), _, _ = gross_z()
        support = list(bits(x[-1]))
        self.assertEqual(len(support), 4)
        witness_mask(support, n, z, x[:-1])
        with self.assertRaisesRegex(ValueError, "not a nontrivial logical"):
            witness_mask(support, n, z, x)
