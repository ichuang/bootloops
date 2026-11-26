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


class NoisyCycleTests(unittest.TestCase):
    def test_corrected_cycle(self):
        from .gross_cycle import cycle, validate_cycle

        _, (_, hx, hz), _, _ = gross_z()
        body = cycle()
        # Drop the early gauge CNOT. A tick schedule still exists, but its
        # signed meter pullback must fail.
        body[1] = [op for op in body[1] if op[0] != "CX" or op[1] != 234]
        with self.assertRaisesRegex(ValueError, "incorrect ideal check"):
            validate_cycle(body, hx, hz)

    def test_record_columns_against_forward_ftqc_audit(self):
        from ftqc_audit.circuit import Circuit
        from ftqc_audit.model import Model
        from ftqc_audit.record_model import RecordModel

        circ = Circuit(
            [[("RX", 0), ("R", 1)], [("CX", 0, 1)], [("M", 0), ("M", 1)]],
            detectors=[],
            observables=[[0, 1]],
            outputs=[],
        )
        native, fast = Model(circ), RecordModel(circ)
        for site, group in enumerate(fast.groups):
            for label, sig in group.items():
                self.assertEqual(native.replay([(site, label)]), (0, sig, 0))
                self.assertEqual(fast.replay([(site, label)]), (0, sig))
        with self.assertRaisesRegex(ValueError, "random"):
            RecordModel(Circuit([[("RX", 0)], [("M", 0)]], detectors=[[0]], outputs=[]))


class CircuitDistanceTests(unittest.TestCase):
    def test_projected_exclusion(self):
        from itertools import combinations
        from random import Random
        from operator import xor
        from functools import reduce
        from ftqc_audit.projected import exclude

        rng = Random(2026)
        for _ in range(50):
            cols = sorted(set(rng.randrange(1, 64) for _ in range(9)))
            for w in range(5):
                sums = (
                    reduce(xor, subset, 0)
                    for r in range(1, w + 1)
                    for subset in combinations(cols, r)
                )
                harmful = any(v and not v & 15 for v in sums)
                ok, _, hit = exclude(cols, 4, w, 100000)
                self.assertEqual(ok, not harmful)
                if hit:
                    self.assertLessEqual(len(hit), w)
                    self.assertTrue(reduce(xor, hit))
                    self.assertEqual(reduce(xor, hit) & 15, 0)

    def test_short_rounds_and_tampering(self):
        from .__main__ import receipt, verify

        s = dict(schema=1, kind="gross_circuit", rounds=2, target=3)
        r = receipt(s)
        self.assertEqual((r["result"]["lower"], r["result"]["upper"]), (2, 2))
        self.assertEqual(r["result"]["status"], "counterexample")
        verify(s, r)
        r["result"]["witness"]["faults"][0][0] += 1
        with self.assertRaises(ValueError):
            verify(s, r)


class SpacetimeTests(unittest.TestCase):
    def test_domains(self):
        from ftqc_audit.spacetime import xor_feasible
        from itertools import product
        from random import Random

        rng = Random(155)
        for _ in range(90):
            groups = [
                set(rng.sample(range(16), rng.randint(1, 4)))
                for _ in range(rng.randint(1, 4))
            ]
            target = rng.randrange(16)
            self.assertEqual(
                xor_feasible(groups, target),
                any(reduce(xor, c) == target for c in product(*groups)),
            )

    def test_native_cover_and_zero_projection_fault(self):
        from ftqc_audit.spacetime import (
            support_enumerator,
            lift_bound,
            enumerate_supports,
        )
        from itertools import combinations
        from random import Random
        import shutil

        if not shutil.which("c++"):
            self.skipTest("C++17 compiler absent")
        rng = Random(12)
        with support_enumerator() as exe:
            for _ in range(15):
                cols = sorted(set(rng.randrange(1, 32) for _ in range(8)))
                patterns, _ = enumerate_supports(exe, cols, 4, 4, 1000000)
                actual = [
                    ids
                    for r in range(1, 5)
                    for ids in combinations(range(len(cols)), r)
                    if reduce(xor, (cols[i] for i in ids)) == 16
                ]
                # Minimum words are irreducible. At most one further site
                # cannot add a nonempty detector-silent/trivial subset.
                d = min(map(len, actual), default=5)
                self.assertEqual(min(map(len, patterns), default=5), d)
                if d >= 2:
                    self.assertEqual(
                        sorted(p for p in patterns if len(p) <= d + 1),
                        sorted(p for p in actual if len(p) <= d + 1),
                    )
            # The three spatial classes sum to a logical. An extra fault
            # with zero spatial projection is needed to cancel its time tag.
            tags = {1: {1}, 2: {2}, 7: {0}, 0: {0, 3}}
            self.assertEqual(lift_bound(exe, tags, 2, 3, 1000)["lower"], 4)
            self.assertEqual(
                lift_bound(exe, tags, 2, 4, 1000)["outcome"], "possible lift"
            )
