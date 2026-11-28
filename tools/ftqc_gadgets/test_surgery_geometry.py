import unittest
from ftqc_gadgets.surgery_design import build
from ftqc_gadgets.surgery import Instrument
from ftqc_gadgets.surgery_geometry import expansion, restricted_cone, skip_tree
from ftqc_gadgets.parsimonious import cellulate


class GeometryTests(unittest.TestCase):
    def test_cone_substitution(self):
        base = dict(
            n=2,
            hx=[],
            hz=[],
            sites=[0, 1],
            edges=[[0, 1]],
            targets=["+ZZ"],
            axis="Z",
            rounds=1,
        )
        for method in ("gauging", "parsimonious", "adapter"):
            out = build(dict(base, method=method, left=[], right=[]))
            Instrument(out["spec"])
            self.assertEqual(out["geometry"]["cone"]["logical_qubits"], 1)

    def test_sparse_cone(self):
        for n in (3, 4, 7):
            cycle = [(i, (i + 1) % n) for i in range(n)]
            g = cellulate(n, cycle)
            self.assertLessEqual(g["max_face"], 5)
            for i, (a, b) in enumerate(cycle):
                self.assertEqual(
                    set(g["edges"][g["original_edges"][i]]),
                    {g["ports"][a], g["ports"][b]},
                )

    def test_pcpp_proof_coordinate(self):
        # The fourth coordinate is a parity proof bit. Its value is quantified.
        p, a = [1, 2, 4, 0], [11, 12]
        e = expansion(p, a, 3, 1000)
        self.assertTrue(e["unit_bound"])
        c = restricted_cone(3, [], [], p, a)
        self.assertEqual(len(c["measured"]), 2)
        with self.assertRaises(ValueError):
            restricted_cone(3, [], [], p, [11, 12, 7], [0])
        self.assertEqual(expansion(p, a, 3, 1)["status"], "inconclusive")

    def test_skiptree(self):
        labels, paths = skip_tree(5, [(0, i) for i in range(1, 5)])
        self.assertEqual(set(labels), set(range(5)))
        self.assertTrue(all(c.bit_count() <= 3 for c in paths))

    def test_physical_surface_seam(self):
        from copy import deepcopy
        from ftqc_extractors.surgery_bench import surface_setting
        from ftqc_extractors.installed_surgery import compile_schedule
        from ftqc_extractors.surgery_layout import transport
        from .surgery_faults import pair_search

        for axis in "ZX":
            setting = surface_setting(3, 3, axis)
            spec, peak = compile_schedule(setting)
            layout = setting["transport"]
            transport(layout["initial"], layout["layers"], spec["ticks"], 1)
            instrument = Instrument(spec)
            self.assertEqual(peak, 35)
            self.assertEqual(pair_search(instrument)["status"], "excluded")
            self.assertEqual(instrument.audit(1, 100000)["status"], "proved")
            wrong = deepcopy(spec)
            # Invert an actual split-frame bit: after the conditional Z,
            # apply Z unconditionally on its whole support (before the X-basis undo).
            index, row = next((i, r) for i, r in enumerate(wrong["ticks"])
                              if any(op[0] == "Z_IF" for op in r))
            phase = [["S", op[1]] for op in row if op[0] == "Z_IF"]
            wrong["ticks"][index+1:index+1] = [phase, phase]
            with self.assertRaises(ValueError):
                Instrument(wrong)
            broken = deepcopy(setting)
            broken["couplers"].pop()
            with self.assertRaises(ValueError):
                compile_schedule(broken)

    def test_surface_patch_bridges(self):
        from .surgery_design import build
        from .surgery import Instrument
        from .surgery_merge import word

        # Two distance-two rotated patches, each with one encoded qubit.
        base = dict(
            n=8,
            hx=[15, 240],
            hz=[3, 12, 48, 192],
            axis="Z",
            rounds=2,
            sites=[0, 2, 4, 6],
            targets=[word(8, 85, "Z")],
            distance_target=2,
        )
        for shape in [
            dict(method="enhanced"),
            dict(
                method="bridge",
                left_size=2,
                left=[[0, 1]],
                right=[[0, 1]],
                links=[[0, 0], [1, 1]],
            ),
        ]:
            out = build(dict(base, **shape))
            Instrument(out["spec"])
            self.assertEqual(out["geometry"]["static"]["lower"], 2)
            self.assertEqual(out["geometry"]["cone"]["logical_qubits"], 1)
