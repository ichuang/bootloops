import unittest
from copy import deepcopy
from .__main__ import receipt, verify
from .algebra import distance_below

SMALL = dict(
    schema=1,
    kind="block_reading",
    memory=dict(n=4, hx=[15], hz=[15]),
    blocks=2,
    pattern=[3],
    axis="Z",
    target=2,
)


class BlockTests(unittest.TestCase):
    def test_reading(self):
        for axis in ("X", "Z"):
            s = dict(SMALL, axis=axis)
            r = receipt(s)
            verify(s, r)
            out = r["result"]
            self.assertEqual((len(out["measured"]), out["survivors"]), (2, 2))
            self.assertIsNone(distance_below((out["width"], out["hx"], out["hz"]), 2))
            self.assertEqual(out["split_basis"], [[1, 1]])
            with self.assertRaises(ValueError):
                receipt(dict(s, expected=[out["measured"][0] ^ out["measured"][1]]))
            bad = deepcopy(r)
            bad["result"]["metachecks"] = []
            with self.assertRaises(ValueError):
                verify(s, bad)

    def test_patterns(self):
        self.assertEqual(
            receipt(dict(SMALL, blocks=3, pattern=[3, 6]))["result"]["survivors"], 2
        )
        with self.assertRaises(ValueError):
            receipt(dict(SMALL, pattern=[3, 3]))


class GraphTests(unittest.TestCase):
    def test_signed_and_routing(self):
        import json
        from pathlib import Path

        s = json.loads(
            (Path(__file__).parent / "examples/fixed-graph.json").read_text()
        )
        for word in ("+XXII", "-YYII", "+ZZII"):
            r = receipt(dict(s, logical=word))
            verify(dict(s, logical=word), r)
            self.assertEqual(r["result"]["readout_sign"], int(word[0] == "-"))
        self.assertEqual(
            receipt(dict(s, routes=[1, 7]))["result"]["obligation"],
            "even-subset routing",
        )
        with self.assertRaises(ValueError):
            receipt(dict(s, checks=["-IIII", "+ZZZZ"]))


class SequenceTests(unittest.TestCase):
    def test_sequence_and_boundaries(self):
        import json
        from pathlib import Path

        s = json.loads(
            (Path(__file__).parent / "examples/three-block-sequence.json").read_text()
        )
        r = receipt(s)
        self.assertEqual(r["result"]["distance_lower"], 2)
        verify(s, r)
        for key, value in [
            ("padding", [0, 2]),
            ("operations", [s["operations"][0]] * 2),
        ]:
            with self.assertRaises(ValueError):
                receipt(dict(s, **{key: value}))
        del s["operations"][0]["metachecks"][:]
        with self.assertRaises(ValueError):
            receipt(s)
        self.assertEqual(
            receipt(dict(SMALL, target=3))["result"]["status"], "counterexample"
        )


class GrossReadingTests(unittest.TestCase):
    def test_both_axes(self):
        from ftqc_codes.codes import CSSCode, load
        from ftqc_faults.audit.record_model import RecordModel
        from .reading_audit import probe
        from .block_reading import reading

        bb = load("bb144")
        for axis in ("X", "Z"):
            mem = bb if axis == "Z" else CSSCode("dual", 144, bb.hz, bb.hx)
            r = reading((144, mem.hx, mem.hz), 2, [3], "Z")
            self.assertEqual(
                (
                    len(r["measured"]),
                    r["survivors"],
                    r["added_data"] + r["added_check_meters"],
                ),
                (12, 12, 216),
            )
            circ, *_ = probe(mem, axis, "X")
            self.assertTrue(RecordModel(circ).locations)


class ExpansionTests(unittest.TestCase):
    def test_expansion_and_soundness(self):
        import json
        from pathlib import Path

        path = Path(__file__).parent / "examples"
        e = json.loads((path / "modular-expansion.json").read_text())
        r = receipt(e)
        verify(e, r)
        self.assertEqual(r["result"]["original"], "1/2")
        self.assertTrue(r["result"]["meets_unit_bound"])
        s = json.loads((path / "intermediate-reading.json").read_text())
        r = receipt(s)
        verify(s, r)
        self.assertEqual(r["result"]["quantum_soundness"], "2")
        self.assertEqual(
            receipt(dict(s, rounds=[1]))["result"]["status"], "counterexample"
        )


class ZXTests(unittest.TestCase):
    def test_homomorphic(self):
        import json
        from pathlib import Path

        s = json.loads(
            (Path(__file__).parent / "examples/homomorphic-surgery.json").read_text()
        )
        r = receipt(s)
        verify(s, r)
        self.assertEqual(r["result"]["split_Z_frame_columns"], [1])
        self.assertEqual(
            set(r["result"]["physical_open_wire_checks"]),
            {"homomorphic", "cone surgery"},
        )


class CompactedDistanceTests(unittest.TestCase):
    def test_low_weight_auxiliary_output(self):
        # Four disjoint logical images, with only the product of all four trivial.
        # The new edge qubit alone supports an undetectable logical Z.
        op = dict(
            aux=SMALL["memory"],
            cycles=[],
            f0=[0],
            f1=[3, 12, 48, 192],
            f2=[1],
            metachecks=[31],
        )
        s = dict(
            schema=1,
            kind="surgery_sequence",
            memory=dict(n=8, hx=[3, 12, 48, 192], hz=[255]),
            operations=[op],
            padding=[2, 2],
            max_degree=8,
            target=2,
        )
        r = receipt(s)
        verify(s, r)
        self.assertEqual(r["result"]["obligation"], "compacted code distance")
        self.assertEqual(r["result"]["witness"], dict(axis="Z", qubits=[8]))


class SpectatorTests(unittest.TestCase):
    def test_one_sector_and_spectators(self):
        from .surgery_sequence import cosystole_witness

        self.assertIsNone(cosystole_witness((3, [], [3, 6]), 3, 100))
        self.assertIsNotNone(distance_below((3, [], [3, 6]), 3))
        s = receipt(SMALL)["result"]["spectator_basis"]
        self.assertEqual(len(s["same_axis"]), 2)
        self.assertTrue(all(s["pairing"]))


class PartialReadingTests(unittest.TestCase):
    def test_spectator_factors(self):
        import json
        from pathlib import Path
        from .algebra import css_dimension

        s = json.loads(
            (Path(__file__).parent / "examples/partial-block-sequence.json").read_text()
        )
        r = receipt(s)
        verify(s, r)
        result = r["result"]
        self.assertEqual(result["distance_lower"], 2)
        self.assertEqual(len(result["measured"][0]), 2)
        self.assertTrue(all(v & 0xF0F0 == 0 for v in result["measured"][0]))
        code = result["compacted"]
        self.assertEqual(css_dimension((code["n"], code["hx"], code["hz"])), 6)
