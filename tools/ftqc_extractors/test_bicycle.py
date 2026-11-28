import copy
import unittest
from .__main__ import receipt, verify
from .algebra import parity


class BicycleMapTests(unittest.TestCase):
    def test_compiler_shift_group(self):
        # bicycle_common::AutomorphismData multiplies exponent pairs modulo six.
        # Compare on every signed logical generator in BootLoops' own basis.
        from ftqc_codes import action_domain, automorphisms as aut
        from ftqc_codes.orbit_minweight import infer_default_automorphisms
        c = action_domain.load("bb144")
        paired = action_domain.default_basis(c)
        actions = [aut.build(c, paired, [{"perm": p}]) for p in infer_default_automorphisms(c)]
        identity = tuple((1 << j, 0, 0) for j in range(12)) + tuple((0, 1 << j, 0) for j in range(12))
        def push(a, frame):
            return tuple(aut.logical_image(a, p) for p in frame)
        seen = set()
        row = identity
        for _ in range(6):
            frame = row
            for _ in range(6):
                seen.add(frame)
                self.assertEqual(push(actions[0], push(actions[1], frame)),
                                 push(actions[1], push(actions[0], frame)))
                frame = push(actions[1], frame)
            self.assertEqual(frame, row)
            row = push(actions[0], row)
        self.assertEqual(row, identity)
        self.assertEqual(len(seen), 36)

    def test_native_settings(self):
        s = dict(schema=1, kind="bicycle_rm15_map")
        r = receipt(s)
        verify(s, r)
        m = r["result"]
        self.assertEqual(len(m["rotations"]), 15)
        self.assertEqual(m["orbit_size"], 36)
        self.assertTrue(any(op["mask_zero_slots"] for op in m["rotations"]))
        for op in m["rotations"]:
            image = sum(
                parity(x, op["physical_z"]) << i for i, (x, z) in enumerate(m["basis"])
            )
            self.assertEqual(
                image, (op["parity"] | 32) ^ sum(1 << i for i in op["mask_zero_slots"])
            )
            self.assertLess(image, 1 << 7)
        for key, value in (("shift", [11, 5]), ("sign", 1), ("mask_zero_slots", [7])):
            other = copy.deepcopy(m)
            other["rotations"][0][key] = value
            with self.assertRaisesRegex(ValueError, "map"):
                receipt(dict(s, candidate=other))
