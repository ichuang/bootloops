import copy
import unittest
from .__main__ import receipt, verify
from .bicycle import mapping
from .algebra import parity


class BicycleMapTests(unittest.TestCase):
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
