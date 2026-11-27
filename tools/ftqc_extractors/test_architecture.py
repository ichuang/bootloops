import unittest
from copy import deepcopy
from .__main__ import receipt, verify
from .algebra import distance_below

SMALL = dict(schema=1, kind="block_reading", memory=dict(n=4,hx=[15],hz=[15]),
             blocks=2, pattern=[3], axis="Z", target=2)


class BlockTests(unittest.TestCase):
    def test_reading(self):
        for axis in ("X", "Z"):
            s = dict(SMALL, axis=axis)
            r = receipt(s)
            verify(s,r)
            out = r["result"]
            self.assertEqual((len(out["measured"]),out["survivors"]), (2,2))
            self.assertIsNone(distance_below((out["width"],out["hx"],out["hz"]),2))
            self.assertEqual(out["split_basis"], [[1,1]])
            with self.assertRaises(ValueError):
                receipt(dict(s,expected=[out["measured"][0]^out["measured"][1]]))
            bad = deepcopy(r)
            bad["result"]["metachecks"] = []
            with self.assertRaises(ValueError):
                verify(s,bad)

    def test_patterns(self):
        self.assertEqual(receipt(dict(SMALL,blocks=3,pattern=[3,6]))["result"]["survivors"],2)
        with self.assertRaises(ValueError):
            receipt(dict(SMALL,pattern=[3,3]))

