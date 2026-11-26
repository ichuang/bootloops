import copy
import unittest
from itertools import combinations

from .contracts import BIT, REP
from .level import level_result


def example():
    vertices = []
    for i in range(3):
        k = str(i)
        vertices.extend(
            [
                {
                    "id": "p" + k,
                    "inputs": [],
                    "outputs": ["block" + k],
                    "gadget": {
                        "op": "prepare0",
                        "sites": [f"{k}/{j}" for j in range(3)],
                    },
                },
                {
                    "id": "d" + k,
                    "inputs": ["block" + k],
                    "outputs": ["bit" + k],
                    "gadget": {"op": "majority", "site": "decoder" + k},
                },
            ]
        )
    vertices.append(
        {
            "id": "v",
            "inputs": ["bit0", "bit1", "bit2"],
            "outputs": ["answer"],
            "gadget": {"op": "vote3", "site": "vote"},
        }
    )
    terms = list(map(list, combinations(["p0", "p1", "p2"], 2))) + [
        [k] for k in ("d0", "d1", "d2", "v")
    ]
    return {
        "kind": "level",
        "inputs": {},
        "outputs": ["answer"],
        "physical_sites": [f"{i}/{j}" for i in range(3) for j in range(3)]
        + ["decoder" + str(i) for i in range(3)]
        + ["vote"],
        "vertices": vertices,
        "outer": {
            "op": "family",
            "sites": [v["id"] for v in vertices],
            "minimal": terms,
        },
        "output_check": "classical",
        "require_closed": True,
    }


class LevelTests(unittest.TestCase):
    def test_majority(self):
        r = level_result(example())
        self.assertEqual(r["bad_faults"]["enumerator"], {"4": 27})
        self.assertEqual(r["output_correctness"]["assignments"], 128)

    def test_hypotheses(self):
        h = {
            "op": "hypothesis",
            "id": "boundary theorem",
            "inputs": [BIT],
            "outputs": [BIT],
            "ideal": "identity-bit",
            "faults": {"op": "family", "sites": ["gate"], "minimal": []},
            "filters": ["filter"],
        }
        s = {
            "kind": "level",
            "inputs": {"x": BIT},
            "outputs": ["y"],
            "physical_sites": ["gate"],
            "vertices": [
                {
                    "id": "u",
                    "inputs": ["x"],
                    "outputs": ["y"],
                    "gadget": {"op": "assume-4.8", "contract": h},
                }
            ],
            "outer": {"op": "family", "sites": ["u"], "minimal": [["u"]]},
            "output_check": "none",
            "require_closed": False,
        }
        self.assertEqual(level_result(s)["status"], "conditional")
        self.assertIsNone(level_result(s)["output_correctness"])
        for key, value, message in [
            ("filters", None, "friendly"),
            ("faults", dict(h["faults"], minimal=[["gate"]]), "output-adjacent"),
        ]:
            broken = copy.deepcopy(s)
            broken["vertices"][0]["gadget"]["contract"][key] = value
            with self.assertRaisesRegex(ValueError, message):
                level_result(broken)
        with self.assertRaisesRegex(ValueError, "undischarged"):
            level_result(dict(s, require_closed=True))
        with self.assertRaisesRegex(ValueError, "trivial"):
            level_result(dict(s, inputs={"x": REP}))

    def test_negative_controls(self):
        s = example()
        with self.assertRaisesRegex(ValueError, "partition"):
            level_result(dict(s, physical_sites=[]))
        s["outer"]["minimal"].remove(["v"])
        with self.assertRaisesRegex(ValueError, "misses a harmful"):
            level_result(s)
        s = example()
        s["vertices"][2]["gadget"]["sites"][0] = "0/0"
        with self.assertRaisesRegex(ValueError, "shared physical"):
            level_result(s)
        s = example()
        s["vertices"][-1]["inputs"][0] = "block0"
        with self.assertRaisesRegex(ValueError, "bundle"):
            level_result(s)
