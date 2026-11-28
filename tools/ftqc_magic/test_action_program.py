import unittest

from ftqc_codes import action_domain, automorphisms
from .action_program import transport


class MagicActionTests(unittest.TestCase):
    def test_signed_rotation_axis(self):
        code = action_domain.load("bb144")
        paired = action_domain.default_basis(code)
        xword = paired["slots"][0]["x"]
        steps = [
            {"layer": [[gate, q] for q, gate in enumerate(xword[1:]) if gate != "I"]}
        ]
        rec = automorphisms.build(code, paired, steps)
        slots = [s["name"] for s in paired["slots"][:5]]
        result = transport(
            dict(kind="rm15-action", code="bb144", slots=slots, action=rec)
        )
        self.assertEqual(len(result["rotations"]), 15)
        for row in result["rotations"]:
            self.assertEqual(row["signed_axis"][0], "-" if row["parity"] & 1 else "+")
        self.assertEqual(result["output_z"][0], "-")
