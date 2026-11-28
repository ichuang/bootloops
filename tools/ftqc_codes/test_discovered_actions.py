"""Replay externally discovered generators and the signed Table 4 GB fold."""

import json
from pathlib import Path
import unittest

from ftqc_codes import action_domain, automorphisms as aut, certify
from ftqc_extractors import action_orbit

EXAMPLES = Path(__file__).with_name("examples")


class ActionTests(unittest.TestCase):
    def test_gb_fold(self):
        code = action_domain.load("bb:9,1,x0+x1,x1+x2+x3+x5+x7+x8")
        paired = action_domain.default_basis(code)
        steps = json.loads((EXAMPLES / "gb18-fold.steps.json").read_text())
        with self.assertRaises(certify.Refusal):
            aut.build(code, paired, steps[:1])
        action = aut.build(code, paired, steps)
        self.assertEqual(aut.logical_image(action, (1, 0, 0)), (1, 2, 0))
        self.assertEqual(aut.logical_image(action, (2, 0, 0)), (2, 1, 1))

    def test_qldpc_hadamard_orbit(self):
        from ftqc_extractors.__main__ import receipt, verify

        code = action_domain.load("steane")
        action = json.loads((EXAMPLES / "steane.qldpc-h.action.json").read_text())
        aut.verify(code, action)
        spec = dict(
            schema=1,
            kind="action_orbit",
            code="steane",
            basis=action["basis"],
            actions=[action],
            seed="+Y",
            targets=["-Y", "+X"],
        )
        rec = receipt(spec)
        verify(spec, rec)
        self.assertEqual(rec["result"]["settings"][0]["word"], [0])
        self.assertFalse(rec["result"]["settings"][1]["reachable"])

    def test_discovered_noncss_and_signed_target(self):
        path = str(EXAMPLES / "five_qubit.stabilizer.json")
        code = action_domain.load(path)
        action = json.loads((EXAMPLES / "five_qubit.autqec.action.json").read_text())
        aut.verify(code, action)
        spec = dict(
            code=path,
            basis=action["basis"],
            actions=[action],
            seed="+X",
            targets=["+Y", "+Z"],
        )
        result = action_orbit.analyze(spec, 6)
        self.assertEqual(result["settings"][0]["word"], [0])
        identity = aut.from_spec(
            code, action["basis"], {"actions": [action], "word": [0, {"inverse": 0}]}
        )
        for p in ((1, 0, 0), (0, 1, 0), (1, 1, 1)):
            self.assertEqual(aut.logical_image(identity, p), p)
        self.assertFalse(result["settings"][1]["reachable"])
