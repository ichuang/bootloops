import unittest

from .contracts import BIT, REP, SCOPE, contract_result


def spec(program):
    return {
        "kind": "contract",
        "scope": SCOPE,
        "program": program,
        "require_friendly": False,
        "require_closed": True,
    }


def premise(name, source, target):
    return {
        "op": "hypothesis",
        "id": name,
        "inputs": [source],
        "outputs": [target],
        "ideal": "declared",
        "faults": {"op": "family", "sites": [name], "minimal": [[name]]},
        "filters": None,
    }


class Contracts(unittest.TestCase):
    def test_real_compositions(self):
        prepare = {"op": "prepare0", "sites": list("abc")}
        one = {"op": "sequence", "left": prepare, "right": {"op": "majority"}}
        two = {"op": "parallel", "left": one, "right": {"op": "bit-wire"}}
        r = contract_result(dict(spec(two), require_friendly=True))
        self.assertEqual(
            (r["status"], r["inputs"], r["outputs"]), ("proved", [BIT], [BIT, BIT])
        )
        self.assertEqual(r["filters"], ["identity"])
        self.assertEqual(r["bad_faults"]["enumerator"], {"2": 3})
        with self.assertRaisesRegex(ValueError, "shared physical"):
            contract_result(spec({"op": "parallel", "left": one, "right": one}))

    def test_conditional_and_join_direction(self):
        tight = dict(REP, bad=[[w] for w in REP["wires"]])
        left = premise("a", REP, tight)
        right = premise("b", REP, BIT)
        s = spec({"op": "sequence", "left": left, "right": right})
        with self.assertRaisesRegex(ValueError, "unproved premises"):
            contract_result(s)
        s["require_closed"] = False
        self.assertEqual(contract_result(s)["status"], "conditional")
        with self.assertRaisesRegex(ValueError, "friendly"):
            contract_result(dict(s, require_friendly=True))
        left["outputs"] = [REP]
        right["inputs"] = [tight]
        with self.assertRaisesRegex(ValueError, "residual"):
            contract_result(s)
        right["inputs"] = [dict(REP, decoder="another")]
        with self.assertRaisesRegex(ValueError, "decoder"):
            contract_result(s)
        right["inputs"] = [dict(REP, wires=["2", "1", "0"])]
        with self.assertRaisesRegex(ValueError, "wires"):
            contract_result(s)
