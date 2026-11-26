import unittest

from .contracts import BIT, REP, SCOPE, contract_result


def spec(program):
    return dict(
        kind="contract",
        scope=SCOPE,
        program=program,
        require_friendly=False,
        require_closed=True,
    )


def premise(name, source, target):
    return dict(
        op="hypothesis",
        id=name,
        inputs=[source],
        outputs=[target],
        ideal="declared",
        faults=dict(op="family", sites=[name], minimal=[[name]]),
        filters=None,
    )


class Contracts(unittest.TestCase):
    def test_real_compositions(self):
        prepare = dict(op="prepare0", sites=list("abc"))
        one = dict(op="sequence", left=prepare, right=dict(op="majority"))
        two = dict(op="parallel", left=one, right=dict(op="bit-wire"))
        r = contract_result(dict(spec(two), require_friendly=True))
        self.assertEqual(
            (r["status"], r["inputs"], r["outputs"]), ("proved", [BIT], [BIT, BIT])
        )
        self.assertEqual(r["filters"], ["identity"])
        self.assertEqual(r["bad_faults"]["enumerator"], {"2": 3})
        with self.assertRaisesRegex(ValueError, "shared physical"):
            contract_result(spec(dict(op="parallel", left=one, right=one)))

    def test_conditional_and_join_direction(self):
        tight = dict(REP, bad=[[w] for w in REP["wires"]])
        left = premise("a", REP, tight)
        right = premise("b", REP, BIT)
        s = spec(dict(op="sequence", left=left, right=right))
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
