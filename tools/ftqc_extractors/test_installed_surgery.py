import copy
import unittest
from ftqc_gadgets.surgery import receipt
from ftqc_gadgets.surgery_build import graph
from ftqc_gadgets.test_surgery import edge
from ftqc_extractors.installed_surgery import analyze, compile_schedule
from ftqc_compose.contracts import derive
from ftqc_extractors.__main__ import evaluate


def setting(g, pos, control=None):
    out = dict(gadget=receipt(g, 0), placement=pos)
    if control is not None:
        out["control"] = control
    return out


def workload():
    g = edge()
    # A cross-block parity joins a surviving HGP logical to an uncoded block.
    word = "+IIIIZZ" + "II" + "Z"
    bridge = graph(
        9,
        g["checks"] and [w + "I" for w in g["checks"]],
        word,
        [(0, 1), (1, 2)],
        rounds=2,
    )
    # X on the first HGP spectator times X on the external wire commutes
    # with the joint Z parity; the original measured Z1 commutes as well.
    follow = graph(9, bridge["checks"], "+IIIIXIXIX", [(0, 1), (1, 2)])
    parts = [
        [setting(g, list(range(8)) + [9, 10])],
        [setting(bridge, list(range(14)))],
        [setting(follow, list(range(14)), [0, 4])],
    ]
    couplers = set()
    for batch in parts:
        for part in batch:
            pos = part["placement"]
            for row in part["gadget"]["spec"]["ticks"]:
                for op in row:
                    if op[0] == "CX":
                        couplers.add(tuple(sorted(pos[i] for i in op[1:])))
    return dict(
        schema=1,
        kind="installed_surgery",
        n=9,
        checks=bridge["checks"],
        targets=[g["targets"][0] + "I", word, follow["targets"][0]],
        couplers=[list(p) for p in sorted(couplers)],
        batches=parts,
        idles=False,
        terminal_syndrome=True,
        target=2,
    )


class InstalledTests(unittest.TestCase):
    def test_controlled_handoff(self):
        s = workload()
        out = analyze(s, 1000000)
        self.assertEqual(out["gadget"]["public_rank"], 3)
        # Physical verification of the joint circuit runs even for ideal component receipts.
        self.assertGreater(out["gadget"]["faults"]["tested"], 0)
        s["target"] = 1
        proof = evaluate(s, "exact", 100000, 60)["gadget"]
        node = dict(
            op="surgery", instance="joint", wires=list("abcdefghi"), receipt=proof
        )
        self.assertEqual(len(derive(node)["outputs"]), 4)
        wrong = copy.deepcopy(proof)
        wrong["spec"]["ticks"].pop()
        with self.assertRaises(ValueError):
            derive(dict(node, receipt=wrong))

    def test_resource_and_bits(self):
        s = workload()
        for key in ("couplers", "parallel", "targets"):
            t = copy.deepcopy(s)
            if key == "couplers":
                t[key].pop()
            elif key == "parallel":
                t["batches"][0] += t["batches"].pop(1)
            else:
                t[key].reverse()
            with self.assertRaises(ValueError):
                compile_schedule(t)
