import cmath
import json
import unittest
from pathlib import Path

from ftqc_gadgets.surgery import Instrument, receipt as gadget_receipt
from ftqc_gadgets.surgery_build import graph
from ftqc_magic.certify import build_receipt
from .__main__ import receipt, verify
from .surgery_products import lifted_product
from .test_pbc import apply, logical_action


def tiny(physical=True):
    s = json.loads((Path(__file__).parent / "examples/bridged-pbc.json").read_text())
    block = s["blocks"][0]
    n, hx, hz = lifted_product([[3]], [[3]], 2)
    assert (n, set(hx), set(hz)) == (4, {15}, {15})
    s.update(
        blocks=[block],
        bridges=[],
        magic_supplies=[0],
        target=1,
        max_degree=16,
        gates=[
            ["H", [0, 0]],
            ["T", [0, 0]],
            ["S", [0, 0]],
            ["S", [0, 0]],
            ["S", [0, 0]],
            ["H", [0, 0]],
        ],
    )
    if physical:
        settings = []
        links = set()
        compiled = receipt(s)["result"]
        for word in sorted({v["logical"] for v in compiled["settings"].values()}):
            size = sum(c != "I" for c in word[1:])
            g = graph(4, block["checks"], word, [(i, i + 1) for i in range(size - 1)])
            inst = Instrument(g)
            settings.append(
                dict(gadget=gadget_receipt(g, 0), placement=list(range(inst.c.n)))
            )
            links.update(
                tuple(sorted(op.qs))
                for _, _, op in inst.c.locations()
                if op.name == "CX"
            )
        s["surgery_installation"] = dict(
            idles=True, couplers=sorted(map(list, links)), settings=settings
        )
    factory = dict(
        kind="factory", catalog="bk_15_1_3_t", p_phys="1/1000", mode="postselect"
    )
    wires = list(range(30, 35))
    attempts = [
        dict(id="lane0/1", start=0, end=15, wires=wires, raw=30, accepted=False),
        dict(id="lane0/2", start=15, end=30, wires=wires, raw=30, accepted=True),
    ]
    s["magic_delivery"] = dict(
        factory=factory,
        receipt=build_receipt(factory),
        attempts=attempts,
        requests=[dict(attempt="lane0/2", block=0, slot=1, member=0, start=32)],
        latency=1,
        couplers=[[i, 30] for i in range(3)],
    )
    return s


class InstalledDeliveryTests(unittest.TestCase):
    def test_swap_and_signed_readout(self):
        s = tiny()
        proof = receipt(s)
        verify(s, proof)
        out = proof["result"]
        d = out["magic_delivery"]
        self.assertEqual(d["transfers"][0]["proof"]["unchanged_logicals"], [0])
        for branch, ledger in zip(out["branches"], d["branches"], strict=True):
            self.assertEqual(ledger["buffer_wire_ticks"], 2)
            self.assertEqual(ledger["steps"][0]["consume_tick"], 38)
            self.assertEqual(ledger["steps"][1]["start"], 38)
            # Contract both actual signed projectors with the delivered T and
            # evaluate the native final signed physical measurement on |0>.
            v = [1 / 2**0.5, 0, cmath.exp(1j * cmath.pi / 4) / 2**0.5, 0]
            raw = branch["outcomes"][0]
            for step, bit in zip(
                branch["program"][1:3], (raw & 1, raw >> 1), strict=True
            ):
                op = logical_action(
                    out["settings"][step["measurement"]], s["blocks"][0]
                )
                pv = apply(v, *op)
                v = [(a + (-1) ** bit * b) / 2 for a, b in zip(v, pv, strict=True)]
            last = branch["program"][-1]
            op = logical_action(out["settings"][last["measurement"]], s["blocks"][0])
            prob = sum(abs(a) ** 2 for a in v)
            self.assertAlmostEqual(prob, 0.25)
            self.assertAlmostEqual(
                sum(
                    a.conjugate() * b for a, b in zip(v, apply(v, *op), strict=True)
                ).real
                / prob,
                2**-0.5,
            )

    def test_rejected_early_reused_and_collision(self):
        for case in ("rejected", "early", "double", "collision"):
            s = tiny(physical=case == "collision")
            d = s["magic_delivery"]
            if case == "rejected":
                d["requests"][0]["attempt"] = "lane0/1"
            elif case == "early":
                d["requests"][0]["start"] = 29
            elif case == "double":
                s["gates"].append(["T", [0, 0]])
                d["requests"] *= 2
            else:
                a = dict(
                    d["attempts"][0],
                    id="other/1",
                    start=38,
                    end=53,
                    wires=[4, 40, 41, 42, 43],
                    raw=40,
                )
                d["attempts"].append(a)
            with self.assertRaises(ValueError):
                receipt(s)
