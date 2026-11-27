import unittest
from ftqc_gadgets.surgery import receipt
from ftqc_gadgets.surgery_design import build
from .surgery_bench import architecture
from .installed_surgery import compile_schedule
from .surgery_layout import transport
from ftqc_gadgets.surgery import Instrument


class LayoutTests(unittest.TestCase):
    def test_movements(self):
        layout=[[0,0],[3,0]]
        tick=[[["CX",0,1]]]
        step=dict(moves=[[1,[1,0]]],move_time=2,gate_time=1)
        self.assertEqual(transport(layout,[step],tick,1)["elapsed"],3)
        with self.assertRaisesRegex(ValueError,"collision"):
            transport(layout,[dict(step,moves=[[0,[3,0]],[1,[0,0]]])],tick,3)
        with self.assertRaisesRegex(ValueError,"radius"):
            transport(layout,[dict(step,moves=[])],tick,1)

    def test_parallel(self):
        targets=None
        for case in ("serial","branched","parallel-parsimonious"):
            s=architecture(case)
            spec,_=compile_schedule(s)
            Instrument(spec)
            if targets is not None:
                self.assertEqual(targets,spec["targets"])
            targets=spec["targets"]
        self.assertEqual(len(targets),3)

    def test_product_and_cyclic(self):
        s=architecture("product-connection")
        self.assertEqual(len(s["targets"]),2)
        Instrument(build(s)["spec"])
        Instrument(compile_schedule(architecture("cyclic-seed"))[0])

    def test_finite_pcpp(self):
        s=architecture("pcpp")
        g=build(dict(s,noise_locations="phenomenological"))
        self.assertTrue(g["geometry"]["relative_expansion"]["unit_bound"])
        self.assertEqual(receipt(g["spec"],1)["faults"]["lower"],2)
        self.assertEqual(receipt(build(s)["spec"],1)["faults"]["upper"],1)

