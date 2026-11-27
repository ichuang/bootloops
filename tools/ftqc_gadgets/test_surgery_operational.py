import unittest
from importlib.util import find_spec
from .surgery_operational import conditional_law,sampled
from .test_surgery import edge

class OperationalTests(unittest.TestCase):
    def test_exclusive_modes(self):
        groups=[[("X",(1,1,0)),("Z",(1,0,1))]]*2
        law=conditional_law(groups,"1/2",20)
        self.assertEqual(law["acceptance"],"1/2")
        self.assertEqual(law["conditional"],[dict(readout=0,residual=0,p="3/4"),dict(readout=1,residual=1,p="1/4")])
        self.assertEqual(conditional_law(groups,"1/2",1)["status"],"inconclusive")

    @unittest.skipUnless(find_spec("stim") and find_spec("numpy"),"optional sampling backend")
    def test_decoding(self):
        result=sampled(edge(),"1/100",40,55,"bp_osd")
        self.assertEqual(result["decoder"],"bp_osd")
        self.assertEqual(result["shots"],40)


    @unittest.skipUnless(find_spec("stim"),"independent Stim fault injection")
    def test_physical_signatures(self):
        from .surgery import Instrument,centralizer
        from .surgery_faults import columns
        from .surgery_stim import compile_probe
        from ftqc_codes.gf2 import dot
        s=edge()
        inst=Instrument(s)
        n=inst.n
        for location,group in enumerate(columns(inst)[0]):
            for label,(d,b,res) in group:
                probe=compile_probe(inst.c,s,inst.initial,injected=(location,label))
                det,obs=probe.compile_detector_sampler().sample(1,separate_observables=True,bit_packed=True)
                self.assertEqual(d,int.from_bytes(det[0].tobytes(),"little"))
                def anticommutes(g):
                    return dot(res,g.z)^dot(res>>n,g.x)
                expected=[0]*len(inst.checks)
                expected.extend(anticommutes(g) for g in centralizer(inst.checks+inst.targets,n))
                for j,t in enumerate(inst.targets):
                    bit=(b>>j)&1
                    expected.extend([bit^anticommutes(t),bit])
                self.assertEqual(sum(v<<j for j,v in enumerate(expected)),int.from_bytes(obs[0].tobytes(),"little"))
