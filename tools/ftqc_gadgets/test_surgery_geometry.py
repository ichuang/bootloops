import unittest
from ftqc_gadgets.surgery_design import build
from ftqc_gadgets.surgery import Instrument
from ftqc_gadgets.surgery_geometry import expansion,restricted_cone,skip_tree
from ftqc_gadgets.parsimonious import cellulate


class GeometryTests(unittest.TestCase):
    def test_cone_substitution(self):
        base=dict(n=2,hx=[],hz=[],sites=[0,1],edges=[[0,1]],targets=["+ZZ"],axis="Z",rounds=1)
        for method in ("gauging","parsimonious","adapter"):
            out=build(dict(base,method=method,left=[],right=[]))
            Instrument(out["spec"])
            self.assertEqual(out["geometry"]["cone"]["logical_qubits"],1)
    def test_sparse_cone(self):
        for n in (3,4,7):
            cycle=[(i,(i+1)%n) for i in range(n)]
            g=cellulate(n,cycle)
            self.assertLessEqual(g["max_face"],5)
            for i,(a,b) in enumerate(cycle):
                self.assertEqual(set(g["edges"][g["original_edges"][i]]),{g["ports"][a],g["ports"][b]})
    def test_pcpp_proof_coordinate(self):
        # The fourth coordinate is a parity proof bit. Its value is quantified.
        p,a=[1,2,4,0],[11,12]
        e=expansion(p,a,3,1000)
        self.assertTrue(e["unit_bound"])
        c=restricted_cone(3,[],[],p,a)
        self.assertEqual(len(c["measured"]),2)
        with self.assertRaises(ValueError):
            restricted_cone(3,[],[],p,[11,12,7],[0])
        self.assertEqual(expansion(p,a,3,1)["status"],"inconclusive")
    def test_skiptree(self):
        labels,paths=skip_tree(5,[(0,i) for i in range(1,5)])
        self.assertEqual(set(labels),set(range(5)))
        self.assertTrue(all(c.bit_count()<=3 for c in paths))


    def test_surface_patch_bridges(self):
        from .surgery_design import build
        from .surgery import Instrument
        from .surgery_merge import word
        # Two distance-two rotated patches, each with one encoded qubit.
        base=dict(n=8,hx=[15,240],hz=[3,12,48,192],axis="Z",rounds=2,
                  sites=[0,2,4,6],targets=[word(8,85,"Z")],distance_target=2)
        for shape in [dict(method="enhanced"),dict(method="bridge",left_size=2,
                      left=[[0,1]],right=[[0,1]],links=[[0,0],[1,1]])]:
            out=build(dict(base,**shape))
            Instrument(out["spec"])
            self.assertEqual(out["geometry"]["static"]["lower"],2)
            self.assertEqual(out["geometry"]["cone"]["logical_qubits"],1)
