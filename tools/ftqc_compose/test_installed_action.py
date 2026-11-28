import unittest
from ftqc_codes import action_domain, automorphisms
from ftqc_compose.__main__ import receipt
from ftqc_uncompute import estimate_reaction_decoder_latency


class InstalledActionTests(unittest.TestCase):
    def test_route_with_spectator(self):
        c = action_domain.load("bb:2,1,x0+x1,x0+x1")
        p = [1,0,3,2]
        a = automorphisms.build(c, action_domain.default_basis(c), [{"perm": p}])
        spec = dict(kind="installed-action", code="bb:2,1,x0+x1,x0+x1", action=a, width=5,
                    ticks=[[["SWAP",0,4],["SWAP",2,3]], [["SWAP",1,4]], [["SWAP",0,4]]],
                    couplers=[[0,4],[1,4],[2,3]])
        res = receipt(spec)["result"]
        self.assertEqual((res["ticks"], res["idle_qubit_ticks"]), (3,7))
        self.assertEqual(res["faults"]["post_gate_and_idle_faults"], 4*15+7*3)
        metrics = dict(num_data_qubits=4, peak_ancilla_qubits=1, rounds_per_reaction=1,
                       ticks_per_round=2, installed_couplers=3)
        cost = estimate_reaction_decoder_latency(2, metrics, installed_action_per_reaction=spec)
        self.assertEqual(cost["total_wall_ticks_ceil"], 10)
        for change in (dict(couplers=[[0,4]]), dict(ticks=spec["ticks"][:-1]),
                       dict(ticks=spec["ticks"]+[[["Z",0]]])):
            with self.assertRaises(ValueError):
                receipt(dict(spec, **change))

