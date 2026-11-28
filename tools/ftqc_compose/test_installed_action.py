import unittest
from ftqc_codes import action_domain, automorphisms
from ftqc_compose.__main__ import receipt
from ftqc_uncompute import estimate_reaction_decoder_latency


class InstalledActionTests(unittest.TestCase):
    def test_single_fault_histogram_against_forward_propagation(self):
        from collections import Counter
        from itertools import product
        from ftqc_faults.action import single_fault_spread

        ticks = [[["H", 0], ["S", 1]], [["CX", 1, 2]], [["CZ", 0, 2]], [["SWAP", 0, 1]]]
        hist = Counter()
        for i, tick in enumerate(ticks):
            active = {q for op in tick for q in op[1:]}
            locations = [op[1:] for op in tick] + [
                [q] for q in range(3) if q not in active
            ]
            for sites in locations:
                for errors in product(range(4), repeat=len(sites)):
                    if not any(errors):
                        continue
                    x = sum((p & 1) << q for p, q in zip(errors, sites))
                    z = sum((p >> 1) << q for p, q in zip(errors, sites))
                    suffix = [{"layer": layer} for layer in ticks[i + 1 :]]
                    x, z, _ = automorphisms.conjugate((x, z, 0), suffix)
                    hist[str(((x | z) & 3).bit_count())] += 1
        self.assertEqual(
            single_fault_spread(3, 2, ticks)["data_weight_histogram"], dict(hist)
        )

    def test_route_with_spectator(self):
        c = action_domain.load("bb:2,1,x0+x1,x0+x1")
        p = [1, 0, 3, 2]
        a = automorphisms.build(c, action_domain.default_basis(c), [{"perm": p}])
        spec = dict(
            kind="installed-action",
            code="bb:2,1,x0+x1,x0+x1",
            action=a,
            width=5,
            ticks=[
                [["SWAP", 0, 4], ["SWAP", 2, 3]],
                [["SWAP", 1, 4]],
                [["SWAP", 0, 4]],
            ],
            couplers=[[0, 4], [1, 4], [2, 3]],
        )
        res = receipt(spec)["result"]
        self.assertEqual((res["ticks"], res["idle_qubit_ticks"]), (3, 7))
        self.assertEqual(res["faults"]["post_gate_and_idle_faults"], 4 * 15 + 7 * 3)
        metrics = dict(
            num_data_qubits=4,
            peak_ancilla_qubits=1,
            rounds_per_reaction=1,
            ticks_per_round=2,
            installed_couplers=3,
        )
        cost = estimate_reaction_decoder_latency(
            2, metrics, installed_action_per_reaction=spec
        )
        self.assertEqual(cost["total_wall_ticks_ceil"], 10)
        for change in (
            dict(couplers=[[0, 4]]),
            dict(ticks=spec["ticks"][:-1]),
            dict(ticks=spec["ticks"] + [[["Z", 0]]]),
        ):
            with self.assertRaises(ValueError):
                receipt(dict(spec, **change))
