import copy
import random
import unittest
from fractions import Fraction
from itertools import combinations

from . import filter_models as fm
from . import filters as fl
from .__main__ import receipt
from .tails import binomial_tail, cone_bound_result

F = Fraction
MARGIN = 1e-12


def model_run(variant, seq, accept, n=2, wmin=2, carry=None, independent=None, pi=(0.03, 0.02), rate=0.01, meas=0.01):
    """The exact model against the calculus, level by level.  Returns (largest ratio of a
    true probability to its bound, list of violated statements)."""
    pi = dict(X=pi[0], Z=pi[1])
    carry = (variant != "independent") if carry is None else carry
    independent = (variant == "independent") if independent is None else independent
    law = fm.leaf_law(n, variant, pi, rate)
    R = lambda x: F(str(x))
    g, m = R(rate), R(meas)
    if variant == "independent":
        rates = dict(gate=rate, meas=meas, idle=rate / 2)
        fresh = 1 - (1 - g) * (1 - g / 2)
        step = dict(fresh=dict(X=fresh, Z=fresh), partner=1 - (1 - g) * (1 - m), clean=((1 - g) * (1 - m)) ** n)
        C = {t: fl.sector(R(pi[t]), [(1, g)]) for t in "XZ"}
    else:
        rates = dict(source=rate, meas=meas)
        fresh = 1 - (1 - g) ** 2 if variant == "overlap" else g
        step = dict(fresh=dict(X=fresh, Z=fresh), partner=1 - (1 - g) ** 2 * (1 - m), clean=((1 - g) ** 2 * (1 - m)) ** n)
        C = {t: fl.sector(R(pi[t]), [(1, 0)], fresh) for t in "XZ"}
    events = fm.increasing_events(n)
    worst, bad = 0.0, []
    for lvl, s in enumerate(seq):
        law, p_acc = fm.step(law, law, s, n, wmin, variant, rates, accept)
        C, alpha = fl.check(C, C, s, n, wmin, step, carry, independent)
        if float(alpha) > p_acc + MARGIN:
            bad.append((lvl, s, "acceptance", float(alpha), p_acc))
        for t in "XZ":
            junk = fm.junk_probability(law, t)
            if junk > float(C[t]["pi"]) + MARGIN:
                bad.append((lvl, t, "junk", junk, float(C[t]["pi"])))
            modes = fl.total(C[t], independent)
            for ev in events:
                true = fm.event_probability(law, t, ev)
                bound = sum(float(w) * fm.iid_probability(n, float(q), ev) for w, q in modes)
                worst = max(worst, true / bound)
                if true > bound + MARGIN:
                    bad.append((lvl, t, sorted(ev), true, bound))
    return worst, bad


class Arithmetic(unittest.TestCase):
    def test_rounding_brackets_the_value(self):
        for x in (F(1, 3), F(2, 7) * F(10) ** -40, F(123456789, 1000), F(10) ** 30 + 1):
            lo, hi = fl.down(x, 12), fl.up(x, 12)
            self.assertLessEqual(lo, x)
            self.assertLessEqual(x, hi)
            self.assertLess(hi - lo, x * F(10) ** -10)
        self.assertEqual(fl.down(F(1, 4), 5), F(1, 4))
        self.assertEqual(fl.down(0), 0)
        self.assertEqual(fl.decimal(F(1, 8)), "0.125")
        self.assertEqual(fl.decimal(F(7)), "7")
        self.assertEqual(fl.decimal(fl.down(F(2, 3), 5)), "0.666666")

    def test_union_of_rates(self):
        self.assertEqual(fl.unite(True, F(1, 2), F(1, 3)), F(2, 3))
        self.assertEqual(fl.unite(False, F(1, 2), F(1, 3)), F(5, 6))
        self.assertEqual(fl.unite(False, F(2, 3), F(2, 3)), 1)
        self.assertGreaterEqual(fl.unite(False, F(1, 9), F(1, 7)), fl.unite(True, F(1, 9), F(1, 7)))

    def test_leaf_contract(self):
        c = fl.leaf(F(1, 100), 6, 3)
        self.assertEqual(c["pi"], 1 - F(99, 100) ** 6)
        self.assertEqual(c["modes"], [(1, 1 - F(99, 100) ** 3)])
        self.assertEqual(c["last"], 0)
        self.assertEqual(fl.leaf(F(1, 100), 0, 0)["pi"], 0)

    def test_pruning_never_lowers_an_event(self):
        rnd = random.Random(5)
        for _ in range(40):
            modes = [(F(rnd.randint(1, 50), 100), F(rnd.randint(1, 400), 1000)) for _ in range(rnd.randint(3, 12))]
            for cap in (1, 2, 4):
                for order in (1, 3):
                    small = fl.prune(modes, cap, order)
                    self.assertLessEqual(len(small), cap)
                    self.assertEqual(sum(w for w, _ in small), sum(w for w, _ in modes))
                    for n, k in ((5, 0), (5, 2), (9, 4)):
                        before = sum(w * binomial_tail(n, k, q) for w, q in modes)
                        after = sum(w * binomial_tail(n, k, q) for w, q in small)
                        self.assertGreaterEqual(after, before)
        self.assertEqual(fl.prune([(F(1, 2), F(1, 5)), (F(1, 4), F(1, 5)), (0, F(1, 2))]), [(F(3, 4), F(1, 5))])
        # the order decides which neighbours are merged: a heavy mode is kept apart from a dirty one
        modes = [(F(1), F(1, 1000)), (F(1, 100), F(2, 1000)), (F(1, 10 ** 6), F(1, 10))]
        self.assertEqual(fl.prune(modes, 2, 4)[0][1], F(2, 1000))

    def test_reset_is_the_stated_inequality(self):
        """P_0 + sum_k w_k (P_k - P_0) never exceeds what the returned modes give, for every
        increasing event of independent coordinates, also when the weights exceed one."""
        rnd = random.Random(11)
        n = 3
        events = fm.increasing_events(n)
        for trial in range(60):
            base = F(rnd.randint(0, 80), 1000)
            inherited = [(F(rnd.randint(1, 90), 100), F(rnd.randint(1, 300), 10 ** rnd.randint(3, 6)))
                         for _ in range(rnd.randint(1, 5))]
            for independent in (True, False):
                modes = fl.reset(inherited, base, n, independent)
                for ev in events:
                    p0 = fm.iid_probability(n, base, ev)
                    stated = p0 + sum(w * (fm.iid_probability(n, fl.unite(independent, base, r), ev) - p0)
                                      for w, r in inherited)
                    got = sum(w * fm.iid_probability(n, q, ev) for w, q in modes)
                    self.assertGreaterEqual(got, stated, (trial, independent, base, inherited))
        self.assertEqual(fl.reset([(F(1, 2), 0)], F(1, 10), 4), [(1, F(1, 10))])
        self.assertEqual(fl.reset([(F(1, 4), F(1, 10))], F(1, 10), 4, True), [(F(3, 4), F(1, 10)), (F(1, 4), F(19, 100))])
        # weights above one: the correction keeps the leading coefficient near one
        heavy = fl.reset([(F(11, 10), F(1, 10 ** 6))], F(1, 100), 5, True)
        self.assertEqual(len(heavy), 1)
        self.assertLess(heavy[0][0], F(1001, 1000))
        self.assertGreaterEqual(heavy[0][0], 1)

    def test_one_check_by_hand(self):
        """A check of two single-mode blocks without junk: every number written out."""
        n, w = 5, 3
        a = dict(X=fl.sector(0, [(1, F(1, 10))]), Z=fl.sector(F(1, 50), [(1, F(1, 20))]))
        step = dict(fresh=dict(X=F(1, 100), Z=F(1, 200)), partner=F(1, 50), clean=F(2, 7))
        c, alpha = fl.check(a, a, "X", n, w, step, False, True)
        nu = 1 - F(9, 10) ** n
        self.assertEqual(alpha, fl.down((1 - nu) ** 2 * F(2, 7)))
        self.assertLess(alpha, (1 - nu) ** 2 * F(2, 7))
        tb = 1 - F(9, 10) * F(49, 50)
        p_w = binomial_tail(n, w - 1, fl.up(1 - F(9, 10) * (1 - tb)))
        self.assertEqual(c["X"]["pi"], fl.up(p_w / alpha))
        inherited = F(1, 10) * tb
        want = [((1 - 1) / alpha, F(1, 100)), (1 / alpha, 1 - F(99, 100) * (1 - inherited))]
        self.assertEqual(c["X"]["modes"], [(fl.up(wt), fl.up(q)) for wt, q in want if wt > 0])
        self.assertEqual(c["Z"]["pi"], fl.up(F(2, 50) / alpha))
        self.assertEqual(c["Z"]["modes"], [(fl.up(1 / alpha), fl.up(1 - F(19, 20) ** 2 * F(199, 200)))])
        # junk on one side leaves the other side's dirt unsquared, with the junk probability as weight
        b = dict(X=fl.sector(F(1, 10), [(1, F(1, 10))]), Z=a["Z"])
        c2, alpha2 = fl.check(a, b, "X", n, None, step, False, True)
        rates = sorted(q for _, q in c2["X"]["modes"])
        self.assertIn(fl.up(1 - F(99, 100) * F(9, 10)), rates)
        self.assertEqual(c2["X"]["pi"], 0)
        with self.assertRaises(ValueError):
            fl.check(a, a, "X", n, w, step, True, True)
        dead, alpha3 = fl.check(a, a, "X", n, w, dict(step, clean=0), False, True)
        self.assertIsNone(dead)
        self.assertEqual(alpha3, 0)

    def test_fold_shares_identical_subtrees(self):
        p = F(1, 1000)
        lf = dict(X=fl.leaf(p, 6, 2), Z=fl.leaf(p, 5, 3), cost=16)
        nodes, root = fl.uniform("XZX")
        step = dict(fresh=dict(X=2 * p, Z=2 * p), partner=2 * p, clean=(1 - p) ** 14, cost=21)
        for nd in nodes:
            nd["step"] = step
        r = fl.fold(7, dict(X=4, Z=3), {i: lf for i in range(8)}, nodes, root)
        self.assertEqual(len(r["accept_lower"]), 3)
        self.assertEqual(len(nodes), 7)
        a0, a1, a2 = r["accept_lower"]
        self.assertEqual(r["cost"], fl.up((2 * fl.up((2 * fl.up((2 * 16 + 21) / a0) + 21) / a1) + 21) / a2))
        # a different leaf somewhere makes its ancestors different
        leaves = {i: lf for i in range(8)}
        leaves[5] = dict(X=fl.leaf(p, 7, 2), Z=fl.leaf(p, 5, 3), cost=16)
        self.assertEqual(len(fl.fold(7, dict(X=4, Z=3), leaves, nodes, root)["accept_lower"]), 5)


class AgainstExactModels(unittest.TestCase):
    def test_independent_locations(self):
        for seq, accept, wmin in (("XZ", "always", 2), ("XX", "adverse", 2), ("ZXZ", "adverse", 3), ("XZX", "clean", 2)):
            worst, bad = model_run("independent", seq, accept, wmin=wmin)
            self.assertEqual(bad, [], (seq, accept, wmin))
            self.assertLess(worst, 1)
        worst, _ = model_run("independent", "XX", "adverse", wmin=3)
        self.assertGreater(worst, 0.9)                         # and the bound is not slack

    def test_shared_sources_with_carried_dirt(self):
        for variant in ("shared", "overlap"):
            for seq, accept in (("XX", "adverse"), ("XZX", "always"), ("ZZ", "adverse")):
                worst, bad = model_run(variant, seq, accept)
                self.assertEqual(bad, [], (variant, seq, accept))
        worst, _ = model_run("shared", "XX", "adverse", wmin=3)
        self.assertGreater(worst, 0.9)

    def test_multiplying_the_latest_dirt_is_wrong_when_sources_are_shared(self):
        """The reason for carrying: without it the model's true probabilities exceed the bound."""
        worst, bad = model_run("shared", "XX", "always", carry=False, independent=False)
        self.assertGreater(worst, 1.5)
        self.assertTrue(any(isinstance(b[2], list) for b in bad))

    def test_conditioning_matters(self):
        """Dropping the division by the acceptance bound is caught by the adverse acceptance."""
        orig = fl.check

        def unconditioned(*a, **k):
            c, alpha = orig(*a, **k)
            return {t: dict(pi=c[t]["pi"] * alpha, modes=[(w * alpha, q) for w, q in c[t]["modes"]], last=c[t]["last"])
                    for t in "XZ"}, alpha

        fl.check = unconditioned
        try:
            _, bad = model_run("independent", "XX", "adverse", wmin=3)
        finally:
            fl.check = orig
        self.assertTrue(bad)


def steane_tree(state, seq):
    from ftqc_gadgets import factories, library
    from ftqc_gadgets.effects import Analysis

    g = library.prep("steane", state, "tree=" + seq)
    return g, factories.Tree(g), Analysis(g)


class AgainstPhysicalTrees(unittest.TestCase):
    def compare(self, state, seq, p, weights=None):
        from ftqc_codes import gf2
        from ftqc_gadgets import factories

        g, tr, an = steane_tree(state, seq)
        profile = tr.profile()
        n = profile["coordinates"]
        exact = factories.accepted_law(g, float(p), tr, an, weights)
        leaves, nodes, root = fl.physical_tree(profile, p)
        worst = 0.0
        for independent in (True, False):
            r = fl.fold(n, profile["wmin"], leaves, nodes, root, independent=independent)
            self.assertFalse(r["vanished"])
            # one bound per level of the balanced tree, against the least exact value
            self.assertLessEqual(float(min(r["accept_lower"])), min(exact["accept"]))
            self.assertLessEqual(float(r["accept_lower"][0]), exact["accept"][0])
            for t in "XZ":
                sec = r["contract"][t]
                modes = fl.total(sec, independent)
                weight = factories.class_weights(exact["single"][t], n)
                part = lambda e: e & exact["classes"][t]
                for k in range(n):
                    true = sum(pr for e, pr in exact["law"].items() if weight[part(e)] > k)
                    bound = float(sec["pi"] + sum(w * binomial_tail(n, k, fl.up(q, 12)) for w, q in modes))
                    self.assertLessEqual(true, bound * (1 + 1e-9), (state, seq, t, "weight above", k))
                    worst = max(worst, true / bound)
                for size in (1, 2, n - 1):
                    for T in combinations(range(n), size):
                        off = gf2.Echelon([exact["single"][t][c] for c in range(n) if c not in T])
                        true = sum(pr for e, pr in exact["law"].items() if part(e) not in off)
                        bound = float(sec["pi"] + sum(w * (1 - (1 - q) ** size) for w, q in modes))
                        self.assertLessEqual(true, bound * (1 + 1e-9), (state, seq, t, "needs", T))
                        worst = max(worst, true / bound)
        return worst, exact

    def test_profile_is_read_as_counted(self):
        _, tr, _ = steane_tree("zero", "XZ")
        profile = tr.profile()
        p = F(1, 1000)
        leaves, nodes, root = fl.physical_tree(profile, p)
        sh = profile["shapes"][0]
        self.assertEqual((sh["X"]["junk"], max(sh["X"]["local"]), sh["Z"]["junk"], max(sh["Z"]["local"])), (6, 2, 5, 3))
        self.assertEqual(len(leaves), 4)
        for lf in leaves.values():
            self.assertEqual(lf["X"], fl.leaf(p, 6, 2))
            self.assertEqual(lf["Z"], fl.leaf(p, 5, 3))
            self.assertEqual(lf["cost"], sh["locations"])
        self.assertEqual([nd["sector"] for nd in nodes], ["X", "X", "Z"])
        two = 1 - (1 - p) ** 2
        for nd in nodes:                                        # per coordinate: a gate, a measurement, an idle
            self.assertEqual(nd["step"], dict(fresh=dict(X=two, Z=two), partner=two, clean=(1 - p) ** 14, cost=21))
        self.assertEqual(root, ("node", 2))
        self.assertEqual(nodes[2]["A"][0], "node")

    def test_steane_trees_under_uniform_faults(self):
        worst, exact = self.compare("zero", "XZ", F(1, 1000))
        self.assertGreater(worst, 0.25)
        self.assertAlmostEqual(sum(exact["law"].values()), 1.0, places=12)
        self.compare("zero", "ZX", F(1, 500))
        self.compare("plus", "XZ", F(1, 1000))
        self.compare("zero", "XX", F(1, 500))

    def test_steane_tree_under_one_sided_faults(self):
        """Every location, when faulty, carries its first (or last) label: nothing in the
        calculus depends on the labels being uniform."""
        for pick in (0, -1):
            weights = lambda loc, pick=pick: {sorted(loc.faults)[pick]: 1.0}
            self.compare("zero", "XZ", F(1, 500), weights)

    def test_exact_law_agrees_with_enumeration_of_single_faults(self):
        """To first order in p the exact law is the sum over single faults of the test tool."""
        from ftqc_gadgets import factories

        g, tr, an = steane_tree("zero", "X")
        p = 1e-7
        exact = factories.accepted_law(g, p, tr, an)
        first = {}
        for loc in an.model.locations:
            for label, eff in loc.faults.items():
                sg = an.sig(*eff)
                if not sg & an.DET:
                    e = sg >> an.shift
                    first[e] = first.get(e, 0.0) + p / len(loc.faults)
        for e, pr in first.items():
            if e:
                self.assertAlmostEqual(exact["law"][e] / pr, 1.0, places=4)
        self.assertEqual(set(exact["law"]) - {0}, {e for e in exact["law"] if e})
        self.assertTrue(set(first) <= set(exact["law"]))


RAW = {"name": "raw outer", "inputs": [], "outputs": [0, 1, 2, 3, 4, 5, 6],
       "layers": [[["NEW", "plus", 0, 1, 2], ["NEW", "zero", 3, 4, 5, 6]],
                  [["CX", 0, 3], ["CX", 1, 4], ["CX", 2, 5]],
                  [["CX", 0, 4], ["CX", 1, 5], ["CX", 2, 6]],
                  [["CX", 0, 6], ["CX", 1, 3], ["CX", 2, 4]]]
                 + 2 * [[["ECa", b] for b in range(7)], [["ECb", b] for b in range(7)]]}
INNER = {k: {s: {"pi": "1/1000000000", "modes": [["1", "1/20000"]]} for s in "XZ"} for k in ("plus", "zero")}


def encoded(**change):
    spec = dict(kind="filter-contract", level="encoded", p="1/100000", inner={"n": 7, "t": 1}, contracts=INNER,
                raw={"program": copy.deepcopy(RAW), "last_layers": 2},
                tree={"n": 7, "wmin": {"X": 4, "Z": 3}, "sequence": "XZ"})
    spec.update(change)
    return spec


class Study(unittest.TestCase):
    def test_physical_from_a_factory_receipt(self):
        from ftqc_gadgets import certify

        g, tr, _ = steane_tree("zero", "XZ")
        rec = certify.build(g, tree=dict(strict=False, light=True))
        spec = dict(kind="filter-contract", level="physical", p="1/1000", factory=rec)
        r = fl.filter_contract_result(spec)
        self.assertEqual(r["status"], "conditional")
        self.assertEqual(r["tree"]["leaves"], 4)
        leaves, nodes, root = fl.physical_tree(tr.profile(), F(1, 1000))
        direct = fl.fold(7, tr.profile()["wmin"], leaves, nodes, root)
        self.assertEqual(r["contract"], fl.as_json({s: fl.summary(direct["contract"][s], 4, True) for s in "XZ"}))
        self.assertEqual(len(r["accept_lower"]), 2)
        self.assertLess(F(r["accept_lower"][0]), 1)
        # the contract is in the form the tail bounds read
        small = {"name": "small", "inputs": ["A"], "outputs": ["A"], "layers": [[["ECa", "A"]], [["ECb", "A"]]]}
        cb = cone_bound_result(dict(kind="cone-bound", noise="independent", p="1/1000", inner={"n": 7, "t": 1},
                                    contracts={"plus": r["contract"], "zero": r["contract"]}, program=small))
        self.assertEqual(cb["status"], "conditional")
        # a receipt that does not verify is refused
        bad = copy.deepcopy(rec)
        bad["gadget"]["circuit"]["ticks"][0] = bad["gadget"]["circuit"]["ticks"][0][:-1]
        with self.assertRaises(ValueError):
            fl.filter_contract_result(dict(spec, factory=bad))
        with self.assertRaises(ValueError):
            fl.filter_contract_result(dict(spec, leaf={"locations": 1, "X": {"junk": 0, "local": 0}, "Z": {"junk": 0, "local": 0}}))

    def test_physical_from_counts(self):
        spec = dict(kind="filter-contract", level="physical", p="1/10000",
                    tree={"n": 23, "wmin": {"X": 8, "Z": 7}, "sequence": "XZXZ"},
                    leaf={"locations": 151, "X": {"junk": 70, "local": 5}, "Z": {"junk": 71, "local": 5}})
        r = fl.filter_contract_result(spec)
        p = F(1, 10000)
        nu = 1 - (1 - p) ** 70 + 1 - (1 - (1 - (1 - p) ** 5)) ** 23
        self.assertEqual(r["accept_lower"][0], fl.decimal(fl.down(fl.down((1 - nu) ** 2 * (1 - p) ** 46), 20)))
        self.assertEqual(r["tree"]["checks"], 15)
        self.assertEqual(len(r["contract"]["X"]["modes"]), 4)
        self.assertEqual(len(fl.filter_contract_result(dict(spec, modes=2))["contract"]["X"]["modes"]), 2)
        rec = receipt(spec)
        self.assertEqual(rec["schema"], "ftqc_compose-1")
        self.assertEqual(rec, receipt(copy.deepcopy(spec)))
        strong = fl.filter_contract_result(dict(spec, p="1/100"))
        self.assertEqual(strong["status"], "no bound")
        for wrong in (dict(spec, tree=dict(spec["tree"], sequence="XY")), dict(spec, modes=0), dict(spec, p="2"),
                      dict(spec, level="logical"), dict(spec, extra=1)):
            with self.assertRaises(ValueError):
                fl.filter_contract_result(wrong)

    def test_banked_example(self):
        import json
        from pathlib import Path

        spec = json.loads((Path(__file__).parent / "examples" / "filter-contract-golay.json").read_text())
        r = receipt(spec)["result"]
        self.assertEqual(r["status"], "conditional")
        self.assertEqual(r["tree"]["leaves"], 64)
        self.assertEqual((r["junk"]["X"], r["junk"]["Z"]), ("1.597e-15", "3.042e-13"))
        self.assertTrue(all(F(9, 10) < F(a) < 1 for a in r["accept_lower"]))

    def test_sources_of_one_column(self):
        rows = fl.shared_sources("XZXZ")
        self.assertEqual([r["sector"] for r in rows], list("XZXZ"))
        self.assertTrue(all(r["rest_disjoint"] and r["last_shared"] for r in rows))
        # without the survivor's correction the old dirt keeps its sources and may not be multiplied
        self.assertFalse(all(r["rest_disjoint"] for r in fl.shared_sources("XX", correct=False)))

    def test_encoded_tree(self):
        r = fl.filter_contract_result(encoded())
        self.assertEqual(r["status"], "conditional")
        self.assertTrue(r["carry_last"])
        self.assertEqual(r["tree"], dict(leaves=4, checks=3, sequence="XZ"))
        for s in "XZ":
            self.assertGreater(F(r["contract"][s]["pi"]), 0)
            self.assertLess(F(r["contract"][s]["pi"]), F(r["raw"]["junk"][s]))       # purification helps
        # one check, from the tail bounds of its three decoding events
        from .tails import event_bounds, sci

        for s, gate in (("X", ["CX", "A", "B"]), ("Z", ["CX", "B", "A"])):
            prog = dict(name="one check", inputs=[], outputs=["A"],
                        layers=[[["NEW", "std", "A", "B"], gate], [["MZ" if s == "X" else "MX", "B"], ["ECa", "A"]],
                                [["ECb", "A"]]])
            probs, _ = event_bounds(dict(kind="cone-bound", noise="independent", p="1/100000", inner={"n": 7, "t": 1},
                                         contracts=INNER, program=prog))
            meas, eca, ecb = probs[0], probs[1], probs[2]
            self.assertEqual(r["step"][s]["partner"], sci(meas))
            self.assertEqual(r["step"][s]["fresh"], {"X": sci(eca), "Z": sci(ecb)})
            self.assertEqual(r["step"][s]["clean_lower"], fl.decimal(fl.down(fl.down(1 - 7 * (meas + eca + ecb)), 12)))
            self.assertGreater(meas, 0)
        self.assertEqual(r["cost_unit"], "raw preparations")
        self.assertGreater(F(r["expected_cost"]), 4)
        self.assertLess(F(r["expected_cost"]), F(41, 10))
        priced = fl.filter_contract_result(encoded(costs={"raw": "1000", "X": "10", "Z": "20"}))
        self.assertGreater(F(priced["expected_cost"]), 4000 + 2 * 10 + 20)
        self.assertEqual(priced["contract"], r["contract"])

    def test_encoded_tree_from_a_factory_receipt(self):
        from ftqc_gadgets import certify, library

        g = library.load("prep:steane,zero,tree=XZ,rounds=1")
        strict = certify.build(g, tree=dict(strict=True, light=True))
        spec = encoded()
        del spec["tree"]
        r = fl.filter_contract_result(dict(spec, factory=strict))
        self.assertEqual(r["tree"]["sequence"], "XZ")
        self.assertEqual(r["contract"], fl.filter_contract_result(encoded())["contract"])
        # classes that trust more of the leaf than its last round do not describe an encoded tree
        loose = certify.build(g, tree=dict(strict=False, light=True))
        with self.assertRaisesRegex(ValueError, "strict"):
            fl.filter_contract_result(dict(spec, factory=loose))

    def test_encoded_refusals(self):
        # the first of the two rounds is not local to a block: its cones reach the blocks it was entangled with
        with self.assertRaisesRegex(ValueError, "shares sources"):
            fl.filter_contract_result(encoded(raw={"program": copy.deepcopy(RAW), "last_layers": 4}))
        # one round only: its cones still hold what the encoding gates spread
        short = copy.deepcopy(RAW)
        short["layers"] = short["layers"][:-2]
        with self.assertRaisesRegex(ValueError, "shares sources"):
            fl.filter_contract_result(encoded(raw={"program": short, "last_layers": 2}))
        idle = copy.deepcopy(RAW)
        idle["layers"].append([])
        with self.assertRaisesRegex(ValueError, "just corrected"):
            fl.filter_contract_result(encoded(raw={"program": idle, "last_layers": 3}))
        nothing = copy.deepcopy(RAW)
        nothing["layers"] = nothing["layers"][:-4]
        with self.assertRaises(ValueError):
            fl.filter_contract_result(encoded(raw={"program": nothing, "last_layers": 2}))
        with self.assertRaisesRegex(ValueError, "7 blocks"):
            fl.filter_contract_result(encoded(tree={"n": 23, "wmin": {"X": 8, "Z": 7}, "sequence": "XZ"}))
        two = copy.deepcopy(RAW)
        two["layers"][0].append(["NEW", "zero2", 9])
        with self.assertRaises(ValueError):
            fl.filter_contract_result(encoded(raw={"program": two, "last_layers": 2}))
        with self.assertRaises(ValueError):
            fl.filter_contract_result(encoded(contracts={"plus": INNER["plus"]}))


if __name__ == "__main__":
    unittest.main()
