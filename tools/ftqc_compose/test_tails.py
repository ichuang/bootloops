import random
import unittest
from fractions import Fraction
from itertools import combinations, product

from . import tails
from .stochastic import stochastic_result
from .tails import binomial_tail, cone_bound_result, local_stochastic_tail, mixture_tail, threshold_enumerator, up

F = Fraction
SMALL = {"name": "small", "inputs": ["A", "B"], "outputs": ["A"],
         "layers": [[["CX", "A", "B"]], [["ECa", "A"], ["NEW", "plus", "C"]], [["ECb", "A"], ["CX", "C", "B"]],
                    [["MZ", "B", "m"], ["CZ", "A", "C"]], [["CXg", "A", "C", "m"]], [["ECa", "C"]], [["ECb", "C"]],
                    [["MX", "C"]]]}


def spec(p="1/100", noise="independent", program=SMALL, n=7, t=1, contracts=None, **more):
    if contracts is None:
        c = {"pi": "0", "modes": [["1", p]]} if noise == "independent" else {"equivalents": 1}
        contracts = {k: {"X": c, "Z": c} for k in ("plus", "zero")}
    return dict(kind="cone-bound", noise=noise, p=p, inner={"n": n, "t": t}, contracts=contracts,
                program=program, **more)


class TailTests(unittest.TestCase):
    def test_rounding_is_upward_and_tight(self):
        for x in (F(1, 3), F(2, 7) * F(10) ** -40, F(123456789, 1000), F(1), F(10) ** 30 + 1):
            y = up(x, 12)
            self.assertGreaterEqual(y, x)
            self.assertLess(y - x, x * F(10) ** -10)
        self.assertEqual(up(0), 0)
        self.assertEqual(up(F(1, 4), 5), F(1, 4))
        self.assertEqual(tails.sci(F(1, 3)), "3.334e-1")
        self.assertEqual(tails.sci(F(99999, 100000)), "1.000e0")
        self.assertEqual(tails.sci(F(5, 1000)), "5.000e-3")
        for x in (F(7, 13) * F(10) ** -30, F(10) ** 9 + 1, F(1, 7)):
            self.assertLess(up(x) - x, x * F(10) ** -48)
        self.assertGreaterEqual(F(tails.sci(F(2, 3) * F(10) ** -17)), F(2, 3) * F(10) ** -17)

    def test_binomial_tail_is_the_definition(self):
        for n, t, q in ((5, 1, F(1, 7)), (7, 2, F(3, 10)), (6, 0, F(1, 2)), (4, 3, F(1, 3))):
            brute = F(0)
            for bits in product((0, 1), repeat=n):
                if sum(bits) > t:
                    pr = F(1)
                    for b in bits:
                        pr *= q if b else 1 - q
                    brute += pr
            self.assertEqual(binomial_tail(n, t, q), brute)
        self.assertEqual(binomial_tail(5, 5, F(1, 2)), 0)
        self.assertEqual(binomial_tail(5, -1, F(1, 2)), 1)
        with self.assertRaises(ValueError):
            binomial_tail(5, 1, F(3, 2))

    def test_threshold_family_matches_explicit_family(self):
        # 4 coordinates, 2 locations each, bad when more than one coordinate is dirty
        n, t, L = 4, 1, 2
        sites = [f"c{c}l{i}" for c in range(n) for i in range(L)]
        minimal = [[f"c{a}l{i}", f"c{b}l{j}"] for a, b in combinations(range(n), 2) for i in range(L) for j in range(L)]
        r = stochastic_result({"kind": "local-stochastic", "p": "1/50",
                               "source": {"kind": "family", "expression": {"op": "family", "sites": sites,
                                                                             "minimal": minimal}}})
        self.assertEqual(r["enumerator"], {"2": threshold_enumerator(n, t, [L] * n)})
        self.assertEqual(F(r["raw_bound"]), local_stochastic_tail(n, t, L, F(1, 50)))
        self.assertEqual(threshold_enumerator(3, 1, [1, 2, 3]), 1 * 2 + 1 * 3 + 2 * 3)
        self.assertEqual(threshold_enumerator(2, 2, [1, 1]), 0)
        # and the same event under independent noise, by enumeration of all fault sets
        p = F(1, 50)
        brute = F(0)
        for bits in product((0, 1), repeat=n * L):
            dirty = sum(1 for c in range(n) if any(bits[c * L:(c + 1) * L]))
            if dirty > t:
                k = sum(bits)
                brute += p ** k * (1 - p) ** (n * L - k)
        self.assertEqual(binomial_tail(n, t, 1 - (1 - p) ** L), brute)
        self.assertLessEqual(brute, local_stochastic_tail(n, t, L, p))

    def test_mixture_is_exact_for_true_mixtures(self):
        # two ancillas, each a genuine mixture of two independent-dirt laws, plus locations
        n, t, base = 6, 1, F(1, 20)
        modes = [(F(3, 4), F(1, 30)), (F(1, 4), F(1, 5))]
        exact = F(0)
        for (w1, q1), (w2, q2) in product(modes, repeat=2):
            exact += w1 * w2 * binomial_tail(n, t, 1 - (1 - base) * (1 - q1) * (1 - q2))
        got, terms, coarse = mixture_tail(n, t, base, [modes, modes])
        self.assertFalse(coarse)
        self.assertEqual(terms, 3)                       # identical ancillas are grouped
        self.assertGreaterEqual(got, exact)
        self.assertLess(got - exact, exact * F(10) ** -40)
        alone = mixture_tail(n, t, base, [])[0]
        self.assertGreaterEqual(alone, binomial_tail(n, t, base))
        self.assertLess(alone - binomial_tail(n, t, base), F(10) ** -45)
        for nn, tt, q in ((7, 1, F(1, 9)), (23, 3, F(1, 1000)), (47, 5, F(3, 10 ** 7)), (9, 6, F(1, 3))):
            fixed = F(tails._tail_fixed(nn, tt, tails._floor(1 - q)), tails.FIX)
            self.assertGreaterEqual(fixed, binomial_tail(nn, tt, q))
            self.assertLess(fixed - binomial_tail(nn, tt, q), binomial_tail(nn, tt, q) * F(10) ** -100)
        # more dirt never lowers the bound
        self.assertGreater(mixture_tail(n, t, base, [modes, modes, modes])[0], got)
        old = tails.MAX_TERMS
        tails.MAX_TERMS = 4
        try:
            many = [[(F(1, 2), F(1, 100 + i)), (F(1, 2), F(1, 50 + i))] for i in range(5)]
            rough, _, did = mixture_tail(n, t, base, many)
        finally:
            tails.MAX_TERMS = old
        fine = mixture_tail(n, t, base, many)[0]
        self.assertTrue(did)
        self.assertGreaterEqual(rough, fine)
        # the coarsened bound stays close when the mass sits in a few terms
        tails.MAX_TERMS = 8
        try:
            lop = [[(F(999, 1000), F(1, 100)), (F(1, 1000), F(1, 20))]] * 12
            rough = mixture_tail(n, t, base, lop)[0]
        finally:
            tails.MAX_TERMS = old
        fine = mixture_tail(n, t, base, lop)[0]
        self.assertGreaterEqual(rough, fine)
        self.assertLess(rough, fine * F(11, 10))

    def test_cone_bound_on_a_small_program(self):
        r = cone_bound_result(spec(p="1/500"))
        self.assertEqual(r["status"], "conditional")
        self.assertEqual(cone_bound_result(spec(p="1/10"))["upper"], "1")      # a sum of bounds is capped
        self.assertEqual(r["cones"]["inner_decoding_events"], 8)
        self.assertEqual(r["terms"], 9)                  # eight decoding events and one delivery
        total = F(r["upper"])
        self.assertTrue(0 < total < 1)
        # lower rate, lower bound; better code, lower bound
        self.assertLess(F(cone_bound_result(spec(p="1/5000"))["upper"]), total)
        self.assertLess(F(cone_bound_result(spec(p="1/500", n=23, t=3))["upper"]), total)
        # junk is charged: once at the delivery, and at each half of a correction
        c = {"pi": "1/1000", "modes": [["1", "1/500"]]}
        r2 = cone_bound_result(spec(p="1/500", contracts={k: {"X": c, "Z": c} for k in ("plus", "zero")}))
        self.assertEqual(F(r2["upper"]) - total, F(2, 1000) + 4 * F(2, 1000))
        ls = cone_bound_result(spec(p="1/500", noise="local-stochastic"))
        self.assertIn("p^|S|", ls["premises"][-1])
        self.assertGreater(F(ls["upper"]), 0)

    def test_against_the_event_formula(self):
        from ftqc_gadgets import cones
        s = spec(p="1/200")
        probs, terms = tails.event_bounds(s)
        flow = cones.Flow(cones.Program.from_json(SMALL))
        p = F(1, 200)
        for ev, pr in zip(flow.events, probs):
            if ev["node"] is None:
                self.assertEqual(pr, 0)
                continue
            L, atoms = flow.weight(flow.cone(ev))
            k = sum(atoms.values())
            want = binomial_tail(7, 1, 1 - (1 - p) ** (L + k))        # each atom is one more source of rate p
            self.assertGreaterEqual(pr, want)
            self.assertLess(pr - want, want * F(10) ** -30)
            lsb = tails.event_bounds(spec(p="1/200", noise="local-stochastic"))[0]
        self.assertTrue(all(a >= b for a, b in zip(lsb, probs)))

    def test_sampled_failures_stay_below_the_bound(self):
        from ftqc_gadgets import cones, library
        prog = cones.Program.from_json(SMALL)
        d = cones.Direct(prog, library.code_of("steane"))
        p = F(1, 100)
        probs, _ = tails.event_bounds(spec(p="1/100"))
        rnd = random.Random(9)
        shots, first = 6000, {}
        by_loc = {}
        for f in d.faults:
            by_loc.setdefault(f[0], []).append(f)
        locs = sorted(by_loc)
        for _ in range(shots):
            items = [rnd.choice(by_loc[i]) for i in locs if rnd.random() < 0.01]
            bad = d.run(items)[0]
            if bad is not None:
                k = d.events[bad]["flow"]
                first[k] = first.get(k, 0) + 1
        self.assertGreater(sum(first.values()), 30)
        for k, cnt in first.items():
            bound = float(probs[k])
            self.assertLessEqual(cnt / shots, bound + 4 * (bound / shots) ** 0.5 + 4 / shots, (k, cnt, bound))
        self.assertLessEqual(sum(first.values()) / shots, float(sum(probs)))

    def test_two_levels(self):
        prog = {"name": "two", "inputs": [], "outputs": ["Q"],
                "layers": [[["NEW", "zero2", "Q"]], [["ECa", "Q"]], [["ECb", "Q"]], [["EC2a", "Q"]], [["EC2b", "Q"]],
                           [["ECa", "Q"]], [["ECb", "Q"]]]}
        c = {"pi": "1/100000", "modes": [["1", "1/1000"]]}
        contracts = {k: {"X": c, "Z": c} for k in ("plus", "zero", "plus2", "zero2")}
        s = spec(p="1/1000", program=prog, contracts=contracts, outer={"n": 7, "t": 1})
        r = cone_bound_result(s)
        self.assertEqual(r["cones"]["outer_decoding_events"], 4)
        self.assertEqual(set(r["by_tag"]), {"ec2", "out", "junk:zero2", "junk:plus2"})
        self.assertEqual(F(r["by_tag"]["junk:zero2"]), 2 * F(2, 100000))
        probs, terms = tails.event_bounds(s)
        outer = [x for t, x in terms if t in ("ec2", "out")]
        self.assertTrue(all(0 < x < F(1, 100) for x in outer))
        # an outer event is a tail over columns of the summed inner failures in its cone
        self.assertLess(max(outer), binomial_tail(7, 1, min(F(1), sum(probs) + F(2, 1000))) * (1 + F(1, 10 ** 20)))
        # each outer term, from the cones and the inner bounds, written out
        from ftqc_gadgets import cones
        flow = cones.Flow(cones.Program.from_json(prog))
        ev_of = {i: x["event"] for i, x in enumerate(flow.sources) if x["level"] == 2 and x["kind"] == "event"}
        mode = [(F(1), F(1, 1000))]
        want = []
        for o in flow.outer:
            mem = cones._bits(flow.cone_of[o["node"]])
            base = sum((probs[ev_of[i]] for i in mem if i in ev_of), F(0))
            atoms = [mode for i in mem if flow.sources[i]["kind"] == "atom"]
            self.assertGreater(base, 0)
            want.append(mixture_tail(7, 1, min(F(1), up(base)), atoms)[0])
        self.assertEqual(outer, want)
        self.assertGreater(len({len([i for i in cones._bits(flow.cone_of[o["node"]]) if i in ev_of])
                                for o in flow.outer}), 1)
        with self.assertRaisesRegex(ValueError, "outer code"):
            cone_bound_result(spec(p="1/1000", program=prog, contracts=contracts))
        with self.assertRaisesRegex(ValueError, "independent noise only"):
            cone_bound_result(dict(s, noise="local-stochastic",
                                   contracts={k: {"X": {"equivalents": 1}, "Z": {"equivalents": 1}} for k in contracts}))

    def test_refusals(self):
        with self.assertRaisesRegex(ValueError, "no contract"):
            cone_bound_result(spec(contracts={}))
        with self.assertRaisesRegex(ValueError, "single-level"):
            cone_bound_result(spec(outer={"n": 7, "t": 1}))
        with self.assertRaises(ValueError):
            cone_bound_result(spec(p="2"))
        with self.assertRaisesRegex(ValueError, "0 <= t < n"):
            cone_bound_result(spec(n=7, t=7))
        with self.assertRaisesRegex(ValueError, "noise is"):
            cone_bound_result(spec(noise="gaussian"))
        with self.assertRaises(ValueError):
            cone_bound_result(spec(contracts={"plus": {"X": {"pi": "0", "modes": []}, "Z": {"pi": "0", "modes": []}}}))
        bad = dict(SMALL, layers=[[["CX", "A", "B"]]] + SMALL["layers"][1:])
        bad["sha256"] = "0" * 64
        with self.assertRaises(ValueError):
            cone_bound_result(spec(program=bad))

    def test_receipt_replays(self):
        from .__main__ import canonical, receipt
        a = receipt(spec())
        self.assertEqual(canonical(a), canonical(receipt(spec())))
        self.assertNotEqual(a["input_sha256"], receipt(spec(p="1/101"))["input_sha256"])


if __name__ == "__main__":
    unittest.main()
