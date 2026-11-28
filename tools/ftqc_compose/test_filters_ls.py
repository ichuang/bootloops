"""The local-stochastic filter calculus against exact laws of small models."""
import unittest
from fractions import Fraction
from itertools import combinations, product

from . import filter_models as fm
from . import filters as fl
from . import filters_ls as ls

F = Fraction


def correlated_leaf(n, pi, rate):
    """A leaf law that is local-stochastic but far from independent: all coordinates are dirty
    together (probability rate^n), junk in a sector comes alone, never with dirt."""
    assert n == 2
    both = rate ** n
    law = {}

    def add(jx, jz, d, pr):
        key = (jx, jz, tuple((d, d, 0) for _ in range(n)))
        law[key] = law.get(key, 0) + pr

    add(0, 0, 1, both)
    add(1, 0, 0, pi["X"])
    add(0, 1, 0, pi["Z"])
    add(0, 0, 0, 1 - both - pi["X"] - pi["Z"])
    return law


def exclusive_leaf(n, pi, rate):
    """Local-stochastic and anticorrelated: exactly one coordinate is dirty, each with
    probability `rate`, so "some coordinate dirty" has the full union weight n * rate."""
    law = {}

    def add(jx, jz, cols, pr):
        law[(jx, jz, cols)] = law.get((jx, jz, cols), 0) + pr

    clean = tuple((0, 0, 0) for _ in range(n))
    for c in range(n):
        add(0, 0, tuple((1, 1, 0) if i == c else (0, 0, 0) for i in range(n)), rate)
    add(1, 0, clean, pi["X"])
    add(0, 1, clean, pi["Z"])
    add(0, 0, clean, 1 - n * rate - pi["X"] - pi["Z"])
    return law


def subsets(n):
    for k in range(1, n + 1):
        for A in combinations(range(n), k):
            yield A


def moment_probability(law, sector, A):
    """P[not junk in the sector, every coordinate of A dirty]."""
    si = 0 if sector == "X" else 1
    return sum(pr for key, pr in law.items() if not key[si] and all(key[2][c][si] for c in A))


def model_run(variant, seq, accept, n=2, wmin=2, carry=None, pi=(F(3, 100), F(2, 100)), rate=F(1, 100),
              meas=F(1, 100), leaf="product", check=None):
    """Exact model against the calculus, level by level: (largest ratio true / bound, violations)."""
    pi = dict(X=pi[0], Z=pi[1])
    carry = (variant != "independent") if carry is None else carry
    g, m = rate, meas
    if variant == "independent":
        rates = dict(gate=g, meas=m, idle=g / 2)
        step = dict(fresh=dict(X=g + g / 2, Z=g + g / 2), partner=g + m, clean=1 - n * (g + m))
        C = {t: fl.sector(pi[t], [(1, g)]) for t in "XZ"}
        law = {"product": lambda: fm.leaf_law(n, variant, pi, g), "correlated": lambda: correlated_leaf(n, pi, g),
               "exclusive": lambda: exclusive_leaf(n, pi, g)}[leaf]()
    else:
        rates = dict(source=g, meas=m)
        fresh = 2 * g if variant == "overlap" else g
        step = dict(fresh=dict(X=fresh, Z=fresh), partner=2 * g + m, clean=1 - n * (2 * g + m))
        C = {t: fl.sector(pi[t], [(1, 0)], fresh) for t in "XZ"}
        law = fm.leaf_law(n, variant, pi, g)
    worst, bad = F(0), []
    for lvl, s in enumerate(seq):
        law, p_acc = fm.step(law, law, s, n, wmin, variant, rates, accept)
        C, alpha = (check or ls.check)(C, C, s, n, wmin, step, carry)
        if C is None:                                       # the weight of rejection reached one: no claim
            break
        if alpha > p_acc:
            bad.append((lvl, s, "acceptance", alpha, p_acc))
        for t in "XZ":
            junk = fm.junk_probability(law, t)
            if junk > C[t]["pi"]:
                bad.append((lvl, t, "junk", junk, C[t]["pi"]))
            modes = ls.total(C[t])
            for A in subsets(n):
                true = moment_probability(law, t, A)
                bound = ls.moment(modes, len(A))
                worst = max(worst, true / bound)
                if true > bound:
                    bad.append((lvl, t, A, true, bound))
    return worst, bad


class Algebra(unittest.TestCase):
    def test_padding_and_unions(self):
        self.assertEqual(ls.pad([(F(1, 4), F(1, 10))]), [(F(1, 4), F(1, 10)), (F(3, 4), 0)])
        self.assertEqual(ls.pad([(F(3, 2), F(1, 10))]), [(F(3, 2), F(1, 10))])
        # union of two processes on disjoint coordinates' sources: sum over splits of A
        a, b = [(F(1, 4), F(1, 10))], [(F(1, 2), F(1, 20)), (F(1, 3), F(1, 5))]
        mom = lambda modes, k: 1 if k == 0 else ls.moment(modes, k)
        for k in (1, 2, 3):
            want = sum(F(len(list(combinations(range(k), j)))) * mom(a, j) * mom(b, k - j) for j in range(k + 1))
            self.assertEqual(ls.moment(ls.united(ls.pad(a), ls.pad(b)), k), want)
        self.assertEqual(ls.leaf(F(1, 100), 6, 3), fl.sector(F(6, 100), [(1, F(3, 100))]))
        self.assertEqual(ls.nontrivial(fl.sector(F(1, 50), [(1, F(1, 100))], F(1, 200)), 4), F(1, 50) + 4 * F(3, 200))

    def test_union_weights_dominate_independent_probabilities(self):
        """With independent sources the weights are union bounds of the true probabilities."""
        p = F(1, 100)
        self.assertGreaterEqual(ls.leaf(p, 6, 3)["pi"], fl.leaf(p, 6, 3)["pi"])
        self.assertGreaterEqual(ls.leaf(p, 6, 3)["modes"][0][1], fl.leaf(p, 6, 3)["modes"][0][1])


class Models(unittest.TestCase):
    def test_bounds_hold_for_every_model_and_acceptance_rule(self):
        for variant in ("independent", "shared", "overlap"):
            for seq in ("X", "XZ", "XX"):
                for accept in ("always", "clean", "adverse"):
                    for wmin in (2, 3):
                        worst, bad = model_run(variant, seq, accept, wmin=wmin, pi=(F(1, 100), F(1, 200)),
                                               rate=F(1, 400), meas=F(1, 400))
                        self.assertEqual(bad, [], (variant, seq, accept, wmin))
                        self.assertLessEqual(worst, 1)
                        self.assertGreater(worst, 0)        # every run reached at least one comparison

    def test_bounds_are_not_slack_by_orders_of_magnitude(self):
        worst, _ = model_run("independent", "X", "always")
        self.assertGreater(worst, F(1, 2))
        worst, _ = model_run("shared", "XX", "adverse", wmin=3)
        self.assertGreater(worst, F(1, 2))

    def test_correlated_local_stochastic_leaf_laws(self):
        """Leaf laws that satisfy the noise condition without independence."""
        for leaf in ("correlated", "exclusive"):
            for seq in ("X", "XZ", "XX"):
                for accept in ("always", "clean", "adverse"):
                    worst, bad = model_run("independent", seq, accept, leaf=leaf, pi=(F(1, 100), F(1, 200)),
                                           rate=F(1, 400), meas=F(1, 400))
                    self.assertEqual(bad, [], (leaf, seq, accept))
        # the calculus for independent faults does not cover the exclusive law: its acceptance
        # bound uses 1 - (1 - q)^n for "some coordinate dirty", and here that event has weight n q
        pi = dict(X=F(1, 100), Z=F(1, 200))
        g = F(1, 400)
        law = exclusive_leaf(2, pi, g)
        zero = dict(gate=F(0), meas=F(0), idle=F(0))
        _, p_acc = fm.step(law, law, "X", 2, 2, "independent", zero, "clean")
        step = dict(fresh=dict(X=F(0), Z=F(0)), partner=F(0), clean=F(1))
        C = {t: fl.sector(pi[t], [(1, g)]) for t in "XZ"}
        _, alpha_independent = fl.check(C, C, "X", 2, 2, step, False, True)
        _, alpha_ls = ls.check(C, C, "X", 2, 2, step, False)
        self.assertGreater(alpha_independent, p_acc)
        self.assertLessEqual(alpha_ls, p_acc)

    def test_planted_defects_are_caught(self):
        orig = ls.check

        def no_retry_factor(*a, **k):                       # forget the division by 1 - rho
            c, alpha = orig(*a, **k)
            return {t: dict(pi=c[t]["pi"] * alpha, modes=[(w * alpha, q) for w, q in c[t]["modes"]], last=c[t]["last"])
                    for t in "XZ"}, alpha

        _, bad = model_run("independent", "XX", "adverse", wmin=3, check=no_retry_factor)
        self.assertTrue(bad)

        def no_one_sided_junk(A, B, s, n, wmin, step, carry=False, independent=False, cap=ls.MAX_MODES):
            zero = {t: dict(B[t], pi=F(0)) for t in "XZ"}    # drop the terms with one junk input
            c, alpha = orig(A, zero, s, n, wmin, step, carry, independent, cap)
            return c, alpha

        _, bad = model_run("independent", "X", "always", wmin=3, pi=(F(1, 5), F(1, 5)), check=no_one_sided_junk)
        self.assertTrue(bad)
        # multiplying the latest dirt when its sources are shared with the partner's side
        _, bad = model_run("shared", "XX", "always", carry=False)
        self.assertTrue(any(isinstance(b[2], tuple) for b in bad))
        with self.assertRaisesRegex(ValueError, "independent"):
            ls.check({}, {}, "X", 2, 2, {}, False, True)


if __name__ == "__main__":
    unittest.main()
