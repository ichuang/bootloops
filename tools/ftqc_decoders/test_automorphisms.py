import unittest

from .decimal_bounds import RefusalError
from .dem import DecodingProblem
from .effective_distance import run_decoder


class EnsembleTests(unittest.TestCase):
    def test_noninvolutive_redundant_checks(self):
        problem = DecodingProblem.from_matrices([3, 6, 5], [7], 3)
        cfg = dict(algorithm="automorphism_ensemble", permutations=[[1, 2, 0]])
        for fault in range(8):
            syndrome = problem.syndrome_of(fault)
            result = run_decoder(problem, syndrome, cfg)
            self.assertEqual(result.stats["converged_members"], 2)
            self.assertEqual(syndrome, problem.syndrome_of(result.error_mask))
            self.assertEqual(result.observable_mask, problem.observable_of(result.error_mask))
            self.assertLessEqual(result.soft_weight, run_decoder(problem, syndrome).soft_weight)

    def test_receipt_transport(self):
        from ftqc_codes import action_domain, automorphisms
        from ftqc_codes.orbit_minweight import infer_default_automorphisms
        c = action_domain.load("bb18")
        a = automorphisms.build(c, action_domain.default_basis(c),
                                [{"perm": infer_default_automorphisms(c)[0]}])
        p = DecodingProblem.from_css_code(c)
        cfg = dict(algorithm="automorphism_ensemble", actions=[a])
        r = run_decoder(p, p.syndrome_of(1), cfg)
        self.assertTrue(r.converged)
        self.assertEqual(r.stats["members"], 2)
        with self.assertRaises(RefusalError):
            run_decoder(DecodingProblem.from_css_code("steane"), 0, cfg)

    def test_unsupported_maps(self):
        cfg = dict(algorithm="automorphism_ensemble", permutations=[[1, 2, 0]])
        for problem in [DecodingProblem.from_matrices([3], [], 3),
                        DecodingProblem.from_matrices([7], [1], 3),
                        DecodingProblem.from_matrices([3, 6], [], 3, priors=["1/4","1/4","1/8"])]:
            with self.assertRaises(RefusalError):
                run_decoder(problem, 0, cfg)

