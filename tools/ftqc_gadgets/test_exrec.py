"""Unit and planted-defect tests for table/causal decoders, zero-postselection Steane gadgets,
coherent channel certificates, and exact malignant-pair enumeration."""

import copy
from fractions import Fraction
from itertools import combinations
import unittest

import ftqc_codes
from ftqc_faults.audit.selftest import dense_run
from ftqc_gadgets import Analysis, certify, library, props, roles
from ftqc_gadgets.gadget import Builder, Decoder, GadgetError


class ExRecAndCoherentTests(unittest.TestCase):
    def test_causal_table_decoder_matches_dense_simulation(self):
        # 3-qubit repetition code: inject an X error via a flagged table decoder mid-circuit
        # before a destructive Z measurement, and cross-check against dense_run.
        rep3 = library.code_of("rep3")
        b = Builder("causal_rep3")
        A = b.input(rep3, "data")
        anc = b.alloc(1)[0]
        b.no_detectors()
        b.tick([("R", anc)])
        b.tick([("CX", anc, A[0])])
        b.tick([("CX", anc, A[1])])
        frec = b.tick([("M", anc)])[0]
        # When frec == 1, apply X on A[0] and A[1] (mask 0b011) right before measurement!
        b.table_decoder("hook_fix", [[frec]], "X", A, {1: 0b011}, at_tick=b.cursor)
        recs = b.measure(A, "Z")
        b.decoder("outcomes", recs, rep3.hz, "record")
        g = b.build(idles=True)
        an = Analysis(g)
        # Fault RX->X on anc before CX(anc, A[0]), CX(anc, A[1]) flips frec and puts X_0 X_1 on A;
        # the causal table decoder cancels X_0 X_1 before M(A), leaving zero logical error.
        anc_loc = next(loc for loc in an.model.locations if loc.op.name == "R" and loc.op.qs == (anc,))
        sg = an.replay([(anc_loc.index, "X")])
        self.assertEqual(an.effect(sg), 0)
        # Independent dense statevector check on the same fault!
        where = {i: (k, j) for i, (k, j, _) in enumerate(g.circuit.locations())}
        for _prob, reported, _state in dense_run(g.circuit, {where[anc_loc.index]: "X"}):
            self.assertEqual(reported[0], 1)
            # Raw Z measurements on A[0], A[1] flipped together with anc; applying X_0 X_1 cancels them.
            corr_a = [reported[1] ^ 1, reported[2] ^ 1, reported[3]]
            self.assertEqual(corr_a[0], corr_a[1])
            self.assertEqual(corr_a[1], corr_a[2])

    def test_zero_postselection_steane_set_and_planted_controls(self):
        code = ftqc_codes.load("steane")
        for spec, role, kind, cs in [
            ("prep:steane,zero,goto,idles", "prep_zero", "prep", 1),
            ("prep:steane,plus,goto,idles", "prep_plus", "prep", 1),
            ("knill:steane,anc=goto,idles", "ec", "ec", 1),
            ("steane_ec:steane,anc=goto,idles", "ec", "ec", 1),
            ("phase_s:steane,anc=goto,idles", "s", "ec", 1),
            ("h:steane,idles", "h", "gate", None),
        ]:
            g = library.load(spec)
            roles.logical_action(g, role, code)
            rec = certify.build(g, sector=True, coherent=[(kind, 1, cs)])
            self.assertTrue(rec["results"][0]["coherent"], spec)
            certify.verify(rec, full=True)

        # Planted negative control 1: bare transversal SDG on Steane fails sector=True at r=1, s=0
        # because SDG maps X_0 Z_1 -> X_0 Z_0 Z_1 (weight 2 in Z).
        b = Builder("bare_sdg")
        A = b.input(code, "data")
        b.no_detectors()
        b.tick([("SDG", q) for q in A])
        b.output(A, code, "data")
        g_bare = b.build(idles=True)
        roles.logical_action(g_bare, "s", code)
        bad_coh = props.coherent_channel(Analysis(g_bare), "gate", 1, sector=True)
        self.assertFalse(bad_coh["coherent"])
        self.assertEqual(bad_coh["holds_through"], 0)

        # Planted negative control 2: postselected purification tree is rejected by coherent_channel.
        g_tree = library.load("prep:steane,zero,tree=XZ")
        tree_coh = props.coherent_channel(Analysis(g_tree), "prep", 1, sector=True)
        self.assertFalse(tree_coh["coherent"])
        self.assertIn("postselected", tree_coh["reason"])

        # Planted negative control 3: table entry with nonzero noiseless output is refused.
        with self.assertRaisesRegex(GadgetError, "syndrome 0"):
            Decoder("bad_tab", [0], [[0]], "X", [0], rule={"type": "table", "entries": {"0": 1}})

    def test_malignant_pairs_matches_bruteforce_replay_and_slices(self):
        # Independent route check: compare signature-algebra malignant_pairs against
        # direct circuit propagation (an.replay) over every unordered location pair of Goto prep.
        g = library.load("prep:steane,zero,goto,idles")
        an = Analysis(g)
        W = an.weigher(sector=True)
        mp = props.malignant_pairs(an, bound=1, sector=True)
        self.assertTrue(mp["t1_holds"])

        locs = an.model.locations
        brute_supp = 0
        brute_dep = Fraction(0)
        for loc_i, loc_j in combinations(locs, 2):
            fi = list(loc_i.faults)
            fj = list(loc_j.faults)
            bij = 0
            for li in fi:
                for lj in fj:
                    sg = an.replay([(loc_i.index, li), (loc_j.index, lj)])
                    if W.weight(an.effect(sg), 1) is None:
                        bij += 1
            if bij > 0:
                brute_supp += 1
                brute_dep += Fraction(bij, len(fi) * len(fj))

        self.assertEqual(mp["malignant_support_pairs"], brute_supp)
        self.assertEqual(Fraction(mp["depolarizing_pair_coeff"]), brute_dep)

        # Slice additivity across 3 slices:
        slices = [props.malignant_pairs(an, bound=1, sector=True, pair_slice=(k, 3)) for k in range(3)]
        self.assertEqual(sum(s["total_location_pairs"] for s in slices), mp["total_location_pairs"])
        self.assertEqual(sum(s["malignant_support_pairs"] for s in slices), mp["malignant_support_pairs"])
        self.assertEqual(sum(Fraction(s["depolarizing_pair_coeff"]) for s in slices), brute_dep)

        # Receipt verification and tamper refusal:
        rec = certify.build(g, sector=True, malignant=1)
        certify.verify(rec, full=True)
        tampered = copy.deepcopy(rec)
        tampered["results"][0]["malignant_support_pairs"] += 1
        with self.assertRaisesRegex(certify.Refusal, "does not recompute"):
            certify.verify(tampered, full=True)
        tampered_w = copy.deepcopy(rec)
        tampered_w["results"][0]["first_malignant_pair"] = [[0, "X"], [0, "X"]]
        with self.assertRaisesRegex(certify.Refusal, "not malignant"):
            certify.verify(tampered_w, full=False)

    def test_exrec_union_malignant_pairs(self):
        emp = props.exrec_malignant_pairs("steane", "prep_zero", ec="knill", anc="goto", bound=1, sector=True)
        self.assertTrue(emp["t1_holds"])
        self.assertEqual(emp["truncations"], [[], [0]])
        # Union over {full, truncated} must be at least the full ExRec malignant count alone.
        full_mp = props.malignant_pairs(Analysis(library.load("exrec:steane,prep_zero")), bound=1, sector=True)
        self.assertGreaterEqual(emp["malignant_support_pairs"], full_mp["malignant_support_pairs"])


if __name__ == "__main__":
    unittest.main()
