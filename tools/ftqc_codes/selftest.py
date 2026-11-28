#!/usr/bin/env python3
"""ftqc_codes battery.  Default tier is stdlib-only and finishes in a few seconds.

    python3 tools/ftqc_codes/selftest.py            # default tier
    python3 tools/ftqc_codes/selftest.py --heavy    # adds QR-47 by SAT, BB72/BB90 targets
    python3 tools/ftqc_codes/selftest.py --lean     # adds the Lean-QEC leg (needs lake + a built Lean-QEC)

Legs that need an external engine (cadical, lake) skip with a named line when
the engine is absent; a skip is not a pass and is counted separately.
Exit 0 only if every check that ran passed.
"""
import copy
import json
import os
import random
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from ftqc_codes import (  # noqa: E402
    basis,
    certify,
    cnf,
    codes,
    constructions,
    gf2,
    interchange,
    lrat,
    minweight,
    sat,
    symmetry,
)

PASS, FAIL, SKIP = [], [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    if not cond:
        print(f"FAIL  {name}  {detail}")


def refuses(name, fn, needle=""):
    try:
        fn()
    except certify.Refusal as e:
        check(name, needle in str(e), f"refused for another reason: {e}")
        return
    check(name, False, "was accepted")


def spec_err_or_refuses(name, fn):
    try:
        fn()
    except (codes.SpecError, certify.Refusal):
        check(name, True)
        return
    check(name, False, "was accepted")


def skip(name, why):
    SKIP.append(name)
    print(f"SKIP  {name}: {why}")


# ---------------------------------------------------------------- 1. linear algebra
def t_gf2():
    rnd = random.Random(1)
    for _ in range(200):
        n = rnd.randint(1, 14)
        rows = [rnd.getrandbits(n) for _ in range(rnd.randint(0, 10))]
        ns = gf2.nullspace(rows, n)
        check("gf2.nullspace orthogonal", all(gf2.dot(a, b) == 0 for a in ns for b in rows))
        check("gf2.rank-nullity", gf2.rank(rows) + len(ns) == n)
        check("gf2.nullspace independent", gf2.rank(ns) == len(ns))
        idx = gf2.independent_rows(rows)
        check("gf2.independent_rows", gf2.rank([rows[i] for i in idx]) == len(idx) == gf2.rank(rows))


# ---------------------------------------------------------------- 2. enumeration vs brute force
def t_minweight():
    rnd = random.Random(2)
    for trial in range(300):
        n = rnd.randint(2, 16)
        k = rnd.randint(1, min(n, 9))
        gens = [rnd.getrandbits(n) for _ in range(k)]
        ebasis = gf2.Echelon(gens).basis
        sub = [b for b in ebasis if rnd.random() < 0.4]
        if rnd.random() < 0.3 and len(ebasis) >= 2:
            sub.append(ebasis[0] ^ ebasis[1])
        want, _ = minweight.brute_min_weight(gens, n, sub)
        got = minweight.min_weight(gens, n, sub)
        if want is None:
            check("minweight empty case", got.get("empty") or got["upper"] is None)
            continue
        ok = got["exact"] and got["lower"] == want and gf2.wt(got["witness"]) == want
        check("minweight equals brute force", ok, f"trial {trial}: {got['lower']} vs {want}")
        if want > 1:
            tgt = rnd.randint(1, want)
            g2 = minweight.min_weight(gens, n, sub, target=tgt)
            check("minweight target stop is a valid bound", tgt <= g2["lower"] <= want,
                  f"trial {trial}")
        capped = minweight.min_weight(gens, n, sub, max_combos=rnd.randint(0, 40))
        check("minweight capped bound never exceeds truth",
              capped["lower"] <= want and (capped["upper"] is None or capped["upper"] >= want),
              f"trial {trial}: {capped['lower']} > {want}")


# ---------------------------------------------------------------- 3. families
EXPECT = {  # spec: (n, k, d, wmin_x, wmin_z, self_dual)
    "steane": (7, 1, 3, 4, 4, True),
    "golay": (23, 1, 7, 8, 8, True),
    "qr:47": (47, 1, 11, 12, 12, True),
    "toric:3": (18, 2, 3, 4, 4, False),
    "toric:4": (32, 2, 4, 4, 4, False),
    "planar:3": (13, 1, 3, 3, 3, False),
    "rotated_surface:3": (9, 1, 3, 2, 2, False),
    "rotated_surface:5": (25, 1, 5, 2, 2, False),
    "bb18": (18, 4, 4, 6, 6, True),
    "bb54": (54, 8, 6, 6, 6, False),
    "gb54": (54, 2, 10, 8, 8, True),
}
RECEIPTS = {}


def t_families():
    for spec, (n, k, d, wx, wz, sd) in EXPECT.items():
        code = codes.load(spec)
        rec = certify.build(code)
        RECEIPTS[spec] = rec
        f = rec["facts"]
        got = (f["n"], f["k"], f["d"]["lower"], f["wmin_x"]["lower"], f["wmin_z"]["lower"],
               f["self_dual"])
        check(f"family {spec} parameters", got == (n, k, d, wx, wz, sd), f"{got}")
        check(f"family {spec} exact", all(f[q]["exact"] for q in ("d", "d_x", "d_z")))
        try:
            certify.verify(rec, full=(n <= 32))
            check(f"family {spec} receipt verifies", True)
        except certify.Refusal as e:
            check(f"family {spec} receipt verifies", False, str(e))
        again = json.loads(json.dumps(rec))
        check(f"family {spec} JSON round trip", again == rec)
    f = RECEIPTS["golay"]["facts"]
    check("golay flags", f["doubly_even_x"] and f["all_ones_logical"] and f["t"] == 3)
    check("rotated surface not doubly even", not RECEIPTS["rotated_surface:3"]["facts"]["doubly_even_x"])


# ---------------------------------------------------------------- 4. planted failures
def t_controls():
    base = RECEIPTS["golay"]
    code = codes.load("golay")

    def tamper(fn):
        r = copy.deepcopy(base)
        fn(r)
        return r

    bad = codes.CSSCode("bad", 4, [0b0011], [0b0110])
    refuses("control: non-commuting checks refused", lambda: certify.build(bad), "do not commute")

    def raise_lower(r):
        r["facts"]["d_x"]["lower"] = 8
        r["facts"]["d_x"]["exact"] = False
    refuses("control: lower bound above witness", lambda: certify.verify(tamper(raise_lower)),
            "exceeds upper")

    def drop_evidence(r):
        r["evidence"] = [e for e in r["evidence"] if e["quantity"] != "d_z"]
    refuses("control: lower bound with no evidence", lambda: certify.verify(tamper(drop_evidence)),
            "not backed")

    def light_witness(r):
        w = r["facts"]["d_x"]["witness"]
        i = w.index("1")
        r["facts"]["d_x"]["witness"] = w[:i] + "0" + w[i + 1:]
        r["facts"]["d_x"]["upper"] = 6
        r["facts"]["d_x"]["lower"] = 6
    refuses("control: forged lighter witness", lambda: certify.verify(tamper(light_witness)),
            "parity check")

    def stab_as_logical(r):
        stab = code.hx[0]
        r["facts"]["d_x"].update(witness=certify._bits(stab, 23), upper=gf2.wt(stab))
        r["facts"]["d_x"]["lower"] = min(r["facts"]["d_x"]["lower"], gf2.wt(stab))
        r["facts"]["d_x"]["exact"] = r["facts"]["d_x"]["lower"] == gf2.wt(stab)
    refuses("control: stabilizer passed off as logical",
            lambda: certify.verify(tamper(stab_as_logical)), "partner")

    def flip_matrix(r):
        row = r["code"]["hx"][0]
        r["code"]["hx"][0] = ("1" if row[0] == "0" else "0") + row[1:]
    refuses("control: matrix edited under a stale hash", lambda: certify.verify(tamper(flip_matrix)),
            "sha256")

    def wrong_k(r):
        r["facts"]["k"] = 2
    refuses("control: wrong k", lambda: certify.verify(tamper(wrong_k)), "fact k")

    def wrong_t(r):
        r["facts"]["t"] = 4
    refuses("control: wrong t", lambda: certify.verify(tamper(wrong_t)), "derived")

    # forged evidence claiming more than is true: caught only by full re-proof
    c2 = codes.load("steane")
    r2 = certify.build(c2, targets=dict(d_x=2))   # stops at a true bound below 3
    check("control: target stop reports a bound, not a value",
          r2["facts"]["d_x"]["lower"] >= 2)

    def forge(r):
        for q in ("d_x", "d_z"):
            e = r["facts"][q]
            e["witness"] = e["partner"] = e["upper"] = None
            e["lower"], e["exact"] = 5, False
        for ev in r["evidence"]:
            if ev["quantity"] in ("d_x", "d_z"):
                ev["claim"] = f"{ev['quantity']} >= 5"
        certify._combine(r["facts"])
    forged = copy.deepcopy(RECEIPTS["steane"])
    forge(forged)
    refuses("control: forged enumeration claim caught by --full",
            lambda: certify.verify(forged, full=True), "re-run proves only")
    refuses("control: another leg contradicting a witness",
            lambda: certify.add_lower_bound(copy.deepcopy(base), "d_x", 9, dict(leg="X")),
            "disagree")

    # spec errors
    for name, d in [("row length", dict(n=3, hx=["11"], hz=["110"])),
                    ("missing field", dict(n=3, hx=["110"])),
                    ("bad symbol", dict(n=3, hx=["1a0"], hz=["110"]))]:
        try:
            codes.CSSCode.from_json(d)
            check(f"control: spec error ({name})", False, "accepted")
        except codes.SpecError:
            check(f"control: spec error ({name})", True)


# ---------------------------------------------------------------- 5. CNF direction and LRAT checker
def t_cnf():
    rnd = random.Random(5)
    for spec in ("steane", "golay", "toric:3", "bb18", "rotated_surface:3"):
        code = codes.load(spec)
        f = RECEIPTS[spec]["facts"]
        for q in sat.QUANTITIES:
            true = f[q]["lower"]
            word = certify._unbits(f[q]["witness"])
            for w in (true, true + 1, code.n):
                F = cnf.for_quantity(code, q, w)
                check(f"cnf {spec} {q}: witness satisfies statement", F.holds(word))
                check(f"cnf {spec} {q}: witness extends to a model",
                      F.satisfied_by(F.extend(word)))
            F = cnf.for_quantity(code, q, true - 1)
            check(f"cnf {spec} {q}: witness excluded below its weight",
                  not F.holds(word) and not F.satisfied_by(F.extend(word)))
            # random words: the induced assignment is a model exactly when the statement holds
            F = cnf.for_quantity(code, q, true + 1)
            ker = gf2.nullspace(F.P, code.n)
            for _ in range(60):
                e = rnd.getrandbits(code.n)
                if rnd.random() < 0.6:          # a word passing the parity rows
                    e = 0
                    for v in ker:
                        if rnd.random() < 2.0 / max(2, len(ker)):
                            e ^= v
                check(f"cnf {spec} {q}: extend agrees with statement",
                      F.satisfied_by(F.extend(e)) == F.holds(e))
    # exhaustive agreement of formula and statement on a tiny code, by brute force over all
    # assignments of e (auxiliaries are functionally determined, so extend is the only model)
    code = codes.load("steane")
    for q in sat.QUANTITIES:
        for w in range(0, 5):
            F = cnf.for_quantity(code, q, w)
            sat_e = [e for e in range(1 << 7) if F.satisfied_by(F.extend(e))]
            want = [e for e in range(1 << 7) if F.holds(e)]
            check(f"cnf steane {q} w={w}: model set equals statement set", sat_e == want)


def t_lrat():
    cnf_text = "p cnf 2 4\n1 2 0\n-1 2 0\n1 -2 0\n-1 -2 0\n"
    good = b"5 2 0 1 2 0\n6 0 5 3 4 0\n"
    try:
        st = lrat.check(cnf_text, good)
        check("lrat accepts a valid refutation", st["steps"] == 2)
    except lrat.ProofError as e:
        check("lrat accepts a valid refutation", False, str(e))

    def rejects(name, proof, text=cnf_text):
        try:
            lrat.check(text, proof)
            check(name, False, "accepted")
        except lrat.ProofError:
            check(name, True)
    rejects("lrat rejects missing hint", b"5 2 0 1 0\n6 0 5 3 4 0\n")
    rejects("lrat rejects non-unit hint", b"5 2 0 3 1 0\n6 0 5 3 4 0\n")
    rejects("lrat rejects RAT hint", b"5 2 0 -1 2 0\n6 0 5 3 4 0\n")
    rejects("lrat rejects proof without empty clause", b"5 2 0 1 2 0\n")
    rejects("lrat rejects use of deleted clause", b"5 2 0 1 2 0\n5 d 3 0\n6 0 5 3 4 0\n")
    rejects("lrat rejects unknown hint", b"5 2 0 1 9 0\n6 0 5 3 4 0\n")
    rejects("lrat rejects non-increasing id", b"3 2 0 1 2 0\n6 0 5 3 4 0\n")
    rejects("lrat rejects fresh variable", b"5 3 0 1 2 0\n6 0 5 3 4 0\n")
    rejects("lrat rejects satisfiable CNF", b"5 0 1 2 0\n", "p cnf 2 2\n1 2 0\n-1 2 0\n")
    rejects("lrat rejects empty proof", b"")
    rejects("lrat rejects header mismatch", good, "p cnf 2 5\n1 2 0\n-1 2 0\n1 -2 0\n-1 -2 0\n")

    def enc(x):
        out = bytearray()
        while True:
            b = x & 0x7F
            x >>= 7
            out.append(b | (0x80 if x else 0))
            if not x:
                return bytes(out)

    def lit(v):
        return enc(2 * abs(v) + (1 if v < 0 else 0))
    binary = (b"a" + enc(5) + lit(2) + b"\0" + lit(1) + lit(2) + b"\0"
              + b"a" + enc(6) + b"\0" + lit(5) + lit(3) + lit(4) + b"\0")
    try:
        lrat.check(cnf_text, binary)
        check("lrat accepts the binary format", True)
    except lrat.ProofError as e:
        check("lrat accepts the binary format", False, str(e))
    rejects("lrat rejects truncated binary", binary[:-2])


# ---------------------------------------------------------------- 5b. BCH leg
def t_bch():
    from ftqc_codes import bch
    rnd = random.Random(7)
    check("bch irreducible test", bch._irreducible(0b10011) and not bch._irreducible(0b10101))
    # the bound never exceeds the true minimum weight, on random cyclic codes
    for n in (7, 15, 21, 31):
        m = bch.mult_order_2(n)
        cosets = sorted({min((z * pow(2, i, n)) % n for i in range(m)) for z in range(n)})
        p, alpha = bch.root_of_unity(n)
        for _ in range(12):
            pick = [c for c in cosets if rnd.random() < 0.45]
            Z = set()
            for c in pick:
                Z |= {(c * pow(2, i, n)) % n for i in range(m)}
            if not Z or len(Z) == n:
                continue
            coeffs = [1]
            for z in sorted(Z):
                a = bch._ppowmod(alpha, z, p)
                nxt = [0] * (len(coeffs) + 1)
                for i, c in enumerate(coeffs):
                    nxt[i + 1] ^= c
                    nxt[i] ^= bch._pmulmod(c, a, p)
                coeffs = nxt
            g = sum(c << i for i, c in enumerate(coeffs))
            cb_basis = codes.cyclic_span(g, n)
            r = bch.code_bound(cb_basis, n)
            check(f"bch n={n}: zero set recovered", r["zeros"] == len(Z)
                  and set(r["zero_coset_leaders"]) == set(pick))
            true = minweight.min_weight(cb_basis, n)["lower"]
            check(f"bch n={n}: bound <= true minimum weight", r["bound"] <= true,
                  f"zeros {pick}: bound {r['bound']} > {true}")
    for spec, d in (("cyclic:31,1,3", 5), ("cyclic:63,1,3,5", 7), ("cyclic:15,1", 3)):
        code = codes.load(spec)
        rec = certify.build(code, legs=())
        bch.run_leg(rec, code)
        check(f"bch {spec}: designed distance {d}", rec["facts"]["d"]["lower"] == d)
        full = certify.build(code, max_combos=100_000)      # finds a witness, may not close
        bch.run_leg(full, code)
        check(f"bch {spec}: bound meets an enumerated witness",
              full["facts"]["d"]["lower"] == d and full["facts"]["d"]["exact"])
    code = codes.load("golay")
    rec = certify.build(code)
    bch.run_leg(rec, code)
    check("bch golay: weaker bound does not lower the proven one",
          rec["facts"]["d"]["lower"] == 7 and
          any(e["leg"] == "D" and e["claim"] == "d_x >= 5" for e in rec["evidence"]))
    code = codes.load("toric:3")
    rec = certify.build(code, legs=())
    bch.run_leg(rec, code)
    check("bch control: non-cyclic code is skipped by name",
          all(e["status"] == "skipped" and "odd length" in e["reason"] or "cyclic" in e["reason"]
              for e in rec["evidence"] if e["leg"] == "D"))
    code = codes.load("rotated_surface:3")
    rec = certify.build(code, legs=())
    bch.run_leg(rec, code)
    check("bch control: odd-length non-cyclic code is skipped by name",
          all(e["status"] == "skipped" for e in rec["evidence"] if e["leg"] == "D"))


# ---------------------------------------------------------------- 5c. cut certificates
def t_cuts():
    from ftqc_codes import cuts
    for spec, d in (("rotated_surface:3", 3), ("rotated_surface:5", 5), ("rotated_surface:9", 9),
                    ("planar:3", 3), ("planar:4", 4)):
        code = codes.load(spec)
        rec = certify.build(code, legs=())
        cuts.run_leg(rec, code)
        check(f"cuts {spec}: d >= {d} from disjoint representatives",
              rec["facts"]["d"]["lower"] == d)
        try:
            lines = certify.verify(rec)
            check(f"cuts {spec}: certificate re-checked in cheap verify",
                  sum("cut certificate re-checked" in l for l in lines) == 2)
        except certify.Refusal as e:
            check(f"cuts {spec}: certificate re-checked in cheap verify", False, str(e))
    # the bound never exceeds the truth on codes where the heuristic finds little
    for spec in ("steane", "golay"):
        code = codes.load(spec)
        rec = copy.deepcopy(RECEIPTS[spec])
        cuts.run_leg(rec, code)
        e = [x for x in rec["evidence"] if x["leg"] == "E" and x["status"] == "proved"]
        check(f"cuts {spec}: bound <= true distance",
              all(int(x["claim"].split(">=")[1]) <= rec["facts"]["d"]["upper"] for x in e))
    code = codes.load("rotated_surface:5")
    rec = certify.build(code, legs=())
    cuts.run_leg(rec, code)

    def tamper(fn):
        r = copy.deepcopy(rec)
        e2 = next(x for x in r["evidence"] if x["leg"] == "E" and x["quantity"] == "d_x")
        fn(e2, r)
        return r
    refuses("control: cut set overlapping another",
            lambda: certify.verify(tamper(lambda e, r: e["cuts"].__setitem__(1, e["cuts"][0]))),
            "overlaps")
    stab_x = gf2.supp(code.hz[0])
    refuses("control: stabilizer used as a cut",
            lambda: certify.verify(tamper(lambda e, r: e["cuts"].__setitem__(0, stab_x))),
            "stabilizer")
    refuses("control: non-logical set used as a cut",
            lambda: certify.verify(tamper(lambda e, r: e["cuts"].__setitem__(0, [0]))),
            "not a logical")

    def inflate(e, r):
        e["claim"] = "d_x >= 6"
        r["facts"]["d_x"]["lower"] = 6
        r["facts"]["d_x"]["upper"] = r["facts"]["d_x"]["witness"] = r["facts"]["d_x"]["partner"] = None
        certify._combine(r["facts"])
    refuses("control: cut certificate claiming more sets than it has",
            lambda: certify.verify(tamper(inflate)), "fewer than it claims")
    code = codes.load("toric:3")
    rec = certify.build(code, legs=())
    cuts.run_leg(rec, code)
    check("cuts control: k = 2 is skipped by name",
          all(x["status"] == "skipped" for x in rec["evidence"] if x["leg"] == "E"))
    refuses("control: cut certificate on a k = 2 code",
            lambda: cuts.verify_cuts(code, "d_x", [[0]]), "k = 1")


# ---------------------------------------------------------------- 5d. concatenation
def t_concat():
    from ftqc_codes import concat
    st, st_rec = codes.load("steane"), RECEIPTS["steane"]
    big = codes.concatenate(st, st)
    rec = certify.build(big, max_combos=4_000_000)
    check("concat steane o steane: enumeration gives [[49,1,9]]",
          rec["facts"]["k"] == 1 and rec["facts"]["d"]["lower"] == 9 and rec["facts"]["d"]["exact"])
    r2 = certify.build(big, legs=())
    concat.run_leg(r2, big, st_rec, st_rec)
    check("concat steane o steane: product bound 3 * 3", r2["facts"]["d"]["lower"] == 9)
    try:
        lines = certify.verify(r2)
        check("concat: structure re-checked in cheap verify",
              sum("concatenation structure re-checked" in l for l in lines) == 2)
    except certify.Refusal as e:
        check("concat: structure re-checked in cheap verify", False, str(e))
    big2 = codes.concatenate(codes.load("toric:3"), st)       # k = 2 outer is allowed
    r3 = certify.build(big2, legs=())
    concat.run_leg(r3, big2, RECEIPTS["toric:3"], st_rec)
    check("concat toric3 o steane: product bound 9, k = 2",
          r3["facts"]["d"]["lower"] == 9 and r3["facts"]["k"] == 2)
    # controls
    r4 = certify.build(big, legs=())
    concat.run_leg(r4, big, RECEIPTS["golay"], st_rec)
    check("concat control: wrong outer code is skipped, no bound",
          r4["facts"]["d"]["lower"] is None
          and all(e["status"] == "skipped" for e in r4["evidence"] if e["leg"] == "F"))
    r5 = certify.build(big, legs=())
    concat.run_leg(r5, big, st_rec, RECEIPTS["toric:3"])
    check("concat control: inner code with k = 2 is skipped",
          all(e["status"] == "skipped" for e in r5["evidence"] if e["leg"] == "F"))
    other = codes.concatenate(st, codes.load("rotated_surface:3"))
    r6 = certify.build(other, legs=())
    concat.run_leg(r6, other, st_rec, st_rec)
    check("concat control: different inner code of other length is skipped",
          all(e["status"] == "skipped" for e in r6["evidence"] if e["leg"] == "F"))

    def swap_outer(r):
        for e in r["evidence"]:
            if e["leg"] == "F":
                e["outer"] = RECEIPTS["rotated_surface:3"]["code"]
    bad = copy.deepcopy(r2)
    swap_outer(bad)
    refuses("control: concat evidence with a substituted outer code",
            lambda: certify.verify(bad), "")

    def inflate(r):
        for e in r["evidence"]:
            if e["leg"] == "F":
                e["claim"] = f"{e['quantity']} >= 12"
        for q in ("d_x", "d_z"):
            r["facts"][q]["lower"] = 12
        certify._combine(r["facts"])
    bad = copy.deepcopy(r2)
    inflate(bad)
    refuses("control: concat claim above the product", lambda: certify.verify(bad), "arithmetic")


# ---------------------------------------------------------------- 5e. hierarchical decoding
def t_hier():
    import itertools
    from ftqc_codes import hier
    st = codes.steane()
    big = codes.concatenate(st, st)
    for sector in "XZ":
        h = hier.Hier(st, st, sector)
        check(f"hier {sector}: radii 1 and 1; checks are the concatenated code's",
              (h.t1, h.t2) == (1, 1) and h.checks == (big.hz if sector == "X" else big.hx))
        stab = gf2.Echelon(big.hx if sector == "X" else big.hz)

        def good(x):
            p = h.decode(h.check_values(x))
            return p is not None and (p ^ x) in stab
        # one arbitrary column (all 128 patterns), one error in each of up to two others
        ok = n = 0
        for a in (0, 4):
            others = [c for c in range(7) if c != a]
            for pat in range(128):
                base = pat << (7 * a)
                cases = [base]
                cases += [base | 1 << (7 * c + q) for c in others for q in range(7)]
                for c1, c2 in itertools.combinations(others, 2):
                    cases += [base | 1 << (7 * c1 + q1) | 1 << (7 * c2 + q2)
                              for q1 in range(7) for q2 in range(7)]
                for x in cases:
                    n += 1
                    ok += good(x)
        check(f"hier {sector}: one arbitrary column (two choices, all 128 patterns) and up to two single errors elsewhere, all corrected",
              ok == n and n == 2 * 128 * (1 + 42 + 15 * 49), f"{ok}/{n}")
        rnd = random.Random(4)
        ok = n = 0
        for _ in range(4000):
            a = rnd.randrange(7)
            x = rnd.randrange(128) << (7 * a)
            for c in range(7):
                if c != a and rnd.random() < 0.8:
                    x |= 1 << (7 * c + rnd.randrange(7))
            n += 1
            ok += good(x)
        check(f"hier {sector}: one arbitrary column and one error in every other column (sampled)", ok == n)
        # two bad columns can defeat it, and weight 4 is where that starts
        fail = None
        for q in itertools.combinations(range(14), 4):
            x = sum(1 << p for p in q)
            if len(h.bad_columns(x)) == 2 and not good(x):
                fail = x
                break
        check(f"hier {sector}: some weight-4 error on two columns is not corrected", fail is not None)
        check(f"hier {sector}: bad columns are counted from the pattern",
              h.bad_columns(0b11 | 1 << 7) == [0] and h.bad_columns(0b1) == [])
        # every error of weight <= 3 has at most one bad column
        ok = all(good(sum(1 << p for p in q)) for q in itertools.combinations(range(49), 3)
                 if rnd.random() < 0.15)
        check(f"hier {sector}: sampled weight-3 errors are corrected", ok)
    cert = hier.certificate(RECEIPTS["golay"], RECEIPTS["steane"])
    check("hier: certificate from receipts (Golay over Steane: t1 = 1, t2 = 3, 161 qubits)",
          cert["n"] == 161 and cert["sectors"]["X"] == dict(t1=1, t2=3) and "Step 3" in cert["proof"])
    refuses("hier: an inner code with k = 2 is refused",
            lambda: hier.certificate(RECEIPTS["steane"], RECEIPTS["toric:3"]), "k = 1")
    inflated = copy.deepcopy(RECEIPTS["steane"])
    for quantity in ("d_x", "wmin_x"):
        inflated["facts"][quantity].update(lower=15, upper=None, exact=False, witness=None, partner=None)
        inflated["evidence"].append(dict(leg="A", quantity=quantity, status="proved", claim=f"{quantity} >= 15"))
    refuses("hier: claimed radii must follow from the actual parity checks",
            lambda: hier.certificate(inflated, RECEIPTS["steane"]), "proved kernel radius")
    tab, u = hier.leader_table(codes.golay().hz, 23)
    check("hier: Golay table has 2048 entries, unique through 3", len(tab) == 2048 and u == 3)
    # radius of a decoder from its own matrix: search, then the BCH bound when the search stops short
    from ftqc_codes import bch
    b63 = bch.cyclic_css(63, [1, 3, 5])
    bound, ev = hier.kernel_radius(b63.hz, 63, 3, cap=200_000)
    check("kernel radius: a search that stops short falls back to the BCH bound on the same kernel",
          bound == 3 and ev["kernel_lower"] == 7 and ev["search_lower"] < 7 and "BCH" in ev["bound_from"], str(ev)[:200])
    bound2, ev2 = hier.kernel_radius(b63.hz, 63, 3, cap=60_000_000)
    check("kernel radius: with a larger budget the search alone proves the same radius",
          bound2 == 3 and "bound_from" not in ev2, str(ev2))
    perm = list(range(63))
    perm[0], perm[5] = perm[5], perm[0]
    moved = [sum(1 << perm[i] for i in gf2.supp(r)) for r in b63.hz]
    bound3, ev3 = hier.kernel_radius(moved, 63, 3, cap=200_000)
    check("kernel radius: the same code in a non-cyclic order gets no such help",
          bound3 < 3 and "bound_from" not in ev3, str(ev3))
    bound4, ev4 = hier.kernel_radius(b63.hz, 63, 4, cap=200_000)
    check("kernel radius: the BCH bound does not prove more than the designed distance gives", bound4 == 3)
    bound5, ev5 = hier.kernel_radius(codes.golay().hz, 23, 3)
    check("kernel radius: Golay by search", bound5 == 3 and "bound_from" not in ev5)


# ---------------------------------------------------------------- 5f. interchange & replayable constructions
def t_interchange_and_constructions():
    # 1. Three baseline controls + read-only cross-check against ftqc_extractors
    lp6 = constructions.lifted_product([[3]], [[3]], 3, name="LP(6,2)")
    lp6_ref = constructions.lifted_product_circulant_ref([[3]], [[3]], 3)
    s_lp6 = certify.structure(lp6)
    check("baseline LP(6,2): n=6, k=2, digest c6860dd462fe, max wt 4",
          (s_lp6["n"], s_lp6["k"]) == (6, 2)
          and lp6.sha256().startswith("c6860dd462fe")
          and max(gf2.wt(r) for r in lp6.hx + lp6.hz) == 4)
    check("baseline LP(6,2): matches independent circulant block expander row-for-row",
          lp6.hx == lp6_ref.hx and lp6.hz == lp6_ref.hz and lp6.sha256() == lp6_ref.sha256())

    rad90 = constructions.radial90()
    s_rad = certify.structure(rad90)
    check("baseline radial90: n=90, k=8, digest b656926c5959, max wt 6",
          (s_rad["n"], s_rad["k"]) == (90, 8)
          and rad90.sha256().startswith("b656926c5959")
          and max(gf2.wt(r) for r in rad90.hx + rad90.hz) == 6)

    bb144 = codes.load("bb144")
    s_bb144 = certify.structure(bb144)
    check("baseline bb144: n=144, k=12, digest 0ae98ac8a6e5, max wt 6",
          (s_bb144["n"], s_bb144["k"]) == (144, 12)
          and bb144.sha256().startswith("0ae98ac8a6e5")
          and max(gf2.wt(r) for r in bb144.hx + bb144.hz) == 6)

    try:
        from ftqc_extractors import algebra as ext_alg
        from ftqc_extractors import surgery_products as ext_sp
        check("cross-check LP(6,2) with ftqc_extractors.surgery_products",
              (lp6.n, lp6.hx, lp6.hz) == ext_sp.lifted_product([[3]], [[3]], 3))
        check("cross-check radial90 with ftqc_extractors.surgery_products",
              (rad90.n, rad90.hx, rad90.hz) == ext_sp.radial90())
        a_seed, b_seed = [0b011, 0b110], [0b0011, 0b0110, 0b1100]
        hgp_c = constructions.hgp(a_seed, 3, b_seed, 4)
        check("cross-check rectangular HGP with codes.hgp and ftqc_extractors.algebra.hgp",
              (hgp_c.n, hgp_c.hx, hgp_c.hz) == ext_alg.hgp(a_seed, 3, b_seed, 4)
              and hgp_c.sha256() == codes.hgp(a_seed, 3, b_seed, 4).sha256())
    except ImportError:
        skip("ftqc_extractors read-only cross-check", "ftqc_extractors not present")

    # 2. Exact small X/Z distances vs direct normalizer-mod-stabilizer brute force
    for nm, c_small in [
        ("LP(6,2)", lp6),
        ("LP_rect(1x2,2x2,L=3)", constructions.lifted_product([[3, 5]], [[3, 6], [5, 3]], 3)),
    ]:
        r_sm = certify.build(c_small)
        want_dx, _ = minweight.brute_min_weight(gf2.nullspace(c_small.hz, c_small.n), c_small.n, c_small.hx)
        want_dz, _ = minweight.brute_min_weight(gf2.nullspace(c_small.hx, c_small.n), c_small.n, c_small.hz)
        check(f"{nm}: certified X/Z distances equal direct normalizer-mod-stabilizer brute force",
              r_sm["facts"]["d_x"]["lower"] == want_dx
              and r_sm["facts"]["d_z"]["lower"] == want_dz
              and r_sm["facts"]["d_x"]["exact"] and r_sm["facts"]["d_z"]["exact"])

    # 3. Rectangular seeds, cancelling coefficients, and shift-reversal detection
    rect_lp = constructions.lifted_product([[3, 5]], [[3, 6, 1], [5, 3, 2]], 5)
    rect_ref = constructions.lifted_product_circulant_ref([[3, 5]], [[3, 6, 1], [5, 3, 2]], 5)
    check("rectangular cyclic LP (1x2 by 2x3, L=5) matches circulant reference",
          rect_lp.n == (2 * 2 + 1 * 3) * 5 and rect_lp.sha256() == rect_ref.sha256())

    canc_lp = constructions.lifted_product(
        [[[0, 1, 2, 1]]], [[[0, 2, 2, 1]]], 4, coeff_encoding="exponents"
    )
    clean_lp = constructions.lifted_product([[0b0101]], [[0b0011]], 4, coeff_encoding="bitmask")
    check("mod-2 repeated exponents cancel and match reduced bitmask",
          canc_lp.sha256() == clean_lp.sha256()
          and canc_lp.meta["construction"]["seed_inputs"]["cancelled_terms"] == 4)

    asym_fwd = constructions.lifted_product([[0b0111]], [[0b1011]], 7, shift_convention="hx_neg_hz_pos")
    asym_rev = constructions.lifted_product([[0b0111]], [[0b1011]], 7, shift_convention="hx_pos_hz_neg")
    asym_rev_ref = constructions.lifted_product_circulant_ref(
        [[0b0111]], [[0b1011]], 7, shift_convention="hx_pos_hz_neg"
    )
    check("asymmetric seed detects reversing cyclic shift while matching circulant reference",
          asym_fwd.sha256() != asym_rev.sha256() and asym_rev.sha256() == asym_rev_ref.sha256())

    # 4. General group lifted product: C_3 agreement + nonabelian S_3 + action validation
    c3_perms = [[(g + h) % 3 for h in range(3)] for g in range(3)]
    glp_c3 = constructions.group_lifted_product(3, c3_perms, c3_perms, [[[0, 1]]], [[[0, 1]]])
    check("group_lifted_product on C_3 matches cyclic LP(6,2) row-for-row",
          glp_c3.sha256() == lp6.sha256())

    import itertools
    s3_elems = list(itertools.permutations(range(3)))
    s3_idx = {p: i for i, p in enumerate(s3_elems)}

    def comp(p, q):
        return tuple(p[q[i]] for i in range(3))

    s3_left = [[s3_idx[comp(g, h)] for h in s3_elems] for g in s3_elems]
    s3_right = [[s3_idx[comp(h, g)] for h in s3_elems] for g in s3_elems]
    glp_s3 = constructions.group_lifted_product(
        6, s3_left, s3_right, [[[0, 1, 2]]], [[[0, 3, 4]]], name="S3_LP"
    )
    check("nonabelian S_3 group_lifted_product commutes and records abelian=False",
          glp_s3.commutes()
          and glp_s3.n == 12
          and glp_s3.meta["construction"]["seed_inputs"]["abelian"] is False
          and s3_left != s3_right)

    bad_right = copy.deepcopy(s3_left)
    refuses("control: inconsistent nonabelian left/right group action refused",
            lambda: constructions.group_lifted_product(6, s3_left, bad_right, [[[0, 1]]], [[[0, 1]]]),
            "inconsistent left/right")
    spec_err_or_refuses("control: ragged protograph rejected",
                        lambda: constructions.lifted_product([[1, 2], [1]], [[1]], 3))
    refuses("control: pre-allocation estimator refuses oversized lift before expansion",
            lambda: constructions.lifted_product([[3]], [[3]], 500, max_n=200),
            "exceeds budget")

    # 5. Budget-exhausted Leg A interval + work count
    capped_rad = certify.build(rad90, max_combos=25)
    ev_dx = next(e for e in capped_rad["evidence"] if e["quantity"] == "d_x")
    check("budget-exhausted Leg A returns interval and work count",
          capped_rad["facts"]["d_x"]["lower"] >= 1
          and "budget" in ev_dx
          and isinstance(ev_dx["combos"], int))

    # 6. Three library-survey exchange gates (banked offline in examples/)
    ex_dir = os.path.join(HERE, "examples")
    bb_imp_path = os.path.join(ex_dir, "qldpc_bb144.import.json")
    lp_imp_path = os.path.join(ex_dir, "qldpc_lp6_2.import.json")
    pp_imp_path = os.path.join(ex_dir, "pp256_q4_j3_l8_p16.import.json")

    bb_imp = codes.load(bb_imp_path)
    check("exchange 1 (qLDPC BB144): byte-for-byte digest match + 12 dual pairs",
          bb_imp.sha256() == bb144.sha256()
          and len(bb_imp.meta["candidate_logicals"]["slot_names"]) == 12)

    lp_imp = codes.load(lp_imp_path)
    check("exchange 2 (qLDPC LP(6,2)): raw 6b6ef87d3693 -> native c6860dd462fe via coordinate map",
          lp_imp.meta["import"]["raw_sha256"].startswith("6b6ef87d3693")
          and lp_imp.sha256() == lp6.sha256()
          and len(lp_imp.meta["candidate_logicals"]["slot_names"]) == 2)

    lp_manifest = json.load(open(lp_imp_path))
    bad_map = copy.deepcopy(lp_manifest)
    bad_map["coordinate_map"]["perm_q"] = [1, 0, 2, 5, 3, 4]
    refuses("exchange 2 control: altered physical coordinate map refused",
            lambda: interchange.load_import_manifest(bad_map, base_dir=ex_dir),
            "sha256")

    bad_pair = copy.deepcopy(lp_manifest)
    bad_pair["candidate_logicals"]["lz"]["supports"].reverse()
    del bad_pair["candidate_logicals"]["lz"]["sha256"]
    refuses("exchange 2 control: swapped logical partner refused",
            lambda: interchange.load_import_manifest(bad_pair, base_dir=ex_dir),
            "pairing")

    # Check-order sensitivity: permuting two rows preserves stabilizer group yet changes sha256
    reordered = codes.CSSCode("LP_reordered", lp6.n, [lp6.hx[1], lp6.hx[0], lp6.hx[2]], lp6.hz)
    check("check-order sensitivity: same stabilizer group, distinct ordered sha256",
          sorted(gf2.Echelon(reordered.hx).basis) == sorted(gf2.Echelon(lp6.hx).basis)
          and reordered.sha256() != lp6.sha256())

    pp_imp = codes.load(pp_imp_path)
    pp_srec = certify.build(pp_imp, legs=())
    certify.verify(pp_srec)
    check("exchange 3 (Okada-Kasai PP256): n=256, k=64, ae0e0d3f792d, wt 12/21 witnesses -> d<=12",
          pp_imp.n == 256
          and pp_srec["facts"]["k"] == 64
          and pp_imp.sha256().startswith("ae0e0d3f792d")
          and max(gf2.wt(r) for r in pp_imp.hx + pp_imp.hz) == 9
          and pp_srec["facts"]["d_x"]["upper"] == 12
          and pp_srec["facts"]["d_z"]["upper"] == 21
          and pp_srec["facts"]["d"]["upper"] == 12
          and pp_srec["facts"]["d"]["lower"] is None)

    non_css = copy.deepcopy(lp_manifest)
    non_css["code_type"] = "non_css"
    refuses("interchange control: non-CSS code refused with symplectic message",
            lambda: interchange.load_import_manifest(non_css, base_dir=ex_dir),
            "symplectic")
    refuses("interchange control: duplicate Matrix Market coordinate refused",
            lambda: interchange.parse_matrix_market(
                "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 1\n1 1 1\n"
            ),
            "duplicate")

    # 7. Construction receipts (ftqc_codes-construction-1)
    lp6_srec = certify.build(lp6)
    lp6_crec = constructions.build_construction_receipt(
        lp6,
        static_receipt=lp6_srec,
        provenance={"paper": "PLANS/high-rate-qec.md", "version": "baseline", "locator": "table 1"},
    )
    check("LP(6,2) construction receipt verifies (reproduced, rate=[2,6])",
          bool(constructions.verify_construction_receipt(lp6_crec, static_receipt=lp6_srec))
          and lp6_crec["payload"]["ancestry"] == "reproduced"
          and lp6_crec["metrics"]["rate"] == [2, 6]
          and lp6_crec["claims"]["proved_lower"]["d"] == 2)

    pp_crec = constructions.build_construction_receipt(pp_imp, static_receipt=pp_srec)
    check("PP256 construction receipt verifies (supplied, paper_claims vs witness_upper vs proved_lower)",
          bool(constructions.verify_construction_receipt(pp_crec, static_receipt=pp_srec))
          and pp_crec["payload"]["ancestry"] == "supplied"
          and pp_crec["metrics"]["rate"] == [64, 256]
          and pp_crec["claims"]["paper_claims"]["d"] == 12
          and pp_crec["claims"]["witness_upper"]["d"] == 12
          and pp_crec["claims"]["proved_lower"]["d"] is None)

    bb_repro = constructions.build_construction_receipt(
        bb_imp,
        reproduction_spec={
            "schema": constructions.SPEC_SCHEMA,
            "kind": "named_family",
            "name": "BB144",
            "family_spec": "bb:12,6,x3+y1+y2,y3+x1+x2",
        },
    )
    check("imported BB144 upgrades ancestry to reproduced when generator matches sha256",
          bool(constructions.verify_construction_receipt(bb_repro))
          and bb_repro["payload"]["ancestry"] == "reproduced")

    bad_crec = copy.deepcopy(pp_crec)
    bad_crec["claims"]["proved_lower"]["d"] = 12
    bad_crec["static_receipt_ref"] = None
    refuses("construction receipt control: proved_lower without static_receipt_ref refused",
            lambda: constructions.verify_construction_receipt(bad_crec),
            "static_receipt_ref")


# ---------------------------------------------------------------- 5g. addressable logical bases & lift symmetries
def t_basis_and_symmetry():
    ex_dir = os.path.join(HERE, "examples")
    lp6 = constructions.lifted_product([[3]], [[3]], 3, name="LP(6,2)")
    b_lp6 = basis.derive_paired_basis(lp6)
    check("LP(6,2) derived basis verifies", bool(basis.verify_basis_receipt(lp6, b_lp6)) and b_lp6["k"] == 2)

    # Imported lower-weight / external candidate bases verify through the same verifier
    for fname, want_k in (
        ("qldpc_bb144.import.json", 12),
        ("qldpc_lp6_2.import.json", 2),
        ("pp256_q4_j3_l8_p16.import.json", 64),
    ):
        c_imp = codes.load(os.path.join(ex_dir, fname))
        b_imp = basis.build_basis_receipt(c_imp, source="imported_candidate")
        check(f"imported candidate basis {fname} verifies (k={want_k})",
              bool(basis.verify_basis_receipt(c_imp, b_imp)) and b_imp["k"] == want_k)

    # 1. Slot swap on LP(6,2): preserves code_sha256, changes basis_sha256 and physical operator for Z on q0
    b_swapped, swap_tr = basis.dress_basis(lp6, b_lp6, permute_slots=[1, 0], new_slot_names=["q0", "q1"])
    z_q0_orig = basis.embed_named_pauli(b_lp6, [("q0", "Z")])
    z_q0_swap = basis.embed_named_pauli(b_swapped, [("q0", "Z")])
    check("LP(6,2) slot swap preserves code_sha256, changes basis_sha256 and physical Z(q0)",
          b_swapped["code_sha256"] == b_lp6["code_sha256"]
          and b_swapped["basis_sha256"] != b_lp6["basis_sha256"]
          and z_q0_orig != z_q0_swap)
    refuses("LP(6,2) slot swap without conversion rejected by verify_logical_target",
            lambda: basis.verify_logical_target(lp6, b_lp6, [("q0", "Z")], z_q0_swap),
            "different logical")

    # Read-only compatibility check with ftqc_gadgets.surgery.validate_basis
    try:
        from ftqc_gadgets import surgery as gad_surg
        exp_b = basis.export_surgery_basis(b_lp6, measured_slots=["q0"])
        checks_g = [
            gad_surg.pauli(basis.pauli_from_xz(lp6.n, r, 0, 0), lp6.n) for r in lp6.hx
        ] + [
            gad_surg.pauli(basis.pauli_from_xz(lp6.n, 0, r, 0), lp6.n) for r in lp6.hz
        ]
        gad_surg.validate_basis(
            exp_b["surgery_basis"],
            checks_g,
            [gad_surg.pauli(z_q0_orig, lp6.n)],
            lp6.n,
        )
        check("exported surgery_basis accepted by ftqc_gadgets.surgery.validate_basis", True)
    except ImportError:
        skip("ftqc_gadgets.surgery.validate_basis check", "ftqc_gadgets not present")

    # 2. Valid stabilizer dressing (including mixed X/Z stabilizers creating Y entries) + tampered dressing refusal
    b_dressed, dress_tr = basis.dress_basis(
        lp6,
        b_lp6,
        dress_x_hx=[[0], [1]],
        dress_x_hz=[[1], []],
        dress_z_hx=[[], [0]],
        dress_z_hz=[[2], [1]],
    )
    check("stabilizer-dressed LP(6,2) basis and transform verify",
          bool(basis.verify_basis_receipt(lp6, b_dressed))
          and basis.verify_basis_transform(lp6, b_lp6, b_dressed, dress_tr))
    z_q0_dressed = basis.embed_named_pauli(b_dressed, [("q0", "Z")])
    check("dressed Z(q0) accepted as equivalent logical target on original basis",
          basis.verify_logical_target(lp6, b_lp6, [("q0", "Z")], z_q0_dressed) == z_q0_orig)

    bad_tr = copy.deepcopy(dress_tr)
    bad_tr["dressings"]["z_hz"][0] = [0]
    refuses("control: wrong stabilizer dressing rejected by verify_basis_transform",
            lambda: basis.verify_basis_transform(lp6, b_lp6, b_dressed, bad_tr),
            "mismatch")

    # 3. Wrong sign in a mixed-Pauli (Y) measurement
    st = codes.steane()
    b_st = basis.derive_paired_basis(st)
    y_true = basis.embed_named_pauli(b_st, [("q0", "Y")])
    y_wrong_sign = ("+" if y_true[0] == "-" else "-") + y_true[1:]
    check("Steane logical Y valid sign accepted",
          basis.verify_logical_target(st, b_st, [("q0", "Y")], y_true) == y_true)
    refuses("control: wrong sign in mixed-Pauli (Y) measurement refused",
            lambda: basis.verify_logical_target(st, b_st, [("q0", "Y")], y_wrong_sign),
            "wrong sign")
    # Also test stabilizer-dressed Y on LP(6,2)
    y_lp_dressed = basis.embed_named_pauli(b_dressed, [("q0", "Y")])
    y_lp_flipped = ("+" if y_lp_dressed[0] == "-" else "-") + y_lp_dressed[1:]
    check("LP(6,2) dressed mixed-Pauli Y accepted against base basis",
          bool(basis.verify_logical_target(lp6, b_lp6, [("q0", "Y")], y_lp_dressed)))
    refuses("control: LP(6,2) dressed mixed-Pauli Y with flipped sign refused",
            lambda: basis.verify_logical_target(lp6, b_lp6, [("q0", "Y")], y_lp_flipped),
            "wrong sign")

    # 4. Duplicated slot refusal
    dup_slots = copy.deepcopy(b_lp6["slots"])
    dup_slots[1]["name"] = dup_slots[0]["name"]
    refuses("control: duplicated slot name in basis refused",
            lambda: basis.build_basis_receipt(lp6, dup_slots),
            "duplicated")
    refuses("control: duplicated measured slot in consumer export refused",
            lambda: basis.export_surgery_basis(b_lp6, measured_slots=["q0", "q0"]),
            "duplicated")
    refuses("control: duplicated slot in embed_named_pauli refused",
            lambda: basis.embed_named_pauli(b_lp6, [("q0", "X"), ("q0", "Z")]),
            "duplicated")

    # 5. Seed-orbit coverage (Zheng et al., arXiv:2607.28605): partial vs full coverage
    z0_supp = basis.parse_pauli(b_lp6["slots"][0]["z"], lp6.n)[1]
    z1_supp = basis.parse_pauli(b_lp6["slots"][1]["z"], lp6.n)[1]
    cov_partial = symmetry.verify_seed_orbit_coverage(lp6, b_lp6, [z0_supp], 3, sector="Z")
    cov_full = symmetry.verify_seed_orbit_coverage(lp6, b_lp6, [z0_supp, z1_supp], 3, sector="Z")
    check("LP(6,2) single invariant seed covers 1/2 slots; two seeds cover 2/2",
          cov_partial["covered_rank"] == 1 and not cov_partial["full_coverage"]
          and cov_full["covered_rank"] == 2 and cov_full["full_coverage"])

    # Zheng et al. Ex. III.6: LP_5^{2x4} = [[100, 26, 4]] (a is 2x4, b = a^T is 4x2)
    a_z100 = [[1, 1, 1, 1], [1, 2, 4, 8]]
    b_z100 = [[a_z100[i][j] for i in range(2)] for j in range(4)]
    lp100 = constructions.lifted_product(a_z100, b_z100, 5, name="LP100_26")
    b_lp100 = basis.derive_paired_basis(lp100)
    left_mask_only = [
        v for v in gf2.nullspace(lp100.hx, lp100.n)
        if (v >> 80) == 0
    ]
    cov_left_all = symmetry.verify_seed_orbit_coverage(lp100, b_lp100, left_mask_only, 5, sector="Z")
    check("Zheng LP_5^{2x4} [[100, 26]]: partial seed-orbit coverage strictly < 26 detected",
          lp100.n == 100 and b_lp100["k"] == 26
          and cov_left_all["covered_rank"] == 22
          and cov_left_all["full_coverage"] is False)

    # 6. Lee odd-order Abelian spectroscopy (arXiv:2607.05386) + fail-closed even-order check
    refuses("spectral basis refuses even cyclic lift L=4 fail-closed",
            lambda: symmetry.spectral_abelian_basis(4, [[[0, 1]]], [[[0, 1]]]),
            "odd-order")
    refuses("spectral basis refuses even Abelian group Z_3 x Z_2 fail-closed",
            lambda: symmetry.spectral_abelian_basis([3, 2], [[[[0, 0]]]], [[[[0, 0]]]]),
            "odd-order")
    # Binary elimination remains available on even lifts
    even_lp = constructions.lifted_product([[3]], [[3]], 4, name="LP_even4")
    check("binary elimination remains available on even lift L=4",
          basis.derive_paired_basis(even_lp)["k"] == 2)

    # Odd cyclic lifts: LP(6,2) (L=3), nontrivial packet C_3 (1+x+x^2 -> k=4), radial90 (L=5 -> k=8), L=15
    c_sp3, b_sp3, rep_sp3 = symmetry.spectral_abelian_basis(3, [[[0, 1]]], [[[0, 1]]], name="LP(6,2)")
    check("spectral basis on LP(6,2) (C_3) matches k=2 and c6860dd462fe",
          c_sp3.sha256() == lp6.sha256() and b_sp3["k"] == 2 and rep_sp3["odd_order_checked"])

    c_pkt3, b_pkt3, rep_pkt3 = symmetry.spectral_abelian_basis(3, [[[0, 1, 2]]], [[[0, 1, 2]]])
    check("spectral basis on C_3 nontrivial degree-2 packet gives k=4 paired slots",
          b_pkt3["k"] == 4 and any(p["degree"] == 2 and p["binary_logical_slots"] == 4 for p in rep_pkt3["packets"]))

    rad_a = [[[2], [3], [4]], [[1], [4], [1]], [[4], [3], [2]]]
    rad_b = [[[2], [2], [0]], [[4], [0], [4]], [[1], [3], [0]]]
    c_rad_sp, b_rad_sp, rep_rad_sp = symmetry.spectral_abelian_basis(5, rad_a, rad_b, name="Radial90")
    check("spectral basis on radial90 (C_5) matches b656926c5959 and k=8",
          c_rad_sp.sha256().startswith("b656926c5959") and b_rad_sp["k"] == 8)

    # Lee T_{2,4}^{15} over Z_15 (n=300, k=66 across 5 cyclotomic packets of degrees 1, 2, 4)
    p15_a = [[[0], [0], [0], [0]], [[0], [1], [2], [3]]]
    p15_b = [[[0], [0]], [[0], [14]], [[0], [13]], [[0], [12]]]
    c_p15, b_p15, rep_p15 = symmetry.spectral_abelian_basis(15, p15_a, p15_b, name="T_2_4_15")
    check("spectral basis on Z_15 2x4 lift T_{2,4}^{15} (n=300, k=66) verifies across degree 1,2,4 packets",
          c_p15.n == 300 and b_p15["k"] == 66 and len(rep_p15["packets"]) == 5)

    # Noncyclic odd Abelian group G = Z_3 x Z_3 (|G| = 9, n = 36, k = 4)
    z3z3_a = [[[[0, 0], [1, 0]], [[0, 0], [0, 1]]]]
    c_z33, b_z33, rep_z33 = symmetry.spectral_abelian_basis([3, 3], z3z3_a, z3z3_a, name="Z3xZ3_LP")
    check("spectral basis on noncyclic Z_3 x Z_3 lift verifies against binary elimination",
          c_z33.n == 36 and b_z33["k"] == 4 == rep_z33["binary_crosscheck_k"])

    # 7. Group-extension lift chain maps (Hirasaki & Lee, arXiv:2607.28621): even vs odd fiber extensions
    src_c, tgt_c, src_b, tgt_b, maps_pb2, cert_pb2 = symmetry.build_cyclic_extension_chain_map(
        [[3]], [[3]], base_lift=3, ext_factor=2, mode="pullback"
    )
    _, _, _, _, _, cert_pf2 = symmetry.build_cyclic_extension_chain_map(
        [[3]], [[3]], base_lift=3, ext_factor=2, mode="pushforward"
    )
    _, _, _, _, _, cert_pb3 = symmetry.build_cyclic_extension_chain_map(
        [[3]], [[3]], base_lift=3, ext_factor=3, mode="pullback"
    )
    _, _, _, _, _, cert_pf3 = symmetry.build_cyclic_extension_chain_map(
        [[3]], [[3]], base_lift=3, ext_factor=3, mode="pushforward"
    )
    check("cyclic extension chain maps Z_3<->Z_6 (rank 1) and Z_3<->Z_9 (rank 2) commute with boundaries and satisfy M_X == M_Z^T",
          cert_pb2["commutes_boundary_1"] and cert_pb2["commutes_boundary_2"] and cert_pb2["rank_mz"] == 1
          and cert_pf2["rank_mz"] == 1
          and cert_pb3["rank_mz"] == 2 and cert_pf3["rank_mz"] == 2)
    bad_f1 = list(maps_pb2[0])
    bad_f1[0] ^= 1
    refuses("control: corrupted lift chain map fails boundary commutation",
            lambda: symmetry.verify_chain_map(src_c, tgt_c, bad_f1, maps_pb2[1], maps_pb2[2]),
            "commutation")


# ---------------------------------------------------------------- 6. SAT leg (needs cadical)
def t_sat(heavy):
    if not sat.find_cadical():
        skip("SAT leg", "cadical not found (PATH or FTQC_CODES_CADICAL)")
        return
    specs = ["steane", "golay", "bb18", "toric:3", "rotated_surface:3"]
    if heavy:
        specs += ["bb54", "gb54", "bb72", "qr:47"]
    for spec in specs:
        code = codes.load(spec)
        rec = certify.build(code, max_combos=3_000_000)
        with tempfile.TemporaryDirectory() as wd:
            sat.run_leg(rec, code, workdir=wd, timeout=1200)
            proved = {e["quantity"]: e for e in rec["evidence"]
                      if e["leg"] == "B" and e["status"] == "proved"}
            check(f"sat {spec}: all four bounds proved with checked LRAT",
                  set(proved) == set(sat.QUANTITIES),
                  str([(e["quantity"], e.get("reason")) for e in rec["evidence"] if e["leg"] == "B"]))
            for q, e in proved.items():
                a = next(x for x in rec["evidence"] if x["leg"] == "A" and x["quantity"] == q)
                check(f"sat {spec} {q}: agrees with enumeration", e["claim"] == a["claim"])
            # one above the truth must come back SAT with a verified word
            q = "d_x"
            d = rec["facts"][q]["lower"]
            r2 = certify.build(code, legs=())
            sat.run_leg(r2, code, targets={q: d + 1}, workdir=wd, timeout=1200)
            ev = next(e for e in r2["evidence"] if e["quantity"] == q)
            check(f"sat {spec}: target d+1 is refuted with a witness",
                  ev["status"] == "refuted" and r2["facts"][q]["upper"] == d)
            try:
                certify.verify(rec)
                check(f"sat {spec}: receipt verifies", True)
            except certify.Refusal as e:
                check(f"sat {spec}: receipt verifies", False, str(e))
    # a corrupted proof must be rejected by the checker
    code = codes.load("golay")
    F = cnf.for_quantity(code, "d_x", 6)
    with tempfile.TemporaryDirectory() as wd:
        res = sat.solve(F, wd, "g", 300)
        check("sat golay raw solve unsat", res[0] == "unsat")
        proof = open(os.path.join(wd, "g.lrat"), "rb").read().splitlines()
        adds = {st[1]: st for st in lrat.steps(b"\n".join(proof) + b"\n") if st[0] == "a"}
        ncl = len(F.clauses)
        final = next(st for st in adds.values() if not st[2])
        used = [h for h in final[3] if h > ncl]
        check("sat control: refutation uses derived lemmas", bool(used))
        victim = str(used[0]).encode() if used else b"-"
        cut = b"\n".join(l for l in proof if l.split()[:1] != [victim]) + b"\n"
        try:
            lrat.check(F.dimacs(), cut)
            check("sat control: proof with a line removed", False, "accepted")
        except lrat.ProofError:
            check("sat control: proof with a line removed", True)
        weaker = cnf.for_quantity(code, "d_x", 7)     # satisfiable formula, same proof
        try:
            lrat.check(weaker.dimacs(), b"\n".join(proof) + b"\n")
            check("sat control: proof replayed against another CNF", False, "accepted")
        except lrat.ProofError:
            check("sat control: proof replayed against another CNF", True)


# ---------------------------------------------------------------- 7. CLI
def t_cli():
    env = dict(os.environ, PYTHONPATH=os.path.dirname(HERE))
    with tempfile.TemporaryDirectory() as wd:
        out = os.path.join(wd, "r.json")

        def run(*args):
            return subprocess.run([sys.executable, "-m", "ftqc_codes", *args], env=env, cwd=wd,
                                  capture_output=True, text=True)
        p = run("certify", "golay", "-o", out)
        check("cli certify exits 0", p.returncode == 0 and "d=7" in p.stdout, p.stderr[-200:])
        p = run("verify", out, "--full")
        check("cli verify --full exits 0", p.returncode == 0, p.stderr[-200:])
        rec = json.load(open(out))
        rec["facts"]["d_x"]["lower"] = 9
        json.dump(rec, open(out, "w"))
        p = run("verify", out)
        check("cli verify refuses a tampered receipt with exit 1",
              p.returncode == 1 and "REFUSED" in p.stderr)
        p = run("certify", "nosuchfamily")
        check("cli unknown family exits 2", p.returncode == 2)
        spec = os.path.join(wd, "bad.json")
        json.dump(dict(name="bad", n=4, hx=["1100"], hz=["0110"]), open(spec, "w"))
        p = run("certify", spec)
        check("cli non-commuting spec exits 1", p.returncode == 1)
        p = run("spec", "steane", "-o", os.path.join(wd, "s.json"))
        p2 = run("certify", os.path.join(wd, "s.json"))
        check("cli spec file round trip", p.returncode == 0 and p2.returncode == 0 and "d=3" in p2.stdout)
        # Tagged construction JSON + CLI construct, basis & verify
        cspec = os.path.join(wd, "lp6.spec.json")
        json.dump({
            "schema": constructions.SPEC_SCHEMA,
            "kind": "lifted_product",
            "name": "LP(6,2)",
            "lift": 3,
            "a": [[3]],
            "b": [[3]],
        }, open(cspec, "w"))
        srec_out = os.path.join(wd, "lp6.static.json")
        crec_out = os.path.join(wd, "lp6.construct.json")
        brec_out = os.path.join(wd, "lp6.basis.json")
        p_s = run("certify", cspec, "-o", srec_out)
        p_c = run("construct", cspec, "--static-receipt", srec_out, "-o", crec_out)
        p_v = run("verify", crec_out, "--static-receipt", srec_out)
        p_b = run("basis", cspec, "--method", "spectral", "-o", brec_out)
        p_vb = run("verify", brec_out, "--code", cspec)
        check("cli tagged construction + construct + spectral basis + verify exits 0",
              p_s.returncode == 0 and p_c.returncode == 0 and p_v.returncode == 0
              and p_b.returncode == 0 and p_vb.returncode == 0
              and "c6860dd462fe" in p_c.stdout, p_b.stderr[-200:] + p_vb.stderr[-200:])
        for ex in ("steane.receipt.json", "golay.receipt.json"):
            ex_rec = json.load(open(os.path.join(HERE, "examples", ex)))
            check(f"banked example {ex} verifies", bool(certify.verify(ex_rec, full=True)))
        import qeccert as _legacy_shim
        from qeccert.codes import CSSCode as _LegacyCSSCode
        check("legacy qeccert shim resolves to ftqc_codes",
              _legacy_shim.VERSION == certify.VERSION and _LegacyCSSCode is codes.CSSCode)


# ---------------------------------------------------------------- 8. High-rate paper families & Leg G orbit search
def t_highrate_families_and_orbit_leg():
    import itertools
    from ftqc_codes import families_highrate, orbit_minweight

    # Row 1: Trivial lift L=1 HGP equality + finite-lift Raveendran LP witnesses/bounds
    c_hgp, c_lp1 = families_highrate.trivial_lift_hgp_equivalence(
        [0b011, 0b110], 3, [0b0011, 0b0110, 0b1100], 4, name="HGP_L1"
    )
    check("trivial lift L=1 LP equals HGP byte-for-byte", c_hgp.sha256() == c_lp1.sha256() and c_hgp.n == 18)

    rav4 = families_highrate.raveendran_lp_instance(4)
    rav5 = families_highrate.raveendran_lp_instance(5)
    rec_rav4 = certify.build(rav4)
    rec_rav5 = certify.build(rav5)
    check("Raveendran finite LP L=4,5 exact receipts verify",
          bool(certify.verify(rec_rav4)) and bool(certify.verify(rec_rav5))
          and rec_rav4["facts"]["d"]["exact"] and rec_rav5["facts"]["d"]["exact"])

    # Row 2: Two-Block Group Algebra (2BGA) over nonabelian S_3 + exact BB reconciliation
    s3_elems = list(itertools.permutations(range(3)))
    s3_idx = {p: i for i, p in enumerate(s3_elems)}
    s3_left = [[s3_idx[tuple(g[h[i]] for i in range(3))] for h in s3_elems] for g in s3_elems]
    s3_right = [[s3_idx[tuple(h[g[i]] for i in range(3))] for h in s3_elems] for g in s3_elems]
    ga_s3 = families_highrate.two_block_ga(
        6, s3_left, s3_right,
        [s3_idx[(0, 1, 2)], s3_idx[(1, 2, 0)], s3_idx[(1, 0, 2)]],
        [s3_idx[(0, 1, 2)], s3_idx[(2, 0, 1)], s3_idx[(0, 2, 1)]],
        name="2BGA_S3",
    )
    check("nonabelian S_3 2BGA commutes and is nonabelian",
          ga_s3.commutes() and ga_s3.n == 12 and ga_s3.meta["two_block_ga"]["abelian"] is False)

    for bb_name, args in (
        ("BB18", (3, 3, ["x0", "x1", "y1"], ["x0", "x2", "y2"])),
        ("GB54", (27, 1, ["x8", "x13", "x15", "x16"], ["x7", "x8", "x10", "x15"])),
        ("BB144", (12, 6, ["x3", "y1", "y2"], ["y3", "x1", "x2"])),
    ):
        ga_c, cmap = families_highrate.reconcile_bb_as_2bga(*args, name=bb_name)
        check(f"2BGA reconciles {bb_name} canonical + inverted-shift coordinate map",
              ga_c.sha256() == codes.load(bb_name.lower()).sha256() and len(cmap["perm_q"]) == ga_c.n)

    # Row 3: Nonabelian balanced product + Quantum Tanner codes (D_8 [[36,8,3]] and S_3 [[288,8]])
    # Left-Cayley checks {g, r*g} are equivariant under RIGHT G-action; right-Cayley checks {g, g*r} under LEFT G-action!
    r_s3 = s3_idx[(1, 2, 0)]
    hx_a = [gf2.from_supp([g, s3_left[r_s3][g]]) for g in range(6)]
    hz_b = [gf2.from_supp([g, s3_right[r_s3][g]]) for g in range(6)]
    bp_s3 = families_highrate.balanced_product(
        6, s3_left, s3_right,
        hx_a, 6, s3_right, s3_right,
        hz_b, 6, s3_left, s3_left,
        name="BP_S3",
    )
    check("nonabelian S_3 balanced product commutes with free 6-element orbits",
          bp_s3.commutes() and bp_s3.n == 12 and bp_s3.meta["balanced_product"]["abelian"] is False)

    # Negative control: non-equivariant H_A in balanced_product refused
    bad_hx_a = list(hx_a)
    bad_hx_a[0] ^= 1
    try:
        families_highrate.balanced_product(
            6, s3_left, s3_right, bad_hx_a, 6, s3_right, s3_right, hz_b, 6, s3_left, s3_left
        )
        check("non-equivariant balanced_product refused", False)
    except certify.Refusal:
        check("non-equivariant balanced_product refused", True)

    qt_d8 = families_highrate.quantum_expanders_d8_36_8_3()
    rec_d8 = certify.build(qt_d8)
    hyps = qt_d8.meta["tanner"]["hypotheses"]
    check("Quantum Tanner D_8 [[36,8,3]] exact parameters + checked/assumed hypotheses",
          rec_d8["facts"]["n"] == 36 and rec_d8["facts"]["k"] == 8
          and rec_d8["facts"]["d"] == {"lower": 3, "upper": 3, "exact": True}
          and hyps["total_non_conjugacy"] == "checked"
          and hyps["spectral_expansion_ramanujan"] == "assumed")

    d8_imp_path = os.path.join(HERE, "examples", "quantum_expanders_d8_36_8_3.import.json")
    if os.path.isfile(d8_imp_path):
        d8_imp = codes.load(d8_imp_path)
        check("banked quantum_expanders_d8_36_8_3.import.json matches native generator",
              d8_imp.sha256() == qt_d8.sha256() and "tanner" in d8_imp.meta)

    # Negative control: TNC violation refused in quantum_tanner_lrcc
    try:
        # Taking gen_b == gen_a guarantees a*e == e*a (TNC violation at g=e)
        d8_l = [[0] * 8 for _ in range(8)]
        # Reuse valid S_3 with overlapping generator!
        t_s3 = s3_idx[(1, 0, 2)]
        families_highrate.quantum_tanner_lrcc(
            6, s3_left, s3_right, [t_s3], [t_s3], [[1]], [[1]], [[1]], [[1]]
        )
        check("TNC violation refused in quantum_tanner_lrcc", False)
    except certify.Refusal:
        check("TNC violation refused in quantum_tanner_lrcc", True)

    qt_s3 = families_highrate.quantum_expanders_s3_288_8()
    st_s3 = certify.structure(qt_s3)
    check("Lifted Quantum Tanner S_3 [[288,8]] commutes with max check weight 12",
          st_s3["n"] == 288 and st_s3["k"] == 8
          and max(gf2.wt(r) for r in qt_s3.hx + qt_s3.hz) == 12
          and qt_s3.meta["paper_claims"]["d_upper"] == 15)

    # Row 4: Pair-Partition CPM & APM high-rate codes + construction receipts
    for fn, want_n, want_k, want_d_claim in (
        (families_highrate.pp184_50, 184, 50, 10),
        (families_highrate.pp256_68, 256, 68, 12),
        (families_highrate.pp576_148, 576, 148, 18),
        (families_highrate.apm1152_580, 1152, 580, 12),
    ):
        c_hr = fn()
        crec_hr = constructions.build_construction_receipt(c_hr)
        constructions.verify_construction_receipt(crec_hr)
        check(f"{c_hr.name} reproduces [[{want_n},{want_k}]] with paper_claims d={want_d_claim} separated from proved_lower",
              crec_hr["metrics"]["n"] == want_n and crec_hr["metrics"]["k"] == want_k
              and crec_hr["payload"]["ancestry"] == "reproduced"
              and crec_hr["claims"]["paper_claims"].get("d") == want_d_claim
              and crec_hr["claims"]["proved_lower"]["d"] is None)

    # Negative controls for cpm_pair_partition: non-CSS refused & broken pair-partition refused
    try:
        families_highrate.cpm_pair_partition(8, [[0, 1]], [[0, 1]], code_type="symplectic")
        check("non-CSS cpm_pair_partition refused", False)
    except certify.Refusal:
        check("non-CSS cpm_pair_partition refused", True)

    try:
        families_highrate.cpm_pair_partition(8, [[0, 1, 2, 3]], [[0, 2, 5, 7]])
        check("broken pair-partition shift array refused", False)
    except certify.Refusal:
        check("broken pair-partition shift array refused", True)

    # Leg G: orbit_minweight exact cross-check against Leg A + cheap & full verification + negative controls
    for spec_name, want_d in (("lp6_2", 2), ("bb18", 4), ("toric:4", 4)):
        c_g = codes.load(spec_name)
        gens_g = orbit_minweight.infer_default_automorphisms(c_g)
        if gens_g is None and spec_name == "toric:4":
            # Toric(L=4) has n = 32 = 16 + 16; shift x by 1 on both 4x4 blocks!
            sh4 = [((i // 4 + 1) % 4) * 4 + (i % 4) for i in range(16)]
            gens_g = [sh4 + [16 + x for x in sh4]]
        rec_a = certify.build(c_g, legs=("A",))
        rec_g = certify.build(c_g, legs=())
        orbit_minweight.run_leg(rec_g, c_g, aut_generators=gens_g)
        check(f"Leg G matches unrestricted Leg A on {spec_name} (d={want_d})",
              rec_g["facts"]["d"] == rec_a["facts"]["d"]
              and rec_g["facts"]["d"]["lower"] == want_d
              and bool(certify.verify(rec_g, full=False))
              and bool(certify.verify(rec_g, full=True)))

    # Negative controls on Leg G: tampered word_to_root and non-automorphism generator refused
    c_bb18 = codes.load("bb18")
    rec_bb18_g = certify.build(c_bb18, legs=())
    orbit_minweight.run_leg(rec_bb18_g, c_bb18)
    bad_rec_g = json.loads(json.dumps(rec_bb18_g))
    for ev in bad_rec_g["evidence"]:
        if ev.get("leg") == "G":
            # Corrupt a non-root word_to_root witness so it lands on the wrong qubit
            orb0 = ev["orbits"][0]
            non_root = [m for m in orb0["members"] if m != orb0["root"]][0]
            orb0["words_to_root"][str(non_root)] = []
            break
    try:
        certify.verify(bad_rec_g, full=False)
        check("tampered Leg G word_to_root witness refused", False)
    except certify.Refusal:
        check("tampered Leg G word_to_root witness refused", True)

    bad_aut = list(range(c_bb18.n))
    bad_aut[0], bad_aut[1] = bad_aut[1], bad_aut[0]
    try:
        orbit_minweight.orbit_min_weight(c_bb18, "d_x", [bad_aut], target=2)
        check("non-automorphism generator refused by Leg G", False)
    except certify.Refusal:
        check("non-automorphism generator refused by Leg G", True)


# ---------------------------------------------------------------- 9. Lean leg (needs lake + Lean-QEC)
def t_lean(enabled):
    try:
        from ftqc_codes import leangen
    except ImportError:
        skip("Lean leg", "leangen module not present")
        return
    code = codes.load("steane")
    src = leangen.render(code, dict(d_x=3, d_z=3), "FtqcCodesSteane")
    check("lean generator emits the theorem", "theorem FtqcCodesSteane_dist" in src
          and "bitvec_sat_translation_correct" in src and "sorry" not in src)
    if not enabled:
        skip("Lean leg build", "not requested (run with --lean)")
        return
    why = leangen.unavailable()
    if why:
        skip("Lean leg build", why)
        return
    for spec, d in (("steane", 3), ("bb18", 4)):
        code = codes.load(spec)
        rec = certify.build(code)
        with tempfile.TemporaryDirectory() as wd:
            leangen.run_leg(rec, code, workdir=wd, timeout=3600)
        evs = [e for e in rec["evidence"] if e["leg"] == "C"]
        check(f"lean {spec}: kernel-checked d >= {d}",
              any(e["status"] == "proved" and e["claim"] == f"d_x >= {d}" for e in evs),
              str([(e["status"], e.get("reason", "")[-300:]) for e in evs]))
    # negative control: a false distance must fail to build
    code = codes.load("steane")
    rec = certify.build(code, legs=())
    with tempfile.TemporaryDirectory() as wd:
        leangen.run_leg(rec, code, targets=dict(d_x=4, d_z=4), workdir=wd, timeout=3600)
    evs = [e for e in rec["evidence"] if e["leg"] == "C"]
    check("lean control: false claim d >= 4 for Steane does not build",
          evs and all(e["status"] != "proved" for e in evs))


def t_code_actions():
    from ftqc_codes import automorphisms as aut
    code = codes.load("steane")
    paired = basis.derive_paired_basis(code)
    steps = [{"layer": [["H", q] for q in range(7)]}]
    receipt = aut.build(code, paired, steps)
    check("transversal H logical Z image", receipt["logical_z"][0]["x_mask"] == 1)
    check("logical H reverses Y sign", aut.logical_image(receipt, (1, 1, 0)) == (1, 1, 1))
    changed = copy.deepcopy(receipt)
    changed["logical_z"][0]["sign"] ^= 1
    refuses("altered logical sign", lambda: aut.verify(code, changed))
    refuses("wrong physical sign", lambda: aut.build(code, paired, [{"layer": [["X", 0]]}]))
    spec_err_or_refuses("fold collision", lambda: aut.build(code, paired, [{"layer": [["CZ", 0, 1], ["S", 1]]}]))
    code = codes.load("bb18")
    paired = basis.derive_paired_basis(code)
    from ftqc_codes.orbit_minweight import infer_default_automorphisms
    perm = infer_default_automorphisms(code)[0]
    action = aut.build(code, paired, [{"perm": perm}])
    reference = symmetry.induced_logical_action(code, paired, perm)
    check("pure permutation agrees with original verifier", all(action[k] == v for k,v in reference.items()))
    inverse = [perm.index(q) for q in range(code.n)]
    identity = aut.build(code, paired, [{"perm": perm}, {"perm": inverse}])
    check("composition with inverse", all(row["x_mask"] == 1 << j for j, row in enumerate(identity["logical_x"])))


def main():
    heavy, lean = "--heavy" in sys.argv, "--lean" in sys.argv
    t0 = time.monotonic()
    for fn, args in (
        (t_gf2, ()),
        (t_minweight, ()),
        (t_families, ()),
        (t_controls, ()),
        (t_cnf, ()),
        (t_lrat, ()),
        (t_bch, ()),
        (t_cuts, ()),
        (t_concat, ()),
        (t_hier, ()),
        (t_interchange_and_constructions, ()),
        (t_basis_and_symmetry, ()),
        (t_code_actions, ()),
        (t_highrate_families_and_orbit_leg, ()),
        (t_sat, (heavy,)),
        (t_cli, ()),
        (t_lean, (lean,)),
    ):
        t = time.monotonic()
        fn(*args)
        print(f"      {fn.__name__:32s} {time.monotonic() - t:6.2f} s")
    print(f"ftqc_codes selftest: {len(PASS)} passed, {len(FAIL)} failed, {len(SKIP)} skipped "
          f"in {time.monotonic() - t0:.1f} s")
    if FAIL:
        print("FAILED: " + "; ".join(sorted(set(FAIL))[:12]))
        return 1
    print("ALL PASS" + (f" ({len(SKIP)} named skips)" if SKIP else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
