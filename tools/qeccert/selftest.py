#!/usr/bin/env python3
"""qeccert battery.  Default tier is stdlib-only and finishes in a few seconds.

    python3 tools/qeccert/selftest.py            # default tier
    python3 tools/qeccert/selftest.py --heavy    # adds QR-47 by SAT, BB72/BB90 targets
    python3 tools/qeccert/selftest.py --lean     # adds the Lean-QEC leg (needs lake + a built Lean-QEC)

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

from qeccert import certify, cnf, codes, gf2, lrat, minweight, sat  # noqa: E402

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
        basis = gf2.Echelon(gens).basis
        sub = [b for b in basis if rnd.random() < 0.4]
        if rnd.random() < 0.3 and len(basis) >= 2:
            sub.append(basis[0] ^ basis[1])
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


# ---------------------------------------------------------------- 6. SAT leg (needs cadical)
def t_sat(heavy):
    if not sat.find_cadical():
        skip("SAT leg", "cadical not found (PATH or QECCERT_CADICAL)")
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
        # remove a lemma that the refutation really uses: a hint of the empty clause
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
            return subprocess.run([sys.executable, "-m", "qeccert", *args], env=env, cwd=wd,
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


# ---------------------------------------------------------------- 8. Lean leg (needs lake + Lean-QEC)
def t_lean(enabled):
    try:
        from qeccert import leangen
    except ImportError:
        skip("Lean leg", "leangen module not present")
        return
    code = codes.load("steane")
    src = leangen.render(code, dict(d_x=3, d_z=3), "QeccertSteane")
    check("lean generator emits the theorem", "theorem QeccertSteane_dist" in src
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


def main():
    heavy, lean = "--heavy" in sys.argv, "--lean" in sys.argv
    t0 = time.monotonic()
    for fn, args in ((t_gf2, ()), (t_minweight, ()), (t_families, ()), (t_controls, ()),
                     (t_cnf, ()), (t_lrat, ()), (t_sat, (heavy,)), (t_cli, ()), (t_lean, (lean,))):
        t = time.monotonic()
        fn(*args)
        print(f"      {fn.__name__:12s} {time.monotonic() - t:6.2f} s")
    print(f"qeccert selftest: {len(PASS)} passed, {len(FAIL)} failed, {len(SKIP)} skipped "
          f"in {time.monotonic() - t0:.1f} s")
    if FAIL:
        print("FAILED: " + "; ".join(sorted(set(FAIL))[:12]))
        return 1
    print("ALL PASS" + (f" ({len(SKIP)} named skips)" if SKIP else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
