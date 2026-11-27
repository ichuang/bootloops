#!/usr/bin/env python3
"""ftqc_audit battery.  Default tier is stdlib-only (plus ftqc_codes from this repository).

    python3 tools/ftqc_audit/selftest.py
    python3 tools/ftqc_audit/selftest.py --heavy     # adds surface d = 5 and Golay preparation at order 2

The SAT leg needs cadical and the Stim leg needs the stim Python package; each skips
with a named line when its engine is absent.  Exit 0 only if every check that ran passed.
"""
import copy
import itertools
import json
import os
import random
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, TOOLS)

from ftqc_faults.audit import analysis, certify, library, stab, stimleg  # noqa: E402
from ftqc_faults.audit.circuit import Circuit, CircuitError  # noqa: E402
from ftqc_faults.audit.model import Model  # noqa: E402

PASS, FAIL, SKIP = [], [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    if not cond:
        print(f"FAIL  {name}  {detail}")


def refuses(name, fn, needle="", exc=certify.Refusal):
    try:
        fn()
    except exc as e:
        check(name, needle in str(e), f"refused for another reason: {e}")
        return
    check(name, False, "was accepted")


def skip(name, why):
    SKIP.append(name)
    print(f"SKIP  {name}: {why}")


# ---------------------------------------------------------------- dense reference simulator
class Dense:
    """State-vector simulation for a handful of qubits (qubit q = bit q of the index)."""

    def __init__(self, n):
        self.n = n
        self.v = [0j] * (1 << n)
        self.v[0] = 1 + 0j

    def h(self, q):
        b = 1 << q
        r = 0.5 ** 0.5
        for i in range(1 << self.n):
            if not i & b:
                a, c = self.v[i], self.v[i | b]
                self.v[i], self.v[i | b] = r * (a + c), r * (a - c)

    def cx(self, c, t):
        bc, bt = 1 << c, 1 << t
        for i in range(1 << self.n):
            if i & bc and not i & bt:
                self.v[i], self.v[i | bt] = self.v[i | bt], self.v[i]

    def pauli(self, q, p):
        b = 1 << q
        if p in "XY":
            for i in range(1 << self.n):
                if not i & b:
                    self.v[i], self.v[i | b] = self.v[i | b], self.v[i]
        if p in "ZY":
            for i in range(1 << self.n):
                if i & b:
                    self.v[i] = -self.v[i]

    def project(self, q, basis, outcome):
        """Project onto the given outcome; returns the squared norm kept."""
        if basis == "X":
            self.h(q)
        b = 1 << q
        norm = 0.0
        for i in range(1 << self.n):
            if bool(i & b) != bool(outcome):
                self.v[i] = 0j
            else:
                norm += abs(self.v[i]) ** 2
        if norm > 1e-12:
            s = norm ** -0.5
            self.v = [a * s for a in self.v]
        if basis == "X":
            self.h(q)
        return norm

    def reduced(self, keep):
        """Density matrix on the qubits in `keep` (list), as a dict of entries."""
        rest = [q for q in range(self.n) if q not in keep]
        rho = {}
        for i in range(1 << self.n):
            if self.v[i] == 0:
                continue
            ki = tuple((i >> q) & 1 for q in keep)
            ri = tuple((i >> q) & 1 for q in rest)
            for j in range(1 << self.n):
                if self.v[j] == 0 or tuple((j >> q) & 1 for q in rest) != ri:
                    continue
                kj = tuple((j >> q) & 1 for q in keep)
                rho[ki, kj] = rho.get((ki, kj), 0j) + self.v[i] * self.v[j].conjugate()
        return rho


def dense_run(circ, faults, forced=None):
    """All branches of a run.  faults: {(tick, op index): label}.  forced: reported
    outcomes to project onto (the believed run, no faults).  Yields
    (probability, reported outcomes, Dense state)."""
    def rec(k, j, sim, prob, reported):
        while k < len(circ.ticks) and j >= len(circ.ticks[k]):
            k, j = k + 1, 0
        if k == len(circ.ticks):
            yield prob, reported, sim
            return
        op = circ.ticks[k][j]
        lab = faults.get((k, j))
        if op.name in ("M", "MX"):
            basis = "Z" if op.name == "M" else "X"
            choices = (0, 1) if forced is None else (forced[len(reported)],)
            for out in choices:
                s2 = Dense(sim.n)
                s2.v = list(sim.v)
                p = s2.project(op.qs[0], basis, out)
                if p < 1e-9:
                    continue
                rep = out ^ (1 if lab == "flip" else 0)
                yield from rec(k, j + 1, s2, prob * p, reported + [rep])
            return
        if op.name == "RX":
            sim.h(op.qs[0])                      # preparations happen on fresh |0>
        elif op.name == "H":
            sim.h(op.qs[0])
        elif op.name in ("S", "SDG"):
            phase = 1j if op.name == "S" else -1j
            for i in range(len(sim.v)):
                if i >> op.qs[0] & 1:
                    sim.v[i] *= phase
        elif op.name == "Z_IF" and reported[op.record]:
            sim.pauli(op.qs[0], "Z")
        elif op.name == "CX":
            sim.cx(*op.qs)
        if lab:
            for q, p in zip(op.qs, lab):
                if p != "I":
                    sim.pauli(q, p)
        yield from rec(k, j + 1, sim, prob, reported)
    yield from rec(0, 0, Dense(circ.n), 1.0, [])


def same_state(r1, r2):
    keys = set(r1) | set(r2)
    return all(abs(r1.get(k, 0) - r2.get(k, 0)) < 1e-7 for k in keys)


def random_circuit(rnd, n):
    ticks = [[("RX" if rnd.random() < 0.5 else "R", q) for q in range(n)]]
    for _ in range(rnd.randint(2, 5)):
        qs = list(range(n))
        rnd.shuffle(qs)
        tick = []
        while qs:
            if len(qs) >= 2 and rnd.random() < 0.6:
                tick.append(("CX", qs.pop(), qs.pop()))
            else:
                q = qs.pop()
                tick.append((rnd.choice(["H", "I", "S", "SDG"]), q))
        ticks.append(tick)
    meas = rnd.sample(range(n), rnd.randint(1, n - 1))
    ticks.append([(rnd.choice(["M", "MX"]), q) for q in meas])
    live = [q for q in range(n) if q not in meas]
    ticks.append([("Z_IF", rnd.choice(live), rnd.randrange(len(meas)))])
    return Circuit(ticks, n)


# ---------------------------------------------------------------- 1. stabilizer simulator
def t_stab():
    ghz = Circuit([[("RX", 0), ("R", 1), ("R", 2)], [("CX", 0, 1)], [("CX", 1, 2)],
                   [("M", 0), ("M", 1), ("M", 2)]])
    sim = stab.Sim(3).run(ghz)
    check("stab GHZ: one random outcome, two detectors",
          [r for _, r in sim.outcomes] == [True, False, False]
          and sim.detectors() == [([0, 1], 0), ([0, 2], 0)])
    bell = Circuit([[("RX", 0), ("R", 1)], [("CX", 0, 1)], [("MX", 0), ("MX", 1)]])
    check("stab Bell in X basis: XX = +1", stab.Sim(2).run(bell).detectors() == [([0, 1], 0)])
    minus = Circuit([[("RX", 0)], [("H", 0)], [("M", 0)]])
    check("stab |+> then H: deterministic 0", stab.Sim(1).run(minus).detectors() == [([0], 0)])
    reuse = Circuit([[("RX", 0), ("R", 1)], [("CX", 0, 1)], [("M", 1)], [("R", 1)], [("CX", 0, 1)],
                     [("M", 1)]])
    check("stab reset and reuse: second outcome repeats the first",
          stab.Sim(2).run(reuse).detectors() == [([0, 1], 0)])


# ---------------------------------------------------------------- 2. model against dense simulation
def t_dense():
    rnd = random.Random(11)
    n_circ = n_sets = 0
    seen = dict(harm_undetected=0, harmless_undetected=0, detected=0, frame_only=0)
    while n_circ < 60:
        n = rnd.randint(2, 4)
        circ = random_circuit(rnd, n)
        model = Model(circ)
        sim = model.sim
        expect = {tuple(d): e for d, e in sim.detectors()}
        n_circ += 1
        locs = model.locations
        pos = {loc.index: (k, j) for loc, (k, j, _) in zip(locs, circ.locations())}
        meas_locs = [l for l in locs if l.op.name in ("M", "MX")]
        for trial in range(8):
            pool = meas_locs if trial >= 5 else locs       # some sets are outcome flips only
            chosen = rnd.sample(pool, rnd.randint(0 if trial < 5 else 1, min(3, len(pool))))
            faults = [(loc.index, rnd.choice(sorted(loc.faults))) for loc in chosen]
            det, obs, syn = model.replay(faults)
            fmap = {pos[i]: lab for i, lab in faults}
            n_sets += 1
            for prob, reported, state in dense_run(circ, fmap):
                fired = 0
                for b, d in enumerate(model.detectors):
                    if (sum(reported[i] for i in d) & 1) != expect[tuple(d)]:
                        fired |= 1 << b
                check("dense: detector pattern equals the model's", fired == det,
                      f"{circ.to_json()['ticks']} {faults}")
                believed = list(dense_run(circ, {}, forced=reported))
                check("dense: reported outcomes are possible iff no detector fires",
                      bool(believed) == (det == 0), f"{circ.to_json()['ticks']} {faults}")
                seen["detected"] += det != 0
                if believed and circ.outputs:
                    same = same_state(state.reduced(circ.outputs), believed[0][2].reduced(circ.outputs))
                    seen["harm_undetected" if syn else "harmless_undetected"] += 1
                    seen["frame_only"] += bool(syn) and all(lab == "flip" for _, lab in faults)
                    check("dense: output state is the believed one iff the syndrome is zero",
                          same == (syn == 0), f"{circ.to_json()['ticks']} {faults} syn={syn}")
    check("dense: enough cases exercised", n_sets >= 300 and all(v >= 20 for v in seen.values()),
          str(seen))


# ---------------------------------------------------------------- 3. linearity and columns
def t_columns():
    rnd = random.Random(3)
    for spec in ("cat:4,verify", "css_prep:steane", "surface:3,1"):
        model = Model(library.load(spec))
        cols = set(model.columns)
        for _ in range(200):
            chosen = rnd.sample(model.locations, rnd.randint(1, 3))
            faults = [(l.index, rnd.choice(sorted(l.faults))) for l in chosen]
            total = [0, 0, 0]
            for f in faults:
                s = model.replay([f])
                check(f"columns {spec}: every single-fault signature is a column or zero",
                      s == (0, 0, 0) or s in cols)
                total = [a ^ b for a, b in zip(total, s)]
            check(f"columns {spec}: effect of a set is the XOR of its members",
                  tuple(total) == model.replay(faults))
        for loc in model.locations:
            sigs = {model.signature(*e) for e in loc.faults.values()} | {(0, 0, 0)}
            closed = all(tuple(a ^ b for a, b in zip(s, t)) in sigs for s in sigs for t in sigs)
            check(f"columns {spec}: a location's signatures are closed under addition",
                  closed or loc.op.name in ("M", "MX"))


# ---------------------------------------------------------------- 4. search against brute force
def t_search():
    rnd = random.Random(4)
    for trial in range(200):
        m = rnd.randint(1, 9)
        cols = sorted({(rnd.getrandbits(3), rnd.getrandbits(2)) for _ in range(m)} - {(0, 0)})
        if not cols:
            continue
        want = None
        for j in range(1, len(cols) + 1):
            for idx in itertools.combinations(range(len(cols)), j):
                d = o = 0
                for i in idx:
                    d ^= cols[i][0]
                    o ^= cols[i][1]
                if d == 0 and o:
                    want = j
                    break
            if want:
                break
        got = analysis.min_undetected(cols, lambda p: any(p), 6)
        if want is not None and want <= 6:
            check("search equals brute force", got["exact"] and got["lower"] == want, f"trial {trial}")
        else:
            check("search reports no witness when none exists", got["upper"] is None and got["lower"] == 7)
        capped = analysis.min_undetected(cols, lambda p: any(p), 6, cap=rnd.randint(1, 12))
        check("search with a cap never exceeds the truth",
              want is None or capped["lower"] <= want, f"trial {trial}")


def brute_strict(model, order):
    """Strict fault tolerance straight from the definition: all location subsets, all labels."""
    for j in range(1, order + 1):
        ok = analysis.weight_table(model, j)
        for locs in itertools.combinations(model.locations, j):
            for labels in itertools.product(*[sorted(l.faults) for l in locs]):
                det, obs, syn = model.replay([(l.index, lab) for l, lab in zip(locs, labels)])
                if det == 0 and syn not in ok:
                    return j - 1
    return order


def t_strict():
    for spec, order in (("cat:3", 2), ("cat:3,verify", 2), ("cat:4", 2), ("cat:4,verify", 2),
                        ("css_prep:steane", 1)):
        model = Model(library.load(spec))
        got = analysis.strict_ft(model, order)
        check(f"strict {spec}: equals the definition by brute force",
              got["holds_through"] == brute_strict(model, order),
              f"{got['holds_through']} vs {brute_strict(model, order)}")


# ---------------------------------------------------------------- 5. known answers
def t_library(heavy):
    rec, model = certify.build(library.load("cat:4"), strict_order=2)
    r = rec["results"]
    check("cat4: a single fault leaves a weight-2 error (the hook)",
          r["strict_ft"]["holds_through"] == 0 and r["strict_ft"]["violation"]["order"] == 1
          and r["fault_distance"]["lower"] == 1 and r["census"]["hook_capable"] >= 1)
    rec, _ = certify.build(library.load("cat:4,verify"), strict_order=2, heavier=(1,))
    check("cat4 with one verification: strict through order 1 at least",
          rec["results"]["strict_ft"]["holds_through"] >= 1)
    h = rec["results"]["heavier_than"][0]
    check("cat4 with one verification: an accepted error heavier than 1 needs 2 faults or more",
          h["lower"] >= 2)
    try:
        certify.verify(rec, full=True)
        check("heavier-than receipt verifies in full", True)
    except certify.Refusal as e:
        check("heavier-than receipt verifies in full", False, str(e))
    raw, _ = certify.build(library.load("cat:4"), max_weight=2, heavier=(1,))
    h = raw["results"]["heavier_than"][0]
    check("cat4 chain: one fault already leaves an accepted error heavier than 1",
          h["exact"] and h["lower"] == 1 and len(h["witness"]) == 1)
    bad = copy.deepcopy(raw)
    bad["results"]["heavier_than"][0]["witness"][0]["fault"] = "Z"
    bad["results"]["heavier_than"][0]["witness"][0]["location"] = 0
    refuses("control: heavier-than witness that is not one", lambda: certify.verify(bad), "does not replay")
    for basis in ("z", "x"):
        good, _ = certify.build(library.load(f"surface:3,2,{basis}"), max_weight=4)
        bad, _ = certify.build(library.load(f"surface:3,2,{basis},bad"), max_weight=4)
        cg, cb = good["results"]["circuit_distance"], bad["results"]["circuit_distance"]
        check(f"surface d=3 {basis}: good schedule has circuit distance 3", cg["exact"] and cg["lower"] == 3)
        check(f"surface d=3 {basis}: hook-aligned schedule has circuit distance 2",
              cb["exact"] and cb["lower"] == 2)
        for rr in (good, bad):
            try:
                certify.verify(rr, full=True)
                check(f"surface d=3 {basis}: receipt verifies in full", True)
            except certify.Refusal as e:
                check(f"surface d=3 {basis}: receipt verifies in full", False, str(e))
    rec, _ = certify.build(library.load("css_prep:steane"), strict_order=1)
    check("steane fan-out preparation is not fault tolerant",
          rec["results"]["strict_ft"]["holds_through"] == 0)
    if heavy:
        rec, _ = certify.build(library.load("surface:5,2"), max_weight=4, cap=6_000_000)
        check("surface d=5: no undetected logical error below weight 5",
              rec["results"]["circuit_distance"]["lower"] == 5)
        rec, _ = certify.build(library.load("surface:5,2,z,bad"), max_weight=4, cap=6_000_000)
        cd = rec["results"]["circuit_distance"]
        check("surface d=5 hook-aligned: circuit distance 3", cd["exact"] and cd["lower"] == 3)


# ---------------------------------------------------------------- 6. receipts and controls
def t_controls():
    circ = library.load("surface:3,2")
    rec, model = certify.build(circ, max_weight=4, strict_order=0)
    again = json.loads(json.dumps(rec))
    check("receipt survives a JSON round trip", again == rec)

    def tamper(fn):
        r = copy.deepcopy(rec)
        fn(r)
        return r
    refuses("control: witness with a fault removed",
            lambda: certify.verify(tamper(lambda r: r["results"]["circuit_distance"]["witness"].pop())),
            "witness")

    def raise_lower(r):
        r["results"]["circuit_distance"]["lower"] = 4
    refuses("control: lower bound above the witness", lambda: certify.verify(tamper(raise_lower)), "exceeds")

    def drop_ev(r):
        r["results"]["circuit_distance"]["evidence"] = []
    refuses("control: lower bound with no evidence", lambda: certify.verify(tamper(drop_ev)), "not backed")

    def edit_circuit(r):
        r["circuit"]["ticks"][1].pop()
    refuses("control: circuit edited under a stale hash", lambda: certify.verify(tamper(edit_circuit)), "")

    def swap_fault(r):
        w = r["results"]["circuit_distance"]["witness"]
        w[0]["location"] = (w[0]["location"] + 7) % r["model"]["locations"]
        w[0]["fault"] = "X"
    refuses("control: witness pointing at another location", lambda: certify.verify(tamper(swap_fault)), "")

    def census(r):
        r["results"]["census"]["hook_capable"] += 1
    refuses("control: census edited", lambda: certify.verify(tamper(census)), "census")

    def forge(r):
        cd = r["results"]["circuit_distance"]
        cd.update(lower=4, upper=None, witness=None, exact=False)
        cd["evidence"] = [dict(leg="A", status="proved", claim=">= 4")]
    refuses("control: forged enumeration claim caught by --full",
            lambda: certify.verify(tamper(forge), full=True), "does not confirm")
    rec2, _ = certify.build(library.load("cat:4"), strict_order=1)

    def fake_violation(r):
        r["results"]["strict_ft"]["violation"]["faults"] = [dict(location=0, tick=0, op="RX 0", fault="Z")]
    r2 = copy.deepcopy(rec2)
    fake_violation(r2)
    refuses("control: strict-FT counterexample that is not one", lambda: certify.verify(r2), "does not replay")
    for phase, sign in (("S", 1), ("SDG", 0)):
        phase_circuit = Circuit([[("RX", 0)], [("S", 0)], [(phase, 0)], [("MX", 0)]])
        check(f"signed ideal S followed by {phase}",
              Model(phase_circuit).sim.parity_form([0]) == sign)
    ticks = [[("R", 0), ("RX", 1)], [("M", 0)], [("Z_IF", 1, 0)], [("MX", 1)]]
    feedback = Circuit(ticks, detectors=[], observables=[[1]])
    check("flipped control reaches measured observable",
          Model(feedback).replay([(2, "flip")]) == (0, 1, 0))
    check("feedback survives native JSON roundtrip",
          Circuit.from_json(feedback.to_json()).canonical() == feedback.canonical())
    refuses("same-tick record cannot control Z",
            lambda: Circuit([ticks[0], ticks[1] + ticks[2]]), "earlier tick", CircuitError)
    refuses("negative native witness location",
            lambda: Model(feedback).replay([(-1, "flip")]), "out of range", ValueError)
    # circuit-level refusals
    for name, fn, needle in [
        ("observables without declared detectors",
         lambda: Circuit([[("R", 0)], [("M", 0)]], observables=[[0]]), "must declare"),
        ("qubit used before preparation", lambda: Circuit([[("H", 0)]]), "before preparation"),
        ("qubit used twice in a tick", lambda: Circuit([[("R", 0), ("R", 1)], [("CX", 0, 1), ("H", 0)]]), "twice"),
        ("unknown operation", lambda: Circuit([[("R", 0)], [("T", 0)]]), "unknown"),
    ]:
        refuses(f"control: {name}", fn, needle, CircuitError)
    bad = Circuit([[("RX", 0)], [("M", 0)]], detectors=[[0]])
    refuses("control: declared detector that is random", lambda: Model(bad), "not deterministic", ValueError)


# ---------------------------------------------------------------- 7. SAT leg
def t_sat():
    from ftqc_codes import sat
    if not sat.find_cadical():
        skip("SAT leg", "cadical not found (PATH or FTQC_CODES_CADICAL)")
        return
    for spec, d in (("surface:3,2", 3), ("surface:3,2,z,bad", 2), ("surface:3,2,x", 3)):
        circ = library.load(spec)
        with tempfile.TemporaryDirectory() as wd:
            rec, _ = certify.build(circ, max_weight=4, legs=("B",), workdir=wd)
            cd = rec["results"]["circuit_distance"]
            check(f"sat {spec}: circuit distance {d} with a checked proof and a replayed witness",
                  cd["exact"] and cd["lower"] == d
                  and (d == 1 or any(e["leg"] == "B" and e["status"] == "proved" for e in cd["evidence"])))
            both, _ = certify.build(circ, max_weight=4, legs=("A", "B"), workdir=wd)
            ev = both["results"]["circuit_distance"]["evidence"]
            check(f"sat {spec}: agrees with enumeration",
                  {e["claim"] for e in ev if e["status"] == "proved"} == {f">= {d}"})


# ---------------------------------------------------------------- 8. Stim leg
def t_stim(heavy):
    why = stimleg.available()
    if why:
        skip("Stim leg", why)
        return
    specs = ["cat:4", "cat:4,verify", "css_prep:steane", "css_prep:golay", "surface:3,2",
             "surface:3,3,x", "surface:3,2,z,bad"] + (["surface:5,2"] if heavy else [])
    for spec in specs:
        circ = library.load(spec)
        model = Model(circ)
        ev = stimleg.compare(circ, model)
        check(f"stim {spec}: same set of fault mechanisms", ev["status"] == "agrees", str(ev))
    rnd = random.Random(8)
    for _ in range(40):
        circ = random_circuit(rnd, rnd.randint(2, 5))
        ev = stimleg.compare(circ, Model(circ))
        check("stim random circuit: discovered detectors accepted, mechanisms equal",
              ev["status"] == "agrees", str(ev))
    # controls: a model with one location's faults removed, and a wrong detector
    circ = library.load("surface:3,2")
    model = Model(circ)
    victim = next(l for l in model.locations if l.op.name == "CX")
    keep = [s for s in model.columns
            if s not in {model.signature(*e) for e in victim.faults.values()}]
    broken = copy.copy(model)
    broken.columns = keep
    check("stim control: a fault model missing mechanisms is reported",
          len(keep) == len(model.columns) or stimleg.compare(circ, broken)["status"] == "DISAGREES")
    wrong = copy.copy(model)
    wrong.detectors = model.detectors[:-1] + [[0]]
    check("stim control: a non-deterministic detector is reported",
          stimleg.compare(circ, wrong)["status"] == "DISAGREES")
    rec, _ = certify.build(circ, max_weight=4)
    ev = stimleg.compare(circ, model, distance_lower=9)
    check("stim control: a bound above Stim's own logical error is reported",
          ev["status"] == "DISAGREES" and "below" in ev["reason"])


# ---------------------------------------------------------------- 9. CLI
def t_cli():
    env = dict(os.environ, PYTHONPATH=TOOLS)
    with tempfile.TemporaryDirectory() as wd:
        out = os.path.join(wd, "r.json")

        def run(*args):
            return subprocess.run([sys.executable, "-m", "ftqc_audit", *args], env=env, cwd=wd,
                                  capture_output=True, text=True)
        p = run("analyze", "surface:3,2", "-o", out)
        check("cli analyze exits 0", p.returncode == 0 and "circuit distance 3" in p.stdout, p.stderr[-300:])
        p = run("verify", out, "--full")
        check("cli verify --full exits 0", p.returncode == 0, p.stderr[-300:])
        rec = json.load(open(out))
        rec["results"]["circuit_distance"]["lower"] = 5
        json.dump(rec, open(out, "w"))
        check("cli verify refuses a tampered receipt with exit 1", run("verify", out).returncode == 1)
        check("cli unknown gadget exits 2", run("analyze", "nosuch:1").returncode == 2)
        cj = os.path.join(wd, "c.json")
        p = run("export", "cat:4,verify", "-o", cj)
        p2 = run("analyze", cj, "--strict", "1")
        check("cli export and re-analyze", p.returncode == 0 and p2.returncode == 0
              and "strict FT through order" in p2.stdout, p2.stderr[-300:])


def main():
    heavy = "--heavy" in sys.argv
    t0 = time.monotonic()
    for fn, args in ((t_stab, ()), (t_dense, ()), (t_columns, ()), (t_search, ()), (t_strict, ()),
                     (t_library, (heavy,)), (t_controls, ()), (t_sat, ()), (t_stim, (heavy,)),
                     (t_cli, ())):
        t = time.monotonic()
        fn(*args)
        print(f"      {fn.__name__:12s} {time.monotonic() - t:6.2f} s")
    print(f"ftqc_audit selftest: {len(PASS)} passed, {len(FAIL)} failed, {len(SKIP)} skipped "
          f"in {time.monotonic() - t0:.1f} s")
    if FAIL:
        print("FAILED: " + "; ".join(sorted(set(FAIL))[:12]))
        return 1
    print("ALL PASS" + (f" ({len(SKIP)} named skips)" if SKIP else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
