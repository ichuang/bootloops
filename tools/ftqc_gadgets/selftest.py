#!/usr/bin/env python3
"""ftqc_gadgets battery.  Stdlib only (plus ftqc_audit and ftqc_codes from this repository).

    python3 tools/ftqc_gadgets/selftest.py
    python3 tools/ftqc_gadgets/selftest.py --heavy    # adds the Golay gadgets at t = 3

Exit 0 only if every check that ran passed.
"""
import copy
import itertools
import math
import json
import os
import random
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from ftqc_faults.audit import analysis as fca  # noqa: E402
from ftqc_faults.audit import library as fclib  # noqa: E402
from ftqc_faults.audit.selftest import Dense, dense_run  # noqa: E402
from ftqc_gadgets import certify, effects, library, props  # noqa: E402
from ftqc_gadgets.effects import Analysis  # noqa: E402
from ftqc_gadgets.gadget import Builder, Code, Gadget, GadgetError  # noqa: E402
from ftqc_codes import gf2  # noqa: E402

PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    if not cond:
        print(f"FAIL  {name}  {detail}")


def refuses(name, fn, needle="", exc=(certify.Refusal, GadgetError)):
    try:
        fn()
    except exc as e:
        check(name, needle in str(e), f"refused for another reason: {e}")
        return
    check(name, False, "was accepted")


# ---------------------------------------------------------------- 1. structure
def t_structure():
    g = library.load("knill:steane,anc=spec")
    check("structure: Knill on Steane has one input block with one reference",
          g.shape() == dict(inputs=[[7, 1]], outputs=[7]) and len(g.decoders) == 2)
    an = Analysis(g)
    check("structure: references and output are in a pure state", an.pure and an.ns == 8)
    g2 = Gadget.from_json(json.loads(json.dumps(g.to_json())))
    check("structure: JSON round trip keeps the hash", g2.sha256() == g.sha256())
    check("structure: same analysis after the round trip", Analysis(g2).stats() == an.stats())
    check("structure: an ideal gadget and its implementation have the same ideal map",
          Analysis(library.load("knill:steane,ideal")).group_id == an.group_id
          and Analysis(library.load("idle:steane")).group_id == an.group_id)
    check("structure: |0_L> and |+_L> preparations have different ideal maps",
          Analysis(library.load("prep:steane,zero")).group_id
          != Analysis(library.load("prep:steane,plus")).group_id)

    def broken_preamble():
        d = g.to_json()
        t = next(t for t in d["circuit"]["ticks"] if any(op[0] == "!CX" for op in t))
        t.pop()
        d["circuit"].pop("sha256")
        d.pop("sha256")
        Analysis(Gadget.from_json(d))
    refuses("structure: a preamble with a gate removed is refused", broken_preamble, "input block")

    def noisy_preamble():
        d = g.to_json()
        t = next(t for t in d["circuit"]["ticks"] if any(op[0] == "!CX" for op in t))
        t[0][0] = "CX"
        d["circuit"].pop("sha256")
        d.pop("sha256")
        Analysis(Gadget.from_json(d))
    refuses("structure: a noisy operation in the preamble is refused", noisy_preamble)

    def random_check():
        b = Builder("bad")
        A = b.input(library.code_of("steane"))
        recs = b.measure(A, "Z")
        b.decoder("d", recs, [0b0000001], "record")
        b.output([], None)
        Analysis(b.build())
    refuses("structure: a decoder check that is random is refused", random_check, "not deterministic")

    def undeclared():
        d = g.to_json()
        d["circuit"]["detectors"] = None
        d["circuit"].pop("sha256")
        d.pop("sha256")
        Gadget.from_json(d)
    refuses("structure: decoders without declared detectors are refused", undeclared, "declare")

    def wrong_code():
        b = Builder("bad")
        A = b.input(library.code_of("steane"))
        b.tick([("CX", A[0], A[1])])
        b.output(A, library.code_of("steane"))
        Analysis(b.build()).weigher()
    refuses("structure: an output block outside its declared code is refused", wrong_code, "declared code")

    # declared detectors of a tree are exactly the discovered ones
    for spec in ("prep:steane,zero,tree=XZ", "prep:steane,plus,tree=ZX", "prep:steane,zero,tree=ZZ"):
        g = library.load(spec)
        d = g.to_json()
        d["circuit"]["detectors"] = None
        d["circuit"].pop("sha256")
        d.pop("sha256")
        a1, a2 = Analysis(g), Analysis(Gadget.from_json(d))
        r1, r2 = gf2.rank(a1.det_masks), gf2.rank(a2.det_masks)
        both = gf2.rank(a1.det_masks + a2.det_masks)
        check(f"structure: declared checks of {spec} span the discovered ones", r1 == r2 == both,
              f"{r1} {r2} {both}")


# ---------------------------------------------------------------- 2. decoder tables
def t_tables():
    for name, uniq, size in (("steane", 1, 8), ("golay", 3, 2048)):
        c = library.code_of(name)
        tab, u = effects.leader_table(c.hz, c.n)
        check(f"tables: {name} has {size} syndromes, unique through weight {uniq}",
              len(tab) == size and u == uniq, f"{len(tab)} {u}")
    rnd = random.Random(5)
    for trial in range(40):
        n = rnd.randint(3, 8)
        rows = [rnd.randrange(1, 1 << n) for _ in range(rnd.randint(1, 4))]
        tab, u = effects.leader_table(rows, n)
        best = {}
        for v in range(1 << n):
            s = sum(1 << i for i, r in enumerate(rows) if gf2.dot(r, v))
            if s not in best or v.bit_count() < best[s][0]:
                best[s] = (v.bit_count(), v)
        ok = set(tab) == set(best) and all(tab[s].bit_count() == best[s][0] for s in tab)
        ok = ok and all(sum(1 << i for i, r in enumerate(rows) if gf2.dot(r, tab[s])) == s for s in tab)
        # unique radius by brute force
        bu = 0
        for w in range(1, n + 1):
            seen, clash = {}, False
            for v in range(1 << n):
                if v.bit_count() <= w:
                    s = sum(1 << i for i, r in enumerate(rows) if gf2.dot(r, v))
                    clash = clash or s in seen
                    seen[s] = v
            if clash:
                break
            bu = w
        check("tables: coset leaders and unique radius agree with brute force", ok and u == bu,
              f"{rows} n={n} u={u} brute={bu}")


# ---------------------------------------------------------------- 3. dense semantics
def iface_vector(circ, state, keep):
    """Pure state on `keep` when every other qubit is in a basis state; None otherwise."""
    basis = {}
    for _, _, q, bas in circ.meas:
        basis[q] = bas
    for q in range(circ.n):
        if q not in keep and basis.get(q) == "X":
            state.h(q)
    rest = sum(1 << q for q in range(circ.n) if q not in keep)
    pat, vec = None, {}
    for i, a in enumerate(state.v):
        if abs(a) < 1e-9:
            continue
        if pat is None:
            pat = i & rest
        elif i & rest != pat:
            return None
        vec[sum(((i >> q) & 1) << k for k, q in enumerate(keep))] = a
    return vec


def overlap(v1, v2):
    return abs(sum(a.conjugate() * v2.get(k, 0) for k, a in v1.items()))


def tiny_bare():
    """One bare input qubit encoded into the 3-qubit code, two syndrome bits read with
    bare ancillas, a syndrome decoder for X errors."""
    b = Builder("tiny_bare")
    A = b.input(None, "in", 1)
    q = b.alloc(2)
    s = b.alloc(2)
    b.no_detectors()
    b.tick([("R", q[0]), ("R", q[1]), ("R", s[0]), ("R", s[1])])
    b.tick([("CX", A[0], q[0])])
    b.tick([("CX", A[0], q[1])])
    data = [A[0], q[0], q[1]]
    b.tick([("CX", data[0], s[0]), ("CX", data[1], s[1])])
    b.tick([("CX", data[1], s[0]), ("CX", data[2], s[1])])
    recs = b.tick([("M", s[0]), ("M", s[1])])
    b.syndrome_decoder("x", [[recs[0]], [recs[1]]], [0b011, 0b110], "X", data)
    b.output(data, None, "out")
    return b.build()


def t_dense(heavy):
    rnd = random.Random(23)
    seen = dict(detected=0, zero=0, nonzero=0, decoder_acted=0, input_error=0, frame_only=0,
                impossible=0)
    gadgets = [library.load(s) for s in ("knill:rep3,anc=raw", "steane_ec:rep3,anc=raw",
                                         "cnot:rep3", "seq_syndrome:rep3",
                                         "prep:rep3,zero,tree=X", "prep:rep3,plus,tree=Z")] + [tiny_bare()]
    for g in gadgets:
        an = Analysis(g)
        circ, m = g.circuit, an.model
        where = {i: (k, j) for i, (k, j, _) in enumerate(circ.locations())}
        in_locs = {i for b in g.inputs for i in b.in_locs}
        floc = [l for l in m.locations if l.index not in in_locs]
        mloc = [l for l in floc if l.op.name in ("M", "MX")]
        cells = [(b, p) for b, blk in enumerate(g.inputs) for p in range(len(blk.qubits))]
        exp_det = [m.sim.parity_form(d) & 1 for d in m.detectors]
        trials = 60 if circ.n <= 10 else 12
        for trial in range(trials):
            pool = mloc if (mloc and trial % 4 == 3) else floc
            chosen = rnd.sample(pool, rnd.randint(0, min(2, len(pool))))
            faults = [(l.index, rnd.choice(sorted(l.faults))) for l in chosen]
            inerr = []
            if cells and rnd.random() < 0.5:
                b, p = rnd.choice(cells)
                inerr = [(b, p, rnd.choice("XYZ"))]
            sg = an.replay(faults, inerr)
            det, dec, _ = an.parts(sg)
            eff = an.effect(sg)
            fmap = {where[i]: lab for i, lab in faults}
            for b, p, pl in inerr:
                fmap[where[g.inputs[b].in_locs[p]]] = pl
            for prob, reported, state in dense_run(circ, fmap):
                fired = 0
                for i, d in enumerate(m.detectors):
                    if (sum(reported[x] for x in d) & 1) != exp_det[i]:
                        fired |= 1 << i
                check("dense: detector pattern equals the model's", fired == det, f"{g.name} {faults}")
                if fired:
                    seen["detected"] += 1
                    continue
                corrected, frame, off, acted = list(reported), [], 0, False
                dec_ok = True
                for d, info in zip(g.decoders, an.dec_info):
                    sy = 0
                    for i, idx in enumerate(d.bits):
                        if (sum(reported[x] for x in idx) & 1) != (m.sim.parity_form(idx) & 1):
                            sy |= 1 << i
                    dec_ok = dec_ok and sy == (dec >> off) & ((1 << len(d.bits)) - 1)
                    off += len(d.bits)
                    for p in gf2.supp(info["table"][sy]):
                        acted = True
                        if d.kind == "record":
                            corrected[d.targets[p]] ^= 1
                        else:
                            frame.append((d.targets[p], d.kind))
                check("dense: decoder check values equal the model's", dec_ok, f"{g.name} {faults}")
                believed = list(dense_run(circ, {}, forced=corrected))
                if not believed:
                    seen["impossible"] += 1
                    continue
                bstate = believed[0][2]
                for q, pl in frame:
                    bstate.pauli(q, pl)
                v1 = iface_vector(circ, state, g.interface)
                v2 = iface_vector(circ, bstate, g.interface)
                if v1 is None or v2 is None:
                    check("dense: interface state is pure", False, g.name)
                    continue
                same = overlap(v1, v2) > 1 - 1e-7
                check("dense: references and outputs are in the believed state iff the effect is zero",
                      same == (eff == 0), f"{g.name} faults={faults} input={inerr} eff={eff}")
                seen["zero" if eff == 0 else "nonzero"] += 1
                seen["decoder_acted"] += acted
                seen["input_error"] += bool(inerr)
                seen["frame_only"] += eff != 0 and bool(faults) and all(l == "flip" for _, l in faults)
    seen.pop("impossible")
    check("dense: enough cases of each kind", all(v >= 20 for v in seen.values()), str(seen))


# ---------------------------------------------------------------- 4. linearity and weights
def t_linear():
    rnd = random.Random(3)
    for spec in ("knill:steane,anc=raw", "cnot:steane", "prep:steane,zero,tree=XZ", "steane_ec:steane"):
        an = Analysis(library.load(spec))
        m = an.model
        in_locs = {i for b in an.g.inputs for i in b.in_locs}
        floc = [l for l in m.locations if l.index not in in_locs]
        cols = set(an.cols)
        ok = True
        for _ in range(60):
            chosen = rnd.sample(floc, rnd.randint(1, 3))
            faults = [(l.index, rnd.choice(sorted(l.faults))) for l in chosen]
            tot = 0
            for i, lab in faults:
                s = an.sig(*m.locations[i].faults[lab])
                ok = ok and (s == 0 or s in cols)
                tot ^= s
            ok = ok and an.replay(faults) == tot
        check(f"linear: replay of a set is the XOR of its members ({spec})", ok)
        ok = True
        for l in floc:
            sigs = {an.sig(*e) for e in l.faults.values()} | {0}
            ok = ok and all((a ^ b) in sigs for a in sigs for b in sigs)
        check(f"linear: each location's effects are closed under addition ({spec})", ok)
    for spec in ("cnot:steane", "knill:steane,anc=spec"):
        an = Analysis(library.load(spec))
        for sector in (False, True):
            ins = an.input_errors(2, sector)
            ok = all(an.replay([], lab) == sg and props._weight(lab, sector) == w
                     for sg, (w, lab) in ins.items())
            check(f"linear: every input error replays to its signature and weight ({spec}, "
                  f"{'sector' if sector else 'pauli'})", ok and len(ins) > 20)
        pairs = props._input_pairs(an, an, 2)
        check(f"linear: input errors of a refinement replay to their signatures ({spec})",
              all(an.replay([], [tuple(x) for x in lab]) == s2 == s1 for (s2, s1), lab in pairs.items()))
    # sums() against brute force, with and without detectors
    for spec, j in (("prep:steane,zero,tree=X", 2), ("cnot:rep3", 3), ("knill:rep3", 2)):
        an = Analysis(library.load(spec))
        targets = {0} | {c & an.DET for c in an.cols[:5]}
        ok = True
        for d in targets:
            fast = sorted(v for v, _ in an.sums(j, d))
            slow = []
            for idx in itertools.combinations(range(len(an.cols)), j):
                v = 0
                for i in idx:
                    v ^= an.cols[i]
                if v & an.DET == d:
                    slow.append(v)
            ok = ok and fast == sorted(slow)
        check(f"linear: accepted subsets of size {j} agree with brute force ({spec})", ok)


def t_weights():
    rnd = random.Random(9)
    for spec, sector in (("cnot:steane", False), ("cnot:steane", True), ("knill:steane,anc=spec", False),
                         ("prep:steane,zero", False), ("fc:cat:4", False)):
        an = Analysis(library.load(spec))
        W = an.weigher(sector)
        base = sum(len(b.refs) for b in an.g.inputs)
        blocks, off = [], base
        for b in an.g.outputs:
            blocks.append(list(range(off, off + len(b.qubits))))
            off += len(b.qubits)
        # brute force: least per-block weight over Paulis of weight <= 1 per block
        best = {}
        per_block = []
        for ip in blocks:
            opts = [(0, 0)]
            for p in ip:
                opts += [(1, an.sx[p]), (1, an.sz[p]), (1, an.sx[p] ^ an.sz[p])]
            if sector:
                opts += [(1, an.sx[p] ^ an.sz[q]) for p in ip for q in ip if p != q]
            per_block.append(opts)
        for combo in itertools.product(*per_block):
            e, w = 0, 0
            for wt, s in combo:
                e ^= s
                w = max(w, wt)
            best[e] = min(best.get(e, 9), w)
        ok = all(W.weight(e, 1) == w for e, w in best.items())
        for _ in range(200):
            e = rnd.getrandbits(an.ns)
            ok = ok and (W.weight(e, 1) is None) == (e not in best)
        check(f"weights: least weight agrees with brute force ({spec}, {'sector' if sector else 'pauli'})", ok)
    an = Analysis(library.load("cnot:steane"))
    W = an.weigher()
    e1 = an.sx[an.mi - 1]
    check("weights: a bound below an already built table is enforced",
          W.weight(e1, 1) == 1 and W.weight(e1, 0) is None)
    refuses("weights: a weight beyond the unique radius is refused",
            lambda: an.weigher().weight(1, 2), "unique", exc=effects.Budget)


# ---------------------------------------------------------------- 5. properties against the definition
def brute_check(an, kind, t, sector=False):
    """The definition, by replay over every location subset, every label, every input Pauli."""
    g, m = an.g, an.model
    W = an.weigher(sector)
    in_locs = {i for b in g.inputs for i in b.in_locs}
    floc = [l for l in m.locations if l.index not in in_locs]
    cells = [(b, p) for b, blk in enumerate(g.inputs) for p in range(len(blk.qubits))]
    holds = -1
    for total in range(t + 1):
        good = True
        for r in range((0 if kind == "prep" else total) + 1):
            s = total - r
            for cs in itertools.combinations(cells, r):
                for paulis in itertools.product("XYZ", repeat=r):
                    inerr = [(b, p, pl) for (b, p), pl in zip(cs, paulis)]
                    for ls in itertools.combinations(floc, s):
                        for labs in itertools.product(*[sorted(l.faults) for l in ls]):
                            sg = an.replay([(l.index, lab) for l, lab in zip(ls, labs)], inerr)
                            if sg & an.DET:
                                continue
                            bound = r + s if kind == "gate" else s      # stated here independently
                            if W.weight(an.effect(sg), bound) is None:
                                good = False
        if not good:
            break
        holds = total
    return holds


def t_props(heavy):
    cases = [("cnot:steane", "gate", 1), ("cnot:steane,leaky", "gate", 1), ("idle:steane", "gate", 1),
             ("knill:steane,anc=spec", "ec", 1), ("knill:steane,anc=raw", "ec", 1),
             ("steane_ec:steane,anc=spec", "ec", 1), ("seq_syndrome:steane", "ec", 1),
             ("prep:steane,zero", "prep", 1), ("prep:steane,zero,spec", "prep", 1),
             ("prep:steane,zero,tree=X", "prep", 1)]
    for spec, kind, t in cases:
        an = Analysis(library.load(spec))
        r = props.check(an, kind, t)
        b = brute_check(an, kind, t)
        check(f"props: {kind} of {spec} agrees with the definition", r["holds_through"] == b,
              f"{r['holds_through']} vs {b}")
        if r["violation"]:
            check(f"props: the witness for {spec} replays as a violation",
                  props.judge(an, kind, r["violation"]))
    known = [("cnot:steane", "gate", 1, 1), ("cnot:steane,leaky", "gate", 1, 0),
             ("knill:steane,anc=spec", "ec", 1, 1), ("knill:steane,anc=raw", "ec", 1, 0),
             ("knill:steane,anc=tree=XZ", "ec", 1, 1), ("knill:steane,anc=tree=X", "ec", 1, 0),
             ("steane_ec:steane,anc=spec", "ec", 1, 1), ("steane_ec:steane,anc=tree=XZ", "ec", 1, 1),
             ("seq_syndrome:steane", "ec", 1, 0),
             ("prep:steane,zero", "prep", 1, 0), ("prep:steane,zero,tree=XZ", "prep", 1, 1),
             ("prep:steane,plus,tree=ZX", "prep", 1, 1), ("cnot:golay", "gate", 2, 2),
             ("knill:golay,anc=spec", "ec", 2, 2), ("knill:golay,anc=raw", "ec", 1, 0)]
    if heavy:
        known += [("cnot:golay", "gate", 3, 3), ("knill:golay,anc=spec", "ec", 3, 3)]
    for spec, kind, t, want in known:
        r = props.check(Analysis(library.load(spec)), kind, t)
        check(f"known: {kind} of {spec} holds through {want}", r["holds_through"] == want,
              f"got {r['holds_through']} {r['note']} {r['violation']}")
    r = props.check(Analysis(library.load("cnot:steane")), "gate", 1, sector=True)
    check("known: transversal CX on Steane, sector weights", r["holds_through"] == 1, str(r["note"]))
    # the error-correction pair
    r = props.close(Analysis(library.load("knill:steane,anc=spec")), 1)
    check("close: Knill with good ancillas, any input, through 1",
          r["holds_through"] == 1 and "does not depend" in r["method"], str(r))
    r = props.close(Analysis(library.load("knill:steane,anc=raw")), 1)
    check("close: Knill with raw ancillas fails at one fault", r["holds_through"] == 0 and r["violation"])
    r = props.close(Analysis(library.load("steane_ec:steane,anc=spec")), 1)
    check("close: Steane-type extraction, shown independent of the input", r["holds_through"] == 1
          and "does not depend" in r["method"], str(r))
    r2 = props.close(Analysis(library.load("seq_syndrome:steane")), 1)
    check("close: sequential extraction corrects X only: decided by enumerating the input span",
          r2["holds_through"] == -1 and "enumerated" in (r2["method"] or ""), str(r2)[:300])
    r = props.close(Analysis(library.load("idle:steane")), 1)
    check("close: a bare wire is not close for any input", r["holds_through"] == -1 and r["violation"])
    refuses("props: a gate property on a preparation is refused",
            lambda: props.check(Analysis(library.load("prep:steane,zero")), "gate", 1), "input")
    # the sequential-extraction counterexample: X on data qubit 2 after the first check
    g = library.load("seq_syndrome:steane")
    an = Analysis(g)
    W = an.weigher()
    code = library.code_of("steane")
    first = gf2.supp(code.hz[0])
    locs = g.circuit.locations()
    last_of_first = [i for i, (_, _, op) in enumerate(locs) if op.name == "CX"][len(first) - 1]
    target = next(p for p in range(7) if not (code.hz[0] >> p) & 1 or p == first[-1])
    bad = [f for f in props.check(an, "ec", 1)["violation"]["faults"]]
    check("known: sequential extraction has a single-fault violation with a replayable witness",
          props.judge(an, "ec", dict(r=0, s=1, input_error=[], faults=bad)))


# ---------------------------------------------------------------- 6. refinement
def brute_effects(an, order):
    """{effect: least number of fault locations} by replay over location subsets and labels."""
    m = an.model
    in_locs = {i for b in an.g.inputs for i in b.in_locs}
    floc = [l for l in m.locations if l.index not in in_locs]
    out = {}
    for j in range(order + 1):
        for ls in itertools.combinations(floc, j):
            for labs in itertools.product(*[sorted(l.faults) for l in ls]):
                sg = an.replay([(l.index, lab) for l, lab in zip(ls, labs)])
                if not sg & an.DET:
                    out.setdefault(an.effect(sg), j)
    return out


def t_refine(heavy):
    pairs = [("fc:cat:4", "fc:cat:4,verify", 2), ("fc:cat:4,verify", "fc:cat:4", 2),
             ("prep:steane,zero,tree=X", "prep:steane,zero,spec", 1),
             ("prep:steane,zero", "prep:steane,zero,spec", 1),
             ("prep:steane,zero,spec", "prep:steane,zero", 1)]
    for s2, s1, order in pairs:
        a2, a1 = Analysis(library.load(s2)), Analysis(library.load(s1))
        r = props.refine(a2, a1, order)
        e2, e1 = brute_effects(a2, order), brute_effects(a1, order)
        want = -1
        for j in range(order + 1):
            if any(c == j and e1.get(e, 99) > j for e, c in e2.items()):
                break
            want = j
        check(f"refine: {s2} against {s1} agrees with brute force", r["holds_through"] == want,
              f"{r['holds_through']} vs {want}")
        if r["violation"]:
            check(f"refine: the unmatched fault set of {s2} replays", props.judge_refine(a2, a1, r["violation"]))
    # strict fault tolerance (ftqc_audit) is refinement of the idealized gadget
    for spec in ("cat:3", "cat:4", "cat:4,verify", "css_prep:steane", "cat:5,verify"):
        circ = fclib.load(spec)
        g = library.from_circuit(circ)
        an = Analysis(g)
        r = props.refine(an, Analysis(library.idealize(g)), 2)
        s = fca.strict_ft(an.model, 2)
        check(f"refine: against the idealized gadget equals ftqc_audit strict_ft ({spec})",
              max(r["holds_through"], 0) == s["holds_through"], f"{r['holds_through']} vs {s['holds_through']}")
    # with inputs
    a_c, a_i = Analysis(library.load("cnot:steane")), Analysis(library.load("cnot:steane,ideal"))
    check("refine: transversal CX and its specification are fault-equivalent, input weight 1",
          props.refine(a_c, a_i, 2, 1)["holds_through"] == 2
          and props.refine(a_i, a_c, 2, 1)["holds_through"] == 2)
    a_l = Analysis(library.load("cnot:steane,leaky"))
    r = props.refine(a_l, a_i, 1, 0)
    check("refine: the leaky CX does not, with a replayable witness",
          r["holds_through"] == 0 and props.judge_refine(a_l, a_i, r["violation"]))
    a_k = Analysis(library.load("knill:steane,anc=tree=XZ"))
    a_ki = Analysis(library.load("knill:steane,ideal"))
    check("refine: Knill with purified ancillas refines ideal correction followed by noisy wires",
          props.refine(a_k, a_ki, 1, 1)["holds_through"] == 1)
    check("refine: a bare wire does not refine ideal correction once the input has an error",
          props.refine(Analysis(library.load("idle:steane")), a_ki, 1, 1)["holds_through"] == 0)
    for s in ("knill:steane,anc=spec", "cnot:steane", "prep:steane,zero,tree=X"):
        a = Analysis(library.load(s))
        b = Analysis(library.load(s))
        check(f"refine: {s} is fault-equivalent to itself",
              props.refine(a, b, 2)["holds_through"] == 2)
    refuses("refine: gadgets with different ideal maps are not comparable",
            lambda: props.refine(Analysis(library.load("prep:steane,zero")),
                                 Analysis(library.load("prep:steane,plus")), 1), "ideal map")
    refuses("refine: gadgets with different blocks are not comparable",
            lambda: props.refine(a_c, a_ki, 1), "different input or output")
    # a specification that is too strong: no wires at all
    g = library.load("prep:steane,zero,tree=XZ")
    b = Builder("no_wires")
    b.no_detectors()
    b.output(library.frag_prep(b, library.code_of("steane"), "zero", "ideal"), library.code_of("steane"))
    r = props.refine(Analysis(g), Analysis(b.build()), 1)
    check("refine: nothing with faults refines a noiseless specification", r["holds_through"] == 0)
    if heavy:
        a2 = Analysis(library.load("prep:golay,zero,tree=XZ"))
        a1 = Analysis(library.load("prep:golay,zero,spec"))
        r = props.refine(a2, a1, 2)
        check("refine: Golay XZ tree refines its specification through order 1 only",
              r["holds_through"] == 1, str(r["holds_through"]))


# ---------------------------------------------------------------- 7. composition and receipts
def t_compose():
    code = "steane"
    target = certify.build(library.load(f"knill:{code},anc=spec"), [("ec", 1)], close=1)
    parts = {}
    for role, seq in (("plus", "ZX"), ("zero", "XZ")):
        parts[role] = certify.build(library.load(f"prep:{code},{role},tree={seq}"),
                                    refines=[(library.prep(code, role, "spec"), 1, 0)])
    comp = certify.compose(target, parts)
    check("compose: Knill with trees inherits both properties through order 1",
          [c["holds_through"] for c in comp["conclusions"]] == [1, 1])
    certify.verify_composition(comp, target, parts)
    direct = props.check(Analysis(library.load(f"knill:{code},anc=tree=XZ")), "ec", 1)
    check("compose: the direct analysis of the composed gadget agrees", direct["holds_through"] >= 1)
    weak = dict(parts)
    weak["plus"] = certify.build(library.load(f"prep:{code},plus,tree=X"),
                                 refines=[(library.prep(code, "plus", "spec"), 1, 0)])
    comp2 = certify.compose(target, weak)
    direct2 = props.check(Analysis(library.load(f"knill:{code},anc=tree=X")), "ec", 1)
    check("compose: a part that refines only through order 0 lowers the conclusion to 0",
          [c["holds_through"] for c in comp2["conclusions"]] == [0, 0] and direct2["holds_through"] == 0)
    refuses("compose: a missing part is refused", lambda: certify.compose(target, {"zero": parts["zero"]}),
            "no refinement receipt")
    swapped = dict(plus=parts["zero"], zero=parts["zero"])
    refuses("compose: a part refined against another specification is refused",
            lambda: certify.compose(target, swapped), "different specification")
    forged = copy.deepcopy(comp)
    forged["conclusions"][0]["holds_through"] = 3
    refuses("compose: an inflated conclusion is refused",
            lambda: certify.verify_composition(forged, target, parts), "conclusions")
    # two levels: a tree whose leaves are specifications, composed with real leaves, then
    # used as the implementation of a part of the error-correction gadget
    sc = library.code_of(code)
    b = Builder("two_spec_leaves_then_X")
    b.no_detectors()
    leaves = b.par(*[(lambda: library.frag_prep(b, sc, "zero", "spec")) for _ in range(2)])
    b.output(library._purify(b, sc, "zero", "X", leaves), sc)
    top = certify.build(b.build(), refines=[(library.prep(code, "zero", "spec"), 1, 0)])
    leaf = certify.build(library.load(f"prep:{code},zero,tree=X"),
                         refines=[(library.prep(code, "zero", "spec"), 1, 0)])
    inner = certify.compose(top, {"zero": leaf})
    check("compose: a refinement is carried through a composition",
          [c["type"] for c in inner["conclusions"]] == ["refine"]
          and inner["conclusions"][0]["holds_through"] == min(top["results"][0]["holds_through"],
                                                              leaf["results"][0]["holds_through"]))
    bundle = dict(composition=inner, target=top, parts={"zero": leaf})
    two = certify.compose(target, dict(plus=parts["plus"], zero=bundle))
    direct3 = props.refine(Analysis(library.load(f"prep:{code},zero,tree=XX")),
                           Analysis(library.prep(code, "zero", "spec")), 1)
    check("compose: a composed implementation can itself be a part; the direct analysis agrees",
          two["conclusions"][0]["holds_through"] == inner["conclusions"][0]["holds_through"]
          <= direct3["holds_through"], f"{two['conclusions']} {direct3['holds_through']}")
    bad_bundle = dict(composition=copy.deepcopy(inner), target=top, parts={"zero": leaf})
    bad_bundle["composition"]["conclusions"][0]["holds_through"] = 3
    refuses("compose: an inflated inner composition is refused",
            lambda: certify.compose(target, dict(plus=parts["plus"], zero=bad_bundle)), "conclusions")
    refuses("compose: a gadget analysed without specification parts is refused",
            lambda: certify.compose(certify.build(library.load("knill:steane,anc=raw"), [("ec", 1)]), parts),
            "no part")


def t_receipts():
    g = library.load("knill:steane,anc=raw")
    rec = certify.build(g, [("ec", 1)], close=1,
                        refines=[(library.load("knill:steane,ideal"), 1, 0)])
    done = certify.verify(rec, full=True)
    check("receipt: builds and verifies in full", len(done) >= 6, str(done))
    rec2 = json.loads(json.dumps(rec))
    check("receipt: survives a JSON round trip", certify.verify(rec2) and rec2 == rec)

    def tamper(fn, full=False):
        r = copy.deepcopy(rec)
        fn(r)
        return lambda: certify.verify(r, full=full)

    def drop_fault(r):
        r["results"][0]["violation"]["faults"] = []
        r["results"][0]["violation"]["s"] = 0
    refuses("receipt: a witness with the fault removed is refused", tamper(drop_fault), "not a violation")

    def harmless(r):
        v = r["results"][0]["violation"]
        loc = v["faults"][0][0]
        v["faults"][0][0] = next(i for i in range(len(g.circuit.locations()))
                                 if i != loc and not _is_violation(g, i, v["faults"][0][1]))
    refuses("receipt: a witness moved to a harmless location is refused", tamper(harmless), "")

    def inflate(r):
        r["results"][0]["holds_through"] = 1
    refuses("receipt: holds_through above the witness is refused", tamper(inflate), "contradicts")

    def no_evidence(r):
        r["results"][0]["violation"] = None
        r["results"][0]["holds_through"] = 1
    refuses("receipt: a positive claim not covered by cases is refused", tamper(no_evidence), "not covered")

    def edit_circuit(r):
        r["gadget"]["circuit"]["ticks"][-1].pop()
    refuses("receipt: an edited circuit under a stale hash is refused", tamper(edit_circuit), "sha256")

    def edit_stats(r):
        r["stats"]["columns"] += 1
    refuses("receipt: edited statistics are refused", tamper(edit_stats), "statistics")

    def forge_cases(r):
        r["results"][0]["per"][0]["cases"] += 1
    refuses("receipt: a forged case count is refused in full verification", tamper(forge_cases, True),
            "re-enumeration")

    def refine_witness(r):
        res = r["results"][2]
        res["violation"]["faults"] = []
    refuses("receipt: a refinement witness that the specification matches is refused",
            tamper(refine_witness), "matched")

    def swap_spec(r):
        r["results"][2]["spec"] = library.load("idle:steane").to_json()
        r["results"][2]["violation"] = None
        r["results"][2]["holds_through"] = 1
    refuses("receipt: a refinement claim against another specification is refused in full",
            tamper(swap_spec, True), "")

    def decoder_edit(r):
        r["gadget"]["decoders"][0]["checks"].pop()
        r["gadget"]["decoders"][0]["bits"].pop()
    refuses("receipt: an edited decoder under a stale hash is refused", tamper(decoder_edit), "sha256")
    refuses("receipt: another schema is refused", lambda: certify.verify(dict(schema="x")), "not a")


def _is_violation(g, loc, label):
    an = Analysis(g)
    try:
        return props.judge(an, "ec", dict(r=0, s=1, input_error=[], faults=[[loc, label]]))
    except ValueError:
        return False


# ---------------------------------------------------------------- 9. locality leg
def _coords(g):
    """qubit -> coordinate for library gadgets: blocks are allocated n qubits at a time."""
    n = len((g.inputs + g.outputs)[0].qubits)
    refs = {q for b in g.inputs for q in b.refs}
    qs = [q for q in range(g.circuit.n) if q not in refs]
    return {q: i % n for i, q in enumerate(qs)}, n


def t_local(heavy):
    rnd = random.Random(17)
    cases = [("cnot:steane", "gate", None, True, False), ("idle:steane", "gate", None, True, False),
             ("knill:steane,anc=spec", "ec", 1, True, True), ("steane_ec:steane,anc=spec", "ec", 1, True, True),
             ("steane_ec:steane,anc=spec,rounds=2", "ec", 1, True, True),
             ("measure:steane", "ec", 1, True, True), ("measure:steane,x", "ec", 1, True, True),
             ("cnot:golay", "gate", None, True, False), ("knill:golay,anc=spec", "ec", 3, True, True),
             ("steane_ec:golay,anc=spec", "ec", 3, True, True), ("measure:golay", "ec", 3, True, True)]
    for spec, kind, through, gate, ec in cases:
        g = library.load(spec)
        an = Analysis(g)
        r = props.local(an)
        ok = r["local"] and r["through"] == through and r["gate"] == gate and r["ec"] == ec and r["sector"]
        check(f"local: {spec} is coordinate-local with the expected radius and properties", ok, str(r)[:300])
        if not r["local"]:
            continue
        t_code = 1 if "steane" in spec else 3
        T = min(through if through is not None else t_code, 1 if "steane" in spec else (3 if heavy else 2))
        e = props.check(an, kind, T)
        check(f"local: enumeration agrees for {spec} through {T}", e["holds_through"] == T, str(e["note"]))
        if g.outputs and ec:
            c = props.close(an, T)
            check(f"local: closeness by enumeration agrees for {spec}", r["close"] and c["holds_through"] == T)
        # the stronger statement: any number of faults inside few coordinates
        coord, n = _coords(g)
        m = an.model
        in_locs = {i for b in g.inputs for i in b.in_locs}
        by_c = {}
        for l in m.locations:
            if l.index not in in_locs:
                cs = {coord[q] for q in l.op.qs}
                if len(cs) == 1:
                    by_c.setdefault(cs.pop(), []).append(l)
        W = an.weigher()
        limit = through if through is not None else t_code
        good = tried = 0
        for _ in range(150):
            k = rnd.randint(1, limit)
            cs = rnd.sample(range(n), k)
            pool = [l for c in cs for l in by_c.get(c, [])]
            chosen = rnd.sample(pool, min(len(pool), rnd.randint(1, 3 * k + 2)))
            faults = [(l.index, rnd.choice(sorted(l.faults))) for l in chosen]
            fc = {coord[m.locations[i].op.qs[0]] for i, _ in faults}
            inerr = [(b, c, rnd.choice("XYZ")) for b in range(len(g.inputs)) for c in cs if rnd.random() < 0.4]
            eff = an.effect(an.replay(faults, inerr))
            bound = len(fc) if kind == "ec" else k
            tried += 1
            good += eff != -1 and W.weight(eff, min(bound, t_code)) is not None if bound <= t_code else eff != -1
        check(f"local: many faults inside at most {limit} coordinates stay within the bound ({spec})",
              good == tried, f"{good}/{tried}")
    for spec, why in (("cnot:steane,leaky", "confined"), ("seq_syndrome:steane", ""),
                      ("knill:steane,anc=raw", ""), ("prep:steane,zero,tree=XZ", "postselects")):
        r = props.local(Analysis(library.load(spec)))
        check(f"local: {spec} is refused", r["local"] is False and why in r["reason"], str(r)[:200])
    # a correction applied one position over: the decoder is triggered at p, the error stays
    sc = library.code_of("steane")
    b = Builder("shifted_correction")
    A = b.input(sc, "data")
    b.no_detectors()
    P = library.frag_prep(b, sc, "plus", "spec")
    b.tick([("CX", a, q) for a, q in zip(A, P)])
    recs = b.measure(P, "Z")
    b.decoder("x_errors", recs, sc.hz, "X", A[1:] + A[:1])
    b.output(A, sc, "data")
    r = props.local(Analysis(b.build()))
    check("local: a correction applied at the neighbouring position is refused",
          r["local"] is False and "not confined to that coordinate" in r["reason"], str(r)[:200])
    b = Builder("measurement with temporary swap")
    data = b.input(sc)
    b.no_detectors()
    for _ in range(2):
        b.tick([("I", data[0])])
        for i, j in ((0, 1), (1, 0), (0, 1)):
            b.tick([("!CX", data[i], data[j])])
    records = b.measure(data, "Z")
    b.decoder("outcomes", records, sc.hz, "record")
    swapped = Analysis(b.build())
    faults = [(7, "X"), (8, "X")]
    check("two faults on physical wire zero flip the decoded Z", all(
        swapped.model.locations[i].op.qs == (data[0],) for i, _ in faults)
        and swapped.effect(swapped.replay(faults)) != 0)
    check("ideal temporary swaps invalidate the physical column certificate",
          not props.local(swapped)["local"])
    # decoder positions numbered differently from the block's
    perm = [3, 0, 1, 2, 6, 4, 5]
    b = Builder("permuted_decoder")
    A = b.input(sc, "data")
    b.no_detectors()
    recs = b.measure(A, "Z")
    rows = [sum(1 << i for i in range(7) if (h >> perm[i]) & 1) for h in sc.hz]
    b.decoder("outcomes", [recs[perm[i]] for i in range(7)], rows, "record")
    g = b.build()
    check("local: the permuted decoder is a correct decoder", props.check(Analysis(g), "ec", 1)["holds_through"] == 1)
    r = props.local(Analysis(g))
    check("local: but its positions are not the block's coordinates, and that is refused",
          r["local"] is False and "two coordinates" in r["reason"] and r.get("input_error"), str(r)[:200])
    ag = Analysis(library.load("steane_ec:golay,anc=spec"))
    r = props.local(ag, dict(x_errors=3, z_errors=1))
    check("local: the smallest decoder radius is the one that counts", r["local"] and r["through"] == 1)
    an = Analysis(library.load("knill:steane,anc=spec"))
    check("local: a radius above the one the table shows is refused", props.local(an, 2)["local"] is False)
    # two blocks, three coordinates of dirt, a code that corrects one: beyond the radius the
    # statement must not be claimed
    r = props.local(an)
    check("local: the certificate never claims more than the decoders' radius", r["through"] == 1)
    # radii from a code receipt; analysis without tables
    import ftqc_codes
    for name, want in (("steane", 1), ("golay", 3)):
        crec = ftqc_codes.build(ftqc_codes.load(name))
        g = library.load(f"knill:{name},anc=spec")
        rr = library.radii(g, crec)
        check(f"local: decoder radii from the {name} code receipt", set(rr.values()) == {want}, str(rr))
        a2 = Analysis(g, tables=False)
        r2 = props.local(a2, rr)
        check(f"local: without decoder tables the certificate is the same ({name})",
              r2["local"] and r2["through"] == want and r2["ec"] and r2["close"]
              and all(d["source"] == "given" for d in r2["decoders"]))
        e = props.check(a2, "ec", 1)
        check(f"local: enumeration without tables stops with a note ({name})",
              e["holds_through"] < 1 and "lookup table" in (e["note"] or ""), str(e)[:200])
        check("local: without a radius and without tables nothing is claimed", props.local(a2)["local"] is False)
    meas = library.measure("steane", "Z")
    rejected = certify.build(meas, local=7, tables=False)
    certify.verify(rejected)
    check("radius seven without tables is disproved for Steane measurement",
          rejected["results"][0]["local"] is False and "proved radius 1" in rejected["results"][0]["reason"])
    pair = Analysis(meas)
    check("two record flips contradict the proposed seven-column tolerance",
          pair.effect(pair.replay([(11, "flip"), (12, "flip")])) != 0)
    valid = certify.build(meas, local=1, tables=False)
    certify.verify(valid)
    decoder = valid["results"][0]["decoders"][0]
    check("the no-table proof establishes the kernel lower bound",
          decoder["proof"] == dict(method="minimum weight of check kernel", kernel_lower=3, empty_kernel=False))
    legacy = copy.deepcopy(valid)
    del legacy["results"][0]["radius_request"]
    del legacy["results"][0]["decoders"][0]["proof"]
    certify.verify(legacy)
    legacy["results"][0]["through"] = legacy["results"][0]["decoders"][0]["radius"] = 7
    refuses("a legacy unproved seven-column receipt cannot verify", lambda: certify.verify(legacy), "recompute")
    forged = copy.deepcopy(valid)
    forged["results"][0]["decoders"][0]["proof"]["kernel_lower"] = 15
    refuses("an inflated kernel proof is recomputed", lambda: certify.verify(forged), "recompute")
    # out of budget, the only other proof is the BCH bound, which needs a cyclic kernel:
    # the Steane code in another qubit order has none, in its usual order it has
    st_ = library.code_of("steane")
    perm = [1, 0, 2, 3, 4, 5, 6]
    moved = Code("steane_moved", 7, [sum(1 << perm[i] for i in gf2.supp(r)) for r in st_.hx],
                 [sum(1 << perm[i] for i in gf2.supp(r)) for r in st_.hz])
    check("a kernel search that runs out of budget cannot certify the proposed radius",
          not props.local(Analysis(library.measure(moved, "Z"), tables=False, cap=0), 1)["local"])
    cyc = props.local(Analysis(meas, tables=False, cap=0), 1)
    check("out of budget on a cyclic kernel, the radius rests on the BCH bound and the proof says so",
          cyc["local"] and "BCH" in cyc["decoders"][0]["proof"].get("bound_from", "")
          and not props.local(Analysis(meas, tables=False, cap=0), 2)["local"], str(cyc)[:300])
    # repeated correction: later rounds do not depend on the input
    a = Analysis(library.load("steane_ec:steane,anc=spec,rounds=2"))
    ind = props.independent_of_input(a)
    check("rounds: the second round's decoders cannot be changed by any input error, the first round's can",
          ind["x_errors_2"]["independent"] and ind["z_errors_2"]["independent"]
          and not ind["x_errors_1"]["independent"] and not ind["z_errors_1"]["independent"], str(ind))
    for b, p, t in [e for v in ind.values() if v["input_error"] for e in v["input_error"]]:
        sg = a.replay([], [(b, p, t)])
        check("rounds: the witness input error does change a first-round check", a.parts(sg)[1] != 0)
    a1 = Analysis(library.load("steane_ec:steane,anc=spec"))
    check("rounds: with one round every decoder depends on the input",
          not any(v["independent"] for v in props.independent_of_input(a1).values()))
    # receipts
    rec = certify.build(library.load("steane_ec:golay,anc=spec"), local=True)
    certify.verify(rec)
    check("local: receipt verifies and its summary states the result",
          "coordinate-local on 23 coordinates" in certify.summary(rec) and "at most 3" in certify.summary(rec))
    bad = copy.deepcopy(rec)
    bad["results"][0]["through"] = 5
    refuses("local: an inflated radius in a receipt is refused", lambda: certify.verify(bad), "recompute")
    bad = copy.deepcopy(rec)
    bad["results"][0]["decoders"][0].update(radius=5, source="given")
    refuses("local: a given radius above the table's is refused", lambda: certify.verify(bad), "recompute")
    rec = certify.build(library.load("cnot:steane,leaky"), local=True)
    certify.verify(rec)
    check("local: a negative certificate carries its fault", rec["results"][0]["local"] is False
          and rec["results"][0].get("fault"))


# ---------------------------------------------------------------- 10. two levels
def t_levels(heavy):
    from ftqc_gadgets import levels
    jobs = [("measure", dict(basis="Z"), 2, 3000), ("measure", dict(basis="X"), 1, 1500),
            ("cnot", {}, 1, 1500), ("steane_ec", dict(anc="spec"), 2 if heavy else 1, 2500)]
    if heavy:
        jobs.append(("knill", dict(anc="spec"), 1, 4000))
    for fam, kw, exh, samples in jobs:
        g = levels.two_level(fam, "steane", "steane", **kw)
        d = levels.Direct(g, "steane", "steane")
        r = d.run(exhaustive=exh, samples=samples, search=40000)
        tag = f"{fam}{''.join(',' + str(v) for v in kw.values())}"
        check(f"two levels, direct ({tag}): with at most one bad column every examined set acts ideally",
              r["failures"] == 0 and r["exhaustive_sets"] > 100 and r["sampled_sets"] > samples // 2
              and (r["t1"], r["t2"]) == (1, 1), str({k: v for k, v in r.items() if k != "beyond"}))
        check(f"two levels, direct ({tag}): two bad columns can fail, and a failing set is found",
              r["beyond"] is not None and r["beyond"]["bad_columns"] == 2, str(r["beyond"]))
        # the outer certificate says the same thing: one bad coordinate, not two
        outer = props.local(Analysis(library.load(
            {"measure": "measure:steane," + kw.get("basis", "Z").lower(), "cnot": "cnot:steane",
             "steane_ec": "steane_ec:steane,anc=spec", "knill": "knill:steane,anc=spec"}[fam])))
        check(f"two levels ({tag}): the outer locality certificate tolerates exactly what the direct test shows",
              outer["local"] and outer["through"] in (None, 1))
    g = levels.two_level("measure", "steane", "steane", basis="Z")
    g2 = Gadget.from_json(json.loads(json.dumps(g.to_json())))
    check("two levels: a gadget with hierarchical decoders survives a JSON round trip",
          g2.sha256() == g.sha256() and g2.decoders[0].rule["type"] == "hier")
    check("two levels: the hierarchical rule is part of the gadget's hash",
          library.measure(ftqc_codes_concat(), "Z").sha256() != g.sha256())
    lite = Code("steane_lite", 7, list(library.code_of("steane").hx), list(library.code_of("steane").hz))
    check("two levels: a plain code record is accepted for either level",
          levels.two_level("measure", lite, lite, basis="Z").decoders[0].rule["sector"] == "X")
    refuses("two levels: a gadget that joins two positions is refused",
            lambda: levels.positions(library.load("cnot:steane,leaky")), "joins two positions", exc=GadgetError)
    # reduction receipts
    outer = certify.build(library.load("steane_ec:steane,anc=spec"), local=True)
    parts = {role: certify.build(library.load(f"prep:steane,{role},tree={seq}"),
                                 refines=[(library.prep("steane", role, "spec"), 1, 0)])
             for role, seq in (("plus", "ZX"), ("zero", "XZ"))}
    ec_spec = certify.build(library.load("steane_ec:steane,anc=spec"), [("ec", 1)], close=1)
    ec = dict(composition=certify.compose(ec_spec, parts), target=ec_spec, parts=parts)
    inner = dict(cnot=certify.build(library.load("cnot:steane"), [("gate", 1)]), ec=ec,
                 measure_z=certify.build(library.load("measure:steane,z"), [("ec", 1)]),
                 measure_x=certify.build(library.load("measure:steane,x"), [("ec", 1)]))
    lv = levels.reduction(outer, inner)
    check("levels: Steane correction over Steane blocks: one bad column, one fault per inner rectangle",
          lv["conclusion"] == dict(bad_columns_tolerated=1, inner_faults_per_rectangle=1)
          and [r["role"] for r in lv["inner"]] == ["cnot", "ec", "measure_x", "measure_z"]
          and len(lv["outer"]["assumed"]) == 2, str(lv["conclusion"]))
    levels.verify_reduction(lv, outer, inner)
    check("levels: summary names the assumptions", "assumed: plus ancilla" in levels.summary(lv))
    refuses("levels: a missing inner receipt is refused",
            lambda: levels.reduction(outer, {k: v for k, v in inner.items() if k != "cnot"}), "no inner receipt")
    refuses("levels: an inner correction analysed with specified ancillas is refused",
            lambda: levels.reduction(outer, dict(inner, ec=ec_spec)), "specified parts")
    refuses("levels: an inner correction without the any-input property is refused",
            lambda: levels.reduction(outer, dict(inner, ec=certify.build(library.load("idle:steane"),
                                                                           [("gate", 1)]))), "does not establish")
    ec_only = certify.build(library.load("steane_ec:steane,anc=spec"), [("ec", 1)])
    refuses("levels: an inner correction with ec but without close is refused",
            lambda: levels.reduction(outer, dict(inner, ec=dict(
                composition=certify.compose(ec_only, parts), target=ec_only, parts=parts))),
            "does not establish")
    refuses("levels: inner receipts about different codes are refused",
            lambda: levels.reduction(outer, dict(inner, cnot=certify.build(library.load("cnot:golay"),
                                                                             [("gate", 1)]))), "same code")
    refuses("levels: an outer receipt without a locality certificate is refused",
            lambda: levels.reduction(ec_spec, inner), "locality")
    refuses("levels: an outer gadget that is not local is refused",
            lambda: levels.reduction(certify.build(library.load("cnot:steane,leaky"), local=True), inner),
            "locality")
    idle = certify.build(library.idle("steane"), [("gate", 1)])
    refuses("one-block idle cannot fill a CNOT role",
            lambda: levels.reduction(outer, dict(inner, cnot=idle)), "arity")
    no_radius = certify.build(library.measure("steane", "Z"), local=7, tables=False)
    refuses("an unsupported outer radius cannot enter a level reduction",
            lambda: levels.reduction(no_radius, inner), "no positive locality")
    swapped = dict(inner, measure_x=inner["measure_z"], measure_z=inner["measure_x"])
    refuses("X and Z measurement roles cannot be exchanged",
            lambda: levels.reduction(outer, swapped), "logical relation")
    def single_gate(name):
        b = Builder("single-qubit gates")
        data = b.input(library.code_of("steane"))
        b.tick([(name, q) for q in data])
        b.output(data, library.code_of("steane"))
        return b.build()
    h = single_gate("H")
    oh = certify.build(h, local=True)
    refuses("transversal H needs an inner H", lambda: levels.reduction(oh, inner), "no inner receipt for h")
    refuses("identity and H have different signed Choi relations",
            lambda: levels.reduction(oh, dict(inner, h=idle)), "logical relation")
    lv_h = levels.reduction(oh, dict(inner, h=certify.build(h, [("gate", 1)])))
    check("a verified logical Hadamard fills the H role", lv_h["inner"][-1]["logical_action"]["operation"] == "h")
    refuses("even an ideal unsupported outer phase gate needs an implementation",
            lambda: levels.reduction(certify.build(single_gate("!S"), local=True), inner), "unsupported outer operation S")
    # Transversal S squared gives a logical Z with the identity's unsigned Choi group.
    b = Builder("signed identity")
    data = b.input(library.code_of("steane"))
    for gate in ("S", "S"):
        b.tick([(gate, q) for q in data])
    b.output(data, library.code_of("steane"))
    oi = certify.build(library.idle("steane"), local=True)
    refuses("a logical Pauli sign cannot be hidden in an idle receipt",
            lambda: levels.reduction(oi, dict(inner, idle=certify.build(b.build(), [("gate", 1)]))), "wrong sign")
    bad = copy.deepcopy(lv)
    bad["conclusion"]["bad_columns_tolerated"] = 2
    refuses("levels: an inflated conclusion is refused", lambda: levels.verify_reduction(bad, outer, inner),
            "conclusion")
    weak = dict(inner, cnot=certify.build(library.load("cnot:steane"), [("gate", 0)]))
    check("levels: a weaker inner receipt lowers the inner radius",
          levels.reduction(outer, weak)["conclusion"]["inner_faults_per_rectangle"] == 0)

def _tree_variant(kind, code_name="steane"):
    """Two raw copies and one X check, with one thing wrong."""
    from ftqc_gadgets.gadget import state_checks
    code = library.code_of(code_name)
    n = code.n
    b = Builder("variant_" + kind)
    b.no_detectors()
    A, B = b.par(lambda: library.frag_prep(b, code, "zero", "raw"),
                 lambda: library.frag_prep(b, code, "plus" if kind == "states" else "zero", "raw"))
    if kind == "inside":
        b.tick([("CX", x, y) for x, y in zip(A, B)])
        b.tick([("CX", A[0], A[1])])
        b.tick([("CX", A[0], A[1])])
    elif kind == "skew":
        b.tick([("CX", x, y) for x, y in zip(A, B[1:] + B[:1])])
    elif kind == "partial":
        b.tick([("CX", x, y) for x, y in zip(A[:-1], B[:-1])])
    else:
        b.tick([("CX", y, x) if kind == "reversed" else ("CX", x, y) for x, y in zip(A, B)])
    if kind == "waits":
        b.tick([("I", y) for y in B])
    recs = b.tick([("M", y) for y in B] + [("I", x) for x in A])
    masks = state_checks(code, "zero", "Z")
    if kind == "fewer":
        masks = masks[:-1]
    if kind == "foreign":
        masks = masks[:-1] + [masks[-1] ^ 1]
    b.checks(recs[:n], masks)
    b.output(A, code)
    return b.build()


def _ec_with_trees(code_name, seq_plus, seq_zero, rounds=1):
    """Steane error correction with purified ancillas, and for every fault location:
    (coordinate or None, quiet?, leaf id or None).  Quiet: an operation of the correction
    itself, a check operation of a tree, or the last idle of a leaf."""
    code = library.code_of(code_name)
    n = code.n
    b = Builder(f"ec_with_trees({code_name})")
    A = b.input(code, "data")
    b.no_detectors()
    mark = {}

    def tree(state, seq):
        before = {id(op) for t in b.ticks for op in t}
        qs = library.frag_prep(b, code, state, "tree=" + seq, rounds)
        for t in b.ticks:
            for op in t:
                if id(op) not in before:
                    mark[id(op)] = state
        return qs
    P, Z = b.par(lambda: tree("plus", seq_plus), lambda: tree("zero", seq_zero))
    b.tick([("CX", a, p) for a, p in zip(A, P)])
    r1 = b.tick([("M", p) for p in P] + [("CX", z, a) for z, a in zip(Z, A)])[:n]
    r2 = b.tick([("MX", z) for z in Z])
    b.decoder("x_errors", r1, code.hz, "X", A)
    b.decoder("z_errors", r2, code.hx, "Z", A)
    b.output(A, code, "data")
    ticks = [list(t) for t in b.ticks]
    g = b.build()
    refs = {q for blk in g.inputs for q in blk.refs}
    wires = [q for q in range(g.circuit.n) if q not in refs]
    pos = {q: i % n for i, q in enumerate(wires)}
    blk = {q: i // n for i, q in enumerate(wires)}
    last_leaf_idle = {}
    info, i = [], 0
    flat = [(k, op) for k, t in enumerate(ticks) for op in t if not op[0].startswith("!")]
    first_cross = {}
    for k, op in flat:
        if op[0] == "CX" and blk[op[1]] != blk[op[2]]:
            for q in op[1:]:
                first_cross.setdefault(blk[q], k)
    for k, op in flat:
        qs = op[1:]
        state = mark.get(id(op))
        cols = {pos[q] for q in qs}
        col = cols.pop() if len(cols) == 1 else None
        if state is None:
            info.append((col, True, None))
        else:
            in_leaf = k < first_cross[blk[qs[0]]]
            quiet = (not in_leaf) or (op[0] == "I" and k == first_cross[blk[qs[0]]] - 1)
            info.append((col, quiet, (state, blk[qs[0]]) if in_leaf else None))
    return g, info


def t_factories(heavy):
    from ftqc_gadgets import factories, levels
    INF = factories.INF
    # ---- the contract against enumeration
    for spec, seq, need in (("prep:steane,zero,raw", "", (1, 1)), ("prep:steane,zero,tree=X", "X", (2, 1)),
                            ("prep:steane,zero,tree=Z", "Z", (1, 2)), ("prep:steane,plus,tree=Z", "Z", (1, 2)),
                            ("prep:steane,zero,tree=XZ", "XZ", (2, 2)), ("prep:steane,plus,tree=ZX", "ZX", (2, 2)),
                            ("prep:steane,zero,tree=XX", "XX", (4, 1)), ("prep:steane,zero,tree=ZXZ", "ZXZ", (2, 4)),
                            ("prep:steane,zero,tree=XZ,rounds=1", "XZ", (2, 2))):
        g = library.load(spec)
        c = factories.certificate(g)
        k = c.get("junk_leaves_needed", {})
        check(f"factory {spec}: recognised, sequence and number of junk leaves needed",
              c["tree"] and c["sequence"] == seq and (k["X"], k["Z"]) == need
              and c["through"] == min(need) - 1 and c["leaves"] == 2 ** len(seq), str(c)[:300])
        tr = factories.Tree(g)
        an = Analysis(g)
        r = factories.test(g, tr, exhaustive=2 if g.circuit.n <= 21 else 1, samples=9000, analysis=an)
        check(f"factory {spec}: the rule agrees with propagation on every examined set",
              r["failures"] == 0 and r["accepted"] > 800 and r["largest_accepted"] >= 8, str(r))
        # the count K against the rule itself: fewest leaves with a junk-capable fault that raise the flag
        junk = {}
        for b_, leaf in tr.leaves.items():
            cls = tr.shapes[leaf["shape"]]["classes"]
            junk[b_] = {t: next((leaf["locs"][i] for i, row in enumerate(cls) if row[t][0] == "junk"), None)
                        for t in "XZ"}
        for t, want in zip("XZ", need):
            best = INF
            blocks = sorted(junk)
            for w in range(1, len(blocks) + 1):
                if best < INF:
                    break
                for sub in itertools.combinations(blocks, w):
                    if all(junk[b_][t] is not None for b_ in sub) and tr.predict([junk[b_][t] for b_ in sub])[t][0]:
                        best = w
                        break
            check(f"factory {spec}: junk leaves needed ({t}) is the least number that raises the flag", best == want,
                  f"{best} vs {want}")
        # the profile a probability calculus reads, and the exact law of the accepted output
        pf = tr.profile()
        kind0 = c["leaf_kinds"][0]
        check(f"factory {spec}: the profile counts every location of the tree once",
              pf["coordinates"] == tr.n and len(pf["leaves"]) == c["leaves"] and len(pf["nodes"]) == c["checks"]
              and all(pf["shapes"][0][t]["junk"] == kind0[f"junk_capable_{t}"]
                      and max(pf["shapes"][0][t]["local"]) == kind0[f"max_local_per_coordinate_{t}"] for t in "XZ")
              and all(nd["gates"] == [1] * tr.n and nd["meas"] == [1] * tr.n for nd in pf["nodes"])
              and sum(pf["shapes"][k]["locations"] for k in pf["leaves"].values())
              + sum(sum(nd[x]) for nd in pf["nodes"] for x in ("gates", "meas", "idles")) == len(an.model.locations)
              and pf["wmin"] == c["wmin"] and tuple(pf["root"]) == tr.root, str(pf)[:300])
        if spec in ("prep:steane,zero,tree=X", "prep:steane,zero,tree=Z", "prep:steane,zero,tree=XZ"):
            p_ = 1e-5
            ex = factories.accepted_law(g, p_, tr, an)
            acc_all = 1.0
            for a_ in ex["accept"]:
                acc_all *= a_
            L_ = len(an.model.locations)
            flt = [(loc.index, an.sig(*eff), 1.0 / len(loc.faults)) for loc in an.model.locations
                   for eff in loc.faults.values()]
            series = {0: (1 - p_) ** L_}
            for _, sg, w_ in flt:
                if not sg & an.DET:
                    series[sg >> an.shift] = series.get(sg >> an.shift, 0.0) + p_ * w_ * (1 - p_) ** (L_ - 1)
            for i_, (l1, s1, w1) in enumerate(flt):
                for l2, s2, w2 in flt[i_ + 1:]:
                    if l1 != l2 and not (s1 ^ s2) & an.DET:
                        e_ = (s1 ^ s2) >> an.shift
                        series[e_] = series.get(e_, 0.0) + p_ * p_ * w1 * w2 * (1 - p_) ** (L_ - 2)
            third = math.comb(L_, 3) * p_ ** 3
            second = sum(v for e_, v in series.items() if e_) - sum(
                p_ * w_ * (1 - p_) ** (L_ - 1) for _, sg, w_ in flt if not sg & an.DET and sg >> an.shift)
            worst = max(abs(ex["law"].get(e_, 0.0) * acc_all - v) for e_, v in series.items())
            check(f"factory {spec}: the exact law of the accepted output agrees with the sum over one and two faults",
                  worst <= third and second > 20 * third and abs(sum(ex["law"].values()) - 1) < 1e-12
                  and set(ex["law"]) >= set(series), f"worst {worst:.3e}, third order {third:.3e}, pairs {second:.3e}")
            cw = factories.class_weights(ex["single"]["X"], tr.n)
            brute = {}
            for pat in range(1 << tr.n):
                e_ = 0
                for q_ in range(tr.n):
                    if (pat >> q_) & 1:
                        e_ ^= ex["single"]["X"][q_]
                brute[e_] = min(brute.get(e_, tr.n), bin(pat).count("1"))
            check(f"factory {spec}: class weights are least numbers of coordinates",
                  cw == brute and len(cw) == 2 ** len(tr.GZ) and max(cw.values()) > 1, str(sorted(cw.values())))
        # the contract implies refinement of the specification through K - 1 (on the Steane
        # code enumeration goes further: every Z pattern on |0_L> is equivalent to one error)
        if len(seq) <= 2 and "rounds" not in spec:
            state = spec.split(",")[1]
            rr = props.refine(an, Analysis(library.prep("steane", state, "spec")), min(need) - 1, 0)
            check(f"factory {spec}: enumeration against the specification confirms the contract's order",
                  rr["holds_through"] == min(need) - 1, str(rr["holds_through"]))
    # on the Golay code the order is exact: enumeration stops where the contract does
    for seq, K in (("X", 1), ("XZ", 2)) if heavy else (("X", 1),):
        g = library.load(f"prep:golay,zero,tree={seq}")
        c = factories.certificate(g)
        rr = props.refine(Analysis(g), Analysis(library.prep("golay", "zero", "spec")), K, 0)
        check(f"factory prep:golay,zero,tree={seq}: enumeration against the specification stops where the "
              f"contract does", c["through"] == K - 1 == rr["holds_through"] and rr["violation"] is not None,
              f"{c['through']} {rr['holds_through']}")
    # ---- classes against the definition, by replay on the raw circuit
    for name in ("steane", "golay") if heavy else ("steane",):
        for state in ("zero", "plus"):
            g = library.prep(name, state, "raw")
            tr = factories.Tree(g)
            an = Analysis(g)
            n = tr.n
            low = (1 << n) - 1
            mask = dict(X=sum(1 << i for i, v in enumerate(an.canon) if not v & low),
                        Z=sum(1 << i for i, v in enumerate(an.canon) if not v >> n))
            ok = True
            for loc in an.model.locations:
                for t, per in (("X", an.sx), ("Z", an.sz)):
                    hits = set()
                    for label in loc.faults:
                        e = an.replay([(loc.index, label)]) >> an.shift & mask[t]
                        if e:
                            cs = [c for c in range(n) if per[c] == e]
                            hits.add(cs[0] if cs else -1)
                    want = "none" if not hits else ("local" if len(hits) == 1 and -1 not in hits else "junk")
                    got = tr.info[loc.index][2][t]
                    ok = ok and got[0] == want and (want != "local" or got[1] == next(iter(hits)))
            check(f"factory classes ({name}, {state}): every location classed as the definition says", ok)
    c = factories.certificate(library.prep("golay", "zero", "raw"))["leaf_kinds"][0]
    check("factory classes: the library's raw Golay preparation (100 locations)",
          (c["locations"], c["junk_capable_X"], c["junk_capable_Z"]) == (100, 66, 65)
          and c["max_local_per_coordinate_X"] <= 3 and c["max_local_per_coordinate_Z"] <= 3, str(c))
    c = factories.certificate(library.prep("steane", "zero", "spec"))
    check("factory: a specification has no junk-capable location and is its own contract",
          c["tree"] and c["junk_leaves_needed"] == dict(X=None, Z=None) and c["through"] is None and c["physical"])
    # ---- strict classes
    g = library.load("prep:steane,zero,tree=XZ,rounds=2")
    cs, cn = factories.certificate(g, strict=True), factories.certificate(g)
    check("factory, strict: only the last idle round of a leaf is local, each idle in its own column",
          cs["tree"] and cs["physical"] and not cn["physical"]
          and cs["leaf_kinds"][0]["junk_capable_X"] == cs["leaf_kinds"][0]["locations"] - 7
          and cs["leaf_kinds"][0]["max_local_per_coordinate_Z"] == 1 and cs["junk_leaves_needed"] == dict(X=2, Z=2))
    r = factories.test(g, factories.Tree(g, strict=True), exhaustive=1, samples=6000)
    check("factory, strict: the rule with demoted classes still agrees with propagation", r["failures"] == 0, str(r))
    c = factories.certificate(library.load("prep:steane,zero,tree=XZ"), strict=True)
    check("factory, strict: a leaf without a final idle round is refused", not c["tree"] and "idle" in c["reason"])
    # ---- what is not a tree
    for kind, needle in (("fewer", "do not span"), ("foreign", "do not span"), ("reversed", "does not complete"),
                         ("inside", "after its leaf"), ("skew", "not a transversal"), ("partial", "every coordinate"),
                         ("waits", "idles between"), ("states", "same state")):
        c = factories.certificate(_tree_variant(kind))
        check(f"factory: a tree with {kind} is refused", not c["tree"] and needle in c["reason"], str(c))
    c = factories.certificate(_tree_variant("ok"))
    check("factory: the unmodified variant is accepted", c["tree"] and c["sequence"] == "X")
    c = factories.certificate(library.load("steane_ec:steane,anc=spec"))
    check("factory: a gadget with inputs is refused", not c["tree"])
    c = factories.certificate(library.load("prep:steane,zero,tree=X"), wmin=5)
    check("factory: a minimum weight that cannot be proved is refused", not c["tree"] and "only" in c["reason"])
    c = factories.certificate(library.load("prep:steane,zero,tree=X"), wmin=2)
    check("factory: a smaller requested minimum weight is used as given", c["tree"] and c["wmin"] == dict(X=2, Z=2))
    from ftqc_codes import bch as _bch
    b63 = _bch.cyclic_css(63, [1, 3, 5])
    c = factories.certificate(library.prep(b63, "zero", "tree=X"), wmin=7, cap=200_000)
    check("factory: on a cyclic code a minimum weight the search cannot reach is taken from the BCH bound",
          c["tree"] and c["wmin"] == dict(X=7, Z=7)
          and any("BCH" in v.get("bound_from", "") for v in c["wmin_proof"].values()), str(c)[:300])
    # ---- probabilities: the recursion against sampling with independent location faults
    from fractions import Fraction
    for spec, p, shots in (("prep:steane,zero,tree=XZ", Fraction(1, 200), 60000),
                           ("prep:steane,plus,tree=ZXZ", Fraction(1, 400), 40000),
                           ("prep:steane,zero,tree=XZ,rounds=1", Fraction(1, 300), 40000)):
        g = library.load(spec)
        tr = factories.Tree(g)
        an = Analysis(g)
        bd = tr.bounds(*tr.iid_rates(p))
        n = tr.n
        low = (1 << n) - 1
        mask = dict(X=sum(1 << i for i, v in enumerate(an.canon) if not v & low),
                    Z=sum(1 << i for i, v in enumerate(an.canon) if not v >> n))
        table = [[an.sig(*eff) for eff in loc.faults.values()] for loc in an.model.locations]
        rnd = random.Random(11)
        pf = float(p)
        acc = 0
        hit = dict(junkX=0, junkZ=0, dirtyX=0, dirtyZ=0)
        for _ in range(shots):
            sg, locs = 0, []
            for i, labs in enumerate(table):
                if rnd.random() < pf:
                    sg ^= rnd.choice(labs)
                    locs.append(i)
            if sg & an.DET:
                continue
            acc += 1
            if locs:
                pred = tr.predict(locs)
                e = sg >> an.shift
                for t in "XZ":
                    hit["junk" + t] += pred[t][0]
                    hit["dirty" + t] += bool(e & mask[t])
        want_acc = 1.0
        for a_ in bd["accept_lower"]:
            want_acc *= float(a_)
        slack = lambda k, bound: hit[k] / acc <= float(bound) + 3 * (float(bound) / acc) ** 0.5 + 3 / acc
        check(f"factory probabilities {spec}: acceptance is at least the product of the bounds",
              acc / shots >= want_acc - 0.01 and acc > shots // 4, f"{acc / shots:.4f} vs {want_acc:.4f}")
        check(f"factory probabilities {spec}: sampled junk and non-trivial rates do not exceed the bounds",
              all(slack(k + t, bd[k][t]) for k in ("junk", "dirty") for t in "XZ")
              and hit["dirtyX"] + hit["dirtyZ"] > 50,
              str({k: round(v / acc, 5) for k, v in hit.items()}) + " vs "
              + str({k: {t: round(float(v), 5) for t, v in bd[k].items()} for k in ("junk", "dirty")}))
        check(f"factory probabilities {spec}: expected leaves at least the number of leaves",
              bd["expected_leaves"] >= len(tr.leaves))
    tr = factories.Tree(library.load("prep:steane,zero,tree=XZXZ"))
    bd = tr.bounds(dict(X=Fraction(1, 100), Z=Fraction(1, 100)), dict(X=Fraction(0), Z=Fraction(0)),
                   Fraction(0), Fraction(0))
    check("factory probabilities: with no dirt the junk bound is the polynomial of the rule (8 x^4 and 64 z^4)",
          bd["junk"] == dict(X=8 * Fraction(1, 100) ** 4, Z=64 * Fraction(1, 100) ** 4)
          and bd["expected_leaves"] == 16, str(bd["junk"]))
    # one X check, every input non-zero, against the formulas written out
    tr1 = factories.Tree(library.load("prep:steane,zero,tree=X"))
    F = Fraction
    jx, jz, ux, uz, sa, sb = F(1, 7), F(1, 11), F(1, 5), F(1, 9), F(1, 13), F(1, 17)
    b1 = tr1.bounds(dict(X=jx, Z=jz), dict(X=ux, Z=uz), sa, sb)
    al = (1 - ux) ** 2 * (1 - sb)
    check("factory probabilities: one check, all four formulas and the retry count",
          b1["accept_lower"] == [al] and b1["junk"] == dict(X=jx * jx / al, Z=2 * jz / al)
          and b1["dirty"] == dict(X=(sa + ux * (ux + sb)) / al, Z=(2 * uz + sa) / al)
          and b1["expected_leaves"] == 2 / al, str(b1))
    rj, rd, ra, rb = tr1.iid_rates(F(1, 100))
    k1 = factories.certificate(library.load("prep:steane,zero,tree=X"))["leaf_kinds"][0]
    cls = tr1.shapes[0]["classes"]
    check("factory probabilities: rates for independent location faults count the right locations",
          rj == {t: 1 - F(99, 100) ** k1[f"junk_capable_{t}"] for t in "XZ"}
          and rd == {t: 1 - F(99, 100) ** sum(1 for r in cls if r[t][0] != "none") for t in "XZ"}
          and ra == rb == 1 - F(99, 100) ** 14 and rj["X"] != rj["Z"] and rj["X"] < rd["X"])
    check("factory probabilities: no bound when a check cannot be shown to accept",
          tr.bounds(dict(X=Fraction(1, 2), Z=Fraction(1, 2)), dict(X=Fraction(1), Z=Fraction(1)),
                    Fraction(0), Fraction(0)) is None)
    # ---- receipts
    g = library.load("prep:steane,zero,tree=XZ,rounds=2")
    rec = certify.build(g, tree=dict(strict=True))
    check("factory receipt: verifies", any("tree" in line for line in certify.verify(rec)))
    check("factory receipt: summary states the contract", "unless this many leaves" in certify.summary(rec))
    for name, edit in (("inflated count", lambda r: r["results"][0]["junk_leaves_needed"].update(X=4)),
                       ("strict flag removed", lambda r: r["results"][0].update(strict=False)),
                       ("physical claimed", lambda r: r["results"][0].update(physical=False)),
                       ("lemma edited", lambda r: r["results"][0].update(lemma="trust me"))):
        bad = copy.deepcopy(rec)
        edit(bad)
        refuses(f"factory receipt: {name} is refused", lambda bad=bad: certify.verify(bad), "does not recompute")
    light = certify.build(g, tree=dict(strict=True, light=True))
    check("factory receipt, light: same certificate without the full fault model, and it verifies",
          light["results"] == rec["results"] and light["stats"].get("light")
          and any("tree" in line for line in certify.verify(light)))
    bad = copy.deepcopy(light)
    bad["results"].append(dict(type="check", kind="prep", t=1, weight="pauli", holds_through=1, per=[],
                               violation=None, note=None))
    refuses("factory receipt, light: another kind of result cannot ride on it",
            lambda: certify.verify(bad), "factory certificate only")
    bad = copy.deepcopy(light)
    bad["results"][0]["leaves"] = 8
    refuses("factory receipt, light: an edited certificate is refused", lambda: certify.verify(bad),
            "does not recompute")
    refuses("factory receipt, light: enumeration cannot be requested with it",
            lambda: certify.build(g, [("prep", 0)], tree=dict(light=True)), "light", exc=GadgetError)
    bad = copy.deepcopy(certify.build(_tree_variant("ok"), tree=True))
    bad["gadget"]["circuit"] = _tree_variant("fewer").to_json()["circuit"]
    refuses("factory receipt: a circuit with a check removed under the old certificate is refused",
            lambda: certify.verify(bad))
    # ---- one level: correction with real trees, any number of faults inside few coordinates
    for name, T, K, trials in (("steane", 1, 2, 3000),) + ((("golay", 3, 2, 8000),) if heavy else ()):
        g, info = _ec_with_trees(name, "ZX", "XZ")
        an = Analysis(g)
        W = an.weigher()
        n = len(g.outputs[0].qubits)
        in_locs = {i for b_ in g.inputs for i in b_.in_locs}
        pool, noisy = {}, []
        k = 0
        for loc in an.model.locations:
            if loc.index in in_locs:
                p = [b_.in_locs.index(loc.index) for b_ in g.inputs if loc.index in b_.in_locs][0]
                for label, eff in loc.faults.items():
                    pool.setdefault(p, []).append((loc.index, label, an.sig(*eff), True))
                k += 1
                continue
            col, quiet, leaf = info[k]
            k += 1
            for label, eff in loc.faults.items():
                item = (loc.index, label, an.sig(*eff), False)
                if quiet:
                    pool.setdefault(col, []).append(item)
                else:
                    noisy.append((leaf, item))
        check(f"correction with trees ({name}): every location has a class", k == len(info))
        rnd = random.Random(5)
        seen = bad = big = 0
        by_det = {}
        for c_ in pool:
            for it in pool[c_]:
                by_det.setdefault((c_, it[2] & an.DET), []).append(it)
        for trial in range(trials):
            cols = rnd.sample(range(n), rnd.randint(1, T))
            items = []
            for c_ in cols:
                items += rnd.sample(pool[c_], rnd.randint(1, min(12, len(pool[c_]))))
            if trial % 3 == 0:                      # fewer than K leaves with a junk-capable fault
                items.append(rnd.choice(noisy)[1])
            if len({it[0] for it in items}) != len(items):
                continue
            sg = 0
            for it in items:
                sg ^= it[2]
            if sg & an.DET:                          # try to make it accepted inside the same columns
                fix = [it for c_ in cols for it in by_det.get((c_, sg & an.DET), [])
                       if it[0] not in {x[0] for x in items}]
                if not fix:
                    continue
                items.append(rnd.choice(fix))
                sg ^= items[-1][2]
            seen += 1
            big = max(big, len(items))
            if W.weight(an.effect(sg), len(cols)) is None:
                bad += 1
        check(f"correction with trees ({name}): faults confined to at most {T} coordinates, however many, "
              f"with fewer than {K} junk leaves, leave the output within that weight",
              bad == 0 and seen > trials // 20 and big > T + 3, f"{bad} of {seen}, largest {big}")
        # and K junk leaves can break it: the same junk fault in two sibling leaves
        leaves = sorted({leaf for leaf, _ in noisy})
        broke = False
        for a_, b_ in zip(leaves[::2], leaves[1::2]):
            la = [it for leaf, it in noisy if leaf == a_]
            lb = [it for leaf, it in noisy if leaf == b_]
            for x, y in zip(la, lb):
                if x[1] == y[1] and not (x[2] ^ y[2]) & an.DET and W.weight(an.effect(x[2] ^ y[2]), T) is None:
                    broke = True
                    break
            if broke:
                break
        check(f"correction with trees ({name}): two leaves with the same junk-capable fault defeat it", broke)
    # ---- the levels receipt without assumed ancillas
    outer = certify.build(library.load("steane_ec:steane,anc=spec"), local=True)
    parts = {role: certify.build(library.load(f"prep:steane,{role},tree={seq}"),
                                 refines=[(library.prep("steane", role, "spec"), 1, 0)])
             for role, seq in (("plus", "ZX"), ("zero", "XZ"))}
    ec_spec = certify.build(library.load("steane_ec:steane,anc=spec"), [("ec", 1)], close=1)
    ec = dict(composition=certify.compose(ec_spec, parts), target=ec_spec, parts=parts)
    inner = dict(cnot=certify.build(library.load("cnot:steane"), [("gate", 1)]), ec=ec,
                 measure_z=certify.build(library.load("measure:steane,z"), [("ec", 1)]),
                 measure_x=certify.build(library.load("measure:steane,x"), [("ec", 1)]),
                 prep_zero=parts["zero"], prep_plus=parts["plus"],
                 idle=certify.build(library.idle("steane"), [("gate", 1)]))
    fac = {role: certify.build(library.load(f"prep:steane,{role},tree={seq},rounds=2"), tree=dict(strict=True))
           for role, seq in (("plus", "ZX"), ("zero", "XZ"))}
    lv = levels.reduction(outer, inner, factories=fac)
    check("levels with factories: no ancilla is assumed, and the junk event is recorded",
          lv["outer"]["assumed"] == [] and [f["role"] for f in lv["factories"]] == ["plus", "zero"]
          and lv["factory_junk_leaves_needed"] == 2 and "idle" in [r["role"] for r in lv["inner"]]
          and lv["conclusion"] == dict(bad_columns_tolerated=1, inner_faults_per_rectangle=1), str(lv)[:400])
    levels.verify_reduction(lv, outer, inner, factories=fac)
    check("levels with factories: summary", "plus ancilla: tree ZX" in levels.summary(lv))
    lv1 = levels.reduction(outer, inner, factories=dict(zero=fac["zero"]))
    check("levels with factories: an ancilla without a factory stays assumed",
          lv1["outer"]["assumed"] == ["plus ancilla as its specification"])
    refuses("levels with factories: a receipt verified without its factories is refused",
            lambda: levels.verify_reduction(lv, outer, inner), "differs")
    loose = certify.build(library.load("prep:steane,zero,tree=XZ,rounds=2"), tree=True)
    refuses("levels with factories: a factory without strict classes is refused",
            lambda: levels.reduction(outer, inner, factories=dict(fac, zero=loose)), "strict")
    refuses("levels with factories: a tree for the other state is refused",
            lambda: levels.reduction(outer, inner, factories=dict(plus=fac["zero"])), "does not prepare")
    refuses("levels with factories: the factory's operations need inner receipts too",
            lambda: levels.reduction(outer, {k: v for k, v in inner.items() if k != "idle"}, factories=fac),
            "no inner receipt for idle")
    refuses("levels with factories: a factory for an ancilla the gadget does not have is refused",
            lambda: levels.reduction(certify.build(library.load("cnot:steane"), local=True), inner,
                                     factories=dict(zero=fac["zero"])), "no such specified ancilla")
    other = certify.build(library.load("prep:golay,zero,tree=XZ,rounds=1"), tree=dict(strict=True))
    refuses("levels with factories: a factory on another code is refused",
            lambda: levels.reduction(outer, inner, factories=dict(zero=other)), "outer code")
    bad = copy.deepcopy(lv)
    bad["factory_junk_leaves_needed"] = 4
    refuses("levels with factories: an inflated junk count is refused",
            lambda: levels.verify_reduction(bad, outer, inner, factories=fac), "factory_junk")

def _random_program(rnd, nblocks=4, depth=10):
    """A valid single-level program: gates, controlled gates, corrections, measurements."""
    from ftqc_gadgets import cones
    live = [f"B{i}" for i in range(nblocks)]
    inputs, layers, bits, fresh, pend = list(live), [], [], [0], {}
    for _ in range(depth):
        layer, free = [], [b for b in live if b not in pend]
        for b in list(pend):
            layer.append(("ECb", b))
            del pend[b]
        rnd.shuffle(free)
        while len(free) >= 2 and rnd.random() < 0.7:
            a, b = free.pop(), free.pop()
            k = rnd.choice(["CX", "CX", "CZ", "CXg", "CZg"])
            if k.endswith("g") and not bits:
                k = k[:2]
            layer.append((k, a, b, rnd.choice(bits)) if k.endswith("g") else (k, a, b))
        for b in free[:]:
            r = rnd.random()
            if r < 0.3:
                layer.append(("ECa", b))
                pend[b] = 1
                free.remove(b)
            elif r < 0.4 and len(live) > 2:
                fresh[0] += 1
                bit = f"m{fresh[0]}"
                layer.append((rnd.choice(["MZ", "MX"]), b, bit))
                live.remove(b)
                free.remove(b)
                new = f"N{fresh[0]}"
                layer.append(("NEW", rnd.choice(["zero", "plus", "std"]), new))
                live.append(new)
                bits.append(bit)                # usable from the next layer on
                bits_new = bit
        layers.append(layer)
        # a record becomes usable one layer after its measurement
    if pend:
        layers.append([("ECb", b) for b in pend])
    # controls must refer to earlier layers only: drop controlled gates whose bit is too new
    seen, clean = set(), []
    for layer in layers:
        row = []
        for op in layer:
            if op[0] in ("CXg", "CZg") and op[3] not in seen:
                op = (op[0][:2], op[1], op[2])
            row.append(op)
        clean.append(row)
        seen.update(op[2] for op in layer if op[0] in ("MZ", "MX"))
    return cones.Program(clean, inputs, live, "random")


def t_cones(heavy):
    from ftqc_gadgets import cones
    steane = library.code_of("steane")

    def what(flow, ev):
        return sorted((s["kind"], repr(s["what"]), s["tick"], s.get("sector")) for s in
                      (flow.sources[i] for i in cones._bits(flow.cone(ev))))
    # ---- known answers
    P = cones.Program([[("CX", "A", "B")], [("MZ", "B")]], ["A", "B"], ["A"], "copy")
    f = cones.Flow(P)
    ev = {(e["kind"], str(e["block"]), e["sector"]): e for e in f.events}
    check("cones: a measurement after a CX sees the gate and its own location, nothing else",
          what(f, ev[("meas", "B", "X")]) == [("class", "('gate', 'CX')", 0, None), ("class", "('meas', 'X')", 1, None)])
    check("cones: the control keeps the gate and its idle in both sectors",
          f.weight(f.cone(ev[("out", "A", "X")])) == (2, {}) and f.weight(f.cone(ev[("out", "A", "Z")])) == (2, {}))
    P = cones.Program([[("CX", "A", "B")], [("ECa", "B")], [("ECb", "B")]], ["A", "B"], ["A", "B"], "sink")
    f = cones.Flow(P)
    ev = {(e["kind"], str(e["block"]), e["sector"]): e for e in f.events}
    after = f.weight(f.cone(ev[("out", "B", "X")])), f.weight(f.cone(ev[("out", "B", "Z")]))
    check("cones: after a full correction a block carries the correction's own sources only",
          after == ((3, {("plus", "X"): 1, ("zero", "X"): 1}), (2, {("zero", "Z"): 1})), str(after))
    check("cones: the correction's event sees what came before",
          f.weight(f.cone(ev[("ec", "B", "X")])) == (3, {("plus", "X"): 1})
          and f.weight(f.cone(ev[("ec", "B", "Z")])) == (4, {("plus", "Z"): 1, ("zero", "Z"): 1}))
    Pstd = cones.Program([[("NEW", "std", "S")]], [], ["S"], "std")
    fs = cones.Flow(Pstd)
    es = {(e["kind"], e["sector"]): e for e in fs.events}
    idle = 1                                            # the fresh block idles in its first tick
    check("cones: a standard block is what a full correction leaves",
          fs.weight(fs.cone(es[("out", "X")])) == (3 + idle, {("plus", "X"): 1, ("zero", "X"): 1})
          and fs.weight(fs.cone(es[("out", "Z")])) == (2 + idle, {("zero", "Z"): 1}))
    P = cones.Program([[("CXg", "A", "B", "m")]], ["A", "B"], ["A", "B"], "g")
    refuses("cones: a control that is not measured yet is refused", lambda: cones.Flow(P), "not measured", exc=GadgetError)
    P = cones.Program([[("MZ", "C", "m")], [("CXg", "A", "B", "m")]], ["A", "B", "C"], ["A", "B"], "g")
    f = cones.Flow(P)
    check("cones: a controlled gate is one class of weight two",
          max(s["weight"] for s in f.sources) == 2 and sum(s["weight"] == 2 for s in f.sources) == 1)
    for name, layers, needle in (("a block used twice", [[("CX", "A", "B"), ("MZ", "A")]], "twice"),
                                 ("an interrupted correction", [[("ECa", "A")], [("CX", "A", "B")]], "middle"),
                                 ("a second half alone", [[("ECb", "A")]], "without its first half"),
                                 ("a dead block", [[("MZ", "A")], [("CX", "A", "B")]], "not live"),
                                 ("an unknown operation", [[("T", "A")]], "unknown")):
        refuses(f"cones: {name} is refused", lambda layers=layers: cones.Flow(cones.Program(layers, ["A", "B"], ["A", "B"])),
                needle, exc=GadgetError)
    refuses("cones: outputs that are not the live blocks are refused",
            lambda: cones.Flow(cones.Program([[("MZ", "A")]], ["A", "B"], ["A", "B"])), "declared outputs", exc=GadgetError)
    # ---- what a schedule spends, and how it follows its steps by corrections
    two = cones.Program([[("NEW", "zero2", "Q"), ("NEW", "and2", "a", "b", "t")], [("ECa", "Q"), ("MZ", "a", "m")],
                         [("ECb", "Q"), ("DROP", "b")], [("EC2a", "Q")], [("EC2b", "Q")], [("ECa", "Q")], [("ECb", "Q")]],
                        [], ["Q", "t"], "ledger")
    led = cones.ledger(cones.Flow(two))
    # live blocks at the end of each layer: Q a b t | Q b t | Q t | Q t ~p | Q t | Q t | Q t
    check("cones: the ledger counts live blocks per layer, corrections and deliveries",
          led == dict(layers=7, block_ticks=4 + 3 + 2 + 3 + 2 + 2 + 2, peak_blocks=4, inner_corrections=2, outer_corrections=1,
                      resources={"and2": 1, "plus2": 1, "zero2": 2}), str(led))
    one = cones.Program([[("NEW", "zero", "A", "B"), ("NEW", "std", "C")], [("CX", "A", "B")], [("MZ", "B"), ("ECa", "A")],
                         [("ECb", "A")]], [], ["A", "C"], "one")
    led1 = cones.ledger(cones.Flow(one))
    check("cones: inner resources are counted per block",
          led1["resources"] == {"std": 1, "zero": 2} and led1["block_ticks"] == 3 + 3 + 2 + 2 and led1["inner_corrections"] == 1,
          str(led1))
    rnd = random.Random(4)
    progs = [_random_program(rnd, rnd.randint(3, 5), rnd.randint(6, 14)) for _ in range(30)] + [two, one]
    check("cones: the count from the layers alone is the ledger of the flow",
          all(cones.count(p) == cones.ledger(cones.Flow(p)) for p in progs))
    strict = cones.Program([[("CX", "A", "B")], [("ECa", "A"), ("ECa", "B")], [("ECb", "A"), ("ECb", "B")], [("MZ", "B")],
                            [("ECa", "A")], [("ECb", "A")]], ["A", "B"], ["A"], "strict")
    pat = cones.correction_pattern(strict)
    # steps: the gate on A and on B (corrected next), the idle tick of A while B is measured (corrected next)
    check("cones: a schedule that corrects after every gate and idle tick has runs of length one only",
          pat["corrects_after_every_step"] and pat["steps"] == 3 and pat["corrected"] == 3 and pat["runs"] == {"1": 3}
          and pat["gates"] == 2 and pat["idle_ticks"] == 1 and pat["longest_run"] == 1
          and pat["gates_not_followed_by_correction"] == 0, str(pat))
    loose = cones.Program([[("CX", "A", "B")], [("CX", "A", "B")], [], [("MZ", "B"), ("ECa", "A")], [("ECb", "A")], []],
                          ["A", "B"], ["A"], "loose")
    pat = cones.correction_pattern(loose)
    # A: gate, gate, idle, then corrected, then one idle tick at the end; B: gate, gate, idle, then measured
    check("cones: uncorrected steps, runs and gates without a correction are counted",
          not pat["corrects_after_every_step"] and pat["steps"] == 7 and pat["uncorrected"] == 4 and pat["corrected"] == 1
          and pat["measured"] == 1 and pat["ended"] == 1 and pat["gates_not_followed_by_correction"] == 4
          and pat["runs"] == {"1": 1, "3": 2} and pat["longest_run"] == 3, str(pat))
    # ---- algebra against the rule
    P = cones.Program([[("CX", "A", "C")], [("CX", "A", "B")], [("CX", "A", "B")]], ["A", "B", "C"],
                      ["A", "B", "C"], "twice")
    f = cones.Flow(P)
    c = cones.contained(P, f, samples=2)
    check("cones: two equal gates cancel in the algebra and not in the rule, which only overestimates",
          not c["outside"] and c["events_where_rule_exceeds_algebra"] >= 1, str(c))
    rnd = random.Random(21)
    progs = [_random_program(rnd, rnd.randint(3, 5), rnd.randint(6, 14)) for _ in range(40)]
    bad = [p for p in progs if cones.contained(p, cones.Flow(p), samples=4)["outside"]]
    check("cones: on forty random programs linear propagation never leaves a cone", not bad
          and sum(len(cones.Flow(p).events) for p in progs) > 400)
    T = cones.Program([[("TOF", 2, "A", "B")], [("ECa", "A")], [("ECb", "A")]], ["A", "B"], ["A", "B"], "tof")
    f = cones.Flow(T)
    ev = {(e["kind"], str(e["block"]), e["sector"]): e for e in f.events}
    check("cones: a non-Clifford location pools both sectors of its blocks",
          f.weight(f.cone(ev[("out", "B", "Z")]))[0] == 2 + 2 and not cones.contained(T, f, 2)["outside"])
    T = cones.Program([[("ECa", "A")], [("ECb", "A")], [("TOF", 1, "A", "B")]], ["A", "B"], ["A", "B"], "tof2")
    f = cones.Flow(T)
    bx = next(e for e in f.events if e["kind"] == "out" and str(e["block"]) == "B" and e["sector"] == "X")
    check("cones: after a non-Clifford location a Z-only source of one block is in the X cone of the other",
          any(f.sources[i]["what"] == ("ecb", "m2") for i in cones._bits(f.cone(bx))))
    M = cones.Program([[("CX", "A", "B")], [("ECa", "A"), ("NEW", "plus", "C")], [("ECb", "A"), ("CX", "C", "B")],
                       [("MZ", "B", "m"), ("CZ", "A", "C")], [("CXg", "A", "C", "m")], [("ECa", "C")], [("ECb", "C")],
                       [("MX", "C")]], ["A", "B"], ["A"], "mixed")
    c = cones.contained(M, cones.Flow(M), samples=4)
    check("cones: where nothing cancels, the algebra reaches every source of every cone",
          not c["outside"] and c["events_where_rule_exceeds_algebra"] == 0 and c["events"] == 8, str(c))
    # ---- the physical test
    small = [("mixed", [[("CX", "A", "B")], [("ECa", "A"), ("NEW", "plus", "C")], [("ECb", "A"), ("CX", "C", "B")],
                        [("MZ", "B", "m"), ("CZ", "A", "C")], [("CXg", "A", "C", "m")], [("ECa", "C")], [("ECb", "C")],
                        [("MX", "C")]], ["A", "B"], ["A"]),
             ("std", [[("NEW", "std", "D"), ("CX", "A", "B")], [("CZ", "D", "A")], [("MZ", "B", "m")],
                      [("CZg", "A", "D", "m")], [("ECa", "D")], [("ECb", "D")]], ["A", "B"], ["A", "D"])]
    for name, layers, ins, outs in small:
        P = cones.Program(layers, ins, outs, name)
        for label, acts in (("on", None), ("off", lambda t, op: False)):
            d = cones.Direct(P, steane, acts)
            r = d.test(exhaustive=1, samples=2500, seed=3)
            check(f"cones, physical ({name}, controls {label}): no decoding goes wrong while every cone has at "
                  f"most t dirty coordinates", r["violations"] == 0 and r["within_hypothesis"] > 900
                  and r["largest_set"] >= 6, str({k: v for k, v in r.items() if k not in ("beyond",)}))
            check(f"cones, physical ({name}, controls {label}): with t + 1 dirty coordinates a decoding fails",
                  r["beyond"] is not None and r["beyond"]["dirty"] == r["t"] + 1)
    # a cone made too small is caught by the physical test
    P = cones.Program(small[0][1], small[0][2], small[0][3], "tampered")
    f = cones.Flow(P)
    k = next(e for e in f.events if e["kind"] == "meas")
    gate = next(i for i in cones._bits(f.cone(k)) if f.sources[i]["what"][0] == "gate")
    f.cone_of[k["node"]] &= ~(1 << gate)
    r = cones.Direct(P, steane, flow=f).test(exhaustive=2, samples=0, search=0)
    check("cones, physical: a cone with one gate removed is caught", r["violations"] > 0)
    if heavy:
        P = cones.Program(small[0][1], small[0][2], small[0][3], "golay")
        r = cones.Direct(P, library.code_of("golay")).test(exhaustive=1, samples=6000, seed=4)
        check("cones, physical (Golay, t = 3): no violation, and a failure with four dirty coordinates",
              r["violations"] == 0 and r["within_hypothesis"] > 4000 and r["beyond"] is not None, str(r)[:300])
    for _ in range(6 if heavy else 3):
        P = _random_program(rnd, 3, 7)
        if any(op[0] == "NEW" and op[1] not in ("zero", "plus", "std") for l in P.layers for op in l):
            continue
        acts = (lambda t, op, m={}: m.setdefault((t, op), rnd.random() < 0.5))
        r = cones.Direct(P, steane, acts).test(exhaustive=1, samples=1200, seed=5, search=300)
        check("cones, physical: a random program with random controls", r["violations"] == 0 and r["sets"] > 500,
              str(r)[:300])
    # ---- two levels
    P2 = cones.Program([[("NEW", "zero2", "Q")], [("ECa", "Q")], [("ECb", "Q")], [("EC2a", "Q")], [("EC2b", "Q")],
                        [("ECa", "Q")], [("ECb", "Q")]], [], ["Q"], "two")
    f = cones.Flow(P2)
    check("cones, two levels: recognised", f.two and len(f.outer) == 4 and len(f.terms) == 3)
    tags = [(o["tag"], o["sector"], len([i for i in cones._bits(f.cone_of[o["node"]])
                                         if f.sources[i]["kind"] == "event"])) for o in f.outer]
    first_x = [i for i, e in enumerate(f.events) if e["kind"] == "ec" and e["sector"] == "X"][0]
    out_x = next(o for o in f.outer if o["tag"] == "out" and o["sector"] == "X")
    mem = {f.sources[i]["event"] for i in cones._bits(f.cone_of[out_x["node"]]) if f.sources[i]["kind"] == "event"}
    check("cones, two levels: outer correction is a sink for the inner events before it",
          first_x not in mem and len(mem) == 3, f"{tags} {sorted(mem)}")
    ec2_x = next(o for o in f.outer if o["tag"] == "ec2" and o["sector"] == "X")
    mem2 = {f.sources[i]["event"] for i in cones._bits(f.cone_of[ec2_x["node"]]) if f.sources[i]["kind"] == "event"}
    check("cones, two levels: and its own event sees them", first_x in mem2)
    refuses("cones: an outer resource in a single-level run is refused",
            lambda: cones.Flow(cones.Program([[("NEW", "zero2", "Q")]], [], ["Q"]), two=False), "outer resource",
            exc=GadgetError)
    refuses("cones: the physical test refuses a two-level program",
            lambda: cones.Direct(P2, steane), "single-level", exc=GadgetError)
    # ---- receipts
    P = cones.Program(small[0][1], small[0][2], small[0][3], "receipt")
    rec = cones.build(P)
    check("cones receipt: verifies, and lists every event", len(cones.verify(rec)) == 3
          and len(rec["inner_events"]) == len(cones.Flow(P).events) and "largest cone" in cones.summary(rec))
    lite = cones.build(P, events=False)
    check("cones receipt: without the event list it still verifies and carries the same hash",
          "inner_events" not in lite and lite["events_sha256"] == rec["events_sha256"] and cones.verify(lite))
    for name, edit in (("an edited weight", lambda r: r["inner_events"][0].update(weight=1)),
                       ("edited statistics", lambda r: r["stats"].update(largest_weight=1)),
                       ("an edited theorem", lambda r: r.update(theorem="trust me")),
                       ("a dropped term", lambda r: r.update(top_level_terms=[{"x": 1}]))):
        bad = copy.deepcopy(rec)
        edit(bad)
        refuses(f"cones receipt: {name} is refused", lambda bad=bad: cones.verify(bad), "does not recompute", exc=GadgetError)
    bad = copy.deepcopy(rec)
    bad["program"]["layers"][1] = [["NEW", "plus", "C"]]
    refuses("cones receipt: an edited program under the old hash is refused", lambda: cones.verify(bad), "sha256",
            exc=GadgetError)
    rt = cones.Program.from_json(json.loads(json.dumps(P2.to_json())))
    check("cones: a two-level program survives a JSON round trip", rt.sha256() == P2.sha256()
          and cones.build(rt)["stats"] == cones.build(P2)["stats"])


def ftqc_codes_concat():
    import ftqc_codes
    from ftqc_codes.codes import concatenate
    return concatenate(ftqc_codes.load("steane"), ftqc_codes.load("steane"))


# ---------------------------------------------------------------- 8. command line
def t_cli():
    env = dict(os.environ, PYTHONPATH=os.path.dirname(HERE))
    with tempfile.TemporaryDirectory() as tmp:
        def run(*args):
            p = subprocess.run([sys.executable, "-m", "ftqc_gadgets", *args], env=env, cwd=tmp,
                               capture_output=True, text=True)
            return p.returncode, p.stdout + p.stderr
        rc, out = run("analyze", "knill:steane,anc=spec", "--check", "ec:1", "--close", "1", "-o", "t.json")
        check("cli: analyze writes a receipt and exits 0", rc == 0 and "holds through 1" in out, out)
        rc, out = run("verify", "t.json", "--full")
        check("cli: verify --full", rc == 0 and "re-enumerated" in out, out)
        rc, out = run("analyze", "knill:steane,anc=raw", "--check", "ec:1")
        check("cli: a failing property exits 1 and prints the fault", rc == 1 and "fails at" in out, out)
        for role, seq in (("plus", "ZX"), ("zero", "XZ")):
            rc, out = run("analyze", f"prep:steane,{role},tree={seq}", "--refines",
                          f"prep:steane,{role},spec:1", "-o", f"{role}.json")
            check(f"cli: refinement receipt for {role}", rc == 0, out)
        rc, out = run("compose", "t.json", "plus=plus.json", "zero=zero.json", "-o", "c.json")
        check("cli: compose", rc == 0 and "holds through 1" in out, out)
        rc, out = run("equiv", "fc:cat:4", "fc:cat:4,verify", "--order", "1")
        check("cli: equiv reports non-equivalence with exit 1", rc == 1 and "unmatched" in out, out)
        rc, out = run("equiv", "cnot:steane", "cnot:steane", "--order", "1")
        check("cli: equiv of a gadget with itself", rc == 0, out)
        rc, out = run("analyze", "nonsense:steane")
        check("cli: unknown gadget exits 2", rc == 2, out)
        rc, out = run("export", "cnot:rep3", "-o", "g.json")
        rc2, out2 = run("analyze", "g.json", "--check", "gate:0")
        check("cli: export and analyze a gadget file", rc == 0 and rc2 == 0, out + out2)


def main():
    import unittest
    suite = unittest.defaultTestLoader.loadTestsFromNames([
        "ftqc_gadgets.test_surgery", "ftqc_gadgets.test_surgery_build",
        "ftqc_gadgets.test_surgery_geometry", "ftqc_gadgets.test_surgery_operational",
        "ftqc_gadgets.test_occupied_transfer",
    ])
    check("signed surgery regressions", unittest.TextTestRunner().run(suite).wasSuccessful())
    heavy = "--heavy" in sys.argv
    t0 = time.time()
    for fn, args in ((t_structure, ()), (t_tables, ()), (t_dense, (heavy,)), (t_linear, ()),
                     (t_weights, ()), (t_props, (heavy,)), (t_refine, (heavy,)), (t_compose, ()),
                     (t_receipts, ()), (t_local, (heavy,)), (t_levels, (heavy,)), (t_factories, (heavy,)),
                     (t_cones, (heavy,)), (t_cli, ())):
        t1 = time.time()
        fn(*args)
        print(f"  {fn.__name__:<14} {time.time() - t1:6.2f} s   ({len(PASS)} passed so far)")
    print(f"ftqc_gadgets selftest: {len(PASS)} passed, {len(FAIL)} failed, {time.time() - t0:.1f} s")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
