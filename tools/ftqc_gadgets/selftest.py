#!/usr/bin/env python3
"""ftqc_gadgets battery.  Stdlib only (plus ftqc_audit and ftqc_codes from this repository).

    python3 tools/ftqc_gadgets/selftest.py
    python3 tools/ftqc_gadgets/selftest.py --heavy    # adds the Golay gadgets at t = 3

Exit 0 only if every check that ran passed.
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
sys.path.insert(0, os.path.dirname(HERE))

from ftqc_audit import analysis as fca  # noqa: E402
from ftqc_audit import library as fclib  # noqa: E402
from ftqc_audit.selftest import Dense, dense_run  # noqa: E402
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
    bad = copy.deepcopy(lv)
    bad["conclusion"]["bad_columns_tolerated"] = 2
    refuses("levels: an inflated conclusion is refused", lambda: levels.verify_reduction(bad, outer, inner),
            "conclusion")
    weak = dict(inner, cnot=certify.build(library.load("cnot:steane"), [("gate", 0)]))
    check("levels: a weaker inner receipt lowers the inner radius",
          levels.reduction(outer, weak)["conclusion"]["inner_faults_per_rectangle"] == 0)


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
    heavy = "--heavy" in sys.argv
    t0 = time.time()
    for fn, args in ((t_structure, ()), (t_tables, ()), (t_dense, (heavy,)), (t_linear, ()),
                     (t_weights, ()), (t_props, (heavy,)), (t_refine, (heavy,)), (t_compose, ()),
                     (t_receipts, ()), (t_local, (heavy,)), (t_levels, (heavy,)), (t_cli, ())):
        t1 = time.time()
        fn(*args)
        print(f"  {fn.__name__:<14} {time.time() - t1:6.2f} s   ({len(PASS)} passed so far)")
    print(f"ftqc_gadgets selftest: {len(PASS)} passed, {len(FAIL)} failed, {time.time() - t0:.1f} s")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
