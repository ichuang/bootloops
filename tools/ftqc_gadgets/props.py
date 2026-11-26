"""Properties of a gadget over "at most s faults and input weight at most r".

check(kind, t)   For every input error of weight r and every accepted set of s
                 faults with r + s <= t, the effect is within weight `bound` of zero
                 on every output block:
                     prep   no inputs,            bound = s
                     gate   bound = r + s         (errors may pass through, none are made worse)
                     ec     bound = s             (errors that came in are gone)
                 "Within weight w of zero" includes the logical action, because the
                 references are part of the state the effect is measured on.
close(s)         For ANY input, every accepted set of j <= s faults leaves every
                 output block within weight j of its code space (logical value not
                 examined).  With `ec` this is the pair of error-correction
                 properties of Aliferis, Gottesman and Preskill.
refine(impl, spec, order, r)
                 Every accepted set of j faults in impl, together with any input error
                 of weight w <= r, w + j <= order, has exactly the effect of some
                 accepted set of at most j faults in spec with the same input error.  Two
                 gadgets that refine each other are fault-equivalent to that order
                 (after Rodatz, Poor and Kissinger, truncated in weight).
"""
from itertools import combinations

from .effects import Budget
from .gadget import GadgetError

BOUND = dict(prep=lambda r, s: s, gate=lambda r, s: r + s, ec=lambda r, s: s)


def check(an, kind, t, sector=False):
    if kind not in BOUND:
        raise GadgetError(f"unknown property {kind!r} (prep, gate or ec)")
    if kind == "prep" and an.g.inputs:
        raise GadgetError("a gadget with inputs is not a preparation")
    if kind != "prep" and not an.g.inputs:
        raise GadgetError(f"property {kind!r} needs an input block")
    W = an.weigher(sector)
    r_max = 0 if kind == "prep" else t
    per, violation, note = [], None, None
    try:
        ins = an.input_errors(r_max, sector)
        for total in range(t + 1):
            for r in range(min(total, r_max) + 1):
                s = total - r
                bound = BOUND[kind](r, s)
                groups = {}
                for sg, (w, lab) in ins.items():
                    if w == r:
                        groups.setdefault(sg & an.DET, []).append((sg, lab))
                cases = bad = worst = 0
                for d, items in groups.items():
                    level = an.levels(s, d)[s]
                    for sg, lab in items:
                        for v in level:
                            w = W.weight(an.effect(v ^ sg), bound)
                            cases += 1
                            if w is None:
                                bad += 1
                                if violation is None:
                                    violation = dict(r=r, s=s, bound=bound, input_error=[list(x) for x in lab],
                                                     faults=an.faults_of(an.find(s, d, v)))
                            elif w > worst:
                                worst = w
                per.append(dict(r=r, s=s, bound=bound, cases=cases, violations=bad, worst_weight=worst))
            if violation:
                break
    except Budget as e:
        note = f"stopped: {e}"
    done = {}
    for p in per:
        done.setdefault(p["r"] + p["s"], []).append(p)
    holds = -1
    for total in range(t + 1):
        rows = done.get(total, [])
        need = min(total, r_max) + 1
        if len(rows) < need or any(p["violations"] for p in rows):
            break
        holds = total
    return dict(kind=kind, t=t, weight="sector" if sector else "pauli", holds_through=holds,
                per=per, violation=violation, note=note)


def judge(an, kind, witness, sector=False):
    """Re-evaluate a claimed violation by replay.  True iff it is one."""
    sg = an.replay([tuple(f) for f in witness["faults"]], [tuple(x) for x in witness["input_error"]])
    if sg & an.DET:
        return False
    r = _weight(witness["input_error"], sector)
    s = len(witness["faults"])
    if r != witness["r"] or s != witness["s"]:
        return False
    return an.weigher(sector).weight(an.effect(sg), BOUND[kind](r, s)) is None


def _weight(labels, sector):
    if not sector:
        return len({(b, p) for b, p, _ in labels})
    x = {(b, p) for b, p, t in labels if t in "XY"}
    z = {(b, p) for b, p, t in labels if t in "ZY"}
    return max(len(x), len(z))


# ---------------------------------------------------------------------- closeness
def close(an, s_max, sector=False, span_cap=1 << 16):
    """For any input: accepted sets of j <= s_max faults leave each output block within
    weight j of its code space."""
    W = an.weigher(sector)
    if any(b["keys"] is None for b in W.blocks):
        raise GadgetError("closeness to the code space needs a code on every output block")
    funcs = [f for b in W.blocks for t in "XZ" for f in b["keys"][t]]

    def code_part(e):
        return tuple((e & f).bit_count() & 1 for f in funcs)

    def block_weight(e, upto):
        worst = 0
        for blk in W.blocks:
            if blk["built"] < upto:
                W._extend(blk, upto)
            if sector:
                hits = [blk["tabs"][t].get(W._key(e, blk["keys"][t])) for t in "XZ"]
            else:
                kx, kz = blk["keys"]["X"], blk["keys"]["Z"]
                hits = [blk["tabs"]["P"].get(W._key(e, kx) | (W._key(e, kz) << len(kx)))]
            if any(h is None or h[0] > upto for h in hits):
                return None
            worst = max([worst] + [h[0] for h in hits])
        return worst

    # The code syndrome of the output after the decoders.  A decoder's correction has
    # a code syndrome that is linear in its check values when (as for syndrome
    # extraction) the checks are the block's own stabilizers; then the whole thing is a
    # linear function of the signature and can be tested on the input columns alone.
    linear, units = True, []
    for off, mask, eff in an.dec:
        if not isinstance(eff, dict):
            linear = False                     # a decoder without a plain table
            break
        width = mask.bit_length()
        unit = [code_part(eff.get(1 << i, 0)) if (1 << i) in eff else None for i in range(width)]
        if any(u is None for u in unit):
            linear = False
            break
        for sy, e in eff.items():
            acc = [0] * len(funcs)
            for i in range(width):
                if (sy >> i) & 1:
                    acc = [a ^ b for a, b in zip(acc, unit[i])]
            if tuple(acc) != code_part(e):
                linear = False
                break
        units.append((off, unit))

    def total_code_part(sg):
        acc = list(code_part(sg >> an.shift))
        d = sg >> an.nd
        for off, unit in units:
            for i, u in enumerate(unit):
                if (d >> (off + i)) & 1:
                    acc = [a ^ b for a, b in zip(acc, u)]
        return acc

    in_cols = [s for row in an.in_sig for pair in row for s in pair]
    independent = linear and all(not (s & an.DET) and not any(total_code_part(s)) for s in in_cols)
    if not linear and any(not isinstance(eff, dict) for _, _, eff in an.dec):
        return dict(holds_through=-1, per=[], violation=None, method=None,
                    note="not decided: a decoder of this gadget has no plain lookup table")
    per, violation, note = [], None, None
    try:
        if independent:
            method = "the code syndrome of the output after the decoders does not depend on the input"
            inputs = [0]
        else:
            from ftqc_codes import gf2
            basis = gf2.Echelon(in_cols).basis
            if 1 << len(basis) > span_cap:
                return dict(holds_through=-1, per=[], violation=None, method=None,
                            note=f"not decided: the output code syndrome depends on the input and "
                                 f"the input errors span 2^{len(basis)} cases")
            method = f"all 2^{len(basis)} input error signatures enumerated"
            inputs = [0]
            for b in basis:
                inputs += [x ^ b for x in inputs]
        for s in range(s_max + 1):
            cases = bad = 0
            by_det = {}
            for sg in inputs:
                by_det.setdefault(sg & an.DET, []).append(sg)
            for d, items in by_det.items():
                level = an.levels(s, d)[s]
                for sg in items:
                    for v in level:
                        e = an.effect(v ^ sg)
                        cases += 1
                        if e == -1 or block_weight(e, s) is None:
                            bad += 1
                            if violation is None:
                                violation = dict(s=s, input_signature=sg,
                                                 faults=an.faults_of(an.find(s, d, v)))
            per.append(dict(s=s, cases=cases, violations=bad))
            if violation:
                break
    except Budget as e:
        note = f"stopped: {e}"
    holds = -1
    for p in per:
        if p["violations"]:
            break
        holds = p["s"]
    return dict(holds_through=holds, per=per, violation=violation, method=method, note=note)


# ---------------------------------------------------------------------- refinement
def _input_pairs(a2, a1, r):
    """Distinct (signature in impl, signature in spec) over input Paulis of weight <= r,
    each with one Pauli that produces it."""
    cells = [(b, p) for b, row in enumerate(a2.in_sig) for p in range(len(row))]
    out = {(0, 0): []}
    for w in range(1, r + 1):
        for sub in combinations(cells, w):
            def rec(i, s2, s1, lab):
                if i == w:
                    out.setdefault((s2, s1), lab)
                    return
                b, p = sub[i]
                x2, z2 = a2.in_sig[b][p]
                x1, z1 = a1.in_sig[b][p]
                rec(i + 1, s2 ^ x2, s1 ^ x1, lab + [[b, p, "X"]])
                rec(i + 1, s2 ^ z2, s1 ^ z1, lab + [[b, p, "Z"]])
                rec(i + 1, s2 ^ x2 ^ z2, s1 ^ x1 ^ z1, lab + [[b, p, "Y"]])
            rec(0, 0, 0, [])
    return out


def comparable(a2, a1):
    if a2.g.shape() != a1.g.shape():
        return "the gadgets have different input or output blocks"
    if a2.group_id != a1.group_id:
        return "the gadgets do not have the same ideal map (stabilizer groups on references and outputs differ)"
    return None


def refine(impl, spec, order, r=0):
    """Budgets are totals: input weight w <= r and j faults with w + j <= order."""
    why = comparable(impl, spec)
    if why:
        raise GadgetError(why)
    per, violation, note = [], None, None
    n_pairs = 0
    try:
        pairs = _input_pairs(impl, spec, min(r, order))
        n_pairs = len(pairs)
        cost1 = {}                                   # spec input signature -> [{effect: least j}, levels done]
        for total in range(order + 1):
            n2 = bad = 0
            for (s2, s1), lab in pairs.items():
                w = len({(b, p) for b, p, _ in lab})
                j = total - w
                if j < 0:
                    continue
                c1 = cost1.setdefault(s1, [{}, 0])
                while c1[1] <= j:
                    for v in spec.levels(c1[1], s1 & spec.DET)[c1[1]]:
                        c1[0].setdefault(spec.effect(v ^ s1), c1[1])
                    c1[1] += 1
                d2 = s2 & impl.DET
                for v in impl.levels(j, d2)[j]:
                    e = impl.effect(v ^ s2)
                    n2 += 1
                    if c1[0].get(e, j + 1) > j:
                        bad += 1
                        if violation is None:
                            violation = dict(order=total, s=j, input_error=lab, effect=e,
                                             faults=impl.faults_of(impl.find(j, d2, v)))
            per.append(dict(order=total, impl_effects=n2, unmatched=bad))
            if violation:
                break
    except Budget as e:
        note = f"stopped: {e}"
    holds = -1
    for p in per:
        if p["unmatched"]:
            break
        holds = p["order"]
    return dict(order=order, input_weight=r, input_errors=n_pairs,
                holds_through=holds, per=per, violation=violation, note=note)


def judge_refine(impl, spec, witness, cap_order=None):
    """Re-evaluate a claimed refinement failure: the impl fault set is accepted, and no
    accepted spec set of at most as many faults has its effect (with the same input error)."""
    lab = [tuple(x) for x in witness["input_error"]]
    s2 = impl.replay([tuple(f) for f in witness["faults"]], lab)
    if s2 & impl.DET:
        return False
    e = impl.effect(s2)
    s1 = spec.replay([], lab)
    j = len(witness["faults"])
    for k in range(j + 1):
        for v in spec.levels(k, s1 & spec.DET)[k]:
            if spec.effect(v ^ s1) == e:
                return False
    return True


# ---------------------------------------------------------------------- locality
LEMMA = ("Call a signature local at coordinate c if no detector fires, each decoder's check "
         "values are zero or those of a single error at position c, and after removing that "
         "position's correction the syndrome is that of a Pauli on the output qubits at position c "
         "(every output block).  Suppose every fault of every location is local at some coordinate, "
         "and every input error at position p is local at p.  Let a fault set and an input error "
         "touch a set C of coordinates.  Signatures add, so each decoder sees the check values of "
         "an error pattern supported on C.  If |C| is at most the radius within which the decoder's "
         "lightest pattern is unique, the decoder returns exactly that pattern, its correction "
         "cancels what the pattern contributed, and what is left is the syndrome of a Pauli "
         "supported on C in each output block.  Hence: the output is within weight |C| of the "
         "ideal state, logical action included, however many faults lie inside the |C| coordinates "
         "(gate).  If in addition what is left of each input error is zero, only the coordinates "
         "of faults remain (ec), and since the output's code syndrome is then a linear function "
         "of the signature that vanishes on input errors, for any input it is the code syndrome "
         "of a Pauli on the fault coordinates (close).  With X and Z counted separately the same "
         "argument applies to each sector when no decoder responds to input errors of both types.")


def local(an, radius=None):
    """Coordinate-locality certificate.  radius: {decoder name: r} or one integer, the
    radius within which each decoder's lightest pattern is unique; by default the value
    found when the decoder's table was built.  Returns a dict; `local` is False with a
    witness when some fault or input error is not confined to one coordinate."""
    from ftqc_codes import gf2
    g = an.g
    sizes = {len(b.qubits) for b in g.inputs + g.outputs} | {d.npos for d in g.decoders}
    if len(sizes) != 1:
        return dict(type="local", local=False, reason="blocks and decoders do not share one set of coordinates")
    n = sizes.pop()
    if an.nd:
        return dict(type="local", local=False, coordinates=n,
                    reason="the gadget postselects; analyse its postselected parts as specifications")
    base = sum(len(b.refs) for b in g.inputs)
    O, OX, OZ = [], [], []
    for c in range(n):
        xs = [an.sx[base + k * n + c] for k in range(len(g.outputs))]
        zs = [an.sz[base + k * n + c] for k in range(len(g.outputs))]
        O.append(gf2.Echelon(xs + zs))
        OX.append(gf2.Echelon(xs))
        OZ.append(gf2.Echelon(zs))
    colmaps = []
    for d, cols in zip(g.decoders, an.dec_cols):
        cm = {}
        for p, v in enumerate(cols):
            if v == 0 or v in cm:
                return dict(type="local", local=False, coordinates=n,
                            reason=f"decoder {d.name} cannot tell single errors apart")
            cm[v] = p
        colmaps.append(cm)
    widths = [(off, (1 << w) - 1) for off, w in an.dec_layout]

    def place(sig, want=None):
        """(coordinate or None, residual, triggered decoder indices) or a reason string."""
        e = sig >> an.shift
        dbits = sig >> an.nd
        c, trig = want, []
        for k, ((off, mask), cm) in enumerate(zip(widths, colmaps)):
            sy = (dbits >> off) & mask
            if not sy:
                continue
            p = cm.get(sy)
            if p is None:
                return f"decoder {g.decoders[k].name} sees check values of no single error"
            if c is not None and p != c:
                return "two coordinates are involved"
            c = p
            trig.append(k)
            e ^= an.dec_per[k][p]
        if e == 0:
            return c, 0, trig
        if c is not None:
            return (c, e, trig) if e in O[c] else "the remaining error is not confined to that coordinate"
        for cc in range(n):
            if e in O[cc]:
                return cc, e, trig
        return "the remaining error is not confined to one coordinate"

    for i, sg in enumerate(an.cols):
        r = place(sg)
        if isinstance(r, str):
            return dict(type="local", local=False, coordinates=n, reason=r,
                        fault=[an.src[i][0], an.src[i][1]])
    gate = ec = True
    seen = [set() for _ in g.decoders]
    sector_gate = True
    for b, row in enumerate(an.in_sig):
        for p, (sx, sz) in enumerate(row):
            for t, sg, ech in (("X", sx, OX), ("Z", sz, OZ), ("Y", sx ^ sz, None)):
                r = place(sg, p)
                if isinstance(r, str):
                    return dict(type="local", local=False, coordinates=n, reason="input error: " + r,
                                input_error=[[b, p, t]])
                _, res, trig = r
                if res:
                    ec = False
                    if ech is not None and res not in ech[p]:
                        sector_gate = False
                if t != "Y":
                    for k in trig:
                        seen[k].add(t)
    sector = all(len(s) < 2 for s in seen) and sector_gate
    # the output's code syndrome is linear in the decoder check values iff the positions'
    # corrections respect every linear relation among the positions' check columns
    try:
        W = an.weigher()
        lin = all(b["keys"] is not None for b in W.blocks)
    except GadgetError:
        W, lin = None, False
    if lin:
        funcs = [f for b in W.blocks for t in "XZ" for f in b["keys"][t]]
        nf = len(funcs)
        for cols, per in zip(an.dec_cols, an.dec_per):
            aug = []
            for v, e in zip(cols, per):
                k = 0
                for i, f in enumerate(funcs):
                    if (e & f).bit_count() & 1:
                        k |= 1 << i
                aug.append(k | (v << nf))
            if gf2.rank(aug) != gf2.rank(cols):
                lin = False
    decs, through = [], None
    for d, info in zip(g.decoders, an.dec_info):
        r = radius.get(d.name) if isinstance(radius, dict) else radius
        src = "given"
        if r is None:
            r, src = info["unique_through"], "table"
        if r is None:
            return dict(type="local", local=False, coordinates=n,
                        reason=f"no radius for decoder {d.name}: its table was not built and none was given")
        if info["unique_through"] is not None and r > info["unique_through"]:
            return dict(type="local", local=False, coordinates=n,
                        reason=f"the radius given for decoder {d.name} exceeds the one its table shows")
        decs.append(dict(name=d.name, radius=int(r), source=src))
        through = int(r) if through is None else min(through, int(r))
    return dict(type="local", local=True, coordinates=n, fault_effects=len(an.cols), decoders=decs,
                through=through, gate=gate, ec=ec, close=bool(ec and lin), sector=sector)


def independent_of_input(an, names=None):
    """Which decoders' check values no input error can change (names: restrict to these).
    For error correction repeated on one block this is the statement that a later round
    is driven by the faults of the rounds before it and not by what the block arrived with."""
    out = {}
    for d, (off, w) in zip(an.g.decoders, an.dec_layout):
        if names is not None and d.name not in names:
            continue
        mask = ((1 << w) - 1) << (an.nd + off)
        hit = None
        for b, row in enumerate(an.in_sig):
            for p, (sx, sz) in enumerate(row):
                for t, sg in (("X", sx), ("Z", sz)):
                    if hit is None and sg & mask:
                        hit = [b, p, t]
        out[d.name] = dict(independent=hit is None, input_error=None if hit is None else [hit])
    return out
