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
    linear, units = (not an.has_causal_dec), []
    for off, mask, eff in an.dec:
        if not linear:
            break
        if not isinstance(eff, dict):
            linear = False                     # a decoder without a plain table
            break
        width = mask.bit_length()
        unit = [code_part(eff.get(1 << i, 0) >> an.shift) if (1 << i) in eff else None for i in range(width)]
        if any(u is None for u in unit):
            linear = False
            break
        for sy, e in eff.items():
            acc = [0] * len(funcs)
            for i in range(width):
                if (sy >> i) & 1:
                    acc = [a ^ b for a, b in zip(acc, unit[i])]
            if tuple(acc) != code_part(e >> an.shift):
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
    from .columns import ColumnError, layout
    try:
        pos, _ = layout(g)
    except ColumnError as err:
        result = dict(type="local", local=False, coordinates=n, reason=str(err))
        for loc in an.model.locations:
            if loc.op is err.op:
                result["fault"] = [loc.index, next(iter(loc.faults))]
                break
        return result
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

    def place(sig, coordinate):
        """Check the entire signature against its assigned physical coordinate."""
        e = sig >> an.shift
        dbits = sig >> an.nd
        c, trig = coordinate, []
        for k, ((off, mask), cm) in enumerate(zip(widths, colmaps)):
            sy = (dbits >> off) & mask
            if not sy:
                continue
            p = cm.get(sy)
            if p is None:
                return f"decoder {g.decoders[k].name} sees check values of no single error"
            if p != c:
                return "two coordinates are involved"
            trig.append(k)
            e ^= an.dec_per[k][p]
        if e == 0:
            return c, 0, trig
        return (c, e, trig) if e in O[c] else "the remaining error is not confined to that coordinate"

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
    # Equal signatures from distinct physical positions must pass separately.
    inputs = {i for b in g.inputs for i in b.in_locs}
    for loc in an.model.locations:
        if loc.index in inputs:
            continue
        coordinate = pos[loc.op.qs[0]]
        for mode, eff in loc.faults.items():
            result = place(an.sig(*eff), coordinate)
            if isinstance(result, str):
                return dict(type="local", local=False, coordinates=n, reason=result,
                            fault=[loc.index, mode])
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
    # The lightest-pattern decoder returns every weight <= r pattern exactly
    # iff no two such patterns share a syndrome, i.e. ker(H) has no nonzero
    # word of weight <= 2r. Re-prove this when a table was omitted. A
    # hierarchical decoder uses a different algorithm and needs its own bound.
    from ftqc_codes.hier import kernel_radius
    kernel_proofs = {}
    for index, (d, info) in enumerate(zip(g.decoders, an.dec_info)):
        r = radius.get(d.name) if isinstance(radius, dict) else radius
        supplied = r is not None
        bound = info["unique_through"]
        proof = "leader table"
        if d.rule is not None:
            # With at most t1 total errors every inner column decodes exactly.
            bound = an.dec[index][2].h.t1
            proof = "hierarchical inner leader tables"
        if not supplied:
            r = bound
        if type(r) is not int or not 0 <= r <= d.npos:
            return dict(type="local", local=False, coordinates=n,
                        reason=f"no valid integer radius for decoder {d.name}")
        evidence = {}
        if bound is None:
            key = (tuple(d.checks), d.npos, r)
            if key not in kernel_proofs:
                kernel_proofs[key] = kernel_radius(d.checks, d.npos, r, an.cap)
            bound, evidence = kernel_proofs[key]
            proof = "minimum weight of check kernel"
        if r > bound:
            return dict(type="local", local=False, coordinates=n,
                        reason=f"radius {r} exceeds proved radius {bound} for decoder {d.name}")
        row = dict(name=d.name, radius=r, source="given" if supplied else "table")
        if supplied or d.rule is not None:
            row["proof"] = dict(method=proof, **evidence)
        decs.append(row)
        through = r if through is None else min(through, r)
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


# ---------------------------------------------------------------------- coherent channels and malignant pairs
COHERENT_LEMMA = (
    "A Clifford gadget with zero postselected detectors defines a trace-preserving quantum "
    "instrument whose classical outcomes select deterministic Pauli-frame corrections.  When "
    "its noiseless Choi state on references and outputs is pure, the ideal map is unitary; when "
    "there are no output blocks, the ideal map is a destructive measurement channel.  If check "
    "holds through order t (and close holds through order s when output blocks carry a code), "
    "and for every fault location any two elementary Pauli fault modes producing the same "
    "classical measurement flips and the same output code syndrome also produce the exact same "
    "post-decoder interface effect, then any coherent linear combination of elementary Pauli "
    "faults at at most t locations acts on every in-filter input with residual error inside the "
    "declared filter and zero off-diagonal logical phase discrepancy (He, Nguyen and Pattison, "
    "Definition 4.8 and Lemma 4.14)."
)

OP_ROLE = {
    "R": "prep_zero",
    "RX": "prep_plus",
    "I": "idle",
    "H": "h",
    "S": "s",
    "SDG": "s",
    "CX": "cnot",
    "M": "measure_z",
    "MX": "measure_x",
}


def coherent_channel(an, kind, t, close_s=None, sector=False):
    """Verify unconditional (zero-postselection) channel simulation and off-diagonal coherence."""
    g = an.g
    if an.nd:
        return dict(
            type="coherent",
            coherent=False,
            kind=kind,
            t=int(t),
            weight="sector" if sector else "pauli",
            reason=f"gadget has {an.nd} postselected detectors",
        )
    if g.outputs and not an.pure:
        return dict(
            type="coherent",
            coherent=False,
            kind=kind,
            t=int(t),
            weight="sector" if sector else "pauli",
            reason="noiseless Choi state on references and outputs is not pure",
        )
    chk = check(an, kind, int(t), sector=sector)
    cls = None
    if close_s is not None and g.outputs:
        cls = close(an, int(close_s), sector=sector)
    W = an.weigher(sector=sector)
    code_funcs = [f for b in W.blocks for tp in "XZ" for f in b["keys"][tp]] if W.blocks and all(
        b["keys"] is not None for b in W.blocks
    ) else []

    def out_code_syn(eff_syn):
        return tuple((eff_syn & f).bit_count() & 1 for f in code_funcs)

    in_locs = {i for b in g.inputs for i in b.in_locs}
    mode_pairs = 0
    same_sector_pairs = 0
    coh_violation = None
    for loc in an.model.locations:
        if loc.index in in_locs:
            continue
        items = list(loc.faults.items())
        for (l1, ef1), (l2, ef2) in combinations(items, 2):
            mode_pairs += 1
            if ef1[2] == ef2[2]:
                e1 = an.effect(an.sig(*ef1))
                e2 = an.effect(an.sig(*ef2))
                if e1 < 0 or e2 < 0 or out_code_syn(e1) == out_code_syn(e2):
                    same_sector_pairs += 1
                    if e1 != e2 and coh_violation is None:
                        coh_violation = dict(location=loc.index, labels=[l1, l2])
    ok = (
        chk["holds_through"] >= int(t)
        and (cls is None or cls["holds_through"] >= int(close_s))
        and coh_violation is None
    )
    out = dict(
        type="coherent",
        coherent=bool(ok),
        kind=kind,
        t=int(t),
        close_s=None if cls is None else int(close_s),
        weight="sector" if sector else "pauli",
        detectors=an.nd,
        pure_choi=bool(an.pure),
        holds_through=chk["holds_through"],
        close_holds_through=None if cls is None else cls["holds_through"],
        close_method=None if cls is None else cls["method"],
        mode_pairs_checked=mode_pairs,
        same_branch_same_syndrome_pairs=same_sector_pairs,
        violation=chk["violation"] or (cls["violation"] if cls else None) or coh_violation,
    )
    if not ok:
        if chk["holds_through"] < int(t):
            out["reason"] = f"check({kind}) holds only through {chk['holds_through']} < {int(t)}"
        elif cls is not None and cls["holds_through"] < int(close_s):
            out["reason"] = f"close holds only through {cls['holds_through']} < {int(close_s)}"
        else:
            out["reason"] = "coherent Pauli modes on the same branch and syndrome differ logically"
    return out


def malignant_pairs(an, bound=1, sector=False, pair_slice=None):
    """Exact 2-location malignant support and rational depolarizing pair enumerator."""
    from fractions import Fraction

    g = an.g
    if an.nd:
        raise GadgetError("malignant_pairs requires a zero-postselection gadget (detectors == 0)")
    W = an.weigher(sector=sector)
    bound = int(bound)
    in_locs = {i for b in g.inputs for i in b.in_locs}
    body_locs = [loc for loc in an.model.locations if loc.index not in in_locs]
    L = len(body_locs)
    if pair_slice is not None:
        sk, sm = int(pair_slice[0]), int(pair_slice[1])
        if sm < 1 or sk < 0 or sk >= sm:
            raise GadgetError(f"invalid pair_slice {pair_slice!r}")
    else:
        sk, sm = 0, 1

    if W.weight(an.effect(0), bound) is None:
        return dict(
            type="malignant_pairs",
            t1_holds=False,
            bound=bound,
            weight="sector" if sector else "pauli",
            violation=dict(faults=[]),
        )

    sig_to_id = {0: 0}
    sigs = [0]
    first_mode_for_loc_sig = []
    loc_counts = []
    loc_nmodes = []
    loc_roles = []
    role_loc_counts = {}
    for loc in body_locs:
        role = OP_ROLE.get(loc.op.name, loc.op.name.lower())
        loc_roles.append(role)
        role_loc_counts[role] = role_loc_counts.get(role, 0) + 1
        cnt = {}
        fmode = {}
        nmodes = len(loc.faults)
        for label, eff in sorted(loc.faults.items()):
            s = an.sig(*eff)
            if W.weight(an.effect(s), bound) is None:
                return dict(
                    type="malignant_pairs",
                    t1_holds=False,
                    bound=bound,
                    weight="sector" if sector else "pauli",
                    violation=dict(faults=[[loc.index, label]]),
                )
            uid = sig_to_id.get(s)
            if uid is None:
                uid = len(sigs)
                sig_to_id[s] = uid
                sigs.append(s)
            cnt[uid] = cnt.get(uid, 0) + 1
            fmode.setdefault(uid, label)
        loc_counts.append(cnt)
        loc_nmodes.append(nmodes)
        first_mode_for_loc_sig.append(fmode)

    U = len(sigs)
    bad_mask = [0] * U
    for u in range(U):
        su = sigs[u]
        for v in range(u + 1, U):
            if W.weight(an.effect(su ^ sigs[v]), bound) is None:
                bad_mask[u] |= 1 << v
                bad_mask[v] |= 1 << u

    loc_sig_mask = [sum(1 << uid for uid in cnt) for cnt in loc_counts]
    bad_for_loc = []
    for cnt in loc_counts:
        bm = 0
        for uid in cnt:
            bm |= bad_mask[uid]
        bad_for_loc.append(bm)

    total_pairs = 0
    malignant_supp = 0
    supp_by_role = {}
    first_witness = None

    for i in range(sk, L, sm):
        bi = bad_for_loc[i]
        ri = loc_roles[i]
        total_pairs += L - 1 - i
        if not bi:
            continue
        for j in range(i + 1, L):
            overlap = bi & loc_sig_mask[j]
            if overlap:
                malignant_supp += 1
                rj = loc_roles[j]
                rkey = f"{ri},{rj}" if ri <= rj else f"{rj},{ri}"
                supp_by_role[rkey] = supp_by_role.get(rkey, 0) + 1
                if first_witness is None:
                    vj = (overlap & -overlap).bit_length() - 1
                    ui = next(u for u in loc_counts[i] if (bad_mask[u] >> vj) & 1)
                    first_witness = [
                        [body_locs[i].index, first_mode_for_loc_sig[i][ui]],
                        [body_locs[j].index, first_mode_for_loc_sig[j][vj]],
                    ]

    dep_by_role = {}
    if sm == 1:
        roles_sorted = sorted(role_loc_counts)
        W_role = {r: [Fraction(0)] * U for r in roles_sorted}
        diag_role = {r: Fraction(0) for r in roles_sorted}
        for i in range(L):
            r = loc_roles[i]
            nm = loc_nmodes[i]
            cnt = loc_counts[i]
            wr = W_role[r]
            for u, mu in cnt.items():
                wr[u] += Fraction(mu, nm)
            d_sum = 0
            for u, mu in cnt.items():
                bu = bad_mask[u]
                if not bu:
                    continue
                for v, mv in cnt.items():
                    if (bu >> v) & 1:
                        d_sum += mu * mv
            if d_sum:
                diag_role[r] += Fraction(d_sum, nm * nm)
        dep_total = Fraction(0)
        for ia, ra in enumerate(roles_sorted):
            wa = W_role[ra]
            quad_aa = Fraction(0)
            for u in range(U):
                if not wa[u] or not bad_mask[u]:
                    continue
                bu = bad_mask[u]
                s_v = Fraction(0)
                for v in range(U):
                    if (bu >> v) & 1 and wa[v]:
                        s_v += wa[v]
                quad_aa += wa[u] * s_v
            val_aa = (quad_aa - diag_role[ra]) / 2
            if val_aa:
                dep_by_role[f"{ra},{ra}"] = str(val_aa)
                dep_total += val_aa
            for rb in roles_sorted[ia + 1:]:
                wb = W_role[rb]
                quad_ab = Fraction(0)
                for u in range(U):
                    if not wa[u] or not bad_mask[u]:
                        continue
                    bu = bad_mask[u]
                    s_v = Fraction(0)
                    for v in range(U):
                        if (bu >> v) & 1 and wb[v]:
                            s_v += wb[v]
                    quad_ab += wa[u] * s_v
                if quad_ab:
                    dep_by_role[f"{ra},{rb}"] = str(quad_ab)
                    dep_total += quad_ab
    else:
        dep_total = Fraction(0)
        for i in range(sk, L, sm):
            bi = bad_for_loc[i]
            if not bi:
                continue
            ri = loc_roles[i]
            nmi = loc_nmodes[i]
            cnti = loc_counts[i]
            for j in range(i + 1, L):
                if not (bi & loc_sig_mask[j]):
                    continue
                rj = loc_roles[j]
                nmj = loc_nmodes[j]
                cntj = loc_counts[j]
                bij = 0
                for u, mu in cnti.items():
                    bu = bad_mask[u]
                    if not bu:
                        continue
                    for v, mv in cntj.items():
                        if (bu >> v) & 1:
                            bij += mu * mv
                term = Fraction(bij, nmi * nmj)
                dep_total += term
                rkey = f"{ri},{rj}" if ri <= rj else f"{rj},{ri}"
                dep_by_role[rkey] = dep_by_role.get(rkey, Fraction(0)) + term
        dep_by_role = {k: str(v) for k, v in sorted(dep_by_role.items()) if v}

    out = dict(
        type="malignant_pairs",
        t1_holds=True,
        bound=bound,
        weight="sector" if sector else "pauli",
        fault_locations=L,
        distinct_signatures=U,
        role_locations=dict(sorted(role_loc_counts.items())),
        total_location_pairs=total_pairs,
        malignant_support_pairs=malignant_supp,
        support_pairs_by_role=dict(sorted(supp_by_role.items())),
        depolarizing_pair_coeff=str(dep_total),
        depolarizing_pairs_by_role=dict(sorted(dep_by_role.items())),
        first_malignant_pair=first_witness,
    )
    if pair_slice is not None:
        out["pair_slice"] = [sk, sm]
    return out


def judge_malignant_pair(an, witness, bound=1, sector=False):
    """Re-evaluate a claimed 2-location malignant pair [(i, lab_i), (j, lab_j)] by direct replay."""
    if len(witness) != 2 or witness[0][0] == witness[1][0]:
        return False
    sg = an.replay([(int(f[0]), f[1]) for f in witness])
    if sg & an.DET:
        return False
    return an.weigher(sector=sector).weight(an.effect(sg), int(bound)) is None


def exrec_malignant_pairs(code, role, ec="knill", anc="goto", bound=1, sector=True, pair_slice=None):
    """Exact union malignant-pair census across all outgoing-wire truncations T of ExRec(role)."""
    from fractions import Fraction
    from . import library
    from .effects import Analysis

    truncs_for_role = {
        "prep_zero": [(), (0,)],
        "prep_plus": [(), (0,)],
        "idle": [(), (0,)],
        "h": [(), (0,)],
        "s": [(), (0,)],
        "cnot": [(), (0,), (1,), (0, 1)],
        "measure_z": [()],
        "measure_x": [()],
    }
    if role not in truncs_for_role:
        raise GadgetError(f"unknown exrec role {role!r}")
    trunc_list = truncs_for_role[role]
    analyses = []
    weighers = []
    loc_maps = []
    for tr in trunc_list:
        gt = library.exrec(code, role, ec=ec, anc=anc, trunc=tr, idles=True)
        ant = Analysis(gt)
        if ant.nd:
            raise GadgetError("exrec_malignant_pairs requires zero-postselection gadgets")
        wt = ant.weigher(sector=sector)
        in_set = {i for b in gt.inputs for i in b.in_locs}
        lmap = {(loc.tick, loc.op.name, loc.op.qs): loc for loc in ant.model.locations if loc.index not in in_set}
        analyses.append(ant)
        weighers.append(wt)
        loc_maps.append(lmap)

    an0 = analyses[0]
    in0 = {i for b in an0.g.inputs for i in b.in_locs}
    body_locs = [loc for loc in an0.model.locations if loc.index not in in0]
    L = len(body_locs)
    bound = int(bound)
    if pair_slice is not None:
        sk, sm = int(pair_slice[0]), int(pair_slice[1])
        if sm < 1 or sk < 0 or sk >= sm:
            raise GadgetError(f"invalid pair_slice {pair_slice!r}")
    else:
        sk, sm = 0, 1

    # Verify every truncation subset T <= Out(role) maps into ExRec_empty(role) and satisfies t=1.
    for t_idx, (tr, ant, wt, lmap) in enumerate(zip(trunc_list, analyses, weighers, loc_maps)):
        if not set(lmap) <= {(loc.tick, loc.op.name, loc.op.qs) for loc in body_locs}:
            raise GadgetError(f"truncation {tr} contains locations outside full ExRec({role})")
        if wt.weight(ant.effect(0), bound) is None:
            return dict(type="exrec_malignant_pairs", role=role, t1_holds=False,
                        truncation=list(tr), violation=dict(faults=[]))

    zero_tup = (0,) * len(trunc_list)
    sig_to_id = {zero_tup: 0}
    sigs = [zero_tup]
    first_mode_for_loc_sig = []
    loc_counts = []
    loc_nmodes = []
    loc_roles = []
    role_loc_counts = {}

    for loc in body_locs:
        key = (loc.tick, loc.op.name, loc.op.qs)
        rname = OP_ROLE.get(loc.op.name, loc.op.name.lower())
        loc_roles.append(rname)
        role_loc_counts[rname] = role_loc_counts.get(rname, 0) + 1
        cnt = {}
        fmode = {}
        nmodes = len(loc.faults)
        for label in sorted(loc.faults):
            tup = []
            for tr, ant, wt, lmap in zip(trunc_list, analyses, weighers, loc_maps):
                tloc = lmap.get(key)
                if tloc is None:
                    tup.append(0)
                else:
                    s_val = ant.sig(*tloc.faults[label])
                    if wt.weight(ant.effect(s_val), bound) is None:
                        return dict(
                            type="exrec_malignant_pairs",
                            role=role,
                            t1_holds=False,
                            truncation=list(tr),
                            violation=dict(faults=[[tloc.index, label]]),
                        )
                    tup.append(s_val)
            tup = tuple(tup)
            uid = sig_to_id.get(tup)
            if uid is None:
                uid = len(sigs)
                sig_to_id[tup] = uid
                sigs.append(tup)
            cnt[uid] = cnt.get(uid, 0) + 1
            fmode.setdefault(uid, label)
        loc_counts.append(cnt)
        loc_nmodes.append(nmodes)
        first_mode_for_loc_sig.append(fmode)

    U = len(sigs)
    bad_mask = [0] * U
    n_trunc = len(trunc_list)
    for u in range(U):
        tu = sigs[u]
        for v in range(u + 1, U):
            tv = sigs[v]
            is_bad = False
            for k in range(n_trunc):
                if weighers[k].weight(analyses[k].effect(tu[k] ^ tv[k]), bound) is None:
                    is_bad = True
                    break
            if is_bad:
                bad_mask[u] |= 1 << v
                bad_mask[v] |= 1 << u

    loc_sig_mask = [sum(1 << uid for uid in cnt) for cnt in loc_counts]
    bad_for_loc = []
    for cnt in loc_counts:
        bm = 0
        for uid in cnt:
            bm |= bad_mask[uid]
        bad_for_loc.append(bm)

    total_pairs = 0
    malignant_supp = 0
    supp_by_role = {}
    first_witness = None
    first_witness_trunc = None

    for i in range(sk, L, sm):
        bi = bad_for_loc[i]
        ri = loc_roles[i]
        total_pairs += L - 1 - i
        if not bi:
            continue
        for j in range(i + 1, L):
            overlap = bi & loc_sig_mask[j]
            if overlap:
                malignant_supp += 1
                rj = loc_roles[j]
                rkey = f"{ri},{rj}" if ri <= rj else f"{rj},{ri}"
                supp_by_role[rkey] = supp_by_role.get(rkey, 0) + 1
                if first_witness is None:
                    vj = (overlap & -overlap).bit_length() - 1
                    ui = next(u for u in loc_counts[i] if (bad_mask[u] >> vj) & 1)
                    tu, tv = sigs[ui], sigs[vj]
                    win_k = next(k for k in range(n_trunc)
                                 if weighers[k].weight(analyses[k].effect(tu[k] ^ tv[k]), bound) is None)
                    ki = (body_locs[i].tick, body_locs[i].op.name, body_locs[i].op.qs)
                    kj = (body_locs[j].tick, body_locs[j].op.name, body_locs[j].op.qs)
                    first_witness = [
                        [loc_maps[win_k][ki].index, first_mode_for_loc_sig[i][ui]],
                        [loc_maps[win_k][kj].index, first_mode_for_loc_sig[j][vj]],
                    ]
                    first_witness_trunc = list(trunc_list[win_k])

    dep_by_role = {}
    if sm == 1:
        roles_sorted = sorted(role_loc_counts)
        W_role = {r: [Fraction(0)] * U for r in roles_sorted}
        diag_role = {r: Fraction(0) for r in roles_sorted}
        for i in range(L):
            r = loc_roles[i]
            nm = loc_nmodes[i]
            cnt = loc_counts[i]
            wr = W_role[r]
            for u, mu in cnt.items():
                wr[u] += Fraction(mu, nm)
            d_sum = 0
            for u, mu in cnt.items():
                bu = bad_mask[u]
                if not bu:
                    continue
                for v, mv in cnt.items():
                    if (bu >> v) & 1:
                        d_sum += mu * mv
            if d_sum:
                diag_role[r] += Fraction(d_sum, nm * nm)
        dep_total = Fraction(0)
        for ia, ra in enumerate(roles_sorted):
            wa = W_role[ra]
            quad_aa = Fraction(0)
            for u in range(U):
                if not wa[u] or not bad_mask[u]:
                    continue
                bu = bad_mask[u]
                s_v = Fraction(0)
                for v in range(U):
                    if (bu >> v) & 1 and wa[v]:
                        s_v += wa[v]
                quad_aa += wa[u] * s_v
            val_aa = (quad_aa - diag_role[ra]) / 2
            if val_aa:
                dep_by_role[f"{ra},{ra}"] = str(val_aa)
                dep_total += val_aa
            for rb in roles_sorted[ia + 1:]:
                wb = W_role[rb]
                quad_ab = Fraction(0)
                for u in range(U):
                    if not wa[u] or not bad_mask[u]:
                        continue
                    bu = bad_mask[u]
                    s_v = Fraction(0)
                    for v in range(U):
                        if (bu >> v) & 1 and wb[v]:
                            s_v += wb[v]
                    quad_ab += wa[u] * s_v
                if quad_ab:
                    dep_by_role[f"{ra},{rb}"] = str(quad_ab)
                    dep_total += quad_ab
    else:
        dep_total = Fraction(0)
        for i in range(sk, L, sm):
            bi = bad_for_loc[i]
            if not bi:
                continue
            ri = loc_roles[i]
            nmi = loc_nmodes[i]
            cnti = loc_counts[i]
            for j in range(i + 1, L):
                if not (bi & loc_sig_mask[j]):
                    continue
                rj = loc_roles[j]
                nmj = loc_nmodes[j]
                cntj = loc_counts[j]
                bij = 0
                for u, mu in cnti.items():
                    bu = bad_mask[u]
                    if not bu:
                        continue
                    for v, mv in cntj.items():
                        if (bu >> v) & 1:
                            bij += mu * mv
                term = Fraction(bij, nmi * nmj)
                dep_total += term
                rkey = f"{ri},{rj}" if ri <= rj else f"{rj},{ri}"
                dep_by_role[rkey] = dep_by_role.get(rkey, Fraction(0)) + term
        dep_by_role = {k: str(v) for k, v in sorted(dep_by_role.items()) if v}

    out = dict(
        type="exrec_malignant_pairs",
        role=role,
        ec=ec,
        anc=anc,
        truncations=[list(t) for t in trunc_list],
        t1_holds=True,
        bound=bound,
        weight="sector" if sector else "pauli",
        qubits=an0.g.circuit.n,
        fault_locations=L,
        distinct_tuple_signatures=U,
        role_locations=dict(sorted(role_loc_counts.items())),
        total_location_pairs=total_pairs,
        malignant_support_pairs=malignant_supp,
        support_pairs_by_role=dict(sorted(supp_by_role.items())),
        depolarizing_pair_coeff=str(dep_total),
        depolarizing_pairs_by_role=dict(sorted(dep_by_role.items())),
        first_malignant_pair=first_witness,
        first_malignant_truncation=first_witness_trunc,
    )
    if pair_slice is not None:
        out["pair_slice"] = [sk, sm]
    return out
