"""Filter calculus under local-stochastic noise: weights of covers in place of probabilities.

The noise condition is Pr[all sources of G faulty] <= product of their weights, for every
finite set G of sources (physical locations of weight p).  Nothing else is assumed: no
independence, and no conditioning.

Covers.  An event is *covered with weight w* when it implies that all sources of some member
of a family of sets are faulty, the weights of the members summing to at most w.  Its
probability is then at most w, the weights of covers of events on disjoint sources multiply,
and the weights of alternatives add.

The filter rule gives every block of a purification tree, in each sector, a junk flag J and
a set D of dirty coordinates as functions of the faulty sources:

    leaf        J = some junk-capable location is faulty;  D_c = some local location at
                coordinate c is faulty.
    check of s  survivor A, partner B, fresh dirt f, step dirt M on the partner's side:
                J' = (J_A and J_B) or (neither, and |D_A + D_B + M| >= wmin);
                D' = f + D_B + M if only A is junk, f + D_A if only B is,
                     f + (D_A & (D_B + M)) if neither;
                other sector:  J' = J_A or J_B,  D' = D_A + D_B + f.

Contract.  A sector's contract is a junk weight pi, modes [(w, q)] of total weight at most
one, and a rate `last`: J is covered with weight pi, and for every nonempty set A of
coordinates the event "not J, and A inside D" is covered with weight
sum_k w_k (q_k + last)^|A|; the empty requirement has weight 1.  A list of total weight
W < 1 is padded with the mode (1 - W, 0).  Then the weight for a union of two dirt processes,
whose sources at different coordinates are disjoint, is the sum over pairs of modes of
w w' (q + q')^|A| (a sum over the ways of splitting A), and the weight for an intersection at
the same coordinates, on disjoint sources, is the sum of w w' (q q')^|A|.  A list of total
weight W > 1 is replaced by (w / W, W q): for |A| >= 1 this does not lower any weight.

Check.  With modes (w_i, q_i) of the survivor's dirt whose sources are disjoint from the
partner's side, modes (y_j, r_j) of the partner's side D_B + M, and the survivor's carried
rate l (zero when nothing is carried), the three cases of the rule in which the survivor is
not junk give

    inherited   (w_i y_j, q_i r_j) neither junk,  (pi_B w_i, q_i) only B,  (pi_A y_j, r_j) only A
    survivor    inherited united with l (carried) and with the fresh rate
    junk        pi_A pi_B + C(n, wmin) sum_ij w_i y_j (q_i + l + r_j)^wmin
    other       junk pi_A + pi_B, modes (w_i w_j, q_i + q_j + fresh)

Each line covers a non-monotone event of the rule by an increasing one: "only B junk and A
inside D_A" lies in (cover of J_B) and (cover of "not J_A, A inside D_A"), and so on.

Retries.  A check that rejects is repeated with fresh inputs, so the delivered block is the
first accepted attempt and attempts use disjoint sources.  An attempt is rejected only if
an incoming pattern of the checked sector is not trivial or a step location that can cause
rejection is faulty, an event covered with weight rho = nu_A + nu_B + rho_step, where
nu = pi + n * sum_k w_k (q_k + last).  An event of the delivered block is covered by the
union over k of (attempts 1..k-1 rejected) and (the event in attempt k), of weight
(weight of the event) / (1 - rho).  This is where the division by an acceptance bound comes
from; no probability is conditioned.

Every weight is exact in rationals and rounded upward, the bound 1 - rho downward.
"""
from fractions import Fraction
from math import comb

from . import filters as fl
from .stochastic import probability
from .supports import fields
from .tails import SIG, sci, up

MAX_MODES = fl.MAX_MODES
ONE = Fraction(1)

RULE = (
    "Local-stochastic filter calculus.  Every quantity is the total weight of a family of fault sets "
    "covering an event of the filter rule.  Weights of covers on disjoint sources multiply, weights of "
    "alternatives add, and the first accepted attempt of a check costs the factor 1 / (1 - rho) with "
    "rho the weight of rejection.  A contract (pi, modes (w, q), last) covers junk with weight pi and "
    "'not junk, and every coordinate of A dirty' with weight sum_k w_k (q_k + last)^|A|.  Check of "
    "sector s with survivor modes (w_i, q_i), partner-side modes (y_j, r_j) and carried rate l: "
    "inherited modes (w_i y_j, q_i r_j), (pi_B w_i, q_i), (pi_A y_j, r_j), united with l and the fresh "
    "rate; junk pi_A pi_B + C(n, wmin) sum_ij w_i y_j (q_i + l + r_j)^wmin; other sector junk "
    "pi_A + pi_B and modes (w_i w_j, q_i + q_j + fresh).")


def pad(modes):
    """Add the mode (1 - W, 0) when the weights sum to W < 1."""
    w = sum((m[0] for m in modes), Fraction(0))
    return list(modes) + ([(1 - w, Fraction(0))] if w < 1 else [])


def cap1(q):
    return min(ONE, Fraction(q))


def normalize(modes):
    """A list of total weight W > 1 becomes (w / W, W q): no weight for |A| >= 1 decreases."""
    modes = [(w, q) for w, q in modes if w > 0]
    total_w = sum((w for w, _ in modes), Fraction(0))
    if total_w <= 1:
        return modes
    return [(w / total_w, cap1(total_w * q)) for w, q in modes]


def united(a, b):
    """Modes of the union of two dirt processes, both lists padded."""
    return [(wa * wb, cap1(qa + qb)) for wa, qa in a for wb, qb in b if wa * wb > 0]


def shifted(modes, rate):
    """Modes of a dirt process united with dirt of one rate."""
    return [(w, cap1(q + rate)) for w, q in pad(modes)]


def moment(modes, k):
    return sum((w * q ** k for w, q in modes), Fraction(0))


def leaf(p, junk, local):
    """Contract of one sector of a raw preparation: `junk` junk-capable locations and at most
    `local` local locations at one coordinate, each of weight p."""
    p = Fraction(p)
    return fl.sector(cap1(junk * p), [(1, cap1(local * p))])


def total(sec):
    """Modes of the whole dirt of a sector, padded."""
    return shifted(sec["modes"], sec["last"])


def nontrivial(sec, n):
    """Weight of a cover of "the sector's pattern is not trivial"."""
    return cap1(sec["pi"] + n * moment(total(sec), 1))


def check(A, B, s, n, wmin, step, carry_last=False, independent=False, cap=MAX_MODES):
    """One accepted check of sector s; same arguments and result as filters.check.
    step["clean"] is 1 - (weight of the step locations that can cause rejection)."""
    if independent:
        raise ValueError("the local-stochastic calculus does not unite rates as independent")
    o = "Z" if s == "X" else "X"
    if not carry_last:
        fold_last = lambda C: {t: dict(pi=C[t]["pi"], modes=total(C[t]), last=Fraction(0)) for t in "XZ"}
        A, B = fold_last(A), fold_last(B)
    a, b = A[s], B[s]
    q_g, q_m = Fraction(step["fresh"][s]), Fraction(step["partner"])
    am = [(w, q) for w, q in a["modes"] if w > 0]                # dirt disjoint from the partner's side
    tb = shifted(b["modes"], b["last"] + q_m)                    # the partner's side
    inherited = [(wi * wj, qi * qj) for wi, qi in am for wj, qj in tb]
    inherited += [(b["pi"] * wi, qi) for wi, qi in am]
    inherited += [(a["pi"] * wj, qj) for wj, qj in tb]
    if wmin is None:
        p_w = Fraction(0)
    else:
        whole = united(shifted(am, a["last"]), tb)
        p_w = comb(n, wmin) * moment(whole, wmin) if wmin <= n else Fraction(0)
    pi_s = a["pi"] * b["pi"] + p_w
    inherited = normalize(fl.prune([(w, q) for w, q in inherited if w > 0], 4 * cap))
    if carry_last:
        modes_s, last_s = shifted(inherited, a["last"]), q_g
    else:
        modes_s, last_s = shifted(inherited, q_g), Fraction(0)
    ao, bo = A[o], B[o]
    q_o = Fraction(step["fresh"][o])
    both = united(pad(ao["modes"]), pad(bo["modes"]))
    if carry_last:
        modes_o, last_o = [(w, cap1(q + ao["last"] + bo["last"])) for w, q in both], q_o
    else:
        modes_o, last_o = [(w, cap1(q + q_o)) for w, q in both], Fraction(0)
    pi_o = ao["pi"] + bo["pi"]
    rho = nontrivial(A[s], n) + nontrivial(B[s], n) + (1 - Fraction(step["clean"]))
    alpha = fl.down(1 - rho)
    if alpha <= 0:
        return None, alpha

    def delivered(pi, modes, last):
        # weights of the first accepted attempt are those of one attempt over 1 - rho; the list is
        # brought back to total weight one by scaling every rate, the carried one included
        modes = fl.prune([(w / alpha, q) for w, q in modes if w > 0], cap)
        scale = max(ONE, sum((w for w, _ in modes), Fraction(0)))
        return dict(pi=cap1(up(pi / alpha)),
                    modes=[(up(w / scale), cap1(up(scale * q))) for w, q in modes],
                    last=cap1(up(scale * last)))

    return {s: delivered(pi_s, modes_s, last_s), o: delivered(pi_o, modes_o, last_o)}, alpha


def physical_tree(profile, p):
    """(leaves, nodes, root) for filters.fold from the profile of a physical tree: every
    location has weight p, and counts of locations add."""
    p = Fraction(p)
    kinds = []
    for sh in profile["shapes"]:
        kinds.append(dict(X=leaf(p, sh["X"]["junk"], max(sh["X"]["local"])),
                          Z=leaf(p, sh["Z"]["junk"], max(sh["Z"]["local"])), cost=sh["locations"]))
    leaves = {int(b): kinds[k] for b, k in profile["leaves"].items()}
    nodes = []
    for nd in profile["nodes"]:
        g, m, i = nd["gates"], nd["meas"], nd["idles"]
        fresh = cap1(max(a + b for a, b in zip(g, i)) * p)
        step = dict(fresh=dict(X=fresh, Z=fresh), partner=cap1(max(a + b for a, b in zip(g, m)) * p),
                    clean=1 - (sum(g) + sum(m)) * p, cost=sum(g) + sum(m) + sum(i))
        nodes.append(dict(sector=nd["sector"], A=tuple(nd["A"]), B=tuple(nd["B"]), step=step))
    return leaves, nodes, tuple(profile["root"])


def summary(sec, modes=4, order=1):
    """The sector as tails.py reads it under local-stochastic noise: junk weight and at most
    `modes` modes of the whole dirt."""
    out = fl.prune([(w, q) for w, q in total(sec) if w > 0 and q > 0], modes, order)
    return dict(junk=sec["pi"], modes=[(up(w), up(q)) for w, q in out] or [(Fraction(0), Fraction(0))])


def as_json(c):
    return {t: dict(junk=str(c[t]["junk"]), modes=[[str(w), str(q)] for w, q in c[t]["modes"]]) for t in "XZ"}


PREMISES_PHYSICAL = [
    "the filter rule for this tree (ftqc_gadgets.factories: the tree is recognised and its leaf "
    "locations classified; the rule is quoted and tested against propagation)",
    "a check accepts when both incoming patterns of its sector are trivial and no step location that can "
    "cause rejection is faulty, and a rejected check is repeated with fresh inputs",
    "Pr[all locations of S faulty] <= p^|S| for every set S of physical locations",
]
PREMISES_ENCODED = [
    "the filter rule for the tree, with decoding events of the level below as its locations",
    "the cone theorem for the raw preparation and for a check (ftqc_gadgets.cones), and the tail bounds over "
    "its cones (cone-bound, local-stochastic)",
    "a check accepts when both incoming patterns of its sector are trivial and no decoding event of the step "
    "fails, and a rejected check is repeated with fresh inputs",
    "the contracts of the ancillas of the level below, as covers on sources disjoint from the tree's locations "
    "and from one another",
    "the sources of different columns are disjoint",
    "Pr[all sources of S faulty] <= the product of their weights for every set S of sources",
]


def _report(n, res, modes, order):
    c = res["contract"]
    out = {s: summary(c[s], modes, order) for s in "XZ"}
    return dict(
        contract=as_json(out),
        junk={s: sci(out[s]["junk"]) for s in "XZ"},
        main_mode={s: [sci(max(out[s]["modes"], key=lambda m: m[0])[0], 6), sci(max(out[s]["modes"], key=lambda m: m[0])[1], 6)]
                   for s in "XZ"},
        not_trivial_upper={s: sci(nontrivial(c[s], n)) for s in "XZ"},
        accept_lower=[fl.decimal(fl.down(a, 20)) for a in res["accept_lower"]],
        expected_cost=fl.decimal(up(res["cost"], 12)),
        checks=[dict(sector=t["sector"], accept_lower=fl.decimal(fl.down(t["accept_lower"], 10)),
                     junk={s: sci(t["contract"][s]["pi"]) for s in "XZ"}) for t in res["trace"]],
    )


def _finish(out, n, res, modes, order):
    if res["vanished"]:
        out.update(status="no bound", reason="the weight of rejection reaches one at some check",
                   accept_lower=[fl.decimal(fl.down(a, 20)) if a > 0 else "0" for a in res["accept_lower"]])
        return out
    out.update(_report(n, res, modes, order))
    return out


def physical_result(spec):
    keys = ["kind", "level", "noise", "p", "factory" if "factory" in spec else "tree"]
    for opt in ("modes", "order", "max_modes", "leaf"):
        if opt in spec:
            keys.append(opt)
    fields(spec, *keys)
    p = probability(spec["p"])
    modes, order, cap = fl._options(spec)
    profile, n, wmin, seq = fl._tree_of(spec)
    if profile is not None:
        if "leaf" in spec:
            raise ValueError("the leaf classes come from the factory receipt")
        leaves, nodes, root = physical_tree(profile, p)
        shape = dict(leaves=len(leaves), checks=len(nodes), leaf_kinds=profile["shapes"])
    else:
        fields(spec["leaf"], "locations", "X", "Z")
        lf = {"cost": fl._positive_int(spec["leaf"]["locations"], "leaf locations")}
        for s in "XZ":
            fields(spec["leaf"][s], "junk", "local")
            j, l = spec["leaf"][s]["junk"], spec["leaf"][s]["local"]
            if type(j) is not int or type(l) is not int or j < 0 or l < 0:
                raise ValueError("leaf counts are nonnegative integers")
            lf[s] = leaf(p, j, l)
        nodes, root = fl.uniform(seq)
        for nd in nodes:
            nd["step"] = dict(fresh=dict(X=cap1(2 * p), Z=cap1(2 * p)), partner=cap1(2 * p), clean=1 - 2 * n * p, cost=3 * n)
        leaves = {i: lf for i in range(2 ** len(seq))}
        shape = dict(leaves=len(leaves), checks=len(nodes), sequence=seq, leaf=spec["leaf"],
                     step="per coordinate one gate, one measurement of the partner, one idle of the survivor")
    res = fl.fold(n, wmin, leaves, nodes, root, False, False, cap, check_fn=check)
    out = dict(status="conditional", rule="local-stochastic filter calculus, first accepted attempt",
               level="physical", noise="local-stochastic", p=str(p), coordinates=n, wmin=wmin, tree=shape,
               uniting="weights of alternatives add",
               rounding=f"weights upward, 1 - (weight of rejection) downward, {SIG} significant digits",
               premises=PREMISES_PHYSICAL, calculus=RULE)
    return _finish(out, n, res, modes, order)


def encoded_result(spec):
    from ftqc_gadgets import cones

    from . import tails

    keys = ["kind", "level", "noise", "p", "inner", "contracts", "raw", "factory" if "factory" in spec else "tree"]
    for opt in ("modes", "order", "max_modes"):
        if opt in spec:
            keys.append(opt)
    fields(spec, *keys)
    p = probability(spec["p"])
    modes, order, cap = fl._options(spec)
    profile, n, wmin, seq = fl._tree_of(spec)
    if profile is not None:
        if not profile["strict"]:
            raise ValueError("an encoded tree needs a factory certified with strict classes: every leaf "
                             "location before the last round counts as junk-capable")
        seq = fl._balanced(profile)
    fields(spec["raw"], "program", "last_layers")
    last_layers = fl._positive_int(spec["raw"]["last_layers"], "last_layers")

    def bounds(program):
        sub = dict(kind="cone-bound", noise="local-stochastic", p=spec["p"], inner=spec["inner"],
                   contracts=spec["contracts"], program=program)
        probs, _ = tails.event_bounds(sub)
        flow = cones.Flow(cones.Program.from_json(program))
        if flow.two:
            raise ValueError("the programs of a filter contract are single-level")
        return flow, probs

    raw = spec["raw"]["program"]
    flow, probs = bounds(raw)
    depth = len(raw["layers"])
    if len(raw["outputs"]) != n:
        raise ValueError(f"the raw preparation delivers {len(raw['outputs'])} blocks, the tree has {n} coordinates")
    ref = cones.Flow(cones.Program.from_json(dict(name="std", inputs=[], outputs=["b"], layers=[[["NEW", "std", "b"]]])))
    std_shape = {}
    for ev in ref.events:
        if ev["kind"] == "out":
            mask = sum(1 << i for i in ref.members_of(ev) if ref.sources[i].get("what", [""])[0] != "idle")
            std_shape[ev["sector"]] = ref.weight(mask)
    pi = dict(X=Fraction(0), Z=Fraction(0))
    last = {b: dict(X=Fraction(0), Z=Fraction(0)) for b in raw["outputs"]}
    seen, own = set(), {b: set() for b in raw["outputs"]}
    for ev, pr in zip(flow.events, probs):
        if ev["kind"] == "out":
            if flow.weight(flow.cone(ev)) != std_shape[ev["sector"]]:
                raise ValueError(f"block {ev['block']} does not leave the raw preparation as a block just corrected")
        elif ev["tick"] >= depth - last_layers:
            if ev["kind"] != "ec" or ev["block"] not in last:
                raise ValueError("the last layers of the raw preparation must be error correction of its outputs")
            last[ev["block"]][ev["sector"]] += pr
            own[ev["block"]] |= set(flow.members_of(ev))
        else:
            for s in ev["sector"]:
                pi[s] += pr
    for b in raw["outputs"]:
        if seen & own[b]:
            raise ValueError(f"the last round of block {b} shares sources with another block: its failures "
                             f"are not local to a coordinate")
        seen |= own[b]
        if not own[b]:
            raise ValueError(f"block {b} is not corrected in the last layers")
    leaf_c = {s: fl.sector(cap1(up(pi[s])), [(1, 0)], up(max(v[s] for v in last.values()))) for s in "XZ"}
    leaf_c["cost"] = Fraction(1)
    steps, step_report = {}, {}
    for s in sorted(set(seq)):
        fl_, pr = bounds(fl._step_program(s))
        fresh, partner, every = dict(X=Fraction(0), Z=Fraction(0)), Fraction(0), Fraction(0)
        for ev, x in zip(fl_.events, pr):
            if ev["kind"] == "out":
                continue
            every += x
            if ev["block"] == "A":
                for t in ev["sector"]:
                    fresh[t] += x
            else:
                partner += x
        clean = fl.down(1 - n * every)
        steps[s] = dict(fresh={t: up(v) for t, v in fresh.items()}, partner=up(partner), clean=clean, cost=Fraction(0))
        step_report[s] = dict(fresh={t: sci(v) for t, v in fresh.items()}, partner=sci(partner),
                              clean_lower=fl.decimal(fl.down(clean, 12)))
    shared = fl.shared_sources(seq)
    if not all(r["rest_disjoint"] for r in shared):
        raise ValueError("older dirt of a survivor shares sources with the partner's side of a check: "
                         "the weights may not be multiplied")
    carry = any(r["last_shared"] for r in shared)
    nodes, root = fl.uniform(seq)
    for nd in nodes:
        nd["step"] = steps[nd["sector"]]
    res = fl.fold(n, wmin, {i: leaf_c for i in range(2 ** len(seq))}, nodes, root, carry, False, cap, check_fn=check)
    out = dict(status="conditional", rule="local-stochastic filter calculus, first accepted attempt",
               level="encoded", noise="local-stochastic", p=str(p), coordinates=n, wmin=wmin, inner=spec["inner"],
               tree=dict(leaves=2 ** len(seq), checks=len(nodes), sequence=seq),
               raw=dict(program_sha256=cones.Program.from_json(raw).sha256(), decoding_events=len(flow.events),
                        junk={s: sci(leaf_c[s]["pi"]) for s in "XZ"}, last_round={s: sci(leaf_c[s]["last"]) for s in "XZ"},
                        last_round_local=True, blocks_end_as_just_corrected=True),
               step=step_report, shared_sources=shared, carry_last=carry, cost_unit="raw preparations",
               uniting="weights of alternatives add",
               rounding=f"weights upward, 1 - (weight of rejection) downward, {SIG} significant digits",
               premises=PREMISES_ENCODED, calculus=RULE)
    return _finish(out, n, res, modes, order)


def filter_contract_result(spec):
    level = spec.get("level")
    if level == "physical":
        return physical_result(spec)
    if level == "encoded":
        return encoded_result(spec)
    raise ValueError("level is physical or encoded")
