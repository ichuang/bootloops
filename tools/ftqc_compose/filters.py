"""Filter calculus: what an accepted purification tree delivers, with conditioning on acceptance.

A purification tree (ftqc_gadgets.factories) joins raw preparations by checks and keeps a
block only when every check on its way accepted.  The filter rule says where the error of an
accepted block can be, in terms of which locations were faulty.  This module turns that
statement into probabilities: the contract that the tail bounds (tails.py, study kind
`cone-bound`) take as a hypothesis for every ancilla kind.

contract    For one sector of an accepted block: a junk probability pi, modes [(w, q)] and a
            rate `last`.  Given acceptance, apart from an event of probability at most pi,
            the block's error is equivalent to one supported on a random set of dirty
            coordinates, and every increasing event E of that set satisfies
                P[not junk, E]  <=  sum over modes of  w * P[E under independent dirt of
                                                             rate q united with `last`].
            Weights are coefficients: conditioning can push their sum slightly above one.
leaf        A raw preparation with J junk-capable locations and at most L local locations at
            one coordinate: pi = 1 - (1 - p)^J, one mode (1, 1 - (1 - p)^L).
check       Survivor A against partner B in sector s (CHECK below): the survivor's old dirt
            survives only where the partner's side is dirty too, so rates multiply; junk
            needs junk on both sides, or an unlikely exceptional event; everything is divided
            by a lower bound on the probability of acceptance.
reset       Fresh dirt of the step is united with the inherited part without letting the
            leading coefficient grow from check to check (RESET below).

Two ways of uniting dirt at one coordinate are kept apart.  Where the sources are disjoint
sets of locations (a physical tree) rates combine as 1 - (1 - a)(1 - b).  Where they may
share sources (a tree of encoded blocks, whose "locations" are decoding events of the level
below) only a + b is used, and the dirt a block received in its latest step is carried
through the next check of its sector without being multiplied by the partner's rate
(`carry_last`): that dirt and the partner's measured syndrome can depend on the same inner
ancillas.  The study kind `filter-contract` decides which case applies from the cones of the
level below instead of taking it on trust.

Arithmetic is exact in rationals, rounded after every check to SIG significant digits: bounds
upward, the acceptance probability downward.
"""
from fractions import Fraction

from .stochastic import probability
from .supports import fields
from .tails import SIG, binomial_tail, sci, up

MAX_MODES = 16

CHECK = (
    "Check of sector s, survivor A and partner B, independent blocks with contracts (pi_A, modes "
    "(w_i, q_i)) and (pi_B, modes (w_j, q_j)); the step dirties a survivor coordinate with rate "
    "q_g and a coordinate on the partner's side with rate q_m.  By the filter rule the survivor "
    "is junk only if both are, or if neither is and the dirty sets of the two sides together "
    "reach the minimum weight w of the trivial patterns: pi = pi_A pi_B + P_W with P_W = "
    "sum_ij w_i w_j P[Bin(n, q_i + q_j + q_m) >= w].  Otherwise its dirt is the step's fresh "
    "dirt united with: the intersection of A's dirt and the partner side's (neither junk), "
    "modes (w_i w_j, q_i (q_j + q_m)); A's dirt (B junk), modes (pi_B w_i, q_i); the partner "
    "side's (A junk), modes (pi_A w_j, q_j + q_m).  In the other sector patterns add: junk "
    "pi_A + pi_B, modes (w_i w_j, q_i + q_j + q_g').  An event together with acceptance is no "
    "more likely than the event, so dividing every coefficient and junk probability by a lower "
    "bound alpha on the probability of acceptance gives the contract given acceptance; the "
    "check accepts when both incoming patterns are trivial and no step location that can "
    "cause rejection is faulty: alpha = (1 - nu_A)(1 - nu_B) * clean, with nu = pi + sum_k w_k "
    "(1 - (1 - q_k)^n).")

RESET = (
    "Dirt g united with c, where g has rate b at every coordinate and c has modes (w_k, r_k), "
    "r_k > 0, of total weight W.  For an increasing event let P_0 be its probability without c "
    "and P_k with c in mode k.  Either the event holds already without c, or it is an increasing "
    "event of c that fails at c = 0, so its probability is at most P_0 + sum_k w_k (P_k - P_0): "
    "modes (1 - W, b) and (w_k, b + r_k) when W <= 1.  When W > 1 the first coefficient is "
    "negative; thinning a vector of rate b + r_1 to rate b leaves it unchanged with probability "
    "at least lam = (b / (b + r_1))^n, so P_0 >= lam P_1 and (1 - W) P_0 <= (1 - W) lam P_1: the "
    "weight of the first mode becomes w_1 - (W - 1) lam.")


def down(x, sig=SIG):
    """Largest decimal with `sig` significant digits that is at most x (x >= 0)."""
    x = Fraction(x)
    if x <= 0:
        return Fraction(0)
    e = ((x.numerator.bit_length() - x.denominator.bit_length()) * 30103) // 100000
    scale = sig + 1 - e
    num = x * (Fraction(10) ** scale)
    return Fraction(num.numerator // num.denominator) / (Fraction(10) ** scale)


def decimal(x):
    """Exact decimal string of a rational whose denominator divides a power of ten."""
    x = Fraction(x)
    m = 0
    while (10 ** m) % x.denominator:
        m += 1
        if m > 4 * SIG:
            raise ValueError("not a terminating decimal")
    digits = str(x.numerator * 10 ** m // x.denominator).zfill(m + 1)
    return digits[:len(digits) - m] + ("." + digits[len(digits) - m:] if m else "")


def unite(independent, *rates):
    """Rate of the union of dirt of the given rates at one coordinate."""
    if not independent:
        return min(Fraction(1), sum((Fraction(r) for r in rates), Fraction(0)))
    clean = Fraction(1)
    for r in rates:
        clean *= 1 - Fraction(r)
    return 1 - clean


def some(q, n):
    """P[at least one of n coordinates dirty] at rate q."""
    return 1 - (1 - Fraction(q)) ** n


def prune(modes, cap=MAX_MODES, order=1):
    """Merge modes of equal rate, then, while there are more than `cap`, move the weight of
    one mode up to the next rate, choosing the move that adds least to an event of the given
    order (one that needs `order` dirty coordinates).  Upward: an increasing event is no less
    likely at a larger rate."""
    acc = {}
    for w, q in modes:
        if w > 0:
            acc[q] = acc.get(q, Fraction(0)) + w
    out = sorted(([w, q] for q, w in acc.items()), key=lambda m: m[1])
    while len(out) > max(1, cap):
        i = min(range(len(out) - 1), key=lambda i: out[i][0] * (out[i + 1][1] ** order - out[i][1] ** order))
        out[i:i + 2] = [[out[i][0] + out[i + 1][0], out[i + 1][1]]]
    return [(w, q) for w, q in out]


def sector(pi, modes, last=0):
    return dict(pi=Fraction(pi), modes=[(Fraction(w), Fraction(q)) for w, q in modes], last=Fraction(last))


def leaf(p, junk, local):
    """Contract of one sector of a raw preparation: `junk` junk-capable locations and at most
    `local` local locations at one coordinate, each faulty with probability at most p."""
    p = Fraction(p)
    return sector(1 - (1 - p) ** junk, [(1, 1 - (1 - p) ** local)])


def total(sec, independent=False):
    """Modes of the whole dirt of a sector: `last` united with every mode."""
    return [(w, unite(independent, q, sec["last"])) for w, q in sec["modes"]]


def nontrivial(sec, n, independent=False):
    """Upper bound of the probability that the sector's pattern is not trivial."""
    return min(Fraction(1), sec["pi"] + sum(w * some(q, n) for w, q in total(sec, independent)))


def reset(inherited, base, n, independent=False, cap=MAX_MODES):
    """Modes of (dirt of rate `base`) united with an inherited part of modes [(w, rate)].
    See RESET."""
    pos = sorted(((w, r) for w, r in inherited if r > 0 and w > 0), key=lambda m: m[1])
    if not pos:
        return [(Fraction(1), Fraction(base))]
    weight = sum(w for w, _ in pos)
    if weight <= 1:
        return prune([(1 - weight, base)] + [(w, unite(independent, base, r)) for w, r in pos], cap)
    if base == 0:
        return prune(pos, cap)
    # the small rates, moved up to the largest of them, take the correction
    small = [m for m in pos if m[1] * 1000 <= base] or pos[:1]
    rest = pos[len(small):]
    w1, r1 = sum(w for w, _ in small), small[-1][1]
    u1 = unite(independent, base, r1)
    lam = down((Fraction(base) / u1) ** n)
    w1 = max(Fraction(0), w1 - (weight - 1) * lam)
    return prune([(w1, u1)] + [(w, unite(independent, base, r)) for w, r in rest], cap)


def check(A, B, s, n, wmin, step, carry_last=False, independent=True, cap=MAX_MODES):
    """One accepted check of sector s.  A, B: {"X": sector, "Z": sector}, each given its own
    acceptance.  step: dict(fresh={X, Z}, partner=, clean=) with the per-coordinate rate of
    fresh dirt on the survivor in each sector, the rate of step faults on the partner's side,
    and a lower bound of the probability that no step location that can cause rejection is
    faulty.  wmin: minimum weight of the trivial patterns of sector s (None: there are none).  Returns (contract,
    lower bound of the probability of acceptance) or (None, bound <= 0).  See CHECK."""
    if carry_last and independent:
        raise ValueError("dirt that is carried because of shared sources cannot be united as independent")
    o = "Z" if s == "X" else "X"
    if not carry_last:
        fold = lambda C: {t: dict(pi=C[t]["pi"], modes=total(C[t], independent), last=Fraction(0)) for t in "XZ"}
        A, B = fold(A), fold(B)
    a, b = A[s], B[s]
    q_g, q_m = Fraction(step["fresh"][s]), Fraction(step["partner"])
    tb = [(w, unite(independent, q, b["last"], q_m)) for w, q in b["modes"]]       # the partner's side
    inherited = [(wi * wj, qi * qj) for wi, qi in a["modes"] for wj, qj in tb]
    inherited += [(b["pi"] * wi, qi) for wi, qi in a["modes"]]
    inherited += [(a["pi"] * wj, qj) for wj, qj in tb]
    p_w = Fraction(0) if wmin is None else sum(
        wi * wj * binomial_tail(n, wmin - 1, up(unite(independent, qi, a["last"], qj)))
        for wi, qi in a["modes"] for wj, qj in tb)
    pi_s = a["pi"] * b["pi"] + p_w
    if carry_last:
        modes_s, last_s = reset(prune(inherited, 4 * cap), a["last"], n, independent, cap), q_g
    else:
        modes_s, last_s = reset(prune(inherited, 4 * cap), q_g, n, independent, cap), Fraction(0)
    ao, bo = A[o], B[o]
    q_o = Fraction(step["fresh"][o])
    if carry_last:
        modes_o = [(wi * wj, unite(independent, qi, qj, ao["last"], bo["last"]))
                   for wi, qi in ao["modes"] for wj, qj in bo["modes"]]
        last_o = q_o
    else:
        modes_o = [(wi * wj, unite(independent, qi, qj, q_o)) for wi, qi in ao["modes"] for wj, qj in bo["modes"]]
        last_o = Fraction(0)
    pi_o = ao["pi"] + bo["pi"]
    alpha = down((1 - nontrivial(A[s], n, independent)) * (1 - nontrivial(B[s], n, independent))
                 * Fraction(step["clean"]))
    if alpha <= 0:
        return None, alpha

    def given(pi, modes, last):
        return dict(pi=min(Fraction(1), up(pi / alpha)),
                    modes=[(up(w), up(q)) for w, q in prune([(w / alpha, q) for w, q in modes], cap)],
                    last=up(last))

    return {s: given(pi_s, modes_s, last_s), o: given(pi_o, modes_o, last_o)}, alpha


def fold(n, wmin, leaves, nodes, root, carry_last=False, independent=True, cap=MAX_MODES, check_fn=None):
    """Contract at the root of a tree.  leaves: {id: dict(X=sector, Z=sector, cost=)};
    nodes: [dict(sector=, A=("leaf"|"node", id), B=..., step=dict(fresh, partner, clean, cost))]
    in an order where children come first.  Identical subtrees are computed once.  Returns
    dict(contract, accept_lower=[one per distinct check, leaves first], cost, vanished)."""
    memo, shape, accs, trace = {}, {}, [], []

    def canon(x):
        return repr(sorted(x.items())) if isinstance(x, dict) else repr(x)

    for b, lf in leaves.items():
        key = ("leaf", canon({t: (lf[t]["pi"], tuple(lf[t]["modes"]), lf[t]["last"]) for t in "XZ"}), lf.get("cost", 0))
        shape[("leaf", b)] = key
        memo.setdefault(key, ({t: lf[t] for t in "XZ"}, Fraction(lf.get("cost", 0))))
    for i, nd in enumerate(nodes):
        ka, kb = shape[tuple(nd["A"])], shape[tuple(nd["B"])]
        st = nd["step"]
        key = ("node", nd["sector"], ka, kb, canon(st["fresh"]), str(st["partner"]), str(st["clean"]), str(st.get("cost", 0)))
        shape[("node", i)] = key
        if key in memo:
            continue
        (ca, ea), (cb, eb) = memo[ka], memo[kb]
        if ca is None or cb is None:
            memo[key] = (None, None)
            continue
        c, alpha = (check_fn or check)(ca, cb, nd["sector"], n, wmin[nd["sector"]], st, carry_last, independent, cap)
        accs.append(alpha)
        if c is None:
            memo[key] = (None, None)
            continue
        memo[key] = (c, up((ea + eb + Fraction(st.get("cost", 0))) / alpha))
        trace.append(dict(sector=nd["sector"], accept_lower=alpha, contract=c))
    c, cost = memo[shape[tuple(root)]]
    return dict(contract=c, accept_lower=accs, cost=cost, vanished=c is None, trace=trace)


def uniform(sequence):
    """Nodes and root of the balanced tree whose level l checks sector sequence[l]: ids of
    the leaves are 0 .. 2^len - 1, the survivor of a check is its first child."""
    cur = [("leaf", i) for i in range(2 ** len(sequence))]
    nodes = []
    for s in sequence:
        nxt = []
        for i in range(0, len(cur), 2):
            nodes.append(dict(sector=s, A=cur[i], B=cur[i + 1]))
            nxt.append(("node", len(nodes) - 1))
        cur = nxt
    return nodes, cur[0]


def physical_tree(profile, p):
    """(leaves, nodes, root) for `fold`, from the profile of a physical purification tree
    (ftqc_gadgets.factories.Tree.profile) when every location is faulty independently with
    probability at most p.  A coordinate takes the largest count over coordinates."""
    p = Fraction(p)
    n = profile["coordinates"]
    keep = 1 - p
    kinds = []
    for sh in profile["shapes"]:
        kinds.append(dict(X=leaf(p, sh["X"]["junk"], max(sh["X"]["local"])),
                          Z=leaf(p, sh["Z"]["junk"], max(sh["Z"]["local"])), cost=sh["locations"]))
    leaves = {("leaf", int(b)): kinds[k] for b, k in profile["leaves"].items()}
    nodes = []
    for nd in profile["nodes"]:
        g, m, i = nd["gates"], nd["meas"], nd["idles"]
        fresh = 1 - keep ** max(a + b for a, b in zip(g, i))
        step = dict(fresh=dict(X=fresh, Z=fresh), partner=1 - keep ** max(a + b for a, b in zip(g, m)),
                    clean=keep ** (sum(g) + sum(m)), cost=sum(g) + sum(m) + sum(i))
        nodes.append(dict(sector=nd["sector"], A=tuple(nd["A"]), B=tuple(nd["B"]), step=step))
    return {b[1]: v for b, v in leaves.items()}, nodes, tuple(profile["root"])


def summary(sec, modes=3, independent=False, order=1):
    """The sector as tails.py reads it: junk probability and at most `modes` modes of the
    whole dirt, merged with least loss for events of the given order."""
    out = prune(total(sec, independent), modes, order)
    return dict(pi=sec["pi"], modes=[(up(w), up(q)) for w, q in out])


def as_json(c):
    return {t: dict(pi=str(c[t]["pi"]), modes=[[str(w), str(q)] for w, q in c[t]["modes"]]) for t in "XZ"}


def event_bound(sec, n, k):
    """What a summary gives for "more than k dirty coordinates, or junk"."""
    return min(Fraction(1), sec["pi"] + sum(w * binomial_tail(n, k, q) for w, q in sec["modes"]))


# ---------------------------------------------------------------------------- the study
PREMISES_PHYSICAL = [
    "the filter rule for this tree (ftqc_gadgets.factories: the tree is recognised and its leaf "
    "locations classified; the rule is quoted and tested against propagation)",
    "faults at different locations are independent, each with probability at most p",
    "the filter calculus (CHECK and RESET, quoted; tested against exact laws of small trees, not proved here)",
]
PREMISES_ENCODED = [
    "the filter rule for the tree, with decoding events of the level below as its locations",
    "the cone theorem for the raw preparation and for a check (ftqc_gadgets.cones), and the tail bounds over "
    "its cones (cone-bound)",
    "the declared contracts of the ancillas of the level below",
    "faults at different physical locations are independent, each with probability at most p",
    "every dirt indicator is an increasing function of independent sources (used for the bound on acceptance)",
    "the filter calculus with shared sources (CHECK and RESET, quoted; tested against exact laws of small "
    "models with shared sources, not proved here)",
]


def _positive_int(x, what):
    if type(x) is not int or x < 1:
        raise ValueError(f"{what} must be a positive integer")
    return x


def _balanced(profile):
    """The check sequence of a profile whose tree is balanced with one sector per level, or
    a ValueError."""
    depth, sector_at = {}, {}

    def level(sub):
        sub = tuple(sub)
        if sub[0] == "leaf":
            return 0
        if sub not in depth:
            nd = profile["nodes"][sub[1]]
            a, b = level(nd["A"]), level(nd["B"])
            if a != b:
                raise ValueError("the tree is not balanced")
            depth[sub] = a + 1
            if sector_at.setdefault(a, nd["sector"]) != nd["sector"]:
                raise ValueError("checks of one level have different sectors")
        return depth[sub]

    top = level(profile["root"])
    if len(profile["leaves"]) != 2 ** top:
        raise ValueError("the tree is not balanced")
    return "".join(sector_at[i] for i in range(top))


def _tree_of(spec):
    """(profile or None, n, wmin, sequence) from a factory receipt or from explicit fields."""
    if "factory" in spec:
        from ftqc_gadgets import certify, factories
        from ftqc_gadgets.gadget import GadgetError

        try:
            profile = factories.profile(spec["factory"])
        except (certify.Refusal, GadgetError, KeyError, TypeError) as e:
            raise ValueError(f"factory receipt: {e}") from None
        return profile, profile["coordinates"], profile["wmin"], None
    fields(spec["tree"], "n", "wmin", "sequence")
    t = spec["tree"]
    n = _positive_int(t["n"], "tree n")
    seq = t["sequence"]
    if not isinstance(seq, str) or not seq or set(seq) - set("XZ"):
        raise ValueError("a check sequence is a nonempty word in X and Z")
    wmin = {}
    for s in "XZ":
        w = t["wmin"].get(s)
        wmin[s] = None if w is None else _positive_int(w, "wmin")
    return None, n, wmin, seq


def _report(n, res, modes, order, independent):
    c = res["contract"]
    out = {s: summary(c[s], modes, independent, order) for s in "XZ"}
    return dict(
        contract=as_json(out),
        junk={s: sci(out[s]["pi"]) for s in "XZ"},
        main_mode={s: [sci(max(out[s]["modes"], key=lambda m: m[0])[0], 6), sci(max(out[s]["modes"], key=lambda m: m[0])[1], 6)]
                   for s in "XZ"},
        not_trivial_upper={s: sci(nontrivial(c[s], n, independent)) for s in "XZ"},
        accept_lower=[decimal(down(a, 20)) for a in res["accept_lower"]],
        expected_cost=decimal(up(res["cost"], 12)),
        checks=[dict(sector=t["sector"], accept_lower=decimal(down(t["accept_lower"], 10)),
                     junk={s: sci(t["contract"][s]["pi"]) for s in "XZ"}) for t in res["trace"]],
    )


def _options(spec):
    modes = spec.get("modes", 4)
    order = spec.get("order", 1)
    cap = spec.get("max_modes", MAX_MODES)
    for v, what in ((modes, "modes"), (order, "order"), (cap, "max_modes")):
        _positive_int(v, what)
    return modes, order, cap


def physical_result(spec):
    keys = ["kind", "level", "p", "factory" if "factory" in spec else "tree"]
    for opt in ("modes", "order", "max_modes", "leaf"):
        if opt in spec:
            keys.append(opt)
    fields(spec, *keys)
    p = probability(spec["p"])
    modes, order, cap = _options(spec)
    profile, n, wmin, seq = _tree_of(spec)
    if profile is not None:
        if "leaf" in spec:
            raise ValueError("the leaf classes come from the factory receipt")
        leaves, nodes, root = physical_tree(profile, p)
        shape = dict(leaves=len(leaves), checks=len(nodes), leaf_kinds=profile["shapes"])
    else:
        fields(spec["leaf"], "locations", "X", "Z")
        lf = {"cost": _positive_int(spec["leaf"]["locations"], "leaf locations")}
        for s in "XZ":
            fields(spec["leaf"][s], "junk", "local")
            j, l = spec["leaf"][s]["junk"], spec["leaf"][s]["local"]
            if type(j) is not int or type(l) is not int or j < 0 or l < 0:
                raise ValueError("leaf counts are nonnegative integers")
            lf[s] = leaf(p, j, l)
        nodes, root = uniform(seq)
        two = 1 - (1 - p) ** 2
        for nd in nodes:
            nd["step"] = dict(fresh=dict(X=two, Z=two), partner=two, clean=(1 - p) ** (2 * n), cost=3 * n)
        leaves = {i: lf for i in range(2 ** len(seq))}
        shape = dict(leaves=len(leaves), checks=len(nodes), sequence=seq, leaf=spec["leaf"],
                     step="per coordinate one gate, one measurement of the partner, one idle of the survivor")
    res = fold(n, wmin, leaves, nodes, root, False, True, cap)
    out = dict(status="conditional", rule="filter calculus with conditioning on acceptance", level="physical",
               p=str(p), coordinates=n, wmin=wmin, tree=shape, uniting="independent sources: 1 - (1 - a)(1 - b)",
               rounding=f"bounds upward, acceptance downward, {SIG} significant digits",
               premises=PREMISES_PHYSICAL, check=CHECK, reset=RESET)
    if res["vanished"]:
        out.update(status="no bound", reason="the lower bound on acceptance is not positive at some check",
                   accept_lower=[decimal(down(a, 20)) if a > 0 else "0" for a in res["accept_lower"]])
        return out
    out.update(_report(n, res, modes, order, True))
    return out


def _std_sources(flow, block):
    return {i for i, s in enumerate(flow.sources)
            if s["blocks"] == [block] and (s.get("origin") == "std" or (s["kind"] == "class" and s["what"][0] == "std"))}


def _step_program(s):
    gate = ["CX", "A", "B"] if s == "X" else ["CX", "B", "A"]
    return dict(name=f"check {s}", inputs=[], outputs=["A"],
                layers=[[["NEW", "std", "A", "B"], gate], [["MZ" if s == "X" else "MX", "B"], ["ECa", "A"]],
                        [["ECb", "A"]]])


def _column_program(seq, correct=True):
    """One column of the whole tree, leaves as blocks just after a correction.  correct=False
    leaves out the correction of the survivor (for testing the analysis)."""
    alive = [f"L{i}" for i in range(2 ** len(seq))]
    layers, plan = [], []
    for lvl, s in enumerate(seq):
        pairs = [(alive[i], alive[i + 1]) for i in range(0, len(alive), 2)]
        first = [["NEW", "std"] + alive] if lvl == 0 else []
        layers.append(first + [["CX", a, b] if s == "X" else ["CX", b, a] for a, b in pairs])
        layers.append([["MZ" if s == "X" else "MX", b] for _, b in pairs] + ([["ECa", a] for a, _ in pairs] if correct else []))
        layers.append([["ECb", a] for a, _ in pairs] if correct else [])
        plan.append(pairs)
        alive = [a for a, _ in pairs]
    return dict(name="column", inputs=[], outputs=alive, layers=layers), plan


def shared_sources(seq, correct=True):
    """Which dirt of a survivor shares sources with the partner's side of its next check, read
    from the cones of one column of the tree.  For every check: `rest` is the dirt the
    survivor's subtree produced before its latest step, `last` the dirt of that step (for a
    leaf: of the correction that ends the raw preparation).  Returns one row per level with
    rest_disjoint (must be true for rates to multiply) and last_shared."""
    from ftqc_gadgets import cones

    prog, plan = _column_program(seq, correct)
    flow = cones.Flow(cones.Program.from_json(prog))
    events = [(ev, set(flow.members_of(ev))) for ev in flow.events if ev["kind"] in ("ec", "meas")]
    below = {f"L{i}": {f"L{i}"} for i in range(2 ** len(seq))}       # blocks of the subtree a block survives
    rows = []
    for lvl, pairs in enumerate(plan):
        s, t0 = seq[lvl], 3 * lvl
        row = dict(level=lvl, sector=s, rest_disjoint=True, last_shared=False)
        for a, b in pairs:
            def cone_of(blocks, pred):
                m = set()
                for ev, mem in events:
                    if ev["block"] in blocks and ev["sector"] == s and pred(ev):
                        m |= mem
                return m
            leaves_a = {x for x in below[a]}
            std_a = set().union(*[_std_sources(flow, x) for x in leaves_a])
            std_b = set().union(*[_std_sources(flow, x) for x in below[b]])
            if lvl == 0:
                rest, last = set(), std_a
            else:
                latest = lambda ev: ev["block"] == a and ev["kind"] == "ec" and ev["tick"] >= t0 - 3
                rest = cone_of(below[a], lambda ev: ev["tick"] < t0 and not latest(ev)) | std_a
                last = cone_of({a}, lambda ev: ev["tick"] < t0 and latest(ev))
            partner = cone_of(below[b], lambda ev: ev["tick"] < t0) | std_b
            partner |= cone_of({b}, lambda ev: ev["kind"] == "meas" and t0 <= ev["tick"] < t0 + 3)
            if rest & partner:
                row["rest_disjoint"] = False
            if last & partner:
                row["last_shared"] = True
        for a, b in pairs:
            below[a] = below[a] | below[b]
        rows.append(row)
    return rows


def encoded_result(spec):
    from ftqc_gadgets import cones

    from . import tails

    keys = ["kind", "level", "p", "inner", "contracts", "raw", "factory" if "factory" in spec else "tree"]
    for opt in ("modes", "order", "max_modes", "costs", "inner_costs"):
        if opt in spec:
            keys.append(opt)
    fields(spec, *keys)
    if "costs" in spec and "inner_costs" in spec:
        raise ValueError("give costs or inner_costs, not both")
    p = probability(spec["p"])
    modes, order, cap = _options(spec)
    profile, n, wmin, seq = _tree_of(spec)
    if profile is not None:
        if not profile["strict"]:
            raise ValueError("an encoded tree needs a factory certified with strict classes: every leaf "
                             "location before the last round counts as junk-capable")
        seq = _balanced(profile)
    fields(spec["raw"], "program", "last_layers")
    last_layers = _positive_int(spec["raw"]["last_layers"], "last_layers")

    def bounds(program):
        sub = dict(kind="cone-bound", noise="independent", p=spec["p"], inner=spec["inner"],
                   contracts=spec["contracts"], program=program)
        probs, _ = tails.event_bounds(sub)
        flow = cones.Flow(cones.Program.from_json(program))
        if flow.two:
            raise ValueError("the programs of a filter contract are single-level")
        return flow, probs

    # the raw preparation: early failures are junk, the last round is local to a block
    raw = spec["raw"]["program"]
    flow, probs = bounds(raw)
    depth = len(raw["layers"])
    if len(raw["outputs"]) != n:
        raise ValueError(f"the raw preparation delivers {len(raw['outputs'])} blocks, the tree has {n} coordinates")
    ref = cones.Flow(cones.Program.from_json(dict(name="std", inputs=[], outputs=["b"], layers=[[["NEW", "std", "b"]]])))
    std_shape = {}                                             # what a block just corrected carries
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
    leaf_c = {s: sector(min(Fraction(1), up(pi[s])), [(1, 0)], up(max(v[s] for v in last.values()))) for s in "XZ"}
    costs = spec.get("costs", {})
    if costs:
        fields(costs, "raw", "X", "Z")
    priced = None
    if "inner_costs" in spec:
        # expected physical locations: the raw preparation and one check, priced by the ledger
        from . import volumes

        fields(spec["inner_costs"], "plus", "zero")
        inner_costs = {k: volumes.cost_of(v, f"inner cost of {k}") for k, v in spec["inner_costs"].items()}
        n1 = spec["inner"]["n"]
        raw_counts, _, raw_total, _ = volumes.program_volume(raw, n1, None, inner_costs)
        costs = dict(raw=str(raw_total))
        priced = dict(raw=dict(counts=raw_counts, locations=decimal(up(raw_total, 15))))
        for s in sorted(set(seq)):
            st_counts, _, st_total, _ = volumes.program_volume(_step_program(s), n1, None, inner_costs)
            costs[s] = str(n * st_total)
            priced[s] = dict(counts=st_counts, locations_per_column=decimal(up(st_total, 15)))
    leaf_c["cost"] = Fraction(str(costs["raw"])) if costs else Fraction(1)
    # one check of each sector
    steps, step_report = {}, {}
    for s in sorted(set(seq)):
        fl, pr = bounds(_step_program(s))
        fresh, partner, every = dict(X=Fraction(0), Z=Fraction(0)), Fraction(0), Fraction(0)
        for ev, x in zip(fl.events, pr):
            if ev["kind"] == "out":
                continue
            every += x
            if ev["block"] == "A":
                for t in ev["sector"]:
                    fresh[t] += x
            else:
                partner += x
        clean = down(1 - n * every)
        steps[s] = dict(fresh={t: up(v) for t, v in fresh.items()}, partner=up(partner), clean=clean,
                        cost=Fraction(str(costs[s])) if costs else Fraction(0))
        step_report[s] = dict(fresh={t: sci(v) for t, v in fresh.items()}, partner=sci(partner), clean_lower=decimal(down(clean, 12)))
    shared = shared_sources(seq)
    if not all(r["rest_disjoint"] for r in shared):
        raise ValueError("older dirt of a survivor shares sources with the partner's side of a check: "
                         "the rates may not be multiplied")
    carry = any(r["last_shared"] for r in shared)
    nodes, root = uniform(seq)
    for nd in nodes:
        nd["step"] = steps[nd["sector"]]
    res = fold(n, wmin, {i: leaf_c for i in range(2 ** len(seq))}, nodes, root, carry, False, cap)
    out = dict(status="conditional", rule="filter calculus with conditioning on acceptance", level="encoded",
               p=str(p), coordinates=n, wmin=wmin, inner=spec["inner"],
               tree=dict(leaves=2 ** len(seq), checks=len(nodes), sequence=seq),
               raw=dict(program_sha256=cones.Program.from_json(raw).sha256(), decoding_events=len(flow.events),
                        junk={s: sci(leaf_c[s]["pi"]) for s in "XZ"}, last_round={s: sci(leaf_c[s]["last"]) for s in "XZ"},
                        last_round_local=True, blocks_end_as_just_corrected=True),
               step=step_report, shared_sources=shared, carry_last=carry,
               cost_unit=("physical locations" if priced else "as declared") if costs else "raw preparations",
               uniting="sources may be shared: a + b",
               rounding=f"bounds upward, acceptance downward, {SIG} significant digits",
               premises=PREMISES_ENCODED, check=CHECK, reset=RESET)
    if priced:                       # only when costs were asked for: earlier receipts recompute unchanged
        out["priced"] = priced
    if res["vanished"]:
        out.update(status="no bound", reason="the lower bound on acceptance is not positive at some check",
                   accept_lower=[decimal(down(a, 20)) if a > 0 else "0" for a in res["accept_lower"]])
        return out
    out.update(_report(n, res, modes, order, False))
    return out


def filter_contract_result(spec):
    """Study kind `filter-contract`: the contract of an accepted purification tree."""
    level = spec.get("level")
    if "noise" in spec:
        if spec["noise"] != "local-stochastic":
            raise ValueError("noise, when given, is local-stochastic; the default is independent faults")
        from . import filters_ls

        return filters_ls.filter_contract_result(spec)
    if level == "physical":
        return physical_result(spec)
    if level == "encoded":
        return encoded_result(spec)
    raise ValueError("level is physical or encoded")


__all__ = ["CHECK", "RESET", "as_json", "check", "decimal", "down", "event_bound", "fold", "leaf", "nontrivial", "physical_tree", "prune",
           "reset", "sci", "sector", "shared_sources", "summary", "total", "uniform", "unite", "filter_contract_result"]
