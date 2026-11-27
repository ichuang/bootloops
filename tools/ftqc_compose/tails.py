"""Tail bounds over fault cones: how likely a decoding event is to see too many dirty coordinates.

A decoding event of a transversal program is wrong only if more than t of its N coordinates
are dirty (ftqc_gadgets.cones).  A coordinate is dirty when one of the L locations the cone
holds at that coordinate is faulty, or an ancilla in the cone carries an error there.

    independent          faults at different locations are independent, each with
                         probability at most p.  A coordinate is dirty with probability
                         q = 1 - (1 - p)^L (1 - q_1) ... (1 - q_k), and the event fails with
                         the binomial tail P[Bin(N, q) > t], evaluated exactly.
    local-stochastic     only Pr[all of S faulty] <= p^|S| is assumed, for every set S of
                         locations.  More than t dirty coordinates needs one faulty location
                         in each of t + 1 coordinates: at most C(N, t+1) (L p)^(t+1).

Ancillas enter through declared contracts.  Under `independent`, a contract of an ancilla
kind and sector is a junk probability pi and modes [(w, q)]: apart from an event of
probability pi, charged where the ancilla is delivered, every increasing event of its dirt
is at most the sum over modes of w times its probability under independent per-coordinate
dirt of rate q.  Several ancillas in one cone multiply out.  Under `local-stochastic` a
contract is a number of location equivalents per coordinate, added to L.  Either way the
contract is a hypothesis of the result, not something this module derives.

With an outer level a column is dirty when one of the inner events of the outer cone fails
(a union bound inside the column), columns are independent, and the outer event fails
with the same kind of tail over the N2 columns.

Arithmetic is exact in rationals except where the result says otherwise: values are rounded
upward to SIG significant digits to keep denominators bounded.  Every rounding is upward
and every function used is increasing in what was rounded.
"""
import hashlib
import json
from fractions import Fraction
from math import comb

from .stochastic import probability
from .supports import fields

SIG = 50
MAX_TERMS = 400
_POWERS = {}


def up(x, sig=SIG):
    """Smallest decimal with `sig` significant digits that is at least x (x >= 0)."""
    x = Fraction(x)
    if x <= 0:
        return Fraction(0)
    # decimal exponent within one of the true one, from bit lengths (no string conversion)
    e = ((x.numerator.bit_length() - x.denominator.bit_length()) * 30103) // 100000
    scale = sig + 1 - e
    num = x * (Fraction(10) ** scale)
    n = -((-num.numerator) // num.denominator)                 # ceiling
    return Fraction(n) / (Fraction(10) ** scale)


def binomial_tail(n, t, q):
    """P[Bin(n, q) > t], exactly."""
    q = Fraction(q)
    if not 0 <= q <= 1:
        raise ValueError("a per-coordinate rate outside [0, 1]")
    if t < 0:
        return Fraction(1)
    if t >= n:
        return Fraction(0)
    a, d = q.numerator, q.denominator                          # integers throughout, one reduction at the end
    b = d - a
    if t + 1 <= n - t:                                         # the shorter of the two sums
        low = sum(comb(n, k) * a ** k * b ** (n - k) for k in range(0, t + 1))
        return Fraction(d ** n - low, d ** n)
    return Fraction(sum(comb(n, k) * a ** k * b ** (n - k) for k in range(t + 1, n + 1)), d ** n)


def threshold_enumerator(n, t, weights):
    """Number of minimal bad supports of "more than t dirty coordinates" when coordinate c
    has weights[c] locations: one location in each of t + 1 coordinates."""
    if len(weights) != n:
        raise ValueError("one weight per coordinate")
    e = [1] + [0] * (t + 1)                                    # elementary symmetric polynomial
    for w in weights:
        for k in range(t + 1, 0, -1):
            e[k] += e[k - 1] * w
    return e[t + 1] if t + 1 <= n else 0


def local_stochastic_tail(n, t, weight, p):
    if t >= n:
        return Fraction(0)
    return min(Fraction(1), comb(n, t + 1) * (weight * Fraction(p)) ** (t + 1))


def combine(base, rates):
    """Probability that a coordinate is dirty: independent sources of the given rates."""
    clean = 1 - Fraction(base)
    for q in rates:
        clean *= 1 - Fraction(q)
    return 1 - clean


FIX = 10 ** 200                 # fixed point of the mixture arithmetic


def _floor(x):
    x = Fraction(x)
    return x.numerator * FIX // x.denominator


def _ceil(x):
    x = Fraction(x)
    return -((-x.numerator * FIX) // x.denominator)


def _tail_fixed(n, t, clean):
    """ceil(FIX * P[Bin(n, q) > t]) for q = 1 - clean / FIX, in integers."""
    if t < 0:
        return FIX
    if t >= n:
        return 0
    a, b = FIX - clean, clean
    d = FIX ** n
    if t + 1 <= n - t:
        num = d - sum(comb(n, k) * a ** k * b ** (n - k) for k in range(0, t + 1))
    else:
        num = sum(comb(n, k) * a ** k * b ** (n - k) for k in range(t + 1, n + 1))
    return -((-num * FIX) // d)


def mixture_tail(n, t, base, atoms):
    """Bound for an event whose coordinates are dirty with rate `base` from locations and
    whose cone holds ancillas with mode lists atoms = [[(w, q), ...], ...].  Returns
    (bound, terms, coarsened).  Internally a mixture is {clean probability: weight} in fixed
    point with 200 digits, clean probabilities rounded downward and weights upward."""
    groups = {}
    for modes in atoms:
        key = tuple((_ceil(w), _floor(1 - Fraction(q))) for w, q in modes)
        groups[key] = groups.get(key, 0) + 1
    terms, coarse = {_floor(1 - Fraction(base)): FIX}, False
    for modes, m in sorted(groups.items()):
        # m identical ancillas: the m-fold product of one mixture, by dynamic programming;
        # the powers of one mixture are kept, since many cones hold the same kind
        chain = _POWERS.setdefault((modes, MAX_TERMS, t), [({FIX: FIX}, False)])
        while len(chain) <= m:
            nxt = {}
            for c0, w0 in chain[-1][0].items():
                for w, c1 in modes:
                    c = c0 * c1 // FIX
                    nxt[c] = nxt.get(c, 0) - ((-w0 * w) // FIX)
            power, did = _coarsen(nxt, t + 1)
            chain.append((power, did or chain[-1][1]))
        power, did = chain[m]
        coarse = coarse or did
        nxt = {}
        for c0, w0 in terms.items():
            for c1, w1 in power.items():
                c = c0 * c1 // FIX
                nxt[c] = nxt.get(c, 0) - ((-w0 * w1) // FIX)
        terms, did = _coarsen(nxt, t + 1)
        coarse = coarse or did
    total = 0
    for clean, w in terms.items():
        total += -((-w * _tail_fixed(n, t, max(0, clean))) // FIX)
    return min(Fraction(1), up(Fraction(total, FIX))), len(terms), coarse


def _coarsen(terms, k=1):
    """When there are many terms, keep the ones that matter most for a tail of order k
    (weight times dirt^k) and merge the rest in runs of neighbouring rates, each run taking
    its dirtiest rate.  Upward throughout."""
    if len(terms) <= MAX_TERMS:
        return terms, False
    keep = MAX_TERMS * 3 // 4
    ranked = sorted(terms, key=lambda c: terms[c] * (FIX - c) ** k, reverse=True)
    merged = {c: terms[c] for c in ranked[:keep]}
    rest = sorted(ranked[keep:])                                # smallest clean = dirtiest first
    step = len(rest) // (MAX_TERMS - keep) + 1
    for i in range(0, len(rest), step):
        chunk = rest[i:i + step]
        merged[chunk[0]] = merged.get(chunk[0], 0) + sum(terms[c] for c in chunk)
    return merged, True


def sci(x, digits=4):
    """Decimal string not below x, with `digits` significant digits."""
    x = Fraction(x)
    if x <= 0:
        return "0"
    e = ((x.numerator.bit_length() - x.denominator.bit_length()) * 30103) // 100000
    while x >= Fraction(10) ** (e + 1):
        e += 1
    while x < Fraction(10) ** e:
        e -= 1
    m = x / Fraction(10) ** e * 10 ** (digits - 1)
    n = -((-m.numerator) // m.denominator)                     # ceiling
    if n == 10 ** digits:
        n, e = 10 ** (digits - 1), e + 1
    s = str(n)
    return (s[0] + "." + s[1:] + f"e{e}") if digits > 1 else f"{s}e{e}"


def _contract(c, where, noise):
    if noise == "local-stochastic":
        fields(c, "equivalents")
        if type(c["equivalents"]) is not int or c["equivalents"] < 0:
            raise ValueError(f"{where}: location equivalents must be a nonnegative integer")
        return dict(equivalents=c["equivalents"])
    fields(c, "pi", "modes")
    modes = c["modes"]
    if not isinstance(modes, list) or not modes:
        raise ValueError(f"{where}: a contract needs at least one mode")
    out = []
    for m in modes:
        if not isinstance(m, list) or len(m) != 2:
            raise ValueError(f"{where}: a mode is [weight, rate]")
        w = Fraction(str(m[0])) if type(m[0]) in (str, int) else None
        if w is None or w < 0:
            raise ValueError(f"{where}: a mode weight must be an exact nonnegative number")
        out.append((w, probability(m[1])))
    return dict(pi=probability(c["pi"]), modes=out)


def _compute(spec):
    from ftqc_gadgets import cones

    keys = ["kind", "noise", "p", "inner", "contracts", "program"]
    for opt in ("outer", "stated"):
        if opt in spec:
            keys.append(opt)
    fields(spec, *keys)
    noise = spec["noise"]
    if noise not in ("independent", "local-stochastic"):
        raise ValueError("noise is independent or local-stochastic")
    p = probability(spec["p"])
    levels = {}
    for name in ("inner", "outer"):
        if name in spec:
            fields(spec[name], "n", "t")
            n, t = spec[name]["n"], spec[name]["t"]
            if type(n) is not int or type(t) is not int or not 0 <= t < n:
                raise ValueError(f"{name}: need integers 0 <= t < n")
            levels[name] = (n, t)
    program = cones.Program.from_json(spec["program"])
    # with an outer code every block is an outer block, whether or not the program delivers
    # outer resources: its outputs are then decoded at both levels
    flow = cones.Flow(program, two=True if "outer" in levels else None)
    check = cones.contained(program, flow, 4, 1)
    if check["outside"]:
        raise ValueError("linear propagation leaves a cone of this program")
    if flow.two and "outer" not in levels:
        raise ValueError("a two-level program needs the outer code's n and t")
    if flow.two and noise != "independent":
        raise ValueError("two levels are bounded under independent noise only")
    if not isinstance(spec["contracts"], dict):
        raise ValueError("contracts: {kind: {X: ..., Z: ...}} or {kind: {fail: ...}}")
    contracts = {}
    for kind, c in spec["contracts"].items():
        if isinstance(c, dict) and set(c) == {"fail"}:
            contracts[kind] = dict(fail=probability(c["fail"]))
            continue
        fields(c, "X", "Z")
        contracts[kind] = {s: _contract(c[s], f"{kind} {s}", noise) for s in "XZ"}

    def need(kind, sector=None):
        if kind not in contracts or (sector and sector not in contracts[kind]):
            raise ValueError(f"no contract for ancilla kind {kind!r}")
        return contracts[kind][sector] if sector else contracts[kind]

    n1, t1 = levels["inner"]
    cache, shapes, coarse = {}, {}, 0

    def inner_bound(mask):
        nonlocal coarse
        L, atoms = flow.weight(mask)
        key = (L, tuple(sorted(atoms.items())))
        if key not in cache:
            if noise == "independent":
                base = 1 - (1 - p) ** L
                lists = [need(k, s)["modes"] for (k, s), m in sorted(atoms.items()) for _ in range(m)]
                b, _, did = mixture_tail(n1, t1, up(base), lists)
                coarse += did
            else:
                weight = L + sum(need(k, s)["equivalents"] * m for (k, s), m in atoms.items())
                b = local_stochastic_tail(n1, t1, weight, p)
            cache[key] = b
        shapes[key] = shapes.get(key, 0) + 1
        return cache[key]

    def junk(kind):
        if noise != "independent":
            return Fraction(0)
        return need(kind, "X")["pi"] + need(kind, "Z")["pi"]

    probs, by_tag = [], {}
    for ev in flow.events:
        if ev["node"] is None:
            if ev["tag"] == "and1":
                pr = need("and1")["fail"] if "fail" in need("and1") else junk("and1")
            else:
                pr = junk(ev["tag"].split(":", 1)[1])
        else:
            pr = inner_bound(flow.cone(ev))
            if ev["kind"] == "ec" and noise == "independent":
                s = ev["sector"]
                pr += need("plus", s)["pi"] + need("zero", s)["pi"]
        probs.append(min(Fraction(1), up(pr)))
    terms = []                                                  # (tag, probability)
    if not flow.two:
        terms = [(ev["tag"], pr) for ev, pr in zip(flow.events, probs)]
    else:
        n2, t2 = levels["outer"]
        ev_of = {i: s["event"] for i, s in enumerate(flow.sources) if s["level"] == 2 and s["kind"] == "event"}
        for o in flow.outer:
            base, lists = Fraction(0), []
            for i in cones._bits(flow.cone_of[o["node"]]):
                s = flow.sources[i]
                if i in ev_of:
                    base += probs[ev_of[i]]
                elif s["kind"] == "atom":
                    lists.append(need(s["what"], s["sector"])["modes"])
            b, _, did = mixture_tail(n2, t2, min(Fraction(1), up(base)), lists)
            coarse += did
            terms.append((o["tag"], b))
        for tm in flow.terms:
            terms.append((tm["tag"], junk(tm["tag"].split(":", 1)[1])))
    total = Fraction(0)
    for tag, pr in terms:
        by_tag[tag] = by_tag.get(tag, Fraction(0)) + pr
        total += pr
    total = min(Fraction(1), up(total))
    cone_events = [pr for ev, pr in zip(flow.events, probs) if ev["node"] is not None]
    digest = hashlib.sha256(json.dumps([str(x) for x in probs] + [[t, str(x)] for t, x in terms]).encode()).hexdigest()
    premises = ["the cone theorem for this program (ftqc_gadgets.cones: hypotheses checked, induction quoted)"]
    if contracts:
        premises.append("the declared ancilla contracts hold" + (
            " and ancillas are independent of the program's locations and of one another" if noise == "independent" else ""))
    if flow.two:
        premises.append("failures of inner events in different columns are independent")
    premises.append("faults at different locations are independent, each with probability at most p"
                    if noise == "independent" else
                    "Pr[all locations of S faulty] <= p^|S| for every set S of physical locations")
    return {
        "status": "conditional",
        "rule": "tail bound over fault cones, summed over decoding events",
        "noise": noise,
        "p": str(p),
        "inner": dict(n=n1, t=t1),
        "outer": dict(n=levels["outer"][0], t=levels["outer"][1]) if flow.two else None,
        "program_sha256": program.sha256(),
        "cones": dict(sources=len(flow.sources), inner_decoding_events=len(cone_events),
                      outer_decoding_events=len(flow.outer), distinct_cone_shapes=len(cache),
                      algebra_settings=check["settings"]),
        "largest_inner_event": sci(max(cone_events, default=Fraction(0))),
        "terms": len(terms),
        "by_tag": {t: sci(v) for t, v in sorted(by_tag.items())},
        "upper": str(total),
        "upper_sci": sci(total),
        "bounds_sha256": digest,
        "rounding": f"upward to {SIG} significant digits",
        "mixtures_coarsened": coarse,
        "premises": premises,
    }, probs, terms


def cone_bound_result(spec):
    """Study kind `cone-bound`: the sum over the decoding events of a program."""
    return _compute(spec)[0]


def event_bounds(spec):
    """For comparisons: (probability bound of every inner event in order, [(tag, bound)] of
    the top-level terms), as Fractions."""
    _, probs, terms = _compute(spec)
    return probs, terms
