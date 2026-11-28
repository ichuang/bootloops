"""Small exact models against which the filter calculus is tested.

`indicator_model` follows the filter rule itself on a few coordinates: every block carries a
junk flag and a set of dirty coordinates per sector, checks combine them by the rule, and the
law of the accepted block is computed by summing over everything.  It has two variants.  In
the independent one a step has a gate, a measurement and an idle location per coordinate.  In
the shared one every block has, at every coordinate, a source bit left by its latest step;
the fresh dirt of a step is its new source (variant "shared") or the new source united with
the survivor's previous one (variant "overlap"), and the partner's side depends on the
survivor's previous source as well, which is how decoding events of an encoded tree depend
on the ancillas of the preceding correction.

Acceptance is a parameter, because the calculus claims its bound for every acceptance event
that contains "everything clean": "always" accepts everything; "clean" accepts only when
neither incoming block nor the partner's side of the step is dirty; "adverse" accepts when nothing is dirty or the survivor's checked sector
ends up dirty or junk, which makes conditioning as harmful as it can be.

Probabilities are floats (the tests compare with a margin) or Fractions.
"""
from itertools import product


def _cols(n, per):
    """All assignments of per-coordinate outcomes: [(tuple of outcomes, probability)]."""
    out = [((), 1)]
    for _ in range(n):
        out = [(a + (x,), pa * px) for a, pa in out for x, px in per]
    return out


def leaf_law(n, variant, pi, rate):
    """Law of a raw block: junk flags of the two sectors (independent, probabilities pi[X],
    pi[Z]); independent variant: per coordinate, dirt in both sectors with probability
    `rate`; otherwise source bits s0, s1 of rate `rate`, dirt s1 ("shared") or s0 or s1
    ("overlap") in both sectors, s1 kept as the block's current source."""
    one = type(rate)(1)
    if variant == "independent":
        per = [((1, 1, 0), rate), ((0, 0, 0), one - rate)]
    else:
        per = []
        for s0, s1 in product((0, 1), repeat=2):
            d = (s0 | s1) if variant == "overlap" else s1
            per.append(((d, d, s1), (rate if s0 else one - rate) * (rate if s1 else one - rate)))
    law = {}
    for cols, pc in _cols(n, per):
        for jx, jz in product((0, 1), repeat=2):
            pr = pc * (pi["X"] if jx else one - pi["X"]) * (pi["Z"] if jz else one - pi["Z"])
            if pr:
                key = (jx, jz, cols)
                law[key] = law.get(key, 0) + pr
    return law


def step(law_a, law_b, s, n, wmin, variant, rates, accept):
    """Law of the survivor of one check of sector s given acceptance, and the probability of
    acceptance.  rates: independent variant dict(gate, meas, idle); otherwise
    dict(source, meas)."""
    si = 0 if s == "X" else 1
    one = type(next(iter(law_a.values())))(1)
    if variant == "independent":
        g, m, i = rates["gate"], rates["meas"], rates["idle"]
        per = [((a, b, c), (g if a else one - g) * (m if b else one - m) * (i if c else one - i))
               for a, b, c in product((0, 1), repeat=3)]
    else:
        g, m = rates["source"], rates["meas"]
        per = [((a, b), (g if a else one - g) * (m if b else one - m)) for a, b in product((0, 1), repeat=2)]
    rand = _cols(n, per)
    out, p_acc = {}, 0
    for (ja_x, ja_z, ca), pa in law_a.items():
        ja = (ja_x, ja_z)
        for (jb_x, jb_z, cb), pb in law_b.items():
            jb = (jb_x, jb_z)
            pab = pa * pb
            for r, pr in rand:
                cols, da, dbp, fresh, allclean = [], 0, 0, 0, not ja[si] and not jb[si]
                new = []
                for c in range(n):
                    if variant == "independent":
                        gate, meas, idle = r[c]
                        f, mm, src = gate | idle, gate | meas, 0
                    else:
                        s_new, m0 = r[c]
                        f = (ca[c][2] | s_new) if variant == "overlap" else s_new
                        mm, src = ca[c][2] | cb[c][2] | m0, s_new
                    xa, xb = ca[c][si], cb[c][si]
                    da += xa
                    dbp += xb | mm
                    fresh += f
                    if xa or xb or mm:
                        allclean = False
                    new.append((f, mm, src, xa, xb))
                # the filter rule
                if ja[si] and jb[si]:
                    js = 1
                elif not ja[si] and not jb[si]:
                    union = sum(1 for (f, mm, src, xa, xb) in new if xa or xb or mm)
                    js = int(union >= wmin)
                else:
                    js = 0
                jo = ja[1 - si] | jb[1 - si]
                dirty_s = 0
                for c, (f, mm, src, xa, xb) in enumerate(new):
                    if ja[si] and not jb[si]:
                        ds = f | xb | mm
                    elif jb[si] and not ja[si]:
                        ds = f | xa
                    else:
                        ds = f | (xa & (xb | mm))
                    do = ca[c][1 - si] | cb[c][1 - si] | f
                    dirty_s += ds
                    cols.append((ds, do, src) if si == 0 else (do, ds, src))
                if accept == "always":
                    acc = True
                elif accept == "clean":
                    acc = allclean
                else:
                    acc = allclean or js or dirty_s > 0
                if not acc:
                    continue
                pr_all = pab * pr
                p_acc += pr_all
                key = ((js, jo) if si == 0 else (jo, js)) + (tuple(cols),)
                out[key] = out.get(key, 0) + pr_all
    return {k: v / p_acc for k, v in out.items()}, p_acc


def increasing_events(n):
    """Every increasing event of n coordinates, as the set of dirty-coordinate masks at which
    it holds (n <= 3)."""
    masks = list(range(1 << n))
    events = []
    for bits in product((0, 1), repeat=len(masks)):
        ok = all(not bits[a] or bits[b] for a in masks for b in masks if a & b == a)
        if ok and any(bits) and not all(bits):
            events.append(frozenset(m for m in masks if bits[m]))
    return events


def event_probability(law, sector, event):
    """P[not junk in the sector, the sector's dirty set in `event`]."""
    si = 0 if sector == "X" else 1
    tot = 0
    for key, pr in law.items():
        if key[si]:
            continue
        mask = sum(1 << c for c, col in enumerate(key[2]) if col[si])
        if mask in event:
            tot += pr
    return tot


def junk_probability(law, sector):
    si = 0 if sector == "X" else 1
    return sum(pr for key, pr in law.items() if key[si])


def iid_probability(n, q, event):
    """Probability of the event when every coordinate is dirty independently with rate q."""
    tot = 0
    for mask in event:
        k = bin(mask).count("1")
        tot += q ** k * (1 - q) ** (n - k)
    return tot
