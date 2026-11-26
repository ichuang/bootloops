"""What faults and input errors do to a gadget, as vectors over GF(2).

Every fault and every input error has a signature with three parts, all linear:
    det   which postselected detectors fire,
    dec   which decoder checks are violated,
    syn   which stabilizers of the ideal joint state of references and outputs are
          violated, counting a Pauli frame made wrong by flipped outcomes.
A signature is packed into one integer; the signature of a set is the XOR.

The decoders are the only non-linear step.  The *effect* of an accepted run is
    effect = syn  XOR  (what each decoder's correction does to syn),
the corrections being looked up from the decoder checks.  Effect zero means: the
state of references and outputs is exactly the one the gadget is believed to
have produced.  Since the references record the input, that includes the logical
action.

`syn` is expressed in a canonical basis of the stabilizer group (reduced echelon
form over the interface: references in input order, then output qubits in block
order), so effects of two gadgets with the same ideal map are comparable bit by bit.
"""
import hashlib
from bisect import bisect_left
from itertools import combinations
from math import comb

from faultcert import stab
from faultcert.model import Model, propagate
from qeccert import gf2

from .gadget import GadgetError


class Budget(Exception):
    pass


def leader_table(checks, npos, cap=5_000_000):
    """Coset leaders of a parity-check matrix: {syndrome: position mask of a
    lightest pattern with that syndrome}, lightest first, ties by position order.
    Returns (table, unique_through): every pattern of weight <= unique_through is its
    own coset's only pattern of that weight or less."""
    col = [sum(1 << i for i, r in enumerate(checks) if (r >> p) & 1) for p in range(npos)]
    need = 1 << gf2.rank(checks)
    tab, w, clash = {0: 0}, 0, None
    while len(tab) < need and w < npos:
        w += 1
        if comb(npos, w) > cap:
            raise Budget(f"decoder table needs patterns of weight {w} on {npos} positions")
        for pos in combinations(range(npos), w):
            s = m = 0
            for p in pos:
                s ^= col[p]
                m |= 1 << p
            if s in tab:
                if clash is None:
                    clash = w
            else:
                tab[s] = m
    if clash is None and len(tab) >= need and w < npos:
        # every syndrome is taken; one more pattern of any weight must clash
        clash = w + 1
    return tab, (w if clash is None else clash - 1)


class Analysis:
    def __init__(self, gadget, cap=3_000_000):
        g = self.g = gadget
        self.cap = cap                         # budget for fault-subset enumeration
        self.small = min(cap, 5_000_000)       # budget for tables held in memory
        m = self.model = Model(g.circuit)
        c = g.circuit
        self._check_inputs()
        # ---- layout of a packed signature
        self.nd = len(m.detectors)
        self.det_masks = list(m.det_masks)
        self.dec_layout, self.dec_masks = [], []
        for d in g.decoders:
            off = len(self.dec_masks)
            for idx in d.bits:
                if m.sim.parity_form(idx) >> 1:
                    raise GadgetError(f"decoder {d.name}: a check is not deterministic in the "
                                      f"noiseless circuit")
                self.dec_masks.append(sum(1 << i for i in idx))
            self.dec_layout.append((off, len(d.checks)))
        self.nc = len(self.dec_masks)
        self.DET = (1 << self.nd) - 1
        self.shift = self.nd + self.nc
        # ---- canonical basis of the stabilizer group on the interface
        pos = {q: i for i, q in enumerate(g.interface)}
        mi = self.mi = len(g.interface)
        basis = []                                     # [vector, combination of generators, pivot]
        for i, (gx, gz, _) in enumerate(m.final):
            v = sum(1 << pos[q] for q in gf2.supp(gx)) | sum(1 << (mi + pos[q]) for q in gf2.supp(gz))
            cmb = 1 << i
            for bv, bc, p in basis:
                if (v >> p) & 1:
                    v ^= bv
                    cmb ^= bc
            if not v:
                raise GadgetError("dependent generators in the output group")
            p = v.bit_length() - 1
            for row in basis:
                if (row[0] >> p) & 1:
                    row[0] ^= v
                    row[1] ^= cmb
            basis.append([v, cmb, p])
        basis.sort(key=lambda r: r[2])
        self.canon = [r[0] for r in basis]
        self.combo = [r[1] for r in basis]
        self.pivot = [r[2] for r in basis]
        self.ns = len(basis)
        h = hashlib.sha256(repr((g.shape(), mi, self.canon)).encode()).hexdigest()
        self.group_id = h
        self.pure = self.ns == mi                      # references and outputs in a pure state
        # ---- per-qubit syndromes on the interface (no frame)
        self.sx = [self._canon(m.syndrome(1 << q, 0)) for q in g.interface]
        self.sz = [self._canon(m.syndrome(0, 1 << q)) for q in g.interface]
        self.ipos = pos
        # ---- columns
        in_locs = {i for b in g.inputs for i in b.in_locs}
        cols = {}
        for loc in m.locations:
            if loc.index in in_locs:
                continue
            for label, eff in loc.faults.items():
                s = self.sig(*eff)
                if s:
                    cols.setdefault(s, (loc.index, label))
        self.cols = sorted(cols)
        self.src = [cols[s] for s in self.cols]
        self.by_det = {}
        for i, s in enumerate(self.cols):
            self.by_det.setdefault(s & self.DET, []).append(i)
        self.in_sig = []                               # per input block, per position: (sig X, sig Z)
        for b in g.inputs:
            row = []
            for i in b.in_locs:
                f = m.locations[i].faults
                row.append((self.sig(*f["X"]), self.sig(*f["Z"])))
            self.in_sig.append(row)
        # ---- decoders
        self.dec = []
        self.dec_info = []
        for d, (off, width) in zip(g.decoders, self.dec_layout):
            tab, uniq = leader_table(d.checks, d.npos)
            per = []
            for p in range(d.npos):
                if d.kind == "record":
                    per.append(self._canon(m.syndrome(0, 0, 1 << d.targets[p])))
                elif d.kind == "X":
                    per.append(self._canon(m.syndrome(1 << d.targets[p], 0)))
                else:
                    per.append(self._canon(m.syndrome(0, 1 << d.targets[p])))
            eff = {}
            for s, lead in tab.items():
                e = 0
                for p in gf2.supp(lead):
                    e ^= per[p]
                eff[s] = e
            self.dec.append((off, (1 << width) - 1, eff))
            self.dec_info.append(dict(name=d.name, kind=d.kind, positions=d.npos,
                                      checks=width, syndromes=len(tab), unique_through=uniq,
                                      table=tab))
        self._effect_cache = {}
        self._sum_cache = {}

    # ------------------------------------------------------------------ inputs
    def _check_inputs(self):
        """The noiseless preamble must leave every input block maximally entangled with
        its references inside the block's code: on block and references together the
        state is pure, on the block alone it is the code's projector."""
        g, c = self.g, self.g.circuit
        if not g.inputs:
            return
        sim = stab.Sim(c.n)
        from faultcert.circuit import Circuit
        pre = Circuit([[("!" + op.name,) + op.qs for op in t] for t in c.ticks[:g.body_start]], c.n)
        if any(op.noisy for t in c.ticks[:g.body_start] for op in t):
            raise GadgetError("the input preamble contains a noisy operation")
        sim.run(pre)
        for b in g.inputs:
            mq = sum(1 << q for q in b.qubits)
            mr = sum(1 << q for q in b.refs)
            both = sim.output_group(mq | mr)
            alone = sim.output_group(mq)
            if len(both) != len(b.qubits) + len(b.refs) or sim.output_group(mr):
                raise GadgetError(f"input block {b.name}: not maximally entangled with its references")
            if len(alone) != len(b.qubits) - len(b.refs):
                raise GadgetError(f"input block {b.name}: wrong number of stabilizers")
            if b.code is not None:
                loc = b.qubits
                try:
                    for r in b.code.hx:
                        if sim._express(sum(1 << loc[p] for p in gf2.supp(r)), 0)[0]:
                            raise AssertionError
                    for r in b.code.hz:
                        if sim._express(0, sum(1 << loc[p] for p in gf2.supp(r)))[0]:
                            raise AssertionError
                except AssertionError:
                    raise GadgetError(f"input block {b.name}: the preamble does not encode "
                                      f"into the declared code") from None
                if gf2.rank(b.code.hx) + gf2.rank(b.code.hz) != len(alone):
                    raise GadgetError(f"input block {b.name}: the preamble does not encode "
                                      f"into the declared code")

    # ------------------------------------------------------------------ signatures
    def _canon(self, raw):
        e = 0
        for i, cmb in enumerate(self.combo):
            if (raw & cmb).bit_count() & 1:
                e |= 1 << i
        return e

    def sig(self, x, z, flips):
        s = 0
        for i, mk in enumerate(self.det_masks):
            if (flips & mk).bit_count() & 1:
                s |= 1 << i
        for i, mk in enumerate(self.dec_masks):
            if (flips & mk).bit_count() & 1:
                s |= 1 << (self.nd + i)
        return s | (self._canon(self.model.syndrome(x, z, flips)) << self.shift)

    def parts(self, sig):
        return sig & self.DET, (sig >> self.nd) & ((1 << self.nc) - 1), sig >> self.shift

    def effect(self, sig):
        """Syndrome after the decoders have acted; -1 if a decoder has no entry."""
        e = self._effect_cache.get(sig)
        if e is None:
            e = sig >> self.shift
            d = sig >> self.nd
            for off, mask, eff in self.dec:
                c = eff.get((d >> off) & mask)
                if c is None:
                    e = -1
                    break
                e ^= c
            if len(self._effect_cache) < 2_000_000:
                self._effect_cache[sig] = e
        return e

    def pauli_syn(self, px, pz):
        """Canonical syndrome of a Pauli given as interface-position masks."""
        e = 0
        for p in gf2.supp(px):
            e ^= self.sx[p]
        for p in gf2.supp(pz):
            e ^= self.sz[p]
        return e

    def replay(self, faults, input_error=()):
        """Signature of a fault set [(location, label)] together with an input error
        [(block, position, Pauli)], recomputed by propagation (not from the columns)."""
        m, c = self.model, self.g.circuit
        in_locs = {i for b in self.g.inputs for i in b.in_locs}
        items = []
        for i, label in faults:
            if i in in_locs:
                raise ValueError(f"location {i} is an input-error location, not a fault location")
            items.append((int(i), label))
        for b, p, pauli in input_error:
            if pauli not in ("X", "Y", "Z"):
                raise ValueError("input errors are X, Y or Z on one position")
            items.append((self.g.inputs[b].in_locs[p], pauli))
        if len({i for i, _ in items}) != len(items):
            raise ValueError("a location appears twice")
        x = z = fl = 0
        for i, label in items:
            loc = m.locations[i]
            if loc.op.name in ("M", "MX"):
                if label != "flip":
                    raise ValueError("measurement locations take the label 'flip'")
                fl ^= loc.faults["flip"][2]
                continue
            if len(label) != len(loc.op.qs) or set(label) - set("IXYZ") or set(label) == {"I"}:
                raise ValueError(f"bad Pauli label {label!r} for location {i}")
            px = sum(1 << q for q, p in zip(loc.op.qs, label) if p in "XY")
            pz = sum(1 << q for q, p in zip(loc.op.qs, label) if p in "ZY")
            fx, fz, ff = propagate(c, loc.tick, px, pz)
            x ^= fx
            z ^= fz
            fl ^= ff
        return self.sig(x, z, fl)

    # ------------------------------------------------------------------ enumeration
    def sums(self, j, det_target=0):
        """Yield (sum, column indices) over the j-subsets of fault columns whose detector
        part equals det_target: choose j - 1 columns, look the last one up."""
        cols, by_det, DET = self.cols, self.by_det, self.DET
        if j == 0:
            if det_target == 0:
                yield 0, ()
            return
        # cost: the (j - 1)-subsets, plus one step per subset found.  Without detectors
        # every j-subset is found; with detectors the number found is counted as it grows.
        if comb(len(cols), j - 1) > self.cap or (not self.nd and comb(len(cols), j) > 40 * self.cap):
            raise Budget(f"{comb(len(cols), j - 1)} subsets of size {j - 1} of {len(cols)} fault "
                         f"effects exceed the budget")
        stack = [(0, 0, ())]
        n = len(cols)
        found, limit = 0, 40 * self.cap
        while stack:
            start, acc, idx = stack.pop()
            if len(idx) == j - 1:
                lst = by_det.get((acc ^ det_target) & DET)
                if lst:
                    for k in range(bisect_left(lst, start), len(lst)):
                        last = lst[k]
                        yield acc ^ cols[last], idx + (last,)
                    found += len(lst)
                    if found > limit:
                        raise Budget(f"more than {limit} accepted subsets of size {j}")
                continue
            for i in range(n - 1, start - 1, -1):
                stack.append((i + 1, acc ^ cols[i], idx + (i,)))

    def levels(self, s, det_target=0):
        """[L_0, ..., L_s]: L_j is the set of sums of accepted fault sets of j columns
        that no smaller accepted set produces."""
        lv, seen = self._sum_cache.setdefault(det_target, ([], set()))
        while len(lv) <= s:
            cur = set()
            for v, _ in self.sums(len(lv), det_target):
                if v not in seen:
                    cur.add(v)
            seen |= cur
            lv.append(cur)
        return lv[:s + 1]

    def find(self, j, det_target, value):
        for v, idx in self.sums(j, det_target):
            if v == value:
                return idx
        return None

    def faults_of(self, idx):
        return [[self.src[i][0], self.src[i][1]] for i in idx]

    # ------------------------------------------------------------------ input errors
    def input_errors(self, r, sector=False):
        """{signature: (weight, [(block, position, Pauli)])} over input errors of weight
        at most r, lightest representative first.  Weight is Pauli weight over all
        input blocks, or with sector=True the larger of the X and Z weights."""
        cells = [(b, p, sx, sz) for b, row in enumerate(self.in_sig) for p, (sx, sz) in enumerate(row)]
        out = {0: (0, [])}
        if not sector:
            for w in range(1, r + 1):
                if comb(len(cells), w) * 3 ** w > self.small:
                    raise Budget(f"too many input errors of weight {w}")
                for sub in combinations(cells, w):
                    def rec(i, acc, lab):
                        if i == w:
                            if acc not in out:
                                out[acc] = (w, lab)
                            return
                        b, p, sx, sz = sub[i]
                        rec(i + 1, acc ^ sx, lab + [(b, p, "X")])
                        rec(i + 1, acc ^ sz, lab + [(b, p, "Z")])
                        rec(i + 1, acc ^ sx ^ sz, lab + [(b, p, "Y")])
                    rec(0, 0, [])
            return out
        one = {"X": {0: (0, [])}, "Z": {0: (0, [])}}
        for t, k in (("X", 2), ("Z", 3)):
            for w in range(1, r + 1):
                for sub in combinations(cells, w):
                    acc = 0
                    for c in sub:
                        acc ^= c[k]
                    one[t].setdefault(acc, (w, [(c[0], c[1], t) for c in sub]))
        if len(one["X"]) * len(one["Z"]) > self.small:
            raise Budget("too many input errors")
        for ax, (wx, lx) in one["X"].items():
            for az, (wz, lz) in one["Z"].items():
                w = max(wx, wz)
                cur = out.get(ax ^ az)
                if cur is None or cur[0] > w:
                    zs = {(b, p) for b, p, _ in lz}
                    xs = {(b, p) for b, p, _ in lx}
                    out[ax ^ az] = (w, [(b, p, "Y" if (b, p) in zs else "X") for b, p, _ in lx]
                                    + [(b, p, "Z") for b, p, _ in lz if (b, p) not in xs])
        return out

    # ------------------------------------------------------------------ output weight
    def weigher(self, sector=False):
        return Weigher(self, sector)

    def stats(self):
        c = self.g.counts()
        return dict(qubits=c["qubits"], ticks=c["ticks"], fault_locations=c["fault_locations"],
                    input_error_locations=c["input_error_locations"], columns=len(self.cols),
                    detectors=self.nd, decoder_checks=self.nc, stabilizers=self.ns,
                    interface=self.mi, pure=self.pure, group_id=self.group_id,
                    decoders=[{k: v for k, v in d.items() if k != "table"} for d in self.dec_info])


class Weigher:
    """Distance of an effect from zero: the least w such that a Pauli of weight at most
    w on every output block (never on a reference) has that effect.  With sector=True
    the weight of a Pauli is the larger of its X weight and Z weight.

    A block with a code is handled through its own syndrome: the effect fixes which of
    the block's stabilizers are violated, a table gives the lightest Pauli on the block
    with that pattern, and after removing one such Pauli per block the effect must be
    zero.  That is exact as long as the lightest Pauli is unique up to the stabilizer
    group, which the table construction checks; beyond that radius it refuses."""

    def __init__(self, an, sector):
        self.an, self.sector = an, sector
        g = an.g
        self.blocks = []
        base = sum(len(b.refs) for b in g.inputs)
        for b in g.outputs:
            ipos = list(range(base, base + len(b.qubits)))
            base += len(b.qubits)
            keys = None
            if b.code is not None:
                keys = dict(X=[self._functional(ipos, 0, r) for r in b.code.hz],   # seen by X errors
                            Z=[self._functional(ipos, r, 0) for r in b.code.hx])
            self.blocks.append(dict(ipos=ipos, keys=keys, tabs={}, built=-1))
        if len(self.blocks) > 1 and any(b["keys"] is None for b in self.blocks):
            raise GadgetError("output weight over several blocks needs a code on each block")
        self.cache = {}

    def _functional(self, ipos, hx_row, hz_row):
        """Mask over canonical syndrome bits whose parity says whether the block
        stabilizer (X on hx_row, Z on hz_row) is violated."""
        an = self.an
        v = sum(1 << ipos[p] for p in gf2.supp(hx_row)) | sum(1 << (an.mi + ipos[p]) for p in gf2.supp(hz_row))
        f = 0
        for i, (cv, p) in enumerate(zip(an.canon, an.pivot)):
            if (v >> p) & 1:
                v ^= cv
                f |= 1 << i
        if v:
            raise GadgetError("an output block is not in its declared code in the noiseless gadget")
        return f

    @staticmethod
    def _key(e, funcs):
        k = 0
        for i, f in enumerate(funcs):
            if (e & f).bit_count() & 1:
                k |= 1 << i
        return k

    def _extend(self, blk, w):
        an = self.an
        ip = blk["ipos"]
        for cur in range(blk["built"] + 1, w + 1):
            if blk["keys"] is None:
                tab = blk["tabs"].setdefault("P", {})
                self._paulis(ip, cur, lambda e, t=tab, c=cur: t.setdefault(e, (c, e)))
            elif self.sector:
                for t, per, funcs in (("X", an.sx, blk["keys"]["X"]), ("Z", an.sz, blk["keys"]["Z"])):
                    tab = blk["tabs"].setdefault(t, {})
                    if comb(len(ip), cur) > an.small:
                        raise Budget(f"weight table at weight {cur}")
                    for sub in combinations(ip, cur):
                        e = 0
                        for p in sub:
                            e ^= per[p]
                        self._store(tab, self._key(e, funcs), cur, e)
            else:
                tab = blk["tabs"].setdefault("P", {})
                kx, kz = blk["keys"]["X"], blk["keys"]["Z"]
                nx = len(kx)
                self._paulis(ip, cur, lambda e, t=tab, c=cur: self._store(
                    t, self._key(e, kx) | (self._key(e, kz) << nx), c, e))
            blk["built"] = cur

    def _paulis(self, ip, w, put):
        an = self.an
        if comb(len(ip), w) * 3 ** w > an.small:
            raise Budget(f"weight table at weight {w}")
        for sub in combinations(ip, w):
            def rec(i, acc):
                if i == w:
                    put(acc)
                    return
                p = sub[i]
                rec(i + 1, acc ^ an.sx[p])
                rec(i + 1, acc ^ an.sz[p])
                rec(i + 1, acc ^ an.sx[p] ^ an.sz[p])
            rec(0, 0)

    @staticmethod
    def _store(tab, key, w, e):
        cur = tab.get(key)
        if cur is None:
            tab[key] = (w, e)
        elif cur[1] != e:
            raise Budget(f"output weight {w} is beyond the radius in which the lightest error "
                         f"of a block is unique")

    def weight(self, e, upto):
        """Least w <= upto with e within weight w of zero on every block; None if there
        is none (a logical error, or an error heavier than upto)."""
        if e == -1:
            return None
        if e == 0:
            return 0
        key = (e, upto)
        if key in self.cache:
            return self.cache[key]
        worst, rest = 0, e
        for blk in self.blocks:
            if blk["built"] < upto:
                self._extend(blk, upto)
            if blk["keys"] is None:
                hit = blk["tabs"]["P"].get(rest)
                if hit is None:
                    worst = None
                    break
                worst = max(worst, hit[0])
                rest ^= hit[1]
            elif self.sector:
                for t in ("X", "Z"):
                    hit = blk["tabs"][t].get(self._key(e, blk["keys"][t]))
                    if hit is None:
                        worst = None
                        break
                    worst = max(worst, hit[0])
                    rest ^= hit[1]
                if worst is None:
                    break
            else:
                kx, kz = blk["keys"]["X"], blk["keys"]["Z"]
                hit = blk["tabs"]["P"].get(self._key(e, kx) | (self._key(e, kz) << len(kx)))
                if hit is None:
                    worst = None
                    break
                worst = max(worst, hit[0])
                rest ^= hit[1]
        if worst is not None and (rest or worst > upto):
            worst = None
        if len(self.cache) < 2_000_000:
            self.cache[key] = worst
        return worst
