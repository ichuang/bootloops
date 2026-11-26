"""Two levels of concatenation.

At the outer level a qubit is an inner code block and an operation is an inner gadget.
What the outer level needs from below is small: an inner gadget that works is an ideal
logical operation, and one that does not is some error on its own outer coordinate.
What it needs of itself is the statement "at most T bad coordinates, hence correct",
which is the coordinate-locality certificate of an outer gadget (props.local): that
certificate counts coordinates, not faults, so it already allows anything at all to
happen inside a bad coordinate.

reduction      puts the two together as a receipt: the outer gadget's locality
               certificate, and for every kind of operation the outer circuit uses, a
               receipt that the inner gadget has the property that makes it a good
               location.  The step from "inner properties" to "a column with few faults
               acts ideally" is the level-reduction argument of the literature (or of a
               monograph's own calculus) and is quoted, not re-proved.
two_level      builds the same gadget directly on the concatenated code, physical qubit
               by physical qubit, with hierarchical decoders (ftqc_codes.hier).
Direct         tests the statement on that physical gadget: fault sets and input errors
               in which at most t2 columns hold more than t1 of them must act ideally.
               Exhaustive at low weight, sampled above, with a search for the first
               failure when two columns are bad.
"""
import random
from itertools import combinations

import ftqc_codes
from ftqc_codes import gf2, hier
from ftqc_codes.codes import concatenate

from . import certify, library, props
from .effects import Analysis
from .gadget import Decoder, Gadget, GadgetError, logical_x

LEVELS = "ftqc_gadgets-levels-1"

STATEMENT = (
    "Replace every qubit of the outer gadget by an inner code block and every operation by the "
    "inner gadget of its receipt (inner error correction where the design puts it).  The physical "
    "qubits of one outer coordinate, across all blocks and ancillas, form a column; the outer "
    "gadget is coordinate-local, so no gate joins two columns.  Call a column good if its inner "
    "gadgets act as ideal logical operations, and bad otherwise.  Whatever happens in a bad "
    "column is an operation on the qubits of one outer coordinate over the course of the gadget, "
    "hence a combination of Pauli faults at that coordinate's outer locations.  By the outer "
    "locality certificate, if at most T2 columns are bad (counting columns in which the input "
    "carries a logical error), the outer gadget has its certified properties.  The inner receipts "
    "certify, for every kind of operation used, the property that level reduction needs (prep, "
    "gate, ec with close, decoded measurement) for input weight r and s faults with r + s <= t1.  "
    "That these properties make a column good when each of its extended rectangles has at most "
    "t1 faults is Aliferis, Gottesman and Preskill (2006), Lemmas 3 to 5; a design without fixed "
    "rectangles must supply its own criterion for a good column.  The receipt records the "
    "hypotheses; it does not re-prove that step.  Outer parts listed as assumed (ancillas given "
    "as specifications) are hypotheses as well: an outer ancilla must refine its specification "
    "with bad columns in place of faults.")

ROLE_OF = {"R": "prep_zero", "RX": "prep_plus", "CX": "cnot", "M": "measure_z", "MX": "measure_x"}
KIND_OF = dict(prep_zero="prep", prep_plus="prep", cnot="gate", ec="ec", measure_z="ec", measure_x="ec")


# ---------------------------------------------------------------------- physical two-level gadgets
def positions(g):
    """qubit -> position within its block, for gadgets whose blocks were allocated whole
    (the library's): the qubits other than references, in order, n at a time.  Refuses a
    gadget in which some noisy operation joins two positions."""
    n = len((g.inputs + g.outputs)[0].qubits)
    refs = {q for b in g.inputs for q in b.refs}
    qs = [q for q in range(g.circuit.n) if q not in refs]
    if len(qs) % n:
        raise GadgetError("qubits do not come in whole blocks")
    pos = {q: i % n for i, q in enumerate(qs)}
    for k, t in enumerate(g.circuit.ticks):
        if k < g.body_start:
            continue
        for op in t:
            if op.noisy and len({pos[q] for q in op.qs}) != 1:
                raise GadgetError(f"operation {op.text()} joins two positions")
    return pos, n


def hierarchical(g, outer, inner):
    """The same gadget with every decoder replaced by the hierarchical one."""
    cc = concatenate(outer, inner)
    decs = []
    for d in g.decoders:
        if d.checks == list(cc.hz):
            sector = "X"
        elif d.checks == list(cc.hx):
            sector = "Z"
        else:
            raise GadgetError(f"decoder {d.name}: its checks are not those of the concatenated code")
        rule = dict(type="hier", sector=sector, outer=outer.to_json(), inner=inner.to_json())
        decs.append(Decoder(d.name, d.checks, d.bits, d.kind, d.targets, rule))
    return Gadget(g.circuit, g.inputs, g.outputs, decs, g.body_start, g.name + "+hier", g.parts)


def two_level(family, outer, inner, **kw):
    """library.<family> on the concatenated code, with hierarchical decoders."""
    outer, inner = library.code_of(outer), library.code_of(inner)
    cc = concatenate(outer, inner)
    return hierarchical(getattr(library, family)(cc, **kw), outer, inner)


class Direct:
    def __init__(self, gadget, outer, inner, cap=3_000_000):
        self.g, self.outer, self.inner = gadget, library.code_of(outer), library.code_of(inner)
        self.m, self.N = self.inner.n, self.outer.n
        self.an = an = Analysis(gadget, cap)
        pos, n = positions(gadget)
        if n != self.m * self.N:
            raise GadgetError("the gadget's blocks are not blocks of the concatenated code")
        in_locs = {i for b in gadget.inputs for i in b.in_locs}
        self.faults = []                       # (column, location, label, signature)
        for loc in an.model.locations:
            if loc.index in in_locs:
                continue
            col = pos[loc.op.qs[0]] // self.m
            for label, eff in loc.faults.items():
                self.faults.append((col, loc.index, label, an.sig(*eff)))
        self.inputs = []                       # (column, (block, position), Pauli, signature)
        for b, row in enumerate(an.in_sig):
            for p, (sx, sz) in enumerate(row):
                for t, sg in (("X", sx), ("Z", sz), ("Y", sx ^ sz)):
                    self.inputs.append((p // self.m, (b, p), t, sg))
        self.hx = hier.Hier(self.outer, self.inner, "X")
        self.hz = hier.Hier(self.outer, self.inner, "Z")
        self.t1 = min(self.hx.t1, self.hz.t1)
        self.t2 = min(self.hx.t2, self.hz.t2)
        W = an.weigher() if gadget.outputs else None
        self.blocks = [] if W is None else W.blocks
        self._ok = {}

    def ideal_ok(self, effect):
        """After an ideal hierarchical decoder on every output block, is the state of
        references and outputs the right one?"""
        if effect == -1:
            return False
        r = self._ok.get(effect)
        if r is None:
            e = effect
            for blk in self.blocks:
                for t, h, per in (("X", self.hx, self.an.sx), ("Z", self.hz, self.an.sz)):
                    vals = 0
                    for i, f in enumerate(blk["keys"][t]):
                        if (effect & f).bit_count() & 1:
                            vals |= 1 << i
                    pat = h.decode(vals)
                    if pat is None:
                        e = -1
                        break
                    for p in gf2.supp(pat):
                        e ^= per[blk["ipos"][p]]
                if e == -1:
                    break
            r = e == 0
            if len(self._ok) < 2_000_000:
                self._ok[effect] = r
        return r

    def verdict(self, items):
        """items: fault and input entries (distinct locations and positions).  Returns
        (number of bad columns, acts ideally?)."""
        dirt, sg = {}, 0
        for it in items:
            dirt[it[0]] = dirt.get(it[0], 0) + 1
            sg ^= it[3]
        bad = sum(1 for v in dirt.values() if v > self.t1)
        return bad, self.ideal_ok(self.an.effect(sg))

    @staticmethod
    def _distinct(items):
        return len({it[1] for it in items}) == len(items)

    def run(self, exhaustive=2, samples=20000, seed=1, search=200000):
        """exhaustive: every set of at most that many faults and input errors.  samples:
        random sets with at most t2 bad columns (each holding 2 to 4 items) and at most
        t1 items in each of a random choice of the other columns.  search: random sets
        with t2 + 1 bad columns, looking for one that does not act ideally."""
        rnd = random.Random(seed)
        pool = self.faults + self.inputs
        by_col = {}
        for it in pool:
            by_col.setdefault(it[0], []).append(it)
        out = dict(t1=self.t1, t2=self.t2, columns=self.N, fault_labels=len(self.faults),
                   input_labels=len(self.inputs), exhaustive_weight=exhaustive, exhaustive_sets=0,
                   sampled_sets=0, failures=0, first_failure=None, beyond=None)

        def note(items, ok, bad):
            if bad <= self.t2 and not ok:
                out["failures"] += 1
                if out["first_failure"] is None:
                    out["first_failure"] = [[it[1], it[2]] for it in items]

        for w in range(exhaustive + 1):
            for items in combinations(pool, w):
                if not self._distinct(items):
                    continue
                bad, ok = self.verdict(items)
                out["exhaustive_sets"] += 1
                note(items, ok, bad)
        cols = sorted(by_col)
        for _ in range(samples):
            items = []
            bads = rnd.sample(cols, rnd.randint(0, self.t2))
            for c in bads:
                items += rnd.sample(by_col[c], rnd.randint(self.t1 + 1, self.t1 + 3))
            for c in cols:
                if c not in bads and rnd.random() < 0.6:
                    items += rnd.sample(by_col[c], rnd.randint(1, self.t1))
            if not self._distinct(items):
                continue
            bad, ok = self.verdict(items)
            out["sampled_sets"] += 1
            note(items, ok, bad)
        tried = 0
        for tried in range(1, search + 1):
            items = []
            for c in rnd.sample(cols, self.t2 + 1):
                items += rnd.sample(by_col[c], self.t1 + 1)
            if not self._distinct(items):
                continue
            bad, ok = self.verdict(items)
            if bad == self.t2 + 1 and not ok:
                out["beyond"] = dict(bad_columns=bad, items=[[it[1], it[2]] for it in items], tries=tried)
                break
        return out


# ---------------------------------------------------------------------- reduction receipts
def supports(receipt, kind):
    """Largest T such that the receipt establishes property `kind` for r + s <= T
    (Pauli weight), from an enumeration, a locality certificate or a composition;
    -1 if it does not.  A refinement of a preparation specification counts as `prep`."""
    best = -1
    if receipt.get("schema") in (certify.COMPOSITION, "gadgetcert-composition-1"):
        for c in receipt["conclusions"]:
            if c["kind"] == kind and c.get("weight") in ("pauli", None):
                best = max(best, c["holds_through"])
            if kind == "prep" and c["type"] == "refine":
                best = max(best, c["holds_through"])
        return best
    for r in receipt.get("results", []):
        if r["type"] == "check" and r["kind"] == kind and r["weight"] == "pauli":
            best = max(best, r["holds_through"])
        elif r["type"] == "close" and kind == "close" and r.get("weight", "pauli") == "pauli":
            best = max(best, r["holds_through"])
        elif r["type"] == "refine" and kind == "prep" and r["input_weight"] == 0:
            best = max(best, r["holds_through"])
        elif r["type"] == "local" and r["local"] and r.get(kind if kind != "prep" else "none"):
            best = max(best, 10 ** 6 if r["through"] is None else r["through"])
    return best


def _block_code(receipt):
    if "gadget" not in receipt:
        return None
    g = Gadget.from_json(receipt["gadget"])
    codes = [b.code for b in g.inputs + g.outputs if b.code is not None]
    return codes[0] if codes else None


def reduction(outer_receipt, inner, cap=3_000_000):
    """outer_receipt: a receipt of the outer gadget that carries a locality certificate.
    inner: {role: receipt} with roles prep_zero, prep_plus, cnot, ec, measure_z,
    measure_x; every role the outer circuit uses is required, `ec` always.  Returns the
    levels receipt; raises Refusal if a hypothesis is missing."""
    certify.verify(outer_receipt, cap=cap)
    loc = [r for r in outer_receipt["results"] if r["type"] == "local"]
    if len(loc) != 1 or not loc[0]["local"]:
        raise certify.Refusal("the outer receipt carries no positive locality certificate")
    loc = loc[0]
    g = Gadget.from_json(outer_receipt["gadget"])
    # noiseless operations belong to parts given as specifications; they are listed as
    # assumptions, not matched with inner gadgets
    used = sorted({ROLE_OF[op.name] for k, t in enumerate(g.circuit.ticks) if k >= g.body_start
                   for op in t if op.noisy and op.name in ROLE_OF} | {"ec"})
    assumed = sorted({f"{p['role']} ancilla as its specification" for p in g.parts
                      if p.get("method") == "spec"})
    rows, t1, code_sha = [], None, None
    for role in used:
        rec = inner.get(role)
        if rec is None:
            raise certify.Refusal(f"no inner receipt for {role}")
        if rec.get("schema") in (certify.COMPOSITION, "gadgetcert-composition-1"):
            raise certify.Refusal(f"{role}: pass the composition as dict(composition=, target=, parts=)")
        if "composition" in rec:
            certify.verify_composition(rec["composition"], rec["target"], rec["parts"], cap)
            evid, base = rec["composition"], rec["target"]
        else:
            certify.verify(rec, cap=cap)
            evid = base = rec
        need = [KIND_OF[role]] + (["close"] if role == "ec" else [])
        got = {k: supports(evid, k) for k in need}
        through = min(got.values())
        if through < 0:
            raise certify.Refusal(f"{role}: the receipt does not establish {' and '.join(need)}")
        bg = Gadget.from_json(base["gadget"])
        if evid is base and any(p.get("method") == "spec" for p in bg.parts):
            raise certify.Refusal(f"{role}: the inner gadget was analysed with specified parts; pass its "
                                  f"composition with real parts as dict(composition=, target=, parts=)")
        if role.startswith("measure") and bg.outputs:
            raise certify.Refusal(f"{role}: a decoded measurement has no output block")
        code = _block_code(base)
        if code is None or len(logical_x(code)) != 1:
            raise certify.Refusal(f"{role}: the inner gadget must act on blocks of a code with k = 1")
        sha = (code.n, tuple(code.hx), tuple(code.hz))
        if code_sha is None:
            code_sha = sha
        elif sha != code_sha:
            raise certify.Refusal(f"{role}: the inner receipts are not all about the same code")
        rows.append(dict(role=role, gadget=bg.name, properties=got, through=through,
                         receipt_sha256=certify.receipt_sha(evid)))
        t1 = through if t1 is None else min(t1, through)
    return certify._clean(dict(
        schema=LEVELS, tool="ftqc_gadgets", tool_version=certify.VERSION,
        outer=dict(gadget=g.name, sha256=g.sha256(), receipt_sha256=certify.receipt_sha(outer_receipt),
                   coordinates=loc["coordinates"], properties=[k for k in ("gate", "ec", "close") if loc[k]],
                   decoders=loc["decoders"], assumed=assumed),
        inner=rows, inner_code=dict(n=code_sha[0]),
        conclusion=dict(bad_columns_tolerated=loc["through"], inner_faults_per_rectangle=min(t1, 10 ** 6)),
        statement=STATEMENT))


def verify_reduction(levels, outer_receipt, inner, cap=3_000_000):
    if levels.get("schema") != LEVELS:
        raise certify.Refusal(f"not a {LEVELS} receipt")
    again = reduction(outer_receipt, inner, cap)
    for k in ("outer", "inner", "inner_code", "conclusion", "statement"):
        if again[k] != levels.get(k):
            raise certify.Refusal(f"levels: {k} differs from what the receipts support")
    return ["outer receipt and locality certificate", "inner receipts", "conclusion recomputed"]


def summary(levels):
    c = levels["conclusion"]
    t2 = "any number of" if c["bad_columns_tolerated"] is None else f"at most {c['bad_columns_tolerated']}"
    lines = [f"{levels['outer']['gadget']} over inner blocks of {levels['inner_code']['n']} qubits: "
             f"{', '.join(levels['outer']['properties'])} with {t2} bad columns; a column is good "
             f"when each inner rectangle has at most {c['inner_faults_per_rectangle']} faults"]
    for a in levels["outer"].get("assumed", []):
        lines.append(f"  assumed: {a}")
    for r in levels["inner"]:
        lines.append(f"  {r['role']}: {r['gadget']} ({', '.join(f'{k} through {v}' for k, v in r['properties'].items())})")
    return "\n".join(lines)
