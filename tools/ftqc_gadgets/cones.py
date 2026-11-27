"""Fault cones of a block-level program over a transversal library.

A program is a list of layers of operations on code blocks (logical qubits):

    ('NEW', kind, b...)   fresh resource.  zero, plus, cat: a verified ancilla used as data;
                          std: a block just after a full error correction; and1: three std
                          blocks carrying an AND resource; any other kind is an outer resource
                          (two levels): std blocks sharing one outer atom per sector
    ('CX', c, t) ('CZ', u, v)            transversal gates
    ('CXg', c, t[, bit]) ('CZg', ...)    the same, classically controlled
    ('MZ', b[, bit]) ('MX', b[, bit])    destructive transversal measurement, decoded
    ('ECa', b) ('ECb', b)                the two halves of Steane's error correction: X
                                         errors against a |+_L> ancilla, then Z against |0_L>
    ('EC2a', b) ('EC2b', b)              outer error correction (two levels)
    ('TOF', n, b...)                     n coordinate-local non-Clifford locations on the blocks
    ('DROP', b)
Blocks not named in a layer idle for that tick.

A *source* is something that can put an error on one coordinate: a location class (one
scheduled operation, standing for one physical location per coordinate) or an atom (the
residual error of a delivered ancilla, per sector).  A *decoding event* is a place where a
bounded-distance decoder is applied.  The cone of an event is the set of sources whose
faults can be part of the pattern its decoder sees.

Three things are computed, by three different routes:

Flow     the cones, by the propagation rules (unions along gates, error correction as a
         sink), written as a graph of values: every operation makes new values from old
         ones and attaches its sources.  With an outer level, the same for outer cones,
         whose sources are the inner events.
exact    for a chosen setting of the classical controls, the cones by linear algebra over
         GF(2): every source contributes symbolic generators, gates add them, and a
         correction adds back exactly what its decoder saw.  Nothing is assumed about
         error correction being a sink; it comes out of the algebra.  The result must be
         contained in the Flow cones.
Direct   the program built qubit by qubit on a code, faults propagated through the real
         circuit, real lookup decoders, corrections carried through the later gates: the
         first decoding that goes wrong must have more than t dirty coordinates in the
         cone Flow gave it.

The statement these support (THEOREM) is proved by induction over the events; the tool
checks its hypotheses on the program and tests its conclusion on physical circuits.
"""
import hashlib
import json
import random
from itertools import combinations

from ftqc_audit.circuit import Circuit
from ftqc_audit.model import Model, propagate
from ftqc_codes import gf2, hier

from .gadget import GadgetError, prep_layers

SCHEMA = "ftqc_gadgets-cones-1"

THEOREM = (
    "Every operation of the program is transversal, so a fault at coordinate c of a location "
    "class, or an error at coordinate c of an ancilla, stays at coordinate c of every block it "
    "reaches.  For each block and sector keep the set of sources that can have reached it: a "
    "location class joins both sectors of the blocks it acts on; CX carries X sources from control "
    "to target and Z sources from target to control; CZ carries X sources of each block to the Z "
    "sector of the other; a classically controlled gate is treated as acting, and as two idle "
    "locations when it does not; a non-Clifford location pools everything on its blocks.  In "
    "Steane's correction of X errors the data's pattern is copied onto a fresh ancilla and "
    "measured; if the measured pattern is decoded exactly, applying it to the data leaves only "
    "what the ancilla, the gate and the measurement contributed, so the step is a sink for the "
    "X sector, and the ancilla's Z errors flow to the data; the Z half is the mirror image.  The "
    "cone of a decoding event is the set of sources of the pattern its decoder sees.  By "
    "induction over the events in time: if every event has at most t coordinates with a faulty "
    "source in its cone, every decoder returns the pattern it faces, every decoded bit and so "
    "every classical control is the ideal one, and each output block ends within t errors per "
    "sector of the ideal state.  One level up, a failed inner event is a fault of the outer "
    "program in one column, and the same statement holds with columns for coordinates.")

ARITY = {"CX": 2, "CXg": 2, "CZ": 2, "CZg": 2, "MZ": 1, "MX": 1, "DROP": 1,
         "ECa": 1, "ECb": 1, "EC2a": 1, "EC2b": 1}
INNER_KINDS = ("zero", "plus", "cat", "std", "and1")


def _bits(mask):
    """Positions of the set bits of a large, sparse mask."""
    out = []
    while mask:
        low = mask & -mask
        out.append(low.bit_length() - 1)
        mask ^= low
    return out


class ConeError(GadgetError):
    pass


def _name(b):
    return tuple(b) if isinstance(b, list) else b


class Program:
    def __init__(self, layers, inputs, outputs, name="program"):
        self.layers = [[tuple(_name(x) for x in op) for op in layer] for layer in layers]
        self.inputs, self.outputs, self.name = list(inputs), list(outputs), name

    def to_json(self):
        return dict(name=self.name, inputs=self.inputs, outputs=self.outputs,
                    layers=[[list(op) for op in layer] for layer in self.layers], sha256=self.sha256())

    def sha256(self):
        body = dict(inputs=self.inputs, outputs=self.outputs, layers=[[list(op) for op in l] for l in self.layers])
        return hashlib.sha256(json.dumps(body, sort_keys=True, default=str).encode()).hexdigest()

    @staticmethod
    def from_json(d):
        p = Program(d["layers"], d["inputs"], d["outputs"], d.get("name", "program"))
        if "sha256" in d and d["sha256"] != p.sha256():
            raise ConeError("sha256 in program file does not match its contents")
        return p

    def blocks_of(self, op):
        k = op[0]
        if k == "NEW":
            return list(op[2:])
        if k == "TOF":
            return list(op[2:])
        if k not in ARITY:
            raise ConeError(f"unknown operation {k!r}")
        if len(op) < 1 + ARITY[k]:
            raise ConeError(f"{k} needs {ARITY[k]} blocks")
        return list(op[1:1 + ARITY[k]])

    def two_level(self):
        return any(op[0] in ("EC2a", "EC2b") or (op[0] == "NEW" and op[1] not in INNER_KINDS)
                   for layer in self.layers for op in layer)


# ---------------------------------------------------------------------- the cones by rule
class Flow:
    """sources: list of dicts (level 1: kind class | atom; level 2: kind event | atom).
    events: inner decoding events, in order; outer: outer decoding events, in order.
    Cones are bit masks over source numbers: cone(event) -> int."""

    def __init__(self, program, two=None):
        self.p = program
        self.two = program.two_level() if two is None else two
        self.sources, self.cone_of = [], []          # node -> mask
        self.events, self.outer, self.terms = [], [], []
        self.S, self.O, self.pend = {}, {}, {}
        self.known, self.tick = set(), 0
        self._ec2_atom = {}
        for b in program.inputs:
            self._blank(b)
        for layer in program.layers:
            self._layer(layer)
            self.tick += 1
        self._finish()

    # -- graph
    def _src(self, **kw):
        kw["tick"] = self.tick
        self.sources.append(kw)
        return len(self.sources) - 1

    def _cls(self, what, blocks, weight=1):
        return self._src(level=1, kind="class", what=what, blocks=list(blocks), weight=weight)

    def _atom(self, anc, sector, block, level=1, origin="new"):
        return self._src(level=level, kind="atom", what=anc, sector=sector, blocks=[block], weight=0,
                         origin=origin)

    def _node(self, preds, srcs=()):
        m = 0
        for p in preds:
            m |= self.cone_of[p]
        for s in srcs:
            m |= 1 << s
        self.cone_of.append(m)
        return len(self.cone_of) - 1

    def _blank(self, b):
        if b in self.S:
            raise ConeError(f"block {b} is already live")
        self.S[b] = {"X": self._node(()), "Z": self._node(())}
        if self.two:
            self.O[b] = {"X": self._node(()), "Z": self._node(())}

    def _event(self, kind, block, sector, node, tag):
        self.events.append(dict(kind=kind, block=block, sector=sector, tick=self.tick, node=node, tag=tag))
        i = len(self.events) - 1
        if not self.two:
            return None
        return self._src(level=2, kind="event", what=tag, event=i, blocks=[block], weight=1)

    def _outer(self, block, sector, node, tag):
        self.outer.append(dict(block=block, sector=sector, tick=self.tick, node=node, tag=tag))

    # -- pieces
    def _post_ec(self, b):
        """A block just after a full error correction: what the two halves leave on it."""
        self._blank(b)
        c1, m1, c2, m2 = (self._cls(("std", w), [b]) for w in ("c1", "m1", "c2", "m2"))
        ax = [self._atom("plus", "X", b, origin="std"), self._atom("zero", "X", b, origin="std")]
        az = [self._atom("zero", "Z", b, origin="std")]
        self.S[b] = {"X": self._node((), [c1, m1, c2] + ax), "Z": self._node((), [c2, m2] + az)}

    def _new(self, kind, names):
        if kind in ("zero", "plus", "cat"):
            for b in names:
                self._blank(b)
                self.S[b] = {s: self._node((), [self._atom(kind, s, b)]) for s in "XZ"}
                e = self._event("junk", b, "XZ", None, "junk:" + kind)
                if self.two:
                    self.O[b] = {s: self._node((), [e]) for s in "XZ"}
        elif kind == "std":
            for b in names:
                self._post_ec(b)
        elif kind == "and1":
            for b in names:
                self._post_ec(b)
            e = self._event("junk", names[0], "XZ", None, "and1")
            if self.two:
                for b in names:
                    self.O[b] = {s: self._node((), [e]) for s in "XZ"}
        else:
            if not self.two:
                raise ConeError(f"outer resource {kind!r} in a single-level program")
            self._new_outer(kind, names)

    def _new_outer(self, kind, names):
        at = {s: self._atom(kind, s, names[0], level=2, origin="outer") for s in "XZ"}
        for b in names:
            self._post_ec(b)
            self.O[b] = {s: self._node((), [at[s]]) for s in "XZ"}
        self.terms.append(dict(tick=self.tick, tag="junk:" + kind, blocks=list(names)))
        return at

    def _gate(self, c, t, kind, weight):
        g = self._cls(("gate", kind), [c, t], weight)
        for store, srcs in ((self.S, [g]), (self.O, [])) if self.two else ((self.S, [g]),):
            a, b = store[c], store[t]
            if kind == "CX":
                new = {(t, "X"): [b["X"], a["X"]], (c, "Z"): [a["Z"], b["Z"]],
                       (c, "X"): [a["X"]], (t, "Z"): [b["Z"]]}
            else:
                new = {(c, "Z"): [a["Z"], b["X"]], (t, "Z"): [b["Z"], a["X"]],
                       (c, "X"): [a["X"]], (t, "X"): [b["X"]]}
            made = {k: self._node(p, srcs) for k, p in new.items()}
            for (blk, s), n in made.items():
                store[blk][s] = n

    def _measure(self, b, sector, tag, keep=False):
        m = self._cls(("meas", sector), [b])
        e = self._event("meas", b, sector, self._node([self.S[b][sector]], [m]), tag)
        if self.two and not keep:
            self._outer(b, sector, self._node([self.O[b][sector]], [e]), tag)
        return e

    def _ec_a(self, b):
        c1, m1 = self._cls(("eca", "c1"), [b]), self._cls(("eca", "m1"), [b])
        ax, az = self._atom("plus", "X", b, origin="eca"), self._atom("plus", "Z", b, origin="eca")
        S = self.S[b]
        e = self._event("ec", b, "X", self._node([S["X"]], [c1, m1, ax]), "ec")
        self.S[b] = {"X": self._node((), [c1, m1, ax]), "Z": self._node([S["Z"]], [c1, az])}
        if self.two:
            self.O[b]["X"] = self._node([self.O[b]["X"]], [e])
        self.pend[b] = "ECb"

    def _ec_b(self, b):
        c2, m2 = self._cls(("ecb", "c2"), [b]), self._cls(("ecb", "m2"), [b])
        ax, az = self._atom("zero", "X", b, origin="ecb"), self._atom("zero", "Z", b, origin="ecb")
        S = self.S[b]
        e = self._event("ec", b, "Z", self._node([S["Z"]], [c2, m2, az]), "ec")
        self.S[b] = {"X": self._node([S["X"]], [c2, ax]), "Z": self._node((), [c2, m2, az])}
        if self.two:
            self.O[b]["Z"] = self._node([self.O[b]["Z"]], [e])
        del self.pend[b]

    def _ec2_a(self, b):
        A = ("~p", b)
        self._ec2_atom[b] = self._new_outer("plus2", [A])
        self._gate(b, A, "CX", 1)
        self.pend[b] = "EC2b"

    def _ec2_b(self, b):
        A = ("~p", b)
        s1 = self._measure(A, "X", "ec2", keep=True)
        self._outer(A, "X", self._node([self.O[A]["X"]], [s1]), "ec2")
        # the data now carries what this ancilla, the gate and the measurement contributed
        self.O[b]["X"] = self._node((), [s1, self._ec2_atom[b]["X"]])
        del self.S[A], self.O[A]
        Z0 = ("~z", b)
        at = self._new_outer("zero2", [Z0])
        self._gate(Z0, b, "CX", 1)
        s2 = self._measure(Z0, "Z", "ec2", keep=True)
        self._outer(Z0, "Z", self._node([self.O[Z0]["Z"]], [s2]), "ec2")
        self.O[b]["Z"] = self._node((), [s2, at["Z"]])
        del self.S[Z0], self.O[Z0], self.pend[b]

    def _tof(self, n, blocks):
        ids = [self._cls(("tof", i), blocks) for i in range(n)]
        for store, srcs in ((self.S, ids), (self.O, [])) if self.two else ((self.S, ids),):
            pool = [store[b][s] for b in blocks for s in "XZ"]
            node = {(b, s): self._node(pool, srcs) for b in blocks for s in "XZ"}
            for (b, s), v in node.items():
                store[b][s] = v

    # -- layers
    def _layer(self, layer):
        used, newbits = set(), []
        for op in layer:
            k = op[0]
            if k == "NEW":
                self._new(op[1], list(op[2:]))
                continue
            bs = self.p.blocks_of(op)
            for b in bs:
                if b in used:
                    raise ConeError(f"tick {self.tick}: block {b} used twice")
                if b not in self.S:
                    raise ConeError(f"tick {self.tick}: block {b} is not live")
                if b in self.pend and k != self.pend[b]:
                    raise ConeError(f"tick {self.tick}: block {b} is in the middle of an error correction")
                used.add(b)
            if k in ("CXg", "CZg") and len(op) > 3 and op[3] not in self.known:
                raise ConeError(f"tick {self.tick}: control {op[3]} is not measured yet")
            if k in ("MZ", "MX") and len(op) > 2:
                if op[2] in self.known or op[2] in newbits:
                    raise ConeError(f"tick {self.tick}: classical record {op[2]} reused")
                newbits.append(op[2])
            if k in ("CX", "CXg", "CZ", "CZg"):
                self._gate(op[1], op[2], k[:2], 2 if k.endswith("g") else 1)
            elif k in ("MZ", "MX"):
                self._measure(op[1], "X" if k == "MZ" else "Z", "meas")
                del self.S[op[1]]
                self.O.pop(op[1], None)
            elif k == "DROP":
                del self.S[op[1]]
                self.O.pop(op[1], None)
            elif k == "ECa":
                self._ec_a(op[1])
            elif k == "ECb":
                if self.pend.get(op[1]) != "ECb":
                    raise ConeError(f"tick {self.tick}: ECb on {op[1]} without its first half")
                self._ec_b(op[1])
            elif k == "EC2a":
                self._ec2_a(op[1])
            elif k == "EC2b":
                if self.pend.get(op[1]) != "EC2b":
                    raise ConeError(f"tick {self.tick}: EC2b on {op[1]} without its first half")
                self._ec2_b(op[1])
            elif k == "TOF":
                self._tof(int(op[1]), bs)
        for b in self.pend:
            if b not in used:
                raise ConeError(f"tick {self.tick}: error correction on {b} interrupted")
        self.known.update(newbits)
        for b in list(self.S):
            if b not in used and not (isinstance(b, tuple) and b and b[0] in ("~p", "~z")):
                g = self._cls(("idle",), [b])
                self.S[b] = {s: self._node([self.S[b][s]], [g]) for s in "XZ"}

    def _finish(self):
        if set(self.p.outputs) != set(self.S):
            raise ConeError("the live blocks at the end are not the declared outputs")
        if self.pend:
            raise ConeError("an error correction is unfinished at the end")
        for b in self.p.outputs:
            for s in "XZ":
                e = self._event("out", b, s, self._node([self.S[b][s]]), "out")
                if self.two:
                    self._outer(b, s, self._node([self.O[b][s]], [e]), "out")

    # -- reading
    def cone(self, event):
        return 0 if event["node"] is None else self.cone_of[event["node"]]

    def members(self, mask):
        return _bits(mask)

    def weight(self, mask):
        """(location classes counted with their weights, {(ancilla, sector): atoms})."""
        L, atoms = 0, {}
        for i in _bits(mask):
            s = self.sources[i]
            if s["kind"] == "atom":
                key = (s["what"], s["sector"])
                atoms[key] = atoms.get(key, 0) + 1
            else:
                L += s["weight"]
        return L, atoms

    def key(self, i):
        """A name for source i that does not depend on numbering."""
        s = self.sources[i]
        return (s["level"], s["kind"], s["tick"], repr(s["what"]), repr(s["blocks"]), s.get("sector"), s.get("event"))


# ---------------------------------------------------------------------- the cones by algebra
def exact(program, flow, acts=None):
    """Inner cones by linear propagation of symbolic generators, for one setting of the
    classical controls: acts(tick, op) -> bool says whether a controlled gate is applied
    (default: all are).  Returns a list of source masks, one per inner decoding event of
    `flow`, in the same order (None for events without a cone).  Sources are numbered as
    in `flow`: the walk creates them in the same order and checks that it does."""
    src_of_gen, x, z = [], {}, {}
    out, state = [], dict(n=0, tick=0)
    level1 = [i for i, s in enumerate(flow.sources) if s["level"] == 1]

    def take(kind, what=None):
        i = level1[state["n"]]
        s = flow.sources[i]
        if s["kind"] != kind or s["tick"] != state["tick"] or (what is not None and s["what"] != what):
            raise ConeError("exact walk and flow disagree on the order of sources")
        state["n"] += 1
        return i

    def gen(src):
        src_of_gen.append(src)
        return 1 << (len(src_of_gen) - 1)

    def emit(mask):
        m = 0
        for g in _bits(mask):
            m |= 1 << src_of_gen[g]
        out.append(m)

    def post_ec(b):
        c = [take("class", ("std", w)) for w in ("c1", "m1", "c2", "m2")]
        ap, a0x, a0z = take("atom", "plus"), take("atom", "zero"), take("atom", "zero")
        x[b] = gen(c[0]) ^ gen(c[1]) ^ gen(c[2]) ^ gen(ap) ^ gen(a0x)
        z[b] = gen(c[2]) ^ gen(c[3]) ^ gen(a0z)

    def gate(c, t, kind, on):
        g = take("class", ("gate", kind))
        if on:
            if kind == "CX":
                x[t] ^= x[c]
                z[c] ^= z[t]
            else:
                zc, zt = z[c] ^ x[t], z[t] ^ x[c]
                z[c], z[t] = zc, zt
        for b in (c, t):
            x[b] ^= gen(g)
            z[b] ^= gen(g)

    def measure(b, sector):
        m = take("class", ("meas", sector))
        emit((x[b] if sector == "X" else z[b]) ^ gen(m))

    def new_outer(names):
        for s in "XZ":
            pass
        for b in names:
            post_ec(b)

    for b in program.inputs:
        x[b] = z[b] = 0
    pend_anc = {}
    for tick, layer in enumerate(program.layers):
        state["tick"] = tick
        used = set()
        for op in layer:
            k = op[0]
            if k == "NEW":
                kind, names = op[1], list(op[2:])
                if kind in ("zero", "plus", "cat"):
                    for b in names:
                        x[b], z[b] = gen(take("atom", kind)), gen(take("atom", kind))
                        out.append(None)
                elif kind == "std":
                    for b in names:
                        post_ec(b)
                elif kind == "and1":
                    for b in names:
                        post_ec(b)
                    out.append(None)
                else:
                    new_outer(names)
                continue
            bs = program.blocks_of(op)
            used.update(bs)
            if k in ("CX", "CZ"):
                gate(op[1], op[2], k, True)
            elif k in ("CXg", "CZg"):
                gate(op[1], op[2], k[:2], True if acts is None else bool(acts(tick, op)))
            elif k in ("MZ", "MX"):
                measure(op[1], "X" if k == "MZ" else "Z")
                del x[op[1]], z[op[1]]
            elif k == "DROP":
                del x[op[1]], z[op[1]]
            elif k == "ECa":
                b = op[1]
                c1, m1 = take("class", ("eca", "c1")), take("class", ("eca", "m1"))
                ax, az = take("atom", "plus"), take("atom", "plus")
                xa = gen(ax) ^ x[b]                 # CX data -> ancilla copies X onto the ancilla
                z[b] ^= gen(az)                     # and the ancilla's Z onto the data
                x[b] ^= gen(c1)
                z[b] ^= gen(c1)
                v = xa ^ gen(c1) ^ gen(m1)
                emit(v)
                x[b] ^= v                           # the decoded pattern applied to the data
            elif k == "ECb":
                b = op[1]
                c2, m2 = take("class", ("ecb", "c2")), take("class", ("ecb", "m2"))
                ax, az = take("atom", "zero"), take("atom", "zero")
                x[b] ^= gen(ax)                     # CX ancilla -> data copies the ancilla's X
                za = gen(az) ^ z[b]                 # and the data's Z onto the ancilla
                x[b] ^= gen(c2)
                z[b] ^= gen(c2)
                v = za ^ gen(c2) ^ gen(m2)
                emit(v)
                z[b] ^= v
            elif k == "EC2a":
                A = ("~p", op[1])
                new_outer([A])
                gate(op[1], A, "CX", True)
            elif k == "EC2b":
                A, Z0 = ("~p", op[1]), ("~z", op[1])
                measure(A, "X")
                del x[A], z[A]
                new_outer([Z0])
                gate(Z0, op[1], "CX", True)
                measure(Z0, "Z")
                del x[Z0], z[Z0]
            elif k == "TOF":
                ids = [take("class", ("tof", i)) for i in range(int(op[1]))]
                pool = 0
                for b in bs:
                    pool |= x[b] | z[b]
                for i in ids:
                    pool |= gen(i)
                for b in bs:
                    x[b] = z[b] = pool
        for b in list(x):
            if b not in used and not (isinstance(b, tuple) and b and b[0] in ("~p", "~z")):
                g = take("class", ("idle",))
                x[b] ^= gen(g)
                z[b] ^= gen(g)
    for b in program.outputs:
        emit(x[b])
        emit(z[b])
    if state["n"] != len(level1) or len(out) != len(flow.events):
        raise ConeError("exact walk and flow disagree on the number of sources or events")
    return out


def contained(program, flow, samples=8, seed=1):
    """exact cones, for all controls on, all off and `samples` random settings, against
    the flow's.  Returns dict(settings, events, outside=[...first few...], strictly_smaller)."""
    rnd = random.Random(seed)
    settings = [("all on", None), ("all off", lambda t, op: False)]
    for i in range(samples):
        memo = {}
        settings.append((f"random {i}", lambda t, op, memo=memo: memo.setdefault((t, op), rnd.random() < 0.5)))
    outside, smaller = [], 0
    union = [0] * len(flow.events)
    for name, acts in settings:
        ex = exact(program, flow, acts)
        for i, (m, ev) in enumerate(zip(ex, flow.events)):
            if m is None:
                continue
            union[i] |= m
            extra = m & ~flow.cone(ev)
            if extra and len(outside) < 5:
                outside.append(dict(setting=name, event=i, kind=ev["kind"], block=str(ev["block"]),
                                    sector=ev["sector"], tick=ev["tick"],
                                    sources=[flow.key(j) for j in _bits(extra)][:5]))
    for m, ev in zip(union, flow.events):
        if ev["node"] is not None and m != flow.cone(ev):
            smaller += 1
    return dict(settings=len(settings), events=sum(1 for e in flow.events if e["node"] is not None),
                outside=outside, events_where_rule_exceeds_algebra=smaller)


# ---------------------------------------------------------------------- the physical test
class Direct:
    """The program on a CSS code, one physical qubit per coordinate.  Ancillas and fresh
    blocks are ideal encoded states followed by one noisy wire per qubit (their atoms).
    Classically controlled gates are applied or not according to `acts`.  Programs with
    non-Clifford locations, and1 or outer operations are refused."""

    def __init__(self, program, code, acts=None, flow=None):
        self.p, self.code, self.n = program, code, code.n
        self.flow = flow or Flow(program, two=False)
        if self.flow.two:
            raise ConeError("the physical test is for single-level programs")
        n = self.n
        ticks, tags, nq = [], [], [0]
        level1 = list(range(len(self.flow.sources)))
        cursor = [0]

        def take():
            cursor[0] += 1
            return level1[cursor[0] - 1]

        def alloc():
            qs = list(range(nq[0], nq[0] + n))
            nq[0] += n
            return qs

        def tick(ops, tag=None):
            """ops: list of op tuples; tag: for noisy ops, function (index in ops) -> tag."""
            ticks.append(list(ops))
            tags.append([None if str(op[0]).startswith("!") else tag(i) for i, op in enumerate(ops)])

        def fresh(state):
            qs = alloc()
            for layer in prep_layers(self.code, state, qs):
                tick([("!" + op[0],) + tuple(op[1:]) for op in layer])
            return qs

        q = {}
        self.events = []                               # dict(records, checks, kind, tick, data, flow index)
        fi = [0]                                       # flow event counter

        def nxt():
            fi[0] += 1
            return fi[0] - 1

        def event(kind, recs_tick, qubits, sector, data, idx):
            self.events.append(dict(kind=kind, tick=recs_tick, qubits=list(qubits), sector=sector,
                                    data=None if data is None else list(data), flow=idx))

        for b in program.inputs:
            q[b] = fresh("zero")
        for t, layer in enumerate(program.layers):
            used = set()
            pre, gates, post = [], [], []              # deliveries, then gates and idles, then measurements
            for op in layer:
                k = op[0]
                if k == "NEW":
                    kind = op[1]
                    for b in op[2:]:
                        if kind in ("zero", "plus"):
                            q[b] = fresh(kind)
                            ax, az = take(), take()
                            pre.append((q[b], ("atom", ax, az)))
                            fi[0] += 1                 # the junk event has no cone
                        elif kind == "std":
                            q[b] = fresh("zero")
                            ids = [take() for _ in range(7)]      # c1 m1 c2 m2, plus_X zero_X, zero_Z
                            xs = [ids[0], ids[1], ids[2], ids[4], ids[5]]
                            zs = [ids[2], ids[3], ids[6]]
                            pre.append((q[b], ("multi", xs, zs)))
                        else:
                            raise ConeError(f"the physical test does not build resource {kind!r}")
                    continue
                bs = program.blocks_of(op)
                used.update(bs)
                if k in ("CX", "CXg", "CZ", "CZg"):
                    g = take()
                    on = True if not k.endswith("g") or acts is None else bool(acts(t, op))
                    gates.append((k[:2], q[op[1]], q[op[2]], g, on))
                elif k in ("MZ", "MX"):
                    post.append((k, q.pop(op[1]), take(), op[1], nxt()))
                elif k == "DROP":
                    q.pop(op[1])
                elif k == "ECa":
                    c1, m1, ax, az = take(), take(), take(), take()
                    A = fresh("plus")
                    pre.append((A, ("atom", ax, az)))
                    gates.append(("CX", q[op[1]], A, c1, True))
                    post.append(("ECa", A, m1, op[1], nxt()))
                elif k == "ECb":
                    c2, m2, ax, az = take(), take(), take(), take()
                    A = fresh("zero")
                    pre.append((A, ("atom", ax, az)))
                    gates.append(("CX", A, q[op[1]], c2, True))
                    post.append(("ECb", A, m2, op[1], nxt()))
                else:
                    raise ConeError(f"the physical test does not build {k}")
            if pre:
                ops, tg = [], []
                for qs, tag in pre:
                    for c, w in enumerate(qs):
                        ops.append(("I", w))
                        tg.append(tag + (c,))
                tick(ops, lambda i, tg=tg: tg[i])
            ops, tg, hs = [], [], []
            for kind, A, B, g, on in gates:
                for c in range(n):
                    if on:
                        ops.append(("CX", A[c], B[c]))
                        tg.append(("class", g, c))
                    else:
                        ops += [("I", A[c]), ("I", B[c])]
                        tg += [("class", g, c), ("class", g, c)]
                if kind == "CZ" and on:
                    hs += B
            idle = [b for b in q if b not in used]
            for b in idle:
                g = take()
                for c, w in enumerate(q[b]):
                    ops.append(("I", w))
                    tg.append(("class", g, c))
            if hs:
                tick([("!H", w) for w in hs])
            if ops:
                tick(ops, lambda i, tg=tg: tg[i])
            if hs:
                tick([("!H", w) for w in hs])
            if post:
                ops, tg = [], []
                for k, qs, m, b, idx in post:
                    for c, w in enumerate(qs):
                        ops.append(("M" if k in ("MZ", "ECa") else "MX", w))
                        tg.append(("class", m, c))
                tick(ops, lambda i, tg=tg: tg[i])
                for k, qs, m, b, idx in post:
                    sector = "X" if k in ("MZ", "ECa") else "Z"
                    event("meas" if k in ("MZ", "MX") else "ec", len(ticks) - 1, qs, sector,
                          q[b] if k in ("ECa", "ECb") else None, idx)
        if cursor[0] != len(level1):
            raise ConeError("the physical build and the flow disagree on the number of sources")
        self.out_blocks = [(b, q[b]) for b in program.outputs]
        for b, qs in self.out_blocks:
            for s in "XZ":
                event("out", None, qs, s, None, nxt())
        if fi[0] != len(self.flow.events):
            raise ConeError("the physical build and the flow disagree on the number of events")
        outs = [w for _, qs in self.out_blocks for w in qs]
        self.circuit = Circuit(ticks, nq[0], detectors=[], outputs=outs, name=program.name)
        self.model = Model(self.circuit)
        flat = [tag for row in tags for tag in row if tag is not None]
        if len(flat) != len(self.model.locations):
            raise ConeError("internal: tags and locations differ")
        self.tags = flat
        rec_of = {m[2]: i for i, m in enumerate(self.circuit.meas)}
        for ev in self.events:
            ev["records"] = None if ev["kind"] == "out" else [rec_of[w] for w in ev["qubits"]]
        self.tab = {"X": hier.leader_table(code.hz, n)[0], "Z": hier.leader_table(code.hx, n)[0]}
        self.radius = {"X": hier.leader_table(code.hz, n)[1], "Z": hier.leader_table(code.hx, n)[1]}
        self.stab = {"X": gf2.Echelon(code.hx), "Z": gf2.Echelon(code.hz)}
        self.checks = {"X": list(code.hz), "Z": list(code.hx)}
        self.t = min(self.radius.values())
        self._corr = {}
        self.faults = []                               # (location, label, (x, z, flips))
        for loc in self.model.locations:
            for label, eff in loc.faults.items():
                self.faults.append((loc.index, label, eff))

    def dirty(self, items):
        """{source: set of coordinates} for a list of faults."""
        d = {}
        for loc, label, _ in items:
            tag = self.tags[loc]
            if tag[0] == "class":
                d.setdefault(tag[1], set()).add(tag[2])
            else:
                hx, hz = label in ("X", "Y"), label in ("Z", "Y")
                if tag[0] == "atom":
                    xs, zs = [tag[1]], [tag[2]]
                else:
                    xs, zs = tag[1], tag[2]
                for s in (xs if hx else []) + (zs if hz else []):
                    d.setdefault(s, set()).add(tag[-1])
        return d

    def cone_dirt(self, event, dirt):
        """Number of coordinates with a faulty source in the event's cone."""
        mask = self.flow.cone(self.flow.events[event["flow"]])
        cs = set()
        for s, where in dirt.items():
            if (mask >> s) & 1:
                cs |= where
        return len(cs)

    def _correction(self, k, p):
        key = (k, p)
        if key not in self._corr:
            ev = self.events[k]
            w = 1 << ev["data"][p]
            self._corr[key] = propagate(self.circuit, ev["tick"], w if ev["sector"] == "X" else 0,
                                        0 if ev["sector"] == "X" else w)
        return self._corr[key]

    def run(self, items):
        """Decode the events in time order with corrections carried forward.  Returns
        (index of the first event decoded wrongly or None, its number of dirty cone
        coordinates, largest number of dirty coordinates over all cones)."""
        x = z = fl = 0
        for _, _, (fx, fz, ff) in items:
            x ^= fx
            z ^= fz
            fl ^= ff
        dirt = self.dirty(items)
        worst = max((self.cone_dirt(ev, dirt) for ev in self.events), default=0)
        for k, ev in enumerate(self.events):
            s = ev["sector"]
            if ev["kind"] == "out":
                frame = x if s == "X" else z
                e = sum(1 << c for c, w in enumerate(ev["qubits"]) if (frame >> w) & 1)
            else:
                e = sum(1 << c for c, r in enumerate(ev["records"]) if (fl >> r) & 1)
            sy = 0
            for i, row in enumerate(self.checks[s]):
                if (row & e).bit_count() & 1:
                    sy |= 1 << i
            pat = self.tab[s].get(sy)
            if pat is None or (pat ^ e) not in self.stab[s]:
                return k, self.cone_dirt(ev, dirt), worst
            if ev["kind"] == "ec":
                for c in _bits(pat):
                    cx, cz, cf = self._correction(k, c)
                    x ^= cx
                    z ^= cz
                    fl ^= cf
        return None, None, worst

    def test(self, exhaustive=1, samples=4000, seed=1, search=4000):
        """Every set of at most `exhaustive` faults; `samples` random sets of three kinds
        (a few faults anywhere; many faults inside at most t coordinates; random sets kept
        when every cone has at most t dirty coordinates); then a search for a failing set
        among those with t + 1 dirty coordinates in some cone.  A violation is a wrong
        decoding whose cone has at most t dirty coordinates."""
        rnd = random.Random(seed)
        t = self.t
        out = dict(t=t, qubits=self.circuit.n, locations=len(self.model.locations), events=len(self.events),
                   fault_labels=len(self.faults), sets=0, within_hypothesis=0, largest_set=0,
                   failures_beyond=0, violations=0, first_violation=None, beyond=None)

        def judge(items):
            if len({f[0] for f in items}) != len(items):
                return None
            bad, d, worst = self.run(items)
            out["sets"] += 1
            if worst <= t:
                out["within_hypothesis"] += 1
                out["largest_set"] = max(out["largest_set"], len(items))
            if bad is not None:
                if d <= t:
                    out["violations"] += 1
                    if out["first_violation"] is None:
                        out["first_violation"] = dict(event=bad, dirty=d, faults=[[f[0], f[1]] for f in items])
                else:
                    out["failures_beyond"] += 1
            return bad, d, worst

        for w in range(exhaustive + 1):
            for items in combinations(self.faults, w):
                judge(items)
        by_coord = {}
        for f in self.faults:
            by_coord.setdefault(self.tags[f[0]][-1], []).append(f)
        for i in range(samples):
            mode = i % 3
            if mode == 0:
                items = rnd.sample(self.faults, rnd.randint(2, t + 3))
            elif mode == 1:
                items = []
                for c in rnd.sample(sorted(by_coord), rnd.randint(1, t)):
                    items += rnd.sample(by_coord[c], min(len(by_coord[c]), rnd.randint(1, 12)))
            else:
                items = rnd.sample(self.faults, rnd.randint(t + 1, 3 * t + 4))
                if len({f[0] for f in items}) == len(items):
                    dirt = self.dirty(items)
                    if max(self.cone_dirt(ev, dirt) for ev in self.events) > t:
                        continue
            judge(items)
        for tries in range(1, search + 1):
            cs = rnd.sample(sorted(by_coord), t + 1)
            items = [rnd.choice(by_coord[c]) for c in cs]
            r = judge(items)
            if r and r[0] is not None and r[1] == t + 1:
                out["beyond"] = dict(dirty=r[1], event=r[0], tries=tries, faults=[[f[0], f[1]] for f in items])
                break
        return out


# ---------------------------------------------------------------------- receipts
def summarize(flow):
    """Per event: what a probability bound needs.  Inner events: location weight and atoms
    by kind; outer events: the inner events and outer atoms in the cone."""
    inner, shapes = [], {}
    for i, ev in enumerate(flow.events):
        if ev["node"] is None:
            inner.append(dict(kind=ev["kind"], tag=ev["tag"], block=str(ev["block"]), tick=ev["tick"]))
            continue
        L, atoms = flow.weight(flow.cone(ev))
        row = dict(kind=ev["kind"], block=str(ev["block"]), sector=ev["sector"], tick=ev["tick"], weight=L,
                   atoms=sorted([k[0], k[1], v] for k, v in atoms.items()))
        inner.append(row)
        key = json.dumps([L, row["atoms"]])
        shapes[key] = shapes.get(key, 0) + 1
    outer = []
    ev_of = {i: s["event"] for i, s in enumerate(flow.sources) if s["level"] == 2 and s["kind"] == "event"}
    for ev in flow.outer:
        mem = _bits(flow.cone_of[ev["node"]])
        atoms = {}
        for j in mem:
            s = flow.sources[j]
            if s["kind"] == "atom":
                atoms[(s["what"], s["sector"])] = atoms.get((s["what"], s["sector"]), 0) + 1
        outer.append(dict(block=str(ev["block"]), sector=ev["sector"], tick=ev["tick"], tag=ev["tag"],
                          inner_events=sorted(ev_of[j] for j in mem if j in ev_of),
                          atoms=sorted([k[0], k[1], v] for k, v in atoms.items())))
    cones = [r for r in inner if "weight" in r]
    return dict(inner_events=inner, outer_events=outer, top_level_terms=flow.terms,
                stats=dict(layers=len(flow.p.layers), sources=len(flow.sources),
                           inner_decoding_events=len(cones), other_inner_events=len(inner) - len(cones),
                           outer_decoding_events=len(outer),
                           largest_weight=max((r["weight"] for r in cones), default=0),
                           largest_atoms=max((sum(a[2] for a in r["atoms"]) for r in cones), default=0),
                           distinct_cone_shapes=len(shapes),
                           largest_outer_cone=max((len(r["inner_events"]) for r in outer), default=0)))


def build(program, samples=8, seed=1, events=True):
    """The cones receipt.  events=False keeps only the statistics and a hash of the event
    list (large programs)."""
    flow = Flow(program)
    summ = summarize(flow)
    chk = contained(program, flow, samples, seed)
    digest = hashlib.sha256(json.dumps([summ["inner_events"], summ["outer_events"]], sort_keys=True).encode()).hexdigest()
    rec = dict(schema=SCHEMA, tool="ftqc_gadgets", program=program.to_json(), two_levels=flow.two,
               stats=summ["stats"], events_sha256=digest, top_level_terms=summ["top_level_terms"],
               algebra=dict(chk, samples=samples, seed=seed), theorem=THEOREM)
    if events:
        rec["inner_events"], rec["outer_events"] = summ["inner_events"], summ["outer_events"]
    return json.loads(json.dumps(rec, default=str))


def verify(receipt):
    """Recompute everything from the program in the receipt.  Raises ConeError."""
    if receipt.get("schema") != SCHEMA:
        raise ConeError(f"not a {SCHEMA} receipt")
    prog = Program.from_json(receipt["program"])
    a = receipt.get("algebra", {})
    again = build(prog, a.get("samples", 8), a.get("seed", 1), events="inner_events" in receipt)
    for k in again:
        if again[k] != receipt.get(k):
            raise ConeError(f"cones: {k} does not recompute")
    if receipt["algebra"]["outside"]:
        raise ConeError("cones: linear propagation reaches an event from outside its cone")
    return ["program hash", "cones recomputed by rule", "linear propagation contained in the cones"]


def summary(receipt):
    s = receipt["stats"]
    lines = [f"{receipt['program']['name']}: {s['layers']} layers, {s['sources']} sources, "
             f"{s['inner_decoding_events']} inner decoding events"
             + (f", {s['outer_decoding_events']} outer" if receipt["two_levels"] else ""),
             f"  largest cone: weight {s['largest_weight']} with {s['largest_atoms']} atoms; "
             f"{s['distinct_cone_shapes']} distinct cone shapes"
             + (f"; largest outer cone {s['largest_outer_cone']} inner events" if receipt["two_levels"] else "")]
    a = receipt["algebra"]
    lines.append(f"  linear propagation under {a['settings']} settings of the controls: "
                 + ("contained in every cone" if not a["outside"] else f"OUTSIDE a cone: {a['outside'][0]}")
                 + f"; the rule exceeds the algebra at {a['events_where_rule_exceeds_algebra']} events")
    return "\n".join(lines)
