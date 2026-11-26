"""Gadgets: Clifford circuits with input blocks, output blocks and decoders.

A gadget is compiled to a ftqc_audit circuit in which nothing is left open:

  input block   Each input block arrives as one half of a noiseless encoded Bell
                pair whose other half is k reference qubits that are never touched
                again (the Choi state of the block).  What the gadget does to an
                unknown codeword is then a property of a state: the joint state of
                references and outputs.  A bare block (no code) gets one reference
                per qubit.
  input errors  One noisy idle on every input qubit, in one tick just before the
                body.  Faults there are the errors the input arrives with; they are
                counted as input weight, never as gadget faults.
  decoder       A set of measured bits, parity checks on them, and what the decoded
                error is applied to: the measurement record (the believed outcomes
                are the reported ones with the decoded error removed) or a Pauli
                frame on output qubits.
  detectors     Parities that are postselected on.  A gadget with decoders must
                declare them; without decoders they may be left to discovery.
"""
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from ftqc_audit.circuit import Circuit, CircuitError  # noqa: E402
from ftqc_codes import gf2  # noqa: E402


class GadgetError(ValueError):
    pass


def logical_x(code):
    """Independent X-type logical operators: a basis of ker(H_Z) modulo rowspace(H_X)."""
    e = gf2.Echelon(code.hx)
    return [v for v in gf2.nullspace(code.hz, code.n) if e.add(v)]


def logical_z(code):
    e = gf2.Echelon(code.hz)
    return [v for v in gf2.nullspace(code.hx, code.n) if e.add(v)]


def _layer(edges):
    """Greedy layering of CX edges, one gate per qubit per tick."""
    busy, out = set(), []
    for c, t in edges:
        k = 0
        while (k, c) in busy or (k, t) in busy:
            k += 1
        busy.add((k, c))
        busy.add((k, t))
        while len(out) <= k:
            out.append([])
        out[k].append(("CX", c, t))
    return out


def prep_layers(code, state, qubits):
    """Layers of a non-fault-tolerant preparation of |0_L> ('zero': +1 eigenstate of
    the X stabilizers and of every Z-type operator commuting with them) or |+_L>
    ('plus', the same with X and Z exchanged) on the given qubits."""
    n = code.n
    if state == "zero":
        rows, piv = gf2.systematic(code.hx, n, list(range(n)))
        first = [("RX" if q in piv else "R", qubits[q]) for q in range(n)]
        edges = [(qubits[p], qubits[t]) for r, p in zip(rows, piv) for t in gf2.supp(r) if t != p]
    elif state == "plus":
        rows, piv = gf2.systematic(code.hz, n, list(range(n)))
        first = [("R" if q in piv else "RX", qubits[q]) for q in range(n)]
        edges = [(qubits[t], qubits[p]) for r, p in zip(rows, piv) for t in gf2.supp(r) if t != p]
    else:
        raise GadgetError(f"unknown state {state!r} (zero or plus)")
    return [first] + _layer(edges)


def state_checks(code, state, basis):
    """Parities of a transversal measurement of |0_L> or |+_L> in the given basis that
    are fixed in the noiseless state, as position masks."""
    if basis == "Z":
        return list(code.hz) + (logical_z(code) if state == "zero" else [])
    return list(code.hx) + (logical_x(code) if state == "plus" else [])


class Code:
    """Minimal CSS code record (n, hx, hz as position masks); either list may be empty."""

    def __init__(self, name, n, hx, hz):
        self.name, self.n, self.hx, self.hz = name, int(n), list(hx), list(hz)


class Block:
    def __init__(self, name, qubits, code=None, refs=(), in_locs=()):
        self.name, self.qubits, self.code = name, list(qubits), code
        self.refs, self.in_locs = list(refs), list(in_locs)

    def to_json(self):
        d = dict(name=self.name, qubits=self.qubits, refs=self.refs, in_locs=self.in_locs)
        if self.code is not None:
            c = self.code
            d["code"] = dict(name=c.name, n=c.n, hx=[gf2.supp(r) for r in c.hx],
                             hz=[gf2.supp(r) for r in c.hz])
        return d

    @staticmethod
    def from_json(d):
        code = None
        if d.get("code"):
            c = d["code"]
            code = Code(c.get("name", "code"), c["n"], [gf2.from_supp(r) for r in c["hx"]],
                        [gf2.from_supp(r) for r in c["hz"]])
        return Block(d["name"], d["qubits"], code, d.get("refs", ()), d.get("in_locs", ()))


class Decoder:
    """A lookup decoder.  The error space has `npos` positions; `checks` are parity
    checks on it (position masks); bits[i] lists the measurements whose parity is the
    value of check i.  The decoded error (a lightest pattern with the observed check
    values) is applied according to `kind`:
        record   targets[p] is a measurement index: that reported outcome is believed
                 to be flipped (transversal measurement followed by classical decoding);
        X, Z     targets[p] is a qubit: a Pauli of that type is applied there as a frame."""

    def __init__(self, name, checks, bits, kind, targets, rule=None):
        self.name, self.checks, self.kind = name, list(checks), kind
        # rule: None for the lightest pattern with the observed check values, or
        # dict(type="hier", sector, outer=code json, inner=code json) for hierarchical
        # decoding of a concatenated code (ftqc_codes.hier), whose checks these must be
        self.rule = rule
        self.bits = [sorted(int(i) for i in b) for b in bits]
        self.targets = [int(t) for t in targets]
        self.npos = len(self.targets)
        if kind not in ("record", "X", "Z"):
            raise GadgetError(f"decoder {name}: unknown kind {kind!r}")
        if len(self.bits) != len(self.checks):
            raise GadgetError(f"decoder {name}: one set of measurements per check")
        for r in self.checks:
            if r >> self.npos:
                raise GadgetError(f"decoder {name}: a check refers to a missing position")

    def to_json(self):
        d = dict(name=self.name, checks=[gf2.supp(r) for r in self.checks], bits=self.bits,
                 kind=self.kind, targets=self.targets)
        if self.rule is not None:
            d["rule"] = self.rule
        return d

    @staticmethod
    def from_json(d):
        return Decoder(d["name"], [gf2.from_supp(r) for r in d["checks"]], d["bits"],
                       d["kind"], d["targets"], d.get("rule"))


class Gadget:
    def __init__(self, circuit, inputs, outputs, decoders, body_start, name=None, parts=None):
        self.circuit, self.inputs, self.outputs = circuit, inputs, outputs
        self.decoders, self.body_start = decoders, body_start
        self.name = name or circuit.name
        self.parts = list(parts or [])          # library fragments used, for composition receipts
        self.interface = [q for b in inputs for q in b.refs] + [q for b in outputs for q in b.qubits]
        if sorted(self.interface) != sorted(circuit.outputs):
            raise GadgetError("circuit outputs must be the references and the output blocks")
        if decoders and circuit.detectors is None:
            raise GadgetError("a gadget with decoders must declare its detectors")
        locs = circuit.locations()
        for b in inputs:
            if len(b.in_locs) != len(b.qubits):
                raise GadgetError(f"input block {b.name}: one input-error location per qubit")
            for q, i in zip(b.qubits, b.in_locs):
                k, _, op = locs[i]
                if op.name != "I" or op.qs != (q,) or k != body_start:
                    raise GadgetError(f"input block {b.name}: location {i} is not the input-error "
                                      f"idle of qubit {q}")
        nm = len(circuit.meas)
        for d in decoders:
            used = [i for b in d.bits for i in b] + (d.targets if d.kind == "record" else [])
            if any(i < 0 or i >= nm for i in used):
                raise GadgetError(f"decoder {d.name} refers to a missing measurement")
            if d.kind != "record" and any(q < 0 or q >= circuit.n for q in d.targets):
                raise GadgetError(f"decoder {d.name} refers to a missing qubit")

    def shape(self):
        """What must agree for two gadgets to be comparable."""
        return dict(inputs=[[len(b.qubits), len(b.refs)] for b in self.inputs],
                    outputs=[len(b.qubits) for b in self.outputs])

    def to_json(self):
        d = dict(name=self.name, circuit=self.circuit.to_json(), body_start=self.body_start,
                 inputs=[b.to_json() for b in self.inputs],
                 outputs=[b.to_json() for b in self.outputs],
                 decoders=[x.to_json() for x in self.decoders], parts=self.parts)
        d["sha256"] = self.sha256()
        return d

    def sha256(self):
        body = dict(circuit=self.circuit.sha256(), body_start=self.body_start,
                    inputs=[b.to_json() for b in self.inputs],
                    outputs=[b.to_json() for b in self.outputs],
                    decoders=[x.to_json() for x in self.decoders])
        return hashlib.sha256(json.dumps(body, sort_keys=True).encode()).hexdigest()

    @staticmethod
    def from_json(d):
        try:
            g = Gadget(Circuit.from_json(d["circuit"]), [Block.from_json(b) for b in d["inputs"]],
                       [Block.from_json(b) for b in d["outputs"]],
                       [Decoder.from_json(x) for x in d["decoders"]], d["body_start"],
                       d.get("name"), d.get("parts"))
        except (KeyError, TypeError, IndexError) as e:
            raise GadgetError(f"malformed gadget: {e}") from None
        if "sha256" in d and d["sha256"] != g.sha256():
            raise GadgetError("sha256 in gadget file does not match its contents")
        return g

    def counts(self):
        c = self.circuit.counts()
        n_in = sum(len(b.qubits) for b in self.inputs)
        c["fault_locations"] = c["locations"] - n_in
        c["input_error_locations"] = n_in
        return c


class Builder:
    """Incremental construction.  Declare inputs first, then the body.  `tick` places
    operations at the cursor and advances it; `par` runs several fragments from the
    same starting tick.  Measurements are referred to by the records `measure`
    returns and are numbered when the gadget is built."""

    def __init__(self, name, ideal=False):
        self.name, self.ideal = name, ideal
        self.ticks, self.cursor, self.nq = [], 0, 0
        self.inputs, self.outputs, self.decoders = [], [], []
        self.detectors = None
        self.body_start = None
        self.parts = []
        self._in_ops = []

    def alloc(self, n):
        qs = list(range(self.nq, self.nq + n))
        self.nq += n
        return qs

    def _emit(self, ops, noisy):
        while len(self.ticks) <= self.cursor:
            self.ticks.append([])
        recs = []
        for op in ops:
            nm = str(op[0])
            if not noisy and not nm.startswith("!"):
                nm = "!" + nm
            rec = [nm] + [int(q) for q in op[1:]]
            self.ticks[self.cursor].append(rec)
            recs.append(rec)
        self.cursor += 1
        return recs

    def _begin(self):
        if self.body_start is None:
            self.cursor = len(self.ticks)
            self.body_start = self.cursor
            for blk in self.inputs:
                self._in_ops.append(self._ops_for(blk))
            flat = [op for ops in self._in_ops for op in ops]
            if flat:
                while len(self.ticks) <= self.cursor:
                    self.ticks.append([])
                self.ticks[self.cursor].extend(flat)
                self.cursor += 1

    def _ops_for(self, blk):
        return [["I", q] for q in blk.qubits]

    def input(self, code=None, name=None, n=None):
        """Declare an input block; returns its qubits."""
        if self.body_start is not None:
            raise GadgetError("declare every input before the body")
        self.cursor = len(self.ticks)
        if code is None:
            qs, refs = self.alloc(n), self.alloc(n)
            self._emit([("RX", r) for r in refs] + [("R", q) for q in qs], False)
            self._emit([("CX", r, q) for r, q in zip(refs, qs)], False)
        else:
            qs = self.alloc(code.n)
            lx = logical_x(code)
            refs = self.alloc(len(lx))
            for layer in prep_layers(code, "zero", qs):
                self._emit(layer, False)
            self._emit([("RX", r) for r in refs], False)
            for r, v in zip(refs, lx):
                for p in gf2.supp(v):
                    self._emit([("CX", r, qs[p])], False)
        self.inputs.append(Block(name or f"in{len(self.inputs)}", qs, code, refs))
        return qs

    def tick(self, ops, noisy=None):
        """noisy=None: noisy unless the builder is ideal; True forces a noisy tick."""
        self._begin()
        return self._emit(ops, (not self.ideal) if noisy is None else noisy)

    def wires(self, qubits):
        """One noisy idle on each qubit: where a specification puts its faults."""
        return self.tick([("I", q) for q in qubits], True)

    def pair_wires(self, pairs):
        """One noisy two-qubit location per pair that is otherwise the identity (a noisy
        CX followed by a noiseless one): where the specification of a transversal
        two-block gate puts its faults."""
        self.tick([("CX", a, c) for a, c in pairs], True)
        self.tick([("CX", a, c) for a, c in pairs], False)

    def layers(self, layers, noisy=None):
        for layer in layers:
            self.tick(layer, noisy)

    def skip(self, k=1):
        self._begin()
        self.cursor += k

    def par(self, *thunks):
        self._begin()
        start, end, out = self.cursor, self.cursor, []
        for f in thunks:
            self.cursor = start
            out.append(f())
            end = max(end, self.cursor)
        self.cursor = end
        return out

    def measure(self, qubits, basis, noisy=None, extra=()):
        """Transversal measurement in one tick; `extra` are other operations for the
        same tick (idles on the survivors).  Returns the measurement records."""
        recs = self.tick([("M" if basis == "Z" else "MX", q) for q in qubits] + list(extra), noisy)
        return recs[:len(qubits)]

    def detector(self, records):
        if self.detectors is None:
            self.detectors = []
        self.detectors.append(list(records))

    def no_detectors(self):
        if self.detectors is None:
            self.detectors = []

    def checks(self, records, masks):
        """Postselect on each parity in `masks` (position masks over `records`)."""
        self.no_detectors()
        for m in masks:
            self.detector([records[p] for p in gf2.supp(m)])

    def decoder(self, name, records, checks, kind="record", targets=()):
        """Decoder for a transversal measurement: position p is records[p]; check i is
        the parity of the records on its support.  kind 'record' corrects the record
        itself; 'X' or 'Z' puts the decoded Pauli on targets[p]."""
        self.no_detectors()
        bits = [[records[p] for p in gf2.supp(r)] for r in checks]
        self.decoders.append((name, list(checks), bits, kind,
                              list(records) if kind == "record" else list(targets)))

    def syndrome_decoder(self, name, bit_records, checks, kind, targets):
        """Decoder fed by syndrome bits: bit_records[i] are the measurements whose
        parity is check i (rows of `checks` are masks over the positions of `targets`)."""
        if kind == "record":
            raise GadgetError("a syndrome decoder applies a Pauli frame (kind X or Z)")
        self.no_detectors()
        self.decoders.append((name, list(checks), [list(b) for b in bit_records], kind, list(targets)))

    def output(self, qubits, code=None, name=None):
        self.outputs.append(Block(name or f"out{len(self.outputs)}", qubits, code))

    def build(self, idles=False):
        self._begin()
        ticks = [list(t) for t in self.ticks]
        if idles:
            lo, hi = {}, {}
            for k, t in enumerate(ticks):
                for op in t:
                    for q in op[1:]:
                        lo.setdefault(q, k)
                        hi[q] = k
            for k in range(self.body_start + 1, len(ticks)):
                used = {q for op in ticks[k] for q in op[1:]}
                for q in range(self.nq):
                    if q in lo and max(lo[q], self.body_start) < k < hi[q] and q not in used:
                        ticks[k].append(["I", q])
        mi, li, meas_of, loc_of = 0, 0, {}, {}
        for t in ticks:
            for op in t:
                nm = op[0].lstrip("!")
                if nm in ("M", "MX"):
                    meas_of[id(op)] = mi
                    mi += 1
                if not op[0].startswith("!"):
                    loc_of[id(op)] = li
                    li += 1
        dets = None if self.detectors is None else [[meas_of[id(r)] for r in d] for d in self.detectors]
        refs = [q for b in self.inputs for q in b.refs]
        outs = [q for b in self.outputs for q in b.qubits]
        try:
            circ = Circuit([[tuple(op) for op in t] for t in ticks], self.nq, detectors=dets,
                           outputs=refs + outs, name=self.name)
        except CircuitError as e:
            raise GadgetError(str(e)) from None
        for blk, ops in zip(self.inputs, self._in_ops):
            blk.in_locs = [loc_of[id(op)] for op in ops]
        decs = [Decoder(nm, chk, [[meas_of[id(r)] for r in b] for b in bits], kind,
                        [meas_of[id(t)] for t in tg] if kind == "record" else tg)
                for nm, chk, bits, kind, tg in self.decoders]
        return Gadget(circ, self.inputs, self.outputs, decs, self.body_start or 0, self.name, self.parts)
