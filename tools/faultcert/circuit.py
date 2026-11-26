"""Clifford gadget circuits as lists of ticks.

Operations (name, qubits...):
    R q      prepare |0>            RX q     prepare |+>
    I q      idle                   H q      Hadamard
    CX c t   controlled-X           M q      measure Z          MX q   measure X

A tick is a list of operations on disjoint qubits.  Measurements are numbered
in order of appearance.  An operation may be marked noiseless (for ideal
encoders or ideal final readout); every other operation is a fault location.

Fault model.  One fault = one location misbehaving:
  after R, RX, I, H        any non-identity Pauli on the qubit
  after CX                 any non-identity Pauli on the two qubits
  at M, MX                 the reported outcome is flipped
The weight of a fault set is the number of locations in it.
"""
import hashlib
import json

ONE_Q = {"R", "RX", "I", "H", "M", "MX"}
ALIASES = {"P0": "R", "PZ": "R", "P": "R", "P+": "RX", "PX": "RX", "PP": "RX",
           "CNOT": "CX", "MZ": "M"}


class CircuitError(ValueError):
    pass


class Op:
    __slots__ = ("name", "qs", "noisy")

    def __init__(self, name, qs, noisy=True):
        self.name, self.qs, self.noisy = name, tuple(qs), noisy

    def text(self):
        return ("" if self.noisy else "!") + self.name + " " + " ".join(map(str, self.qs))


class Circuit:
    def __init__(self, ticks, n=None, detectors=None, observables=None, outputs=None, name="circuit"):
        """ticks: list of lists of (name, q...) tuples or Op; a name starting
        with '!' marks a noiseless operation.
        detectors / observables: lists of measurement-index lists.  detectors=None
        means "discover them" (every deterministic outcome parity); that is only
        allowed without observables, i.e. for preparations and other gadgets where
        every deterministic parity is a legitimate check.
        outputs: qubits that carry the gadget's output; default: qubits not measured
        after their last preparation."""
        self.name = name
        self.ticks = []
        top = -1
        for k, tick in enumerate(ticks):
            seen, row = set(), []
            for op in tick:
                if not isinstance(op, Op):
                    nm = str(op[0])
                    noisy = not nm.startswith("!")
                    nm = nm.lstrip("!").upper()
                    nm = ALIASES.get(nm, nm)
                    op = Op(nm, [int(q) for q in op[1:]], noisy)
                if op.name == "CX":
                    if len(op.qs) != 2 or op.qs[0] == op.qs[1]:
                        raise CircuitError(f"tick {k}: CX needs two distinct qubits")
                elif op.name in ONE_Q:
                    if len(op.qs) != 1:
                        raise CircuitError(f"tick {k}: {op.name} takes one qubit")
                else:
                    raise CircuitError(f"tick {k}: unknown operation {op.name!r}")
                for q in op.qs:
                    if q < 0:
                        raise CircuitError("negative qubit index")
                    if q in seen:
                        raise CircuitError(f"tick {k}: qubit {q} used twice")
                    seen.add(q)
                    top = max(top, q)
                row.append(op)
            self.ticks.append(row)
        self.n = top + 1 if n is None else int(n)
        if self.n <= top:
            raise CircuitError("n is smaller than the largest qubit index")
        # every qubit must be prepared before anything else touches it
        prepared = set()
        self.meas = []                     # (tick, op index, qubit, basis)
        last_prep, last_meas = {}, {}
        for k, tick in enumerate(self.ticks):
            for j, op in enumerate(tick):
                if op.name in ("R", "RX"):
                    prepared.add(op.qs[0])
                    last_prep[op.qs[0]] = k
                else:
                    for q in op.qs:
                        if q not in prepared:
                            raise CircuitError(f"tick {k}: qubit {q} used before preparation")
                if op.name in ("M", "MX"):
                    last_meas[op.qs[0]] = k
                    self.meas.append((k, j, op.qs[0], "Z" if op.name == "M" else "X"))
        if outputs is None:
            outputs = [q for q in sorted(prepared)
                       if q not in last_meas or last_meas[q] < last_prep[q]]
        self.outputs = sorted(int(q) for q in outputs)
        nm = len(self.meas)
        for group in (detectors or []), (observables or []):
            for d in group:
                if any(i < 0 or i >= nm for i in d):
                    raise CircuitError("detector or observable refers to a missing measurement")
        if observables and detectors is None:
            # A logical observable is itself a deterministic parity, so discovery would
            # list it (or a detector containing it) and nothing could flip it unseen.
            raise CircuitError("a circuit that declares observables must declare its detectors")
        self.detectors = None if detectors is None else [sorted(set(d)) for d in detectors]
        self.observables = [sorted(set(d)) for d in (observables or [])]

    # ------------------------------------------------------------------
    def locations(self):
        """[(tick, op index, Op)] for the noisy operations, in circuit order."""
        return [(k, j, op) for k, tick in enumerate(self.ticks) for j, op in enumerate(tick) if op.noisy]

    def canonical(self):
        lines = [f"n {self.n}", "outputs " + " ".join(map(str, self.outputs))]
        for tick in self.ticks:
            lines.append("TICK")
            lines += [op.text() for op in tick]
        if self.detectors is not None:
            lines += ["DET " + " ".join(map(str, d)) for d in self.detectors]
        lines += ["OBS " + " ".join(map(str, d)) for d in self.observables]
        return "\n".join(lines) + "\n"

    def sha256(self):
        return hashlib.sha256(self.canonical().encode()).hexdigest()

    def to_json(self):
        return dict(name=self.name, n=self.n, outputs=self.outputs,
                    ticks=[[[("" if op.noisy else "!") + op.name, *op.qs] for op in t]
                           for t in self.ticks],
                    detectors=self.detectors, observables=self.observables,
                    sha256=self.sha256())

    @staticmethod
    def from_json(d):
        try:
            c = Circuit(d["ticks"], d.get("n"), d.get("detectors"), d.get("observables"),
                        d.get("outputs"), d.get("name", "circuit"))
        except (KeyError, TypeError, IndexError) as e:
            raise CircuitError(f"malformed circuit: {e}") from None
        if "sha256" in d and d["sha256"] != c.sha256():
            raise CircuitError("sha256 in circuit file does not match its contents")
        return c

    def counts(self):
        out = {}
        for tick in self.ticks:
            for op in tick:
                out[op.name] = out.get(op.name, 0) + 1
        return dict(qubits=self.n, ticks=len(self.ticks), measurements=len(self.meas),
                    locations=len(self.locations()), ops=out)


def load(path):
    with open(path) as f:
        return Circuit.from_json(json.load(f))
