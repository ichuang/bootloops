"""The fault model of a circuit: what each fault does, as vectors over GF(2).

For Pauli faults in a Clifford circuit the effect of a fault set is the XOR of
the effects of its members.  The effect of one fault is
    flips     which measurement outcomes change,
    (x, z)    the residual Pauli left on the qubits at the end,
from which the model derives
    det       which detectors fire,
    obs       which logical observables flip,
    syn       which stabilizers of the output state are violated, counting both the
              residual and a Pauli frame made wrong by flipped outcomes (zero iff
              the output state is the one the reported outcomes imply).

`Model.columns` lists the distinct non-trivial (det, obs, syn) signatures that a
single location can produce, over all its non-identity Paulis.  A set of j
locations produces a sum of at most j columns and conversely, because the
signatures a location can produce are closed under addition.  Every question
about "at most j faults" is therefore a question about sums of at most j columns.
"""
from . import stab


def propagate(circuit, tick, x, z, flips=0):
    """Push a Pauli (x, z) present after tick `tick` to the end of the circuit.
    Returns (x, z, flips)."""
    mi = sum(1 for m in circuit.meas if m[0] <= tick)       # next measurement index
    for k in range(tick + 1, len(circuit.ticks)):
        for op in circuit.ticks[k]:
            nm = op.name
            if nm == "CX":
                c, t = op.qs
                if (x >> c) & 1:
                    x ^= 1 << t
                if (z >> t) & 1:
                    z ^= 1 << c
                continue
            q = op.qs[0]
            b = 1 << q
            if nm == "H":
                xb, zb = x & b, z & b
                if bool(xb) != bool(zb):
                    x ^= b
                    z ^= b
            elif nm in ("R", "RX"):
                x &= ~b
                z &= ~b
            elif nm == "M":
                if x & b:
                    flips ^= 1 << mi
                z &= ~b                  # Z on a Z eigenstate is a phase
                mi += 1
            elif nm == "MX":
                if z & b:
                    flips ^= 1 << mi
                x &= ~b
                mi += 1
    return x, z, flips


PAULI_1 = {1: "X", 2: "Z", 3: "Y"}


class Location:
    __slots__ = ("index", "tick", "op", "faults")

    def __init__(self, index, tick, op):
        self.index, self.tick, self.op = index, tick, op
        self.faults = {}         # label -> (x, z, flips) at the end of the circuit


class Model:
    def __init__(self, circuit):
        self.circuit = circuit
        sim = stab.Sim(circuit.n).run(circuit)
        self.sim = sim
        found = sim.detectors()
        if circuit.detectors is None:
            self.detectors = [d for d, _ in found]
        else:
            self.detectors = circuit.detectors
            for d in self.detectors:
                if sim.parity_form(d) >> 1:
                    raise ValueError(f"declared detector {d} is not deterministic")
        self.det_masks = [sum(1 << i for i in d) for d in self.detectors]
        self.obs_masks = [sum(1 << i for i in o) for o in circuit.observables]
        for o in circuit.observables:
            if sim.parity_form(o) >> 1:
                raise ValueError(f"observable {o} is not deterministic in the noiseless circuit")
        self.out_mask = sum(1 << q for q in circuit.outputs)
        # stabilizer group of the output state, with the outcomes each sign depends on
        self.final = [(g.x, g.z, sim.frame_mask(g.s)) for g in sim.output_group(self.out_mask)]
        self.locations = []
        meas_index = {(m[0], m[1]): i for i, m in enumerate(circuit.meas)}
        for idx, (k, j, op) in enumerate(circuit.locations()):
            loc = Location(idx, k, op)
            if op.name in ("M", "MX"):
                loc.faults["flip"] = (0, 0, 1 << meas_index[(k, j)])
            else:
                gens = []
                for q in op.qs:
                    gens.append(("X", q, propagate(circuit, k, 1 << q, 0)))
                    gens.append(("Z", q, propagate(circuit, k, 0, 1 << q)))
                for m in range(1, 1 << len(gens)):
                    x = z = fl = 0
                    per_q = {}
                    for i, (p, q, (gx, gz, gf)) in enumerate(gens):
                        if (m >> i) & 1:
                            x ^= gx
                            z ^= gz
                            fl ^= gf
                            per_q[q] = per_q.get(q, 0) | (1 if p == "X" else 2)
                    label = "".join(PAULI_1[per_q[q]] if q in per_q else "I" for q in op.qs)
                    loc.faults[label] = (x, z, fl)
            self.locations.append(loc)
        self._columns()

    # -- signatures -----------------------------------------------------
    def signature(self, x, z, flips):
        det = 0
        for i, m in enumerate(self.det_masks):
            if (flips & m).bit_count() & 1:
                det |= 1 << i
        obs = 0
        for i, m in enumerate(self.obs_masks):
            if (flips & m).bit_count() & 1:
                obs |= 1 << i
        return det, obs, self.syndrome(x, z, flips)

    def syndrome(self, x, z, flips=0):
        """Which output-state stabilizers the actual output violates, relative to the
        state the reported outcomes imply.  Bit i is the commutation of the residual
        with generator i, plus the parity of the flipped outcomes that the sign of
        generator i depends on (a wrong Pauli frame is an error like any other).
        Zero iff the output state is the believed one."""
        x &= self.out_mask
        z &= self.out_mask
        s = 0
        for i, (gx, gz, fm) in enumerate(self.final):
            if ((x & gz).bit_count() + (z & gx).bit_count() + (flips & fm).bit_count()) & 1:
                s |= 1 << i
        return s

    def _columns(self):
        cols = {}
        for loc in self.locations:
            for label, eff in loc.faults.items():
                sig = self.signature(*eff)
                if sig == (0, 0, 0):
                    continue
                cols.setdefault(sig, (loc.index, label))
        self.columns = sorted(cols)               # distinct signatures
        self.column_source = [cols[s] for s in self.columns]   # one (location, label) each

    # -- direct replay of a fault set (independent of the column bookkeeping) -------
    def replay(self, faults):
        """faults: [(location index, label)], distinct locations.  Returns the signature."""
        if len({i for i, _ in faults}) != len(faults):
            raise ValueError("a location appears twice in a fault set")
        x = z = fl = 0
        for i, label in faults:
            loc = self.locations[i]
            if loc.op.name in ("M", "MX"):
                if label != "flip":
                    raise ValueError("measurement locations take the label 'flip'")
                fx, fz, ff = loc.faults["flip"]
            else:
                if len(label) != len(loc.op.qs) or set(label) - set("IXYZ") or set(label) == {"I"}:
                    raise ValueError(f"bad Pauli label {label!r} for location {i}")
                px = sum(1 << q for q, p in zip(loc.op.qs, label) if p in "XY")
                pz = sum(1 << q for q, p in zip(loc.op.qs, label) if p in "ZY")
                fx, fz, ff = propagate(self.circuit, loc.tick, px, pz)
            x ^= fx
            z ^= fz
            fl ^= ff
        return self.signature(x, z, fl)

    def describe(self, faults):
        out = []
        for i, label in faults:
            loc = self.locations[i]
            out.append(dict(location=i, tick=loc.tick, op=loc.op.text(), fault=label))
        return out

    def stats(self):
        return dict(locations=len(self.locations),
                    faults=sum(len(l.faults) for l in self.locations),
                    columns=len(self.columns), detectors=len(self.detectors),
                    observables=len(self.obs_masks), outputs=len(self.circuit.outputs))
