"""Batched adjoint fault signatures for CSS circuits with record-only endpoints."""

from .model import propagate, PAULI_1


class RecordModel:
    def __init__(self, circuit):
        if circuit.outputs or circuit.detectors is None:
            raise ValueError(
                "explicit record detectors and empty quantum output required"
            )
        self.circuit = circuit
        self.groups = []
        self.locations = circuit.locations()
        self.nd = len(circuit.detectors)
        self.mask = (1 << self.nd) - 1
        self.detectors = circuit.detectors
        self.parities = [*circuit.detectors, *circuit.observables]
        records = [0] * len(circuit.meas)
        for j, ids in enumerate(self.parities):
            for k in ids:
                records[k] |= 1 << j
        self.record_signatures = records
        self.meas = {(t, i): j for j, (t, i, _, _) in enumerate(circuit.meas)}
        sx, sz = [0] * circuit.n, [0] * circuit.n
        for t in reversed(range(len(circuit.ticks))):
            for i in reversed(range(len(circuit.ticks[t]))):
                op = circuit.ticks[t][i]
                g, q = op.name, op.qs[0]
                if g in ("M", "MX"):
                    if sx[q] or sz[q]:
                        raise ValueError("measured wire reused without reset")
                    v = records[self.meas[t, i]]
                    if op.noisy:
                        self.groups.append({"flip": v})
                    sx[q], sz[q] = (v, 0) if g == "M" else (0, v)
                    continue
                if op.noisy:
                    group = {}
                    for m in range(1, 1 << (2 * len(op.qs))):
                        value = 0
                        label = ""
                        for k, site in enumerate(op.qs):
                            p = (m >> (2 * k)) & 3
                            label += PAULI_1.get(p, "I")
                            if p & 1:
                                value ^= sx[site]
                            if p & 2:
                                value ^= sz[site]
                        group[label] = value
                    self.groups.append(group)
                if g == "CX":
                    target = op.qs[1]
                    sx[q] ^= sx[target]
                    sz[target] ^= sz[q]
                elif g in ("R", "RX"):
                    if sz[q] if g == "R" else sx[q]:
                        raise ValueError("random declared parity")
                    sx[q] = sz[q] = 0
                elif g != "I":
                    raise ValueError("unsupported operation in CSS record model")
        self.groups.reverse()
        if any(sx) or any(sz):
            raise ValueError("unclosed input parity")

    @property
    def columns(self):
        return {
            (v & self.mask, v >> self.nd, 0)
            for group in self.groups
            for v in group.values()
            if v
        }

    def replay(self, faults):
        """Cross-check the chosen fault set through faultcert's forward transporter."""
        if len({i for i, _ in faults}) != len(faults):
            raise ValueError("repeated fault site")
        raw = expected = 0
        for loc, label in faults:
            if (
                type(loc) is not int
                or not 0 <= loc < len(self.locations)
                or label not in self.groups[loc]
            ):
                raise ValueError("invalid physical witness")
            t, i, op = self.locations[loc]
            x = z = rec = 0
            if label == "flip":
                rec = 1 << self.meas[t, i]
            else:
                for q, p in zip(op.qs, label):
                    if p in "XY":
                        x |= 1 << q
                    if p in "ZY":
                        z |= 1 << q
            _, _, f = propagate(self.circuit, t, x, z, rec)
            raw ^= f
            expected ^= self.groups[loc][label]
        actual = 0
        while raw:
            b = raw & -raw
            actual ^= self.record_signatures[b.bit_length() - 1]
            raw -= b
        if actual != expected:
            raise ValueError("forward and adjoint witness signatures disagree")
        return actual & self.mask, actual >> self.nd
