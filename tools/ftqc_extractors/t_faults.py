"""Independent forward and adjoint propagation of Pauli fault signatures."""

from collections import defaultdict
from functools import reduce, cached_property
from operator import xor
from .algebra import bits
from .t_schedule import TSchedule


class FaultCircuit(TSchedule):
    @cached_property
    def forward(self):
        # Pack independent single-fault experiments into bit planes.
        x, z, rec = defaultdict(int), defaultdict(int), {}
        offset = 0
        for t, (g, q, *tail) in enumerate(self.ops):
            modes = self.modes[t]
            flags = [
                sum(1 << (offset + j) for j, mode in enumerate(modes) if mode >> i & 1)
                for i in range(4)
            ]
            offset += len(modes)
            if g in ("MZ", "MX"):
                rec[t] = (x if g == "MZ" else z)[q] ^ flags[0]
                x[q] = z[q] = 0
                continue
            if g in ("RX", "RZ"):
                x[q] = flags[0] if g == "RZ" else 0
                z[q] = flags[0] if g == "RX" else 0
                continue
            if g == "CX":
                other = tail[0]
                x[other] ^= x[q]
                z[q] ^= z[other]
                x[other] ^= flags[2]
                z[other] ^= flags[3]
            if g == "H":
                x[q], z[q] = z[q], x[q]
            if g in ("S", "SDG"):
                z[q] ^= x[q]
            x[q] ^= flags[0]
            z[q] ^= flags[1]
        signatures = [0] * offset
        for i, (px, pz, r) in enumerate(self.flows):
            row = (
                reduce(xor, (z[j] for j in bits(px)), 0)
                ^ reduce(xor, (x[j] for j in bits(pz)), 0)
                ^ reduce(xor, (rec[j] for j in bits(r)), 0)
            )
            for j in bits(row):
                signatures[j] |= 1 << i
        out = []
        offset = 0
        for modes in self.modes:
            out.append(list(zip(modes, signatures[offset : offset + len(modes)])))
            offset += len(modes)
        return out

    @cached_property
    def backward(self):
        # Propagate the output tests instead of the input fault patterns.
        sx, sz, rec = defaultdict(int), defaultdict(int), defaultdict(int)
        for i, (px, pz, r) in enumerate(self.flows):
            for j in bits(px):
                sz[j] |= 1 << i
            for j in bits(pz):
                sx[j] |= 1 << i
            for j in bits(r):
                rec[j] |= 1 << i
        out = [None] * len(self.ops)
        for t in reversed(range(len(self.ops))):
            g, q, *tail = self.ops[t]
            if g[0] == "M":
                out[t] = [(1, rec[t])]
                sx[q] = rec[t] if g == "MZ" else 0
                sz[q] = rec[t] if g == "MX" else 0
            elif g[0] == "R":
                out[t] = [(1, sx[q] if g == "RZ" else sz[q])]
                sx[q] = sz[q] = 0
            else:
                v = [sx[q], sz[q]]
                if tail:
                    v.extend([sx[tail[0]], sz[tail[0]]])
                out[t] = [
                    (mode, reduce(xor, (v[j] for j in bits(mode)), 0))
                    for mode in self.modes[t]
                ]
                if g == "CX":
                    other = tail[0]
                    sx[q] ^= sx[other]
                    sz[other] ^= sz[q]
                if g == "H":
                    sx[q], sz[q] = sz[q], sx[q]
                if g in ("S", "SDG"):
                    sx[q] ^= sz[q]
        return out

    def events(self, reverse=False):
        return self.backward if reverse else self.forward

    def signature(self, faults):
        return reduce(xor, (dict(self.backward[t])[m] for t, m in faults.items()), 0)


def pair_witness(schedule, reverse=False):
    """Two distinct locations are silent together iff their syndromes agree."""
    seen = {}
    for t, group in enumerate(schedule.events(reverse)):
        for mode, value in group:
            syndrome = value & schedule.mask
            for old, loc, m in seen.get(syndrome, []):
                if old != value:
                    faults = [[loc, m], [t, mode]]
                    if schedule.signature(dict(faults)) != old ^ value:
                        raise ValueError("fault replay mismatch")
                    return faults
        for mode, value in group:
            bucket = seen.setdefault(value & schedule.mask, [])
            if len(bucket) < 2 and all(record[0] != value for record in bucket):
                bucket.append((value, t, mode))
    return None
