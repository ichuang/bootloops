"""Build a deterministic syndrome recovery, retaining a collision witness."""
from ftqc_codes.gf2 import Echelon, dot


def syndrome(word, checks):
    return sum(dot(word, h) << i for i, h in enumerate(checks))


class Recovery:
    def __init__(self, hx, hz):
        self.hx, self.hz = hx, hz
        self.x_stabs, self.z_stabs = Echelon(hx), Echelon(hz)
        self.table = {}
        self.origins = {}

    def observe(self, history, x, z, origin):
        key = (tuple(history), syndrome(x, self.hz), syndrome(z, self.hx))
        value = (self.x_stabs.reduce(x), self.z_stabs.reduce(z))
        if key in self.table and self.table[key] != value:
            return dict(record=key, first=self.origins[key], second=origin,
                        first_residual=self.table[key], second_residual=value)
        if key not in self.table:
            self.table[key] = value
            self.origins[key] = origin
        return None

    def rows(self):
        return [[list(history), sx, sz, x, z]
                for (history, sx, sz), (x, z) in sorted(self.table.items())]

