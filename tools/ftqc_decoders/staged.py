"""Build a deterministic syndrome recovery, retaining a collision witness."""

from ftqc_codes.gf2 import Echelon, dot


def syndrome(word, checks):
    return sum(dot(word, h) << i for i, h in enumerate(checks))


class Recovery:
    def __init__(self, hx, hz, supplied=None):
        self.hx, self.hz = hx, hz
        self.x_stabs, self.z_stabs = Echelon(hx), Echelon(hz)
        self.table = {}
        self.origins = {}
        self.supplied = None
        if supplied is not None:
            self.supplied = {}
            width = max(h.bit_length() for h in [*hx, *hz])
            for row in supplied:
                if not isinstance(row, list) or len(row) != 5:
                    raise ValueError("invalid recovery row")
                history, sx, sz, x, z = row
                if (
                    not isinstance(history, list)
                    or any(type(v) is not int or v < 0 for v in history)
                    or any(
                        type(v) is not int or not 0 <= v < 1 << bits
                        for v, bits in (
                            (sx, len(hz)),
                            (sz, len(hx)),
                            (x, width),
                            (z, width),
                        )
                    )
                ):
                    raise ValueError("invalid recovery record or correction")
                key = tuple(history), sx, sz
                if key in self.supplied:
                    raise ValueError("duplicate recovery record")
                self.supplied[key] = (x, z)

    def decode(self, history, sx, sz):
        """Perform the correction from observed bits; default identity off-table."""
        table = self.table if self.supplied is None else self.supplied
        return table.get((tuple(history), sx, sz), (0, 0))

    def observe(self, history, x, z, origin):
        key = (tuple(history), syndrome(x, self.hz), syndrome(z, self.hx))
        value = (self.x_stabs.reduce(x), self.z_stabs.reduce(z))
        if self.supplied is not None:
            cx, cz = self.decode(*key)
            residual = self.x_stabs.reduce(x ^ cx), self.z_stabs.reduce(z ^ cz)
            if any(residual):
                return dict(
                    record=key, origin=origin, correction=[cx, cz], residual=residual
                )
        if key in self.table and self.table[key] != value:
            return dict(
                record=key,
                first=self.origins[key],
                second=origin,
                first_residual=self.table[key],
                second_residual=value,
            )
        if key not in self.table:
            self.table[key] = value
            self.origins[key] = origin
        return None

    def rows(self):
        return [
            [list(history), sx, sz, x, z]
            for (history, sx, sz), (x, z) in sorted(
                (self.table if self.supplied is None else self.supplied).items()
            )
        ]
