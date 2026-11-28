"""GF(2) linear algebra on Python ints (bit j of a row = coordinate j). Stdlib only."""


def wt(v):
    return v.bit_count()


def dot(a, b):
    return (a & b).bit_count() & 1


def rows_from_lists(rows):
    """[[0,1,...],...] -> list of ints."""
    out = []
    for r in rows:
        v = 0
        for j, x in enumerate(r):
            if x not in (0, 1):
                raise ValueError("matrix entries must be 0 or 1")
            if x:
                v |= 1 << j
        out.append(v)
    return out


def rows_to_lists(rows, n):
    return [[(r >> j) & 1 for j in range(n)] for r in rows]


def supp(v):
    out, j = [], 0
    while v:
        if v & 1:
            out.append(j)
        v >>= 1
        j += 1
    return out


def from_supp(s):
    v = 0
    for j in s:
        v |= 1 << j
    return v


class Echelon:
    """Reduced echelon basis of a row space, with membership and reduction."""

    def __init__(self, rows=()):
        self.basis = []   # reduced rows
        self.pivots = []  # pivot bit of each row (its highest set bit)
        self.picked = []  # indices of the input rows that were independent
        self._seen = 0
        for r in rows:
            self.add(r)

    def reduce(self, v):
        for b, p in zip(self.basis, self.pivots):
            if (v >> p) & 1:
                v ^= b
        return v

    def add(self, r):
        idx = self._seen
        self._seen += 1
        r = self.reduce(r)
        if not r:
            return False
        p = r.bit_length() - 1
        for i in range(len(self.basis)):
            if (self.basis[i] >> p) & 1:
                self.basis[i] ^= r
        self.basis.append(r)
        self.pivots.append(p)
        self.picked.append(idx)
        return True

    def __contains__(self, v):
        return self.reduce(v) == 0

    @property
    def rank(self):
        return len(self.basis)


def rank(rows):
    return Echelon(rows).rank


def independent_rows(rows):
    """Indices of a maximal independent subset, in input order."""
    return Echelon(rows).picked


def nullspace(rows, n):
    """Basis of {v in GF(2)^n : v.r = 0 for every row r}."""
    e = Echelon(rows)
    piv = set(e.pivots)
    out = []
    for f in range(n):
        if f in piv:
            continue
        v = 1 << f
        for b, p in zip(e.basis, e.pivots):
            if (b >> f) & 1:
                v |= 1 << p
        out.append(v)
    return out


def systematic(gens, n, cols):
    """Row-reduce `gens` choosing pivots only among `cols` (in that order of
    preference).  Returns (rows, pivot_cols): rows[i] has a 1 in pivot_cols[i]
    and 0 in the other pivot columns.  len(rows) is the rank of gens restricted
    to cols; generators with no pivot in cols are dropped."""
    rows = list(gens)
    out, piv = [], []
    for c in cols:
        m = 1 << c
        k = next((i for i, r in enumerate(rows) if r & m), None)
        if k is None:
            continue
        r = rows.pop(k)
        rows = [x ^ r if x & m else x for x in rows]
        out = [x ^ r if x & m else x for x in out]
        out.append(r)
        piv.append(c)
    return out, piv
