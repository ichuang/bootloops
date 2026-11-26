"""Leg B, part 1: the CNF whose unsatisfiability proves a minimum-weight bound.

Statement encoded, for parity rows P, detector rows D and a weight cap w:

    exists e in GF(2)^n :  P e = 0,  D_j . e = 1 for some j,  wt(e) <= w.

For a logical-operator distance, P is the check matrix the operator must
commute with and D is a basis of the kernel of the other check matrix, so
"D_j . e = 1 for some j" says e is not a stabilizer.  For the lightest
stabilizer, P is a basis of the kernel of the generators and D is the set of
unit vectors, so the condition says e is nonzero.

Only one direction of the encoding carries the proof: every e satisfying the
statement must extend to a satisfying assignment.  `extend` computes that
assignment explicitly, and the battery checks it against every clause, so the
direction is tested without a SAT solver.  UNSAT then means no such e exists.

Variables 1..n are the bits of e.  Parities use chained XOR definitions
(Tseitin), the weight cap uses the sequential counter of Sinz (2005).
"""
from . import gf2


class Formula:
    def __init__(self, n, P, D, w):
        if w < 0:
            raise ValueError("weight cap must be >= 0")
        self.n, self.w = n, w
        self.P, self.D = list(P), list(D)
        self.clauses = []
        self.nvars = n
        self._defs = []   # (var, a, b): var = a xor b, a and b literals > 0
        self._det = []    # literal carrying D_j . e
        self._cnt = {}    # (i, j) -> var s_{i,j}, i in 1..n-1, j in 1..w
        self._build()

    def _new(self):
        self.nvars += 1
        return self.nvars

    def _xor_chain(self, row):
        """Literal equal to the parity of e over supp(row); None for an empty row."""
        s = [j + 1 for j in gf2.supp(row)]
        if not s:
            return None
        acc = s[0]
        for x in s[1:]:
            t = self._new()
            self._defs.append((t, acc, x))
            self.clauses += [[-t, acc, x], [-t, -acc, -x], [t, -acc, x], [t, acc, -x]]
            acc = t
        return acc

    def _build(self):
        n, w = self.n, self.w
        for row in self.P:
            lit = self._xor_chain(row)
            if lit is not None:
                self.clauses.append([-lit])
        det = []
        for row in self.D:
            lit = self._xor_chain(row)
            if lit is not None:
                det.append(lit)
        self._det = det
        self.clauses.append(list(det))      # empty clause if D gives nothing
        # at most w of x_1..x_n
        if w == 0:
            self.clauses += [[-i] for i in range(1, n + 1)]
            return
        if w >= n:
            return
        s = self._cnt
        for i in range(1, n):
            for j in range(1, w + 1):
                s[i, j] = self._new()
        self.clauses.append([-1, s[1, 1]])
        for j in range(2, w + 1):
            self.clauses.append([-s[1, j]])
        for i in range(2, n):
            self.clauses.append([-i, s[i, 1]])
            self.clauses.append([-s[i - 1, 1], s[i, 1]])
            for j in range(2, w + 1):
                self.clauses.append([-i, -s[i - 1, j - 1], s[i, j]])
                self.clauses.append([-s[i - 1, j], s[i, j]])
            self.clauses.append([-i, -s[i - 1, w]])
        self.clauses.append([-n, -s[n - 1, w]])

    # -- the proof-carrying direction ---------------------------------------
    def extend(self, e):
        """Full assignment (list indexed by variable, entry 0 unused) induced
        by the error vector e (an int)."""
        val = [False] * (self.nvars + 1)
        for j in range(self.n):
            val[j + 1] = bool((e >> j) & 1)
        for t, a, b in self._defs:
            val[t] = val[a] != val[b]
        run = 0
        for i in range(1, self.n):
            run += val[i]
            for j in range(1, self.w + 1):
                if (i, j) in self._cnt:
                    val[self._cnt[i, j]] = run >= j
        return val

    def satisfied_by(self, val):
        return all(any(val[l] if l > 0 else not val[-l] for l in c) for c in self.clauses)

    def holds(self, e):
        """The statement itself, evaluated directly on e."""
        return (all(gf2.dot(e, r) == 0 for r in self.P)
                and any(gf2.dot(e, r) for r in self.D)
                and gf2.wt(e) <= self.w)

    def dimacs(self):
        out = [f"p cnf {self.nvars} {len(self.clauses)}"]
        out += [" ".join(map(str, c)) + " 0" if c else "0" for c in self.clauses]
        return "\n".join(out) + "\n"


def for_quantity(code, q, w):
    """Formula asserting that quantity q (d_x, d_z, wmin_x, wmin_z) is <= w."""
    n = code.n
    units = [1 << j for j in range(n)]
    if q == "d_x":
        return Formula(n, code.hz, gf2.nullspace(code.hx, n), w)
    if q == "d_z":
        return Formula(n, code.hx, gf2.nullspace(code.hz, n), w)
    if q == "wmin_x":
        return Formula(n, gf2.nullspace(code.hx, n), units, w)
    if q == "wmin_z":
        return Formula(n, gf2.nullspace(code.hz, n), units, w)
    raise ValueError(q)
