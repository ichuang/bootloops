"""Symbolic stabilizer simulation of the noiseless circuit.

The state is a full set of n stabilizer generators.  Each generator is a Pauli
(x, z bitmasks) with a sign that is an affine GF(2) form over "random-outcome
variables": bit 0 of the form is the constant, bit i >= 1 is the variable
introduced by the i-th measurement whose outcome was random.  (Symbolic
stabilizer states in this sense are the device of Chen et al., arXiv:2501.14380.)

Two things are read off:
  * every measurement outcome as an affine form in earlier random outcomes; an
    outcome that introduces no new variable is deterministic given the earlier
    ones, which yields a detector (a parity of outcomes with a fixed value);
  * the final generators, which define when a residual Pauli is harmless.
"""


def _g_sum(x1, z1, x2, z2):
    """Exponent of i (mod 4) picked up when multiplying Pauli 1 by Pauli 2
    (Aaronson-Gottesman convention, summed over qubits)."""
    y1, xo, zo = x1 & z1, x1 & ~z1, z1 & ~x1
    plus = (y1 & z2 & ~x2) | (xo & z2 & x2) | (zo & x2 & ~z2)
    minus = (y1 & x2 & ~z2) | (xo & z2 & ~x2) | (zo & x2 & z2)
    return (plus.bit_count() - minus.bit_count()) % 4


class Gen:
    __slots__ = ("x", "z", "s")

    def __init__(self, x, z, s=0):
        self.x, self.z, self.s = x, z, s

    def mul(self, other):
        """self <- other * self (both Hermitian and commuting, so the phase is real)."""
        e = _g_sum(other.x, other.z, self.x, self.z)
        if e % 2:
            raise AssertionError("product of anticommuting generators")
        self.s ^= other.s ^ (1 if e == 2 else 0)
        self.x ^= other.x
        self.z ^= other.z


def _anti(g, px, pz):
    return ((g.x & pz).bit_count() + (g.z & px).bit_count()) & 1


class Sim:
    def __init__(self, n):
        self.n = n
        self.gens = [Gen(0, 1 << q) for q in range(n)]      # |0...0>
        self.nvars = 0
        self.var_meas = {}        # variable index -> measurement index that introduced it
        self.outcomes = []        # affine form of each measurement outcome

    # -- Clifford gates -------------------------------------------------
    def h(self, q):
        b = 1 << q
        for g in self.gens:
            xb, zb = g.x & b, g.z & b
            if xb and zb:
                g.s ^= 1
            if bool(xb) != bool(zb):
                g.x ^= b
                g.z ^= b

    def cx(self, c, t):
        bc, bt = 1 << c, 1 << t
        for g in self.gens:
            xc, zc, xt, zt = g.x & bc, g.z & bc, g.x & bt, g.z & bt
            if xc and zt and (bool(xt) == bool(zc)):
                g.s ^= 1
            if xc:
                g.x ^= bt
            if zt:
                g.z ^= bc

    # -- measurement of a Pauli (px, pz), sign + ------------------------
    def _express(self, px, pz):
        """(sign form, mask of generators) of the group element equal to the Pauli
        (px, pz); the Pauli must commute with every generator."""
        rows = [(g.x, g.z, 1 << i) for i, g in enumerate(self.gens)]
        tx, tz, used = px, pz, 0
        n = self.n
        for col in range(2 * n):
            def bit(x, z):
                return (x >> col) & 1 if col < n else (z >> (col - n)) & 1
            piv = next((i for i, r in enumerate(rows) if bit(r[0], r[1])), None)
            if piv is None:
                continue
            rx, rz, rm = rows.pop(piv)
            rows = [(x ^ rx, z ^ rz, m ^ rm) if bit(x, z) else (x, z, m) for x, z, m in rows]
            if bit(tx, tz):
                tx, tz, used = tx ^ rx, tz ^ rz, used ^ rm
        if tx or tz:
            raise AssertionError("Pauli commutes with the group but is not in it")
        acc = Gen(0, 0, 0)
        for i, g in enumerate(self.gens):
            if (used >> i) & 1:
                acc.mul(g)
        assert acc.x == px and acc.z == pz
        return acc.s, used

    def _measure(self, px, pz, reset):
        anti = [i for i, g in enumerate(self.gens) if _anti(g, px, pz)]
        if anti:
            g0 = self.gens[anti[0]]
            for i in anti[1:]:
                self.gens[i].mul(g0)
            if reset:
                form = 0
            else:
                self.nvars += 1
                form = 1 << self.nvars
            self.gens[anti[0]] = Gen(px, pz, form)
            return form, True
        form, used = self._express(px, pz)
        if reset:
            # change basis so that the Pauli is itself a generator, strip it from the
            # other generators, then set its sign to + (the reset)
            i0 = (used & -used).bit_length() - 1
            me = Gen(px, pz, form)
            self.gens[i0] = me
            b = px | pz
            for i, g in enumerate(self.gens):
                if i != i0 and ((g.x | g.z) & b):
                    g.mul(me)
            me.s = 0
        return form, False

    def run(self, circuit):
        for tick in circuit.ticks:
            for op in tick:
                nm, q = op.name, op.qs[0]
                if nm == "H":
                    self.h(q)
                elif nm in ("S", "SDG"):
                    for g in self.gens:
                        b = 1 << q
                        if g.x & b:
                            if bool(g.z & b) != (nm == "SDG"):
                                g.s ^= 1
                            g.z ^= b
                elif nm == "Z_IF":
                    control, _ = self.outcomes[op.record]
                    for g in self.gens:
                        if g.x >> q & 1:
                            g.s ^= control
                elif nm == "CX":
                    self.cx(op.qs[0], op.qs[1])
                elif nm == "R":
                    self._measure(0, 1 << q, True)
                elif nm == "RX":
                    self._measure(1 << q, 0, True)
                elif nm in ("M", "MX"):
                    px, pz = (0, 1 << q) if nm == "M" else (1 << q, 0)
                    form, random = self._measure(px, pz, False)
                    if random:
                        self.var_meas[self.nvars] = len(self.outcomes)
                    self.outcomes.append((form, random))
        return self

    # -- results --------------------------------------------------------
    def detectors(self):
        """[(measurement indices, expected parity)]: one per deterministic outcome."""
        out = []
        for j, (form, random) in enumerate(self.outcomes):
            if random:
                continue
            idx = [j]
            for v in range(1, self.nvars + 1):
                if (form >> v) & 1:
                    idx.append(self.var_meas[v])
            out.append((sorted(idx), form & 1))
        return out

    def parity_form(self, meas_indices):
        """Affine form of a parity of outcomes; deterministic iff no variable survives."""
        f = 0
        for j in meas_indices:
            f ^= self.outcomes[j][0]
        return f

    def output_group(self, out_mask):
        """Generators (x, z, sign form) of the subgroup of the final stabilizer group
        supported on the output qubits: the stabilizer group of the output state when
        the other qubits are in a product state with it."""
        rows = [Gen(g.x, g.z, g.s) for g in self.gens]
        rest = ((1 << self.n) - 1) & ~out_mask
        for q in range(self.n):
            if not (rest >> q) & 1:
                continue
            b = 1 << q
            for part in ("x", "z"):
                piv = next((i for i, g in enumerate(rows) if getattr(g, part) & b), None)
                if piv is None:
                    continue
                p = rows.pop(piv)
                for g in rows:
                    if getattr(g, part) & b:
                        g.mul(p)
        return [g for g in rows if g.x or g.z]

    def frame_mask(self, form):
        """Measurement indices whose reported outcome enters a sign form."""
        m = 0
        for v in range(1, self.nvars + 1):
            if (form >> v) & 1:
                m |= 1 << self.var_meas[v]
        return m
