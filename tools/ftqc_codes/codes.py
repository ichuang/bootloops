"""CSS code specifications and family constructors.

A CSS code is a pair of binary matrices (H_X, H_Z) with H_X H_Z^T = 0.  Rows
of H_X are supports of X-type stabilizer generators, rows of H_Z of Z-type.
X-type logical operators are the vectors of ker(H_Z) outside rowspace(H_X);
their least weight is d_X.  Z-type: ker(H_X) outside rowspace(H_Z), d_Z.
"""
import hashlib
import json
import os

from . import gf2


class SpecError(ValueError):
    pass


class CSSCode:
    def __init__(self, name, n, hx, hz, meta=None):
        self.name, self.n = name, int(n)
        self.hx, self.hz = list(hx), list(hz)
        self.meta = dict(meta or {})
        if self.n <= 0:
            raise SpecError("n must be positive")
        for r in self.hx + self.hz:
            if r < 0 or r >> self.n:
                raise SpecError("row has support outside 0..n-1")
        if not self.hx or not self.hz:
            raise SpecError("both H_X and H_Z need at least one row")

    # -- canonical text and hash -------------------------------------------
    def _bits(self, r):
        return "".join("1" if (r >> j) & 1 else "0" for j in range(self.n))

    def canonical(self):
        return (f"n={self.n}\nHX\n" + "\n".join(self._bits(r) for r in self.hx)
                + "\nHZ\n" + "\n".join(self._bits(r) for r in self.hz) + "\n")

    def sha256(self):
        return hashlib.sha256(self.canonical().encode()).hexdigest()

    def to_json(self):
        return dict(name=self.name, n=self.n, hx=[self._bits(r) for r in self.hx],
                    hz=[self._bits(r) for r in self.hz], meta=self.meta,
                    sha256=self.sha256())

    @staticmethod
    def from_json(d):
        try:
            n = int(d["n"])

            def rows(key):
                out = []
                for s in d[key]:
                    if isinstance(s, str):
                        if len(s) != n or set(s) - {"0", "1"}:
                            raise SpecError(f"{key} row is not a 0/1 string of length n")
                        out.append(int(s[::-1], 2))
                    else:
                        if len(s) != n:
                            raise SpecError(f"{key} row has wrong length")
                        out.append(gf2.rows_from_lists([s])[0])
                return out
            code = CSSCode(d.get("name", "unnamed"), n, rows("hx"), rows("hz"), d.get("meta"))
        except KeyError as e:
            raise SpecError(f"spec lacks field {e}") from None
        if "sha256" in d and d["sha256"] != code.sha256():
            raise SpecError("sha256 in spec does not match its matrices")
        return code

    def commutes(self):
        return all(gf2.dot(a, b) == 0 for a in self.hx for b in self.hz)


# ---------------------------------------------------------------- families
def _is_prime(p):
    return p > 1 and all(p % q for q in range(2, int(p ** 0.5) + 1))


def cyclic_span(word, n):
    mask = (1 << n) - 1
    return [((word << s) | (word >> (n - s))) & mask for s in range(n)]


def qr(p):
    """Self-dual CSS code [[p,1,d]] from the binary quadratic-residue code of
    prime length p = 7 mod 8: H_X = H_Z = a basis of the even-like QR code
    (dimension (p-1)/2), which is contained in its dual."""
    p = int(p)
    if p % 8 != 7 or not _is_prime(p):
        raise SpecError("qr family needs a prime p = 7 mod 8")
    Q = sorted({(i * i) % p for i in range(1, p)})
    N = [i for i in range(1, p) if i not in Q]
    for S in (Q, N):
        for const in (0, 1):
            e = sum(1 << r for r in S) | const
            C = gf2.Echelon(cyclic_span(e, p))
            if C.rank != (p + 1) // 2:
                continue
            Cd = gf2.nullspace(C.basis, p)
            if all(h in C for h in Cd) and all(gf2.wt(h) % 4 == 0 for h in gf2.Echelon(Cd).basis):
                rows = sorted(gf2.Echelon(Cd).basis)
                return CSSCode(f"QR{p}", p, rows, rows, dict(family="qr", p=p))
    raise SpecError("no doubly-even self-orthogonal QR subcode found")


def steane():
    c = qr(7)
    c.name = "Steane"
    return c


def golay():
    c = qr(23)
    c.name = "Golay"
    return c


def _shift(n, k):
    """n x n cyclic shift matrix as a list of rows: row i has a 1 at i+k."""
    return [1 << ((i + k) % n) for i in range(n)]


def _kron(A, na, B, nb):
    """Kronecker product of row lists (A: rows over na columns, B over nb)."""
    out = []
    for a in A:
        for b in B:
            v = 0
            for i in gf2.supp(a):
                v |= b << (i * nb)
            out.append(v)
    return out


def _eye(n):
    return [1 << i for i in range(n)]


def _add(A, B):
    return [a ^ b for a, b in zip(A, B)]


def _transpose(A, ncols):
    out = [0] * ncols
    for i, a in enumerate(A):
        for j in gf2.supp(a):
            out[j] |= 1 << i
    return out


def _hstack(A, B, na):
    return [a | (b << na) for a, b in zip(A, B)]


def bb(l, m, A, B, name=None):
    """Bivariate bicycle code (Bravyi et al. 2024).  A and B are lists of
    monomials such as ["x3", "y1", "y2"]; x = S_l (x) I_m, y = I_l (x) S_m.
    H_X = [A | B], H_Z = [B^T | A^T], n = 2 l m."""
    l, m = int(l), int(m)
    N = l * m

    def poly(terms):
        M = [0] * N
        for t in terms:
            t = t.strip().lower()
            if not t or t[0] not in "xy" or not t[1:].isdigit():
                raise SpecError(f"bad monomial {t!r}; use x<k> or y<k>")
            k = int(t[1:])
            mono = _kron(_shift(l, k), l, _eye(m), m) if t[0] == "x" else \
                _kron(_eye(l), l, _shift(m, k), m)
            M = _add(M, mono)
        return M
    Am, Bm = poly(A), poly(B)
    hx = _hstack(Am, Bm, N)
    hz = _hstack(_transpose(Bm, N), _transpose(Am, N), N)
    return CSSCode(name or f"BB(l={l},m={m})", 2 * N, hx, hz,
                   dict(family="bb", l=l, m=m, A=list(A), B=list(B)))


def hgp(H1, n1, H2, n2, name=None):
    """Hypergraph product of classical checks H1 (m1 x n1) and H2 (m2 x n2):
    H_X = [H1 (x) I_n2 | I_m1 (x) H2^T],  H_Z = [I_n1 (x) H2 | H1^T (x) I_m2]."""
    m1, m2 = len(H1), len(H2)
    left = n1 * n2
    hx = _hstack(_kron(H1, n1, _eye(n2), n2), _kron(_eye(m1), m1, _transpose(H2, n2), m2), left)
    hz = _hstack(_kron(_eye(n1), n1, H2, n2), _kron(_transpose(H1, n1), m1, _eye(m2), m2), left)
    return CSSCode(name or "HGP", left + m1 * m2, hx, hz, dict(family="hgp"))


def toric(L):
    L = int(L)
    if L < 2:
        raise SpecError("toric family needs L >= 2")
    H = [(1 << i) | (1 << ((i + 1) % L)) for i in range(L)]
    c = hgp(H, L, H, L, f"Toric(L={L})")
    c.meta = dict(family="toric", L=L)
    return c


def planar(L):
    """Unrotated planar surface code [[L^2 + (L-1)^2, 1, L]]."""
    L = int(L)
    if L < 2:
        raise SpecError("planar family needs L >= 2")
    H = [(1 << i) | (1 << (i + 1)) for i in range(L - 1)]
    c = hgp(H, L, H, L, f"Planar(L={L})")
    c.meta = dict(family="planar", L=L)
    return c


def rotated_plaquettes(d):
    """Yield (basis, corner coordinate, occupied corners) in syndrome-record order."""
    if d < 2:
        raise SpecError("rotated_surface family needs d >= 2")
    for i in range(d + 1):
        for j in range(d + 1):
            corners = {"NW": (i-1, j-1), "NE": (i-1, j), "SW": (i, j-1), "SE": (i, j)}
            cells = {k: r*d+c for k, (r,c) in corners.items() if 0 <= r < d and 0 <= c < d}
            axis = "Z" if (i+j) % 2 else "X"
            if len(cells) == 4 or len(cells) == 2 and ((i in (0, d)) == (axis == "X")):
                yield axis, (i, j), cells


def rotated_surface(d):
    """Rotated surface code [[d^2, 1, d]]: data qubits on a d x d grid, bulk
    plaquettes in a checkerboard, weight-2 X checks on the top and bottom
    edges and weight-2 Z checks on the left and right edges."""
    d = int(d)
    hx, hz = [], []
    for axis, _, cells in rotated_plaquettes(d):
        (hx if axis == "X" else hz).append(gf2.from_supp(cells.values()))
    return CSSCode(f"RotatedSurface(d={d})", d * d, hx, hz,
                   dict(family="rotated_surface", d=d))


def concatenate(outer, inner, name=None):
    """Concatenation of an outer CSS code with an inner CSS code that has k = 1."""
    from . import concat
    hx, hz = concat.expected_stabilizers(outer, inner)
    return CSSCode(name or f"{outer.name} o {inner.name}", outer.n * inner.n, hx, hz,
                   dict(family="concat", outer=outer.name, inner=inner.name))


def _radial90():
    from . import constructions
    return constructions.radial90()


def _lp6_2():
    from . import constructions
    return constructions.lifted_product([[3]], [[3]], 3, name="LP(6,2)")


def _pp184():
    from . import families_highrate
    return families_highrate.pp184_50()


def _pp256():
    from . import families_highrate
    return families_highrate.pp256_68()


def _pp576():
    from . import families_highrate
    return families_highrate.pp576_148()


def _apm1152():
    from . import families_highrate
    return families_highrate.apm1152_580()


def _qt_d8_36():
    from . import families_highrate
    return families_highrate.quantum_expanders_d8_36_8_3()


def _qt_s3_288():
    from . import families_highrate
    return families_highrate.quantum_expanders_s3_288_8()


FAMILIES = {
    "steane": (steane, "Steane [[7,1,3]]"),
    "golay": (golay, "Golay [[23,1,7]]"),
    "qr": (qr, "qr:P, self-dual quadratic-residue CSS code, prime P = 7 mod 8"),
    "toric": (toric, "toric:L, [[2L^2,2,L]]"),
    "planar": (planar, "planar:L, unrotated [[L^2+(L-1)^2,1,L]]"),
    "rotated_surface": (rotated_surface, "rotated_surface:D, [[D^2,1,D]]"),
    "cyclic": (None, "cyclic:N,l1,l2,..., self-dual CSS code from the cyclic code with zero "
                     "cosets l1, l2, ... (BCH codes: cyclic:127,1,3,5,7,9)"),
    "bb": (None, "bb:L,M,A,B with A,B like x3+y1+y2, e.g. bb:6,6,x3+y1+y2,y3+x1+x2"),
}

# Named instances.  The parameters of the bivariate and generalized bicycle
# codes are those of the Lean-QEC examples (LeanQEC/Stabilizer/Examples/BB/).
NAMED = {
    "bb18": lambda: bb(3, 3, ["x0", "x1", "y1"], ["x0", "x2", "y2"], "BB18"),
    "bb54": lambda: bb(3, 9, ["x0", "y2", "y4"], ["x1", "x2", "y3"], "BB54"),
    "gb54": lambda: bb(27, 1, ["x8", "x13", "x15", "x16"], ["x7", "x8", "x10", "x15"], "GB54"),
    "bb70": lambda: bb(7, 5, ["y2", "x2", "x3", "x4"], ["y1", "x1", "x3"], "BB70"),
    "bb72": lambda: bb(6, 6, ["x3", "y1", "y2"], ["y3", "x1", "x2"], "BB72"),
    "bb90": lambda: bb(15, 3, ["x9", "y1", "y2"], ["y0", "x2", "x7"], "BB90"),
    "bb108": lambda: bb(9, 6, ["x3", "y1", "y2"], ["y3", "x1", "x2"], "BB108"),
    "bb144": lambda: bb(12, 6, ["x3", "y1", "y2"], ["y3", "x1", "x2"], "BB144"),
    "radial90": _radial90,
    "lp6_2": _lp6_2,
    "pp184": _pp184,
    "pp256": _pp256,
    "pp576": _pp576,
    "apm1152": _apm1152,
    "qt_d8_36": _qt_d8_36,
    "qt_s3_288": _qt_s3_288,
}


def load(spec):
    """spec: path to a JSON code / import / construction file, or 'family[:args]'."""
    if spec.endswith(".json"):
        base_dir = os.path.dirname(os.path.abspath(spec))
        with open(spec) as f:
            d = json.load(f)
        if isinstance(d, dict) and d.get("schema") == "ftqc_codes-import-1":
            from . import interchange
            return interchange.load_import_manifest(d, base_dir=base_dir)
        if isinstance(d, dict) and d.get("schema") == "ftqc_codes-construction-spec-1":
            from . import constructions
            return constructions.from_spec(d, base_dir=base_dir)
        return CSSCode.from_json(d)
    fam, _, arg = spec.partition(":")
    fam = fam.lower()
    if fam in NAMED and not arg:
        return NAMED[fam]()
    if fam == "bb":
        parts = arg.split(",")
        if len(parts) != 4:
            raise SpecError("bb needs L,M,A,B")
        return bb(parts[0], parts[1], parts[2].split("+"), parts[3].split("+"))
    if fam == "cyclic":
        from . import bch
        parts = arg.split(",")
        if len(parts) < 2:
            raise SpecError("cyclic needs N,leader[,leader...]")
        return bch.cyclic_css(parts[0], parts[1:])
    if fam not in FAMILIES:
        raise SpecError(f"unknown family {fam!r}; known: {sorted(FAMILIES) + sorted(NAMED)}")
    fn = FAMILIES[fam][0]
    try:
        return fn(*([arg] if arg else []))
    except TypeError:
        raise SpecError(f"wrong arguments for {fam}: {FAMILIES[fam][1]}") from None
