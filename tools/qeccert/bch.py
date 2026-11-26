"""Leg D: the BCH bound for cyclic codes, with the zeros computed exactly.

For a binary cyclic code C of odd length n with generator polynomial g, let
alpha be a primitive n-th root of unity in GF(2^m) and Z = {j : g(alpha^j) = 0}.
BCH bound: if Z contains b, b+s, ..., b+(delta-2)s (mod n) for some step s
coprime to n, then every nonzero word of C has weight at least delta.

What this module computes exactly: that C is cyclic, its generator
polynomial, the field, an element of order n, and the zero set.  What it takes
from the literature: the BCH theorem itself (Bose and Ray-Chaudhuri 1960,
Hocquenghem 1959; the step-s form follows by replacing alpha with alpha^s).
A lower bound from this leg is therefore "exact zeros + cited theorem", and
the evidence entry says so.

Applied to a CSS code: d_x >= bound(ker H_Z), d_z >= bound(ker H_X),
wmin_x >= bound(rowspace H_X), wmin_z >= bound(rowspace H_Z).
"""
from math import gcd

from . import certify, gf2


class NotCyclic(Exception):
    pass


# ---------------------------------------------------------------- GF(2)[x]
def _deg(p):
    return p.bit_length() - 1


def _pmod(a, b):
    db = _deg(b)
    while a and _deg(a) >= db:
        a ^= b << (_deg(a) - db)
    return a


def _pmulmod(a, b, mod):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if _deg(a) >= _deg(mod):
            a ^= mod
    return r


def _ppowmod(a, e, mod):
    r = 1
    while e:
        if e & 1:
            r = _pmulmod(r, a, mod)
        a = _pmulmod(a, a, mod)
        e >>= 1
    return r


def _pgcd(a, b):
    while b:
        a, b = b, _pmod(a, b)
    return a


def _prime_factors(n):
    out, q = [], 2
    while q * q <= n:
        if n % q == 0:
            out.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        out.append(n)
    return out


def _irreducible(p):
    """Rabin's test for p in GF(2)[x]."""
    m = _deg(p)
    x = 2
    if _ppowmod(x, 1 << m, p) != _pmod(x, p):
        return False
    for q in _prime_factors(m):
        if _pgcd(_ppowmod(x, 1 << (m // q), p) ^ x, p) != 1:
            return False
    return True


def find_irreducible(m):
    if m == 1:
        return 0b11
    for mid in range(1, 1 << (m - 1)):
        p = (1 << m) | (mid << 1) | 1
        if _irreducible(p):
            return p
    raise AssertionError("no irreducible polynomial found")


def mult_order_2(n):
    m, v = 1, 2 % n
    while v != 1:
        v = (v * 2) % n
        m += 1
    return m


def root_of_unity(n):
    """(field polynomial, alpha) with alpha of exact order n in GF(2^m)."""
    if n % 2 == 0 or n < 3:
        raise NotCyclic("BCH leg needs odd length n >= 3")
    m = mult_order_2(n)
    p = find_irreducible(m)
    cof = ((1 << m) - 1) // n
    for beta in range(2, 1 << m):
        a = _ppowmod(beta, cof, p)
        if a != 1 and all(_ppowmod(a, n // q, p) != 1 for q in _prime_factors(n)):
            assert _ppowmod(a, n, p) == 1
            return p, a
    raise AssertionError("no element of order n found")


# ---------------------------------------------------------------- cyclic codes
def generator_poly(basis, n):
    """Generator polynomial of span(basis) if that code is cyclic."""
    ech = gf2.Echelon(basis)
    if ech.rank == 0:
        raise NotCyclic("zero code")
    mask = (1 << n) - 1
    for b in ech.basis:
        if (((b << 1) | (b >> (n - 1))) & mask) not in ech:
            raise NotCyclic("code is not cyclic in the given coordinate order")
    g = min(ech.basis, key=lambda b: b.bit_length())
    if _deg(g) != n - ech.rank:
        raise NotCyclic("minimal-degree word does not generate the code")
    if _pmod((1 << n) | 1, g) != 0:
        raise NotCyclic("generator does not divide x^n - 1")
    return g


def zeros(g, n, field=None):
    p, alpha = field or root_of_unity(n)
    Z, aj = [], 1
    for j in range(n):
        # Horner evaluation of g at alpha^j in GF(2^m)
        acc = 0
        for i in range(_deg(g), -1, -1):
            acc = _pmulmod(acc, aj, p) ^ ((g >> i) & 1)
        if acc == 0:
            Z.append(j)
        aj = _pmulmod(aj, alpha, p)
    if len(Z) != _deg(g):
        raise AssertionError("zero count differs from deg g")
    return Z, p, alpha


def bch_bound(Z, n, max_steps=None):
    """Best delta over all starts b and steps s coprime to n, with the run."""
    zs = set(Z)
    best = (1, 0, 1, 0)        # delta, start, step, run length
    steps = [s for s in range(1, n) if gcd(s, n) == 1]
    if max_steps:
        steps = steps[:max_steps]
    for s in steps:
        for b in Z:
            if (b - s) % n in zs:
                continue        # not the start of a run
            run, j = 0, b
            while j in zs and run < n:
                run += 1
                j = (j + s) % n
            if run + 1 > best[0]:
                best = (run + 1, b, s, run)
    return best


def code_bound(basis, n):
    g = generator_poly(basis, n)
    Z, p, alpha = zeros(g, n)
    delta, b, s, run = bch_bound(Z, n, max_steps=None if n <= 2000 else 1)
    leaders = sorted({min((z * pow(2, i, n)) % n for i in range(mult_order_2(n))) for z in Z})
    return dict(bound=delta, generator=hex(g), field_poly=hex(p), alpha=hex(alpha),
                zero_coset_leaders=leaders, zeros=len(Z), run_start=b, run_step=s, run_length=run)


def run_leg(receipt, code):
    n = code.n
    spaces = {"d_x": gf2.nullspace(code.hz, n), "d_z": gf2.nullspace(code.hx, n),
              "wmin_x": code.hx, "wmin_z": code.hz}
    for q, basis in spaces.items():
        ev = dict(leg="D", quantity=q,
                  method="BCH bound: zeros computed exactly, theorem cited (not machine-checked)")
        if receipt["facts"][q]["upper"] is None and receipt["facts"]["k"] == 0 and q in ("d_x", "d_z"):
            ev.update(status="skipped", reason="k = 0")
            receipt["evidence"].append(ev)
            continue
        try:
            r = code_bound(basis, n)
        except NotCyclic as e:
            ev.update(status="skipped", reason=str(e))
            receipt["evidence"].append(ev)
            continue
        b = r.pop("bound")
        ev.update(status="proved", claim=f"{q} >= {b}", **r)
        certify.add_lower_bound(receipt, q, b, ev)
    return receipt


def cyclic_css(n, leaders, name=None):
    """Self-dual CSS code from the cyclic code C whose zeros are the cyclotomic
    cosets of `leaders` (needs Z and -Z disjoint, so that C^perp is inside C):
    H_X = H_Z = a basis of C^perp, k = n - 2|Z|."""
    from .codes import CSSCode, SpecError, cyclic_span
    n = int(n)
    p, alpha = root_of_unity(n)
    m = mult_order_2(n)
    Z = set()
    for s in leaders:
        Z |= {(int(s) * pow(2, i, n)) % n for i in range(m)}
    if Z & {(-z) % n for z in Z} or 0 in Z:
        raise SpecError("zero set meets its negative: C^perp is not inside C")
    # g(x) = prod (x + alpha^z), coefficients in GF(2^m), must land in GF(2)
    coeffs = [1]
    for z in sorted(Z):
        a = _ppowmod(alpha, z, p)
        nxt = [0] * (len(coeffs) + 1)
        for i, c in enumerate(coeffs):
            nxt[i + 1] ^= c
            nxt[i] ^= _pmulmod(c, a, p)
        coeffs = nxt
    if any(c not in (0, 1) for c in coeffs):
        raise AssertionError("generator polynomial is not binary")
    g = sum(c << i for i, c in enumerate(coeffs))
    C = gf2.Echelon(cyclic_span(g, n))
    Cd = gf2.nullspace(C.basis, n)
    if not all(h in C for h in Cd):
        raise SpecError("C^perp is not inside C")
    rows = sorted(gf2.Echelon(Cd).basis)
    return CSSCode(name or f"Cyclic(n={n},zeros={sorted(int(s) for s in leaders)})", n, rows, rows,
                   dict(family="cyclic", n=n, leaders=sorted(int(s) for s in leaders)))
