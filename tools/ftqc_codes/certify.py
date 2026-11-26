"""Code-facts receipts: build them (leg A here; legs B and C attach evidence)
and re-verify them.

A receipt separates three kinds of statement about a quantity q:
  upper bound   a witness word, checked by dot products only;
  lower bound   a proof by some leg, recorded in "evidence" with what ran;
  value         reported only when a proven lower bound meets a witnessed upper bound.
Nothing is reported as proven that no leg proved.
"""
import platform
import sys

from . import gf2, minweight
from .codes import CSSCode

SCHEMA = "ftqc_codes-receipt-1"
LEGACY_SCHEMAS = ("qeccert-receipt-1",)
VERSION = "0.1.0"


class Refusal(Exception):
    """The input is not a CSS code, or a receipt fails verification."""


# The four minimum-weight quantities.  Each is  min wt (span(G) \ span(E)).
def _quantities(code):
    n = code.n
    return {
        # X-type logicals: ker(H_Z) outside rowspace(H_X)
        "d_x": (gf2.nullspace(code.hz, n), code.hx),
        "d_z": (gf2.nullspace(code.hx, n), code.hz),
        # lightest nonzero stabilizer of each type
        "wmin_x": (code.hx, []),
        "wmin_z": (code.hz, []),
    }


def _bits(v, n):
    return "".join("1" if (v >> j) & 1 else "0" for j in range(n))


def _unbits(s):
    return int(s[::-1], 2) if s else 0


def _partner(word, dual_basis):
    """A vector of dual_basis (or a sum of two) with odd overlap with word."""
    for z in dual_basis:
        if gf2.dot(z, word):
            return z
    return None


def structure(code):
    if not code.commutes():
        raise Refusal("H_X and H_Z do not commute: some X check overlaps some Z check oddly")
    n = code.n
    ex, ez = gf2.Echelon(code.hx), gf2.Echelon(code.hz)
    k = n - ex.rank - ez.rank
    self_dual = ex.rank == ez.rank and all(r in ez for r in ex.basis)

    def doubly_even(rows):
        b = gf2.Echelon(rows).basis
        return all(gf2.wt(r) % 4 == 0 for r in b) and \
            all(gf2.wt(r & s) % 2 == 0 for r in b for s in b)
    return dict(n=n, k=k, rank_hx=ex.rank, rank_hz=ez.rank,
                rows_hx=len(code.hx), rows_hz=len(code.hz),
                self_dual=self_dual,
                doubly_even_x=doubly_even(code.hx), doubly_even_z=doubly_even(code.hz),
                all_ones_logical=_all_ones_logical(code, ex, ez))


def _all_ones_logical(code, ex, ez):
    """True when 1^n is an X-type and a Z-type logical (transversal X and Z act
    as logical operators), as the self-dual constructions in the manuscripts need."""
    ones = (1 << code.n) - 1
    in_kz = all(gf2.dot(ones, h) == 0 for h in code.hz)
    in_kx = all(gf2.dot(ones, h) == 0 for h in code.hx)
    return bool(in_kz and in_kx and ones not in ex and ones not in ez)


def leg_a(code, targets=None, max_combos=5_000_000):
    """Exact enumeration.  targets: {"d_x": 7, ...} stops each search once the
    claimed value is proven as a lower bound; absent keys run to exactness."""
    targets = targets or {}
    out, evidence = {}, []
    n = code.n
    for q, (gens, excl) in _quantities(code).items():
        entry = dict(lower=None, upper=None, exact=False, witness=None, partner=None)
        ev = dict(leg="A", quantity=q, method="information-set enumeration (Brouwer-Zimmermann)")
        r = minweight.min_weight(gens, n, excl, target=targets.get(q), max_combos=max_combos)
        if r.get("empty"):
            ev.update(status="empty", reason="no such word exists (k = 0 or no stabilizers)")
            evidence.append(ev)
            out[q] = entry
            continue
        entry.update(lower=r["lower"], upper=r["upper"], exact=r["exact"])
        if r["witness"] is not None:
            entry["witness"] = _bits(r["witness"], n)
            if q in ("d_x", "d_z"):
                dual = gf2.nullspace(code.hx if q == "d_x" else code.hz, n)
                p = _partner(r["witness"], dual)
                entry["partner"] = _bits(p, n) if p is not None else None
        ev.update(status="proved", claim=f"{q} >= {r['lower']}", depth=r["depth"],
                  combos=r["combos"], info_sets=r["info_sets"])
        if r["budget"]:
            ev["budget"] = r["budget"]
        evidence.append(ev)
        out[q] = entry
    return out, evidence


def _combine(facts):
    """Derived quantities from d_x, d_z."""
    dx, dz = facts["d_x"], facts["d_z"]
    lo = None if dx["lower"] is None or dz["lower"] is None else min(dx["lower"], dz["lower"])
    ups = [u for u in (dx["upper"], dz["upper"]) if u is not None]
    up = min(ups) if ups else None
    facts["d"] = dict(lower=lo, upper=up, exact=(lo is not None and lo == up))
    facts["t"] = None if lo is None else (lo - 1) // 2
    return facts


def build(code, targets=None, max_combos=5_000_000, legs=("A",)):
    facts = structure(code)
    evidence = []
    if "A" in legs:
        mw, evidence = leg_a(code, targets, max_combos)
    else:
        mw = {q: dict(lower=None, upper=None, exact=False, witness=None, partner=None)
              for q in _quantities(code)}
    facts.update(mw)
    _combine(facts)
    return dict(schema=SCHEMA, code=code.to_json(), facts=facts, evidence=evidence,
                tool=dict(name="ftqc_codes", version=VERSION,
                          python=platform.python_version(), argv=sys.argv[1:]))


def add_lower_bound(receipt, quantity, bound, evidence):
    """Attach a lower bound proved by another leg.  A bound above a witnessed
    upper bound is a contradiction between legs and is refused."""
    e = receipt["facts"][quantity]
    if e["upper"] is not None and bound > e["upper"]:
        raise Refusal(f"leg {evidence.get('leg')} claims {quantity} >= {bound} but a witness "
                      f"of weight {e['upper']} is on record: the legs disagree")
    if e["lower"] is None or bound > e["lower"]:
        e["lower"] = bound
    e["exact"] = e["upper"] is not None and e["lower"] == e["upper"]
    receipt["evidence"].append(evidence)
    _combine(receipt["facts"])
    return receipt


# ---------------------------------------------------------------- verification
def verify(receipt, full=False, max_combos=5_000_000):
    """Re-check a receipt.  Cheap mode: hash, commutation, ranks, flags, every
    upper-bound witness, and internal consistency.  Full mode also re-proves
    every leg-A lower bound.  Returns a list of check lines; raises Refusal on
    the first failure."""
    if receipt.get("schema") not in (SCHEMA, *LEGACY_SCHEMAS):
        raise Refusal("unknown receipt schema")
    try:
        code = CSSCode.from_json(receipt["code"])
    except ValueError as e:
        raise Refusal(f"code block invalid: {e}") from None
    if "sha256" not in receipt["code"]:
        raise Refusal("code block carries no sha256")
    lines = [f"code sha256 {code.sha256()[:16]} matches matrices"]
    f = receipt["facts"]
    s = structure(code)
    for key, val in s.items():
        if f.get(key) != val:
            raise Refusal(f"fact {key}: receipt says {f.get(key)!r}, recomputed {val!r}")
    lines.append(f"structure recomputed: n={s['n']} k={s['k']} ranks {s['rank_hx']},{s['rank_hz']}")
    n = code.n
    quantities = _quantities(code)
    for q in quantities:
        e = f[q]
        lo, up = e["lower"], e["upper"]
        if up is not None:
            if not e["witness"]:
                raise Refusal(f"{q}: upper bound {up} without a witness")
            w = _unbits(e["witness"])
            if len(e["witness"]) != n or gf2.wt(w) != up or w == 0:
                raise Refusal(f"{q}: witness does not have weight {up}")
            if q in ("d_x", "d_z"):
                checks = code.hz if q == "d_x" else code.hx
                other = code.hx if q == "d_x" else code.hz
                if any(gf2.dot(w, h) for h in checks):
                    raise Refusal(f"{q}: witness fails a parity check")
                if not e["partner"]:
                    raise Refusal(f"{q}: witness has no partner proving it is not a stabilizer")
                p = _unbits(e["partner"])
                if any(gf2.dot(p, h) for h in other) or gf2.dot(p, w) != 1:
                    raise Refusal(f"{q}: partner does not certify a nontrivial logical")
            else:
                rows = code.hx if q == "wmin_x" else code.hz
                if w not in gf2.Echelon(rows):
                    raise Refusal(f"{q}: witness is not a stabilizer")
            lines.append(f"{q} <= {up}: witness verified")
        if lo is not None and up is not None and lo > up:
            raise Refusal(f"{q}: lower bound {lo} exceeds upper bound {up}")
        if e["exact"] != (lo is not None and lo == up):
            raise Refusal(f"{q}: 'exact' flag inconsistent with its bounds")
        proved = [ev for ev in receipt["evidence"]
                  if ev.get("quantity") == q and ev.get("status") == "proved"]
        claimed = max([int(ev["claim"].split(">=")[1]) for ev in proved], default=None)
        if lo is not None and (claimed is None or claimed < lo):
            raise Refusal(f"{q}: lower bound {lo} is not backed by any evidence entry")
        for ev in proved:
            if ev.get("leg") == "E":      # cut certificates are cheap to re-check
                from . import cuts
                m = cuts.verify_cuts(code, q, ev.get("cuts") or [])
                if m < int(ev["claim"].split(">=")[1]):
                    raise Refusal(f"{q}: cut certificate has {m} sets, fewer than it claims")
                lines.append(f"{q} >= {m}: cut certificate re-checked")
            if ev.get("leg") == "F":      # concatenation: structure and arithmetic
                from . import concat
                do, di = concat.recheck(code, q, ev)
                lines.append(f"{q} >= {do} * {di}: concatenation structure re-checked")
        if full and lo is not None:
            gens, excl = quantities[q]
            a_claims = [int(ev["claim"].split(">=")[1]) for ev in proved if ev.get("leg") == "A"]
            if a_claims:
                tgt = max(a_claims)
                r = minweight.min_weight(gens, n, excl, target=tgt, max_combos=max_combos)
                if r["lower"] < tgt:
                    raise Refusal(f"{q}: leg A re-run proves only {r['lower']}, receipt says {tgt}")
                lines.append(f"{q} >= {tgt}: leg A re-proved ({r['combos']} combinations)")
    want = dict(f)
    _combine(want)
    d = want["d"]
    if f["d"] != d or f["t"] != want["t"]:
        raise Refusal("derived d or t inconsistent with d_x, d_z")
    lines.append(f"d in [{d['lower']}, {d['upper']}], t = {f['t']}")
    return lines


def summary(receipt):
    f = receipt["facts"]

    def rng(e):
        if e["lower"] is None and e["upper"] is None:
            return "n/a"
        if e["exact"]:
            return str(e["lower"])
        return f"[{e['lower'] if e['lower'] is not None else '?'}, " \
               f"{e['upper'] if e['upper'] is not None else '?'}]"
    legs = sorted({ev["leg"] for ev in receipt["evidence"] if ev.get("status") == "proved"})
    return (f"{receipt['code']['name']}: n={f['n']} k={f['k']} d={rng(f['d'])} "
            f"(d_x={rng(f['d_x'])}, d_z={rng(f['d_z'])}) t={f['t']} "
            f"wmin_x={rng(f['wmin_x'])} wmin_z={rng(f['wmin_z'])} "
            f"self_dual={f['self_dual']} legs={','.join(legs) or 'none'}")
