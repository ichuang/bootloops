"""Receipts for circuit-level fault properties, and their re-verification.

As in ftqc_codes, a number is an interval: the upper end carries a witness (a fault
set that `verify` replays through the circuit), the lower end names the leg that
proved it.  A strict-fault-tolerance failure carries its counterexample.
"""
import os
import platform
import sys

from . import analysis
from .circuit import Circuit
from .model import Model

SCHEMA = "ftqc_audit-receipt-1"
VERSION = "0.2.0"
NOISE = ("one fault = one location: any non-identity Pauli after R, RX, I, H, S, SDG, Z_IF (one qubit) or "
         "CX (two qubits), or a flipped outcome at M, MX; weight = number of locations")


class Refusal(Exception):
    pass


def _ftqc_codes():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if here not in sys.path:
        sys.path.insert(0, here)
    from ftqc_codes import cnf, sat
    return cnf, sat


def _problem(model, kind):
    """Reduced columns (det, payload) for a distance question, the harm predicate,
    and for each reduced column one model column that realizes it."""
    index = {}
    for i, (d, o, s) in enumerate(model.columns):
        pay = (o,) if kind == "circuit_distance" else (o, s)
        if d == 0 and not any(pay):
            continue
        index.setdefault((d,) + pay, i)
    cols = sorted(index)
    return cols, (lambda p: any(p)), [index[c] for c in cols]


def _faults(model, col_indices):
    """Fault set realizing a sum of model columns: one location each; two columns from
    the same location are merged into that location's combined Pauli."""
    by_loc = {}
    for ci in col_indices:
        loc, label = model.column_source[ci]
        by_loc.setdefault(loc, []).append(label)
    out = []
    for loc, labels in sorted(by_loc.items()):
        if len(labels) == 1:
            out.append((loc, labels[0]))
            continue
        acc = [0] * len(labels[0])
        code = {"I": 0, "X": 1, "Z": 2, "Y": 3}
        for lab in labels:
            acc = [a ^ code[c] for a, c in zip(acc, lab)]
        lab = "".join("IXZY"[a] for a in acc)
        if set(lab) != {"I"}:
            out.append((loc, lab))
    return out


def distance(model, kind, max_weight, legs=("A",), cap=3_000_000, workdir=None, sat_timeout=900):
    """kind: 'circuit_distance' (harm = a logical observable flips) or
    'fault_distance' (harm = an observable flips or the output state is wrong)."""
    cols, harmful, back = _problem(model, kind)
    res = dict(lower=None, upper=None, exact=False, witness=None, evidence=[])
    if kind == "circuit_distance" and not model.obs_masks:
        res["evidence"].append(dict(leg="-", status="skipped", reason="circuit declares no observable"))
        return res
    if not cols:
        res["evidence"].append(dict(leg="-", status="empty",
                                    reason="no fault can cause harm without firing a detector alone"))
    if "A" in legs:
        r = analysis.min_undetected(cols, harmful, max_weight, cap)
        ev = dict(leg="A", method="exact enumeration, meet in the middle on detector signatures",
                  columns=len(cols))
        if r["witness"] is not None:
            faults = _faults(model, [back[i] for i in r["witness"]])
            _set_upper(model, kind, res, faults)
        if r["lower"] is not None:
            res["lower"] = r["lower"]
            ev.update(status="proved", claim=f">= {r['lower']}")
        else:
            ev.update(status="empty", reason=r["note"])
        if r["note"]:
            ev["note"] = r["note"]
        res["evidence"].append(ev)
    if "B" in legs:
        _leg_b(model, kind, cols, back, res, max_weight, workdir, sat_timeout)
    res["exact"] = res["lower"] is not None and res["lower"] == res["upper"]
    return res


def _set_upper(model, kind, res, faults):
    det, obs, syn = model.replay(faults)
    harm = obs if kind == "circuit_distance" else (obs or syn)
    if det or not harm:
        raise Refusal("internal: a candidate witness does not replay as an undetected harmful fault set")
    if res["upper"] is None or len(faults) < res["upper"]:
        res["upper"] = len(faults)
        res["witness"] = model.describe(faults)
    if res["lower"] is not None and res["upper"] < res["lower"]:
        raise Refusal(f"a replayed fault set of weight {res['upper']} contradicts the proven "
                      f"lower bound {res['lower']}: the legs disagree")


def _leg_b(model, kind, cols, back, res, max_weight, workdir, timeout):
    import tempfile
    cnf, sat = _ftqc_codes()
    base = dict(leg="B", method="SAT on the fault matrix (ftqc_codes.cnf), LRAT checked by ftqc_codes.lrat",
                solver=sat.solver_version())
    if not cols:
        return
    m = len(cols)
    ndet = max((c[0].bit_length() for c in cols), default=0)
    P = [sum(1 << i for i, c in enumerate(cols) if (c[0] >> b) & 1) for b in range(ndet)]
    D = []
    for part in range(1, len(cols[0])):
        nb = max((c[part].bit_length() for c in cols), default=0)
        D += [sum(1 << i for i, c in enumerate(cols) if (c[part] >> b) & 1) for b in range(nb)]
    tmp = None
    if workdir is None:
        tmp = tempfile.mkdtemp(prefix="ftqc_audit-sat-")
        workdir = tmp
    os.makedirs(workdir, exist_ok=True)
    try:
        # prove the bound just below the witnessed upper bound, or walk upward
        if res["upper"] is not None:
            todo = [res["upper"] - 1]
        else:
            todo = list(range(res["lower"] or 1, max_weight + 1))
        for w in todo:
            if w < 0:
                continue
            F = cnf.Formula(m, P, D, w)
            out = sat.solve(F, workdir, f"{model.circuit.sha256()[:12]}_{kind}_le{w}", timeout)
            if out[0] == "unsat":
                ev = dict(base, status="proved", claim=f">= {w + 1}", **out[1])
                res["evidence"].append(ev)
                res["lower"] = max(res["lower"] or 0, w + 1)
                if res["upper"] is not None and res["lower"] > res["upper"]:
                    raise Refusal("SAT proves a bound above a replayed witness: the legs disagree")
            elif out[0] == "sat":
                word = out[1]
                faults = _faults(model, [back[i] for i in range(m) if (word >> i) & 1])
                _set_upper(model, kind, res, faults)
                res["evidence"].append(dict(base, status="witness",
                                            reason=f"found a fault set of weight {len(faults)}", **out[2]))
                break
            else:
                res["evidence"].append(dict(base, status="skipped", reason=out[1]))
                break
    finally:
        if tmp:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)


def build(circuit, max_weight=4, strict_order=0, legs=("A",), cap=3_000_000, workdir=None,
          sat_timeout=900, want=("circuit_distance", "fault_distance", "census"), heavier=()):
    model = Model(circuit)
    rec = dict(schema=SCHEMA, circuit=circuit.to_json(), noise_model=NOISE,
               counts=circuit.counts(), model=model.stats(),
               detectors=model.detectors, results={},
               tool=dict(name="ftqc_audit", version=VERSION, python=platform.python_version()))
    for kind in ("circuit_distance", "fault_distance"):
        if kind in want:
            rec["results"][kind] = distance(model, kind, max_weight, legs, cap, workdir, sat_timeout)
    if strict_order:
        r = analysis.strict_ft(model, strict_order, cap)
        if r["violation"]:
            faults = _faults(model, r["violation"]["columns"])
            r["violation"] = dict(order=r["violation"]["order"], faults=model.describe(faults))
        rec["results"]["strict_ft"] = r
    if heavier:
        rows = []
        for w in heavier:
            r = analysis.min_heavier(model, w, max_weight, cap)
            row = dict(weight=w, lower=r["lower"], upper=r["upper"], exact=r["exact"], note=r["note"],
                       witness=None)
            if r["witness"] is not None:
                row["witness"] = model.describe(_faults(model, r["witness"]))
            rows.append(row)
        rec["results"]["heavier_than"] = rows
    if "census" in want:
        rec["results"]["census"] = analysis.census(model)
    return rec, model


def verify(receipt, full=False, cap=3_000_000):
    if receipt.get("schema") not in (SCHEMA, "faultcert-receipt-1"):
        raise Refusal("unknown receipt schema")
    try:
        circuit = Circuit.from_json(receipt["circuit"])
    except ValueError as e:
        raise Refusal(f"circuit block invalid: {e}") from None
    if "sha256" not in receipt["circuit"]:
        raise Refusal("circuit block carries no sha256")
    model = Model(circuit)
    lines = [f"circuit sha256 {circuit.sha256()[:16]} matches"]
    if receipt["model"] != model.stats() or receipt["counts"] != circuit.counts():
        raise Refusal("model statistics differ from the recomputed ones")
    if receipt["detectors"] != model.detectors:
        raise Refusal("detector list differs from the recomputed one")
    lines.append("fault model recomputed: " + ", ".join(f"{k}={v}" for k, v in model.stats().items()))
    for kind in ("circuit_distance", "fault_distance"):
        r = receipt["results"].get(kind)
        if not r:
            continue
        lo, up = r["lower"], r["upper"]
        if up is not None:
            faults = [(w["location"], w["fault"]) for w in r["witness"] or []]
            if len(faults) != up:
                raise Refusal(f"{kind}: witness does not have {up} faults")
            try:
                det, obs, syn = model.replay(faults)
            except ValueError as e:
                raise Refusal(f"{kind}: witness invalid: {e}") from None
            harm = obs if kind == "circuit_distance" else (obs or syn)
            if det or not harm:
                raise Refusal(f"{kind}: witness fires a detector or does no harm")
            lines.append(f"{kind} <= {up}: witness replayed")
        if lo is not None and up is not None and lo > up:
            raise Refusal(f"{kind}: lower bound {lo} exceeds upper bound {up}")
        if r["exact"] != (lo is not None and lo == up):
            raise Refusal(f"{kind}: 'exact' flag inconsistent")
        claims = [int(e["claim"].split(">=")[1]) for e in r["evidence"] if e.get("status") == "proved"]
        if lo is not None and (not claims or max(claims) < lo):
            raise Refusal(f"{kind}: lower bound {lo} is not backed by any evidence entry")
        if full and lo is not None and lo > 1:
            cols, harmful, _ = _problem(model, kind)
            again = analysis.min_undetected(cols, harmful, lo - 1, cap)
            if again["witness"] is not None or again["lower"] < lo:
                raise Refusal(f"{kind}: re-enumeration does not confirm >= {lo}")
            lines.append(f"{kind} >= {lo}: re-proved by enumeration")
    s = receipt["results"].get("strict_ft")
    if s and s.get("violation"):
        v = s["violation"]
        faults = [(w["location"], w["fault"]) for w in v["faults"]]
        det, obs, syn = model.replay(faults)
        if det or len(faults) > v["order"] or syn in analysis.weight_table(model, v["order"]):
            raise Refusal("strict_ft: the recorded violation does not replay as one")
        lines.append(f"strict fault tolerance fails at order {v['order']}: counterexample replayed")
    if s and full:
        again = analysis.strict_ft(model, max([p["order"] for p in s["per_order"]], default=0), cap)
        if again["holds_through"] != s["holds_through"] or again["per_order"] != s["per_order"]:
            raise Refusal("strict_ft: re-enumeration gives different counts")
        lines.append(f"strict fault tolerance through order {s['holds_through']}: re-enumerated")
    for row in receipt["results"].get("heavier_than", []):
        if row["upper"] is not None:
            faults = [(w["location"], w["fault"]) for w in row["witness"] or []]
            det, obs, syn = model.replay(faults)
            if det or len(faults) != row["upper"] or syn in analysis.weight_table(model, row["weight"]):
                raise Refusal(f"heavier_than {row['weight']}: witness does not replay as claimed")
            lines.append(f"accepted error heavier than {row['weight']} with {row['upper']} faults: replayed")
        if full and row["lower"] is not None and row["lower"] > 1:
            again = analysis.min_heavier(model, row["weight"], row["lower"] - 1, cap)
            if again["witness"] is not None or again["lower"] < row["lower"]:
                raise Refusal(f"heavier_than {row['weight']}: re-enumeration does not confirm the bound")
    c = receipt["results"].get("census")
    if c and c != analysis.census(model):
        raise Refusal("census differs from the recomputed one")
    return lines


def summary(receipt):
    def rng(r):
        if r is None or (r["lower"] is None and r["upper"] is None):
            return "n/a"
        if r["exact"]:
            return str(r["lower"])
        return f"[{r['lower']}, {r['upper'] if r['upper'] is not None else '?'}]"
    res = receipt["results"]
    parts = [f"{receipt['circuit']['name']}: {receipt['model']['locations']} locations, "
             f"{receipt['model']['detectors']} detectors"]
    if "circuit_distance" in res:
        parts.append(f"circuit distance {rng(res['circuit_distance'])}")
    if "fault_distance" in res:
        parts.append(f"fault distance {rng(res['fault_distance'])}")
    if "strict_ft" in res:
        s = res["strict_ft"]
        parts.append(f"strict FT through order {s['holds_through']}"
                     + (f" (fails at {s['violation']['order']})" if s["violation"] else ""))
    if "census" in res:
        parts.append(f"hook-capable locations {res['census']['hook_capable']}")
    return "; ".join(parts)
