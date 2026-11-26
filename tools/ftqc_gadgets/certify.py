"""Receipts for gadget properties, refinements and compositions.

A receipt records the gadget (circuit, blocks, decoders, with a hash), the size of
its fault model, and a list of results.  Every negative result carries a witness (an
input error and a fault set) that `verify` replays by propagation.  Positive results
are re-enumerated by `verify(full=True)`.

A composition receipt states what follows when preparations that a gadget was
analysed with as specifications are replaced by implementations that refine those
specifications; see `compose`.
"""
import hashlib
import json

from . import library, props
from .effects import Analysis
from .gadget import Gadget, GadgetError

VERSION = "0.1.0"
SCHEMA = "ftqc_gadgets-receipt-1"
COMPOSITION = "ftqc_gadgets-composition-1"


class Refusal(Exception):
    pass


def _clean(d):
    return json.loads(json.dumps(d))


def build(gadget, checks=(), close=None, refines=(), sector=False, cap=3_000_000, analysis=None,
          local=None, tables=True):
    """local: True, an integer radius or {decoder: radius} adds the coordinate-locality
    certificate (props.local).  tables=False skips decoder tables (large codes; only the
    locality certificate is then available).
    checks: [(kind, t)] or [(kind, t, sector)].  close: s or None.  refines: [(spec gadget, order, input weight)]."""
    an = analysis or Analysis(gadget, cap, tables)
    results = []
    for item in checks:
        kind, t = item[0], item[1]
        r = props.check(an, kind, int(t), item[2] if len(item) > 2 else sector)
        r["type"] = "check"
        results.append(r)
    if close is not None:
        r = props.close(an, int(close), sector)
        r.update(type="close", s=int(close), weight="sector" if sector else "pauli")
        results.append(r)
    for spec, order, rin in refines:
        r = props.refine(an, Analysis(spec, cap), int(order), int(rin))
        r.update(type="refine", spec=spec.to_json())
        results.append(r)
    if local is not None and local is not False:
        r = props.local(an, None if local is True else local)
        r["lemma"] = props.LEMMA
        results.append(r)
    return _clean(dict(schema=SCHEMA, tool="ftqc_gadgets", tool_version=VERSION,
                       gadget=gadget.to_json(), stats=an.stats(), tables=bool(tables), results=results))


def receipt_sha(receipt):
    return hashlib.sha256(json.dumps(receipt, sort_keys=True).encode()).hexdigest()


def _first_bad(per, key):
    return next((p for p in per if p[key]), None)


def verify(receipt, full=False, cap=3_000_000):
    """Raises Refusal; returns the list of things that were checked."""
    if receipt.get("schema") not in (SCHEMA, "gadgetcert-receipt-1"):
        raise Refusal(f"not a {SCHEMA} receipt")
    try:
        g = Gadget.from_json(receipt["gadget"])
        an = Analysis(g, cap, receipt.get("tables", True))
    except (GadgetError, ValueError) as e:
        raise Refusal(f"gadget: {e}") from None
    done = ["gadget hash, input preamble, decoder checks, detectors"]
    if _clean(an.stats()) != receipt.get("stats"):
        raise Refusal("the fault-model statistics differ from the recomputed ones")
    done.append("fault-model statistics")
    for i, r in enumerate(receipt.get("results", [])):
        tag = f"result {i} ({r.get('type')})"
        typ = r.get("type")
        sector = r.get("weight") == "sector"
        v = r.get("violation")
        try:
            if typ == "check":
                if v is not None:
                    if not props.judge(an, r["kind"], v, sector):
                        raise Refusal(f"{tag}: the witness is not a violation")
                    if r["holds_through"] >= v["r"] + v["s"]:
                        raise Refusal(f"{tag}: holds_through contradicts the witness")
                    done.append(f"{tag}: witness replayed")
                elif r["holds_through"] >= 0 and not r.get("per"):
                    raise Refusal(f"{tag}: a positive claim without evidence")
                rows = {(p["r"], p["s"]) for p in r.get("per", []) if not p["violations"]}
                rmax = 0 if r["kind"] == "prep" else r["t"]
                for tot in range(r["holds_through"] + 1):
                    for rr in range(min(tot, rmax) + 1):
                        if (rr, tot - rr) not in rows:
                            raise Refusal(f"{tag}: holds_through {r['holds_through']} is not "
                                          f"covered by the recorded cases")
                if full:
                    again = _clean(props.check(an, r["kind"], r["t"], sector))
                    for k in ("holds_through", "per", "note"):
                        if again[k] != r[k]:
                            raise Refusal(f"{tag}: re-enumeration gives a different {k}")
                    done.append(f"{tag}: re-enumerated")
            elif typ == "close":
                if v is not None:
                    sg = an.replay([tuple(f) for f in v["faults"]]) ^ v["input_signature"]
                    again = _clean(props.close(an, v["s"], sector))
                    if sg & an.DET or again["violation"] is None or again["holds_through"] >= v["s"]:
                        raise Refusal(f"{tag}: the witness is not a violation")
                    done.append(f"{tag}: witness order re-enumerated")
                if r["holds_through"] >= (v["s"] if v else r["s"] + 1):
                    raise Refusal(f"{tag}: holds_through exceeds what was examined")
                if full:
                    again = _clean(props.close(an, r["s"], sector))
                    for k in ("holds_through", "per", "method"):
                        if again[k] != r[k]:
                            raise Refusal(f"{tag}: re-enumeration gives a different {k}")
                    done.append(f"{tag}: re-enumerated")
            elif typ == "refine":
                spec = Analysis(Gadget.from_json(r["spec"]), cap)
                why = props.comparable(an, spec)
                if why:
                    raise Refusal(f"{tag}: {why}")
                if v is not None:
                    if not props.judge_refine(an, spec, v):
                        raise Refusal(f"{tag}: the witness is matched by the specification")
                    if r["holds_through"] >= v["order"]:
                        raise Refusal(f"{tag}: holds_through contradicts the witness")
                    done.append(f"{tag}: witness replayed against the specification")
                if r["holds_through"] > r["order"] or len(
                        [p for p in r.get("per", []) if not p["unmatched"]]) < r["holds_through"] + 1:
                    raise Refusal(f"{tag}: holds_through is not covered by the recorded orders")
                if full:
                    again = _clean(props.refine(an, spec, r["order"], r["input_weight"]))
                    for k in ("holds_through", "per"):
                        if again[k] != r[k]:
                            raise Refusal(f"{tag}: re-enumeration gives a different {k}")
                    done.append(f"{tag}: re-enumerated")
            elif typ == "local":
                given = {d["name"]: d["radius"] for d in r.get("decoders", []) if d.get("source") == "given"}
                again = _clean(props.local(an, given or None))
                again["lemma"] = props.LEMMA
                if again != r:
                    raise Refusal(f"{tag}: the locality certificate does not recompute")
                done.append(f"{tag}: every fault and input error placed again"
                            + (" (decoder radii are taken as given)" if given else ""))
            else:
                raise Refusal(f"{tag}: unknown result type")
        except (KeyError, TypeError, ValueError, IndexError) as e:
            raise Refusal(f"{tag}: malformed ({e})") from None
    return done


# ---------------------------------------------------------------------- composition
THEOREM = ("Let G be a gadget in which each listed part is a preparation fragment used as its "
           "specification (noiseless preparation, then one noisy wire per qubit), and let G' be G "
           "with each part replaced by an implementation: a preparation without inputs that touches "
           "the rest of G only through its output qubits, all of whose detectors are postselected. "
           "If each implementation refines its specification through order w, then every accepted "
           "fault set of j <= w faults in G' has the effect of an accepted fault set of at most j "
           "faults in G with the same input error.  Hence a property that holds in G for all input "
           "weights r and fault counts s with r + s <= T holds in G' for r + s <= min(T, w), provided "
           "its bound does not decrease when faults are added (prep, gate, ec and close all qualify), and "
           "a refinement of G against any specification through order T is a refinement of G' "
           "through order min(T, w).")


def compose(target, parts, cap=3_000_000):
    """target: receipt of G (built with spec parts).  parts: {role: refinement receipt of
    the implementation against library.prep(code, role, 'spec')}.  Returns a composition
    receipt; raises Refusal when a part is missing or does not match."""
    verify(target, cap=cap)
    g = Gadget.from_json(target["gadget"])
    need = [p for p in g.parts if p.get("method") == "spec"]
    if not need:
        raise Refusal("the target gadget has no part that was analysed as a specification")
    rows, w = [], None
    for p in need:
        role = p["role"]
        rec = parts.get(role)
        if rec is None:
            raise Refusal(f"no refinement receipt for part {role!r}")
        want = library.prep(library.code_of(_code_json(g, p)), role, "spec")
        if "composition" in rec:
            # the implementation is itself a composed gadget: a bundle
            # dict(composition=..., target=..., parts=...) whose conclusions include a
            # refinement against this part's specification
            verify_composition(rec["composition"], rec["target"], rec["parts"], cap)
            comp = rec["composition"]
            hit = [c for c in comp["conclusions"] if c["type"] == "refine"
                   and c.get("specification_sha256") == want.sha256()]
            if len(hit) != 1:
                raise Refusal(f"part {role!r}: the composition does not conclude a refinement of "
                              f"this part's specification")
            inner = Gadget.from_json(rec["target"]["gadget"])
            if inner.inputs or inner.decoders:
                raise Refusal(f"part {role!r}: the implementation must be a preparation without decoders")
            through = hit[0]["holds_through"]
            if through < 0:
                raise Refusal(f"part {role!r}: the refinement does not hold at any order")
            rows.append(dict(role=role, implementation=comp["target"] + " with " + ", ".join(
                                 f"{q['role']} := {q['implementation']}" for q in comp["parts"]),
                             implementation_sha256=comp["target_sha256"],
                             specification_sha256=want.sha256(), refines_through=through,
                             receipt_sha256=receipt_sha(comp)))
            w = through if w is None else min(w, through)
            continue
        verify(rec, cap=cap)
        ref = [r for r in rec["results"] if r.get("type") == "refine"]
        if len(ref) != 1:
            raise Refusal(f"part {role!r}: the receipt must hold exactly one refinement")
        impl = Gadget.from_json(rec["gadget"])
        if impl.inputs or impl.decoders:
            raise Refusal(f"part {role!r}: the implementation must be a preparation without decoders")
        spec = Gadget.from_json(ref[0]["spec"])
        if spec.sha256() != want.sha256():
            raise Refusal(f"part {role!r}: the refinement is against a different specification")
        if ref[0]["input_weight"] != 0 or ref[0]["holds_through"] < 0:
            raise Refusal(f"part {role!r}: the refinement does not hold at any order")
        rows.append(dict(role=role, implementation=impl.name, implementation_sha256=impl.sha256(),
                         specification_sha256=spec.sha256(), refines_through=ref[0]["holds_through"],
                         receipt_sha256=receipt_sha(rec)))
        w = ref[0]["holds_through"] if w is None else min(w, ref[0]["holds_through"])
    uniq = []
    for row in rows:
        if row not in uniq:
            uniq.append(row)
    rows = uniq
    concl = []
    for r in target["results"]:
        if r.get("type") in ("check", "close"):
            concl.append(dict(type=r["type"], kind=r.get("kind", "close"), weight=r.get("weight"),
                              target_holds_through=r["holds_through"],
                              holds_through=min(r["holds_through"], w)))
        elif r.get("type") == "refine":
            concl.append(dict(type="refine", kind="refines " + r["spec"]["name"], weight=None,
                              specification_sha256=r["spec"]["sha256"],
                              target_holds_through=r["holds_through"],
                              holds_through=min(r["holds_through"], w)))
    return _clean(dict(schema=COMPOSITION, tool="ftqc_gadgets", tool_version=VERSION,
                       target=g.name, target_sha256=g.sha256(), target_receipt_sha256=receipt_sha(target),
                       parts=rows, conclusions=concl, theorem=THEOREM))


def _code_json(g, part):
    for b in g.inputs + g.outputs:
        if b.code is not None and b.code.name == part["code"]:
            return b.code
    raise Refusal(f"part {part['role']!r}: its code {part['code']!r} is not a block code of the gadget")


def verify_composition(comp, target, parts, cap=3_000_000):
    if comp.get("schema") not in (COMPOSITION, "gadgetcert-composition-1"):
        raise Refusal(f"not a {COMPOSITION} receipt")
    again = compose(target, parts, cap)
    for k in ("target_sha256", "target_receipt_sha256", "parts", "conclusions", "theorem"):
        if again[k] != comp.get(k):
            raise Refusal(f"composition: {k} differs from what the receipts support")
    return ["target receipt", "part receipts", "conclusions recomputed"]


def summary(receipt):
    if receipt.get("schema") in (COMPOSITION, "gadgetcert-composition-1"):
        out = [f"{receipt['target']} with " + ", ".join(
            f"{p['role']} := {p['implementation']} (refines its specification through order "
            f"{p['refines_through']})" for p in receipt["parts"])]
        for c in receipt["conclusions"]:
            out.append(f"  {c['kind']}: holds through {c['holds_through']} "
                       f"(with specifications: {c['target_holds_through']})")
        return "\n".join(out)
    st = receipt["stats"]
    out = [f"{receipt['gadget']['name']}: {st['qubits']} qubits, {st['fault_locations']} fault "
           f"locations, {st['columns']} distinct fault effects, {st['detectors']} detectors, "
           f"{st['decoder_checks']} decoder checks"]
    for r in receipt["results"]:
        if r["type"] == "local":
            if not r["local"]:
                out.append(f"  not coordinate-local: {r['reason']}"
                           + (f" (fault {r['fault']})" if r.get("fault") else "")
                           + (f" (input {r['input_error']})" if r.get("input_error") else ""))
            else:
                thr = "any number of" if r["through"] is None else f"at most {r['through']}"
                props_ = [k for k in ("gate", "ec", "close") if r[k]]
                out.append(f"  coordinate-local on {r['coordinates']} coordinates: {', '.join(props_)} "
                           f"hold when faults and input errors touch {thr} coordinates"
                           + ("; also with X and Z counted separately" if r["sector"] else ""))
            continue
        if r["type"] == "check":
            line = (f"  {r['kind']} ({r['weight']} weight), r + s <= {r['t']}: holds through "
                    f"{r['holds_through']}")
            if r["violation"]:
                v = r["violation"]
                line += f"; fails at r = {v['r']}, s = {v['s']}: input {v['input_error']}, faults {v['faults']}"
        elif r["type"] == "close":
            line = f"  close to the code space for any input, s <= {r['s']}: holds through {r['holds_through']}"
            if r.get("method"):
                line += f" ({r['method']})"
        else:
            line = (f"  refines {r['spec']['name']} (input weight <= {r['input_weight']}), order <= "
                    f"{r['order']}: holds through {r['holds_through']}")
            if r["violation"]:
                v = r["violation"]
                line += f"; unmatched at order {v['order']}: faults {v['faults']}"
        if r.get("note"):
            line += f" [{r['note']}]"
        out.append(line)
    return "\n".join(out)
