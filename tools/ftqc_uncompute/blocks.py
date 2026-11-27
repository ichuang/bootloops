"""Arithmetic programs on code blocks: deletions, meaning, and binding to a schedule.

A block program is a list of logical operations on named wires (one wire is one logical
qubit, one code block):

    ["ZERO", b]             fresh wire holding 0
    ["CX", c, t]            controlled NOT
    ["XF", b]               logical X (a Pauli-frame update: no physical location)
    ["AND", u, v, t]        t := u AND v on a fresh wire, by consuming an AND resource
    ["ERASE", t, u, v, z]   measurement deletion of t, claimed valid because t = (u AND v) XOR z;
                            z is a wire, null, or ["const", 0 | 1]
    ["DROP", b]             retire a wire claimed to hold 0
    ["CUT", b, name]        hint: cut here (cuts.py); never needed for soundness
    ["PROBE", name]         name this point of the program for use in the relation
    ["LEMMA", ...]          ignored (annotation of another verifier)

Three things are certified.

deletions  The program is translated to a KickMix circuit (AND as CCX on a zero wire, ERASE as
           HMR with the fix-ups CZ u v, Z z, NEG; DROP as a reset) and verified with cuts:
           every reset wire is zero, every outcome-dependent phase cancels, outputs do not
           depend on outcomes, no phase is left.  So every measurement branch is the same
           permutation of basis states.
meaning    An integer relation between wire values at entry, at exit and at probes, proved
           for all inputs by `cuts.prove_relation`, and a list of wires proved to leave as
           they entered.  The relation is a statement; what it means (for instance that an
           output register holds a sum) is read off it by the reader, and `reading` in the
           specification is carried into the result unchecked for that purpose.
binding    Given the layers of a scheduled program (the format of ftqc_gadgets.cones), the
           layers without their error-correction steps are the lowering of the program:
           the same operations on every wire and every classical bit in the same order, with
           every controlled correction after the measurement it depends on.  Two words with
           equal projections on every wire and bit are the same circuit.  The lowering
           replaces AND and ERASE by gadgets whose identities are checked exactly here
           (`gadget_identities`), including the Pauli-frame terms: X on t when both
           measurement outcomes of an AND are 1, Z on z when the outcome of an ERASE is 1.
           Those frame terms and XF have no location in the layers; they belong to the
           classical control, which is assumed perfect.

100% Python standard library, exact integer and rational arithmetic only.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from typing import Any

from .cuts import CutResult, prove_equal, prove_relation, verify_with_cuts
from .gf2 import ANF, ANF_ONE, ANF_ZERO, anf_var, anf_xor
from .hmr import _SymbolicOverflow
from .ir_kickmix import KickMixCircuit, KickMixOp

SCHEMA = "ftqc_uncompute-blocks-1"
EC_OPS = ("ECa", "ECb", "EC2a", "EC2b")
DEFAULT_KINDS = {"zero": ["zero2", "zero", "std"], "and": ["and2", "and1"]}


class BlockError(ValueError):
    pass


def _ops(program: dict[str, Any]) -> list[list[Any]]:
    if not isinstance(program, dict) or not {"inputs", "outputs", "ops"} <= set(program):
        raise BlockError("a block program has inputs, outputs and ops")
    out = []
    for op in program["ops"]:
        if not isinstance(op, (list, tuple)) or not op or not isinstance(op[0], str):
            raise BlockError(f"not an operation: {op!r}")
        out.append(list(op))
    return out


def sha256(program: dict[str, Any]) -> str:
    return hashlib.sha256(json.dumps(program, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def to_kickmix(program: dict[str, Any]) -> dict[str, Any]:
    """The program as a KickMix circuit.  Returns dict(circuit, inputs={wire: qubit},
    outputs={wire: qubit}, zero=[qubits that enter as 0], probes={name: {wire: qubit}},
    counts)."""
    live: dict[str, int] = {}
    ops: list[KickMixOp] = []
    nq = nb = 0
    counts = {k: 0 for k in ("ZERO", "CX", "XF", "AND", "ERASE", "DROP")}
    probes: dict[str, dict[str, int]] = {}

    def fresh(w: str) -> int:
        nonlocal nq
        if not isinstance(w, str):
            raise BlockError(f"a wire name is a string, got {w!r}")
        if w in live:
            raise BlockError(f"wire {w} is already live")
        live[w] = nq
        nq += 1
        return live[w]

    def need(w: Any) -> int:
        if w not in live:
            raise BlockError(f"wire {w!r} is not live")
        return live[w]

    inputs = {w: fresh(w) for w in program["inputs"]}
    for i, op in enumerate(_ops(program)):
        k, a = op[0], op[1:]
        if k in counts:
            counts[k] += 1
        if k == "ZERO" and len(a) == 1:
            fresh(a[0])
        elif k == "CX" and len(a) == 2:
            if a[0] == a[1]:
                raise BlockError(f"op {i}: CX on one wire")
            ops.append(KickMixOp("CX", (need(a[0]), need(a[1]))))
        elif k == "XF" and len(a) == 1:
            ops.append(KickMixOp("X", (need(a[0]),)))
        elif k == "AND" and len(a) == 3:
            u, v = need(a[0]), need(a[1])
            if u == v:
                raise BlockError(f"op {i}: AND of a wire with itself")
            ops.append(KickMixOp("CCX", (u, v, fresh(a[2]))))
        elif k == "ERASE" and len(a) == 4:
            t, u, v, z = need(a[0]), need(a[1]), need(a[2]), a[3]
            if len({t, u, v}) < 3:
                raise BlockError(f"op {i}: ERASE on repeated wires")
            bit = nb
            nb += 1
            ops.append(KickMixOp("HMR", (t,), (bit,)))
            ops.append(KickMixOp("CZ", (u, v), condition=bit))
            if isinstance(z, str):
                if need(z) == t:
                    raise BlockError(f"op {i}: ERASE of a wire against itself")
                ops.append(KickMixOp("Z", (need(z),), condition=bit))
            elif isinstance(z, (list, tuple)) and len(z) == 2 and z[0] == "const" and z[1] in (0, 1):
                if z[1]:
                    ops.append(KickMixOp("NEG", condition=bit))
            elif z is not None:
                raise BlockError(f"op {i}: the last operand of ERASE is a wire, null or ['const', bit]")
            del live[a[0]]
        elif k == "DROP" and len(a) == 1:
            ops.append(KickMixOp("R", (need(a[0]),)))
            del live[a[0]]
        elif k == "CUT" and len(a) == 2:
            ops.append(KickMixOp("INFO_TAG", info_tag=f"@cut:q{need(a[0])}"))
        elif k == "PROBE" and len(a) == 1:
            if a[0] in probes or a[0] in ("in", "out") or not isinstance(a[0], str):
                raise BlockError(f"op {i}: probe name {a[0]!r} is taken")
            probes[a[0]] = dict(live)
            ops.append(KickMixOp("INFO_TAG", info_tag=f"@probe:{a[0]}"))
        elif k == "LEMMA":
            pass
        else:
            raise BlockError(f"op {i}: unknown operation {op!r}")
    if set(live) != set(program["outputs"]) or len(set(program["outputs"])) != len(program["outputs"]):
        raise BlockError("the live wires at the end are not the declared outputs")
    outputs = {w: live[w] for w in program["outputs"]}
    in_q = set(inputs.values())
    circuit = KickMixCircuit(tuple(ops), num_qubits=nq, num_bits=nb)
    return dict(circuit=circuit, inputs=inputs, outputs=outputs, probes=probes, counts=counts,
                zero=[q for q in range(nq) if q not in in_q])


# ---------------------------------------------------------------------------- meaning
def _value(result: CutResult, km: dict[str, Any], factor: dict[str, Any]) -> ANF:
    if not isinstance(factor, dict) or not {"at", "xor"} <= set(factor) or set(factor) - {"at", "xor", "not"}:
        raise BlockError(f"a factor is {{'at': point, 'xor': [wires]}}, got {factor!r}")
    at = factor["at"]
    acc: ANF = ANF_ONE if factor.get("not") else ANF_ZERO
    for w in factor["xor"]:
        if at == "in":
            if w not in km["inputs"]:
                raise BlockError(f"{w!r} is not an input")
            f = anf_var(result.input_var[km["inputs"][w]])
        elif at == "out":
            if w not in km["outputs"]:
                raise BlockError(f"{w!r} is not an output")
            f = result.outputs[km["outputs"][w]]
        else:
            if at not in km["probes"] or w not in km["probes"][at]:
                raise BlockError(f"wire {w!r} is not live at probe {at!r}")
            f = result.probes[at].get(km["probes"][at][w], ANF_ZERO)
        acc = anf_xor(acc, f)
    return acc


def _run_primitives(prims: list[tuple[Any, ...]], state: dict[tuple[tuple[str, int], ...], Fraction],
                    outcomes: dict[Any, int]) -> dict[tuple[tuple[str, int], ...], Fraction]:
    """Exact run of lowered operations on a superposition of basis states of named wires, for
    one choice of measurement outcomes.  Amplitudes are not normalised after a measurement."""
    half = Fraction(1, 2)
    for op in prims:
        k = op[0]
        nxt: dict[tuple[tuple[str, int], ...], Fraction] = {}

        def put(d: dict[Any, int], amp: Fraction) -> None:
            key = tuple(sorted(d.items(), key=repr))
            nxt[key] = nxt.get(key, Fraction(0)) + amp

        for key, amp in state.items():
            d = dict(key)
            if k == "NEW" and op[1] == "zero":
                d[op[2]] = 0
                put(d, amp)
            elif k == "NEW":                                      # the AND resource (1/2) sum |a, b, ab>
                for x in (0, 1):
                    for y in (0, 1):
                        e = dict(d)
                        e[op[2]], e[op[3]], e[op[4]] = x, y, x & y
                        put(e, amp * half)
            elif k == "CX":
                d[op[2]] ^= d[op[1]]
                put(d, amp)
            elif k == "MZ":
                if d.pop(op[1]) == outcomes[op[2]]:
                    put(d, amp)
            elif k == "MX":
                put(d, amp * (-1) ** (outcomes[op[2]] & d.pop(op[1])))
            elif k == "CXg":
                if outcomes[op[3]]:
                    d[op[2]] ^= d[op[1]]
                put(d, amp)
            elif k == "CZg":
                put(d, amp * (-1) ** (outcomes[op[3]] & d[op[1]] & d[op[2]]))
            else:
                raise BlockError(f"no rule for {op!r}")
        state = {key: amp for key, amp in nxt.items() if amp}
    return state


def gadget_identities() -> dict[str, Any]:
    """Exact check of the two measurement gadgets, as `lower` writes them, on every basis
    input and for every outcome.  A branch that maps every basis input to the same multiple
    of its target is the target isometry on all superpositions.  The Pauli-frame terms are
    applied here as the rule states them: X on t when both outcomes of an AND are 1; Z on z
    (or a sign, for a constant) when the outcome of an ERASE is 1."""
    report: dict[str, Any] = {}
    prims = lower(dict(inputs=["u", "v"], outputs=["u", "v", "t"], ops=[["AND", "u", "v", "t"]]))
    bits = [op[2] for op in prims if op[0] == "MZ"]
    for frame in (True, False):
        consts, ok = set(), True
        for u in (0, 1):
            for v in (0, 1):
                for m in ((0, 0), (0, 1), (1, 0), (1, 1)):
                    out = _run_primitives(prims, {(("u", u), ("v", v)): Fraction(1)}, dict(zip(bits, m)))
                    if frame and m == (1, 1):
                        out = {tuple(sorted({**dict(k), "t": dict(k)["t"] ^ 1}.items(), key=repr)): a for k, a in out.items()}
                    want = tuple(sorted({"u": u, "v": v, "t": u & v}.items(), key=repr))
                    if set(out) != {want}:
                        ok = False
                    else:
                        consts.add(out[want])
        report["and_with_frame" if frame else "and_without_frame"] = bool(ok and len(consts) == 1)
    for name, z in (("erase_wire", "z"), ("erase_none", None), ("erase_const1", ["const", 1])):
        wires = ["u", "v"] + (["z"] if z == "z" else [])
        prims = lower(dict(inputs=wires + ["t"], outputs=wires, ops=[["ERASE", "t", "u", "v", z]]))
        bit = next(op[2] for op in prims if op[0] == "MX")
        for frame in (True, False):
            signs: dict[int, set[Fraction]] = {0: set(), 1: set()}
            ok = True
            for vals in ((u, v, zz) for u in (0, 1) for v in (0, 1) for zz in ((0, 1) if z == "z" else (0,))):
                u, v, zz = vals
                zval = zz if z == "z" else (1 if z else 0)
                start = {"u": u, "v": v, "t": (u & v) ^ zval}
                if z == "z":
                    start["z"] = zz
                for r in (0, 1):
                    out = _run_primitives(prims, {tuple(sorted(start.items(), key=repr)): Fraction(1)}, {bit: r})
                    want = tuple(sorted({w: start[w] for w in wires}.items(), key=repr))
                    if set(out) != {want}:
                        ok = False
                        continue
                    amp = out[want]
                    if frame and z == "z":
                        amp *= (-1) ** (r & zz)
                    signs[r].add(amp)
            report[f"{name}_{'with' if frame else 'without'}_frame"] = bool(ok and all(len(x) == 1 for x in signs.values()))
    report["holds"] = bool(report["and_with_frame"] and not report["and_without_frame"]
                           and report["erase_wire_with_frame"] and not report["erase_wire_without_frame"]
                           and report["erase_none_with_frame"] and report["erase_const1_with_frame"])
    return report


# ---------------------------------------------------------------------------- binding
def lower(program: dict[str, Any]) -> list[tuple[Any, ...]]:
    """Primitive operations of the program, with canonical names for the ancillas and bits of
    its gadgets: ("~a", t, k) and ("~b", t, k) for the k-th AND onto wire t, ("m", ancilla) for
    their measured bits, ("r", t, k) for the k-th ERASE of wire t."""
    out: list[tuple[Any, ...]] = []
    n_and: dict[str, int] = {}
    n_erase: dict[str, int] = {}
    for op in _ops(program):
        k, a = op[0], op[1:]
        if k == "ZERO":
            out.append(("NEW", "zero", a[0]))
        elif k == "CX":
            out.append(("CX", a[0], a[1]))
        elif k == "DROP":
            out.append(("DROP", a[0]))
        elif k == "AND":
            u, v, t = a
            i = n_and.get(t, 0)
            n_and[t] = i + 1
            x, y = ("~a", t, i), ("~b", t, i)
            out += [("NEW", "and", x, y, t), ("CX", u, x), ("CX", v, y), ("MZ", x, ("m", x)), ("MZ", y, ("m", y)),
                    ("CXg", v, t, ("m", x)), ("CXg", u, t, ("m", y))]
        elif k == "ERASE":
            t, u, v, _ = a
            i = n_erase.get(t, 0)
            n_erase[t] = i + 1
            out += [("MX", t, ("r", t, i)), ("CZg", u, v, ("r", t, i))]
    return out


def _projections(seq: list[tuple[int, tuple[Any, ...]]]) -> dict[Any, list[tuple[Any, ...]]]:
    proj: dict[Any, list[tuple[Any, ...]]] = {}
    for _, op in seq:
        k = op[0]
        wires = list(op[2:]) if k == "NEW" else list(op[1:])
        for w in wires:
            proj.setdefault(w, []).append(op)
    return proj


def bind(program: dict[str, Any], layers: list[list[list[Any]]], kinds: dict[str, list[str]] | None = None) -> dict[str, Any]:
    """Whether the layers, without error correction, are the lowering of the program."""
    kinds = kinds or DEFAULT_KINDS
    want = [(0, op) for op in lower(program)]
    got: list[tuple[int, tuple[Any, ...]]] = []
    rename: dict[Any, Any] = {}
    bit_of: dict[Any, Any] = {}
    made_at: dict[Any, int] = {}
    n_and: dict[Any, int] = {}
    n_erase: dict[Any, int] = {}
    skipped = 0
    problems: list[str] = []

    def nm(b: Any) -> Any:
        b = tuple(b) if isinstance(b, list) else b
        return rename.get(b, b)

    for li, layer in enumerate(layers):
        touched: dict[Any, int] = {}
        ordered = sorted(layer, key=lambda o: o[0] != "NEW")
        for op in ordered:
            k = op[0]
            if k in EC_OPS:
                skipped += 1
                continue
            if k == "NEW":
                kind, blocks = op[1], [tuple(b) if isinstance(b, list) else b for b in op[2:]]
                if len(blocks) == 1 and kind in kinds["zero"]:
                    got.append((li, ("NEW", "zero", blocks[0])))
                elif len(blocks) == 3 and kind in kinds["and"]:
                    t = blocks[2]
                    i = n_and.get(t, 0)
                    n_and[t] = i + 1
                    rename[blocks[0]], rename[blocks[1]] = ("~a", t, i), ("~b", t, i)
                    got.append((li, ("NEW", "and", ("~a", t, i), ("~b", t, i), t)))
                else:
                    problems.append(f"layer {li}: resource {op!r} is not a zero wire or an AND resource")
                continue
            if k in ("CX", "DROP"):
                new = (k,) + tuple(nm(b) for b in op[1:])
            elif k in ("MZ", "MX"):
                b = nm(op[1])
                if len(op) != 3:
                    problems.append(f"layer {li}: {op!r} does not record its outcome")
                    continue
                if k == "MZ":
                    canon = ("m", b)
                else:
                    i = n_erase.get(b, 0)
                    n_erase[b] = i + 1
                    canon = ("r", b, i)
                if op[2] in bit_of:
                    problems.append(f"layer {li}: bit {op[2]!r} is written twice")
                bit_of[op[2]] = canon
                made_at[canon] = li
                new = (k, b, canon)
            elif k in ("CXg", "CZg"):
                if len(op) != 4 or op[3] not in bit_of:
                    problems.append(f"layer {li}: {op!r} is controlled by a bit that was not measured before")
                    continue
                canon = bit_of[op[3]]
                if made_at[canon] >= li:
                    problems.append(f"layer {li}: {op!r} is in the layer of the measurement it depends on")
                new = (k, nm(op[1]), nm(op[2]), canon)
            else:
                problems.append(f"layer {li}: operation {op!r} has no counterpart in a block program")
                continue
            for w in new[1:3] if k in ("CX", "CXg", "CZg") else new[1:2]:
                touched[w] = touched.get(w, 0) + 1
                if touched[w] > 1:
                    problems.append(f"layer {li}: wire {w!r} is used twice")
            got.append((li, new))
    pw, pg = _projections(want), _projections(got)
    for w in sorted(set(pw) | set(pg), key=repr):
        if pw.get(w) != pg.get(w):
            a, b = pw.get(w, []), pg.get(w, [])
            i = next((j for j, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
            problems.append(f"wire {w!r}: operation {i} is {b[i] if i < len(b) else None!r} in the layers and "
                            f"{a[i] if i < len(a) else None!r} in the program")
            if len(problems) > 20:
                break
    return dict(bound=not problems, primitive_operations=len(want), layers=len(layers),
                error_correction_steps_skipped=skipped, wires_and_bits_compared=len(set(pw) | set(pg)),
                problems=problems[:20])


# ---------------------------------------------------------------------------- the certificate
def certify(spec: dict[str, Any]) -> dict[str, Any]:
    """Result for {"kind": "block-program", "program": ..., "spec": {"relation": [[coef,
    [factor, ...]], ...], "preserved": [wires], "reading": text}, "layers": [...] (optional),
    "kinds": {"zero": [...], "and": [...]} (optional), "max_monomials": n (optional)}.
    A factor is {"at": "in" | "out" | probe name, "xor": [wires], "not": bool (optional)}."""
    if not isinstance(spec, dict) or spec.get("kind") != "block-program":
        raise BlockError("kind must be block-program")
    extra = set(spec) - {"kind", "program", "spec", "layers", "kinds", "max_monomials"}
    if extra or "program" not in spec or "spec" not in spec:
        raise BlockError(f"expected program and spec; unknown fields {sorted(extra)}")
    program, sp = spec["program"], spec["spec"]
    if not isinstance(sp, dict) or set(sp) - {"relation", "preserved", "reading"} or "relation" not in sp:
        raise BlockError("spec has relation, and optionally preserved and reading")
    budget = spec.get("max_monomials", 8192)
    if type(budget) is not int or budget < 16:
        raise BlockError("max_monomials is an integer of at least 16")
    km = to_kickmix(program)
    res = verify_with_cuts(km["circuit"], inputs=sorted(km["inputs"].values()), outputs=sorted(km["outputs"].values()),
                           max_monomials=budget)
    out: dict[str, Any] = dict(
        program_sha256=sha256(program), wires=dict(inputs=len(km["inputs"]), outputs=len(km["outputs"]),
                                                   total=km["circuit"].num_qubits),
        operations=km["counts"], deletions=res.to_dict(), reading=sp.get("reading"),
        limits=["classical control, including Pauli-frame updates (XF, and the frame terms of the AND and ERASE "
                "gadgets), is assumed perfect",
                "the AND resource is assumed ideal here; its quality is a separate contract",
                "the reading of the relation is not checked"])
    relation = preserved = None
    if res.valid:
        terms = []
        for item in sp["relation"]:
            if not isinstance(item, (list, tuple)) or len(item) != 2 or type(item[0]) is not int:
                raise BlockError(f"a relation term is [integer, [factors]], got {item!r}")
            terms.append((item[0], tuple(_value(res, km, f) for f in item[1])))
        relation = prove_relation(res, terms, max_monomials=budget).to_dict()
        preserved = {}
        for w in sp.get("preserved", []):
            if w not in km["inputs"] or w not in km["outputs"]:
                raise BlockError(f"{w!r} is not both an input and an output")
            try:
                preserved[w] = prove_equal(res, res.outputs[km["outputs"][w]],
                                           anf_var(res.input_var[km["inputs"][w]]), budget)
            except _SymbolicOverflow:
                preserved[w] = False
    out["relation"] = relation
    out["preserved"] = None if preserved is None else dict(wires=len(preserved), all=all(preserved.values()),
                                                           failed=sorted(w for w, ok in preserved.items() if not ok))
    out["gadgets"] = gadget_identities()
    binding = None
    if "layers" in spec:
        binding = bind(program, spec["layers"], spec.get("kinds"))
    out["binding"] = binding
    ok = (res.valid and res.phase_free and relation is not None and relation["proved"] and out["preserved"]["all"]
          and out["gadgets"]["holds"] and (binding is None or binding["bound"]))
    out["status"] = "proved" if ok else "refused"
    return out


def receipt(spec: dict[str, Any]) -> dict[str, Any]:
    digest = hashlib.sha256(json.dumps(spec, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {"schema": SCHEMA, "input_sha256": digest, "result": certify(spec)}


def verify_receipt(spec: dict[str, Any], stored: dict[str, Any]) -> None:
    again = receipt(spec)
    if json.dumps(again, sort_keys=True) != json.dumps(stored, sort_keys=True):
        raise BlockError("receipt differs from recomputed evidence")


def simulate(program: dict[str, Any], inputs: dict[str, int]) -> dict[str, int]:
    """The program on one basis input (a classical run: phases are not seen).  Raises
    BlockError when an ERASE or DROP claim fails on this input."""
    v = dict(inputs)
    for i, op in enumerate(_ops(program)):
        k, a = op[0], op[1:]
        if k == "ZERO":
            v[a[0]] = 0
        elif k == "CX":
            v[a[1]] ^= v[a[0]]
        elif k == "XF":
            v[a[0]] ^= 1
        elif k == "AND":
            v[a[2]] = v[a[0]] & v[a[1]]
        elif k == "ERASE":
            z = a[3]
            zv = 0 if z is None else (z[1] if isinstance(z, (list, tuple)) else v[z])
            if v[a[0]] != (v[a[1]] & v[a[2]]) ^ zv:
                raise BlockError(f"op {i}: the erased wire is not (u AND v) XOR z on this input")
            del v[a[0]]
        elif k == "DROP":
            if v[a[0]]:
                raise BlockError(f"op {i}: the dropped wire is 1 on this input")
            del v[a[0]]
    return v


# ---------------------------------------------------------------------------- arithmetic programs
def _carry_chain(ops: list[list[Any]], A: list[str] | None, B: list[str], C: list[str], n: int, const: int = 0) -> None:
    """Carries of B + A (or B + const) into C[1..n]; afterwards B[i] holds b_i xor a_i and
    C[i], for i < n, the carry into bit i xor a_i."""
    for i in range(n):
        bit = (const >> i) & 1
        if A is None:
            if bit:
                ops += [["XF", B[i]], ["XF", C[i]]]
        else:
            ops += [["CX", A[i], B[i]], ["CX", A[i], C[i]]]
        ops.append(["AND", C[i], B[i], C[i + 1]])
        if A is None:
            if bit:
                ops.append(["XF", C[i + 1]])
        else:
            ops.append(["CX", A[i], C[i + 1]])


def adder(n: int, family: str = "qq", const: int | None = None) -> dict[str, Any]:
    """Ripple-carry adder with measured carry deletion.  qq: B += A; constant: B += const;
    controlled: B += H * A.  The overflow leaves in C{n} (controlled: OUT).  Wires A{i}, B{i},
    C{i}, H."""
    if n < 1 or family not in ("qq", "constant", "controlled"):
        raise BlockError("adder: n >= 1 and family qq, constant or controlled")
    A = None if family == "constant" else [f"A{i}" for i in range(n)]
    B = [f"B{i}" for i in range(n)]
    C = [f"C{i}" for i in range(n + 1)]
    if family == "constant" and const is None:
        const = (1 << n) - 1
    const = const or 0
    ops: list[list[Any]] = [["ZERO", C[0]]]
    _carry_chain(ops, A, B, C, n, const)
    if family == "controlled":
        ops.append(["AND", "H", C[n], "OUT"])
        for i in range(n - 1, -1, -1):
            ops += [["ERASE", C[i + 1], C[i], B[i], A[i]], ["CX", A[i], B[i]], ["AND", "H", C[i], f"T{i}"],
                    ["CX", f"T{i}", B[i]], ["ERASE", f"T{i}", "H", C[i], None], ["CX", A[i], C[i]]]
    else:
        for i in range(n - 1, -1, -1):
            bit = (const >> i) & 1
            if i + 1 < n:
                ops.append(["ERASE", C[i + 1], C[i], B[i], A[i] if A else ["const", bit]])
            if A is None:
                if bit:
                    ops.append(["XF", C[i]])
            else:
                ops.append(["CX", A[i], C[i]])
            ops.append(["CX", C[i], B[i]])
    ops.append(["DROP", C[0]])
    ins = B + (A or []) + (["H"] if family == "controlled" else [])
    outs = ins + (["OUT"] if family == "controlled" else [C[n]])
    return dict(name=f"{family}_adder_{n}", inputs=ins, outputs=outs, ops=ops)


def subtractor(n: int) -> dict[str, Any]:
    """B := (B - A) mod 2^n, with the borrow [A > B] in C{n}: complement, add, complement."""
    p = adder(n, "qq")
    B = [f"B{i}" for i in range(n)]
    ops = [["XF", b] for b in B] + p["ops"] + [["XF", b] for b in B]
    return dict(name=f"subtractor_{n}", inputs=p["inputs"], outputs=p["outputs"], ops=ops)


def comparator(n: int) -> dict[str, Any]:
    """F := [A > B]: the top carry of A + not(B), copied out, with the chain deleted."""
    A = [f"A{i}" for i in range(n)]
    B = [f"B{i}" for i in range(n)]
    C = [f"C{i}" for i in range(n + 1)]
    ops: list[list[Any]] = [["XF", b] for b in B] + [["ZERO", C[0]]]
    _carry_chain(ops, A, B, C, n)
    ops += [["ZERO", "F"], ["CX", C[n], "F"]]
    for i in range(n - 1, -1, -1):
        ops += [["ERASE", C[i + 1], C[i], B[i], A[i]], ["CX", A[i], C[i]], ["CX", A[i], B[i]]]
    ops += [["DROP", C[0]]] + [["XF", b] for b in B]
    return dict(name=f"comparator_{n}", inputs=B + A, outputs=B + A + ["F"], ops=ops)


def multiplier(n: int) -> dict[str, Any]:
    """P := A * B on 2n fresh wires: one controlled addition of B per bit of A."""
    A = [f"A{i}" for i in range(n)]
    B = [f"B{i}" for i in range(n)]
    P = [f"P{k}" for k in range(2 * n)]
    ops: list[list[Any]] = [["ZERO", P[k]] for k in range(n)]
    for i in range(n):
        C = [f"C{i}_{j}" for j in range(n + 1)]
        W = P[i:i + n]
        ops.append(["ZERO", C[0]])
        _carry_chain(ops, B, W, C, n)
        ops.append(["AND", A[i], C[n], P[i + n]])
        for j in range(n - 1, -1, -1):
            T = f"T{i}_{j}"
            ops += [["ERASE", C[j + 1], C[j], W[j], B[j]], ["CX", B[j], W[j]], ["AND", A[i], C[j], T],
                    ["CX", T, W[j]], ["ERASE", T, A[i], C[j], None], ["CX", B[j], C[j]]]
        ops.append(["DROP", C[0]])
    return dict(name=f"multiplier_{n}", inputs=A + B, outputs=A + B + P, ops=ops)


def arithmetic_spec(program: dict[str, Any], kind: str, n: int, family: str = "qq", const: int | None = None) -> dict[str, Any]:
    """The block-program specification of an adder (kind add), subtractor (sub), comparator
    (cmp) or multiplier (mul) that uses the wire names of this module: A{i}, B{i}, C{i}, H,
    OUT, F, P{k}.  For a comparator a probe is placed before the wire F is made; the sum bits
    of A + not(B) are B{i} xor C{i} xor A{i} there."""
    A = [f"A{i}" for i in range(n)]
    B = [f"B{i}" for i in range(n)]
    w = lambda at, *ws: {"at": at, "xor": list(ws)}
    ops = [list(o) for o in program["ops"]]
    rel: list[list[Any]] = []
    preserved: list[str] = []
    if kind == "add":
        top = "OUT" if family == "controlled" else f"C{n}"
        rel += [[1 << i, [w("in", B[i])]] for i in range(n)] + [[-(1 << i), [w("out", B[i])]] for i in range(n)]
        rel.append([-(1 << n), [w("out", top)]])
        if family == "qq":
            rel += [[1 << i, [w("in", A[i])]] for i in range(n)]
            preserved, reading = A, f"B_out + 2^{n} {top} = B + A"
        elif family == "constant":
            c = (1 << n) - 1 if const is None else const
            rel.append([c, []])
            reading = f"B_out + 2^{n} {top} = B + {c}"
        elif family == "controlled":
            rel += [[1 << i, [w("in", "H"), w("in", A[i])]] for i in range(n)]
            preserved, reading = A + ["H"], f"B_out + 2^{n} {top} = B + H * A"
        else:
            raise BlockError("family is qq, constant or controlled")
    elif kind == "sub":
        rel += [[1 << i, [w("in", B[i])]] for i in range(n)] + [[-(1 << i), [w("in", A[i])]] for i in range(n)]
        rel += [[-(1 << i), [w("out", B[i])]] for i in range(n)] + [[1 << n, [w("out", f"C{n}")]]]
        preserved = A
        reading = f"B_out - 2^{n} C{n} = B - A with B_out below 2^{n}: B_out = (B - A) mod 2^{n} and C{n} = [A > B]"
    elif kind == "cmp":
        at = next((i for i, o in enumerate(ops) if o[:2] == ["ZERO", "F"]), None)
        if at is None:
            raise BlockError("a comparator makes the wire F")
        ops.insert(at, ["PROBE", "mid"])
        rel += [[1 << i, [w("in", A[i])]] for i in range(n)] + [[-(1 << i), [w("in", B[i])]] for i in range(n)]
        rel.append([(1 << n) - 1, []])
        rel += [[-(1 << i), [w("mid", B[i], f"C{i}", A[i])]] for i in range(n)]
        rel.append([-(1 << n), [w("out", "F")]])
        preserved = A + B
        reading = (f"A + (2^{n} - 1 - B) = S + 2^{n} F with S a sum of {n} bits of weights 2^i, "
                   f"so S is below 2^{n} and F = [A > B]")
    elif kind == "mul":
        rel += [[1 << (i + j), [w("in", A[i]), w("in", B[j])]] for i in range(n) for j in range(n)]
        rel += [[-(1 << k), [w("out", f"P{k}")]] for k in range(2 * n)]
        preserved, reading = A + B, "P = A * B"
    else:
        raise BlockError("kind is add, sub, cmp or mul")
    prog = dict(name=program.get("name", kind), inputs=list(program["inputs"]), outputs=list(program["outputs"]), ops=ops)
    return dict(kind="block-program", program=prog, spec=dict(relation=rel, preserved=preserved, reading=reading))
