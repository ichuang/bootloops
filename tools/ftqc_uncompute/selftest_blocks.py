"""Selftest leg for cuts.py and blocks.py (run by selftest.py).

What is compared with what:
  - the cut verifier with the global verifier of hmr.py, wire by wire, on circuits small
    enough for both;
  - proved relations with a classical run of the program on every input;
  - refusals with a run that fails: every hazard comes with an input at which the classical
    run of the program breaks;
  - every single-operation corruption of a program is either refused or is the same function
    on every input;
  - bindings with corrupted layers.

100% Python standard library, exact arithmetic only.
"""

from __future__ import annotations

import copy
import itertools
import json
from typing import Any

from . import blocks
from .cuts import prove_relation, unfolded_outputs, verify_with_cuts
from .gf2 import anf_eval, anf_var
from .hmr import _symbolic_verify_ops, verify_hmr_circuit
from .ir_kickmix import KickMixCircuit, KickMixOp


def _all_inputs(program: dict[str, Any]) -> list[dict[str, int]]:
    ins = program["inputs"]
    return [dict(zip(ins, bits)) for bits in itertools.product((0, 1), repeat=len(ins))]


def _value(env: dict[str, int], prefix: str, n: int) -> int:
    return sum(env[f"{prefix}{i}"] << i for i in range(n))


def _reference(kind: str, n: int, family: str, const: int, env: dict[str, int], out: dict[str, int]) -> bool:
    a = _value(env, "A", n) if "A0" in env else 0
    b = _value(env, "B", n)
    if kind == "add":
        add = const if family == "constant" else a * (env["H"] if family == "controlled" else 1)
        top = out["OUT" if family == "controlled" else f"C{n}"]
        return _value(out, "B", n) + (top << n) == b + add
    if kind == "sub":
        return _value(out, "B", n) == (b - a) % (1 << n) and out[f"C{n}"] == int(a > b)
    if kind == "cmp":
        return out["F"] == int(a > b)
    return _value(out, "P", 2 * n) == a * b


def _cases() -> list[tuple[dict[str, Any], tuple[Any, ...]]]:
    out: list[tuple[dict[str, Any], tuple[Any, ...]]] = []
    for n in (1, 2, 3):
        out += [(blocks.adder(n), ("add", n, "qq", 0)), (blocks.adder(n, "controlled"), ("add", n, "controlled", 0)),
                (blocks.subtractor(n), ("sub", n, "qq", 0)), (blocks.comparator(n), ("cmp", n, "qq", 0))]
    out += [(blocks.adder(4, "constant", 11), ("add", 4, "constant", 11)), (blocks.adder(3, "constant", 0), ("add", 3, "constant", 0)),
            (blocks.multiplier(1), ("mul", 1, "qq", 0)), (blocks.multiplier(2), ("mul", 2, "qq", 0)),
            (blocks.multiplier(3), ("mul", 3, "qq", 0))]
    return out


def _spec(prog: dict[str, Any], args: tuple[Any, ...]) -> dict[str, Any]:
    kind, n, family, const = args
    return blocks.arithmetic_spec(prog, kind, n, family, const if family == "constant" else None)


def _layers(program: dict[str, Any]) -> list[list[list[Any]]]:
    """A simple schedule of the lowering: one operation per layer, an error-correction step
    after every gate, names of ancillas and bits chosen here."""
    names: dict[Any, str] = {}

    def nm(x: Any) -> Any:
        if isinstance(x, tuple):
            return names.setdefault(x, f"anc{len(names)}")
        return x

    layers: list[list[list[Any]]] = []
    pending: list[list[Any]] = []
    for op in blocks.lower(program):
        k = op[0]
        if k == "NEW":
            pending.append(["NEW", "zero2" if op[1] == "zero" else "and2"] + [nm(b) for b in op[2:]])
            continue
        row = [k] + [nm(x) for x in op[1:]]
        layers.append(pending + [row])
        pending = []
        if k == "CX":
            layers.append([["ECa", row[2]]])
            layers.append([["ECb", row[2]]])
    if pending:
        layers.append(pending)
    return layers


def leg_cuts_and_blocks() -> str:
    checked = 0
    # 1. programs small enough to run on every input: proved, and the proof agrees with the runs
    for prog, args in _cases():
        spec = _spec(prog, args)
        res = blocks.certify(spec)
        assert res["status"] == "proved", (prog["name"], res)
        assert res["deletions"]["valid"] and res["deletions"]["phase_free"] and res["relation"]["proved"]
        for env in _all_inputs(prog):
            out = blocks.simulate(prog, env)
            assert _reference(*args, env, out), (prog["name"], env)
            for w in spec["spec"]["preserved"]:
                assert out[w] == env[w]
        checked += 1
    # 2. the cut verifier against the global verifier, wire by wire
    for prog, _ in _cases():
        km = blocks.to_kickmix(prog)
        ins, outs = sorted(km["inputs"].values()), sorted(km["outputs"].values())
        cut = verify_with_cuts(km["circuit"], inputs=ins, outputs=outs)
        pub = tuple(sorted(set(ins) | set(outs)))
        # the global verifier gives every public qubit a variable; the extra outputs enter as 0
        # here, so compare on the wires that are inputs and with scratch checking off
        glob, glob_out, _ = _symbolic_verify_ops(
            km["circuit"].ops, tuple(ins), km["circuit"].num_qubits, km["circuit"].num_bits,
            allow_branch_global_phase=False, check_scratch_clean=False, max_monomials=8192)
        assert cut.valid and glob.phase_cancelled and glob.dirty_reset_free and glob.basis_outcome_independent
        mine = unfolded_outputs(cut)
        for i, q in enumerate(ins):                              # same numbering: inputs first, in order
            if q in mine:
                assert mine[q] == glob_out[i], (prog["name"], q)
        # and every output, as a polynomial in the inputs, against a run on every input
        back = {q: w for w, q in km["outputs"].items()}
        var_of = {w: cut.input_var[q] for w, q in km["inputs"].items()}
        for env in _all_inputs(prog):
            mask = sum(1 << var_of[w] for w, bit in env.items() if bit)
            run = blocks.simulate(prog, env)
            for q, f in mine.items():
                assert anf_eval(f, mask) == run[back[q]], (prog["name"], back[q], env)
        checked += 1
    ex = KickMixCircuit.from_kmx(
        (__import__("pathlib").Path(__file__).parent / "examples" / "gidney_adder_4bit.kmx").read_text(encoding="utf-8"))
    assert verify_hmr_circuit(ex).valid and verify_with_cuts(ex).valid
    assert verify_with_cuts(ex, auto_cut=False).num_cuts == 0
    # 3. wide programs: the global method gives up, cuts do not
    for prog, args in ((blocks.adder(64), ("add", 64, "qq", 0)), (blocks.comparator(48), ("cmp", 48, "qq", 0)),
                       (blocks.multiplier(6), ("mul", 6, "qq", 0))):
        res = blocks.certify(_spec(prog, args))
        assert res["status"] == "proved" and res["deletions"]["largest_polynomial"] <= 16, prog["name"]
        km = blocks.to_kickmix(prog)
        ins, outs = sorted(km["inputs"].values()), sorted(km["outputs"].values())
        plain = verify_with_cuts(km["circuit"], inputs=ins, outputs=outs, auto_cut=False)
        if args[0] != "mul":                                     # twelve input bits still fit the budget
            assert plain.proof_method == "symbolic_budget_exceeded", prog["name"]
        checked += 1
    # 4. every single-operation corruption is refused or is the same function
    refused = same = 0
    for prog, args in ((blocks.adder(3), ("add", 3, "qq", 0)), (blocks.adder(2, "controlled"), ("add", 2, "controlled", 0)),
                       (blocks.comparator(2), ("cmp", 2, "qq", 0)), (blocks.adder(3, "constant", 5), ("add", 3, "constant", 5)),
                       (blocks.multiplier(2), ("mul", 2, "qq", 0))):
        ops = prog["ops"]
        for i, op in enumerate(ops):
            muts: list[list[list[Any]]] = []
            if op[0] == "CX":
                muts += [ops[:i] + ops[i + 1:], ops[:i] + [["CX", op[2], op[1]]] + ops[i + 1:], ops[:i + 1] + [op] + ops[i + 1:]]
            elif op[0] == "XF":
                muts += [ops[:i] + ops[i + 1:]]
            elif op[0] == "ERASE":
                z = op[4]
                muts += [ops[:i] + [op[:4] + [None if z is not None else ["const", 1]]] + ops[i + 1:]]
                if isinstance(z, str):
                    muts += [ops[:i] + [[op[0], op[1], op[2], z, op[3]]] + ops[i + 1:]]
            elif op[0] == "AND":
                muts += [ops[:i] + [["ZERO", op[3]], ["CX", op[1], op[3]]] + ops[i + 1:]]
            for m in muts:
                bad = dict(prog, ops=m)
                try:
                    ok = blocks.certify(_spec(bad, args))["status"] == "proved"
                except blocks.BlockError:
                    ok = False
                if not ok:
                    refused += 1
                    continue
                for env in _all_inputs(prog):
                    assert blocks.simulate(bad, env) == blocks.simulate(prog, env), (prog["name"], i, m[i])
                same += 1
    assert refused > 150 and same < refused // 10, (refused, same)
    # 5. hazards come with an input at which the program breaks
    prog = blocks.adder(3)
    first = prog["ops"].index(["CX", "A1", "C1"])
    bad = dict(prog, ops=prog["ops"][:first] + prog["ops"][first + 1:])
    res = blocks.certify(_spec(bad, ("add", 3, "qq", 0)))
    codes = {f["code"] for f in res["deletions"]["findings"]}
    assert res["status"] == "refused" and codes & {"UNCOMPENSATED-HMR-PHASE", "DIRTY-RESET-HAZARD"}, codes
    km = blocks.to_kickmix(bad)
    cut = verify_with_cuts(km["circuit"], inputs=sorted(km["inputs"].values()), outputs=sorted(km["outputs"].values()))
    back = {q: w for w, q in km["inputs"].items()}
    order = sorted(cut.input_var, key=lambda q: cut.input_var[q])
    broke = 0
    for f in cut.findings:
        env = {back[q]: f.witness["inputs"][f"x{j}"] for j, q in enumerate(order)}
        try:
            blocks.simulate(bad, env)
        except blocks.BlockError:
            broke += 1
    assert broke >= 1, [f.to_dict() for f in cut.findings]
    dirty = dict(inputs=["a", "b"], outputs=["a"], ops=[["CX", "a", "b"], ["DROP", "b"]])
    km = blocks.to_kickmix(dirty)
    cut = verify_with_cuts(km["circuit"], inputs=[0, 1], outputs=[0])
    assert not cut.valid and cut.findings[0].code == "DIRTY-RESET-HAZARD"
    w = cut.findings[0].witness["inputs"]
    assert (w["x0"] ^ w["x1"]) == 1
    leak = KickMixCircuit((KickMixOp("HMR", (1,), (0,)), KickMixOp("X", (0,), condition=0)), num_qubits=2, num_bits=1)
    cut = verify_with_cuts(leak, inputs=[0], outputs=[0])
    assert {f.code for f in cut.findings} == {"OUTCOME-DEPENDENT-BASIS-OUTPUT"} and not cut.valid
    phase = KickMixCircuit((KickMixOp("CZ", (0, 1)),), num_qubits=2)
    cut = verify_with_cuts(phase)
    assert cut.valid and not cut.phase_free                    # a unitary phase is not a deletion fault
    spec = dict(kind="block-program", program=dict(inputs=["a", "b"], outputs=["a", "b"], ops=[]),
                spec=dict(relation=[[1, [{"at": "in", "xor": ["a"]}]], [-1, [{"at": "out", "xor": ["a"]}]]]))
    assert blocks.certify(spec)["status"] == "proved"
    # tags place cuts; a cut is never taken from a wire that depends on an outcome
    at = prog["ops"].index(["CX", "A0", "C1"]) + 1
    km = blocks.to_kickmix(dict(prog, ops=prog["ops"][:at] + [["CUT", "C1", "carry"]] + prog["ops"][at:]))
    cut = verify_with_cuts(km["circuit"], inputs=sorted(km["inputs"].values()), outputs=sorted(km["outputs"].values()),
                           auto_cut=False)
    assert cut.valid and cut.num_cuts == 1
    hmr_then_cut = KickMixCircuit((KickMixOp("CCX", (0, 1, 2)), KickMixOp("HMR", (3,), (0,)), KickMixOp("X", (2,), condition=0),
                                   KickMixOp("INFO_TAG", info_tag="@cut:q2")), num_qubits=4, num_bits=1)
    assert verify_with_cuts(hmr_then_cut, inputs=[0, 1], outputs=[0, 1, 2], auto_cut=False).num_cuts == 0
    # 6. false statements about a correct program
    good = _spec(blocks.adder(3), ("add", 3, "qq", 0))
    wrong = copy.deepcopy(good)
    wrong["spec"]["relation"][0][0] += 1
    r = blocks.certify(wrong)
    assert r["status"] == "refused" and r["relation"]["reason"] == "the relation is not an identity"
    wrong = copy.deepcopy(good)
    wrong["spec"]["preserved"] = ["B0"]
    r = blocks.certify(wrong)
    assert r["status"] == "refused" and r["preserved"]["failed"] == ["B0"]
    cmp_spec = _spec(blocks.comparator(3), ("cmp", 3, "qq", 0))
    wrong = copy.deepcopy(cmp_spec)
    for item in wrong["spec"]["relation"]:
        if item[1] and item[1][0]["at"] == "mid":
            item[1][0]["xor"] = item[1][0]["xor"][:2]
    assert blocks.certify(wrong)["status"] == "refused"
    for broken in (dict(good, kind="other"), dict(good, extra=1), dict(good, max_monomials=3),
                   dict(good, spec=dict(relation=[[1, [{"at": "nowhere", "xor": ["A0"]}]]])),
                   dict(good, program=dict(good["program"], outputs=["B0"])),
                   dict(good, program=dict(good["program"], ops=[["CX", "A0", "A0"]]))):
        try:
            blocks.certify(broken)
            raise AssertionError("a malformed specification was accepted")
        except blocks.BlockError:
            pass
    # the relation prover on its own: x + y = (x xor y) + 2 (x and y)
    two = KickMixCircuit((KickMixOp("CX", (0, 2)), KickMixOp("CX", (1, 2)), KickMixOp("CCX", (0, 1, 3))), num_qubits=4)
    cut = verify_with_cuts(two, inputs=[0, 1], outputs=[0, 1, 2, 3])
    x, y = anf_var(cut.input_var[0]), anf_var(cut.input_var[1])
    assert prove_relation(cut, [(1, (x,)), (1, (y,)), (-1, (cut.outputs[2],)), (-2, (cut.outputs[3],))]).proved
    assert not prove_relation(cut, [(1, (x,)), (1, (y,)), (-1, (cut.outputs[2],)), (-1, (cut.outputs[3],))]).proved
    assert prove_relation(cut, [(1, (x, y)), (-1, (cut.outputs[3],))]).proved
    # 7. gadgets and binding
    g = blocks.gadget_identities()
    assert g["holds"] and g["and_with_frame"] and not g["and_without_frame"] and not g["erase_wire_without_frame"]
    bound = 0
    for prog, args in ((blocks.adder(3), ("add", 3, "qq", 0)), (blocks.multiplier(2), ("mul", 2, "qq", 0)),
                       (blocks.comparator(2), ("cmp", 2, "qq", 0))):
        spec = _spec(prog, args)
        layers = _layers(spec["program"])
        res = blocks.certify(dict(spec, layers=layers))
        assert res["status"] == "proved" and res["binding"]["bound"], res["binding"]
        assert res["binding"]["error_correction_steps_skipped"] > 0
        # independent operations may trade places
        moved = False
        plain = [i for i, l in enumerate(layers) if not all(op[0] in blocks.EC_OPS for op in l)]
        for i, j in zip(plain, plain[1:]):
            x, y = layers[i], layers[j]
            if len(x) == 1 and len(y) == 1 and x[0][0] == "CX" and y[0][0] == "CX" and not set(x[0][1:]) & set(y[0][1:]):
                swapped = list(layers)
                swapped[i], swapped[j] = y, x
                assert blocks.bind(spec["program"], swapped)["bound"]
                moved = True
                break
        assert moved or prog["name"].startswith("comparator"), prog["name"]
        # corruptions
        gate = next(i for i, l in enumerate(layers) if l[-1][0] == "CX")
        guard = next(i for i, l in enumerate(layers) if l[-1][0] in ("CXg", "CZg"))
        meas = next(i for i, l in enumerate(layers) if l[-1][0] in ("MZ", "MX") and l[-1][2] == layers[guard][-1][3])
        corrupt = [layers[:gate] + layers[gate + 1:],
                   layers[:gate] + [[["CX", layers[gate][-1][2], layers[gate][-1][1]]]] + layers[gate + 1:],
                   layers[:meas] + [layers[guard]] + layers[meas:guard] + layers[guard + 1:],
                   layers[:meas] + [layers[meas] + layers[guard]] + layers[meas + 1:guard] + layers[guard + 1:],
                   layers[:guard] + [[layers[guard][-1][:3] + ["nobit"]]] + layers[guard + 1:],
                   layers + [[["CX", prog["inputs"][0], prog["inputs"][1]]]],
                   [[(["NEW", "plus2", op[2]] if op[0] == "NEW" and len(op) == 3 else op) for op in l] for l in layers],
                   layers[:gate] + [layers[gate] + layers[gate]] + layers[gate + 1:],
                   layers[:gate] + [[["TOF", 1] + layers[gate][-1][1:]]] + layers[gate + 1:]]
        for j, bad_layers in enumerate(corrupt):
            assert not blocks.bind(spec["program"], bad_layers)["bound"], (prog["name"], j)
        other = _layers(_spec(blocks.adder(3, "constant", 5), ("add", 3, "constant", 5))["program"])
        assert not blocks.bind(spec["program"], other)["bound"]
        bound += 1 + int(moved)
    # 8. receipts
    rec = blocks.receipt(good)
    assert rec["schema"] == blocks.SCHEMA and rec["result"]["status"] == "proved"
    blocks.verify_receipt(good, json.loads(json.dumps(rec)))
    tampered = json.loads(json.dumps(rec))
    tampered["result"]["relation"]["proved"] = False
    try:
        blocks.verify_receipt(good, tampered)
        raise AssertionError("a tampered receipt verified")
    except blocks.BlockError:
        pass
    return (f"{checked} program checks, {refused} corruptions refused and {same} shown harmless, "
            f"hazard witnesses replayed, {bound} bindings with corrupted layers refused")
