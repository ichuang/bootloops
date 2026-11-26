"""Reference gadgets and the fragments they are made of.

Gadget specs (CODE is a qeccert spec without commas, e.g. steane, golay, or rep3):

  prep:CODE,STATE[,METHOD]     STATE zero | plus.  METHOD:
                                 raw        information-set preparation, not fault tolerant
                                 tree=SEQ   2^len(SEQ) raw copies purified pairwise; SEQ is a
                                            word in X and Z (X: CX kept->checked, measure the
                                            checked copy in Z; Z: CX checked->kept, measure in X)
                                 spec       noiseless preparation followed by one noisy wire
                                            per qubit (the specification of a preparation)
  cnot:CODE[,leaky]            transversal CX between two input blocks; `leaky` adds a pair
                               of cancelling CX inside the control block (negative control)
  idle:CODE                    one noisy wire per qubit on an input block
  knill:CODE[,anc=METHOD]      teleportation error correction: |+_L>,|0_L> ancillas by METHOD,
                               CX to make the encoded Bell pair, Bell measurement of the input
                               against one half, two record decoders
  steane_ec:CODE[,anc=METHOD]  syndrome extraction with |+_L> and |0_L> ancillas and two
                               Pauli-frame decoders on the data block
  seq_syndrome:CODE            the Z checks measured one after another, each with one bare
                               ancilla, and a Pauli-frame decoder for X errors (negative
                               control: an X fault between two checks is miscorrected)
  fc:SPEC                      a faultcert library circuit as a gadget without inputs
Suffixes: `,ideal` makes every operation noiseless and adds noisy wires on the outputs (the
specification of the gadget; for cnot one noisy two-qubit location per qubit pair, since one
fault of a transversal gate may leave an error on both blocks); `,idles` inserts idle locations.
"""
from faultcert import library as fclib
import qeccert
from qeccert import gf2

from .gadget import Builder, Block, Code, Gadget, GadgetError, prep_layers, state_checks

T_KNOWN = {"steane": 1, "golay": 3, "rep3": 1}


def code_of(spec):
    if not isinstance(spec, str):
        return spec
    if spec == "rep3":
        return Code("rep3", 3, [], [0b011, 0b110])
    try:
        return qeccert.load(spec)
    except Exception as e:
        raise GadgetError(f"unknown code {spec!r}: {e}") from None


# ---------------------------------------------------------------------- fragments
def frag_prep(b, code, state, method="raw"):
    """Emit a preparation of |0_L> or |+_L> at the builder's cursor; returns its qubits."""
    qs = b.alloc(code.n)
    if method == "raw":
        b.layers(prep_layers(code, state, qs))
    elif method == "ideal":
        b.layers(prep_layers(code, state, qs), False)
    elif method == "spec":
        b.layers(prep_layers(code, state, qs), False)
        b.wires(qs)
    elif method.startswith("tree="):
        seq = method[5:].upper()
        if not seq or set(seq) - set("XZ"):
            raise GadgetError(f"tree sequence must be a word in X and Z, got {seq!r}")
        blocks = b.par(*[(lambda: frag_prep(b, code, state, "raw")) for _ in range(2 ** len(seq))])
        qs = _purify(b, code, state, seq, blocks)
    else:
        raise GadgetError(f"unknown preparation method {method!r}")
    b.parts.append(dict(role=state, code=code.name, method=method))
    return qs


def _purify(b, code, state, seq, blocks):
    alive = blocks
    for s in seq:
        nxt, gate = [], []
        for i in range(0, len(alive), 2):
            A, B = alive[i], alive[i + 1]
            gate += [("CX", a, c) if s == "X" else ("CX", c, a) for a, c in zip(A, B)]
            nxt.append(A)
        b.tick(gate)
        recs = b.tick([("M" if s == "X" else "MX", c) for i in range(1, len(alive), 2) for c in alive[i]]
                      + [("I", a) for A in nxt for a in A])
        n = code.n
        masks = state_checks(code, state, "Z" if s == "X" else "X")
        for k in range(len(alive) // 2):
            b.checks(recs[k * n:(k + 1) * n], masks)
        alive = nxt
    return alive[0]


# ---------------------------------------------------------------------- gadgets
def prep(code, state="zero", method="raw", idles=False):
    code = code_of(code)
    b = Builder(f"prep({code.name},{state},{method})")
    b.no_detectors()
    qs = frag_prep(b, code, state, method)
    b.output(qs, code)
    return b.build(idles)


def cnot(code, leaky=False, ideal=False, idles=False):
    code = code_of(code)
    b = Builder(f"cnot({code.name}" + (",leaky" if leaky else "") + (",ideal" if ideal else "") + ")",
                ideal)
    A, B = b.input(code, "control"), b.input(code, "target")
    if leaky:
        b.tick([("CX", A[0], A[1])])
        b.tick([("CX", A[0], A[1])])
    b.tick([("CX", a, c) for a, c in zip(A, B)])
    if ideal:
        b.pair_wires(list(zip(A, B)))
    b.output(A, code, "control")
    b.output(B, code, "target")
    return b.build(idles)


def idle(code, idles=False):
    code = code_of(code)
    b = Builder(f"idle({code.name})")
    A = b.input(code, "data")
    b.wires(A)
    b.output(A, code, "data")
    return b.build(idles)


def knill(code, anc="raw", ideal=False, idles=False):
    code = code_of(code)
    b = Builder(f"knill({code.name},anc={anc}" + (",ideal" if ideal else "") + ")", ideal)
    A = b.input(code, "data")
    b.no_detectors()
    m = "ideal" if ideal else anc
    Bq, Cq = b.par(lambda: frag_prep(b, code, "plus", m), lambda: frag_prep(b, code, "zero", m))
    b.tick([("CX", x, y) for x, y in zip(Bq, Cq)])
    b.tick([("CX", x, y) for x, y in zip(A, Bq)])
    recs = b.tick([("MX", q) for q in A] + [("M", q) for q in Bq])
    n = code.n
    b.decoder("data_x_outcomes", recs[:n], code.hx, "record")
    b.decoder("ancilla_z_outcomes", recs[n:], code.hz, "record")
    if ideal:
        b.wires(Cq)
    b.output(Cq, code, "data")
    return b.build(idles)


def steane_ec(code, anc="raw", ideal=False, idles=False):
    code = code_of(code)
    b = Builder(f"steane_ec({code.name},anc={anc}" + (",ideal" if ideal else "") + ")", ideal)
    A = b.input(code, "data")
    b.no_detectors()
    m = "ideal" if ideal else anc
    P, Z = b.par(lambda: frag_prep(b, code, "plus", m), lambda: frag_prep(b, code, "zero", m))
    b.tick([("CX", a, p) for a, p in zip(A, P)])          # X errors of the data reach P
    r1 = b.tick([("M", p) for p in P] + [("CX", z, a) for z, a in zip(Z, A)])[:code.n]
    r2 = b.tick([("MX", z) for z in Z])
    b.decoder("x_errors", r1, code.hz, "X", A)
    b.decoder("z_errors", r2, code.hx, "Z", A)
    if ideal:
        b.wires(A)
    b.output(A, code, "data")
    return b.build(idles)


def seq_syndrome(code, idles=False):
    code = code_of(code)
    b = Builder(f"seq_syndrome({code.name})")
    A = b.input(code, "data")
    b.no_detectors()
    bits = []
    for row in code.hz:
        a = b.alloc(1)[0]
        b.tick([("R", a)])
        for p in gf2.supp(row):
            b.tick([("CX", A[p], a)])
        bits.append(b.tick([("M", a)]))
    b.syndrome_decoder("x_errors", bits, code.hz, "X", A)
    b.output(A, code, "data")
    return b.build(idles)


def from_circuit(circ, code=None):
    """A faultcert circuit as a gadget without inputs (one output block)."""
    blk = Block("out", circ.outputs, code)
    return Gadget(circ, [], [blk], [], 0, circ.name)


def idealize(g):
    """The specification of a gadget: the same circuit with every operation noiseless
    (input-error idles excepted), followed by one noisy wire on every output qubit."""
    from faultcert.circuit import Circuit
    c = g.circuit
    ticks = []
    for k, t in enumerate(c.ticks):
        row = []
        for op in t:
            keep = g.inputs and k == g.body_start and op.name == "I"
            row.append((op.name if keep else "!" + op.name,) + op.qs)
        ticks.append(row)
    outs = [q for b in g.outputs for q in b.qubits]
    ticks.append([("I", q) for q in outs])
    dets = [list(d) for d in c.detectors] if c.detectors is not None else None
    c2 = Circuit(ticks, c.n, detectors=dets, outputs=list(c.outputs), name=c.name + "_spec")
    # input-error locations keep their order: they are the first noisy operations
    ins, i = [], 0
    for b in g.inputs:
        ins.append(Block(b.name, b.qubits, b.code, b.refs, list(range(i, i + len(b.qubits)))))
        i += len(b.qubits)
    return Gadget(c2, ins, [Block(b.name, b.qubits, b.code) for b in g.outputs],
                  list(g.decoders), g.body_start, g.name + "_spec")


def load(spec):
    if spec.endswith(".json"):
        import json
        with open(spec) as f:
            return Gadget.from_json(json.load(f))
    fam, _, arg = spec.partition(":")
    if fam == "fc":
        return from_circuit(fclib.load(arg))
    parts = [p for p in arg.split(",") if p]
    if not parts:
        raise GadgetError(f"missing code in {spec!r}")
    flags = set(parts[1:])
    idles, ideal = "idles" in flags, "ideal" in flags
    flags -= {"idles", "ideal"}
    anc = next((f[4:] for f in flags if f.startswith("anc=")), "raw")
    flags = {f for f in flags if not f.startswith("anc=")}
    try:
        if fam == "prep":
            rest = [p for p in parts[1:] if p not in ("idles", "ideal")]
            state = rest[0] if rest else "zero"
            method = rest[1] if len(rest) > 1 else "raw"
            return prep(parts[0], state, method, idles)
        if fam == "cnot" and flags <= {"leaky"}:
            return cnot(parts[0], "leaky" in flags, ideal, idles)
        if fam == "idle" and not flags:
            return idle(parts[0], idles)
        if fam == "knill" and not flags:
            return knill(parts[0], anc, ideal, idles)
        if fam == "seq_syndrome" and not flags:
            return seq_syndrome(parts[0], idles)
        if fam == "steane_ec" and not flags:
            return steane_ec(parts[0], anc, ideal, idles)
    except IndexError:
        raise GadgetError(f"missing arguments in {spec!r}") from None
    raise GadgetError(f"unknown gadget {spec!r}; see gadgetcert.library")
