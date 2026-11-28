"""Reference gadgets and the fragments they are made of.

Gadget specs (CODE is a ftqc_codes spec without commas, e.g. steane, golay, or rep3):

  prep:CODE,STATE[,METHOD]     STATE zero | plus.  METHOD:
                                 raw        information-set preparation, not fault tolerant
                                 goto       zero-postselection flagged Goto preparation with
                                            4-tick parallel dual checks and conditional hook
                                            table decoder (Steane [[7,1,3]])
                                 tree=SEQ   2^len(SEQ) raw copies purified pairwise; SEQ is a
                                            word in X and Z (X: CX kept->checked, measure the
                                            checked copy in Z; Z: CX checked->kept, measure in X)
                                 spec       noiseless preparation followed by one noisy wire
                                            per qubit (the specification of a preparation)
                               `,rounds=K` adds K ticks of idles after every raw preparation
  cnot:CODE[,leaky]            transversal CX between two input blocks; `leaky` adds a pair
                               of cancelling CX inside the control block (negative control)
  h:CODE                       transversal H on an input block
  phase_s:steane[,anc=goto]    zero-postselection cross-sector-compensated logical S on Steane
  idle:CODE                    one noisy wire per qubit on an input block
  knill:CODE[,anc=METHOD]      teleportation error correction: |+_L>,|0_L> ancillas by METHOD,
                               CX to make the encoded Bell pair, Bell measurement of the input
                               against one half, two record decoders
  steane_ec:CODE[,anc=METHOD][,rounds=K]
                               syndrome extraction with |+_L> and |0_L> ancillas and two
                               Pauli-frame decoders on the data block
  measure:CODE[,x|z]           destructive transversal measurement with a record decoder (no
                               output block: correct means the decoded outcome is right)
  exrec:CODE,ROLE[,ec=knill|steane_ec][,anc=METHOD][,trunc=0|1|0+1]
                               extended rectangle (leading ECs + role gadget + un-truncated
                               trailing ECs) for ROLE in prep_zero, prep_plus, idle, h, s,
                               cnot, measure_z, measure_x
  seq_syndrome:CODE            the Z checks measured one after another, each with one bare
                               ancilla, and a Pauli-frame decoder for X errors (negative
                               control: an X fault between two checks is miscorrected)
  fc:SPEC                      a ftqc_faults.audit library circuit as a gadget without inputs
Suffixes: `,ideal` makes every operation noiseless and adds noisy wires on the outputs (the
specification of the gadget; for cnot one noisy two-qubit location per qubit pair, since one
fault of a transversal gate may leave an error on both blocks); `,idles` inserts idle locations.
"""
from ftqc_faults.audit import library as fclib
import ftqc_codes
from ftqc_codes import gf2

from .gadget import Builder, Block, Code, Gadget, GadgetError, prep_layers, state_checks

T_KNOWN = {"steane": 1, "golay": 3, "rep3": 1}


def code_of(spec):
    if not isinstance(spec, str):
        return spec
    if spec == "rep3":
        return Code("rep3", 3, [], [0b011, 0b110])
    try:
        return ftqc_codes.load(spec)
    except Exception as e:
        raise GadgetError(f"unknown code {spec!r}: {e}") from None


# ---------------------------------------------------------------------- fragments
_GOTO_CACHE = {}


def _konig_color_4(edges):
    n_edges = len(edges)
    colors = [-1] * n_edges
    used_u, used_v = {}, {}

    def dfs(idx):
        if idx == n_edges:
            return True
        u, v = edges[idx]
        for c in range(4):
            if not used_u.get((u, c), False) and not used_v.get((v, c), False):
                colors[idx] = c
                used_u[(u, c)] = used_v[(v, c)] = True
                if dfs(idx + 1):
                    return True
                used_u[(u, c)] = used_v[(v, c)] = False
        return False

    if not dfs(0):
        raise GadgetError("bipartite 4-coloring of Steane checks failed")
    layers = [[] for _ in range(4)]
    for (u, v), c in zip(edges, colors):
        layers[c].append((u, v))
    return layers


def _min_coset_rep(mask, stab_rows, n=7):
    space = [0]
    for r in stab_rows:
        space += [x ^ r for x in space]
    best = (n + 1, mask)
    for s in space:
        c = mask ^ s
        cand = (c.bit_count(), c)
        if cand < best:
            best = cand
    return best[1]


def _goto_steane_spec(code, state):
    from itertools import permutations
    from ftqc_faults.audit.model import Model
    from ftqc_faults.zx import make_steane_7_1_3
    from ftqc_faults.zx.synthesize import synth_steane_goto_prep

    key = (code.n, tuple(code.hx), tuple(code.hz), state)
    cached = _GOTO_CACHE.get(key)
    if cached is not None:
        return cached
    if code.n != 7 or len(code.hx) != 3 or len(code.hz) != 3:
        raise GadgetError("method 'goto' requires the [[7,1,3]] Steane code")
    zxcode = make_steane_7_1_3()
    want = gf2.Echelon(code.hx)
    perm = None
    for cand in permutations(range(7)):
        mapped = [sum(1 << cand[q] for q in range(7) if (r >> q) & 1) for r in zxcode.hx_rows]
        if all(v in want for v in mapped) and gf2.rank(mapped) == 3:
            perm = cand
            break
    if perm is None:
        raise GadgetError("code is not coordinate-isomorphic to Steane [[7,1,3]]")

    def map_q(q):
        return perm[q] if q < 7 else 7

    goto = synth_steane_goto_prep()
    prep_ops = []
    cx_pairs = []
    meas_info = None
    for op in goto.ops:
        if op.kind in ("PREP_0", "PREP_PLUS"):
            mq = map_q(op.qubits[0])
            nm = (("R" if op.kind == "PREP_0" else "RX") if state == "zero"
                  else ("RX" if op.kind == "PREP_0" else "R"))
            prep_ops.append((nm, mq))
        elif op.kind == "CX":
            c, t = map_q(op.qubits[0]), map_q(op.qubits[1])
            if state == "plus":
                c, t = t, c
            cx_pairs.append((c, t))
        elif op.kind == "MEAS_Z":
            meas_info = ("M" if state == "zero" else "MX", map_q(op.qubits[0]))
    ticks_cx, last_t = [], {}
    for c, t in cx_pairs:
        earliest = max(last_t.get(c, -1), last_t.get(t, -1)) + 1
        while len(ticks_cx) <= earliest:
            ticks_cx.append([])
        ticks_cx[earliest].append((c, t))
        last_t[c] = last_t[t] = earliest

    rows_dual = list(code.hz) if state == "zero" else list(code.hx)
    rows_same = list(code.hx) if state == "zero" else list(code.hz)
    check_layers = _konig_color_4([(p, a_idx) for a_idx, r in enumerate(rows_dual) for p in gf2.supp(r)])

    # Probe the 8-qubit Goto core to extract the exact conditional hook table on flag==1.
    pb = Builder(f"probe_goto_{state}")
    pb.no_detectors()
    pqs = pb.alloc(8)
    pb.tick([(nm, pqs[mq]) for nm, mq in prep_ops])
    for layer in ticks_cx:
        pb.tick([("CX", pqs[c], pqs[t]) for c, t in layer])
    pb.tick([(meas_info[0], pqs[meas_info[1]])])
    pb.output(pqs[:7], code, "out")
    pg = pb.build(idles=True)
    pm = Model(pg.circuit)
    htab = {0: 0}
    for loc in pm.locations:
        for _lbl, (fx, fz, ff) in loc.faults.items():
            if not (ff & 1):
                continue
            src = fx if state == "zero" else fz
            raw_err = sum((1 << p) for p, q in enumerate(pqs[:7]) if (src >> q) & 1)
            crep = _min_coset_rep(raw_err, rows_same, 7)
            syn = sum(((crep & r).bit_count() & 1) << i for i, r in enumerate(rows_dual))
            if syn in htab and htab[syn] != crep:
                raise GadgetError(f"Goto hook table collision at syndrome {syn}")
            htab[syn] = crep
    for q in range(7):
        syn = sum((((1 << q) & r).bit_count() & 1) << i for i, r in enumerate(rows_dual))
        htab.setdefault(syn, 1 << q)

    spec = (prep_ops, ticks_cx, meas_info, check_layers, htab)
    _GOTO_CACHE[key] = spec
    return spec


def _emit_goto_prep_raw(b, code, state):
    if state not in ("zero", "plus"):
        raise GadgetError(f"unknown state {state!r} (zero or plus)")
    prep_ops, ticks_cx, meas_info, check_layers, htab = _goto_steane_spec(code, state)
    qs = b.alloc(8)
    b.tick([(nm, qs[mq]) for nm, mq in prep_ops])
    for layer in ticks_cx:
        b.tick([("CX", qs[c], qs[t]) for c, t in layer])
    flag_rec = b.tick([(meas_info[0], qs[meas_info[1]])])[0]
    data_qs = qs[:7]
    ancs = b.alloc(3)
    b.tick([("R" if state == "zero" else "RX", a) for a in ancs])
    for layer in check_layers:
        if state == "zero":
            b.tick([("CX", data_qs[p], ancs[a_idx]) for p, a_idx in layer])
        else:
            b.tick([("CX", ancs[a_idx], data_qs[p]) for p, a_idx in layer])
    syn_recs = b.tick([("M" if state == "zero" else "MX", a) for a in ancs])
    return data_qs, flag_rec, syn_recs, htab


def frag_prep(b, code, state, method="raw", rounds=0):
    """Emit a preparation of |0_L> or |+_L> at the builder's cursor; returns its qubits.
    rounds: noisy idle ticks on every qubit after each raw preparation (one level up, an
    idle is an inner error correction: a raw outer preparation ends with such rounds)."""
    if method == "goto":
        if rounds:
            raise GadgetError("rounds is for raw outer preparations, not goto")
        b.no_detectors()
        data_qs, flag_rec, syn_recs, htab = _emit_goto_prep_raw(b, code, state)
        kind_main = "X" if state == "zero" else "Z"
        entries = {(1 | (s << 1)): htab[s] for s in range(8) if htab[s]}
        tag = len(b.decoders)
        b.table_decoder(f"goto_hook_{kind_main.lower()}_{tag}",
                        [[flag_rec], [syn_recs[0]], [syn_recs[1]], [syn_recs[2]]],
                        kind_main, data_qs, entries, at_tick=b.cursor)
        b.parts.append(dict(role=state, code=code.name, method=method))
        return data_qs
    qs = b.alloc(code.n)
    if method == "raw":
        b.layers(prep_layers(code, state, qs))
        for _ in range(rounds):
            b.wires(qs)
    elif method == "ideal":
        b.layers(prep_layers(code, state, qs), False)
    elif method == "spec":
        b.layers(prep_layers(code, state, qs), False)
        b.wires(qs)
    elif method.startswith("tree="):
        seq = method[5:].upper()
        if not seq or set(seq) - set("XZ"):
            raise GadgetError(f"tree sequence must be a word in X and Z, got {seq!r}")
        blocks = b.par(*[(lambda: frag_prep(b, code, state, "raw", rounds)) for _ in range(2 ** len(seq))])
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
def prep(code, state="zero", method="raw", idles=False, rounds=0):
    code = code_of(code)
    b = Builder(f"prep({code.name},{state},{method}" + (f",rounds={rounds}" if rounds else "") + ")")
    b.no_detectors()
    qs = frag_prep(b, code, state, method, rounds)
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


def steane_ec(code, anc="raw", ideal=False, idles=False, rounds=1):
    """rounds > 1 repeats the extraction; each later round decodes its check values
    relative to the frame the previous round left (the parity of the two rounds'
    outcomes), as a frame-tracking implementation does."""
    code = code_of(code)
    b = Builder(f"steane_ec({code.name},anc={anc}" + (",ideal" if ideal else "")
                + (f",rounds={rounds}" if rounds != 1 else "") + ")", ideal)
    A = b.input(code, "data")
    b.no_detectors()
    m = "ideal" if ideal else anc
    last = None
    for k in range(rounds):
        P, Z = b.par(lambda: frag_prep(b, code, "plus", m), lambda: frag_prep(b, code, "zero", m))
        b.tick([("CX", a, p) for a, p in zip(A, P)])          # X errors of the data reach P
        r1 = b.tick([("M", p) for p in P] + [("CX", z, a) for z, a in zip(Z, A)])[:code.n]
        r2 = b.tick([("MX", z) for z in Z])
        if rounds == 1:
            b.decoder("x_errors", r1, code.hz, "X", A)
            b.decoder("z_errors", r2, code.hx, "Z", A)
        else:
            for nm, recs, old, H, kind in (("x_errors", r1, last and last[0], code.hz, "X"),
                                           ("z_errors", r2, last and last[1], code.hx, "Z")):
                bits = [[recs[p] for p in gf2.supp(row)] + ([old[p] for p in gf2.supp(row)] if old else [])
                        for row in H]
                b.syndrome_decoder(f"{nm}_{k + 1}", bits, H, kind, A)
        last = (r1, r2)
    if ideal:
        b.wires(A)
    b.output(A, code, "data")
    return b.build(idles)


def measure(code, basis="Z", idles=False):
    """Destructive transversal measurement of an input block, decoded.  No output block:
    the effect is zero exactly when the decoded logical outcomes are the right ones."""
    code = code_of(code)
    basis = basis.upper()
    if basis not in ("X", "Z"):
        raise GadgetError("measurement basis is X or Z")
    b = Builder(f"measure({code.name},{basis})")
    A = b.input(code, "data")
    b.no_detectors()
    recs = b.measure(A, basis)
    b.decoder("outcomes", recs, code.hz if basis == "Z" else code.hx, "record")
    return b.build(idles)


def h_gate(code, ideal=False, idles=False):
    code = code_of(code)
    b = Builder(f"h({code.name}" + (",ideal" if ideal else "") + ")", ideal)
    A = b.input(code, "data")
    b.no_detectors()
    b.tick([("H", q) for q in A])
    if ideal:
        b.wires(A)
    b.output(A, code, "data")
    return b.build(idles)


def frag_phase_s(b, code, A, anc="goto", causal=False):
    if anc != "goto":
        raise GadgetError("phase_s requires anc='goto' on Steane [[7,1,3]]")
    lead_hx, lead_hz = {0: 0}, {0: 0}
    for q in range(7):
        shx = sum((((1 << q) & r).bit_count() & 1) << i for i, r in enumerate(code.hx))
        shz = sum((((1 << q) & r).bit_count() & 1) << i for i, r in enumerate(code.hz))
        lead_hx[shx] = 1 << q
        lead_hz[shz] = 1 << q

    def syn_hx(mask):
        return sum((((mask & r).bit_count() & 1) << i) for i, r in enumerate(code.hx))

    def syn_hz(mask):
        return sum((((mask & r).bit_count() & 1) << i) for i, r in enumerate(code.hz))

    (B1, fp1, sp1, htab_p1), (C1, fz1, sz1, htab_z1), (B2, fp2, sp2, htab_p2) = b.par(
        lambda: _emit_goto_prep_raw(b, code, "plus"),
        lambda: _emit_goto_prep_raw(b, code, "zero"),
        lambda: _emit_goto_prep_raw(b, code, "plus"),
    )
    b.tick([("CX", a, x) for a, x in zip(A, B1)])
    rz1 = b.measure(B1, "Z", extra=[("I", q) for q in A + C1 + B2])
    b.tick([("CX", x, a) for x, a in zip(C1, A)])
    rx1 = b.measure(C1, "X", extra=[("SDG", q) for q in A] + [("I", q) for q in B2])
    b.tick([("CX", a, x) for a, x in zip(A, B2)])
    rz2 = b.measure(B2, "Z", extra=[("I", q) for q in A])

    hz1_bits = [[rz1[p] for p in gf2.supp(r)] for r in code.hz]
    hx1_bits = [[rx1[p] for p in gf2.supp(r)] for r in code.hx]
    hz2_bits = [[rz2[p] for p in gf2.supp(r)] for r in code.hz]
    at = b.cursor if causal else None
    pfx = f"s{len(b.decoders)}_"

    tab_x = {}
    for val in range(128):
        xh = htab_z1[(val >> 1) & 7] if (val & 1) else 0
        tab_x[val] = xh ^ lead_hz[((val >> 4) & 7) ^ syn_hz(xh)]
    b.table_decoder(pfx + "x_post_sdg",
                    [[fz1], [sz1[0]], [sz1[1]], [sz1[2]], hz2_bits[0], hz2_bits[1], hz2_bits[2]],
                    "X", A, tab_x, at_tick=at)

    tab_z1 = {}
    for val in range(128):
        zp = htab_p1[(val >> 1) & 7] if (val & 1) else 0
        tab_z1[val] = zp ^ lead_hx[((val >> 4) & 7) ^ syn_hx(zp)]
    b.table_decoder(pfx + "z_pre_sdg",
                    [[fp1], [sp1[0]], [sp1[1]], [sp1[2]], hx1_bits[0], hx1_bits[1], hx1_bits[2]],
                    "Z", A, tab_z1, at_tick=at)

    tab_b2h = {(1 | (s << 1)): htab_p2[s] for s in range(8) if htab_p2[s]}
    b.table_decoder(pfx + "b2_hook_z", [[fp2], [sp2[0]], [sp2[1]], [sp2[2]]],
                    "Z", A, tab_b2h, at_tick=at)

    tab_xz = {}
    for val in range(1 << 13):
        fz_b = val & 1
        sz_v = (val >> 1) & 7
        hz1_v = (val >> 4) & 7
        hx1_v = (val >> 7) & 7
        hz2_v = (val >> 10) & 7
        if fz_b:
            tab_xz[val] = htab_z1[sz_v]
        elif hz1_v or hx1_v:
            tab_xz[val] = lead_hz[hz2_v]
    b.table_decoder(pfx + "z_from_x_across_sdg",
                    [[fz1], [sz1[0]], [sz1[1]], [sz1[2]],
                     hz1_bits[0], hz1_bits[1], hz1_bits[2],
                     hx1_bits[0], hx1_bits[1], hx1_bits[2],
                     hz2_bits[0], hz2_bits[1], hz2_bits[2]],
                    "Z", A, tab_xz, at_tick=at)
    return A


def phase_s(code, anc="goto", ideal=False, idles=False):
    code = code_of(code)
    if ideal:
        b = Builder(f"phase_s({code.name},ideal)", True)
        A = b.input(code, "data")
        b.no_detectors()
        b.tick([("SDG", q) for q in A])
        b.wires(A)
        b.output(A, code, "data")
        return b.build(idles)
    b = Builder(f"phase_s({code.name},anc={anc})")
    A = b.input(code, "data")
    b.no_detectors()
    frag_phase_s(b, code, A, anc=anc, causal=False)
    b.output(A, code, "data")
    return b.build(idles)


def frag_ec(b, code, A, ec="knill", anc="goto"):
    n = code.n
    tag = len(b.decoders)
    if ec == "knill":
        Bq, Cq = b.par(lambda: frag_prep(b, code, "plus", anc),
                       lambda: frag_prep(b, code, "zero", anc))
        b.tick([("CX", x, y) for x, y in zip(Bq, Cq)])
        b.tick([("CX", x, y) for x, y in zip(A, Bq)])
        recs = b.tick([("MX", q) for q in A] + [("M", q) for q in Bq])
        b.decoder(f"knill_ax_{tag}", recs[:n], code.hx, "record")
        b.decoder(f"knill_bz_{tag}", recs[n:], code.hz, "record")
        return Cq
    if ec == "steane_ec":
        P, Z = b.par(lambda: frag_prep(b, code, "plus", anc),
                     lambda: frag_prep(b, code, "zero", anc))
        b.tick([("CX", a, p) for a, p in zip(A, P)])
        r1 = b.tick([("M", p) for p in P] + [("CX", z, a) for z, a in zip(Z, A)])[:n]
        r2 = b.tick([("MX", z) for z in Z])
        b.decoder(f"steane_x_{tag}", r1, code.hz, "X", A, at_tick=b.cursor)
        b.decoder(f"steane_z_{tag}", r2, code.hx, "Z", A, at_tick=b.cursor)
        return A
    raise GadgetError(f"unknown ec method {ec!r} (knill or steane_ec)")


def exrec(code, role, ec="knill", anc="goto", trunc=(), idles=True):
    """Extended rectangle ExRec_T(role) with leading ECs and un-truncated trailing ECs."""
    code = code_of(code)
    tset = {int(x) for x in trunc}
    valid_outs = {
        "prep_zero": {0}, "prep_plus": {0}, "idle": {0}, "h": {0}, "s": {0},
        "cnot": {0, 1}, "measure_z": set(), "measure_x": set(),
    }
    if role not in valid_outs:
        raise GadgetError(f"unknown exrec role {role!r}")
    if not (tset <= valid_outs[role]):
        raise GadgetError(f"invalid truncation indices {sorted(tset)} for role {role!r}")
    tsuf = f",trunc={'+'.join(str(x) for x in sorted(tset))}" if tset else ""
    b = Builder(f"exrec({code.name},{role},ec={ec},anc={anc}{tsuf})")
    b.no_detectors()
    if role in ("prep_zero", "prep_plus"):
        qs = frag_prep(b, code, "zero" if role == "prep_zero" else "plus", anc)
        if 0 not in tset:
            qs = frag_ec(b, code, qs, ec, anc)
        b.output(qs, code, "out")
    elif role in ("idle", "h", "s"):
        A = b.input(code, "data")
        A = frag_ec(b, code, A, ec, anc)
        if role == "idle":
            b.wires(A)
        elif role == "h":
            b.tick([("H", q) for q in A])
        else:
            A = frag_phase_s(b, code, A, anc=anc, causal=(0 not in tset))
        if 0 not in tset:
            A = frag_ec(b, code, A, ec, anc)
        b.output(A, code, "data")
    elif role == "cnot":
        A, B = b.input(code, "control"), b.input(code, "target")
        A, B = b.par(lambda: frag_ec(b, code, A, ec, anc),
                     lambda: frag_ec(b, code, B, ec, anc))
        b.tick([("CX", a, c) for a, c in zip(A, B)])
        if 0 not in tset and 1 not in tset:
            A, B = b.par(lambda: frag_ec(b, code, A, ec, anc),
                         lambda: frag_ec(b, code, B, ec, anc))
        elif 0 not in tset:
            A = frag_ec(b, code, A, ec, anc)
        elif 1 not in tset:
            ec_q = (2 * (code.n + 4)) if anc == "goto" else (2 * code.n)
            _, B = b.par(lambda: b.alloc(ec_q), lambda: frag_ec(b, code, B, ec, anc))
        b.output(A, code, "control")
        b.output(B, code, "target")
    else:
        A = b.input(code, "data")
        A = frag_ec(b, code, A, ec, anc)
        basis = "Z" if role == "measure_z" else "X"
        recs = b.measure(A, basis)
        b.decoder("outcomes", recs, code.hz if basis == "Z" else code.hx, "record")
    return b.build(idles)


def radii(gadget, code_receipt):
    """Radius of each decoder of a gadget from an ftqc_codes receipt of its block code: a
    decoder whose checks are H_Z separates X patterns up to (m - 1) // 2, m the smaller of
    the proven lower bounds on d_x and on the lightest X stabilizer (two patterns with the
    same check values differ by a nonzero element of ker H_Z, a logical or a stabilizer);
    likewise for H_X with d_z."""
    from ftqc_codes import certify as cc
    cc.verify(code_receipt)
    code = ftqc_codes.CSSCode.from_json(code_receipt["code"])
    f = code_receipt["facts"]
    lo = {q: f[q]["lower"] for q in ("d_x", "d_z", "wmin_x", "wmin_z")}
    out = {}
    for d in gadget.decoders:
        if d.checks == list(code.hz):
            out[d.name] = (min(lo["d_x"], lo["wmin_x"]) - 1) // 2
        elif d.checks == list(code.hx):
            out[d.name] = (min(lo["d_z"], lo["wmin_z"]) - 1) // 2
        else:
            raise GadgetError(f"decoder {d.name}: its checks are not H_X or H_Z of the receipt's code")
    return out


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
    """A ftqc_faults.audit circuit as a gadget without inputs (one output block)."""
    blk = Block("out", circ.outputs, code)
    return Gadget(circ, [], [blk], [], 0, circ.name)


def idealize(g):
    """The specification of a gadget: the same circuit with every operation noiseless
    (input-error idles excepted), followed by one noisy wire on every output qubit."""
    from ftqc_faults.audit.circuit import Circuit
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
    rounds = next((f[7:] for f in flags if f.startswith("rounds=")), None)
    if rounds is not None and not rounds.isdigit():
        raise GadgetError(f"rounds takes a number in {spec!r}")
    rounds = None if rounds is None else int(rounds)
    flags = {f for f in flags if not f.startswith("anc=") and not f.startswith("rounds=")}
    ec_opt = next((f[3:] for f in flags if f.startswith("ec=")), "knill")
    trunc_opt = next((f[6:] for f in flags if f.startswith("trunc=")), None)
    flags = {f for f in flags if not f.startswith("ec=") and not f.startswith("trunc=")}
    try:
        if fam == "prep":
            rest = [p for p in parts[1:] if p not in ("idles", "ideal") and not p.startswith("rounds=")]
            state = rest[0] if rest else "zero"
            method = rest[1] if len(rest) > 1 else "raw"
            return prep(parts[0], state, method, idles, rounds or 0)
        if fam == "cnot" and flags <= {"leaky"}:
            return cnot(parts[0], "leaky" in flags, ideal, idles)
        if fam == "h" and not flags:
            return h_gate(parts[0], ideal, idles)
        if fam == "phase_s" and not flags:
            return phase_s(parts[0], anc if anc != "raw" else "goto", ideal, idles)
        if fam == "idle" and not flags:
            return idle(parts[0], idles)
        if fam == "knill" and not flags:
            return knill(parts[0], anc, ideal, idles)
        if fam == "seq_syndrome" and not flags:
            return seq_syndrome(parts[0], idles)
        if fam == "steane_ec" and not flags:
            return steane_ec(parts[0], anc, ideal, idles, 1 if rounds is None else rounds)
        if fam == "measure" and flags <= {"x", "z", "X", "Z"}:
            return measure(parts[0], (sorted(flags) or ["Z"])[0], idles)
        if fam == "exrec":
            rest = [p for p in parts[1:] if p not in ("idles", "ideal")
                    and not p.startswith("anc=") and not p.startswith("ec=") and not p.startswith("trunc=")]
            if len(rest) != 1:
                raise GadgetError(f"exrec takes CODE,ROLE[,ec=...][,anc=...][,trunc=...], got {spec!r}")
            tr = () if not trunc_opt else tuple(int(x) for x in trunc_opt.split("+") if x != "")
            use_idles = True if "noidles" not in flags else False
            return exrec(parts[0], rest[0], ec=ec_opt, anc=(anc if anc != "raw" else "goto"),
                         trunc=tr, idles=use_idles)
    except IndexError:
        raise GadgetError(f"missing arguments in {spec!r}") from None
    raise GadgetError(f"unknown gadget {spec!r}; see ftqc_gadgets.library")
