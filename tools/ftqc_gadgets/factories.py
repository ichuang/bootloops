"""Purification trees as ancilla factories: a structural certificate of the filter rule.

A factory here is a preparation without inputs made of raw preparations of one
stabilizer state (the leaves) joined by checks.  A check of sector X applies a
transversal CX from a survivor block to a partner block, measures the partner in Z and
accepts when the outcome has the parities of the ideal state; a check of sector Z is
the mirror image (CX partner -> survivor, partner measured in X).

Enumeration (props.refine) decides such a gadget against its specification for small
fault counts.  The statement a level of concatenation needs is of another kind: any
number of faults that stay on their own coordinate are harmless, and only faults that
can spread are counted.  That statement follows from the shape of the circuit:

classes   every leaf location is, per sector, harmless, local at one coordinate (all
          its faults leave a pattern that is trivial or equivalent to an error at that
          coordinate) or junk-capable (anything else).
tree      the circuit is recognised as leaves and checks, operation by operation; the
          detectors of every check must span exactly the fixed parities of the state.
rule      the filter rule (LEMMA) then says, for every accepted fault set, where the
          output error of each sector can be: on the coordinates of the faulty local
          and step locations, unless enough leaves hold a junk-capable fault.

`certificate` returns the contract; `Tree.predict` evaluates the rule on one fault set
and `test` compares it with the propagated effect of that set, which is how the rule
itself is tested.
"""
import random
from itertools import combinations

from ftqc_audit.circuit import Circuit
from ftqc_codes import gf2, minweight

from .effects import Analysis
from .gadget import Block, Gadget, GadgetError

INF = 10 ** 9

LEMMA = (
    "Every leaf prepares the same stabilizer state; S_X and S_Z are the X-type and Z-type "
    "patterns that act trivially on it.  A leaf location is junk-capable in a sector if some "
    "fault there leaves a pattern not equivalent, modulo S, to an error on one coordinate, and "
    "local at c if every fault there leaves a pattern trivial or equivalent to one at c.  For a "
    "subtree T and a sector s define a flag J_s(T) and a set D_s(T) of coordinates.  Leaf: J_s "
    "is true when a junk-capable location is faulty, D_s is the set of coordinates of faulty "
    "local locations.  Check of sector s on survivor A and partner B: let F be the coordinates "
    "of faulty gate locations and of faulty idles of the survivor, and D_B+ the union of D_s(B) "
    "with the coordinates of faulty gate and measurement locations.  The transversal gate copies "
    "s patterns from A to B, and the check accepts exactly when the copied pattern plus that of "
    "B and of the step is in S_s (the detectors span the fixed parities of the state, whose "
    "annihilator is S_s).  So on acceptance the pattern A arrived with is equivalent to the one "
    "seen on B's side.  Hence J_s = J_s(A) and J_s(B); if J_s(A) is false, D_s = F + D_s(A); "
    "otherwise if J_s(B) is false, D_s = F + D_B+; and if both are false and |D_s(A) + D_B+| is "
    "below the minimum weight of S_s, the two descriptions agree and D_s = F + (D_s(A) & D_B+).  "
    "In the other sector patterns add: J = J(A) or J(B), D = D(A) + D(B) + F.  By induction, "
    "for every accepted fault set and every sector whose flag is false at the root, the output "
    "pattern is equivalent modulo S to one supported on D_s(root), which is contained in the "
    "coordinates of the faulty local and step locations, however many faults there are.  The "
    "flag of sector s can be true only if at least K_s leaves each hold a faulty junk-capable "
    "location, K_s computed by: leaf 1 (or none, if the leaf has no such location); check of "
    "sector s, K(A) + K(B); other sector, min(K(A), K(B)).  So with fewer than K_s such leaves "
    "in each sector the factory's output is its specification (ideal state, one noisy wire per "
    "qubit) with wire faults only on coordinates of faulty local and step locations.")


class TreeError(GadgetError):
    def __init__(self, message, where=None):
        super().__init__(message)
        self.where = where


def _span_equal(a, b):
    ea, eb = gf2.Echelon(a), gf2.Echelon(b)
    return len(ea.basis) == len(eb.basis) and all(v in eb for v in ea.basis)


class Tree:
    """The parsed factory.  Attributes: n, leaves (block -> leaf record), nodes (list of
    checks in circuit order), root (subtree id of the output), info (location index ->
    record used by predict)."""

    def __init__(self, g, strict=False, wmin=None, cap=3_000_000):
        if g.inputs or g.decoders or len(g.outputs) != 1:
            raise TreeError("a factory is a preparation: no inputs, no decoders, one output block")
        c = g.circuit
        if c.detectors is None:
            raise TreeError("a factory must declare its detectors")
        self.g, self.strict = g, strict
        out = g.outputs[0]
        n = self.n = len(out.qubits)
        if c.n % n:
            raise TreeError("qubits do not come in whole blocks")
        nb = c.n // n
        if out.qubits != list(range(out.qubits[0], out.qubits[0] + n)) or out.qubits[0] % n:
            raise TreeError("the output block is not one of the allocated blocks, in order")
        loc_index = {(k, j): i for i, (k, j, _) in enumerate(c.locations())}
        meas_index = {(m[0], m[1]): i for i, m in enumerate(c.meas)}
        phase = ["leaf"] * nb
        current = [("leaf", b) for b in range(nb)]
        leaf_ops = [[] for _ in range(nb)]
        pending = {}
        self.nodes = []
        for k, tick in enumerate(c.ticks):
            cross, single = {}, []
            for j, op in enumerate(tick):
                if op.name == "Z_IF":
                    raise TreeError("feedback is not part of a purification tree", (k, j))
                bs = [q // n for q in op.qs]
                if len(set(bs)) == 2:
                    if op.name != "CX" or op.qs[0] % n != op.qs[1] % n:
                        raise TreeError(f"{op.text()} joins two blocks and is not a transversal CX", (k, j))
                    cross.setdefault((bs[0], bs[1]), []).append((k, j, op))
                else:
                    single.append((bs[0], k, j, op))
            for (cb, tb), ops in cross.items():
                if sorted(op.qs[0] % n for _, _, op in ops) != list(range(n)):
                    raise TreeError(f"the gate between blocks {cb} and {tb} does not cover every coordinate",
                                    ops[0][:2])
                for b in (cb, tb):
                    if phase[b] == "dead" or b in pending:
                        raise TreeError(f"block {b} enters a second gate before its check is complete",
                                        ops[0][:2])
                    phase[b] = "alive"
                pend = dict(ctrl=cb, tgt=tb, gates=ops, idles={cb: [], tb: []})
                pending[cb] = pending[tb] = pend
            measured = {}
            crossing = {b for pair in cross for b in pair}
            for b, k_, j, op in single:
                if b in crossing:
                    raise TreeError(f"{op.text()} shares a tick with its block's transversal gate", (k_, j))
                if phase[b] == "leaf":
                    leaf_ops[b].append((k_, j, op))
                elif phase[b] == "dead":
                    raise TreeError(f"{op.text()} acts on a measured block", (k_, j))
                elif op.name in ("M", "MX"):
                    measured.setdefault(b, []).append((k_, j, op))
                elif op.name == "I":
                    if b in pending:
                        pending[b]["idles"][b].append((k_, j, op))
                    else:
                        self.nodes[current[b][1]]["idles"].append((k_, j, op))
                else:
                    raise TreeError(f"{op.text()} acts inside a block after its leaf", (k_, j))
            for b, ops in measured.items():
                pend = pending.get(b)
                names = {op.name for _, _, op in ops}
                if pend is None or len(ops) != n or len(names) != 1:
                    raise TreeError(f"block {b} is not measured whole, right after a gate", ops[0][:2])
                if b == pend["tgt"] and names == {"M"}:
                    sector, a = "X", pend["ctrl"]
                elif b == pend["ctrl"] and names == {"MX"}:
                    sector, a = "Z", pend["tgt"]
                else:
                    raise TreeError(f"the measurement of block {b} does not complete a check", ops[0][:2])
                if pend["idles"][b]:
                    raise TreeError(f"the partner block {b} idles between gate and measurement",
                                    pend["idles"][b][0][:2])
                if a in measured:
                    raise TreeError("both blocks of a check are measured", ops[0][:2])
                recs = {meas_index[(kk, jj)]: op.qs[0] % n for kk, jj, op in ops}
                self.nodes.append(dict(sector=sector, A=a, B=b, childA=current[a], childB=current[b],
                                       gates=pend["gates"], meas=ops, idles=list(pend["idles"][a]),
                                       records=recs))
                current[a] = ("node", len(self.nodes) - 1)
                phase[b] = "dead"
                del pending[a], pending[b]
        if pending:
            raise TreeError("a transversal gate is not followed by the measurement of one block")
        # a block that no operation touches is not part of the factory
        alive = [b for b in range(nb) if phase[b] != "dead" and (leaf_ops[b] or phase[b] != "leaf")]
        if alive != [out.qubits[0] // n]:
            raise TreeError("exactly the output block must be left unmeasured")
        self.root = current[alive[0]]
        self.loc_index = loc_index
        self._leaves(leaf_ops, out, cap)
        self._detectors(c)
        self._wmin(wmin, cap)
        self._info()

    # ------------------------------------------------------------------ leaves
    def _leaves(self, leaf_ops, out, cap):
        n, c = self.n, self.g.circuit
        shapes, self.leaves, self.shapes = {}, {}, []
        for b, ops in enumerate(leaf_ops):
            if not ops:
                if any(b in (nd["A"], nd["B"]) for nd in self.nodes):
                    raise TreeError(f"block {b} has no preparation")
                continue
            if any(op.name in ("M", "MX") for _, _, op in ops):
                raise TreeError(f"block {b} is measured inside its leaf", ops[0][:2])
            k0 = ops[0][0]
            key = tuple(sorted((k - k0, j_rank, op.name, op.noisy, tuple(q % n for q in op.qs))
                               for j_rank, (k, j, op) in enumerate(ops)))
            # the rank keeps the order of operations within a tick, which fixes location numbers
            if key not in shapes:
                shapes[key] = len(self.shapes)
                self.shapes.append(self._analyse(ops, k0, out.code, cap))
            locs = [self.loc_index[(k, j)] for k, j, op in ops if op.noisy]
            self.leaves[b] = dict(shape=shapes[key], locs=locs, ops=ops)
        first = self.shapes[0]
        for s in self.shapes[1:]:
            if s["canon"] != first["canon"]:
                raise TreeError("the leaves do not all prepare the same state")
        self.GX, self.GZ = first["GX"], first["GZ"]
        if out.code is not None:
            ex, ez = gf2.Echelon(self.GX), gf2.Echelon(self.GZ)
            if any(r not in ex for r in out.code.hx) or any(r not in ez for r in out.code.hz):
                raise TreeError("the prepared state is not in the declared code")

    def _analyse(self, ops, k0, code, cap):
        n = self.n
        depth = max(k for k, _, _ in ops) - k0 + 1
        ticks = [[] for _ in range(depth)]
        for k, j, op in ops:
            ticks[k - k0].append((("" if op.noisy else "!") + op.name,) + tuple(q % n for q in op.qs))
        try:
            circ = Circuit(ticks, n, detectors=[], outputs=list(range(n)), name="leaf")
            la = Analysis(Gadget(circ, [], [Block("out", list(range(n)))], [], 0, "leaf"), cap)
        except (ValueError, GadgetError) as e:
            raise TreeError(f"leaf: {e}") from None
        low = (1 << n) - 1
        GX = [v for v in la.canon if not v >> n]
        GZ = [v >> n for v in la.canon if not v & low]
        if len(GX) + len(GZ) != n or la.ns != n:
            raise TreeError("a leaf does not prepare a CSS stabilizer state")
        mask = dict(X=sum(1 << i for i, v in enumerate(la.canon) if not v & low),
                    Z=sum(1 << i for i, v in enumerate(la.canon) if not v >> n))
        per = dict(X={}, Z={})
        for q in range(n):
            if la.sx[q]:
                per["X"].setdefault(la.sx[q], q)
            if la.sz[q]:
                per["Z"].setdefault(la.sz[q], q)
        last = [op for k, _, op in ops if k - k0 == depth - 1]
        round_ok = (sorted(op.qs[0] % n for op in last if op.name == "I" and op.noisy) == list(range(n))
                    and len(last) == n)
        if self.strict and not round_ok:
            raise TreeError("strict classes need every leaf to end with one noisy idle on each qubit")
        classes, physical, single = [], True, []
        for loc in la.model.locations:
            in_round = round_ok and loc.tick == depth - 1
            row, lab = {}, {}
            for t in "XZ":
                where, junk = set(), False
                for label, eff in loc.faults.items():
                    part = la.sig(*eff) & mask[t]
                    if not part:
                        continue
                    cc = per[t].get(part)
                    if cc is None:
                        junk = True
                        continue
                    where.add(cc)
                    lab.setdefault((t, cc), label)
                if junk or len(where) > 1 or (self.strict and not in_round):
                    row[t] = ("junk", None)
                elif where:
                    cc = where.pop()
                    row[t] = ("local", cc)
                    if any(q != cc for q in loc.op.qs):
                        physical = False
                else:
                    row[t] = ("none", None)
            classes.append(row)
            single.append(lab)
        summ = dict(locations=len(classes))
        for t in "XZ":
            cnt = {}
            for row in classes:
                if row[t][0] == "local":
                    cnt[row[t][1]] = cnt.get(row[t][1], 0) + 1
            summ[f"junk_capable_{t}"] = sum(1 for row in classes if row[t][0] == "junk")
            summ[f"max_local_per_coordinate_{t}"] = max(cnt.values(), default=0)
        # single[i][(sector, c)]: a label of leaf location i that leaves exactly an error at c
        return dict(canon=list(la.canon), GX=GX, GZ=GZ, classes=classes, physical=physical, single=single,
                    summary=summ, round=round_ok)

    # ------------------------------------------------------------------ checks
    def _detectors(self, c):
        owner = {}
        for i, nd in enumerate(self.nodes):
            for r in nd["records"]:
                owner[r] = i
        masks = [[] for _ in self.nodes]
        for d in c.detectors:
            who = {owner.get(r) for r in d}
            if len(who) != 1 or None in who:
                raise TreeError("a detector is not a parity of one check's measurement")
            i = who.pop()
            masks[i].append(sum(1 << self.nodes[i]["records"][r] for r in d))
        for i, nd in enumerate(self.nodes):
            want = self.GZ if nd["sector"] == "X" else self.GX
            if not _span_equal(masks[i], want):
                raise TreeError(f"check {i} ({nd['sector']}): its detectors do not span the fixed "
                                f"parities of the state", nd["meas"][0][:2])

    def _wmin(self, wmin, cap):
        """Lower bounds on the minimum weight of S_X and S_Z (the X-type and Z-type
        stabilizer patterns of the state), proved by search; a request caps the search."""
        self.wmin, self.wmin_proof = {}, {}
        for t, G in (("X", self.GX), ("Z", self.GZ)):
            want = wmin.get(t) if isinstance(wmin, dict) else wmin
            if not G:
                self.wmin[t], self.wmin_proof[t] = INF, "no pattern acts trivially"
                continue
            if want is not None and (type(want) is not int or want < 1):
                raise TreeError("a minimum weight must be a positive integer")
            res = minweight.min_weight(G, self.n, target=want, max_combos=cap)
            lower, how = res["lower"], "search"
            if want is not None and lower < want:
                # the search stopped short: a cyclic group may still have its BCH bound
                from ftqc_codes import bch
                try:
                    r = bch.code_bound(G, self.n)
                    if r["bound"] > lower:
                        lower, how = r["bound"], "BCH bound: zeros computed exactly, theorem cited"
                except (bch.NotCyclic, ValueError):
                    pass
            if want is not None:
                if lower < want:
                    raise TreeError(f"minimum weight {want} of S_{t} was requested but only {lower} is proved")
                lower = want
            self.wmin[t] = lower
            self.wmin_proof[t] = dict(lower=lower if how != "search" else res["lower"],
                                      exact=bool(res.get("exact")) and how == "search")
            if how != "search":
                self.wmin_proof[t].update(search_lower=res["lower"], bound_from=how)

    # ------------------------------------------------------------------ rule
    def _info(self):
        self.info = {}
        for b, leaf in self.leaves.items():
            cls = self.shapes[leaf["shape"]]["classes"]
            for i, loc in enumerate(leaf["locs"]):
                self.info[loc] = ("leaf", b, cls[i])
        n = self.n
        for i, nd in enumerate(self.nodes):
            for kind in ("gates", "meas", "idles"):
                for k, j, op in nd[kind]:
                    if op.noisy:
                        self.info[self.loc_index[(k, j)]] = (kind, i, op.qs[0] % n)

    def predict(self, locations):
        """{sector: (junk flag, set of coordinates or None)} at the root for a set of
        faulty locations (labels do not matter: a location's class covers all of them)."""
        leaf_f, step = {}, {}
        for loc in locations:
            rec = self.info[loc]
            if rec[0] == "leaf":
                leaf_f.setdefault(rec[1], []).append(rec[2])
            else:
                step.setdefault(rec[1], dict(gates=set(), meas=set(), idles=set()))[rec[0]].add(rec[2])
        memo = {}

        def st(sub):
            if sub in memo:
                return memo[sub]
            kind, i = sub
            if kind == "leaf":
                rows = leaf_f.get(i, [])
                res = {t: (any(r[t][0] == "junk" for r in rows),
                           frozenset(r[t][1] for r in rows if r[t][0] == "local")) for t in "XZ"}
            else:
                nd = self.nodes[i]
                a, b = st(nd["childA"]), st(nd["childB"])
                f = step.get(i, dict(gates=set(), meas=set(), idles=set()))
                FA = frozenset(f["gates"] | f["idles"])
                res = {}
                for t in "XZ":
                    (ja, da), (jb, db) = a[t], b[t]
                    if t == nd["sector"]:
                        dbp = db | f["gates"] | f["meas"]
                        if not ja and not jb and len(da | dbp) < self.wmin[t]:
                            res[t] = (False, FA | (da & dbp))
                        elif not ja:
                            res[t] = (False, FA | da)
                        elif not jb:
                            res[t] = (False, FA | dbp)
                        else:
                            res[t] = (True, frozenset())
                    else:
                        res[t] = (ja or jb, da | db | FA)
            memo[sub] = res
            return res

        # iterative order: nodes were created children first
        for i in range(len(self.nodes)):
            st(("node", i))
        return st(self.root)

    def needed(self):
        """{sector: number of leaves that must hold a junk-capable fault for the flag to
        be true at the root} (INF: never)."""
        memo = {}
        for b, leaf in self.leaves.items():
            cls = self.shapes[leaf["shape"]]["classes"]
            memo[("leaf", b)] = {t: 1 if any(r[t][0] == "junk" for r in cls) else INF for t in "XZ"}
        for i, nd in enumerate(self.nodes):
            a, b = memo[nd["childA"]], memo[nd["childB"]]
            memo[("node", i)] = {t: min(INF, a[t] + b[t]) if t == nd["sector"] else min(a[t], b[t])
                                 for t in "XZ"}
        return memo[self.root]

    def sequence(self):
        """Sectors of the checks along the output block, first to last."""
        seq, sub = [], self.root
        while sub[0] == "node":
            seq.append(self.nodes[sub[1]]["sector"])
            sub = self.nodes[sub[1]]["childA"]
        return "".join(reversed(seq))

    def physical(self):
        return all(self.shapes[leaf["shape"]]["physical"] for leaf in self.leaves.values())


def certificate(g, strict=False, wmin=None, cap=3_000_000):
    """The factory contract of a gadget, or dict(tree=False, reason=...)."""
    out = dict(type="tree", strict=bool(strict))
    if wmin is not None:
        out["wmin_request"] = wmin
    try:
        tr = Tree(g, strict, wmin, cap)
    except TreeError as e:
        out.update(tree=False, reason=str(e))
        if e.where is not None:
            out["where"] = list(e.where)
        return out
    need = tr.needed()
    fin = {t: (None if v >= INF else v) for t, v in need.items()}
    least = min(need.values())
    out.update(tree=True, coordinates=tr.n, leaves=len(tr.leaves), checks=len(tr.nodes),
               sequence=tr.sequence(), stabilizers=dict(X=len(tr.GX), Z=len(tr.GZ)),
               wmin={t: (None if v >= INF else v) for t, v in tr.wmin.items()}, wmin_proof=tr.wmin_proof,
               leaf_kinds=[s["summary"] for s in tr.shapes],
               physical=tr.physical(), junk_leaves_needed=fin,
               through=None if least >= INF else least - 1)
    return out


def test(g, tree=None, exhaustive=1, samples=20000, seed=1, analysis=None, cap=3_000_000, big=40):
    """Compare the rule with propagation.  Every set of at most `exhaustive` locations
    (all labels), then `samples` random sets: a few faults, many faults on local and step
    locations only, and sets repaired to be accepted by adding one fault that cancels the
    detectors.  A failure is an accepted set for which a sector is not flagged and its
    part of the propagated effect is not a pattern on the predicted coordinates."""
    tr = tree or Tree(g, cap=cap)
    an = analysis or Analysis(g, cap)
    n = tr.n
    low = (1 << n) - 1
    mask = dict(X=sum(1 << i for i, v in enumerate(an.canon) if not v & low),
                Z=sum(1 << i for i, v in enumerate(an.canon) if not v >> n))
    per = dict(X=an.sx, Z=an.sz)
    faults = []                                    # (location, label, signature)
    for loc in an.model.locations:
        for label, eff in loc.faults.items():
            faults.append((loc.index, label, an.sig(*eff)))
    by_det = {}
    for f in faults:
        by_det.setdefault(f[2] & an.DET, []).append(f)
    out = dict(exhaustive_order=exhaustive, sets=0, accepted=0, flagged=0, failures=0,
               first_failure=None, largest_accepted=0)

    def judge(items):
        if len({f[0] for f in items}) != len(items):
            return
        sg = 0
        for f in items:
            sg ^= f[2]
        out["sets"] += 1
        if sg & an.DET:
            return
        out["accepted"] += 1
        out["largest_accepted"] = max(out["largest_accepted"], len(items))
        e = sg >> an.shift
        pred = tr.predict([f[0] for f in items])
        for t in "XZ":
            junk, dirt = pred[t]
            if junk:
                out["flagged"] += 1
                continue
            if (e & mask[t]) not in gf2.Echelon([per[t][c] for c in dirt]):
                out["failures"] += 1
                if out["first_failure"] is None:
                    out["first_failure"] = dict(sector=t, predicted=sorted(dirt),
                                                faults=[[f[0], f[1]] for f in items])

    for w in range(exhaustive + 1):
        for items in combinations(faults, w):
            judge(items)
    rnd = random.Random(seed)
    look = {(f[0], f[1]): f for f in faults}
    blocks = sorted(tr.leaves)
    first = [nd for nd in tr.nodes if nd["childA"][0] == "leaf" and nd["childB"][0] == "leaf"]

    def mirrored():
        """The same fault in a run of 2^k neighbouring leaves: what a check cannot see."""
        m = 1 << rnd.randint(1, max(1, len(blocks).bit_length() - 1))
        j = rnd.randrange(max(1, len(blocks) // m))
        group = blocks[j * m:(j + 1) * m]
        i = rnd.randrange(len(tr.leaves[group[0]]["locs"]))
        loc0 = tr.leaves[group[0]]["locs"][i]
        label = rnd.choice(sorted(an.model.locations[loc0].faults))
        return [look[(tr.leaves[b]["locs"][i], label)] for b in group
                if len(tr.leaves[b]["locs"]) > i and (tr.leaves[b]["locs"][i], label) in look]

    def split():
        """Single-coordinate errors on the two sides of a first check that add up to a
        pattern acting trivially: accepted, and the survivor keeps its own share."""
        nd = rnd.choice(first)
        t = nd["sector"]
        G = tr.GX if t == "X" else tr.GZ
        word = 0
        for v in rnd.sample(G, min(len(G), rnd.randint(1, 2))):
            word ^= v
        items = []
        for c in gf2.supp(word):
            b = rnd.choice((nd["A"], nd["B"]))
            leaf = tr.leaves[b]
            sg = tr.shapes[leaf["shape"]]["single"]
            cand = [(leaf["locs"][i], lab[(t, c)]) for i, lab in enumerate(sg) if (t, c) in lab]
            if not cand:
                return []
            items.append(look[rnd.choice(cand)])
        return items

    quiet = [f for f in faults if tr.info[f[0]][0] != "leaf"
             or all(tr.info[f[0]][2][t][0] != "junk" for t in "XZ")]
    for i in range(samples):
        mode = i % 6
        if mode == 4 and len(blocks) > 1:
            items = mirrored() + (rnd.sample(quiet, rnd.randint(0, 3)) if quiet else [])
        elif mode == 5 and first:
            items = split() + (rnd.sample(quiet, rnd.randint(0, 1)) if quiet and rnd.random() < 0.3 else [])
        elif mode == 0:
            items = rnd.sample(faults, rnd.randint(2, 5))
        elif mode == 1 and quiet:
            items = rnd.sample(quiet, min(len(quiet), rnd.randint(2, big)))
        elif mode == 2:
            items = rnd.sample(faults, rnd.randint(1, 4)) + (
                rnd.sample(quiet, min(len(quiet), rnd.randint(0, big // 2))) if quiet else [])
        else:
            items = rnd.sample(faults, rnd.randint(2, 6))
        sg = 0
        for f in items:
            sg ^= f[2]
        if sg & an.DET and mode in (1, 2, 3):
            fix = by_det.get(sg & an.DET)
            if fix:
                items = items + [rnd.choice(fix)]
        judge(items)
    return out
