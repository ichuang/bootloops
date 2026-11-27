"""Signed projective Pauli instruments on an arbitrary stabilizer-code input."""

from functools import reduce
from itertools import combinations, product
from operator import xor
from types import SimpleNamespace

from ftqc_codes.gf2 import Echelon, dot, nullspace, supp
from ftqc_faults.audit.circuit import Circuit
from ftqc_faults.audit.model import propagate
from ftqc_faults.audit.stab import Gen, Sim


def pauli(word, n):
    if not isinstance(word, str) or len(word) != n + 1 or word[0] not in "+-" or set(word[1:]) - set("IXYZ"):
        raise ValueError("expected a signed Pauli, qubit zero first")
    x = sum(1 << i for i, c in enumerate(word[1:]) if c in "XY")
    z = sum(1 << i for i, c in enumerate(word[1:]) if c in "ZY")
    return Gen(x, z, int(word[0] == "-"))


def commute(a, b):
    return dot(a.x, b.z) == dot(a.z, b.x)


def centralizer(gens, n):
    mask = (1 << n) - 1
    return [Gen(v & mask, v >> n) for v in nullspace([g.z | g.x << n for g in gens], 2*n)]


class Group:
    """Echelon form with exact Pauli phases and symbolic outcome signs."""
    def __init__(self, width, gens=()):
        self.width = width
        self.rows = {}
        for g in gens:
            self.add(g)

    def reduce(self, g):
        r = Gen(g.x, g.z, g.s)
        while r.x or r.z:
            pivot = (r.x | r.z << self.width).bit_length()
            if pivot not in self.rows:
                break
            r.mul(self.rows[pivot])
        return r

    def add(self, g):
        r = self.reduce(g)
        if r.x or r.z:
            self.rows[(r.x | r.z << self.width).bit_length()] = r
            return True
        if r.s:
            raise ValueError("inconsistent stabilizer or wrong branch sign")
        return False


def paired(g, shift):
    # Transposition of each Y on the reference changes its sign.
    return Gen(g.x | g.x << shift, g.z | g.z << shift, dot(g.x, g.z))


def record_mask(ids, count):
    if not isinstance(ids, list) or any(type(i) != int or not 0 <= i < count for i in ids) or len(ids) != len(set(ids)):
        raise ValueError("invalid record indices")
    return sum(1 << i for i in ids)


def validate_basis(basis,checks,targets,n):
    """Check the named logical interface before and after the projective read."""
    if set(basis)!={"input","surviving"}:
        raise ValueError("expected input and surviving logical names")
    pairs={}
    previous=[]
    for item in basis["input"]:
        if set(item)!={"name","x","z"} or not isinstance(item["name"],str) or not item["name"] or item["name"] in pairs:
            raise ValueError("invalid logical pair name")
        a,b=(pauli(item[key],n) for key in ("x","z"))
        if commute(a,b) or any(not commute(v,g) for v in (a,b) for g in checks+previous):
            raise ValueError("invalid signed symplectic logical basis")
        pairs[item["name"]]=(a,b)
        previous.extend([a,b])
    survivors=basis["surviving"]
    k=n-len(Group(n,checks).rows)
    if len(pairs)!=k or not isinstance(survivors,list) or len(set(survivors))!=len(survivors) or len(survivors)!=k-len(targets):
        raise ValueError("logical interface has incorrect rank")
    if any(name not in pairs or any(not commute(v,t) for v in pairs[name] for t in targets) for name in survivors):
        raise ValueError("named survivor is disturbed by measurement")


class Instrument:
    def __init__(self, spec):
        if set(spec)-{"basis"} != {"n","checks","targets","ticks","readouts","detectors","terminal_syndrome"}:
            raise ValueError("invalid instrument fields")
        self.n = n = spec["n"]
        if type(n) != int or n < 1 or type(spec["terminal_syndrome"]) != bool:
            raise ValueError("invalid memory or boundary")
        self.checks = checks = [pauli(w, n) for w in spec["checks"]]
        self.targets = targets = [pauli(w, n) for w in spec["targets"]]
        generators = checks + targets
        if not targets or any(not commute(a,b) for a,b in combinations(generators,2)):
            raise ValueError("expected commuting checks and logical targets")
        group = Group(n, checks)
        if not all(group.add(t) for t in targets):
            raise ValueError("dependent or trivial logical targets")
        if "basis" in spec:
            validate_basis(spec["basis"],checks,targets,n)
        self.boundary = spec["terminal_syndrome"]
        self.residual = Echelon([g.x | g.z << n for g in generators])
        canonical={"I","CX","Z_IF","S","SDG","H","RX","R","MX","M","MZ"}
        for row in spec["ticks"]:
            for op in row:
                if not isinstance(op,(list,tuple)) or not op or not isinstance(op[0],str) or op[0].removeprefix("!") not in canonical or any(type(q)!=int or q<0 for q in op[1:]):
                    raise ValueError("invalid physical operation or wire")
        self.c = c = Circuit([[("!R", q) for q in range(n)], *spec["ticks"]], detectors=[])
        live = set(range(n))
        for tick in c.ticks[1:]:
            for op in tick:
                if op.name in ("R","RX"):
                    if op.qs[0] in live:
                        raise ValueError("reset of a live wire")
                    live.add(op.qs[0])
                elif not set(op.qs) <= live:
                    raise ValueError("operation on a released ancilla")
                if op.name in ("MX","M"):
                    if op.qs[0] < n:
                        raise ValueError("destructive data readout")
                    live.remove(op.qs[0])
        if live != set(range(n)):
            raise ValueError("unreleased workspace")
        self.detectors = [record_mask(d,len(c.meas)) for d in spec["detectors"]]
        if len(spec["readouts"]) != len(targets):
            raise ValueError("one signed readout per target")
        shift = c.n
        # Untouched physical reference block. Initial encoded Bell state is
        # fixed by both codes and by all paired centralizer generators.
        width = shift + n
        initial = Group(width)
        for g in checks:
            initial.add(g)
            initial.add(Gen(g.x << shift,g.z << shift,g.s ^ dot(g.x,g.z)))
        for g in centralizer(checks,n):
            initial.add(paired(g,shift))
        for q in range(n,shift):
            initial.add(Gen(0,1 << q))
        if len(initial.rows) != width:
            raise ValueError("incomplete Choi preparation")
        self.initial = list(initial.rows.values())
        self.readouts = []
        for read in spec["readouts"]:
            if set(read) != {"records","sign"} or type(read["sign"]) != int or read["sign"] not in (0,1):
                raise ValueError("invalid signed readout")
            self.readouts.append(record_mask(read["records"],len(c.meas)))
        if width >= 64:
            try:
                import stim  # noqa: F401 - optional accelerated tableau
            except ImportError:
                pass
            else:
                from .surgery_stim import compile_probe
                compile_probe(c,spec,self.initial)
                return
        sim = Sim(width)
        sim.gens = [Gen(g.x,g.z,g.s) for g in self.initial]
        sim.run(SimpleNamespace(ticks=c.ticks[1:]))
        final = Group(width,sim.gens)
        def require(g):
            # Group.reduce multiplies commuting Paulis only.
            if any(not commute(g,h) for h in sim.gens):
                raise ValueError("ideal output leaks or disturbs logical information")
            r = final.reduce(g)
            if r.x or r.z or r.s:
                raise ValueError("ideal output or reported sign differs")
        for g in checks:
            require(g)
            require(Gen(g.x << shift,g.z << shift,g.s ^ dot(g.x,g.z)))
        for g in centralizer(generators,n):
            require(paired(g,shift))
        variables = Echelon()
        for g, read in zip(targets,spec["readouts"],strict=True):
            if set(read) != {"records","sign"} or type(read["sign"]) != int or read["sign"] not in (0,1):
                raise ValueError("invalid signed readout")
            form = sim.parity_form(read["records"]) ^ read["sign"]
            if not variables.add(form >> 1):
                raise ValueError("ideal public outcomes lack full rank")
            require(Gen(g.x,g.z,g.s ^ form))
            require(Gen(g.x << shift,g.z << shift,g.s ^ form ^ dot(g.x,g.z)))
        if any(sim.parity_form(supp(mask)) for mask in self.detectors):
            raise ValueError("nondeterministic or negative acceptance parity")

    def response(self, x, z, rec):
        mask = (1 << self.n)-1
        residual = (x & mask) | (z & mask) << self.n
        det = sum(dot(rec,d) << i for i,d in enumerate(self.detectors))
        if self.boundary:
            for i,g in enumerate(self.checks,len(self.detectors)):
                det |= (dot(g.x,z) ^ dot(g.z,x)) << i
        read = sum(dot(rec,r) << i for i,r in enumerate(self.readouts))
        return det, read, self.residual.reduce(residual)

    def locations(self):
        mi = {(k,j):i for i,(k,j,_,_) in enumerate(self.c.meas)}
        for k,j,op in self.c.locations():
            if op.name in ("M","MX"):
                modes = [("flip",0,0,1 << mi[k,j])]
            else:
                modes = []
                for labels in product("IXZY", repeat=len(op.qs)):
                    if all(c == "I" for c in labels):
                        continue
                    x = sum(1 << q for c,q in zip(labels,op.qs) if c in "XY")
                    z = sum(1 << q for c,q in zip(labels,op.qs) if c in "ZY")
                    modes.append(("".join(labels),x,z,0))
            yield [(label,self.response(*propagate(self.c,k,x,z,r))) for label,x,z,r in modes]

    def audit(self, order, cap):
        if type(order) != int or type(cap) != int or order < 0 or cap < 1:
            raise ValueError("invalid fault search budget")
        from .surgery_faults import columns
        locs = columns(self)[0] if order else []
        tested = 0
        for w in range(1,min(order,len(locs))+1):
            for sites in combinations(range(len(locs)),w):
                for modes in product(*(locs[i] for i in sites)):
                    if tested == cap:
                        return dict(status="inconclusive",lower=w,tested=tested)
                    tested += 1
                    det,bit,res = [reduce(xor,items,0) for items in zip(*(v for _,v in modes))]
                    if not det and (bit or res):
                        return dict(status="counterexample",lower=w,upper=w,tested=tested,
                                    faults=[[i,label] for i,(label,_) in zip(sites,modes)],
                                    readout_flip=bit,residual=res)
        return dict(status="proved",lower=order+1,tested=tested)


def receipt(spec, order=1, cap=1000000):
    instrument = Instrument(spec)
    return dict(schema="ftqc_gadgets-surgery-1",spec=spec,order=order,cap=cap,
                circuit_sha256=instrument.c.sha256(),ideal="signed projective instrument",
                public_rank=len(instrument.targets),faults=instrument.audit(order,cap))


def verify(rec):
    rebuilt = receipt(rec["spec"],rec["order"],rec["cap"])
    import json
    if json.dumps(rebuilt,sort_keys=True) != json.dumps(rec,sort_keys=True):
        raise ValueError("surgery receipt differs from replay")
    return rebuilt

