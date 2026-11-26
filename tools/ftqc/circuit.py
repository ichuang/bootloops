"""Nondestructive Z measurements on an arbitrary HGP input codespace."""
from itertools import combinations, product
from algebra import (bits, echelon, hgp, kernel, parity, reduce_row, seed, Limit)


class Schedule:
    def __init__(self,spec):
        a,na = seed(spec['a'])
        b,nb = seed(spec['b'])
        self.n,x,z = hgp(a,na,b,nb)
        self.ops = spec['ops']
        m = spec['measured_z']
        if not 0<m<1<<self.n or any(parity(m,r) for r in x) or not reduce_row(m,echelon(z)):
            raise ValueError("target is not a nontrivial Z logical")
        live = set(range(self.n))
        self.records = []
        self.modes = []
        for t,op in enumerate(self.ops):
            g,*qs = op
            if g not in ('RZ','RX','MZ','MX','CX','H','I') or len(qs)!=(2 if g=='CX' else 1) or len(set(qs))!=len(qs) or any(type(q)!=int or not 0<=q<100000 for q in qs):
                raise ValueError("invalid scheduled operation")
            q=qs[0]
            if g.startswith('R'):
                if q in live or q<self.n:
                    raise ValueError("reset of live qubit")
                live.add(q)
            elif any(q not in live for q in qs):
                raise ValueError("operation on inactive qubit")
            if g.startswith('M'):
                if q<self.n:
                    raise ValueError("data measurement at open output boundary")
                self.records.append(t)
                live.remove(q)
            self.modes.append(range(1,16 if g=='CX' else 2 if g[0] in 'RM' else 4))
        if live!=set(range(self.n)):
            raise ValueError("unclosed ancilla")
        def recmask(ids):
            if len(set(ids))!=len(ids) or any(type(i)!=int or not 0<=i<len(self.records) for i in ids):
                raise ValueError("bad measurement record")
            return sum(1<<self.records[i] for i in ids)
        # Each tuple is an output Pauli X,Z and a parity of measurement records.
        self.det = [(r,0,0) for r in x]+[(0,r,0) for r in z]
        self.det += [(0,0,recmask(ids)) for ids in spec['detectors']]
        mx = kernel(z+[m],self.n)
        mz = kernel(x,self.n)
        self.obs = [(0,0,recmask(spec['readout']))]+[(r,0,0) for r in mx]+[(0,r,0) for r in mz]
        self.flows = self.det+self.obs
        for i,flow in enumerate(self.flows):
            ox,oz,phase = self.pull(flow)
            ex,ez = (0,0) if i<len(self.det) else (0,m) if i==len(self.det) else flow[:2]
            if phase or reduce_row(ox^ex,echelon(x)) or reduce_row(oz^ez,echelon(z)):
                raise ValueError("ideal detector, readout, or logical flow fails")
        self.mask = (1<<len(self.det))-1

    def pull(self,flow,faults=None):
        """Heisenberg pullback, including phase and destructive measurement."""
        x,z,rec=flow
        phase=0
        for t in reversed(range(len(self.ops))):
            g,q,*tail=self.ops[t]
            bit=1<<q
            mode=(faults or {}).get(t,0)
            if g.startswith('M'):
                if rec>>t&1:
                    phase ^= bool(mode)
                    if g=='MX': x ^= bit
                    else: z ^= bit
                continue
            fx,fz=self.pauli(t,mode)
            phase ^= parity(fx,z)^parity(fz,x)
            if g=='CX':
                target=1<<tail[0]
                if x&bit: x ^= target
                if z&target: z ^= bit
            elif g=='H':
                px,pz=bool(x&bit),bool(z&bit)
                phase ^= px and pz
                if px!=pz: x ^= bit; z ^= bit
            elif g.startswith('R'):
                if (z if g=='RX' else x)&bit:
                    raise ValueError("flow has random reset boundary")
                x &= ~bit
                z &= ~bit
        return x,z,phase

    def pauli(self,t,mode):
        g,q,*tail=self.ops[t]
        if g=='RX': return 0,mode<<q
        if g=='RZ': return mode<<q,0
        x,z=(mode&1)<<q,((mode>>1)&1)<<q
        if tail:
            x |= ((mode>>2)&1)<<tail[0]
            z |= ((mode>>3)&1)<<tail[0]
        return x,z

    def effect(self,faults):
        """Forward propagation independent of the reverse flow check."""
        x=z=rec=0
        for t,(g,q,*tail) in enumerate(self.ops):
            bit=1<<q
            mode=faults.get(t,0)
            if g.startswith('M'):
                flip=bool((z if g=='MX' else x)&bit)^bool(mode)
                rec |= int(flip)<<t
                x &= ~bit
                z &= ~bit
                continue
            if g=='CX':
                target=1<<tail[0]
                if x&bit: x ^= target
                if z&target: z ^= bit
            if g=='H' and bool(x&bit)!=bool(z&bit):
                x ^= bit
                z ^= bit
            fx,fz=self.pauli(t,mode)
            x ^= fx
            z ^= fz
        return sum((parity(x,b)^parity(z,a)^parity(rec,r))<<i for i,(a,b,r) in enumerate(self.flows))

    def signature(self,faults):
        return sum(int(self.pull(flow,faults)[2])<<i for i,flow in enumerate(self.flows))

    def events(self,reverse=False):
        evaluate=self.signature if reverse else self.effect
        return [[(mode,evaluate({t:mode})) for mode in modes] for t,modes in enumerate(self.modes)]


def search(schedule,target,budget,reverse=False):
    events=schedule.events(reverse)
    tested=0
    for w in range(1,min(target,len(events)+1)):
        for locations in combinations(range(len(events)),w):
            for choices in product(*(events[i] for i in locations)):
                tested+=1
                if tested>budget:
                    return dict(status='inconclusive',lower=w,reason='fault enumeration budget',tested=tested-1)
                effect=0
                for _,s in choices: effect ^= s
                if effect and not effect&schedule.mask:
                    faults=[[t,choice[0]] for t,choice in zip(locations,choices)]
                    if schedule.signature(dict(faults))!=effect:
                        raise ValueError("forward/reverse fault disagreement")
                    return dict(status='counterexample',lower=w,upper=w,faults=faults,flipped=[i-len(schedule.det) for i in bits(effect)],tested=tested)
    return dict(status='proved',lower=target,tested=tested)
