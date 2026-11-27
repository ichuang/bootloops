"""Reverse propagation of complete detector, public-bit and residual signatures."""
from itertools import combinations,product
from functools import reduce
from operator import xor
from math import comb


def columns(inst):
    c=inst.c
    nd=len(inst.detectors)+(len(inst.checks) if inst.boundary else 0)
    nr=len(inst.readouts)
    def pack(x,z,r):
        d,b,p=inst.response(x,z,r)
        return d|(b<<nd)|(p<<(nd+nr))
    def unpack(s):
        return s&((1<<nd)-1),(s>>nd)&((1<<nr)-1),s>>(nd+nr)
    x=[pack(1<<q,0,0) for q in range(c.n)]
    z=[pack(0,1<<q,0) for q in range(c.n)]
    rec=[pack(0,0,1<<i) for i in range(len(c.meas))]
    mi=len(rec)
    result=[]
    for tick in reversed(c.ticks[1:]):
        for op in reversed(tick):
            q=op.qs[0]
            name=op.name
            if name in ("M","MX"):
                mi-=1
                modes=[("flip",unpack(rec[mi]))]
                x[q],z[q]=(rec[mi],0) if name=="M" else (0,rec[mi])
            else:
                modes=[]
                if op.noisy:
                    for labels in product("IXZY",repeat=len(op.qs)):
                        if any(p!="I" for p in labels):
                            value=0
                            for p,i in zip(labels,op.qs):
                                if p in "XY":
                                    value^=x[i]
                                if p in "ZY":
                                    value^=z[i]
                            modes.append(("".join(labels),unpack(value)))
                if name=="CX":
                    t=op.qs[1]
                    x[q]^=x[t]
                    z[t]^=z[q]
                elif name=="H":
                    x[q],z[q]=z[q],x[q]
                elif name in ("S","SDG"):
                    x[q]^=z[q]
                elif name=="Z_IF":
                    rec[op.record]^=z[q]
                elif name in ("R","RX"):
                    x[q]=z[q]=0
            if op.noisy:
                result.append(modes)
    incoming=[[("X",unpack(x[i])),("Z",unpack(z[i])),("Y",unpack(x[i]^z[i]))] for i in range(inst.n)]
    return result[::-1],incoming


def effects(groups,weight):
    for sites in combinations(range(len(groups)),weight):
        for modes in product(*(groups[i] for i in sites)):
            sig=tuple(reduce(xor,terms,0) for terms in zip(*(m[1] for m in modes))) if modes else (0,0,0)
            yield [[i,m[0]] for i,m in zip(sites,modes)],sig


def residual_profile(inst,r,s,cap):
    """Joint enumeration by input-qubit weight and faulty-location weight."""
    if any(type(i)!=int or i<0 for i in (r,s,cap)):
        raise ValueError("invalid incoming, internal or enumeration budget")
    faults,inputs=columns(inst)
    if r>len(inputs) or s>len(faults):
        raise ValueError("weight exceeds number of sites")
    tested=0
    profile=[]
    for iw in range(r+1):
        for fw in range(s+1):
            accepted=wrong=max_weight=0
            worst=None
            for incoming,a in effects(inputs,iw):
                for internal,b in effects(faults,fw):
                    if tested>=cap:
                        return dict(status="inconclusive",tested=tested,completed=profile)
                    tested+=1
                    d,bit,res=(u^v for u,v in zip(a,b))
                    if d:
                        continue
                    accepted+=1
                    wrong+=bool(bit)
                    support=(res|(res>>inst.n))&((1<<inst.n)-1)
                    if support.bit_count()>max_weight:
                        max_weight=support.bit_count()
                        worst=dict(incoming=incoming,internal=internal,residual=res)
            profile.append(dict(input_weight=iw,location_weight=fw,accepted=accepted,
                                wrong_public=wrong,residual_weight_upper=max_weight,witness=worst))
    return dict(status="proved",tested=tested,completed=profile)

