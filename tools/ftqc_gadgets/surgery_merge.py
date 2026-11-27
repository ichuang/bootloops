"""Compile a CSS cone to repeated merged-code syndrome extraction, then split once."""
from ftqc_codes.gf2 import supp
from .surgery_build import combine,independent_section
from .surgery_geometry import restricted_cone


def word(n,mask,axis):
    return "+" + "".join(axis if mask>>i & 1 else "I" for i in range(n))


def merge(n,hx,hz,ports,aux,rounds,axis="Z",faces=None,phenomenological=False):
    if type(rounds)!=int or rounds<1 or axis not in ("X","Z"):
        raise ValueError("invalid merged syndrome rounds or axis")
    dual=axis=="X"
    a,b=(hz,hx) if dual else (hx,hz)
    cone=restricted_cone(n,a,b,ports,aux,faces)
    ticks=[]
    def gate(op,*args):
        name="!"+op if phenomenological and op not in ("M","MX") else op
        ticks.append([[name,*args]])
    def basis(ordinary,conjugate):
        return conjugate if dual else ordinary
    def cx(c,t):
        gate("CX",t,c) if dual else gate("CX",c,t)
    e=len(aux)
    meters=cone["n"]
    for i in range(e):
        gate(basis("RX","R"),n+i)
    nr=0
    detectors=[]
    prior=None
    for _ in range(rounds):
        if phenomenological:
            ticks.append([["I",q] for q in range(n+e)])
        currents=[]
        for j,(is_x,row) in enumerate([(False,r) for r in cone["hz"]]+[(True,r) for r in cone["hx"]]):
            meter=meters+j
            gate(basis("RX","R") if is_x else basis("R","RX"),meter)
            for q in supp(row):
                if is_x:
                    cx(meter,q)
                else:
                    cx(q,meter)
            gate(basis("MX","M") if is_x else basis("M","MX"),meter)
            currents.append(nr)
            nr+=1
        if prior is None:
            first=currents
            detectors.extend([r] for j,r in enumerate(currents) if j<len(b) or j>=len(cone["hz"]))
            for z,u in cone["relations"]:
                detectors.append([first[i] for i in supp(z)]+[first[len(b)+i] for i in supp(u)])
        else:
            detectors.extend([a,b] for a,b in zip(prior,currents))
        prior=currents
    start=nr
    for i in range(e):
        gate(basis("MX","M"),n+i)
    # The last measured cycle and the final split parity must agree.
    for i,f in enumerate(cone["faces"]):
        detectors.append([prior[len(cone["hz"])+len(a)+i]]+[start+j for j in supp(f)])
    section=independent_section(aux,len(ports))
    for i in range(e):
        data=combine(ports,sum((r>>i&1)<<j for j,r in enumerate(section)))
        for q in supp(data):
            if dual:
                gate("H",q)
            gate("Z_IF",q,start+i)
            if dual:
                gate("H",q)
    checks=[word(n,r,"X") for r in hx]+[word(n,r,"Z") for r in hz]
    instrument=dict(n=n,checks=checks,targets=[word(n,r,axis) for r in cone["measured"]],
                    ticks=ticks,readouts=[dict(records=[first[len(b)+j] for j in supp(u)],sign=0) for u in cone["selectors"]],
                    detectors=detectors,terminal_syndrome=True)
    return cone,instrument

