"""Finite tensor and cyclic lifts for fast surgery protocols."""
from itertools import combinations
from ftqc_codes.gf2 import Echelon,dot,nullspace,supp
from ftqc_gadgets.surgery_geometry import transpose,check_rows,restricted_cone
from ftqc_gadgets.surgery_build import combine
from .algebra import hgp


def product_sequence(c,nc,d,nd,words,target,padding,degree):
    """Chang et al. Eq. (42), tensor the graph-to-classical-code map with id_D."""
    check_rows(c,nc)
    check_rows(d,nd)
    mc,md=len(c),len(d)
    n,hx,hz=hgp(c,nc,transpose(d,nd),md)
    ops=[]
    for word in words:
        if type(word)!=int or not 0<word<1<<nc or any(dot(r,word) for r in c):
            raise ValueError("logical row is outside the classical kernel")
        sites=supp(word)
        edges=list(combinations(range(len(sites)),2))
        incidence=[(1<<a)|(1<<b) for a,b in edges]
        # Each original check restricts to an even subset. Pair its terminals.
        dressing=[0]*len(edges)
        index={p:i for i,p in enumerate(edges)}
        for i,row in enumerate(c):
            hits=[j for j,q in enumerate(sites) if row>>q&1]
            for pair in zip(hits[::2],hits[1::2]):
                dressing[index[pair]] |= 1<<i
        cycles=nullspace(transpose(incidence,len(sites)),len(edges))
        av,ax,az=hgp(incidence,len(sites),transpose(d,nd),md)
        ax += [sum(1<<(len(sites)*md+e*nd+j) for e in supp(cy)) for cy in cycles for j in range(nd)]
        f1=[1<<(q*md+j) for q in sites for j in range(md)]
        f1 += [sum(1<<(nc*md+i*nd+j) for i in supp(row)) for row in dressing for j in range(nd)]
        f0=[sum(1<<(i*md+j) for i in supp(row)) for row in dressing for j in range(md)]
        f0 += [0]*(len(cycles)*nd)
        f2=[1<<(q*nd+j) for q in sites for j in range(nd)]
        ops.append(dict(aux=dict(n=av,hx=ax,hz=az),f0=f0,f1=f1,f2=f2,
                        cycles=nullspace(transpose(ax,av),len(ax)),
                        metachecks=[col|(a<<len(hz)) for col,a in zip(f2,az)]))
    return dict(schema=1,kind="surgery_sequence",memory=dict(n=n,hx=hx,hz=hz),operations=ops,
                target=target,padding=padding,max_degree=degree)


def lifted_product(a,b,lift):
    """C1 tensor D0 + C0 tensor D1, each ring entry is a cyclic-polynomial mask."""
    if type(lift)!=int or lift<1 or not a or not b or not a[0] or not b[0]:
        raise ValueError("empty lifted product")
    m,n,r,k=len(a),len(a[0]),len(b),len(b[0])
    for matrix,cols in ((a,n),(b,k)):
        if any(len(row)!=cols for row in matrix):
            raise ValueError("ragged protograph")
        for row in matrix:
            check_rows(row,lift)
    width=(n*r+m*k)*lift
    def orbit(cell,shift,coeff,reverse=False):
        return sum(1<<(cell*lift+(shift+(-e if reverse else e))%lift) for e in supp(coeff))
    hx=[sum(orbit(u*r+j,t,a[i][u],True) for u in range(n)) |
        sum(orbit(n*r+i*k+v,t,b[j][v],True) for v in range(k))
        for i in range(m) for j in range(r) for t in range(lift)]
    hz=[sum(orbit(u*r+j,t,b[j][v]) for j in range(r)) |
        sum(orbit(n*r+i*k+v,t,a[i][u]) for i in range(m))
        for u in range(n) for v in range(k) for t in range(lift)]
    if any(dot(x,z) for x in hx for z in hz):
        raise ValueError("lifted square fails")
    return width,hx,hz


def radial90():
    return lifted_product([[4,8,16],[2,16,2],[16,8,4]],
                          [[4,4,1],[16,1,16],[2,8,1]],5)


def shift_word(word,n,lift,amount=1):
    if n%lift or word<0 or word>>n:
        raise ValueError("invalid cyclic orbit dimensions")
    return sum(1<<(q//lift*lift+(q%lift+amount)%lift) for q in supp(word))


def invariant_logicals(checks,opposite,n,lift):
    """Orbit norm over an odd lift; independent invariant homology classes."""
    basis=Echelon(checks)
    selected=[]
    for v in nullspace(opposite,n):
        orbit=0
        for j in range(lift):
            orbit ^= shift_word(v,n,lift,j)
        if basis.add(orbit):
            selected.append(orbit)
    return selected


def radial_design():
    """A rank-one projected chain map in the completely fine radial block.

    f1 = v w^T, f0=f2=0 with invariant Z and X representatives. This is an
    explicit dense representative of the homotopy class, before optimization.
    """
    n,hx,hz=radial90()
    zs=invariant_logicals(hz,hx,n,5)
    xs=invariant_logicals(hx,hz,n,5)
    v=min(zs,key=int.bit_count)
    w=min((u for u in xs if dot(v,u)),key=int.bit_count)
    ports=[v if w>>q&1 else 0 for q in range(n)]
    c=restricted_cone(n,hx,hz,ports,hx)
    if c["logical_qubits"]!=7 or len(c["measured"])!=1:
        raise ValueError("wrong radial logical rank")
    from ftqc_gadgets.surgery_merge import word
    return dict(n=n,hx=hx,hz=hz,method="linear",ports=ports,aux=hx,axis="Z",
                targets=[word(n,v,"Z")],rounds=1,distance_target=2)

