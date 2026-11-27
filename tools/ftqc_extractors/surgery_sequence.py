"""Finite mapping-cone and padded-sequence hypotheses of arXiv:2510.14895."""
from ftqc_codes.gf2 import Echelon
from .algebra import bits, css_dimension, distance_below, kernel, transpose
from .block_reading import combine, memory, quotient, rows


def maps(columns, domain, codomain):
    if len(columns)!=domain:
        raise ValueError("chain map domain mismatch")
    return rows(columns,codomain)


def sparse(matrices, limit):
    if type(limit) is not int or limit<1:
        raise ValueError("a positive finite locality bound is required")
    for matrix,width in matrices:
        if any(r.bit_count()>limit for r in matrix+transpose(matrix,width)):
            raise ValueError("sparsity bound exceeded")


def append_cone(code, op, offset, x, z, locality):
    """Append a cone summand. Maps are columns in original memory coordinates."""
    if set(op)!={"aux","cycles","f0","f1","f2","metachecks"}:
        raise ValueError("missing or unknown surgery fields")
    n,hx,hz=code
    v,ax,az=memory(op["aux"])
    f0=maps(op["f0"],len(ax),len(hx))
    f1=maps(op["f1"],v,n)
    f2=maps(op["f2"],len(az),len(hz))
    cycles=rows(op["cycles"],len(ax))
    at=transpose(ax,v)
    xt=transpose(hx,n)
    if any(r.bit_count()>1 for r in transpose(f1,n)):
        raise ValueError("overlapping vertex map")
    if any(combine(xt,col)!=combine(f0,a) for col,a in zip(f1,at)):
        raise ValueError("chain square at degree one")
    if any(combine(f1,a)!=combine(hz,col) for a,col in zip(az,f2)):
        raise ValueError("chain square at degree two")
    if (len(cycles)!=len(ax)-Echelon(ax).rank or Echelon(cycles).rank!=len(cycles)
            or any(combine(ax,c) for c in cycles)):
        raise ValueError("cycle checks must form a complete independent basis")
    sparse([(f0,len(hx)),(f1,n),(f2,len(hz)),(ax,v),(az,v),(cycles,len(ax))],locality)
    for i,extra in enumerate(transpose(f0,len(hx))):
        x[i] |= extra<<offset
    x.extend(c<<offset for c in cycles)
    first=len(z)
    z.extend(col|(edge<<offset) for col,edge in zip(f1,at))
    # The local schedule numbers old checks first, followed by its own vertices.
    meta=[col|(a<<len(hz)) for col,a in zip(f2,az)]
    if op["metachecks"]!=meta or any(combine(z[:len(hz)]+z[first:],m) for m in meta):
        raise ValueError("metachecks missing or inconsistent with the cone")
    logicals=[combine(f1,vv) for vv in quotient(ax,az,v)]
    return offset+len(ax),logicals


def padding(s):
    p=s["padding"]
    if not isinstance(p,list) or len(p)!=2 or any(type(x) is not int or x<s["target"] for x in p):
        raise ValueError("sequence requires target syndrome rounds at both boundaries")


def analyze(s,budget):
    code=memory(s["memory"])
    n,hx,hz=code
    padding(s)
    ops=s["operations"]
    if not isinstance(ops,list) or not ops:
        raise ValueError("empty surgery sequence")
    d=s["target"]
    sparse([(hx,n),(hz,n)],s["max_degree"])
    bad=distance_below(code,d,budget)
    if bad:
        return dict(status="counterexample",obligation="base code",witness=bad)
    x,z=list(hx),list(hz)
    width=n
    measured=[]
    logical_basis=Echelon(hz)
    for i,op in enumerate(ops):
        av,ax,az=memory(op["aux"])
        width,images=append_cone(code,op,width,x,z,s["max_degree"])
        # Check the full auxiliary CSS distance (sufficient for the required X cosystole).
        bad=distance_below((av,ax,az),d,budget)
        if bad:
            return dict(status="counterexample",obligation="auxiliary CSS distance",index=i,witness=bad)
        for v in images:
            if not logical_basis.add(v):
                raise ValueError("dependent or trivial induced logical measurement")
        measured.append(images)
    if css_dimension((width,x,z))!=css_dimension(code)-sum(map(len,measured)):
        raise ValueError("unexpected compacted logical space")
    bad=distance_below((width,x,z),d,budget)
    if bad:
        return dict(status="counterexample",obligation="compacted code distance",witness=bad)
    return dict(status="proved",theorem="2510.14895 Theorem 5.3",model="phenomenological Z surgery",
                distance_lower=d,padding=s["padding"],merge_rounds=[1]*len(ops),measured=measured,
                compacted=dict(n=width,hx=x,hz=z))


def full_operation(code,blocks,pattern):
    """Column representation of one transversal pattern for the general cone verifier."""
    from .block_reading import copies
    n,hx,hz=code
    def replicate(width):
        return [copies(1<<i,pattern,width) for i in range(width)]
    f2=replicate(len(hz))
    return dict(aux=dict(n=n,hx=hx,hz=hz),cycles=kernel(transpose(hx,n),len(hx)),
                f0=replicate(len(hx)),f1=replicate(n),f2=f2,
                metachecks=[col|(h<<(blocks*len(hz))) for col,h in zip(f2,hz)])

