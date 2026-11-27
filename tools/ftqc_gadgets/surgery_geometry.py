"""Finite CSS mapping cones and the SkipTree repetition adapter."""
from collections import deque
from itertools import combinations
from fractions import Fraction
from ftqc_codes.gf2 import Echelon, dot, nullspace
from .surgery_build import combine, independent_section


def transpose(rows,width):
    return [sum(((a>>j)&1)<<i for i,a in enumerate(rows)) for j in range(width)]


def check_rows(rows,width):
    if type(width)!=int or width<0 or any(type(r)!=int or not 0<=r<1<<width for r in rows):
        raise ValueError("matrix entry outside declared width")


def restricted_cone(n,hx,hz,ports,aux,faces=None):
    """Build actual deformed checks; solve the chain square for each data X check."""
    v,e=len(ports),len(aux)
    for rows,width in [(hx,n),(hz,n),(ports,n),(aux,v)]:
        check_rows(rows,width)
    if any(dot(x,z) for x in hx for z in hz):
        raise ValueError("CSS checks anticommute")
    incidence=transpose(aux,v)
    if faces is None:
        faces=nullspace(incidence,e)
    check_rows(faces,e)
    if Echelon(faces).rank!=e-Echelon(aux).rank or any(combine(aux,f) for f in faces):
        raise ValueError("incomplete or nonclosed auxiliary faces")
    section=independent_section(aux,v)
    lifts=[]
    for x in hx:
        overlaps=sum(dot(x,p)<<j for j,p in enumerate(ports))
        row=combine(section,overlaps)
        if combine(aux,row)!=overlaps:
            raise ValueError("no dressing for memory check")
        lifts.append(row)
    mx=[x|(a<<n) for x,a in zip(hx,lifts)]+[f<<n for f in faces]
    mz=hz+[p|(i<<n) for p,i in zip(ports,incidence)]
    code=Echelon(hz)
    selectors=[]
    dependencies=[]
    for u in nullspace(aux,v):
        if code.add(combine(ports,u)):
            selectors.append(u)
    relations=nullspace(transpose(hz+[combine(ports,u) for u in nullspace(aux,v)],n),
                        len(hz)+len(nullspace(aux,v)))
    kernel=nullspace(aux,v)
    for r in relations:
        u=combine(kernel,r>>len(hz))
        if u:
            dependencies.append([r&((1<<len(hz))-1),u])
    before=n-Echelon(hx).rank-Echelon(hz).rank
    after=n+e-Echelon(mx).rank-Echelon(mz).rank
    if before-after!=len(selectors) or any(dot(x,z) for x in mx for z in mz):
        raise ValueError("invalid merged code")
    return dict(n=n+e,hx=mx,hz=mz,logical_qubits=after,selectors=selectors,
                measured=[combine(ports,u) for u in selectors],relations=dependencies,
                dressing=lifts,faces=faces,
                max_check=max(map(int.bit_count,mx+mz),default=0),
                max_qubit_degree=max(map(int.bit_count,transpose(mx+mz,n+e)),default=0))


def distance(n,hx,hz,target,budget):
    """Static CSS distance interval for the actual finite deformed code."""
    spaces=[(hz,Echelon(hx),"X"),(hx,Echelon(hz),"Z")]
    tested=0
    for w in range(1,min(n,target-1)+1):
        for sites in combinations(range(n),w):
            mask=sum(1<<i for i in sites)
            for opposite,stabs,axis in spaces:
                if tested>=budget:
                    return dict(status="inconclusive",lower=w,tested=tested)
                tested+=1
                if mask not in stabs and all(not dot(mask,r) for r in opposite):
                    return dict(status="counterexample",lower=w,upper=w,axis=axis,support=mask,tested=tested)
    return dict(status="proved",lower=target,tested=tested)


def expansion(ports,aux,tau,budget):
    """Definition 2.14 at degree one, including arbitrary proof coordinates."""
    v=len(ports)
    check_rows(aux,v)
    if type(tau)!=int or tau<1 or any(type(p)!=int or p<0 for p in ports):
        raise ValueError("invalid relative map or target")
    image=Echelon(combine(ports,u) for u in nullspace(aux,v))
    rank=image.rank
    if 1<<(v+rank)>budget:
        return dict(status="inconclusive",reason="relative expansion enumeration")
    codewords=[combine(image.basis,s) for s in range(1<<rank)]
    best=None
    witness=None
    for x in range(1<<v):
        fx=combine(ports,x)
        den=min(tau,min((fx^c).bit_count() for c in codewords))
        if den:
            ratio=Fraction(sum(dot(x,h) for h in aux),den)
            if best is None or ratio<best:
                best,witness=ratio,x
    return dict(status="proved",ratio=str(best) if best is not None else None,
                minimizer=witness,unit_bound=best is None or best>=1)


def skip_tree(size,edges):
    """Alg. 1 of Swaroop et al.; return cyclic labels and original edge paths."""
    adj=[[] for _ in range(size)]
    for i,(a,b) in enumerate(edges):
        if not 0<=a<size or not 0<=b<size or a==b:
            raise ValueError("invalid edge")
        adj[a].append((b,i))
        adj[b].append((a,i))
    children=[[] for _ in adj]
    path={0:0}
    queue=deque([0])
    while queue:
        a=queue.popleft()
        for b,i in adj[a]:
            if b not in path:
                children[a].append(b)
                path[b]=path[a]^(1<<i)
                queue.append(b)
    if len(path)!=size:
        raise ValueError("disconnected adapter ports")
    order=[]
    stack=[(0,True,False)]
    while stack:
        a,first,emit=stack.pop()
        if emit:
            order.append(a)
            continue
        if first:
            order.append(a)
        else:
            stack.append((a,first,True))
        stack.extend((b,not first,False) for b in reversed(children[a]))
    paths=[path[a]^path[b] for a,b in zip(order,order[1:]+order[:1])]
    incidence=[(1<<a)|(1<<b) for a,b in edges]
    for i,t in enumerate(paths):
        if combine(incidence,t)!=(1<<order[i])^(1<<order[(i+1)%size]) or t.bit_count()>3:
            raise ValueError("SkipTree path identity failed")
    if any(r.bit_count()>2 for r in transpose(paths,len(edges))):
        raise ValueError("SkipTree congestion bound failed")
    return order,paths


def adapter(size,left,right):
    """Eq. (19) on two equal-size connected port sets; measure one joint parity."""
    if type(size)!=int or size<1:
        raise ValueError("invalid port count")
    l,lp=skip_tree(size,left)
    r,rp=skip_tree(size,right)
    edges=left+[(a,size+b) for a,b in zip(l,r)]+[(size+a,size+b) for a,b in right]
    rows=[(1<<a)|(1<<b) for a,b in edges]
    start=len(left)
    faces=nullspace(transpose(rows[:start],size*2),start)
    rightstart=start+size
    faces.extend(c<<rightstart for c in nullspace(transpose(rows[rightstart:],size*2),len(right)))
    for i in range(size-1):
        faces.append(lp[i]|(3<<(start+i))|(rp[i]<<rightstart))
    if any(combine(rows,f) for f in faces) or Echelon(faces).rank!=len(edges)-2*size+1:
        raise ValueError("invalid repetition-adapter cycles")
    return dict(edges=edges,faces=faces,matching=list(zip(l,r)),left_paths=lp,right_paths=rp)

