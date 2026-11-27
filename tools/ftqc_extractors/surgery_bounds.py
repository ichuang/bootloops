"""Finite modular expansion, Cartesian thickening, and intermediate-time surgery."""
from collections import deque
from fractions import Fraction
from functools import reduce
from operator import or_

from .algebra import css_dimension, distance_below, parity, span, transpose, Limit
from .block_reading import memory, rows
from .surgery_sequence import append_cone, padding, sparse, cosystole_witness


def expansion(s,budget):
    n,layers=s["vertices"],s["layers"]
    if type(n) is not int or n<1 or type(layers) is not int or layers<1:
        raise ValueError("invalid hypergraph dimensions")
    edges=rows(s["edges"],n)
    generators=rows(s["components"],n)
    if not generators or any(not e for e in edges) or any(parity(e,w) for e in edges for w in generators):
        raise ValueError("components must be orthogonal to incidence")
    W=span(generators,budget)
    U=reduce(or_,W,0)
    if not U:
        raise ValueError("empty distinguished port set")
    def minimum(width,inc):
        if width>=budget.bit_length() or (1<<width)*len(W)>budget:
            raise Limit("modular expansion enumeration")
        best=None
        for v in range(1<<width):
            den=min(s["target"],min(((v&U)^w).bit_count() for w in W))
            if den:
                ratio=Fraction(sum(parity(v,e) for e in inc),den)
                if best is None or ratio<best:
                    best,shore=ratio,v
        if best is None:
            raise ValueError("expansion denominator vanishes for all cuts")
        return best,shore
    a,shore=minimum(n,edges)
    if n*layers>budget.bit_length():
        raise Limit("thickened hypergraph enumeration")
    thick=[e<<(i*n) for i in range(layers) for e in edges]
    thick += [(1<<i)|(1<<(n+i)) for i in range((layers-1)*n)]
    b,other=minimum(n*layers,thick)
    return dict(status="proved",model="exact finite modular expansion, ports in layer zero",
                original=str(a),original_cut=shore,thickened=str(b),thickened_cut=other,
                vertices=n*layers,hyperedges=thick,meets_unit_bound=b>=1)


def soundness(code,budget):
    """Minimum over syndromes of n |s| / (m * coset-leader weight)."""
    n,hx,hz=code
    m=len(hx)+len(hz)
    x=transpose(hz,n)
    z=[v<<len(hz) for v in transpose(hx,n)]
    flips=set(x+z+[a^b for a,b in zip(x,z)])-{0}
    dist={0:0}
    queue=deque([0])
    while queue:
        s=queue.popleft()
        budget-=len(flips)
        if budget<0:
            raise Limit("quantum soundness syndrome budget")
        for f in flips:
            if s^f not in dist:
                dist[s^f]=dist[s]+1
                queue.append(s^f)
    if not flips:
        raise ValueError("soundness requires nontrivial checks")
    return min(Fraction(n*s.bit_count(),m*w) for s,w in dist.items() if w)


def intermediate(s,budget):
    padding(s)
    code=memory(s["memory"])
    n,hx,hz=code
    d=s["target"]
    sparse([(hx,n),(hz,n)],s["max_degree"])
    rho=soundness(code,budget)
    if rho<Fraction(2*n,len(hx)+len(hz)):
        return dict(status="counterexample",obligation="Theorem 6.5 soundness",rho=str(rho))
    bad=distance_below(code,d,budget)
    if bad:
        return dict(status="counterexample",obligation="memory distance",witness=bad)
    rounds=s["rounds"]
    if not rounds or len(rounds)!=len(s["operations"]):
        raise ValueError("one positive round count per surgery")
    measured=[]
    for op,r in zip(s["operations"],rounds):
        if type(r) is not int or r<1:
            raise ValueError("invalid round count")
        ax=memory(op["aux"])
        bad=cosystole_witness(ax,(d+r-1)//r,budget)
        if bad:
            return dict(status="counterexample",obligation="auxiliary distance times rounds",witness=bad)
        x,z=list(hx),list(hz)
        width,images=append_cone(code,op,n,x,z,s["max_degree"])
        merged=(width,x,z)
        if css_dimension(merged)!=css_dimension(code)-len(images):
            raise ValueError("unexpected deformed logical dimension")
        bad=distance_below(merged,d,budget)
        if bad:
            return dict(status="counterexample",obligation="deformed code distance",witness=bad)
        measured.append(images)
    return dict(status="proved",theorem="2510.14895 Theorem 6.5",model="phenomenological",
                distance_lower=d,padding=s["padding"],merge_rounds=rounds,measured=measured,
                quantum_soundness=str(rho))

