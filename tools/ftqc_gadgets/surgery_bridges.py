"""Finite graph enhancement and matching bridges for CSS measurements."""
from itertools import combinations
from ftqc_codes.gf2 import Echelon,nullspace
from .surgery_geometry import transpose


def enhanced(sites,opposite,budget):
    """Pair restricted Tanner incidences, then add edges across deficient cuts."""
    n=len(sites)
    if not n:
        raise ValueError("empty logical support")
    edges=set()
    for check in opposite:
        hits=[j for j,q in enumerate(sites) if check>>q&1]
        if len(hits)%2:
            raise ValueError("support anticommutes with a check")
        edges.update(zip(hits[::2],hits[1::2]))
    # A path ensures a connected induced apparatus even on reducible supports.
    edges.update((i,i+1) for i in range(n-1))
    spent=0
    while True:
        cut_found=False
        for mask in range(1,1<<(n-1)):
            spent+=1
            if spent>budget:
                return dict(edges=sorted(edges),expansion_status="inconclusive",tested=budget)
            weight=min(mask.bit_count(),n-mask.bit_count())
            crossing=[p for p in combinations(range(n),2) if ((mask>>p[0])^(mask>>p[1]))&1]
            if sum(p in edges for p in crossing)<weight:
                edges.add(next(p for p in crossing if p not in edges))
                cut_found=True
                break
        if not cut_found:
            return dict(edges=sorted(edges),expansion_status="unit lower bound",tested=spent)


def bridge(left_size,right_size,left,right,links):
    """Two induced graphs joined by bridge qubits. Local cycles plus B-1 bridge relations."""
    if min(left_size,right_size)<1 or not links or len(set(a for a,b in links))!=len(links) or len(set(b for a,b in links))!=len(links):
        raise ValueError("bridge links must be a nonempty matching")
    local=[]
    for edges,n,offset in ((left,left_size,0),(right,right_size,left_size)):
        if any(not (0<=a<n and 0<=b<n and a!=b) for a,b in edges):
            raise ValueError("invalid component edge")
        rows=[(1<<a)|(1<<b) for a,b in edges]
        if len(set(rows))!=len(rows) or Echelon(rows).rank!=n-1:
            raise ValueError("bridge endpoints require connected simple components")
        local.extend([(offset+a,offset+b) for a,b in edges])
    if any(not (0<=a<left_size and 0<=b<right_size) for a,b in links):
        raise ValueError("invalid bridge endpoint")
    edges=local+[(a,left_size+b) for a,b in links]
    incidence=[(1<<a)|(1<<b) for a,b in edges]
    faces=nullspace(transpose(incidence,left_size+right_size),len(edges))
    projected=[f>>len(local) for f in faces]
    if Echelon(projected).rank!=len(links)-1 or any(v.bit_count()%2 for v in projected):
        raise ValueError("bridge restrictions fail to generate the even subspace")
    return dict(edges=edges,faces=faces,bridge_restrictions=projected)

