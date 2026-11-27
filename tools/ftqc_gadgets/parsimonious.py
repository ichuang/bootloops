"""Finite cellulated cones for graph surgery, after Yuan et al., Sec. IV.

Represent compressed prefix trees by sets of incidence leaves. Between two
bit orders related by disjoint adjacent swaps, collapse the changed splits
in both trees to their common hierarchy. Glue the two mapping cylinders and
cellulate original edges along their routed paths. Verify the resulting
binary boundaries and homology independently.
"""
from ftqc_codes.gf2 import Echelon
from .surgery_build import combine


def hierarchy(keys,order):
    layers={frozenset(range(len(keys)))}
    groups=[list(range(len(keys)))]
    for bit in order:
        next_groups=[]
        for group in groups:
            for digit in (0,1):
                subset=[j for j in group if keys[j]>>bit & 1 == digit]
                if subset:
                    next_groups.append(subset)
                    layers.add(frozenset(subset))
        groups=next_groups
    return layers


def parents(groups):
    root=max(groups,key=len)
    return {c:min((s for s in groups if c<s),key=len) for c in groups if c!=root}


def cellulate(n,graph):
    if n<2 or not graph:
        raise ValueError("a connected graph of at least two vertices is required")
    base=[(1<<a)|(1<<b) for a,b in graph if type(a)==int and type(b)==int and 0<=a<n and 0<=b<n and a!=b]
    if len(base)!=len(graph) or len(set(base))!=len(base) or Echelon(base).rank!=n-1:
        raise ValueError("expected a simple connected graph")
    # Subdivide once: left labels are vertices and right labels are edge IDs.
    incidences=[(v,e) for e,pair in enumerate(graph) for v in pair]
    left=max(1,(n-1).bit_length())
    right=max(1,(len(graph)-1).bit_length())
    keys=[(e<<left)|v for v,e in incidences]
    order=list(range(left+right))
    orders=[order.copy()]
    for parity in range(left+right):
        changed=False
        for i in range(parity%2,len(order)-1,2):
            if order[i]<left<=order[i+1]:
                order[i],order[i+1]=order[i+1],order[i]
                changed=True
        if changed:
            orders.append(order.copy())
    if order[:right]!=list(range(left,left+right)):
        raise ValueError("bit-order transposition failed")
    trees=[hierarchy(keys,o) for o in orders]
    nodes=0
    edges=[]
    faces=[]
    def edge(a,b):
        if a==b:
            raise ValueError("degenerate cone edge")
        edges.append((a,b))
        return len(edges)-1
    def tree(groups):
        nonlocal nodes
        ids={c:nodes+i for i,c in enumerate(sorted(groups,key=lambda s:(len(s),tuple(sorted(s)))))}
        nodes+=len(ids)
        parent=parents(groups)
        up={c:edge(ids[c],ids[p]) for c,p in parent.items()}
        return ids,parent,up
    levels=[tree(t) for t in trees]
    joins=[]
    for i in range(len(trees)-1):
        shared=trees[i]&trees[i+1]
        ids,_,up=tree(shared)
        pair=[]
        for t,level in ((trees[i],levels[i]),(trees[i+1],levels[i+1])):
            local,parent,local_up=level
            collapse={c:min((s for s in shared if c<=s),key=len) for c in t}
            vertical={c:edge(local[c],ids[collapse[c]]) for c in t}
            for c,p in parent.items():
                side=[local_up[c],vertical[c],vertical[p]]
                if collapse[c]!=collapse[p]:
                    side.append(up[collapse[c]])
                faces.append(side)
            pair.append(vertical)
        joins.append(pair)
    ports=[]
    paths=[]
    for v in range(n):
        c=frozenset(j for j,(a,_) in enumerate(incidences) if a==v)
        ports.append(levels[0][0][c])
    for j,(v,e) in enumerate(incidences):
        leaf=frozenset([j])
        trail=[]
        c=leaf
        while levels[0][0][c]!=ports[v]:
            trail.append(levels[0][2][c])
            c=levels[0][1][c]
        trail.reverse()
        for a,b in joins:
            trail.extend((a[leaf],b[leaf]))
        terminal=frozenset((2*e,2*e+1))
        c=leaf
        while c!=terminal:
            trail.append(levels[-1][2][c])
            c=levels[-1][1][c]
        paths.append(trail)
    embedded=[]
    for i,(a,b) in enumerate(graph):
        route=paths[2*i]+paths[2*i+1][::-1]
        # The two paths can share a suffix on the final tree.
        trail=[]
        for eid in route:
            if trail and trail[-1]==eid:
                trail.pop()
            else:
                trail.append(eid)
        vertices=[ports[a]]
        for eid in trail:
            u,v=edges[eid]
            if vertices[-1] not in (u,v):
                raise ValueError("broken deformation path")
            vertices.append(u^v^vertices[-1])
        if vertices[-1]!=ports[b] or len(set(vertices))!=len(vertices):
            raise ValueError("invalid deformation endpoints")
        lo,hi=0,len(trail)
        chord=edge(vertices[lo],vertices[hi])
        embedded.append(chord)
        while hi-lo>2:
            nxt=edge(vertices[lo+1],vertices[hi])
            faces.append([chord,trail[lo],nxt])
            chord=nxt
            lo+=1
            if hi-lo>2:
                nxt=edge(vertices[lo],vertices[hi-1])
                faces.append([chord,trail[hi-1],nxt])
                chord=nxt
                hi-=1
        faces.append([chord,*trail[lo:hi]])
    rows=[(1<<a)|(1<<b) for a,b in edges]
    cycles=[sum(1<<e for e in f) for f in faces]
    rank=Echelon(rows).rank
    if rank!=nodes-1 or any(combine(rows,f) for f in cycles) or Echelon(cycles).rank!=len(cycles) or len(cycles)!=len(edges)-rank:
        raise ValueError("cellulation failed exact homology")
    return dict(vertices=nodes,edges=edges,faces=cycles,ports=ports,original_edges=embedded,
                max_face=max(map(len,faces)),layers=len(levels))

