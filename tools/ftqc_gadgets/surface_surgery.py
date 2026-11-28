"""Nearest-neighbor merge and split of two mirrored rotated planar patches.

A shared boundary meter measures the two half checks in successive subrounds
when split, and their product in one subround when merged. The complementary
seam checks multiply to the joint logical Pauli. Qubit zero is the low bit.
"""
from functools import reduce
from operator import xor

from ftqc_codes.codes import rotated_plaquettes
from ftqc_codes.gf2 import dot, from_supp, supp, systematic
from .surgery_merge import word


def patch_pair(d, rounds, axis):
    if type(d) is not int or d < 3 or d % 2 == 0:
        raise ValueError("expected an odd patch distance >= 3")
    if axis not in ("X", "Z") or type(rounds) is not int or rounds < 1:
        raise ValueError("expected X or Z and positive merge and split rounds")
    normal = 0 if axis == "Z" else 1
    n = 2*d*d
    points = []  # doubled row and column; meters are at odd coordinates
    def reflect(point):
        reflected = list(point)
        reflected[normal] = 4*d-2 - point[normal]
        return tuple(reflected)

    data = [(2*r, 2*c) for r in range(d) for c in range(d)]
    points.extend(data + [reflect(p) for p in data])
    ids = {point: q for q, point in enumerate(points)}
    halves = [{}, {}]
    split = {}
    merged = {}

    def check(kind, center, sites):
        key = (kind, from_supp(sites))
        if center not in ids:
            ids[center] = len(points)
            points.append(center)
        return key, ids[center]

    for block in (0, 1):
        for kind, (i, j), cells in rotated_plaquettes(d):
            center = (2*i-1, 2*j-1)
            sites = [q + block*d*d for q in cells.values()]
            key, meter = check(kind, reflect(center) if block else center, sites)
            split[key] = meter
            if center[normal] == 2*d-1:
                halves[block][meter] = key
            else:
                merged[key] = meter

    # On alternating seam plaquettes, the old half checks become one full
    # check. All other seam checks measure new gauge degrees of freedom.
    seam = []
    for j in range(d+1):
        center = (2*d-1, 2*j-1)
        if normal == 1:
            center = center[::-1]
        sites = [q for q in range(n) if
                 all(abs(a-b) == 1 for a,b in zip(points[q], center))]
        meter = ids.get(center)
        if meter in halves[0]:
            kind, mask = halves[0][meter]
            key = kind, mask | halves[1][meter][1]
        else:
            if ("Z" if (d+j) % 2 else "X") != axis:
                continue  # unused corner at the outer boundary
            key, meter = check(axis, center, sites)
            seam.append(key)
        merged[key] = meter

    ticks, detectors = [], []
    record_count = 0

    def emit(ops, live):
        busy = {q for op in ops for q in (op[1:2] if op[0] == "Z_IF" else op[1:])}
        ticks.append(ops + [["I", q] for q in sorted(live - busy)])

    def subround(checks):
        nonlocal record_count
        live = set(range(n)) | set(checks.values())
        emit([["RX" if k == "X" else "R", m] for (k,_),m in checks.items()], live)
        orders = {"X": [(-1,-1),(-1,1),(1,-1),(1,1)],
                  "Z": [(-1,-1),(1,-1),(-1,1),(1,1)]}
        for step in range(4):
            ops = []
            for (kind, mask), meter in checks.items():
                i,j = points[meter]
                di,dj = orders[kind][step]
                q = ids.get((i+di,j+dj))
                if q is not None and mask >> q & 1:
                    ops.append(["CX", meter, q] if kind == "X" else ["CX", q, meter])
            emit(ops, live)
        emit([["MX" if k == "X" else "M", m] for (k,_),m in checks.items()], live)
        records = {k: record_count+i for i,k in enumerate(checks)}
        record_count += len(records)
        return records

    def split_round():
        delayed = set(halves[1].values())
        first = subround({k:m for k,m in split.items() if k not in delayed})
        return first | subround({k:m for k,m in split.items() if k in delayed})

    before = split_round()
    detectors.extend([[r] for r in before.values()])
    previous = subround(merged)
    detectors.extend([[before[k], previous[k]] for k in sorted(before.keys() & previous.keys())])
    for m, (k, left) in halves[0].items():
        right = halves[1][m][1]
        detectors.append([before[k,left], before[k,right], previous[k,left|right]])
    for _ in range(1, rounds):
        current = subround(merged)
        detectors.extend([[previous[k], current[k]] for k in current])
        previous = current
    read = [previous[k] for k in seam]

    after = split_round()
    detectors.extend([[previous[k], after[k]] for k in sorted(previous.keys() & after.keys())])
    for m, (k,left) in halves[0].items():
        right = halves[1][m][1]
        detectors.append([previous[k,left|right], after[k,left], after[k,right]])
    for _ in range(1, rounds):
        current = split_round()
        detectors.extend([[after[k], current[k]] for k in current])
        after = current

    # Solve the split syndrome using only the measured seam generators.
    # Each pivot flips exactly one left half and its matching right half.
    left_keys = list(halves[0].values())
    h = len(left_keys)
    augmented = [sum(dot(mask,l) << i for i,(_,l) in enumerate(left_keys)) | mask << h
                 for _,mask in seam]
    solutions, pivots = systematic(augmented, n+h, range(h))
    if len(pivots) != h:
        raise ValueError("seam checks do not span split gauge syndromes")
    live = set(range(n))
    for solution, pivot in zip(solutions, pivots):
        sites = supp(solution >> h)
        if axis == "X":
            emit([["H",q] for q in sites], live)
        emit([["Z_IF",q,after[left_keys[pivot]]] for q in sites], live)
        if axis == "X":
            emit([["H",q] for q in sites], live)
    joint = reduce(xor, (mask for _,mask in seam))
    spec = dict(n=n, checks=[word(n,mask,k) for k,mask in split],
                targets=[word(n,joint,axis)], ticks=ticks,
                detectors=detectors, readouts=[dict(records=read, sign=0)], terminal_syndrome=True)
    # Rotate the doubled checkerboard onto the integer square coupling grid.
    coordinates = [[(r+c)//2,(r-c)//2] for r,c in points]
    return spec, coordinates
