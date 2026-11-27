"""Full CSS block reading, Cowtan et al. Definition 3.6. Rows are bit masks."""
from functools import reduce
from operator import xor

from ftqc_codes.gf2 import Echelon
from .algebra import bits, css_dimension, distance_below, kernel, parity, transpose


def combine(rows, coefficients):
    return reduce(xor, (rows[i] for i in bits(coefficients)), 0)


def rows(values, width):
    if not isinstance(values, list) or any(type(v) is not int or v < 0 or v.bit_length() > width for v in values):
        raise ValueError("invalid matrix rows or column width")
    return values


def memory(spec):
    if set(spec) != {"n", "hx", "hz"} or type(spec["n"]) is not int or spec["n"] < 1:
        raise ValueError("memory requires n, hx and hz")
    n = spec["n"]
    code = n, rows(spec["hx"], n), rows(spec["hz"], n)
    if css_dimension(code) < 1:
        raise ValueError("memory encodes no logical qubits")
    return code


def quotient(checks, stabilizers, width):
    basis = Echelon(stabilizers)
    return [v for v in kernel(checks, width) if basis.add(v)]


def copies(word, pattern, width):
    return sum(word << (b * width) for b in bits(pattern))


def split_section(hx, n):
    """Basis of allowed X outcomes and data Z corrections on pivot coordinates."""
    columns = transpose(hx, n)
    pivots = Echelon(columns).picked
    return [[columns[j], 1 << j] for j in pivots]


def row_solver(matrix):
    """Carry coefficients while eliminating rows to lift arbitrary row-space elements."""
    basis={}
    for i,row in enumerate(matrix):
        coeff=1<<i
        for pivot,(v,c) in sorted(basis.items(),reverse=True):
            if row>>pivot&1:
                row^=v
                coeff^=c
        if row:
            basis[row.bit_length()-1]=(row,coeff)
    def solve(word):
        result=0
        for pivot,(v,c) in sorted(basis.items(),reverse=True):
            if word>>pivot&1:
                word^=v
                result^=c
        if word:
            raise ValueError("surviving logical has no lift")
        return result
    return solve


def reading(code, blocks, pattern, axis):
    if axis not in ("X", "Z") or type(blocks) is not int or blocks < 2:
        raise ValueError("invalid reading axis or block count")
    n, hx, hz = code
    if axis == "X":
        hx, hz = hz, hx
    rows(pattern, blocks)
    if not pattern or any(p.bit_count() < 2 for p in pattern) or Echelon(pattern).rank != len(pattern):
        raise ValueError("pattern rows must be independent with at least two participating blocks")
    data = blocks * n
    mx = len(hx)
    columns = transpose(hx, n)
    extend = [copies(1, p, mx) if mx else 0 for p in transpose(pattern, blocks)]
    x = [(h << (b*n)) | (extend[b] << (data+i))
         for b in range(blocks) for i, h in enumerate(hx)]
    old_z = [h << (b*n) for b in range(blocks) for h in hz]
    ports = [copies(1 << q, p, n) | (columns[q] << (data+r*mx))
             for r, p in enumerate(pattern) for q in range(n)]
    z = old_z + ports
    width = data + len(pattern)*mx
    logicals = quotient(hx, hz, n)
    measured = [copies(v, p, n) for p in pattern for v in logicals]
    readouts = [v << (r*n) for r in range(len(pattern)) for v in logicals]
    if [combine(ports, v) for v in readouts] != measured:
        raise ValueError("port readout does not equal requested physical logical")
    # Shorten to data to detect ANY additional projection, including spectator loss.
    expected = Echelon(old_z + measured)
    from .algebra import echelon
    short = [v for v in echelon(z).values() if v.bit_length() <= data]
    if len(short) != expected.rank or any(v not in expected for v in short):
        raise ValueError("unexpected data projection")
    k = css_dimension((width, x, z))
    if k != len(logicals)*(blocks-len(pattern)):
        raise ValueError("unexpected surviving logical dimension")
    metas = [copies(1 << i, p, len(hz)) | (h << (len(old_z)+r*n))
             for r,p in enumerate(pattern) for i,h in enumerate(hz)]
    if any(combine(z,m) for m in metas):
        raise ValueError("invalid port metacheck")
    old_x=[h<<(b*n) for b in range(blocks) for h in hx]
    z_survivors=quotient(old_x,old_z+measured,data)
    x_survivors=quotient(old_z+measured,old_x,data)
    solver=row_solver(hx)
    lifts=[]
    for v in x_survivors:
        lifted=v
        for r,p in enumerate(pattern):
            overlap=reduce(xor,((v>>(b*n))&((1<<n)-1) for b in bits(p)),0)
            lifted |= solver(overlap)<<(data+r*mx)
        if any(parity(lifted,h) for h in z):
            raise ValueError("spectator lift anticommutes with a merged check")
        lifts.append(lifted)
    pairing=[sum(parity(a,b)<<i for i,b in enumerate(z_survivors)) for a in lifts]
    if len(lifts)!=k or len(z_survivors)!=k or Echelon(pairing).rank!=k:
        raise ValueError("incomplete spectator logical algebra")
    return dict(width=width, hx=x if axis=="Z" else z, hz=z if axis=="Z" else x,
                axis=axis, measured=measured, port_readouts=readouts, metachecks=metas,
                survivors=k,spectator_basis=dict(same_axis=z_survivors,dual_axis_lifts=lifts,pairing=pairing), split_basis=split_section(hx,n), split_blocks=pattern,
                added_data=width-data, added_check_meters=len(ports))


def analyze(spec, budget):
    code = memory(spec["memory"])
    result = reading(code, spec["blocks"], spec["pattern"], spec["axis"])
    if "expected" in spec and spec["expected"] != result["measured"]:
        raise ValueError("measured logicals differ from expected ordered representatives")
    witness = distance_below(code, spec["target"], budget)
    if witness:
        return dict(status="counterexample", model="memory-distance", witness=witness)
    return dict(status="proved", model="ideal full block reading, static merged CSS distance",
                distance_lower=spec["target"], theorem="2510.14895 Lemma 3.8", **result)

