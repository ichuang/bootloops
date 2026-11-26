"""Single-site product-state injection (Liu et al., Eq. 6), all-plus postselection."""
import hashlib
import json
from qeccert import gf2
from .algebra import bits, echelon, hgp, hgp_parameters, seed
from .magic import branches
from .t_schedule import TSchedule, pull_ideal


def source_record(s,reverse=False):
    if s['resource']!='in-situ-all-plus': raise ValueError('unsupported source')
    gate=TSchedule(s)
    a,na=seed(s['a'])
    b,nb=seed(s['b'])
    n,hx,hz=hgp(a,na,b,nb)
    X,Z=s['logical_x'],s['logical_z']
    common=X&Z
    if common.bit_count()!=1: raise ValueError('source requires a single overlapping site')
    raw=common.bit_length()-1
    plus=X^common
    # V embeds the raw qubit with fixed |+>,|0> factors. X_L V=V X and Z_L V=V Z.
    # Projection on all + checks commutes with both logicals, hence V* Pi V=cI.
    # Evaluate c on raw |0>. The rank below counts independent random syndromes.
    rank=gf2.rank if reverse else lambda rows:len(echelon(rows))
    exponent=rank([r & ~plus for r in hx])+rank([r & plus for r in hz])
    ops=[('RT' if i==raw else 'RX' if plus>>i&1 else 'RZ',i) for i in range(n)]
    for axis,rows in [('X',hx),('Z',hz)]:
        for row in rows:
            if not row: continue
            block=[('R'+axis,n)]
            for q in bits(row):
                block.append(('CX',n,q) if axis=='X' else ('CX',q,n))
            block.append(('M'+axis,n))
            expected=(row,0,0) if axis=='X' else (0,row,0)
            if pull_ideal(block,(0,0,1<<(len(block)-1)))!=expected:
                raise ValueError('source syndrome measures wrong check')
            ops+=block
    # Raw Z is exactly Z_L on the initialized product because all other Z_L
    # factors have eigenvalue +1. It therefore survives accepted projection.
    rows=branches(resource_z=1)
    if any(r['implements_T'] for r in rows) or any((r['columns'][1][1]-r['columns'][0][1])%8!=5 for r in rows):
        raise ValueError('raw fault did not yield Z times T')
    return dict(status='counterexample' if s['target']>1 else 'proved',exact_distance=1,
                resource=s['resource'],base=hgp_parameters(a,na,b,nb),
                source_locations=len(ops),consumer_locations=len(gate.ops),raw_site=raw,
                projection='all measured stabilizers +1',ideal_acceptance=f'2^-{exponent}',
                schedule_sha256=hashlib.sha256(json.dumps(ops,separators=(',',':')).encode()).hexdigest(),
                source_fault=dict(location=raw,after='RT',pauli='Z',accepted_error='logical Z on resource'),
                faulty_branches=rows,output='Z_L T_L',source_guarantee='ideal exact input on acceptance; no stochastic noise bound')
