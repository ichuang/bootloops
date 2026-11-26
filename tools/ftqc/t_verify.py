"""Distance evidence for the T output and both subsequent signed-X readout arms."""
from .circuit import search
from .magic import branches, check_readout_classifier
from .t_faults import FaultCircuit, pair_witness
from .proof import solve, replay_sat, verify_witness


def analyze(s,backend,budget,seconds,reverse=False,stored=None):
    if s['resource']!='ideal-encoded-T':
        raise ValueError("conditional gate proof requires an ideal encoded T input")
    check_readout_classifier()
    bound=s['target']
    results={}
    for name,arm in [('gate',None),('read_X',0),('read_Y',1)]:
        circuit=FaultCircuit(s,arm)
        if backend=='sat':
            if stored is None:
                evidence=solve(circuit,bound,seconds)
            else:
                evidence=stored['cases'][name]['evidence']
                replay_sat(circuit,bound,evidence)
        else:
            evidence=search(circuit,bound,budget,reverse)
        row=dict(locations=len(circuit.ops),alternatives=sum(map(len,circuit.modes)),
                 physical_sign_by_x=circuit.read_sign,evidence=evidence)
        if evidence['status']=='proved' and bound==2:
            w=pair_witness(circuit) if stored is None else stored['cases'][name].get('weight_two')
            if w is not None:
                if len(w)!=2: raise ValueError("incorrect matching upper witness")
                verify_witness(circuit,3,w)
                row.update(exact_distance=2,weight_two=w)
        results[name]=row
    statuses=[r['evidence']['status'] for r in results.values()]
    return dict(status='counterexample' if 'counterexample' in statuses else 'inconclusive' if 'inconclusive' in statuses else 'proved',
                ideal_branches=branches(),resource=s['resource'],cases=results)
