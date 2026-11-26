"""Proof-producing circuit search using qeccert's CNF, CaDiCaL and LRAT."""
from itertools import combinations
from pathlib import Path
from tempfile import TemporaryDirectory
from qeccert import cnf, lrat, sat
from algebra import bits, Limit


def formula(schedule,target,reverse=False):
    events=schedule.events(reverse)
    labels=[]
    columns=[]
    groups=[]
    for t,choices in enumerate(events):
        group=[]
        for mode,signature in choices:
            labels.append((t,mode))
            columns.append(signature)
            group.append(len(labels))
        groups.append(group)
    if len(labels)*min(len(labels),target)>2000000:
        raise Limit("circuit SAT cardinality budget")
    rows=[sum(((s>>i)&1)<<j for j,s in enumerate(columns)) for i in range(len(schedule.flows))]
    boundary=len(schedule.det)
    problem=cnf.Formula(len(columns),rows[:boundary],rows[boundary:],target-1)
    for group in groups:
        problem.clauses.extend([[-a,-b] for a,b in combinations(group,2)])
    return problem,labels


def solve(schedule,target,seconds):
    problem,labels=formula(schedule,target)
    with TemporaryDirectory() as tmp:
        result=sat.solve(problem,tmp,'circuit',seconds)
        if result[0]=='unknown':
            return dict(status='inconclusive',reason=result[1])
        if result[0]=='unsat':
            proof=(Path(tmp)/'circuit.lrat').read_text()
            return dict(status='proved',lower=target,lrat=proof)
        faults=[list(labels[j]) for j in bits(result[1])]
        verify_witness(schedule,target,faults)
        return dict(status='counterexample',upper=len(faults),faults=faults)


def verify_witness(schedule,target,faults):
    sites=dict(faults)
    if not 0<len(sites)==len(faults)<target:
        raise ValueError("invalid fault weight")
    for t,m in faults:
        if type(t)!=int or type(m)!=int or not 0<=t<len(schedule.ops) or m not in schedule.modes[t]:
            raise ValueError("invalid fault location or mode")
    value=schedule.signature(sites)
    if not value or value&schedule.mask:
        raise ValueError("fault witness failed independent propagation")


def replay_sat(schedule,target,result):
    if result['status']=='proved':
        problem,_=formula(schedule,target,True)
        lrat.check(problem.dimacs(),result['lrat'].encode())
        if set(result)!={'status','lower','lrat'} or result['lower']!=target:
            raise ValueError("invalid lower bound")
    elif result['status']=='counterexample':
        verify_witness(schedule,target,result['faults'])
        if set(result)!={'status','upper','faults'} or result['upper']!=len(result['faults']):
            raise ValueError("invalid upper bound")
    else:
        raise ValueError("no SAT evidence to replay")
