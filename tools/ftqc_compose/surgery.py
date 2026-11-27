"""A checked accepted-Pauli leaf for a standalone or placed projective measurement."""
from itertools import combinations
from math import comb

from ftqc_gadgets.surgery import verify
from .contracts import theorem, Violation
from .supports import Family, fields


def leaf(node):
    fields(node,"op","instance","receipt","wires")
    r=verify(node["receipt"])
    s=r["spec"]
    wires=node["wires"]
    if not isinstance(wires,list) or len(wires)!=s["n"] or len(set(wires))!=len(wires) or not all(isinstance(w,str) and w for w in wires):
        raise ValueError("expected named physical data wires in order")
    instance=node["instance"]
    if not isinstance(instance,str) or not instance:
        raise ValueError("expected instance name")
    if r["faults"]["status"]=="counterexample":
        raise Violation(r["faults"])
    if r["faults"]["status"]!="proved":
        raise ValueError("surgery leaf requires a closed fault exclusion")
    import json
    def port(checks):
        return dict(kind="Q",wires=wires,code=json.dumps(checks),basis="physical signed Pauli",
                    decoder="clean codespace, recorded projective branch",bad=[[q] for q in wires])
    outputs=[port(s["checks"])]
    for i in range(len(s["targets"])):
        w=f"{instance}/bit/{i}"
        outputs.append(dict(kind="C",wires=[w],code="bit",basis="Z",decoder="identity",bad=[[w]]))
    sites=[f"{instance}/loc/{i}" for i in range(sum(not op[0].startswith("!") for row in s["ticks"] for op in row))]
    d=r["faults"]["lower"]
    if comb(len(sites),d)>1000000:
        raise ValueError("surgery support family exceeds enumeration budget")
    return theorem([port(s["checks"])],outputs,dict(project=s["targets"]),Family(sites,combinations(sites,d)),
                   None,dict(circuit=r["circuit_sha256"],readouts=s["readouts"],terminal_syndrome=s["terminal_syndrome"],
                             branches="residual identity modulo checks and measured targets on accepted branch"))

