"""Bind compiled signed PBC requests to individually certified physical settings."""
from ftqc_gadgets.surgery import Group,commute,pauli,receipt
from .installed_surgery import compile_schedule


def bind(compiled,source,budget):
    menu=source["surgery_installation"]
    if set(menu)!={"couplers","settings","idles"}:
        raise ValueError("expected a fixed physical surgery inventory")
    blocks=source["blocks"]
    sizes=[b["n"] for b in blocks]
    n=sum(sizes)
    offsets=[sum(sizes[:j]) for j in range(len(sizes))]
    checks=[w[0]+"I"*base+w[1:]+"I"*(n-base-size)
            for b,size,base in zip(blocks,sizes,offsets) for w in b["checks"]]
    stabilizers=Group(n,[pauli(w,n) for w in checks])
    options=[]
    for part in menu["settings"]:
        spec=part["gadget"]["spec"]
        if len(spec["targets"])!=1 or set(part)!={"gadget","placement"}:
            raise ValueError("PBC settings must each produce one public parity")
        w=spec["targets"][0]
        chars=["I"]*n
        for i,c in enumerate(w[1:]):
            q=part["placement"][i]
            if type(q)!=int or not 0<=q<n:
                raise ValueError("invalid PBC data placement")
            chars[q]=c
        target=w[0]+"".join(chars)
        physical,peak=compile_schedule(dict(n=n,checks=checks,batches=[[part]],targets=[target],
            couplers=menu["couplers"],idles=menu["idles"],terminal_syndrome=True))
        proof=receipt(physical,source["target"]-1,budget)
        if proof["faults"]["status"]!="proved":
            raise ValueError("PBC surgery exclusion incomplete: "+str(proof["faults"]))
        options.append((pauli(target,n),proof,peak))
    links={tuple(sorted(p)) for p in menu["couplers"]}
    from collections import Counter
    if max(Counter(q for p in links for q in p).values(),default=0)>source["max_degree"]:
        raise ValueError("installed surgery degree exceeds PBC maximum")
    selected={}
    for key,reading in compiled["settings"].items():
        parts=iter(reading["logical"][1:])
        chars=["I"]*n
        for i in reading["blocks"]:
            sl=slice(offsets[i],offsets[i]+sizes[i])
            chars[sl]=[next(parts) for _ in range(sizes[i])]
        want=pauli(reading["logical"][0]+"".join(chars),n)
        matches=[]
        for index,(actual,proof,peak) in enumerate(options):
            if not commute(want,actual):
                continue
            from ftqc_faults.audit.stab import Gen
            delta=Gen(want.x,want.z,want.s)
            delta.mul(actual)
            rem=stabilizers.reduce(delta)
            if rem.x or rem.z:
                continue
            matches.append((len(proof["spec"]["ticks"]),index,rem.s,peak))
        if not matches:
            raise ValueError("no verified installed surgery for "+key)
        ticks,index,flip,peak=min(matches)
        proof=options[index][1]
        selected[key]=dict(setting=index,public_xor=flip,ticks=ticks,peak_live=peak,
                           circuit=proof["circuit_sha256"],faults=proof["faults"])
    for branch in compiled["branches"]:
        time=0
        for step in branch["program"]:
            if "measurement" in step:
                step["surgery_start_tick"]=time
                time+=selected[step["measurement"]]["ticks"]
        branch["surgery_ticks"]=time
    return dict(selected=selected,couplers=[list(p) for p in sorted(links)],
                scope="per-setting exclusion at ideal codespace boundaries; serialized measurement ticks")

