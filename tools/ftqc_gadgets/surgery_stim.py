"""Independent Stim boundary checks and circuit-level noise for signed surgery."""
from ftqc_codes.gf2 import dot,supp
from .surgery import centralizer,paired,pauli,record_mask
from ftqc_faults.audit.stab import Gen


def compile_probe(c,spec,initial,noise=0.0):
    import stim
    if not 0<=noise<1:
        raise ValueError("invalid physical error probability")
    n,shift=spec["n"],c.n
    width=n+shift
    def ps(g):
        p=stim.PauliString(width)
        for q in supp(g.x|g.z):
            p[q]=(1 if g.x>>q&1 else 0)^(3 if g.z>>q&1 else 0)
        p.sign=-1 if g.s else 1
        return p
    out=stim.Tableau.from_stabilizers([ps(g) for g in initial]).to_circuit()
    for tick in c.ticks[1:]:
        for op in tick:
            name=op.name
            qs=list(op.qs)
            if name=="Z_IF":
                out.append("CZ",[stim.target_rec(op.record-out.num_measurements),qs[0]])
            elif name in ("MX","M"):
                out.append(name,qs,noise if op.noisy else 0)
                continue
            else:
                out.append({"SDG":"S_DAG"}.get(name,name),qs)
            if noise and op.noisy:
                out.append("DEPOLARIZE2" if name=="CX" else "DEPOLARIZE1",qs,noise)
        out.append("TICK")
    physical_records=out.num_measurements
    count=0
    conditions=[]
    def condition(g,records,detector):
        nonlocal count
        ids=list(records)
        if g is not None:
            if g.x|g.z:
                ids.append(out.num_measurements)
                terms=[]
                for q in supp(g.x|g.z):
                    kind=(g.x>>q&1,g.z>>q&1)
                    factory={(1,0):stim.target_x,(0,1):stim.target_z,(1,1):stim.target_y}[kind]
                    if terms:
                        terms.append(stim.target_combiner())
                    terms.append(factory(q,invert=bool(g.s and not terms)))
                out.append("MPP",terms)
            elif g.s:
                raise ValueError("inconsistent identity sign")
        conditions.append(ids)
        rec=[stim.target_rec(i-out.num_measurements) for i in ids]
        if detector:
            out.append("DETECTOR",rec)
        else:
            out.append("OBSERVABLE_INCLUDE",rec,count)
            count+=1
    for ids in spec["detectors"]:
        record_mask(ids,physical_records)
        condition(None,ids,True)
    checks=[pauli(w,n) for w in spec["checks"]]
    targets=[pauli(w,n) for w in spec["targets"]]
    for g in checks:
        condition(g,[],spec["terminal_syndrome"])
        condition(Gen(g.x<<shift,g.z<<shift,g.s^dot(g.x,g.z)),[],False)
    for g in centralizer(checks+targets,n):
        condition(paired(g,shift),[],False)
    for target,read in zip(targets,spec["readouts"],strict=True):
        sign=read["sign"]^target.s
        condition(Gen(target.x,target.z,sign),read["records"],False)
        condition(Gen(target.x<<shift,target.z<<shift,sign^dot(target.x,target.z)),read["records"],False)
    sample=out.reference_sample()
    if any(sum(int(sample[i]) for i in ids)%2 for ids in conditions):
        raise ValueError("reference branch has wrong sign")
    # DEM construction checks determinism on every branch, including private records.
    out.detector_error_model()
    return out

