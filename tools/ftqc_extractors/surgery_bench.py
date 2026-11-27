"""Reproducible finite architecture studies and source generation."""
import argparse
import json
from pathlib import Path
from ftqc_gadgets.surgery import Instrument,receipt
from ftqc_gadgets.surgery_design import build
from ftqc_gadgets.surgery_merge import word
from .algebra import hgp,seed
from .surgery_export import export
from .surgery_products import radial_design,product_sequence


def installed(gadgets,n,checks,target,idles=False):
    parts=[]
    links=set()
    targets=[]
    for spec,pos,control in gadgets:
        local=receipt(spec,0)
        setting=dict(gadget=local,placement=pos)
        if control is not None:
            setting["control"]=control
        parts.append([setting])
        for row in spec["ticks"]:
            for name,*args in row:
                if name=="CX":
                    links.add(tuple(sorted(pos[i] for i in args)))
        for w in spec["targets"]:
            out=["I"]*n
            for j,a in enumerate(w[1:]):
                out[pos[j]]=a
            targets.append(w[0]+"".join(out))
    return dict(schema=1,kind="installed_surgery",n=n,checks=checks,targets=targets,
                couplers=[list(p) for p in sorted(links)],batches=parts,idles=idles,
                terminal_syndrome=True,target=target)


def two_block(variant="adapter"):
    source=json.loads((Path(__file__).parent/"examples/edge-extractor.json").read_text())
    first=export(source)
    n,x,z=hgp(*seed(source["a"]),*seed(source["b"]))
    hx=x+[v<<n for v in x]
    hz=z+[v<<n for v in z]
    common=dict(n=2*n,hx=hx,hz=hz,axis="Z",rounds=2,targets=[word(2*n,48|(5<<n),"Z")],
                sites=[4,5,8,10])
    if variant=="adapter":
        shape=dict(common,method="adapter",left=[[0,1]],right=[[0,1]])
    else:
        shape=dict(common,method="parsimonious",edges=[[0,1],[1,2],[2,3]])
    bridge=build(shape)["spec"]
    follow=build(dict(n=n,hx=x,hz=z,axis="X",sites=[4,6],edges=[[0,1]],
                      rounds=2,targets=[word(n,80,"X")],method="gauging"))["spec"]
    def count(s):
        return 1+max(q for row in s["ticks"] for op in row for q in (op[1:2] if op[0]=="Z_IF" else op[1:]))
    return installed([(first,list(range(n))+[2*n,2*n+1],None),
                      (bridge,list(range(count(bridge))),None),
                      (follow,list(range(n,2*n))+list(range(2*n,2*n+count(follow)-n)),[0,[12,13]])],
                     2*n,bridge["checks"],2)


def architecture(case):
    from .surgery_layout import product_connection,translated_settings
    from .surgery_products import lifted_product
    if case in ("shuttle-packed","shuttle-serial"):
        from ftqc_gadgets.surgery_parallel import branch
        from .installed_surgery import compile_schedule
        g=branch(1,[],["+Z"],2)
        s=installed([(g,[0,2],None),(g,[1,3],None)],2,[],2)
        if case=="shuttle-packed":
            s["batches"]=[s["batches"][0]+s["batches"][1]]
        spec,_=compile_schedule(s)
        layers=[]
        for i,row in enumerate(spec["ticks"]):
            moves=[[2,[1,0]],[3,[1,3]]] if i==0 else []
            if i==len(spec["ticks"])-1:
                moves=[[2,[3,0]],[3,[3,3]]]
            layers.append(dict(moves=moves,move_time=2 if moves else 0,gate_time=1))
        s["transport"]=dict(initial=[[0,0],[0,3],[3,0],[3,3]],layers=layers,gate_range=1)
        return s
    if case=="product-connection":
        return product_connection([[5]],[[5]],[[1]],[[3]],4,2)
    if case=="pcpp":
        return dict(method="linear",n=8,hx=[15,240],hz=[15,240],axis="Z",ports=[1,2,4,0],
                    aux=[11,12],rounds=2,distance_target=2,targets=[word(8,i,"Z") for i in (5,6)])
    if case=="cyclic-seed":
        n,x,z=lifted_product([[5]],[[5]],4)
        from ftqc_codes.gf2 import nullspace
        for seed_word in nullspace(x,n):
            try:
                designs=translated_settings(n,x,z,4,seed_word,[0,1],2)
            except ValueError:
                continue
            parts=[build(d)["spec"] for d in designs]
            return installed([(p,list(range(Instrument(p).c.n)),None) for p in parts],n,parts[0]["checks"],2)
        raise ValueError("no independent translated seeds")
    # Two inequivalent small CSS blocks. Both joint measurements touch data 0 and 3.
    from ftqc_gadgets.surgery_parallel import branch
    from ftqc_gadgets.surgery_build import graph
    checks=[word(7,r,"Z") for r in (7,24)]
    targets=[word(7,r,"Z") for r in (9,43)]
    if case=="branched":
        parts=[branch(7,checks,targets,2)]
    elif case=="parallel-parsimonious":
        parts=[]
        for m in (9,43):
            sites=[j for j in range(7) if m>>j&1]
            parts.append(build(dict(n=7,hx=[],hz=[7,24],targets=[word(7,m,"Z")],axis="Z",rounds=2,
                         method="parsimonious",sites=sites,edges=[[i,i+1] for i in range(len(sites)-1)]))["spec"])
    else:
        parts=[graph(7,checks,t,[(j,j+1) for j in range(t.count("Z")-1)],"homomorphic",2) for t in targets]
    parts.append(graph(7,checks,word(7,59,"X"),[(j,j+1) for j in range(4)],"homomorphic",2))
    return installed([(p,list(range(Instrument(p).c.n)),None) for p in parts],7,checks,2)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("case",choices=["two-block","two-block-parsimonious","radial90","radial90-reference","hgp-padded","shuttle-packed","shuttle-serial","product-connection","cyclic-seed","pcpp","branched","serial","parallel-parsimonious"])
    ap.add_argument("-o",required=True,type=Path)
    args=ap.parse_args()
    if args.case.startswith("two-block"):
        result=two_block("adapter" if args.case=="two-block" else "parsimonious")
    elif args.case=="hgp-padded":
        result=product_sequence([7],3,[3,3],2,[3,5],2,[2,2],12)
    elif not args.case.startswith("radial"):
        result=architecture(args.case)
    else:
        result=radial_design()
        if args.case.endswith("reference"):
            sites=[i for i,c in enumerate(result["targets"][0][1:]) if c=="Z"]
            result={k:v for k,v in result.items() if k not in ("ports","aux")}
            result.update(method="gauging",sites=sites,edges=[[i,i+1] for i in range(len(sites)-1)],rounds=10)
    args.o.write_text(json.dumps(result,indent=1)+"\n")


if __name__=="__main__":
    main()

