"""Physical addresses for extractor blocks, bridges and adapters (He et al. §4)."""
from .algebra import bits, kernel, parity
from .extractor import incidence, matching
from .fixed_graph import analyze, pauli, signed_basis


def word(n,x,z,sign=0):
    return ("-" if sign else "+")+"".join("IXZY"[((x>>q)&1)+2*((z>>q)&1)] for q in range(n))


def physical(block,x,z,sign):
    """Embed a signed logical Pauli, including the i phase of each logical Y."""
    pairs=block["basis"]
    px=pz=0
    phase=2*sign+(x&z).bit_count()
    for axis,mask in enumerate((x,z)):
        for i in bits(mask):
            p=pauli(pairs[i][axis],block["n"])
            phase+=(p.x&p.z).bit_count()+2*(p.s+parity(pz,p.x))
            px ^= p.x
            pz ^= p.z
    phase=(phase-(px&pz).bit_count())%4
    if phase%2:
        raise ValueError("non-Hermitian physical Pauli")
    return word(block["n"],px,pz,phase//2)


def local(b,target,degree):
    if set(b)!={"n","checks","vertices","edges","ports","routes","cycles","basis"}:
        raise ValueError("unknown or missing block architecture fields")
    n=b["n"]
    gens=[pauli(w,n) for w in b["checks"]]
    group=signed_basis(gens,n)
    pairs=b["basis"]
    if not pairs or len(pairs)!=n-group.rank or any(len(pair)!=2 for pair in pairs):
        raise ValueError("expected a complete ordered logical X/Z basis")
    for i,pair in enumerate(pairs):
        for a,w in enumerate(pair):
            g=pauli(w,n)
            if any(parity(g.x,h.z)!=parity(h.x,g.z) for h in gens):
                raise ValueError("logical anticommutes with memory")
            for j,other in enumerate(pairs):
                for k,v in enumerate(other):
                    h=pauli(v,n)
                    if (parity(g.x,h.z)^parity(g.z,h.x))!=int(i==j and k!=a):
                        raise ValueError("invalid paired logical basis")
    return {**{k:v for k,v in b.items() if k!="basis"},"target":target,"max_degree":degree}


def join(a,b,links):
    """Disjoint union plus a matching, with one new cycle per extra bridge edge."""
    av,ae,an=a["vertices"],len(a["edges"]),a["n"]
    bv,be,bn=b["vertices"],len(b["edges"]),b["n"]
    if not links or any(len(pair)!=2 for pair in links):
        raise ValueError("bridge requires a nonempty matching")
    left,right=zip(*links)
    if len(set(left))!=len(left) or len(set(right))!=len(right):
        raise ValueError("bridge endpoints overlap")
    if any(u not in a["ports"] or v not in b["ports"] or min(u,v)<0 for u,v in links):
        raise ValueError("bridge endpoints must be ports")
    edges=a["edges"]+[[av+u,av+v] for u,v in b["edges"]]+[[u,av+v] for u,v in links]
    cycles=a["cycles"]+[h<<ae for h in b["cycles"]]
    for i,(u,v) in enumerate(links[1:],1):
        ca=matching(av,a["edges"],(1<<left[0])^(1<<u))
        cb=matching(bv,b["edges"],(1<<right[0])^(1<<v))
        cycles.append(ca|(cb<<ae)|((1|(1<<i))<<(ae+be)))
    return dict(n=an+bn,checks=[w+"I"*bn for w in a["checks"]]+[w[0]+"I"*an+w[1:] for w in b["checks"]],
                vertices=av+bv,edges=edges,cycles=cycles,
                ports=a["ports"]+[av+v if v>=0 else -1 for v in b["ports"]],
                routes=a["routes"]+[h<<ae for h in b["routes"]],
                target=a["target"],max_degree=a["max_degree"])


class Hardware:
    def __init__(self,s,budget):
        self.budget=budget
        self.blocks=s["blocks"]
        self.layouts={}
        self.cache={}
        self.wires={}
        self.labels={}
        locals_=[local(b,s["target"],s["max_degree"]) for b in self.blocks]
        for i,b in enumerate(locals_):
            self.layouts[i,]=b
        for link in s["bridges"]:
            i,j=link["blocks"]
            if not 0<=i<j<len(locals_) or (i,j) in self.layouts:
                raise ValueError("expected unique ordered bridge endpoints")
            self.layouts[i,j]=join(locals_[i],locals_[j],link["matching"])
        for ids,b in self.layouts.items():
            # Select one supported logical to verify the universal wiring hypotheses.
            logical=self.blocks[ids[0]]["basis"][0][1]+"I"*(b["n"]-locals_[ids[0]]["n"])
            cert=self.select(ids,logical)
            def labels(kind,key):
                return [(i,kind,j) for i in ids for j in range(len(locals_[i][key]))]
            data=[(i,"data",q) for i in ids for q in range(locals_[i]["n"])]
            data+=labels("edge","edges")
            checks=labels("memory","checks")+[(i,"vertex",v) for i in ids for v in range(locals_[i]["vertices"])]
            checks+=labels("cycle","cycles")
            data += [(ids,"bridge-edge",i) for i in range(b["n"]+len(b["edges"])-len(data))]
            checks += [(ids,"bridge-cycle",i) for i in range(len(cert["physical_check_couplers"])-len(checks))]
            self.labels[ids]=(data,checks)
            for check,mask in zip(checks,cert["physical_check_couplers"],strict=True):
                self.wires.setdefault(check,set()).update(data[q] for q in bits(mask))
        degrees={}
        for neighbors in self.wires.values():
            for q in neighbors:
                degrees[q]=degrees.get(q,0)+1
        if max([*degrees.values(),*(len(v) for v in self.wires.values())])>s["max_degree"]:
            raise ValueError("aggregate degree of installed block map exceeds bound")

    def select(self,ids,logical):
        if ids not in self.layouts:
            raise ValueError("requested measurement has no installed bridge/adapter")
        key=ids,logical
        if key not in self.cache:
            report=analyze(dict(self.layouts[ids],logical=logical), self.budget)
            if report["status"]!="proved":
                raise ValueError(f"uncertified measurement wiring: {report}")
            self.cache[key]=report
        return self.cache[key]

