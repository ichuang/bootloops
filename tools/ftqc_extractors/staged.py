"""Check an explicit fixed assignment of data, cat couplers and native triples."""
from ftqc_codes.codes import steane
from ftqc_codes.gf2 import supp
from ftqc_gadgets import staged


def install(s):
    record = s["receipt"]
    staged.verify(record)
    result = record["result"]
    if result["status"] != "proved":
        raise ValueError("installation needs a proved staged gadget")
    spec = record["spec"]
    w = result["workload"]
    sites = s["sites"]
    if (not isinstance(sites,list) or len(sites) != w["data"]+w["peak_cat_legs"]
            or any(type(q) is not int or q<0 for q in sites) or len(set(sites))!=len(sites)):
        raise ValueError("expected injective data and cat placement")
    def inventory(items, size):
        if not isinstance(items, list) or any(not isinstance(e,list) or len(e)!=size or
                any(type(v) is not int or v<0 for v in e) or len(set(e))!=size for e in items):
            raise ValueError("invalid coupler or native triple")
        return {tuple(sorted(e)) for e in items}
    couplers = inventory(s["couplers"],2)
    triples = inventory(s["ccz_triples"],3)
    checks = spec.get("checks", [])
    if spec["operation"] == "pieceable-ccz":
        checks = [h << (7*b) for b in range(3) for h in steane().hz]
    used = {tuple(sorted((sites[q], sites[w["data"]+j]))) for h in checks for j,q in enumerate(supp(h))}
    needed = {tuple(sorted(sites[q] for q in gate)) for piece in spec.get("pieces", []) for gate in piece}
    if used-couplers or needed-triples:
        raise ValueError("missing installed check coupler or native CCZ triple")
    return dict(status="proved", stage=record["sha256"], used_couplers=sorted(map(list,used)),
                native_triples=sorted(map(list,needed)), placement=sites,
                scope="cat-data connectivity and native triples; supplied seeds and terminal EC")

