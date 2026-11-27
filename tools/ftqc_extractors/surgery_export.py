"""Export a legacy HGP measurement as a self-contained Gadgets instrument."""
from .algebra import hgp, seed
from .__main__ import validate


def export(spec):
    validate(spec)
    if spec["kind"] != "circuit":
        raise ValueError("expected an HGP measurement circuit")
    n,hx,hz = hgp(*seed(spec["a"]),*seed(spec["b"]))
    def word(mask,axis):
        return "+"+"".join(axis if mask >> i & 1 else "I" for i in range(n))
    return dict(n=n,checks=[word(x,"X") for x in hx]+[word(z,"Z") for z in hz],
                targets=[word(spec["measured_z"],"Z")],
                ticks=[[["R" if g=="RZ" else g,*qs]] for g,*qs in spec["ops"]],
                readouts=[dict(records=spec["readout"],sign=0)],
                detectors=spec["detectors"],terminal_syndrome=True)

