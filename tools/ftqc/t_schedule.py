"""Two HGP blocks coupled by verified cat-state Pauli meters."""
from functools import reduce
from operator import xor
from .algebra import bits, kernel, hgp, seed, parity
from .magic import check_classifier, frame_observable


def mask(indices):
    return reduce(xor,(1<<i for i in indices),0)


class Meter:
    def __init__(self, first_aux):
        self.first_aux=first_aux
        self.ops=[]
        self.det=[]
        self.parts=[]

    def measure(self,x,z,rounds):
        """Postselect verified GHZ, couple one cat site per operand, read X parity."""
        sites=list(bits(x|z))
        cat=list(range(self.first_aux,self.first_aux+len(sites)))
        meter=cat[-1]+1
        first=None
        for _ in range(rounds):
            start=len(self.ops)
            self.ops.append(('RX',cat[0]))
            for i in range(1,len(cat)):
                self.ops.extend([('RZ',cat[i]),('CX',cat[i-1],cat[i])])
            flags=[]
            for u,v in zip(cat,cat[1:]):
                self.ops.extend([('RZ',meter),('CX',u,meter),('CX',v,meter),('MZ',meter)])
                flags.append(len(self.ops)-1)
                self.det.append((0,0,1<<(len(self.ops)-1)))
            for a,q in zip(cat,sites):
                if x>>q&1 and z>>q&1:
                    self.ops.extend([('SDG',q),('CX',a,q),('S',q)])
                elif x>>q&1:
                    self.ops.append(('CX',a,q))
                else:
                    self.ops.extend([('H',q),('CX',a,q),('H',q)])
            rec=mask(range(len(self.ops),len(self.ops)+len(cat)))
            self.ops.extend(('MX',a) for a in cat)
            self.parts.append((start,len(self.ops),x,z,rec,flags))
            if first is None: first=rec
            else: self.det.append((0,0,first^rec))
        return first


def pull_ideal(ops,flow):
    """Exact adjoint in the convention i^phase X^x Z^z."""
    x,z,rec=flow
    phase=0
    for i in reversed(range(len(ops))):
        g,q,*tail=ops[i]
        bit=1<<q
        if g.startswith('M'):
            if (x|z)&bit: raise ValueError("retired qubit has output flow")
            if rec>>i&1:
                if g=='MX': x ^= bit
                else: z ^= bit
        elif g.startswith('R'):
            if (z if g=='RX' else x)&bit: raise ValueError("random ideal flow")
            x &= ~bit
            z &= ~bit
        elif g=='CX':
            t=1<<tail[0]
            if x&bit: x ^= t
            if z&t: z ^= bit
        elif g=='H':
            if x&z&bit: phase+=2
            if bool(x&bit)!=bool(z&bit): x^=bit; z^=bit
        elif g in ('S','SDG'):
            if x&bit:
                phase+=1 if g=='SDG' else -1
                z^=bit
        elif g!='I': raise ValueError("unsupported physical operation")
    return x,z,phase%4


def verify_meters(meter):
    """Check record projector, ideal flags, and identity on its full commutant."""
    seen=set()
    for start,end,x,z,rec,flags in meter.parts:
        ops=tuple(meter.ops[start:end])
        key=(ops,x,z,rec>>start,tuple(i-start for i in flags))
        if key in seen: continue
        seen.add(key)
        if pull_ideal(ops,(0,0,rec>>start))!=(x,z,(x&z).bit_count()%4):
            raise ValueError("incorrect measured Pauli or phase")
        for i in flags:
            if pull_ideal(ops,(0,0,1<<(i-start)))!=(0,0,0):
                raise ValueError("non-deterministic acceptance flag")
        # Outside the coupled support the schedule acts identically.
        units=[(1<<q,0) for q in bits(x|z)]+[(0,1<<q) for q in bits(x|z)]
        anti=next((a,b) for a,b in units if parity(a,z)^parity(b,x))
        for a,b in units:
            if parity(a,z)^parity(b,x): a^=anti[0]; b^=anti[1]
            if pull_ideal(ops,(a,b,0))!=(a,b,0):
                raise ValueError("ideal meter disturbs commuting operator")


class TSchedule:
    def __init__(self,spec,arm=None):
        n,hx,hz=hgp(*seed(spec['a']),*seed(spec['b']))
        X,Z=spec['logical_x'],spec['logical_z']
        if any(type(w)!=int or not 0<w<1<<n for w in (X,Z)):
            raise ValueError("logical port support out of range")
        if not parity(X,Z) or any(parity(X,r) for r in hz) or any(parity(Z,r) for r in hx):
            raise ValueError("logical ports must be a normalizer pair")
        rounds=spec['rounds']
        if type(rounds)!=int or not 1<=rounds<=32: raise ValueError("invalid meter repetitions")
        check_classifier()
        m=Meter(2*n)
        p=m.measure(0,Z|(Z<<n),rounds)
        x=m.measure(X<<n,0,rounds)
        for op in spec.get('trailing_identity',[]):
            if len(op)!=3 or op[0]!='CX' or not all(type(q)==int and 0<=q<2*n for q in op[1:]) or op[1]==op[2]:
                raise ValueError("invalid identity control")
            m.ops.append(tuple(op))
        extra=spec.get('trailing_identity',[])
        for q in {q for _,a,b in extra for q in (a,b)}:
            for px,pz in [(1<<q,0),(0,1<<q)]:
                if pull_ideal(extra,(px,pz,0))!=(px,pz,0): raise ValueError("tail is not identity")
        for signed_p in (0,1):
            for signed_x in (0,1): frame_observable(signed_p,signed_x)
        # For an actual reported p=1 use the physical Y string; the known sign
        # converts the product of site Ys to logical i X_L Z_L.
        self.y_sign=(((X&Z).bit_count()-1)//2)%2
        self.read_sign=None
        if arm is not None:
            read=m.measure(X,Z if arm else 0,rounds)
            self.read_sign=[frame_observable(arm,b)[1]^(self.y_sign if arm else 0) for b in (0,1)]
        verify_meters(m)
        self.ops=m.ops
        self.parts=m.parts
        self.det=list(m.det)
        for shift in (0,n):
            self.det.extend((row<<shift,0,0) for row in hx)
            self.det.extend((0,row<<shift,0) for row in hz)
        self.mask=(1<<len(self.det))-1
        self.obs=[(0,0,p),(0,Z,0),(X,0,x)]
        # Require identity on the full spectator algebra in BOTH blocks.
        for shift in (0,n):
            self.obs.extend((r<<shift,0,0) for r in kernel(hz+[Z],n))
            self.obs.extend((0,r<<shift,0) for r in kernel(hx+[X],n))
        if arm is not None:
            # Sufficient branch-channel equality: p is correct, final corrected
            # Pauli trivial, and signed follow-on measurement equals ideal.
            self.obs.append((0,0,read^x))
        self.flows=self.det+self.obs
        self.modes=[range(1,16 if g=='CX' else 2 if g[0] in 'RM' else 4) for g,*_ in self.ops]
        self.n=n
        self.arm=arm
