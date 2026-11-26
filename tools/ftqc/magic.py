"""Exact branch maps. An exponent j denotes exp(i*pi*j/4), modulo 8.
Each postmeasurement data branch is a monomial 2x2 matrix divided by 2.
The columns below are its action on |0> and |1>. This suffices on entangled inputs.
"""
from itertools import product


def branches(dp=0, dx=0, data_x=0, data_z=0, resource_z=0, correction=True):
    results = []
    for p, x in product(range(2), repeat=2):
        columns = []
        for data in (0,1):
            ancilla = p ^ data
            phase = ancilla + 4 * ((x ^ resource_z) * ancilla + data_z * data)
            out = data ^ data_x
            phase += (2 * (p ^ dp) * correction + 4 * (x ^ dx)) * out
            columns.append([out,phase%8])
        relative = (columns[1][1]-columns[0][1])%8
        results.append(dict(p=p,x=x,columns=columns,probability="1/4",
                            implements_T=not data_x and relative==1))
    return results


def damaged(dp,dx,lx,lz):
    return not all(row['implements_T'] for row in branches(dp,dx,lx,lz))


def check_classifier():
    # Compare the SAT predicates with exact branch matrices on every joint label.
    for dp,dx,lx,lz in product(range(2),repeat=4):
        if damaged(dp,dx,lx,lz)!=bool(dp or lx or dx^lz):
            raise ValueError("branch predicate disagrees with exact logical channel")


def frame_observable(p,x):
    # F = S^p Z^x. F-dagger X F = (-1)^(p+x) Y when p=1, else (-1)^x X.
    axis='Y' if p else 'X'
    sign=p^x
    for b in (0,1):
        # Compare F-dagger X F and the signed Pauli on both basis vectors.
        phase=((2*p+4*x)*(2*b-1))%8
        physical=(4*sign+(2-4*b if axis=='Y' else 0))%8
        if physical!=phase:
            raise ValueError("signed readout disagrees with Clifford frame")
    return axis,sign
