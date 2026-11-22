"""Check the bilinear projection of an AHH-type quartic pullback.

The selected Hamiltonian is a specialization of their moment-map
prescription, not a claim that they print this distinguished Hamiltonian.
Time derivatives below are scaled by i: i*d/ds = Hamiltonian rate.
"""
import json
from pathlib import Path
import sympy as sp


def check(q, d):
    F = sp.Matrix(q, d, lambda i, a: sp.Symbol(f'f_{i}_{a}'))
    G = sp.Matrix(q, d, lambda i, a: sp.Symbol(f'g_{i}_{a}'))
    y = sp.symbols(f'y_0:{d}')
    Y = sp.diag(*y)
    M = G.T * F
    H = sp.expand(sp.trace(Y * M * M))
    C = Y*M + M*Y
    Frate = H.diff(G)
    Grate = -H.diff(F)
    expand = lambda A: A.applyfunc(sp.expand)
    assert expand(Frate-F*C) == sp.zeros(q, d)
    assert expand(Grate+G*C.T) == sp.zeros(q, d)
    X = -M
    Xrate = -(Grate.T*F+G.T*Frate)
    assert expand(Xrate-(Y*X**2-X**2*Y)) == sp.zeros(d)
    I = lambda i,a: F[i,a]*G[i,a]
    powers = sum(y[a]*sum(I(i,a) for i in range(q))**2
                 for a in range(d))
    powers += sum((y[a]+y[b])*sum(I(i,a)*I(i,b) for i in range(q))
                  for a in range(d) for b in range(a+1,d))
    rectangles = sum((y[a]+y[b])*(
        F[i,a]*F[j,b]*G[i,b]*G[j,a]
        +F[i,b]*F[j,a]*G[i,a]*G[j,b])
        for a in range(d) for b in range(a+1,d)
        for i in range(q) for j in range(i+1,q))
    assert sp.expand(H-powers-rectangles) == 0
    return dict(q=q,d=d,canonical_rates=True,bilinear_closed_rate=True,
                power_and_rectangle_decomposition=True)


if __name__ == '__main__':
    results = {'quartic_pullback': [check(2,3),check(3,3)],
               'scope': 'Symbolic parent projection. The separate all-size '
                        'block and phase check is ahh_general_multiwave_map.py. '
                        'No canonical equivalence or pure-BA dictionary for '
                        'the projected entries is asserted.'}
    target = Path(__file__).parent/'results'/'ahh_quartic_to_cubic.json'
    target.write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))
