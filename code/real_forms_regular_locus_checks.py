"""Exact checks supporting the signed real forms and node-smoothing argument.

Run from the repository root with .venv/bin/python. The local determinant
checks illustrate the general Schur-complement calculation; they do not
replace the compactness and connectedness argument in the manuscript.
"""
from itertools import product
import sympy as sp

u, v, w = sp.symbols('u v w', complex=True)
bu, bv, bw, coupling = sp.symbols('bu bv bw coupling', real=True)
invariant = []
for signs in product((1, -1), repeat=3):
    modes = (u, v, w)
    partners = tuple(e*sp.conjugate(a) for e, a in zip(signs, modes))
    rates, partner_rates = [], []
    for k, b in enumerate((bu, bv, bw)):
        j, h = (k+1) % 3, (k+2) % 3
        rates.append(-sp.I*(b*modes[k]+coupling*partners[j]*partners[h]))
        partner_rates.append(sp.I*(b*partners[k]+coupling*modes[j]*modes[h]))
    ok = all(sp.simplify(db-e*sp.conjugate(da)) == 0
             for db, da, e in zip(partner_rates, rates, signs))
    assert ok == (sp.prod(signs) == 1)
    if ok:
        invariant.append(signs)
print('Invariant signed real forms:', invariant)

eps = sp.symbols('eps')
for sizes in ((2, 2, 2), (1, 2, 3)):
    slopes = [s for n, s in zip(sizes, (1, -1, 0)) for _ in range(n)]
    d = sum(sizes)
    intercepts = [sp.Integer(2)**k for k in range(d)]
    # All off-group products are nonzero; no within-group perturbation.
    A = sp.Matrix(d, d, lambda i, j:
                  0 if slopes[i] == slopes[j] else 1+2*i+3*j)
    checked = 0
    for i in range(d):
        for j in range(i+1, d):
            if slopes[i] == slopes[j]:
                continue
            nu = (intercepts[j]-intercepts[i])/sp.Integer(slopes[i]-slopes[j])
            lam = intercepts[i]+slopes[i]*nu
            diagonal = [lam-intercepts[k]-slopes[k]*nu for k in range(d)]
            assert all(diagonal[k] != 0 for k in range(d) if k not in (i, j))
            outside = sp.prod(diagonal[k] for k in range(d) if k not in (i, j))
            determinant = (sp.diag(*diagonal)-eps*A).det(method='domain-ge')
            coefficient = sp.Poly(determinant, eps).coeff_monomial(eps**2)
            assert coefficient == -A[i, j]*A[j, i]*outside != 0
            checked += 1
    print(sizes, ': nonzero node-smoothing coefficients at', checked, 'intersections')
