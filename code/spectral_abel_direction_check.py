"""Exact small checks of trace-flow Abel directions, with no BA evaluation.

See the manuscript hierarchy theorem for the normalization and spanning argument.
Time is geometric: delta = i d/dt. Residue vectors use the eigenline
transition convention +lambda**(r-1)/(nu-a); dualizing reverses all signs.
"""
from pathlib import Path
import json
import sympy as s


def power_coefficients(X, D, n):
    """Coefficients of (X+a D)^n without expanding symbolic matrices."""
    coeff = [s.eye(X.rows)]
    for _ in range(n):
        nxt = [s.zeros(X.rows) for _ in range(len(coeff) + 1)]
        for k, C in enumerate(coeff):
            nxt[k] += C * X
            nxt[k + 1] += C * D
        coeff = nxt
    return coeff


def lax_checks():
    X = s.Matrix([
        [-2, 0, 1, 2, 1], [0, -6, 3, -1, 2],
        [2, 1, -4, 0, 3], [-1, 2, 0, -10, -1],
        [-2, -3, -1, -2, 1],
    ])
    D = s.diag(1, 1, -1, -1, 0)
    nu = s.Symbol('nu')
    pencil = X + nu * D
    for r, m in [(3, 0), (4, 0), (5, 0), (5, 1)]:
        A = power_coefficients(X, D, r - 1)
        T = -sum((A[k] * nu ** (k - 1 - m)
                  for k in range(m + 1, r)), s.zeros(5))
        assert all(s.expand(v) == 0 for v in
                   pencil * T - T * pencil - (D * A[m] - A[m] * D))
    return {'sizes': [2, 2, 1], 'trace_and_mixed_lax_checks': 4}


def abel_checks():
    lam, nu, a = s.symbols('lambda nu a')
    # Generic monic degree-five F(lambda,nu). Newton's identities imply
    # these first two coefficients for every pencil, not just the witness.
    t, td, x2, xd, d2 = s.symbols('trX trD trX2 trXD trD2')
    e1 = t + nu * td
    h2 = ((t + nu * td)**2 + x2 + 2 * nu * xd + nu**2 * d2) / 2
    c3, c4, c5 = s.symbols('c3 c4 c5')
    F = lam**5 - e1 * lam**4 + (e1**2 - h2) * lam**3 + c3*lam**2 + c4*lam + c5
    # A holomorphic adjoint basis for the regular (2,2,1) plane model.
    numerators = [s.Integer(1), lam - nu, lam + nu, lam**2 - nu**2]

    def direction(r):
        # Lagrange identity: sum R(lambda_j)/F_lambda(lambda_j)
        # is the lambda^4 coefficient of R mod F.
        return s.Matrix([
            s.expand(s.rem(lam**(r-1)*Q, F, lam)).coeff(lam, 4)
            for Q in numerators
        ])

    U3, U4, U5 = [direction(r).subs(nu, 0) for r in (3, 4, 5)]
    assert U3 == s.Matrix([0, 0, 0, 1])
    assert U4 == s.Matrix([0, 1, 1, t])
    assert U5 == s.Matrix([1, t, t, (t*t+x2)/2])
    # Coefficient of a in tr((X+aD)^5)/5: differentiate the direction
    # of the shifted flow at a=0, with the curve itself fixed.
    U51 = direction(5).diff(nu).subs({nu: 0, td: 0, d2: 4})
    assert U51 == s.Matrix([0, -1, 1, xd])
    A = s.Matrix.hstack(U3, U4, U5, U51)
    assert s.expand(A.det()) == 2

    # Independent residue calculation at one exact unramified fibre.
    roots = [-5, -2, 1, 3, 7]
    poly = s.prod(lam-v for v in roots)
    fp = s.diff(poly, lam)
    for r in (3, 4, 5):
        for Q in numerators:
            R = (lam**(r-1)*Q).subs(nu, 0)
            residue_sum = sum(R.subs(lam,v)/fp.subs(lam,v) for v in roots)
            assert residue_sum == s.rem(R, poly, lam).expand().coeff(lam, 4)

    # A bound valid for every unshifted trace: at nu=0 the middle two
    # numerators agree, so every such direction has equal middle entries.
    assert numerators[1].subs(nu,0) == numerators[2].subs(nu,0)
    return {
        'holomorphic_basis_numerators': [str(q) for q in numerators],
        'U3': [str(v) for v in U3], 'U4': [str(v) for v in U4],
        'U5': [str(v) for v in U5],
        'U_coefficient_a_H5': [str(v) for v in U51],
        'unshifted_trace_direction_rank': 3,
        'with_mixed_quartic_direction_rank': 4,
        'determinant': str(s.expand(A.det())),
        'exact_fibre_residue_checks': 12,
    }


if __name__ == '__main__':
    report = {'lax': lax_checks(), 'abel': abel_checks()}
    destination = Path(__file__).parent / 'results/spectral_abel_direction_check.json'
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
