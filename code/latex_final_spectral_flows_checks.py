"""Exact checks of the all-size trace-flow formulas in the manuscript.

Symbolic scalar blocks verify the canonical bracket and quartic expansion;
rational larger blocks check the same identities with all three sizes > 1.
The matrix commutator argument in the manuscript proves the general result.
"""
import sympy as s


def check(L, M, N, symbolic=False):
    count = 0

    def matrix(name, rows, cols):
        nonlocal count
        if symbolic:
            return s.Matrix(rows, cols, lambda i, j: s.Symbol(f'{name}_{i}_{j}'))
        entries = [s.Rational((count+k)*7 % 17-8, 3) for k in range(rows*cols)]
        count += rows*cols
        return s.Matrix(rows, cols, entries)

    u, us = matrix('u', L, M), matrix('us', L, M)
    v, vs = matrix('v', M, N), matrix('vs', M, N)
    w, ws = matrix('w', N, L), matrix('ws', N, L)
    diagonal = (list(s.symbols(f'b0:{L+M+N}')) if symbolic
                else [s.Integer(i*i+1) for i in range(L+M+N)])
    A, B, C = s.diag(*diagonal[:L]), s.diag(*diagonal[L:L+M]), s.diag(*diagonal[L+M:])
    X = s.BlockMatrix([[A, s.sqrt(2)*u, ws.T],
                       [s.sqrt(2)*us.T, B, v], [-w, -vs.T, C]]).as_explicit()
    D = s.diag(*([1]*L+[-1]*M+[0]*N))
    Y = X-s.diag(*diagonal)
    fourth = (2*s.trace((u*us.T)**2)+s.trace((ws.T*w)**2)/2
              +s.trace((v*vs.T)**2)/2-2*s.trace(u*us.T*ws.T*w)
              -2*s.trace(us.T*u*v*vs.T)+s.trace(w*ws.T*vs.T*v))
    assert s.expand(s.trace(Y**4)/4-fourth) == 0

    nu = s.Symbol('nu')
    pencil = X+nu*D
    for r in (2, 3, 4, 5):
        gradient = X**(r-1)
        T = (gradient-pencil**(r-1)).applyfunc(lambda e: s.cancel(e/nu))
        assert all(s.denom(e).has(nu) is False for e in T)
        assert all(s.expand(e) == 0 for e in
                   pencil*T-T*pencil-(D*gradient-gradient*D))
    if symbolic:
        q, p = list(u)+list(v)+list(w), list(us)+list(vs)+list(ws)
        H4 = s.expand(s.trace(X**4)/4)
        def velocity(f):
            return s.expand(sum(s.diff(f,a)*s.diff(H4,b)-s.diff(f,b)*s.diff(H4,a)
                                for a,b in zip(q,p)))
        assert all(s.expand(e) == 0 for e in
                   X.applyfunc(velocity)-(D*X**3-X**3*D))
    print(f'{L},{M},{N}: quartic mode expansion and trace Lax flows passed')


if __name__ == '__main__':
    check(1, 1, 1, symbolic=True)
    check(2, 2, 2)
    check(2, 3, 2)
