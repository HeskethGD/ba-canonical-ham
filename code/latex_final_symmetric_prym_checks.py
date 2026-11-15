"""Exact fold checks in the three-matrix notation of the manuscript.

The genus and divisor arguments are proved in the manuscript. These checks
verify preservation of the constraints and both matrix dualities, retaining
arbitrary symbolic endpoint phases and the central phase.
"""
import sympy as s


def check(k):
    def symmetric(prefix):
        return s.Matrix(k, k, lambda i, j: s.Symbol(
            f'{prefix}_{min(i, j)}_{max(i, j)}'))

    U, Us = symmetric('u'), symmetric('us')
    V = s.Matrix(s.symbols(f'v0:{k}'))
    Vs = s.Matrix(s.symbols(f'vs0:{k}'))
    B = s.diag(*s.symbols(f'b0:{k}'))
    c, nu = s.symbols('central_phase nu')
    zero = s.zeros(k)
    X = s.BlockMatrix([
        [-2*B, s.sqrt(2)*U, Vs],
        [s.sqrt(2)*Us, -2*B, V],
        [-V.T, -Vs.T, s.Matrix([[c]])],
    ]).as_explicit()
    D = s.diag(*([1]*k + [-1]*k + [0]))
    S = s.BlockMatrix([
        [zero, s.eye(k), s.zeros(k, 1)],
        [s.eye(k), zero, s.zeros(k, 1)],
        [s.zeros(1, k), s.zeros(1, k), s.Matrix([[-1]])],
    ]).as_explicit()
    T = s.BlockMatrix([
        [-B, zero, Vs], [zero, B, -V],
        [-V.T, Vs.T, s.Matrix([[-nu]])],
    ]).as_explicit()
    assert S*S == s.eye(2*k+1)
    assert X == S*X.T*S
    assert S*D*S == -D
    assert T == -S*T.subs(nu, -nu).T*S
    # Differentiate the full Lax equation: the folded tangent remains folded.
    dX = ((X+nu*D)*T-T*(X+nu*D)).applyfunc(s.expand)
    assert dX == S*dX.T*S
    assert all(nu not in entry.free_symbols for entry in dX)
    assert dX[:k, :k] == zero
    assert dX[k:2*k, k:2*k] == zero
    assert dX[2*k, 2*k] == 0
    genus = k*k
    quotient = k*(k-1)//2
    ramification = 2*k+2
    assert 2*genus-2 == 2*(2*quotient-2)+ramification
    assert genus-quotient == k*(k+1)//2
    print(f'K={k}: fold, pencil/time dualities, genus and Prym count passed')


if __name__ == '__main__':
    for size in (1, 2, 3):
        check(size)
