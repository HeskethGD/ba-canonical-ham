"""General L,M,N map from AHH Part II (1990), equation (2.1).

This checks an equation-level solution transfer, not a canonical coordinate
identification with the separate Moser variables F,G.

Use their r=d=L+M+N, degree m=1 and source loop variable lambda=nu.
Their source evolution is
    d_s Lhat = [Bhat,Lhat], Bhat=(Lhat**2/nu)_+.
Write L=X+nu D with D=diag(Id_L,-Id_M,0_N) and the mode blocks X
from the manuscript. Define Kph=diag(-3 B_L,3 B_M,0_N).
Then our solutions are
    L(nu,t)=exp(-i t Kph) Lhat(nu,i t) exp(i t Kph).
Indeed B=(L**2/nu)_+=nu D**2+D X+X D and
    T_ours=B-Kph-nu Id_d,
so i d_t L=[L,T_ours] exactly.

The source initial polynomial is Lhat(nu,0)=X(initial)+nu D. Its
reconstructed off-diagonal blocks give u,v,w and their independent
partners with the fixed sqrt(2), transpose and sign factors below.
The time/phase conversion restores our original cubic Hamiltonian.

All phase constants and all modes are symbolic. The checks include unequal
group sizes and a case with L,M,N all greater than one. Since the displayed
block identities are dimension-independent, these tests certify the
implementation of the general formula, rather than establish its scope
by extrapolation from small cases.
"""

from pathlib import Path
import json
import sympy as sp


def check(group_sizes):
    L, M, N = group_sizes
    d = L + M + N
    nu = sp.Symbol('nu')

    def array(name, rows, cols):
        return sp.Matrix(rows, cols,
                         lambda i, j: sp.Symbol(f'{name}_{i+1}_{j+1}'))

    u, ub = array('u', L, M), array('ubar', L, M)
    v, vb = array('v', M, N), array('vbar', M, N)
    w, wb = array('w', N, L), array('wbar', N, L)
    bL = sp.symbols(f'bL_1:{L+1}')
    bM = sp.symbols(f'bM_1:{M+1}')
    bN = sp.symbols(f'bN_1:{N+1}')
    BL, BM, BN = sp.diag(*bL), sp.diag(*bM), sp.diag(*bN)
    X = sp.BlockMatrix([
        [-2*BL, sp.sqrt(2)*u, wb.T],
        [sp.sqrt(2)*ub.T, -2*BM, v],
        [-w, -vb.T, BN],
    ]).as_explicit()
    D = sp.diag(*([1]*L+[-1]*M+[0]*N))
    Kph = sp.diag(*([-3*b for b in bL]+[3*b for b in bM]+[0]*N))
    pencil = X + nu*D
    B = nu*D**2 + D*X + X*D
    Tours = sp.BlockMatrix([
        [-BL, sp.zeros(L, M), wb.T],
        [sp.zeros(M, L), BM, -v],
        [-w, vb.T, -nu*sp.eye(N)],
    ]).as_explicit()

    zero = sp.zeros(d)
    expand = lambda A: A.applyfunc(sp.expand)
    assert expand(pencil**2 - (nu*B + X**2)) == zero
    assert expand(Tours - B + Kph + nu*sp.eye(d)) == zero
    source_rate = expand(B*pencil-pencil*B)
    assert source_rate == expand(D*X**2-X**2*D)
    target_rate = expand(-source_rate + Kph*X-X*Kph)
    assert target_rate == expand(pencil*Tours-Tours*pencil)

    H = sum((bL[i]+bM[j])*u[i,j]*ub[i,j]
            for i in range(L) for j in range(M))
    H += sum((bM[j]+bN[k])*v[j,k]*vb[j,k]
             for j in range(M) for k in range(N))
    H += sum((bN[k]+bL[i])*w[k,i]*wb[k,i]
             for k in range(N) for i in range(L))
    H += sp.sqrt(2)*sum(u[i,j]*v[j,k]*w[k,i]+ub[i,j]*vb[j,k]*wb[k,i]
                        for i in range(L) for j in range(M) for k in range(N))
    hamilton_rate = sp.zeros(d)
    for a, abar in ((u,ub),(v,vb),(w,wb)):
        for x, xb in zip(a, abar):
            hamilton_rate += X.diff(x)*sp.diff(H,xb)
            hamilton_rate -= X.diff(xb)*sp.diff(H,x)
    assert expand(hamilton_rate-target_rate) == zero

    H3 = sp.trace(X**3)/3
    charges = 3*sp.trace(BL*(wb.T*w-u*ub.T))
    charges += 3*sp.trace(BM*(v*vb.T-ub.T*u))
    constant = -sp.Rational(8,3)*sp.trace(BL**3)
    constant -= sp.Rational(8,3)*sp.trace(BM**3)
    constant += sp.trace(BN**3)/3
    assert sp.expand(H+H3-charges-constant) == 0
    # Directly verify the general mode bracket becomes the frozen matrix bracket.
    basic = []
    for a, abar in ((u,ub),(v,vb),(w,wb)):
        basic.extend(zip(a,abar))
    for h in range(d):
        for i in range(d):
            for j in range(d):
                for k in range(d):
                    bracket = sum(
                        X[h,i].diff(x)*X[j,k].diff(xb)
                        - X[h,i].diff(xb)*X[j,k].diff(x)
                        for x,xb in basic
                    )
                    expected = (D[h,h]-D[i,i]) if (h==k and i==j) else 0
                    assert sp.expand(bracket-expected) == 0
    return {'L_M_N': list(group_sizes),
            'source_polynomial_flow': True,
            'time_and_phase_map': True,
            'all_canonical_mode_equations': True,
            'cubic_hamiltonian_trace_and_charges': True,
            'canonical_to_frozen_bracket': True}


if __name__ == '__main__':
    result = {
        'checks': [check(sizes) for sizes in ((1,2,2),(2,3,1),(2,2,2))],
        'source': 'AHH 1990 Part II, eqs. (2.1)-(2.2), P(z,nu^-1)=z^2/nu',
        'scope': 'Explicit all-size polynomial-Lax solution transfer. '
                 'Does not assert an F,G canonical reduction or BA regularity '
                 'at every exceptional spectral level.',
    }
    out = Path(__file__).parent/'results'/'ahh_general_multiwave_map.json'
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
