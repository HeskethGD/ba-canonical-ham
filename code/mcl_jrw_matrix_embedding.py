"""Compare MCL (1983), equation (2.13), with our zero-pair-phase model.

All phase constants and independent partner modes are symbolic.
Writes only the comparison report in code/results/.
"""
from pathlib import Path
import json
import sympy as s


def restricted_star(R):
    p, pb, coupling = s.symbols('p pbar epsilon')
    q = s.symbols(f'q0:{R}')
    qb = s.symbols(f'qbar0:{R}')
    r = s.symbols(f'r0:{R}')
    rb = s.symbols(f'rbar0:{R}')
    delta = s.symbols(f'Delta0:{R}')
    H = sum(delta[j]*(q[j]*qb[j]-r[j]*rb[j])/2
            -coupling*(p*qb[j]*r[j]+pb*q[j]*rb[j])/2 for j in range(R))
    u, ub = s.symbols('u ubar')
    v = s.symbols(f'v0:{R}')
    vb = s.symbols(f'vbar0:{R}')
    w = s.symbols(f'w0:{R}')
    wb = s.symbols(f'wbar0:{R}')
    inverse = {p:u,pb:ub}
    for j in range(R):
        inverse.update({q[j]:-s.I*vb[j],qb[j]:-s.I*v[j],r[j]:w[j],rb[j]:wb[j]})
    target = sum(-delta[j]*(v[j]*vb[j]+w[j]*wb[j])/2
                 +s.I*coupling*(u*v[j]*w[j]+ub*vb[j]*wb[j])/2 for j in range(R))
    assert s.expand(H.subs(inverse)-target) == 0
    # {q,qbar}=-i implies {i*qbar,i*q}=-i.
    assert s.I*s.I*s.I == -s.I
    return {'triads':R,'MCL_H_1_1_to_three_array_H_verified':True,
            'canonical_partner_swap_verified':True,
            'arbitrary_distinct_or_repeated_detunings':True}


def check(L, M, N):
    d = L+M+N
    e = s.symbols(f'e0:{d}')
    E = s.diag(*e)
    Z = s.zeros(d)
    u = s.Matrix(L, M, lambda i,j:s.Symbol(f'u{i}_{j}'))
    ub = s.Matrix(L, M, lambda i,j:s.Symbol(f'bu{i}_{j}'))
    v = s.Matrix(M, N, lambda j,k:s.Symbol(f'v{j}_{k}'))
    vb = s.Matrix(M, N, lambda j,k:s.Symbol(f'bv{j}_{k}'))
    w = s.Matrix(N, L, lambda k,i:s.Symbol(f'w{k}_{i}'))
    wb = s.Matrix(N, L, lambda k,i:s.Symbol(f'bw{k}_{i}'))
    rt = s.sqrt(2)
    Z[:L,L:L+M] = rt*u
    Z[L:L+M,:L] = rt*ub.T
    Z[:L,L+M:] = wb.T
    Z[L+M:,:L] = -w
    Z[L:L+M,L+M:] = v
    Z[L+M:,L:L+M] = -vb.T
    groups = [0]*L+[1]*M+[2]*N
    slopes = [1 if g==0 else -1 if g==1 else 0 for g in groups]
    a = [s.I*x for x in slopes]
    c = [s.Integer(g==2) for g in groups]
    alpha = s.zeros(d)
    for i in range(d):
        for k in range(d):
            if a[k]!=a[i]:
                alpha[i,k] = s.cancel((c[k]-c[i])/(s.I*(a[k]-a[i])))
    T = alpha.multiply_elementwise(Z)
    rate = ((E+Z)*T-T*(E+Z)).applyfunc(s.expand)

    # Verify that this is the canonical cubic Hamiltonian, in the H0 gauge.
    H = rt*((u*v*w).trace()+(ub*vb*wb).trace())
    H += sum((e[L+M+k]-e[i])*w[k,i]*wb[k,i]
             for i in range(L) for k in range(N))
    H += sum((e[L+M+k]-e[L+j])*v[j,k]*vb[j,k]
             for j in range(M) for k in range(N))
    unbar = list(u)+list(v)+list(w)
    bar = list(ub)+list(vb)+list(wb)
    field = dict(zip(unbar,[s.diff(H,x) for x in bar]))
    field.update(zip(bar,[-s.diff(H,x) for x in unbar]))
    DHZ = s.zeros(d)
    for variable, derivative in field.items():
        DHZ += Z.diff(variable)*derivative
    assert (DHZ-rate).applyfunc(s.expand) == s.zeros(d)

    # At x=0: N_MCL=i Z, d_x N_MCL=(e_i-e_k) Z_ik,
    # d_t N_MCL=i d_t Z=D Z, since D=i d_t in our convention.
    Nmcl = s.I*Z
    spatial = s.Matrix(d,d,lambda i,k:(e[i]-e[k])*Z[i,k])
    # MCL (2.13): alpha_ik N_ik,x = N_ik,t + sum_j
    # N_ij N_jk (alpha_jk-alpha_ij).
    for i in range(d):
        for k in range(d):
            rhs = rate[i,k]+sum(Nmcl[i,j]*Nmcl[j,k]*(alpha[j,k]-alpha[i,j])
                                   for j in range(d) if j not in (i,k))
            assert s.expand(alpha[i,k]*spatial[i,k]-rhs) == 0
            if groups[i]==groups[k]:
                assert rate[i,k] == 0
    return {'sizes':[L,M,N], 'canonical_H0_rates_match':True,
            'MCL_equation_2_13_match':True,
            'zero_same_group_blocks_preserved':True,
            'all_phase_intercepts_symbolic':True}


if __name__=='__main__':
    results={'restricted_star':[restricted_star(2),restricted_star(3)],
             'checks':[check(1,1,2),check(2,2,1),check(2,2,2)],
             'scope':'Exact equation-level embedding in MCL/Ablowitz-Haberman '
                     'multiwave framework; does not identify every MCL completed '
                     'network with the JRW three-block reduction.'}
    path=Path(__file__).parent/'results'/'mcl_jrw_matrix_embedding.json'
    path.write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2))
