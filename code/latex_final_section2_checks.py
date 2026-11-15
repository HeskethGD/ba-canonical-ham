"""Check Section 2 in its L,M,N mode notation, including the JRW map.

Partner-symbol suffixes denote independent Hamiltonian partners.
Run from the repository root with .venv/bin/python.
"""
import json
import sympy as s


def check(L, M, N):
    def matrix(name, rows, cols):
        return s.Matrix(rows, cols, lambda i, j: s.Symbol(f'{name}_{i}_{j}'))

    u, us = matrix('u', L, M), matrix('us', L, M)
    v, vs = matrix('v', M, N), matrix('vs', M, N)
    w, ws = matrix('w', N, L), matrix('ws', N, L)
    q, p = list(u) + list(v) + list(w), list(us) + list(vs) + list(ws)
    bL = s.symbols(f'bL0:{L}')
    bM = s.symbols(f'bM0:{M}')
    c = s.symbols(f'c0:{N}')
    diag = lambda entries: s.diag(*entries)
    D = diag([1] * L + [-1] * M + [0] * N)
    X = s.BlockMatrix([
        [diag([b / 2 for b in bL]), s.sqrt(2) * u, ws.T],
        [s.sqrt(2) * us.T, diag([b / 2 for b in bM]), v],
        [-w, -vs.T, diag(c)],
    ]).as_explicit()
    H = -sum((bL[i] / 2 - c[n]) * w[n, i] * ws[n, i]
             for i in range(L) for n in range(N))
    H -= sum((bM[j] / 2 - c[n]) * v[j, n] * vs[j, n]
             for j in range(M) for n in range(N))
    H += s.sqrt(2) * sum(
        u[i, j] * v[j, n] * w[n, i]
        + us[i, j] * vs[j, n] * ws[n, i]
        for i in range(L) for j in range(M) for n in range(N))

    def pb(f, h):
        return s.expand(sum(s.diff(f, a) * s.diff(h, b)
                            - s.diff(f, b) * s.diff(h, a)
                            for a, b in zip(q, p)))

    def zero(A):
        return all(s.expand(e) == 0 for e in A)

    nu, mu = s.symbols('nu mu')
    pencil = X + nu * D
    T = s.BlockMatrix([
        [s.zeros(L), s.zeros(L, M), ws.T],
        [s.zeros(M, L), s.zeros(M), -v],
        [-w, vs.T, -nu * s.eye(N)],
    ]).as_explicit()
    dotX = X.applyfunc(lambda e: pb(e, H))
    assert zero(dotX - pencil * T + T * pencil)

    QL = [sum(w[n, i] * ws[n, i] for n in range(N))
          - sum(u[i, j] * us[i, j] for j in range(M)) for i in range(L)]
    QM = [sum(v[j, n] * vs[j, n] for n in range(N))
          - sum(u[i, j] * us[i, j] for i in range(L)) for j in range(M)]
    QN = [sum(w[n, i] * ws[n, i] for i in range(L))
          - sum(v[j, n] * vs[j, n] for j in range(M)) for n in range(N)]
    assert s.expand(sum(QL) - sum(QM) - sum(QN)) == 0
    assert all(pb(charge, H) == 0 for charge in QL + QM + QN)
    Htrace = -s.trace(X ** 3) / 3
    Htrace -= sum(bL[i] * QL[i] for i in range(L))
    Htrace -= sum(bM[j] * QM[j] for j in range(M))
    Htrace += sum(b ** 3 for b in bL + bM) / 24
    Htrace += sum(a ** 3 for a in c) / 3
    assert s.expand(H - Htrace) == 0

    # Symmetric endpoint-sum convention used in the revised manuscript.
    betaL = s.symbols(f'betaL0:{L}')
    betaM = s.symbols(f'betaM0:{M}')
    betaN = s.symbols(f'betaN0:{N}')
    phases = {bL[i]: -4 * betaL[i] for i in range(L)}
    phases.update({bM[j]: -4 * betaM[j] for j in range(M)})
    phases.update({c[n]: betaN[n] for n in range(N)})
    H0 = H.subs(phases)
    Hsym = s.sqrt(2) * sum(
        u[i, j] * v[j, n] * w[n, i]
        + us[i, j] * vs[j, n] * ws[n, i]
        for i in range(L) for j in range(M) for n in range(N))
    Hsym += sum((betaL[i] + betaM[j]) * u[i, j] * us[i, j]
                for i in range(L) for j in range(M))
    Hsym += sum((betaM[j] + betaN[n]) * v[j, n] * vs[j, n]
                for j in range(M) for n in range(N))
    Hsym += sum((betaN[n] + betaL[i]) * w[n, i] * ws[n, i]
                for n in range(N) for i in range(L))
    gauge_H = -sum(betaL[i] * QL[i] for i in range(L))
    gauge_H -= sum(betaM[j] * QM[j] for j in range(M))
    assert s.expand(Hsym - H0 - gauge_H) == 0
    A = diag(list(betaL) + [-b for b in betaM] + [0] * N)
    Xsym = X.subs(phases)
    Lsym = Xsym + nu * D
    Tsym = T - A
    assert zero(Xsym.applyfunc(lambda e: pb(e, Hsym)) - Lsym * Tsym + Tsym * Lsym)
    BL, BM, BN = diag(betaL), diag(betaM), diag(betaN)
    matrix_rates = [
        (u, BL*u + u*BM + s.sqrt(2)*ws.T*vs.T),
        (us.T, -BM*us.T - us.T*BL - s.sqrt(2)*v*w),
        (v, BM*v + v*BN + s.sqrt(2)*us.T*ws.T),
        (vs.T, -BN*vs.T - vs.T*BM - s.sqrt(2)*w*u),
        (w, BN*w + w*BL + s.sqrt(2)*vs.T*us.T),
        (ws.T, -BL*ws.T - ws.T*BN - s.sqrt(2)*u*v),
    ]
    for modes, target in matrix_rates:
        assert zero(modes.applyfunc(lambda e: pb(e, Hsym)) - target)
    assert zero(Xsym.applyfunc(lambda e: pb(e, gauge_H)) - A * Xsym + Xsym * A)
    symmetric_trace_H = -s.trace(Xsym ** 3) / 3
    symmetric_trace_H += 3 * sum(betaL[i] * QL[i] for i in range(L))
    symmetric_trace_H += 3 * sum(betaM[j] * QM[j] for j in range(M))
    symmetric_trace_H -= s.Rational(8, 3) * sum(b ** 3 for b in betaL + betaM)
    symmetric_trace_H += sum(b ** 3 for b in betaN) / 3
    assert s.expand(Hsym - symmetric_trace_H) == 0
    all_trace_H = -s.trace(Xsym ** 3) / 3
    all_trace_H += 3 * s.trace(BL * (ws.T*w - u*us.T))
    all_trace_H += 3 * s.trace(BM * (v*vs.T - us.T*u))
    all_trace_H -= s.Rational(8, 3) * (s.trace(BL**3) + s.trace(BM**3))
    all_trace_H += s.trace(BN**3) / 3
    assert s.expand(Hsym - all_trace_H) == 0
    rates = [betaL[i] + betaM[j] for i in range(L) for j in range(M)]
    rates += [-betaM[j] for j in range(M) for n in range(N)]
    rates += [-betaL[i] for n in range(N) for i in range(L)]
    for unstarred, starred, rate in zip(q, p, rates):
        assert s.expand(pb(unstarred, Hsym) - pb(unstarred, H0) - rate * unstarred) == 0
        assert s.expand(pb(starred, Hsym) - pb(starred, H0) + rate * starred) == 0

    # Full ordered entry-bracket check for the new orientation of w.
    for a in range(X.rows):
        for b in range(X.cols):
            for k in range(X.rows):
                for r in range(X.cols):
                    expected = D[a, a] - D[b, b] if a == r and b == k else 0
                    assert pb(X[a, b], X[k, r]) == expected

    # JRW complex coordinates and the constant diagonal extension.
    theta = diag([1] * L + [-1] * M)
    tau = w.T.col_join(vs)
    dual = ws.row_join(v.T)
    varpi = s.BlockMatrix([
        [s.zeros(L), -s.I * s.sqrt(2) * us],
        [-s.I * s.sqrt(2) * u.T, s.zeros(M)],
    ]).as_explicit()
    Z = diag([b / (2 * s.I) for b in bL + bM])
    E = diag([a / s.I for a in c])
    HJ = s.I * s.trace((varpi - Z) * tau * dual + tau * E * dual)
    assert s.expand(HJ - H) == 0
    J = s.BlockMatrix([[Z + varpi, tau], [dual, E]]).as_explicit()
    J += mu * diag([-1] * L + [1] * M + [0] * N)
    G = diag([1] * (L + M) + [s.I] * N)
    assert zero(s.I * J - G * pencil.subs(nu, -s.I * mu).T * G.inv())
    Y = s.BlockMatrix([
        [s.zeros(L + M), s.I * theta * tau],
        [s.I * dual * theta, s.I * mu * s.eye(N)],
    ]).as_explicit()
    assert zero(J.applyfunc(lambda e: pb(e, H)) - Y * J + J * Y)

    if N == 1:
        witness = dict.fromkeys(list(u) + list(us), 0)
        witness.update(dict.fromkeys(list(v) + list(vs) + list(w) + list(ws), 1))
        lam = s.Symbol('lambda')
        for i in range(L):
            for j in range(M):
                # Determinant differentiation at the witness avoids expanding F.
                A = (lam * s.eye(X.rows) - pencil).subs({
                    lam: (bL[i] + bM[j]) / 4,
                    nu: (bM[j] - bL[i]) / 4,
                }).subs(witness)
                expected = s.sqrt(2)
                expected *= s.prod((bL[i] - bL[r]) / 2 for r in range(L) if r != i)
                expected *= s.prod((bM[j] - bM[r]) / 2 for r in range(M) if r != j)
                for a in range(L):
                    for k in range(M):
                        # dA/du_(a,k) has its only entry -sqrt(2) at (a,L+k).
                        derivative = -s.sqrt(2) * A.cofactor(a, L + k)
                        target = expected if (a, k) == (i, j) else 0
                        assert s.simplify(s.expand(derivative - target)) == 0, (
                            L, M, N, i, j, a, k, s.factor(derivative), s.factor(target))
                        # The displayed symmetric minor is this same identity
                        # after b_old = -4 beta, with no time or curve change.
                        symmetric_target = s.sqrt(2) if (a, k) == (i, j) else s.S.Zero
                        symmetric_target *= s.prod(2 * (betaL[r] - betaL[i])
                                                   for r in range(L) if r != i)
                        symmetric_target *= s.prod(2 * (betaM[r] - betaM[j])
                                                   for r in range(M) if r != j)
                        assert s.simplify(s.expand(derivative.subs(phases) - symmetric_target)) == 0
    return {'sizes': [L, M, N], 'canonical_Lax_equations': True,
            'charges_and_trace_identity': True, 'matrix_bracket': True,
            'JRW_Hamiltonian_pencil_and_time_matrix': True,
            'symmetric_Hamiltonian_and_Lax_pair': True,
            'six_matrix_equations': True,
            'symmetric_trace_identity': True,
            'all_trace_Hamiltonian_identity': True,
            'canonical_time_gauge_and_diagonal_similarity': True,
            'independence_minor': 'verified' if N == 1 else 'not tested for N > 1'}


if __name__ == '__main__':
    print(json.dumps([check(*sizes) for sizes in
                      [(1, 1, 1), (2, 2, 1), (2, 2, 2)]], indent=2))
