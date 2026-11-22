# Independent check of the al-function example in Section 5.4 (branch-point evaluations).
# Run from the repository root: .venv/bin/python code/al_function_checks.py
# No sigma functions are used. On the al family every mode is algebraic in the moving divisor:
#   Psi(B_k)^2 = (-1)^g P(z,a_k),  P(z,x) = prod_i (x - x_i),
# and the z_g-flow of the divisor is the Dubrovin flow x_i' = y_i / prod_{j != i} (x_i - x_j).
import sympy as sp, mpmath as mp
mp.mp.dps = 40
s2 = sp.sqrt(2)

def al_state(g, branch, alpha, div_x):
    """Modes of (5.x) at the 2g+1 branch points, groups L=g, M=g+1, and their z_g-derivatives."""
    L, M = g, g + 1; K = L + M
    a = [sp.Rational(v) for v in branch]
    X = sp.symbols('x1:%d' % (g + 1)); Y = sp.symbols('y1:%d' % (g + 1))
    f = lambda t: 4 * sp.prod([t - ak for ak in a])
    dX = [Y[i] / sp.prod([X[i] - X[j] for j in range(g) if j != i]) for i in range(g)]
    dY = [sp.diff(f(X[i]), X[i]) * dX[i] / (2 * Y[i]) for i in range(g)]
    D = lambda e: sum(sp.diff(e, X[i]) * dX[i] + sp.diff(e, Y[i]) * dY[i] for i in range(g))
    RL = lambda t: sp.prod([t - a[i] for i in range(L)])
    RM = lambda t: sp.prod([t - a[j] for j in range(L, K)])
    t = sp.symbols('t')
    kap = [1 / (2 * sp.diff(RL(t), t).subs(t, a[i])) for i in range(L)] + \
          [(a[j] - alpha) / (2 * sp.diff(RM(t), t).subs(t, a[j])) for j in range(L, K)]
    # moments (5.x): sum kappa a^r = 0 (r <= g-2), 1/2 (r = g-1), in both groups
    for grp in (range(L), range(L, K)):
        assert [sum(kap[k] * a[k] ** r for k in grp) for r in range(g)] == [0] * (g - 1) + [sp.Rational(1, 2)]
    eps = (-1) ** g; e = sp.sqrt(eps)                      # Psi(B_k) = e * phi_k, e^2 = (-1)^g
    phi = [sp.sqrt(sp.prod([a[k] - X[i] for i in range(g)])) for k in range(K)]
    xi = [e * p for p in phi]                              # c_k = 1
    xib = [eps * kap[k] * e * phi[k] for k in range(K)]    # sheet partner coincides at a branch point
    W = lambda p, q: p * D(q) - D(p) * q
    C = [[s2 / (a[j] - a[i]) * W(xib[i], xib[j]) for j in range(L, K)] for i in range(L)]
    Cb = [[-s2 / (a[j] - a[i]) * W(xi[i], xi[j]) for j in range(L, K)] for i in range(L)]
    state = xi + xib + [c for row in C for c in row] + [c for row in Cb for c in row]
    dstate = [D(s) for s in state]
    pot = [D(D(p)) / p - a[k] for k, p in enumerate(phi)]  # common Schroedinger potential
    vals = {X[i]: div_x[i] for i in range(g)}
    vals.update({Y[i]: mp.sqrt(4 * mp.fprod([div_x[i] - mp.mpf(str(v)) for v in branch])) for i in range(g)})
    num = lambda es: [mp.mpmathify(sp.N(ex.subs(vals), 40)) for ex in es]
    return L, M, num(state), num(dstate), num(pot), [float(k) for k in kap]

def rhs_prime(z, L, M):
    """z_g-derivative (= i d/dt) required by Hamilton's equations of H0 with all b_k = 0."""
    K = L + M; r2 = mp.sqrt(2)
    x, xb = z[:K], z[K:2 * K]
    C = [[z[2 * K + i * M + j] for j in range(M)] for i in range(L)]
    Cb = [[z[2 * K + L * M + i * M + j] for j in range(M)] for i in range(L)]
    dx = [r2 * sum(Cb[i][j] * xb[L + j] for j in range(M)) for i in range(L)] + \
         [r2 * sum(Cb[i][j] * xb[i] for i in range(L)) for j in range(M)]
    dxb = [-r2 * sum(C[i][j] * x[L + j] for j in range(M)) for i in range(L)] + \
          [-r2 * sum(C[i][j] * x[i] for i in range(L)) for j in range(M)]
    dC = [r2 * xb[i] * xb[L + j] for i in range(L) for j in range(M)]
    dCb = [-r2 * x[i] * x[L + j] for i in range(L) for j in range(M)]
    return dx + dxb + dC + dCb

def discriminant_roots(z, L, M):
    lam, nu = sp.symbols('lam nu'); K = L + M; d = K + 1
    cv = lambda c: sp.Rational(float(c.real)) + sp.I * sp.Rational(float(c.imag))
    Mx = sp.zeros(d, d)
    for i in range(L):
        Mx[i, i] = nu; Mx[i, d - 1] = cv(z[K + i]); Mx[d - 1, i] = -cv(z[i])
        for j in range(M):
            Mx[i, L + j] = cv(mp.sqrt(2) * z[2 * K + i * M + j]); Mx[L + j, i] = cv(mp.sqrt(2) * z[2 * K + L * M + i * M + j])
    for j in range(M):
        Mx[L + j, L + j] = -nu; Mx[L + j, d - 1] = cv(z[L + j]); Mx[d - 1, L + j] = -cv(z[K + L + j])
    F = (lam * sp.eye(d) - Mx).det(method='berkowitz')
    Dp = sp.Poly(sp.resultant(F, sp.diff(F, lam), lam), nu, extension=sp.I)
    mp.mp.dps = 80
    co = [mp.mpc(sp.N(sp.re(c), 80), sp.N(sp.im(c), 80)) for c in Dp.all_coeffs()]
    roots = mp.polyroots(co, maxsteps=8000, extraprec=2000)
    mp.mp.dps = 40
    sep = [min(abs(roots[i] - roots[j]) for j in range(len(roots)) if j != i) for i in range(len(roots))]
    return Dp.degree(), sep

cases = [(1, [-1, 0, 2], sp.Rational(1, 3), [mp.mpc('0.3', '0.4')]),
         (2, [-2, -1, 0, 1, 3], sp.Rational(3, 8), [mp.mpc('0.3', '0.4'), mp.mpc('-1.3', '0.2')]),
         (3, [-3, -2, -1, 0, 1, 2, 4], sp.Rational(3, 8), [mp.mpc('0.3', '0.4'), mp.mpc('-1.3', '0.2'), mp.mpc('2.6', '-0.7')])]
for g, branch, alpha, div_x in cases:
    L, M, z, dz, pot, kap = al_state(g, branch, alpha, div_x)
    need = rhs_prime(z, L, M)
    res = max(abs(p - q) / (1 + abs(p)) for p, q in zip(dz, need))
    spread = max(abs(p - pot[0]) for p in pot)
    assert mp.isfinite(res) and res < mp.mpf('1e-25'), (g, 'Hamilton residual', res)
    assert mp.isfinite(spread) and spread < mp.mpf('1e-25'), (g, 'potential spread', spread)
    print('g=%d, L=%d, M=%d: Hamilton residual with all phases zero %.1e | spread of Schroedinger potential over the %d branch points %.1e'
          % (g, L, M, res, 2 * g + 1, spread))
    if g <= 2:
        deg, sep = discriminant_roots(z, L, M)
        assert deg == 2*g-2+2*(L+M+1), (g, 'discriminant degree', deg)
        assert min(sep) > mp.mpf('1e-8'), (g, 'unresolved discriminant roots', min(sep))
        print('     discriminant of F in lambda: degree %d (a genus-%d curve with no affine singularity predicts %d), smallest root separation %.2e'
              % (deg, g, 2 * g - 2 + 2 * (L + M + 1), float(min(sep))))
