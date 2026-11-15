# Independent checks of the higher-flow theorem (Section 4.7: Lemma 4.7, Theorem 4.8).
# Run from the repository root: .venv/bin/python code/higher_flow_evolution_diagnostic.py
# No sigma or theta functions are used. The theorem predicts, for the flow of H_{r,nu0}:
#   (a) the Abel direction sum_alpha v_alpha U_{n_alpha}, and (b) constant mode phases b_k^{(r,nu0)},
# both from a meromorphic function f on the curve. These predictions are computed here from the
# curve alone and compared with direct integration of the matrix flow  i X' = [D, (X+nu0 D)^(r-1)].
from pathlib import Path
import numpy as np
rng = np.random.default_rng(7)

def flow(X, D, field, T, n=4000):
    """RK4 for i X' = field(X) over (possibly complex) time T."""
    h = T / n
    f = lambda Y: -1j * field(Y)
    for _ in range(n):
        k1 = f(X); k2 = f(X + h/2*k1); k3 = f(X + h/2*k2); k4 = f(X + h*k3)
        X = X + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return X

def higher_field(D, r, nu0):
    return lambda X: D @ np.linalg.matrix_power(X + nu0*D, r-1) - np.linalg.matrix_power(X + nu0*D, r-1) @ D

def h0_field(L, M):
    K = L + M
    def field(X):                                   # [L0, T0] of (4.x), evaluated at nu = 0
        T0 = np.zeros_like(X)
        T0[:L, K] = X[:L, K]; T0[L:K, K] = -X[L:K, K]; T0[K, :L] = X[K, :L]; T0[K, L:K] = -X[K, L:K]
        return X @ T0 - T0 @ X
    return field

def gauge(X, L, M, theta):
    """c_k -> c_k exp(-i theta_k): xi_k -> exp(-i theta_k) xi_k, partners inverse, pair modes exp(+i(theta_i+theta_j))."""
    s = np.concatenate([np.exp(1j*theta[:L]), np.exp(-1j*theta[L:]), [1]])
    return (s[:, None] * X) / s[None, :]

def scaled(A, B): return (np.abs(A - B) / (1 + np.maximum(np.abs(A), np.abs(B)))).max()

def laurent(g, n_lo, n_hi, centre=0, radius=0.05, N=256):
    """Laurent coefficients c_n, n_lo <= n <= n_hi, of g about centre by the trapezoid rule on a circle."""
    th = 2*np.pi*np.arange(N)/N; w = radius*np.exp(1j*th)
    vals = np.array([g(centre + wi) for wi in w])
    return {n: (vals * w**(-n)).mean() for n in range(n_lo, n_hi + 1)}

# ---------------------------------------------------------------- genus one, generic data, L = M = 1
print('Genus one: generic initial data of the 3x3 pencil (smooth plane cubic)')
L = M = 1; K = 2; D = np.diag([1., -1., 0.]).astype(complex)
X0 = 0.6*(rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3)))
dvals = [1, -1, 0]
def branch(k):
    """lambda as a function of w = 1/nu on the infinity branch with slope dvals[k]."""
    def lam(w):
        ev = np.linalg.eigvals(D + w*X0)
        return ev[np.argmin(np.abs(ev - dvals[k]))] / w
    return lam
for r, nu0 in [(3, 0.0), (4, 0.0), (4, 0.4 - 0.3j), (5, 0.2j)]:
    m = r - 2
    mons = [(a, b) for a in range(m + 1) for b in range(m + 1 - a) if a < 3]      # lambda^a nu^b modulo the cubic
    q = lambda lam, nu: lam**(r-1) / (nu - nu0)
    rows, rhs = [], []
    co = {}
    for k in range(3):
        lam = branch(k)
        cm = [laurent(lambda w, a=a, b=b: lam(w)**a * w**(-b), -m, 0) for a, b in mons]
        cq = laurent(lambda w: q(lam(w), 1/w), -m, 0)
        co[k] = (cm, cq)
        if k < 2:                                   # outer points: cancel the principal part of q
            for n in range(-m, 0):
                rows.append([c[n] for c in cm]); rhs.append(-cq[n])
        else:                                       # centre: only a simple pole, and zero constant term
            for n in list(range(-m, -1)) + [0]:
                rows.append([c[n] for c in cm]); rhs.append(0)
    coef = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)[0]
    v1 = sum(c*cm[-1] for c, cm in zip(coef, co[2][0]))          # f = v1*nu + O(1/nu) = -v1/s + O(s)
    b = np.array([-dvals[k]*(sum(c*cm[0] for c, cm in zip(coef, co[k][0])) + co[k][1][0]) for k in range(2)])
    t = 0.25
    direct = flow(X0, D, higher_field(D, r, nu0), t)
    predicted = gauge(flow(X0, D, h0_field(L, M), v1*t), L, M, t*b)
    print('  r=%d nu0=%s: v1=%.6f%+.6fi | flow of H_r vs (H0 flow for time v1*t, then the predicted gauge): %.1e | without the gauge: %.1e'
          % (r, nu0, v1.real, v1.imag, scaled(direct, predicted), scaled(direct, flow(X0, D, h0_field(L, M), v1*t))))

# ---------------------------------------------------------------- genus two, stored hyperelliptic example
print('Genus two: stored hyperelliptic realization, L = M = 2 (plane image has two nodes)')
d = np.load(Path(__file__).resolve().parent/'results/hyperelliptic_sigma_ba.npz'); bst = d['phases']; z0 = d['BA'][0]
L = M = 2; K = 4; s2 = np.sqrt(2)
x_, xb_, C_, Cb_ = z0[:K], z0[K:2*K], z0[2*K:2*K+L*M].reshape(L, M), z0[2*K+L*M:].reshape(L, M)
X0 = np.vstack([np.hstack([np.diag(-bst[:L]), s2*C_, xb_[:L, None]]),
                np.hstack([s2*Cb_.T, np.diag(-bst[L:]), x_[L:, None]]),
                np.hstack([-x_[:L][None, :], -xb_[L:][None, :], [[0]]])]).astype(complex)
D = np.diag([1., 1., -1., -1., 0.]).astype(complex)
a = np.array([0.1, 2.3, 0.5, 3.0]); f = lambda t: 4*t*(t*t - 1)*(t*t - 4)
eta = np.sqrt(f(a)) * np.array([1, 1, -1, -1])                    # ordinates of the pole points of nu
raw = np.array([1/(2*(a[0]-a[1])), 1/(2*(a[1]-a[0])), 1/(2*(a[2]-a[3])), 1/(2*(a[3]-a[2]))])
hL = lambda x, y: sum(raw[j]*(y + eta[j])/(x - a[j]) for j in (0, 1))
hM = lambda x, y: sum(raw[j]*(y + eta[j])/(x - a[j]) for j in (2, 3))
nu_ = lambda x, y: -(hL(x, y) + hM(x, y))/2
lam_ = lambda x, y: (hM(x, y) - hL(x, y))/2
dk = [1, 1, -1, -1]
# H0 check of the matrix form against the stored trajectory
print('  H0 matrix flow vs stored BA values at the fourth stored time: %.1e' % (
      scaled(flow(X0, D, h0_field(L, M), 6.0*4/40, 8000)[K, :K], np.concatenate([-d['BA'][4][:L], -d['BA'][4][K+L:2*K]]))))

def normalization(r, nu0):
    """f = (A(x) + B(x) y) / R(x)^m with deg A, deg B <= 4m-1; returns v1, v3 and the phases b_k."""
    m = r - 2; R = lambda x: np.prod([x - ak for ak in a], axis=0)
    q = lambda x, y: lam_(x, y)**(r-1) / (nu_(x, y) - nu0)
    nb = 4*m
    def local(k, sign):
        yk = lambda x: sign*eta[k]*np.sqrt(f(x)/f(a[k]))
        basisA = [laurent(lambda x, p=p: x**p / R(x)**m, -m, 0, a[k], 0.03) for p in range(nb)]
        basisB = [laurent(lambda x, p=p: x**p * yk(x) / R(x)**m, -m, 0, a[k], 0.03) for p in range(nb)]
        return basisA + basisB, laurent(lambda x: q(x, yk(x)), -m, 0, a[k], 0.03)
    rows, rhs, keep = [], [], {}
    for k in range(K):
        bas, cq = local(k, +1); keep[k] = (bas, cq)
        for n in range(-m, 0):
            rows.append([c[n] for c in bas]); rhs.append(-cq[n])
        bas2, _ = local(k, -1)                                   # sheet partner: f must be regular there
        for n in range(-m, 0):
            rows.append([c[n] for c in bas2]); rhs.append(0)
    coef = np.linalg.solve(np.array(rows), np.array(rhs))
    B = coef[nb:]
    beta1 = B[nb-1]; beta2 = B[nb-2] + m*a.sum()*B[nb-1]         # B/R^m = beta1/x + beta2/x^2 + ...
    v3, v1 = -2*beta1, -2*beta2                                  # y = 2 s^-5 (1 + O(s^4)), x = s^-2
    b = np.array([-dk[k]*(sum(c*cm[0] for c, cm in zip(coef, keep[k][0])) + keep[k][1][0]) for k in range(K)])
    return v1, v3, b

def divisor_velocity(X, field):
    """Velocities of the two divisor abscissae, from the mode powers xi_k xibar_k = kappa_k P(z, a_k)."""
    def P_coeffs(Y):
        I = np.array([-Y[k, K]*Y[K, k] for k in range(K)]) / raw
        return np.linalg.lstsq(np.vstack([a, np.ones(K)]).T, I - a**2, rcond=None)
    (c1, c0), res = P_coeffs(X)[0], P_coeffs(X)[1]
    e = 1e-6; Xd = -1j*field(X)
    dc = (P_coeffs(X + e*Xd)[0] - P_coeffs(X - e*Xd)[0]) / (2*e)
    roots = np.roots([1, c1, c0])
    return roots, np.array([-(dc[0]*x + dc[1]) / (2*x + c1) for x in roots]), res

roots, vel0, res = divisor_velocity(X0, h0_field(L, M))
print('  mode powers lie on one monic quadratic in a_k (four points, two unknowns): residual %.1e' % (res[0] if len(res) else 0))
pred = {}
for r, nu0 in [(3, 0.0), (4, 0.0), (4, 0.35), (5, 0.0)]:
    v1, v3, b = normalization(r, nu0); pred[(r, nu0)] = (v1, v3, b)
    _, vel, _ = divisor_velocity(X0, higher_field(D, r, nu0))
    ratio = vel / vel0
    expect = np.array([v1 - roots[1]*v3, v1 - roots[0]*v3])
    print('  r=%d nu0=%.2f: v1=%.5f%+.5fi v3=%.5f%+.5fi | divisor velocity relative to H0, measured vs predicted v1 - x_other*v3: %.1e'
          % (r, nu0, v1.real, v1.imag, v3.real, v3.imag, np.abs(ratio - expect).max() / np.abs(expect).max()))
# consistency between flows, which also tests the phases
t = 0.02
v1, v3, b = pred[(3, 0.0)]                                        # the third trace moves along U_1 only
direct = flow(X0, D, higher_field(D, 3, 0.0), t, 6000); comp = flow(X0, D, h0_field(L, M), v1*t, 6000)
print('  flow (3, 0.0) vs H0 for time v1*t with gauge t*b: %.1e | without the gauge: %.1e'
      % (scaled(direct, gauge(comp, L, M, t*b)), scaled(direct, comp)))
for Bk, Ck in [((4, 0.0), (4, 0.35)), ((4, 0.0), (5, 0.0)), ((4, 0.35), (5, 0.0))]:
    (v1B, v3B, bB), (v1C, v3C, bC) = pred[Bk], pred[Ck]
    beta = v3C / v3B; alpha = v1C - beta*v1B                      # direction of C = alpha*U_1 + beta*(direction of B)
    direct = flow(X0, D, higher_field(D, *Ck), t, 6000)
    comp = flow(flow(X0, D, higher_field(D, *Bk), beta*t, 20000), D, h0_field(L, M), alpha*t, 20000)
    print('  flow %s vs [H0 for %.4f t, flow %s for %.4f t] with gauge t(b^C - beta b^B): %.1e | without the gauge: %.1e'
          % (Ck, alpha.real, Bk, beta.real, scaled(direct, gauge(comp, L, M, t*(bC - beta*bB))), scaled(direct, comp)))

# ---------------------------------------------------------------- fifth trace again, normalization at 40 digits
# The double-precision linear solve for r = 5 (24 unknowns, monomials up to x^11) limits the phases to
# about 1e-8. Repeating the interpolation in mpmath removes that limit.
import mpmath as mp
mp.mp.dps = 40
def normalization_mp(r, nu0, N=160, radius=mp.mpf('0.03')):
    m = r - 2; nb = 4*m
    A_ = [mp.mpf(str(v)) for v in a]; fm = lambda t: 4*t*(t*t - 1)*(t*t - 4)
    sg = [1, 1, -1, -1]; ETA = [sg[k]*mp.sqrt(fm(A_[k])) for k in range(K)]
    RAW = [1/(2*(A_[0]-A_[1])), 1/(2*(A_[1]-A_[0])), 1/(2*(A_[2]-A_[3])), 1/(2*(A_[3]-A_[2]))]
    hLm = lambda x, y: sum(RAW[j]*(y + ETA[j])/(x - A_[j]) for j in (0, 1))
    hMm = lambda x, y: sum(RAW[j]*(y + ETA[j])/(x - A_[j]) for j in (2, 3))
    Rm = lambda x: mp.fprod([x - ak for ak in A_])
    q = lambda x, y: ((hMm(x, y) - hLm(x, y))/2)**(r-1) / (-(hLm(x, y) + hMm(x, y))/2 - nu0)
    ws = [radius*mp.expjpi(mp.mpf(2*n)/N) for n in range(N)]
    def lau(vals, n): return mp.fsum(v*w**(-n) for v, w in zip(vals, ws))/N
    rows, rhs, keep = [], [], {}
    for k in range(K):
        for sign in (1, -1):
            xs = [A_[k] + w for w in ws]
            ys = [sign*ETA[k]*mp.sqrt(fm(x)/fm(A_[k])) for x in xs]
            Rs = [Rm(x)**m for x in xs]
            cols = [[x**p/R for x, R in zip(xs, Rs)] for p in range(nb)] + [[x**p*y/R for x, y, R in zip(xs, ys, Rs)] for p in range(nb)]
            qv = [q(x, y) for x, y in zip(xs, ys)]
            for n in range(-m, 0):
                rows.append([lau(c, n) for c in cols]); rhs.append(-lau(qv, n) if sign == 1 else 0)
            if sign == 1: keep[k] = ([lau(c, 0) for c in cols], lau(qv, 0))
    coef = mp.lu_solve(mp.matrix(rows), mp.matrix(rhs))
    B = [coef[nb + p] for p in range(nb)]
    v3 = -2*B[nb-1]; v1 = -2*(B[nb-2] + m*sum(A_)*B[nb-1])
    b = [-dk[k]*(sum(coef[p]*keep[k][0][p] for p in range(2*nb)) + keep[k][1]) for k in range(K)]
    return complex(v1), complex(v3), np.array([complex(x) for x in b])
t = 0.02
predmp = {key: normalization_mp(*key) for key in [(4, 0.0), (5, 0.0)]}
for key in predmp:
    print('  %s: change in (v1, v3, phases) between double and 40-digit interpolation: %.1e, %.1e, %.1e'
          % (key, abs(predmp[key][0]-pred[key][0]), abs(predmp[key][1]-pred[key][1]), np.abs(predmp[key][2]-pred[key][2]).max()))
(v1B, v3B, bB), (v1C, v3C, bC) = predmp[(4, 0.0)], predmp[(5, 0.0)]
beta = v3C / v3B; alpha = v1C - beta*v1B
direct = flow(X0, D, higher_field(D, 5, 0.0), t, 20000)
comp = flow(flow(X0, D, higher_field(D, 4, 0.0), beta*t, 40000), D, h0_field(L, M), alpha*t, 40000)
print('  fifth trace vs [H0, fourth trace] with 40-digit directions and phases: %.1e' % scaled(direct, gauge(comp, L, M, t*(bC - beta*bB))))
