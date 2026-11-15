"""Separate RK4 truncation and precision effects using recorded BA data.

Uses the saved initial modes and phase constants; it does not recompute
sigma values or alter the published validation records.
"""
from pathlib import Path
import json
import math
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent


def discrepancy(a, b):
    return float(np.max(abs(a-b)/(1+np.maximum(abs(a), abs(b)))))


def integrate(initial, phases, L, M, duration, steps, samples, dtype):
    K = L+M
    P = L*M
    phases = phases.astype(dtype)
    coupling = dtype(math.sqrt(2))  # Same fixed coefficient as the main run.
    h = np.longdouble(duration)/steps if dtype == np.clongdouble else duration/steps
    z = initial.astype(dtype)

    def rhs(v):
        u, partner = v[:K], v[K:2*K]
        C = v[2*K:2*K+P].reshape(L, M)
        Cbar = v[2*K+P:].reshape(L, M)
        du, dp = phases*u, -phases*partner
        du[:L] += coupling*Cbar@partner[L:]
        du[L:] += coupling*Cbar.T@partner[:L]
        dp[:L] -= coupling*C@u[L:]
        dp[L:] -= coupling*C.T@u[:L]
        return -1j*np.concatenate((du, dp,
            (coupling*partner[:L, None]*partner[None, L:]).ravel(),
            (-coupling*u[:L, None]*u[None, L:]).ravel()))

    trajectory = [z.copy()]
    stride = steps//(samples-1)
    assert steps % (samples-1) == 0
    for step in range(1, steps+1):
        k1 = rhs(z)
        k2 = rhs(z+h*k1/2)
        k3 = rhs(z+h*k2/2)
        k4 = rhs(z+h*k3)
        z += h*(k1+2*k2+2*k3+k4)/6
        if step % stride == 0:
            trajectory.append(z.copy())
    return np.array(trajectory)


def main():
    data = np.load(ROOT/'results/trigonal_sigma_ba.npz')
    reference = data['BA']
    duration = float(data['analytic_t'][-1])
    sizes = [1200, 2400, 6000, 12000, 24000, 48000]
    results = {}
    trajectories = {}
    for label, dtype in [('double', np.complex128), ('extended', np.clongdouble)]:
        rows = []
        for steps in sizes:
            trajectory = integrate(reference[0], data['phases'], 3, 3,
                                   duration, steps, len(reference), dtype)
            trajectories[label, steps] = trajectory
            row = dict(steps=steps, BA_discrepancy=discrepancy(trajectory, reference))
            if (label, steps//2) in trajectories:
                row['step_halving_difference'] = discrepancy(
                    trajectory, trajectories[label, steps//2])
            rows.append(row)
            print(label, row, flush=True)
        results[label] = rows
    results['precision_comparison'] = [dict(steps=n, discrepancy=discrepancy(
        trajectories['double', n], trajectories['extended', n])) for n in sizes]
    results['extended_mantissa_bits'] = np.finfo(np.longdouble).nmant
    results['double_mantissa_bits'] = np.finfo(np.float64).nmant
    # Extrapolation removes the leading h^4 truncation term.
    results['extended_Richardson_BA_discrepancy'] = discrepancy(
        (16*trajectories['extended', 48000]-trajectories['extended', 24000])/15,
        reference)
    path = ROOT/'results/trigonal_convergence_diagnostic.json'
    path.write_text(json.dumps(results, indent=2)+'\n')
    print('Wrote', path, flush=True)


def multiprecision_check(example='trigonal'):
    """Use genuinely wider arithmetic when NumPy long double is binary64."""
    from mpmath import mp
    ctx = mp.clone()
    ctx.dps = 35
    data = np.load(ROOT/f'results/{example}_sigma_ba.npz')
    reference = data['BA']
    L = M = 2 if example == 'hyperelliptic' else 3
    K, P = L+M, L*M
    # Promote the same stored binary64 coefficients and initial values.
    # This isolates integration arithmetic, not parameter reconstruction.
    phases = [ctx.mpc(complex(z)) for z in data['phases']]
    initial = [ctx.mpc(complex(z)) for z in reference[0]]
    coupling = ctx.mpf(math.sqrt(2))

    def rhs(z):
        u, partner = z[:K], z[K:2*K]
        C, Cbar = z[2*K:2*K+P], z[2*K+P:]
        du = [phases[k]*u[k] for k in range(K)]
        dp = [-phases[k]*partner[k] for k in range(K)]
        dc, dcb = [], []
        for i in range(L):
            for j in range(M):
                k = M*i+j
                du[i] += coupling*Cbar[k]*partner[L+j]
                du[L+j] += coupling*Cbar[k]*partner[i]
                dp[i] -= coupling*C[k]*u[L+j]
                dp[L+j] -= coupling*C[k]*u[i]
                dc.append(coupling*partner[i]*partner[L+j])
                dcb.append(-coupling*u[i]*u[L+j])
        return [-ctx.j*a for a in du+dp+dc+dcb]

    rows, trajectories = [], {}
    for steps in [6000, 12000, 24000]:
        z = initial.copy()
        h = ctx.mpf(str(float(data['analytic_t'][-1])))/steps
        out = [np.array([complex(v) for v in z])]
        for step in range(1, steps+1):
            k1 = rhs(z)
            k2 = rhs([v+h*a/2 for v, a in zip(z, k1)])
            k3 = rhs([v+h*a/2 for v, a in zip(z, k2)])
            k4 = rhs([v+h*a for v, a in zip(z, k3)])
            z = [v+h*(a+2*b+2*c+d)/6 for v, a, b, c, d in zip(z,k1,k2,k3,k4)]
            if step % (steps//(len(reference)-1)) == 0:
                out.append(np.array([complex(v) for v in z]))
        trajectory = np.array(out)
        trajectories[steps] = trajectory
        row = dict(steps=steps, BA_discrepancy=discrepancy(trajectory, reference))
        if steps//2 in trajectories:
            row['step_halving_difference'] = discrepancy(trajectory, trajectories[steps//2])
        rows.append(row)
        print('mpmath 35 digits', row, flush=True)
    path = ROOT/f'results/{example}_convergence_diagnostic.json'
    results = json.loads(path.read_text()) if path.exists() else {}
    results['multiprecision_dps'] = ctx.dps
    results['multiprecision'] = rows
    results['multiprecision_Richardson_BA_discrepancy'] = discrepancy(
        (16*trajectories[24000]-trajectories[12000])/15, reference)
    results['multiprecision_successive_difference_ratio'] = (
        rows[1]['step_halving_difference']/rows[2]['step_halving_difference'])
    results['precision_note'] = ('NumPy extended and double have identical mantissa sizes '
        'on this host. Only the mpmath run tests higher integration precision. '
        'It promotes the same saved coefficients and initial data exactly.')
    path.write_text(json.dumps(results, indent=2)+'\n')


if __name__ == '__main__':
    if '--multiprecision' in sys.argv:
        example = ('hyperelliptic' if '--hyperelliptic' in sys.argv else 'trigonal')
        multiprecision_check(example)
    else:
        main()
