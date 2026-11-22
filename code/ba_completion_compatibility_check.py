"""Genus-two finite completion and deliberately incompatible BA weights.

Complete one prescribed point/weight in each group by three new points.
Exact rational interpolation preserves the prescribed weights and avoids
zero new weights. Independent public Generapy BA jets then check all
canonical rates. Perturbing either group's weights makes the elementary
rates require time-dependent phase coefficients.
"""
from pathlib import Path
import importlib.metadata
import json
import sympy as sp
from mpmath import mp
from generapy import Curve, baker_akhiezer_jet


def complete(nodes, prescribed):
    matrix = sp.Matrix([[1]*3, nodes[1:]])
    target = sp.Matrix([-prescribed, sp.Rational(1, 2)-prescribed*nodes[0]])
    particular, parameters = matrix.gauss_jordan_solve(target)
    parameter = parameters[0]
    # The null vector has no zero coordinate, so only finitely many values
    # are forbidden. Choose deterministically without relying on a grid fit.
    for integer in range(1, 8):
        solution = particular.subs(parameter, sp.Rational(integer, 13))
        if all(solution):
            weights = [prescribed, *solution]
            assert sum(weights) == 0
            assert sum(a*k for a, k in zip(nodes, weights)) == sp.Rational(1, 2)
            assert all(matrix[:, [i, j]].det() != 0 for i in range(3) for j in range(i+1, 3))
            return weights
    raise AssertionError('No nonzero point on the completion affine line')


def check():
    ctx = mp.clone(); ctx.dps = 36
    groups = [[sp.Rational(x, 3) for x in (9, 11, 13, 15)],
              [sp.Rational(x, 3) for x in (10, 12, 14, 16)]]
    exact = complete(groups[0], sp.Rational(2, 7))+complete(groups[1], -sp.Rational(3, 8))
    a = [ctx.mpf(str(x.p))/int(x.q) for group in groups for x in group]
    weights = [ctx.mpf(str(x.p))/int(x.q) for x in exact]
    f = lambda x: 4*x*(x*x-1)*(x*x-4)
    y = [ctx.sqrt(f(x))*(1 if k < 4 else -1) for k, x in enumerate(a)]
    curve = Curve({(0, 2): 1, (5, 0): -4, (3, 0): 20, (1, 0): -16}, ctx=ctx)
    first, second = curve.periods_kind_1(), curve.periods_kind_2()
    char = curve.riemann_constant().characteristic
    fixed = []
    for x, ordinate in zip(a, y):
        fixed.append([(curve.abel_map_kind_1((x, sign*ordinate)).value,
                       curve.abel_map_kind_2((x, sign*ordinate)).value) for sign in (1, -1)])
    times = [ctx.zero, ctx.mpc('.12', '.03'), ctx.mpc('.25', '-.04')]
    jets = []
    for t in times:
        z = [ctx.mpc('.17', '.08'), ctx.mpc('.31', '.12')+t]
        jets.append([[baker_akhiezer_jet(z, shift=A, second_kind=R,
                      omega=first.omega, tau=first.tau, kappa=second.kappa,
                      characteristic=char, direction=[0, 1], order=2,
                      normalization='spectral', ctx=ctx) for A, R in pair] for pair in fixed])

    def residuals(kap):
        errors, effective = [], []
        phase = [-sum(kap[j]*(y[j]-y[i])/(a[j]-a[i])
                      for j in (range(4, 8) if i < 4 else range(4))) for i in range(8)]
        for raw in jets:
            u = [pair[0] for pair in raw]
            v = [[kap[k]*entry for entry in pair[1]] for k, pair in enumerate(raw)]
            C, Cbar, dC, dCbar = {}, {}, {}, {}
            for i in range(4):
                for j in range(4, 8):
                    scale = ctx.sqrt(2)/(a[j]-a[i])
                    C[i, j] = scale*(v[i][0]*v[j][1]-v[i][1]*v[j][0])
                    Cbar[i, j] = -scale*(u[i][0]*u[j][1]-u[i][1]*u[j][0])
                    dC[i, j] = scale*(v[i][0]*v[j][2]-v[i][2]*v[j][0])
                    dCbar[i, j] = -scale*(u[i][0]*u[j][2]-u[i][2]*u[j][0])
            state_errors, b_needed = [], []
            scaled = lambda p, q: abs(p-q)/(1+max(abs(p), abs(q)))
            for k in range(8):
                opposite = range(4, 8) if k < 4 else range(4)
                key = lambda j: (k, j) if k < 4 else (j, k)
                mix = ctx.sqrt(2)*sum(v[j][0]*Cbar[key(j)] for j in opposite)
                dual_mix = ctx.sqrt(2)*sum(u[j][0]*C[key(j)] for j in opposite)
                state_errors += [scaled(u[k][1], phase[k]*u[k][0]+mix),
                                 scaled(v[k][1], -phase[k]*v[k][0]-dual_mix)]
                b_needed.append((u[k][1]-mix)/u[k][0])
            for i in range(4):
                for j in range(4, 8):
                    state_errors += [scaled(dC[i, j], ctx.sqrt(2)*v[i][0]*v[j][0]),
                                     scaled(dCbar[i, j], -ctx.sqrt(2)*u[i][0]*u[j][0])]
            errors.append(max(state_errors)); effective.append(b_needed)
        phase_variation = max(abs(effective[n][k]-effective[0][k])
                              for n in range(1, len(times)) for k in range(8))
        return max(errors), phase_variation

    good_error, good_variation = residuals(weights)
    assert good_error < ctx.mpf('1e-24'), good_error
    assert good_variation < ctx.mpf('1e-24'), good_variation
    rejected = []
    for index in (0, 4):
        bad = weights.copy(); bad[index] += ctx.mpf('0.05')
        error, variation = residuals(bad)
        assert error > ctx.mpf('1e-6') and variation > ctx.mpf('1e-6'), (index, error, variation)
        rejected.append({'perturbed_group': 'L' if index == 0 else 'M',
                         'mode_rate_residual': str(error),
                         'required_phase_variation': str(variation)})
    return {'genus': 2, 'L_M_N': [4, 4, 1], 'canonical_pairs': 24,
            'dps': ctx.dps, 'generapy_version': importlib.metadata.version('generapy'),
            'prescribed_weights': ['2/7', '-3/8'],
            'completed_weights': list(map(str, exact)),
            'all_new_weights_nonzero': True, 'exact_moment_equations': True,
            'all_48_mode_rate_error': str(good_error),
            'compatible_phase_variation': str(good_variation),
            'incompatible_cases': rejected,
            'scope': 'A positive-genus completion and two negative fixtures; not a numerical proof of necessity.'}


if __name__ == '__main__':
    report = check()
    (Path(__file__).parent/'results/ba_completion_compatibility_check.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)
