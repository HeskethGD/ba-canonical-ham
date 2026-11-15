"""Audit the paper's elementary BA normalization against Generapy.

Uses shifted genus-one through genus-three curves so that lambda_(2g) is
nonzero. No changes to Generapy; first- and second-kind lifts share paths.
"""
from pathlib import Path
import json
import sympy as s
from mpmath import mp
from generapy import Curve, kleinian_p
from generapy.kleinian import kleinian_baker_akhiezer


def check(genus):
    ctx = mp.clone()
    ctx.dps = 32
    x = s.Symbol('x')
    f = s.Poly(4*s.prod(x-s.Rational(1, 5)-j
                       for j in range(-genus, genus+1)), x)
    coefficients = [ctx.mpf(str(f.nth(i))) for i in range(2*genus+2)]
    curve = Curve({(0, 2): 1, **{(i, 0): -a
                  for i, a in enumerate(coefficients)}}, ctx=ctx)
    first, second = curve.periods_kind_1(), curve.periods_kind_2()
    characteristic = curve.riemann_constant().characteristic
    a = ctx.mpf(genus+2)
    y = ctx.sqrt(ctx.polyval(list(reversed(coefficients)), a))
    points = []
    for ordinate in (y, -y):
        target = (a, ordinate)
        points.append((curve.abel_map_kind_1(target).value,
                       curve.abel_map_kind_2(target).value))
    z = [ctx.mpf('0.13') + ctx.mpf('0.02')*i for i in range(genus)]

    def psi(t, sheet):
        argument = z[:-1]+[t]
        abel, integral = points[sheet]
        return kleinian_baker_akhiezer(
            argument, abel, integral, first.omega, first.tau,
            second.kappa, characteristic, ctx=ctx)

    def polynomial(t):
        argument = z[:-1]+[t]
        return a**genus-sum(kleinian_p(
            argument, first.omega, first.tau, second.kappa,
            (genus-1, j), characteristic, ctx=ctx)*a**j
            for j in range(genus))

    t = z[-1]
    P, J = polynomial(t), ctx.diff(polynomial, t)
    p, ps = psi(t, 0), (-1)**genus*psi(t, 1)
    dp = ctx.diff(lambda q: psi(q, 0), t)
    dps = (-1)**genus*ctx.diff(lambda q: psi(q, 1), t)
    potential = (2*kleinian_p(z, first.omega, first.tau, second.kappa,
                             (genus-1, genus-1), characteristic, ctx=ctx)
                 +coefficients[2*genus]/4)
    pairs = {
        'sheet_product': (p*ps, P),
        'minus_sheet_rate': (dp/p, (J-y)/(2*P)),
        'plus_sheet_rate': (dps/ps, (J+y)/(2*P)),
        'wronskian': (p*dps-dp*ps, y),
        'potential': (ctx.diff(lambda q: psi(q, 0), t, 2)/p, a+potential),
    }
    errors = {name: abs(got-want)/max(1, abs(got), abs(want))
              for name, (got, want) in pairs.items()}
    for name, error in errors.items():
        assert error < ctx.mpf('1e-24'), (genus, name, ctx.nstr(error))
    report = {'genus': genus, 'curve_coefficients': [str(v) for v in coefficients],
              'raw_sheet_product_sign': (-1)**genus,
              'scaled_errors': {k: ctx.nstr(v, 8) for k, v in errors.items()}}
    print(json.dumps(report), flush=True)
    return report


if __name__ == '__main__':
    reports = [check(g) for g in (1, 2, 3)]
    target = Path(__file__).with_name('results')/'latex_final_ba_conventions.json'
    target.write_text(json.dumps(reports, indent=2)+'\n')
