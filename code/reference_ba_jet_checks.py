"""Check the paper's public-BA adapter against differentiated sigma ratios.

Includes ordinary/dual signs, nonfinal Abel coordinates, numerator zeros
and regular-chart failures. No paper result files are overwritten.
"""
from mpmath import mp
from generapy import kleinian_sigma
from reference_sigma_ba import GenerapySigma, ReferenceSigmaBA


def fixture(genus, coordinate, odd):
    ctx = mp.clone()
    ctx.dps = 40
    characteristic = [[ctx.zero]*genus for _ in range(2)]
    if odd:
        characteristic[0][0] = characteristic[1][0] = ctx.mpf('.5')
    sigma = GenerapySigma(ctx, ctx.eye(genus)/2, ctx.j*ctx.eye(genus),
                          ctx.eye(genus)/10, characteristic, coordinate)
    abel = tuple(ctx.mpf('.37') + j/ctx.mpf(10) for j in range(genus))
    z0 = tuple(ctx.mpf('.2') + j/ctx.mpf(20) for j in range(genus))
    rate = ctx.mpc('.4', '.07')
    return sigma, abel, rate, z0


def check_point(sigma, abel, rate, z0, t):
    ctx = sigma.ctx
    mode = ReferenceSigmaBA(sigma, abel, rate, z0)

    def scalar_sigma(u):
        return kleinian_sigma(u, sigma.omega, sigma.tau, sigma.kappa,
                              sigma.characteristic, normalization='theta', ctx=ctx)

    jets = mode.jets(t)
    shorter = mode.jets(t, order=2)
    maximum = ctx.zero
    for side, sign in enumerate((1, -1)):
        normalizer = scalar_sigma([a+sign*z for a, z in zip(abel, z0)])/scalar_sigma(z0)

        def direct(q):
            z = list(z0)
            z[sigma.coordinate] += q
            return (ctx.exp(sign*rate*q)
                    * scalar_sigma([a+sign*v for a, v in zip(abel, z)])
                    / scalar_sigma(z) / normalizer)

        for n in range(4):
            expected = ctx.diff(direct, t, n)
            got = jets[side][n]
            assert ctx.isfinite(got), (side, n, t)
            residual = abs(got-expected)/max(1, abs(got), abs(expected))
            assert residual < ctx.mpf('1e-32'), (side, n, t, ctx.nstr(residual))
            maximum = max(maximum, residual)
            if n <= 2:
                assert abs(shorter[side][n]-got)/max(1, abs(got)) < ctx.mpf('1e-32')
    return maximum, jets


def require_singular(mode, t):
    try:
        mode.jets(t)
    except ZeroDivisionError:
        return
    raise AssertionError('A singular BA chart was accepted')


def main():
    for genus, coordinate, odd in ((1, -1, True), (2, 0, False), (2, -1, True)):
        sigma, abel, rate, z0 = fixture(genus, coordinate, odd)
        ctx = sigma.ctx
        maximum = ctx.zero
        for t in (ctx.zero, ctx.mpc('.17', '-.08')):
            residual, jets = check_point(sigma, abel, rate, z0, t)
            maximum = max(maximum, residual)
            if not t:
                assert all(abs(jet[0]-1) < ctx.mpf('1e-35') for jet in jets)
        print(f'PASS: genus {genus}, coordinate {sigma.coordinate}, odd={odd}; '
              f'maximum derivative residual {ctx.nstr(maximum, 5)}', flush=True)

    for genus in (1, 2):
        sigma, abel, rate, original_z0 = fixture(genus, 0, True)
        ctx = sigma.ctx
        for side, sign in enumerate((1, -1)):
            # Arrange A+sign*z=0 while varying only the first coordinate.
            z0 = (original_z0[0],) + tuple(-sign*a for a in abel[1:])
            zero = -sign*abel[0]-z0[0]
            for offset in (ctx.zero, ctx.mpf('1e-12'), -ctx.mpf('1e-12')):
                residual, jets = check_point(sigma, abel, rate, z0, zero+offset)
                if not offset:
                    assert abs(jets[side][0]) < ctx.mpf('1e-35')
                    assert abs(jets[side][1]) > ctx.mpf('1e-5')
            print(f'PASS: genus {genus} {("ordinary", "dual")[side]} '
                  'numerator zero and both sides', flush=True)

    sigma, abel, rate, z0 = fixture(1, -1, True)
    ctx = sigma.ctx
    require_singular(ReferenceSigmaBA(sigma, abel, rate, z0), -z0[0])
    require_singular(ReferenceSigmaBA(sigma, abel, rate, [ctx.zero]), 0)
    require_singular(ReferenceSigmaBA(sigma, abel, rate, abel), 0)
    require_singular(ReferenceSigmaBA(sigma, abel, rate, tuple(-a for a in abel)), 0)
    print('PASS: moving denominator and singular reference charts rejected', flush=True)


if __name__ == '__main__':
    main()
