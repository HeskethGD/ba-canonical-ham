"""Paper conventions adapted to Generapy's public sigma and BA jets.

No curve-specific sigma scale or sigma-sharp denominator is used. The
supplied second-kind rate and differential direction must be compatible.
"""
import numpy as np
from generapy import baker_akhiezer_jet, kleinian_sigma_directional_jet


class GenerapySigma:
    """Sigma derivatives needed separately for the pair-mode potential."""

    def __init__(self, ctx, omega, tau, kappa, characteristic, coordinate=-1):
        self.ctx = ctx
        self.omega = omega
        self.tau = tau
        self.kappa = kappa
        self.characteristic = characteristic
        self.genus = omega.rows
        self.coordinate = coordinate % self.genus
        self.direction = tuple(int(j == self.coordinate) for j in range(self.genus))

    def jet(self, u, order=2):
        u = [self.ctx.mpc(z.item() if isinstance(z, np.generic) else z) for z in u]
        # Match the public BA evaluator's outer guard, so its denominator
        # requests can reuse the same precision-aware theta cache entry.
        with self.ctx.extraprec(10):
            return kleinian_sigma_directional_jet(
                u, self.omega, self.tau, self.kappa, order,
                self.characteristic, normalization='theta',
                direction=self.direction, ctx=self.ctx)


class ReferenceSigmaBA:
    """Ordinary/dual BA jets in the paper's geometric Abel time.

    For z = z0 + t*e_coordinate, the ordinary and dual numerators are
    sigma(A+z) and sigma(A-z), with exponential rates +rate and -rate.
    The network later converts geometric time to -i*physical_time.

    Generapy uses sigma(shift-z)*exp(-R.(z-z0)). For the ordinary mode
    choose shift=-A, R_coordinate=-rate; for the dual choose shift=A,
    R_coordinate=rate. Sigma parity cancels in the ordinary reference
    ratio. This requires the half-integral characteristic used by the
    paper's examples. Other R components are zero because only one
    coordinate evolves; this is not a full second-kind integral vector
    for evolution in arbitrary Abel directions.
    """

    def __init__(self, sigma, abel, rate, z0):
        self.sigma = sigma
        ctx = sigma.ctx
        if any(not ctx.isint(2*a) for row in sigma.characteristic for a in row):
            raise ValueError('The BA pair requires a half-integral characteristic')
        self.z0 = tuple(ctx.mpc(z) for z in z0)
        abel = tuple(ctx.mpc(a) for a in abel)
        rate = ctx.mpc(rate)
        self.data = []
        for sign in (1, -1):
            shift = tuple(-sign*a for a in abel)
            second_kind = tuple(-sign*rate if j == sigma.coordinate else ctx.zero
                                for j in range(sigma.genus))
            self.data.append((shift, second_kind))

    def jets(self, t, order=3):
        if order not in (2, 3):
            raise ValueError('BA jets require order two or three')
        sigma = self.sigma
        z = list(self.z0)
        z[sigma.coordinate] += sigma.ctx.mpc(t)
        return tuple(baker_akhiezer_jet(
            z, shift=shift, second_kind=second_kind,
            omega=sigma.omega, tau=sigma.tau, kappa=sigma.kappa,
            characteristic=sigma.characteristic, direction=sigma.direction,
            order=order, normalization='reference', reference=self.z0,
            ctx=sigma.ctx) for shift, second_kind in self.data)
