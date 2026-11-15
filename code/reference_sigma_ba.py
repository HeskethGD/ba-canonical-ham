"""Reference-normalized BA evaluations using existing sigma data.

No curve-specific sigma scale or sigma-sharp denominator is used. The
supplied second-kind rate and differential direction must be compatible.
"""
import numpy as np
from generapy.kleinian import kleinian_sigma,kleinian_sigma_jet


class GenerapySigma:
    """Existing Generapy sigma and its jet along one Baker coordinate."""

    def __init__(self, ctx, omega, tau, kappa, characteristic, coordinate=-1):
        self.ctx=ctx;self.omega=omega;self.tau=tau;self.kappa=kappa
        self.characteristic=characteristic;self.genus=omega.rows
        self.coordinate=coordinate%self.genus

    def jet(self,u,order=2):
        ctx=self.ctx
        u=[ctx.mpc(z.item() if isinstance(z,np.generic) else z) for z in u]
        if order==0:
            return (kleinian_sigma(u,self.omega,self.tau,self.kappa,
                                   self.characteristic,normalization='theta',ctx=ctx),)
        values=kleinian_sigma_jet(u,self.omega,self.tau,self.kappa,order,
                                  self.characteristic,normalization='theta',ctx=ctx)
        result=[]
        for n in range(order+1):
            key=[0]*self.genus;key[self.coordinate]=n
            result.append(values[tuple(key)])
        return tuple(result)


class ReferenceSigmaBA:
    """Ordinary/dual BA pair and derivatives through order three, normalized at z0.

    The argument t here is the geometric Abel displacement, so
    z=z0+t*e_coordinate. The network uses displacement -i*physical_time. ``rate`` is the compatible regularized
    second-kind integral in that direction; it is fixed per spectral point.
    """

    def __init__(self,sigma,abel,rate,z0,reference_denominator=None):
        self.sigma=sigma;ctx=sigma.ctx
        self.abel=tuple(ctx.mpc(a) for a in abel)
        self.rate=ctx.mpc(rate);self.z0=tuple(ctx.mpc(z) for z in z0)
        denominator=(sigma.jet(self.z0,0)[0] if reference_denominator is None
                     else reference_denominator)
        self.fref=sigma.jet([a+z for a,z in zip(self.abel,self.z0)],0)[0]/denominator
        self.gref=sigma.jet([a-z for a,z in zip(self.abel,self.z0)],0)[0]/denominator
        if any(v==0 or not ctx.isfinite(v) for v in (self.fref,self.gref,denominator)):
            raise ValueError('Singular BA reference chart')

    def jets(self,t,denominator=None,order=3):
        if order not in (2,3):raise ValueError('BA jets require order two or three')
        z=list(self.z0);z[self.sigma.coordinate]+=self.sigma.ctx.mpc(t)
        d=self.sigma.jet(z,order) if denominator is None else denominator
        logd=d[1]/d[0];logd2=d[2]/d[0]-logd**2
        result=[]
        for sign,ref in ((1,self.fref),(-1,self.gref)):
            arg=[a+sign*zz for a,zz in zip(self.abel,z)]
            n=self.sigma.jet(arg,order)
            value=self.sigma.ctx.exp(sign*self.rate*t)*(n[0]/d[0])/ref
            rate=sign*n[1]/n[0]-logd+sign*self.rate
            rateprime=n[2]/n[0]-(n[1]/n[0])**2-logd2
            jet=[value,value*rate,value*(rate*rate+rateprime)]
            if order==3:
                logn3=n[3]/n[0]-3*n[2]*n[1]/n[0]**2+2*(n[1]/n[0])**3
                logd3=d[3]/d[0]-3*d[2]*d[1]/d[0]**2+2*logd**3
                rate2=sign*logn3-logd3
                jet.append(value*(rate**3+3*rate*rateprime+rate2))
            result.append(tuple(jet))
        return tuple(result)
