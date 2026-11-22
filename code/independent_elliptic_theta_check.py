"""Cross-check public Generapy genus-one theta jets with mpmath.jtheta.

This independent validation backend is not used by the paper's production
BA evaluators or plots. Both the theta derivatives and BA-like ratios are
compared, including a point close to a theta zero and period shifts.
"""
from pathlib import Path
import json
import mpmath as mp
from generapy import rtheta_jet


def check():
    mp.mp.dps=55
    tau=mp.mpc('.23','1.17');nome=mp.exp(mp.pi*mp.j*tau)
    characteristic=([mp.mpf('.5')],[mp.mpf('.5')])
    points=[mp.mpc('.17','.04'),mp.mpc('.43','-.12'),mp.mpc('1e-7','2e-7')]
    points += [points[0]+1,points[0]+tau]
    errors=[];ratios=[]
    for z in points:
        jet=rtheta_jet([z],[[tau]],3,characteristic,ctx=mp.mp)
        # theta[1/2,1/2](z|tau)=-JacobiTheta1(pi*z,exp(pi*i*tau)).
        reference=[-mp.pi**k*mp.jtheta(1,mp.pi*z,nome,k) for k in range(4)]
        errors.extend(abs(jet[(k,)]-reference[k])/(1+abs(reference[k])) for k in range(4))
        ratios.append(abs(jet[(1,)]/jet[(0,)]-reference[1]/reference[0])/(1+abs(reference[1]/reference[0])))
    def theta(z):return rtheta_jet([z],[[tau]],1,characteristic,ctx=mp.mp)
    def jacobi(z):return [-mp.jtheta(1,mp.pi*z,nome),-mp.pi*mp.jtheta(1,mp.pi*z,nome,1)]
    for z in points[:2]:
        for shift in (mp.mpc('.14','.08'),mp.mpc('-.25','.19')):
            A=theta(z+shift);B=theta(z);a=jacobi(z+shift);b=jacobi(z)
            ratios.append(abs(A[(0,)]/B[(0,)]-a[0]/b[0])/(1+abs(a[0]/b[0])))
            x=A[(1,)]/A[(0,)]-B[(1,)]/B[(0,)];y=a[1]/a[0]-b[1]/b[0]
            ratios.append(abs(x-y)/(1+abs(y)))
    assert max(errors+ratios)<mp.mpf('1e-45')
    return {'dps':55,'independent_backend':'mpmath.jtheta, analytic Jacobi derivatives',
            'maximum_scaled_jet_error':mp.nstr(max(errors),14),
            'maximum_scaled_ratio_and_log_rate_error':mp.nstr(max(ratios),14),
            'production_backend':'Generapy rtheta_jet; unchanged'}


if __name__=='__main__':
    data=check();(Path(__file__).parent/'results/independent_elliptic_theta_check.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))
