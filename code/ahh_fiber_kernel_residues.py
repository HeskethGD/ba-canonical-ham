"""Check the complete-fiber residue identity with elliptic spectral functions."""
import json
from pathlib import Path
import mpmath as mp


def check():
    mp.mp.dps=55
    tau=mp.mpc('.23','1.17')
    nome=mp.exp(mp.pi*mp.j*tau)
    th=lambda z:mp.jtheta(1,mp.pi*z,nome)
    thp=lambda z:mp.pi*mp.jtheta(1,mp.pi*z,nome,1)
    thpp=lambda z:mp.pi**2*mp.jtheta(1,mp.pi*z,nome,2)
    e=lambda z:thp(z)/th(z)
    ep=lambda z:thpp(z)/th(z)-e(z)**2
    phi=lambda z,a:thp(0)*th(z+a)/(th(z)*th(a))
    p=[mp.mpc('.13','.14'),mp.mpc('.41','.29'),mp.mpc('.73','.53')]
    rho=list(map(mp.mpf,['.7','-.2','-.5']))
    kappa=list(map(mp.mpf,['.1','.6','-.7']))
    alpha=mp.mpc('.31','.22')
    lc=mp.mpc('.17','-.09');wc=mp.mpc('-.12','.07')
    lam=lambda u:lc+sum(rho[a]*e(u-p[a]) for a in range(3))
    lamp=lambda u:sum(rho[a]*ep(u-p[a]) for a in range(3))
    w=lambda u:wc+sum(kappa[a]*e(u-p[a]) for a in range(3))

    def normalize(u):
        b=mp.im(u)/mp.im(tau)
        a=mp.re(u)-mp.re(tau)*b
        return u-mp.floor(a)-mp.floor(b)*tau

    roots=[]
    for j in range(6):
        for k in range(6):
            seed=(mp.mpf(j)+mp.mpf('.37'))/6+tau*(mp.mpf(k)+mp.mpf('.43'))/6
            try:
                root=normalize(mp.findroot(lambda u:lam(u)-alpha,seed,
                                           df=lamp,solver='newton',maxsteps=60))
            except (ValueError,ZeroDivisionError):
                continue
            if abs(lam(root)-alpha)>mp.mpf('1e-40'):
                continue
            if all(abs(root-old)>mp.mpf('1e-35') for old in roots):
                roots.append(root)
            if len(roots)==3:
                break
        if len(roots)>=2:
            break
    # Abel's theorem fixes the last point of a degree-three fiber once two
    # are known. Refine it with the defining equation rather than accept
    # the Abel relation as a numerical root test.
    if len(roots)==2:
        seed=normalize(sum(p)-sum(roots))
        root=normalize(mp.findroot(lambda u:lam(u)-alpha,seed,
                                   df=lamp,solver='newton',maxsteps=60))
        assert abs(lam(root)-alpha)<mp.mpf('1e-40')
        assert all(abs(root-old)>mp.mpf('1e-35') for old in roots)
        roots.append(root)
    assert len(roots)==3,roots
    weights=[w(u)/lamp(u) for u in roots]
    y=[kappa[a]/rho[a] for a in range(3)]
    constants=[]
    for a in range(3):
        l0=lc+sum(rho[b]*e(p[a]-p[b]) for b in range(3) if b!=a)
        w0=wc+sum(kappa[b]*e(p[a]-p[b]) for b in range(3) if b!=a)
        constants.append(w0-y[a]*(l0-alpha))
    errors=[]
    for z in [mp.mpc('.11','.06'),mp.mpc('.24','-.04'),mp.mpc('.003','.002')]:
        for a in range(3):
            for b in range(3):
                value=sum(weights[i]*phi(z,p[a]-roots[i])*phi(z,roots[i]-p[b])
                          for i in range(3))
                expected=((y[a]-y[b])*phi(z,p[a]-p[b]) if a!=b
                          else constants[a]/rho[a])
                errors.append(abs(value-expected)/(1+abs(expected)))
    assert max(errors)<mp.mpf('1e-40'),max(errors)
    assert abs(sum(weights))<mp.mpf('1e-40')
    return {'dps':mp.mp.dps,'fiber_degree':3,'complete_regular_fiber':True,
            'scaled_residue_identity_error':mp.nstr(max(errors),12),
            'derived_zero_sum_weight_error':mp.nstr(abs(sum(weights)),12),
            'diagonal_and_all_off_diagonal_entries':True,
            'scope':'Elliptic verification of an arbitrary-curve residue theorem. '
                    'The curve functions and weights are constructed together; '
                    'the weights are not freely chosen.'}


if __name__=='__main__':
    result=check()
    target=Path(__file__).parent/'results'/'ahh_fiber_kernel_residues.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
