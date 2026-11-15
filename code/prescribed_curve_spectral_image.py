"""Check the prescribed-curve spectral map using recorded BA realizations.

No sigma reevaluation is needed: the initial modes and curve points are
read from the current records. This is a numerical check of the algebraic
map, not a substitute for the eigenvector and birationality proofs.
"""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent
rng=np.random.default_rng(20261114)
def unpack(value):return complex(float(value[0]),float(value[1]))
def pencil(z,phases,L,M,nu,rho=0):
    K=L+M;count=L*M
    u,ubar=z[:K],z[K:2*K]
    C=z[2*K:2*K+count].reshape(L,M)
    Cbar=z[2*K+count:].reshape(L,M)
    return np.block([
        [np.diag(rho-phases[:L]+nu),np.sqrt(2)*C,ubar[:L,None]],
        [np.sqrt(2)*Cbar.T,np.diag(rho-phases[L:]-nu),u[L:,None]],
        [-u[None,:L],-ubar[None,L:],np.array([[rho]])]])
records=[]
for kind in ['hyperelliptic','trigonal']:
    geometry=json.loads((ROOT/'results'/f'{kind}_sigma_ba_geometry_24.json').read_text())
    data=np.load(ROOT/'results'/f'{kind}_sigma_ba.npz')
    report=json.loads((ROOT/'results'/f'{kind}_sigma_ba.json').read_text())
    L,M=report['group_sizes'];K=L+M
    x=np.array([unpack(p['x']) for p in geometry['points']])
    y=np.array([unpack(p['y']) for p in geometry['points']])
    if kind=='hyperelliptic':
        weights=np.r_[1/(2*(x[0]-x[1])),1/(2*(x[1]-x[0])),
                      1/(2*(x[2]-x[3])),1/(2*(x[3]-x[2]))]
        def h(group,a,b):return sum(weights[k]*(b+y[k])/(a-x[k]) for k in group)
        def ordinate(a):return np.sqrt(complex(4*a*(a*a-1)*(a*a-4)))*rng.choice([-1,1])
    else:
        residues=[]
        for group in [slice(0,L),slice(L,K)]:
            residues.extend(np.linalg.solve(np.array([np.ones(len(x[group])),x[group],y[group]]),[0,0,1]))
        weights=np.array(residues)
        def h(group,a,b):return sum(weights[k]*(b*b+y[k]*b+y[k]**2)/(a-x[k]) for k in group)
        def ordinate(a):return complex(a**4+2*a**3+3*a*a+5*a+7)**(1/3)*np.exp(2j*np.pi*rng.integers(3)/3)
    phase=np.array([-h(range(L,K),x[i],y[i]) for i in range(L)]
                   +[h(range(L),x[j],y[j]) for j in range(L,K)])
    phase_error=float(np.max(abs(phase-data['phases'])))
    residuals=[]
    # Initial, interior and final BA states test the same polynomial map.
    for state in data['BA'][[0,len(data['BA'])//2,-1]]:
        for _ in range(12):
            a=rng.normal()+1j*rng.normal();b=ordinate(a)
            hL=h(range(L),a,b);hM=h(range(L,K),a,b)
            nu=-(hL+hM)/2;lam=(hM-hL)/2
            matrix=lam*np.eye(K+1)-pencil(state,phase,L,M,nu)
            sv=np.linalg.svd(matrix,compute_uv=False)
            residuals.append(float(sv[-1]/sv[0]))
    maximum=max(residuals)
    assert phase_error<1e-11 and maximum<1e-10,(kind,phase_error,maximum)
    records.append(dict(case=kind,genus=report['genus'],sizes=[L,M],
                        phase_reconstruction_error=phase_error,
                        tested_curve_points_and_states=len(residuals),
                        maximum_scaled_singularity_residual=maximum,
                        expected_total_excess_singularity_defect=L*M-report['genus']))
    print(kind,'PASS: phase error',phase_error,'spectral-map residual',maximum)
(ROOT/'results/prescribed_curve_spectral_image.json').write_text(json.dumps(records,indent=2)+'\n')
