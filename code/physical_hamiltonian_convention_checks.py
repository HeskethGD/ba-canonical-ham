"""Check the -i canonical convention against independent matrix identities.

The local matrix reference uses geometric unit-bracket rates. Its
Hamiltonian, assembly and finite differences independently check the
physical -i convention used in the paper.
"""
from pathlib import Path
import json
from support import matrix_reference
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
old=vars(matrix_reference)
rng=np.random.default_rng(20261114)
def random(shape):return rng.normal(size=shape)+1j*rng.normal(size=shape)
def comm(A,B):return A@B-B@A
records=[]
for L,M,N in [(1,1,1),(2,2,1),(2,3,2),(3,1,2)]:
    for real_form in [False,True]:
        phases=[rng.normal(size=n) if real_form else random(n) for n in (L,M,N)]
        arrays=[]
        for shape in [(L,M),(M,N),(N,L)]:
            a=random(shape);arrays.extend([a,a.conj() if real_form else random(shape)])
        args=(L,M,N,*phases,*arrays)
        X,D=old['build'](*args)
        rates=tuple(-1j*a for a in old['rates'](*args))
        dX,_=old['build'](L,M,N,*[a*0 for a in phases],*rates)
        nu=.37+.21j;Pencil=X+nu*D
        T=old['T'](L,M,N,phases[0],phases[1],*arrays[2:],nu)
        lax=float(np.max(abs(1j*dX-comm(Pencil,T))))
        gradient=[];eps=1e-5
        for k in range(6):
            g=np.empty_like(arrays[k])
            for idx in np.ndindex(g.shape):
                plus=[a.copy() for a in arrays];minus=[a.copy() for a in arrays]
                plus[k][idx]+=eps;minus[k][idx]-=eps
                g[idx]=(old['H'](L,M,N,*phases,*plus)-old['H'](L,M,N,*phases,*minus))/(2*eps)
            gradient.append(g)
        expected=[(-1j if k%2==0 else 1j)*gradient[k+1 if k%2==0 else k-1] for k in range(6)]
        ham=max(float(np.max(abs(a-b))) for a,b in zip(rates,expected))
        reality=(max(float(np.max(abs(rates[k+1]-rates[k].conj()))) for k in (0,2,4)) if real_form else None)
        r=4;Tr=(np.linalg.matrix_power(X,r-1)-np.linalg.matrix_power(Pencil,r-1))/nu
        higher=float(np.max(abs(comm(D,np.linalg.matrix_power(X,r-1))-comm(Pencil,Tr))))
        assert lax<1e-10 and ham<1e-7 and higher<1e-9
        if real_form:assert reality<1e-12 and abs(old['H'](*args).imag)<1e-10
        records.append(dict(sizes=[L,M,N],real_form=real_form,hamiltonian_error=ham,
                            lax_error=lax,higher_flow_error=higher,reality_error=reality))
path=ROOT/'code/results/physical_hamiltonian_convention_checks.json'
path.write_text(json.dumps(records,indent=2)+'\n')
print('PASS:',len(records),'matrix/conjugacy cases; maximum Lax error',max(a['lax_error'] for a in records))
