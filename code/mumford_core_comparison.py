"""Check the Mumford polynomial core against our BA spectral flow.

Source: Fernandes--Vanhaecke, math-ph/0011051, Section 3.1, (4)--(6).
These are literature brackets, not a new Poisson construction. Polynomial
u,v,w here are NOT individual BA modes or their full products U,V.
"""
from itertools import combinations
from pathlib import Path
import json
import sympy as s


def build(g,phi):
    x,y=s.symbols('x y')
    uc=s.symbols(f'u0:{g}');vc=s.symbols(f'v0:{g}');wc=s.symbols(f'w0:{g+1}')
    u=x**g+sum(uc[r]*x**r for r in range(g))
    v=sum(vc[r]*x**r for r in range(g))
    w=x**(g+1)+sum(wc[r]*x**r for r in range(g+1))
    uy,vy,wy=[p.subs(x,y) for p in [u,v,w]]
    phiy=phi.subs(x,y)
    quot=lambda value:s.cancel(value/(x-y))
    Buv=quot(u*phiy-uy*phi)
    Buw=-2*quot(v*phiy-vy*phi)
    Bvw=quot(w*phiy-wy*phi)-u*phiy
    Bww=2*(v*phiy-vy*phi)
    coords=list(uc+vc+wc);n=len(coords);P=s.zeros(n)
    def putblock(expr,left,right,offsetl,offsetr):
        poly=s.Poly(expr,x,y)
        for i in range(left):
            for j in range(right):
                value=poly.coeff_monomial(x**i*y**j)
                P[offsetl+i,offsetr+j]=value
                P[offsetr+j,offsetl+i]=-value
    putblock(Buv,g,g,0,g);putblock(Buw,g,g+1,0,2*g)
    putblock(Bvw,g,g+1,g,2*g)
    ww=s.Poly(Bww,x,y)
    for i in range(g+1):
        for j in range(g+1):P[2*g+i,2*g+j]=ww.coeff_monomial(x**i*y**j)
    return coords,P,u,v,w


def check(g,root_adapted):
    x=s.Symbol('x')
    nodes=list(map(s.Integer,range(g)))
    phi=s.prod(x-node for node in nodes) if root_adapted else s.Integer(1)
    coords,P,u,v,w=build(g,phi)
    for i,j,k in combinations(range(len(coords)),3):
        jac=sum(P[i,l]*s.diff(P[j,k],coords[l])+P[j,l]*s.diff(P[k,i],coords[l])
                +P[k,l]*s.diff(P[i,j],coords[l]) for l in range(len(coords)))
        assert s.expand(jac)==0
    def grad(h):return s.Matrix([s.diff(h,c) for c in coords])
    def pb(f,h):return s.expand((grad(f).T*P*grad(h))[0])
    H=s.Poly(s.expand(u*w+v*v),x)
    coefficients=[H.nth(r) for r in range(2*g+1)]
    for f,h in combinations(coefficients,2):assert pb(f,h)==0
    if root_adapted:
        casimirs=[H.nth(2*g)]+[H.eval(node) for node in nodes]
    else:casimirs=[H.nth(r) for r in range(g,2*g+1)]
    for c in casimirs:assert (P*grad(c)).applyfunc(s.expand)==s.zeros(len(coords),1)
    quotient=s.Poly(s.div(H.as_expr(),phi,x)[0],x)
    generator=quotient.nth(g-1)
    q=w.coeff(x,g)-u.coeff(x,g-1)
    assert s.expand(pb(u,generator)-2*v)==0
    assert s.expand(pb(v,generator)-((x+q)*u-w))==0
    assert s.expand(pb(w,generator)+2*(x+q)*v)==0
    sample={c:s.Integer(j+2) for j,c in enumerate(coords)}
    rank=P.subs(sample).rank()
    assert rank==2*g
    report=dict(genus=g,phi=str(phi),phase_dimension=len(coords),generic_rank=rank,
                jacobi_verified=True,all_spectral_coefficients_in_involution=True,
                casimirs_count=len(casimirs),same_BA_polynomial_flow=True,
                Hamiltonian=str(generator))
    if g==2:
        report['spectral_integrals']=list(map(str,coefficients))
        report['casimirs']=list(map(str,casimirs))
        report['coordinates']=list(map(str,coords))
        report['poisson_matrix']=[[str(entry) for entry in row] for row in P.tolist()]
        if root_adapted:
            # Canonical power coordinates at the roots 0,1 of phi.
            I=[u.subs(x,node) for node in nodes]
            momentum=[-v.subs(x,node)/(s.diff(phi,x).subs(x,node)*I[i])
                      for i,node in enumerate(nodes)]
            for i in range(2):
                for j in range(2):
                    assert s.cancel(pb(I[i],momentum[j])-(1 if i==j else 0))==0
                    assert s.cancel(pb(I[i],I[j]))==0
                    assert s.cancel(pb(momentum[i],momentum[j]))==0
            report['canonical_power_coordinates']=list(map(str,I))
            report['canonical_power_momenta']=list(map(str,momentum))
    return report


if __name__=='__main__':
    results=[]
    for g in range(1,4):
        for adapted in [False,True]:
            print('Mumford core:',g,'root-adapted:',adapted,flush=True)
            results.append(check(g,adapted))
    path=Path(__file__).parent/'results'/'mumford_core_comparison.json'
    path.write_text(json.dumps(dict(checks=results),indent=2)+'\n')
    print(json.dumps([{k:v for k,v in row.items() if k not in ['poisson_matrix','spectral_integrals']}
                      for row in results],indent=2))
