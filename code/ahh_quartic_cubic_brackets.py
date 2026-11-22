"""Compare canonical parent, induced Lie bracket, and frozen cubic bracket.

Brackets here omit the common physical factor -i.
"""
import json
from pathlib import Path
import sympy as sp


def check():
    q,d=2,3
    F=sp.Matrix(q,d,lambda i,a:sp.Symbol(f'f{i}{a}'))
    G=sp.Matrix(q,d,lambda i,a:sp.Symbol(f'g{i}{a}'))
    X=-G.T*F
    pb=lambda a,b:sp.expand(sum(sp.diff(a,f)*sp.diff(b,g)-sp.diff(a,g)*sp.diff(b,f)
                               for f,g in zip(F,G)))
    delta=lambda a,b:sp.Integer(a==b)
    for a in range(d):
        for b in range(d):
            for c in range(d):
                for e in range(d):
                    expected=delta(a,e)*X[c,b]-delta(b,c)*X[a,e]
                    assert sp.expand(pb(X[a,b],X[c,e])-expected)==0
    Z=sp.Matrix(d,d,lambda a,b:sp.Symbol(f'x{a}{b}'))
    y=sp.symbols('y0:3');D=sp.diag(*y)
    H2=sp.trace(D*Z**2);H3=sp.trace(Z**3)/3
    grad=lambda h:h.diff(Z).T
    lie=grad(H2)*Z-Z*grad(H2)
    frozen=D*grad(H3)-grad(H3)*D
    expand=lambda A:A.applyfunc(sp.expand)
    assert expand(lie-frozen)==sp.zeros(d)
    assert expand(lie-(D*Z**2-Z**2*D))==sp.zeros(d)
    shift=sp.Symbol('shift')
    for a in range(d):
        for b in range(d):
            for c in range(d):
                for e in range(d):
                    p1=delta(a,e)*Z[c,b]-delta(b,c)*Z[a,e]
                    shifted=delta(a,e)*(Z+shift*D)[c,b]-delta(b,c)*(Z+shift*D)[a,e]
                    p0=(y[a]-y[b])*delta(a,e)*delta(b,c)
                    assert sp.expand(shifted-p1+shift*p0)==0
    root=sp.sqrt(2)
    u,ub=Z[0,1]/root,Z[1,0]/root
    v,vb=Z[1,2],-Z[2,1]
    w,wb=-Z[2,0],Z[0,2]
    def p1(a,b):
        return sp.expand(sum(sp.diff(a,Z[i,j])*sp.diff(b,Z[k,l])*
                             (delta(i,l)*Z[k,j]-delta(j,k)*Z[i,l])
                             for i in range(d) for j in range(d)
                             for k in range(d) for l in range(d)))
    assert p1(u,v)==-wb/root
    assert p1(v,w)==root*ub
    assert p1(w,u)==-vb/root
    assert sp.expand(p1(u,ub)-(Z[1,1]-Z[0,0])/2)==0
    result={'parent_induced_Lie_bracket':True,
            'same_flow_quadratic_Lie_and_cubic_frozen_Hamiltonians':True,
            'Poisson_pencil_from_constant_shift':True,
            'explicit_three_mode_cross_brackets':True,
            'convention':'Common physical factor -i suppressed.'}
    return result


if __name__=='__main__':
    result=check()
    target=Path(__file__).parent/'results'/'ahh_quartic_cubic_brackets.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
