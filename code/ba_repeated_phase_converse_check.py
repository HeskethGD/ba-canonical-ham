"""Reconstruct a repeated-phase (1,1,3) launch through active elliptic BA data.

Start with rational canonical modes, not prescribed curve samples. The
third block moment has rank two and a one-dimensional inactive eigenspace.
A constant canonical basis reduces the pencil to (1,1,2); its discriminant
has a forced square and a smooth quartic normalization. Recover both BA
divisors, use public Generapy theta jets, and lift the reconstructed flow
back to every original mode. No claim about all repeated-phase profiles
is inferred from this one test.
"""
from pathlib import Path
from functools import lru_cache
import importlib.metadata
import json
import sympy as sp
from mpmath import mp
from generapy import Curve, rtheta_jet


def check(dps=40):
    ctx=mp.clone();ctx.dps=dps
    z,nu=sp.symbols('lambda nu');rat=sp.Rational
    b=rat(1,3);bL=rat(1,7);bM=-rat(2,9)
    # All six original incident modes are nonzero, despite one inactive
    # linear combination. This is direct rational canonical launch data.
    X0=sp.Matrix([[-2*bL,(7+2*sp.I)/10,rat(17,40),-rat(7,40),-rat(1,40)],
                  [(4-3*sp.I)/10,-2*bM,rat(21,20),-rat(9,20),-rat(1,20)],
                  [-rat(17,10),-rat(7,10),b,0,0],
                  [-rat(11,5),-rat(9,10),0,b,0],
                  [-rat(19,10),-rat(4,5),0,0,b]])
    Dfull=sp.diag(1,-1,0,0,0)
    F=X0[:2,2:];G=X0[2:,:2]
    moment=-G[:,0]*F[0,:]+G[:,1]*F[1,:]
    assert F.rank()==G.rank()==moment.rank()==2
    assert (F*G).det()!=0
    assert len(F.nullspace())==1
    eigen=moment.eigenvects()
    active=[vector for value,multiplicity,vectors in eigen if value!=0 for vector in vectors]
    inactive=F.nullspace()[0]
    assert moment*inactive==sp.zeros(3,1)
    S=sp.Matrix.hstack(*active,inactive);T=sp.diag(sp.eye(2),S)
    adapted=sp.simplify(T.inv()*X0*T)
    assert adapted[:4,4]==sp.zeros(4,1) and adapted[4,:4]==sp.zeros(1,4)
    XS=adapted[:4,:4];D=sp.diag(1,-1,0,0)
    M=z*sp.eye(4)-XS-nu*D
    FS=sp.Poly(sp.expand(M.det()),nu)
    aa,bb,cc=[FS.nth(k) for k in (2,1,0)]
    normalized=sp.Poly(sp.cancel((bb**2-4*aa*cc)/(z-b)**2),z)
    assert normalized.degree()==4 and sp.gcd(normalized,normalized.diff()).degree()==0
    assert sp.expand((z*sp.eye(5)-X0-nu*Dfull).det()-(z-b)*FS.as_expr())==0
    A2=sp.cancel(aa/(z-b)**2);B1=sp.cancel(bb/(z-b))
    assert not A2.has(z)
    anchor=sp.Matrix([0,0,1,1]);adj=M.adjugate()
    section=sp.expand((anchor.T*adj*anchor)[0])
    resultant=sp.factor(sp.resultant(FS.as_expr(),section,nu))
    factors=sp.factor_list(resultant,z)[1]
    divisor_poly=sp.Poly(sp.prod(factor for factor,power in factors if sp.degree(factor,z)>1),z)
    assert divisor_poly.degree()==2
    cv=lambda a:ctx.mpc(str(sp.N(sp.re(a),dps+8)),str(sp.N(sp.im(a),dps+8)))
    coeff=lambda a:[cv(c) for c in sp.Poly(a,z).all_coeffs()]
    av,bv,ccv=map(coeff,(aa,bb,cc));fv=coeff(normalized.as_expr());b1v=coeff(B1)
    pol=lambda c,x:ctx.polyval(c,x)
    deriv=lambda c:[(len(c)-1-k)*a for k,a in enumerate(c[:-1])]
    bc=cv(b);a2=cv(A2);X=ctx.matrix([[cv(a) for a in row] for row in XS.tolist()])
    def adjugate(lam,n):
        m=lam*ctx.eye(4)-X-n*ctx.diag([1,-1,0,0])
        return ctx.matrix([[(-1)**(i+j)*ctx.det(ctx.matrix([[m[r,c] for c in range(4) if c!=i] for r in range(4) if r!=j])) for j in range(4)] for i in range(4)])
    branches=ctx.polyroots(fv,maxsteps=1000,extraprec=80)
    e=min(branches,key=abs)
    # x=1/(lambda-e), Y=y*x^2 turns the smooth quartic into an odd cubic.
    odd=[ctx.zero]*4;dd=fv[:]
    for k in range(1,5):
        dd=deriv(dd);odd[4-k]=pol(dd,e)/ctx.factorial(k)
    f=lambda x:ctx.polyval(list(reversed(odd)),x)
    fp=lambda x:sum(k*odd[k]*x**(k-1) for k in range(1,4))
    def point(lam,y):
        x=1/(lam-e);want=y*x*x;root=ctx.sqrt(f(x))
        return (x,min((root,-root),key=lambda r:abs(r-want)))
    ordinary=[];dual=[];root_errors=[]
    for lam in ctx.polyroots(coeff(divisor_poly.as_expr()),maxsteps=1000,extraprec=80):
        yy=ctx.sqrt(pol(fv,lam))
        for y in (yy,-yy):
            n=(-pol(b1v,lam)+y)/(2*a2*(lam-bc));A=adjugate(lam,n)
            col=A*ctx.matrix([0,0,1,1]);row=ctx.matrix([[0,0,1,1]])*A
            cr=max(abs(a) for a in col);rr=max(abs(a) for a in row)
            root_errors.append(min(cr,rr)/(1+max(cr,rr)))
            (ordinary if cr<rr else dual).append(point(lam,y))
    assert len(ordinary)==len(dual)==2
    R=[];weights=[];regular=[]
    yb=ctx.sqrt(pol(fv,bc));xr=1/(bc-e)
    for yy in (yb,-yb):
        rlam=(-pol(b1v,bc)+yy)/(2*a2)
        n0=(-pol(deriv(b1v),bc)+pol(deriv(fv),bc)/(2*yy))/(2*a2)
        R.append(point(bc,yy));weights.append(rlam*xr*xr)
        regular.append(n0-rlam*xr)
    # Reorder punctures to match the moment-diagonalized canonical basis.
    # nu ~ r_lambda/(lambda-b), so the moment eigenvalue is r_lambda.
    expected=[cv((-XS[2+k,0]*XS[0,2+k]+XS[2+k,1]*XS[1,2+k])) for k in range(2)]
    residues=[w/(xr*xr) for w in weights]
    if abs(residues[0]-expected[0])>abs(residues[1]-expected[0]):
        R.reverse();weights.reverse();regular.reverse()
    for k in range(2):assert abs(weights[k]/(xr*xr)-expected[k])<ctx.mpf('1e-28')
    # Their exact odd-model limits are x=0, Y=+/-sqrt(leading quartic coefficient).
    qy=ctx.sqrt(odd[0]);Q=[(ctx.zero,-qy),(ctx.zero,qy)]
    pv=lambda c,x:ctx.polyval(list(reversed(c)),x)
    def roots_poly(points):
        result=[ctx.one]
        for x,y in points:
            new=[ctx.zero]*(len(result)+1)
            for k,a in enumerate(result):new[k]-=x*a;new[k+1]+=a
            result=new
        return result
    def interpolate(points):
        return list(ctx.lu_solve(ctx.matrix([[x**j for j in range(len(points))] for x,y in points]),ctx.matrix([y for x,y in points])))
    qD=roots_poly(ordinary);qDual=roots_poly(dual)
    def basis(points,q):
        p=interpolate(points);fn=lambda P:(P[1]+pv(p,P[0]))/pv(q,P[0]);at=[fn(P) for P in R]
        return lambda P,k:(fn(P)-at[1-k])/(at[k]-at[1-k])
    psi0=basis(ordinary,qD);dual0=basis(dual,qDual)
    beta_scale=-weights[0]*R[0][1]/pv(qD,xr)
    beta=lambda P:beta_scale*pv(qD,P[0])/((P[0]-xr)**2*P[1])
    for k in range(2):assert abs(beta_scale*pv(qD,xr)/R[k][1]+weights[k])<ctx.mpf('1e-28')
    def kernel0(P,Q):
        p=interpolate(ordinary+[Q])
        def raw(P):
            x,y=P
            if abs(x-Q[0])<ctx.mpf('1e-25'):
                assert abs(y+Q[1])<ctx.mpf('1e-25')
                return (fp(x)/(2*y)+sum(k*p[k]*x**(k-1) for k in range(1,len(p))))/pv(qD,x)
            return (y+pv(p,x))/(pv(qD,x)*(x-Q[0]))
        return (raw(P)-sum(psi0(P,k)*raw(R[k]) for k in range(2)))*pv(qD,Q[0])/(2*Q[1]*beta(Q))
    initial=[X[0,1]/ctx.sqrt(2),X[1,0]/ctx.sqrt(2)]+[X[1,k+2] for k in range(2)]+[-X[k+2,1] for k in range(2)]+[-X[k+2,0] for k in range(2)]+[X[0,k+2] for k in range(2)]
    ci=sum(initial[6:8]);cj=sum(initial[2:4]);ki=ci*sum(initial[8:10]);kj=cj*sum(initial[4:6])
    recovered=[ctx.sqrt(2)*ki*kj/(ci*cj)*kernel0(Q[0],Q[1]),-ctx.sqrt(2)*ci*cj*kernel0(Q[1],Q[0])]
    recovered += [cj*psi0(Q[1],k) for k in range(2)]+[kj/cj*dual0(Q[1],k) for k in range(2)]
    recovered += [ci*dual0(Q[0],k) for k in range(2)]+[ki/ci*psi0(Q[0],k) for k in range(2)]
    initial_error=max(abs(a-b)/(1+abs(b)) for a,b in zip(recovered,initial))
    print('Repeated-phase algebraic recovery',ctx.nstr(initial_error,8),flush=True)
    assert initial_error<ctx.mpf('1e-25'),initial_error
    curve=Curve({(j,0):-a for j,a in enumerate(odd)}|{(0,2):ctx.one},ctx=ctx)
    periods=curve.periods_kind_1();inv=(2*periods.omega)**-1;tau=periods.tau;char=curve.riemann_constant().characteristic
    @lru_cache(None)
    def abel(P):return (inv*ctx.matrix(curve.abel_map_kind_1(P).value))[0]
    ar=[abel(P) for P in R];z0=sum(ar)-sum(abel(P) for P in ordinary)
    ER=[inv[0,0]/P[1] for P in R];EdR=[-inv[0,0]*fp(P[0])/(2*P[1]**3) for P in R]
    direction=sum(weights[k]*ER[k] for k in range(2))
    @lru_cache(None)
    def jet(z):return rtheta_jet([z],tau,1,char,ctx=ctx)
    th=lambda z:jet(z)[(0,)];dt=lambda z:jet(z)[(1,)]
    Xi=lambda P:sum(weights[k]*ER[k]*dt(abel(P)-ar[k])/th(abel(P)-ar[k]) for k in range(2))
    finite=[-weights[k]*EdR[k]/(2*ER[k])+weights[1-k]*ER[1-k]*dt(ar[k]-ar[1-k])/th(ar[k]-ar[1-k])+regular[k] for k in range(2)]
    bl,bm=cv(bL),cv(bM)
    def evaluate(t):
        zz=z0+ctx.j*t*direction
        def element(base,delta,rate):
            value=base*ctx.exp(-ctx.j*t*rate)*th(zz+delta)/th(zz)*th(z0)/th(z0+delta)
            irate=rate-direction*dt(zz+delta)/th(zz+delta)+direction*dt(zz)/th(zz)
            return value,-ctx.j*irate*value
        def psi(P,k,dual=False,phase=0):
            delta=abel(P)-ar[k];rate=Xi(P)-finite[k]
            return element(dual0(P,k) if dual else psi0(P,k),-delta if dual else delta,(-rate if dual else rate)+phase)
        def omega(P,Q,phase):return element(kernel0(P,Q),abel(P)-abel(Q),Xi(P)-Xi(Q)+phase)
        pairs=[tuple(ctx.sqrt(2)*ki*kj/(ci*cj)*a for a in omega(Q[0],Q[1],bl+bm)),tuple(-ctx.sqrt(2)*ci*cj*a for a in omega(Q[1],Q[0],-bl-bm))]
        pairs += [tuple(cj*a for a in psi(Q[1],k,phase=-bm)) for k in range(2)]
        pairs += [tuple(kj/cj*a for a in psi(Q[1],k,True,phase=bm)) for k in range(2)]
        pairs += [tuple(ci*a for a in psi(Q[0],k,True,phase=-bl)) for k in range(2)]
        pairs += [tuple(ki/ci*a for a in psi(Q[0],k,phase=bl)) for k in range(2)]
        return list(zip(*pairs))
    Tmp=ctx.matrix([[cv(a) for a in row] for row in T.tolist()]);Tmpinv=Tmp**-1
    theta_errors=[];zero_error=ctx.zero
    def matrices(values,rates):
        def assemble(v,velocity=False):
            u,ub=v[:2];vv=v[2:4];vb=v[4:6];w=v[6:8];wb=v[8:10]
            return ctx.matrix([[0 if velocity else -2*bl,ctx.sqrt(2)*u,*wb,0],
              [ctx.sqrt(2)*ub,0 if velocity else -2*bm,*vv,0],
              [-w[0],-vb[0],0 if velocity else bc,0,0],[-w[1],-vb[1],0,0 if velocity else bc,0],
              [0,0,0,0,0 if velocity else bc]])
        return Tmp*assemble(values)*Tmpinv,Tmp*assemble(rates,True)*Tmpinv
    def rhs(Y):
        # Direct Hamilton equations in all seven original canonical pairs,
        # independently of the spectral eigenvectors and theta formulas.
        root=ctx.sqrt(2);u=Y[0,1]/root;ub=Y[1,0]/root
        v=[Y[1,k+2] for k in range(3)];vb=[-Y[k+2,1] for k in range(3)]
        w=[-Y[k+2,0] for k in range(3)];wb=[Y[0,k+2] for k in range(3)]
        out=ctx.matrix(5)
        out[0,1]=-ctx.j*root*((bl+bm)*u+root*sum(vb[k]*wb[k] for k in range(3)))
        out[1,0]=ctx.j*root*((bl+bm)*ub+root*sum(v[k]*w[k] for k in range(3)))
        for k in range(3):
            out[1,k+2]=-ctx.j*((bm+bc)*v[k]+root*ub*wb[k])
            out[k+2,1]=-ctx.j*((bm+bc)*vb[k]+root*u*w[k])
            out[k+2,0]=ctx.j*((bc+bl)*w[k]+root*ub*vb[k])
            out[0,k+2]=ctx.j*((bc+bl)*wb[k]+root*u*v[k])
        return out
    for t in (ctx.zero,ctx.mpf('.07'),ctx.mpf('.2')):
        values,rates=evaluate(t);Y,dY=matrices(values,rates);expected=rhs(Y)
        error=max(abs(a-b)/(1+abs(b)) for a,b in zip(dY,expected));theta_errors.append(error)
        if not t:zero_error=max(abs(Y[i,j]-cv(X0[i,j]))/(1+abs(cv(X0[i,j]))) for i in range(5) for j in range(5))
    assert max(theta_errors)<ctx.mpf('1e-22'),theta_errors
    assert zero_error<ctx.mpf('1e-25'),zero_error
    # Independent matrix RK4 on the original unreduced canonical coordinates.
    exact0=ctx.matrix([[cv(a) for a in row] for row in X0.tolist()])
    def integrate(steps):
        Y=exact0.copy();h=ctx.mpf('.2')/steps
        for _ in range(steps):
            k1=rhs(Y);k2=rhs(Y+h*k1/2);k3=rhs(Y+h*k2/2);k4=rhs(Y+h*k3)
            Y+=h*(k1+2*k2+2*k3+k4)/6
        return Y
    analytic=matrices(*evaluate(ctx.mpf('.2')))[0]
    errors=[]
    for n in (32,64):
        Y=integrate(n);errors.append(max(abs(a-b)/(1+abs(b)) for a,b in zip(Y,analytic)))
    assert errors[1]<ctx.mpf('1e-10') and 12<errors[0]/errors[1]<20,errors
    return {'dps':dps,'L_M_N':[1,1,3],'active_L_M_N':[1,1,2],'normalization_genus':1,
            'inactive_directions':1,'canonical_pairs':7,'all_original_mode_rates':True,
            'launch_origin':'Direct rational canonical modes, not BA samples',
            'source_matrix':[[str(a) for a in row] for row in X0.tolist()],
            'constant_basis':[[str(a) for a in row] for row in T.tolist()],
            'moment_eigenvalues':list(map(str,moment.eigenvals())),
            'forced_linear_factor':str(z-b),'forced_discriminant_square':str((z-b)**2),
            'normalized_quartic':str(normalized.as_expr()),'quartic_squarefree':True,
            'algebraic_initial_recovery_error':str(initial_error),'restored_initial_error':str(zero_error),
            'restored_mode_rate_error':str(max(theta_errors)),
            'rk4_steps':[32,64],'rk4_duration':'.2','rk4_errors':list(map(str,errors)),
            'rk4_refinement_ratio':str(errors[0]/errors[1]),
            'generapy_backend':'Curve periods/Abel maps and public rtheta_jet',
            'generapy_version':importlib.metadata.version('generapy'),
            'scope':'One smooth active elliptic normalization with a genuine inactive direction; not every phase profile.'}


if __name__=='__main__':
    report=check()
    (Path(__file__).parent/'results/ba_repeated_phase_converse_check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
