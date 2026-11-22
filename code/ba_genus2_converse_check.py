"""Reconstruct a generic (L,M,N)=(1,1,2) solution from canonical modes.

The launch data are prescribed directly, not sampled from BA functions.
An exact characteristic polynomial gives a genus-two hyperelliptic model;
adjugate zeroes determine both BA divisors. Generapy supplies periods,
Abel maps and theta jets only after the spectral reconstruction.
"""
from pathlib import Path
from functools import lru_cache
import argparse
import importlib.metadata
import json
import numpy as np
import sympy as sp
from mpmath import mp
from generapy import Curve,rtheta_jet

ROOT=Path(__file__).resolve().parent


def build(dps=40):
    ctx=mp.clone();ctx.dps=dps
    z,nu=sp.symbols('lambda nu');ii=sp.I
    rat=lambda a,b=1:sp.Rational(a,b)
    bl=rat(1,7);bm=-rat(2,9);bn=[rat(1,3),-rat(1,5)]
    u=(7+2*ii)/(10*sp.sqrt(2));ub=(4-3*ii)/(10*sp.sqrt(2))
    v=[(6+ii)/10,(-3+5*ii)/10];vb=[(8-2*ii)/10,(2+4*ii)/10]
    w=[(5-4*ii)/10,(-2+3*ii)/10];wb=[(9+2*ii)/10,(4-ii)/10]
    initial_exact=[u,ub]+v+vb+w+wb
    XS=sp.Matrix([[-2*bl,sp.sqrt(2)*u,*wb],
                  [sp.sqrt(2)*ub,-2*bm,*v],
                  [-w[0],-vb[0],bn[0],0],[-w[1],-vb[1],0,bn[1]]])
    DM=sp.diag(1,-1,0,0);MS=z*sp.eye(4)-XS-nu*DM
    FS=sp.Poly(sp.expand(MS.det()),nu)
    aa,bb,cc=[sp.expand(FS.nth(k)) for k in (2,1,0)]
    disc=sp.Poly(sp.expand(bb**2-4*aa*cc),z)
    assert disc.degree()==6 and sp.gcd(disc,disc.diff()).degree()==0
    anchor=sp.Matrix([0,0,1,1]);adjS=MS.adjugate();ss=sp.Poly(sp.expand((anchor.T*adjS*anchor)[0]),nu)
    ls=sp.expand(aa*ss.nth(1)-ss.nth(2)*bb);cs=sp.expand(aa*ss.nth(0)-ss.nth(2)*cc)
    resultant=sp.factor(sp.resultant(FS.as_expr(),ss.as_expr(),nu))
    print('Exact spectral discriminant degree',disc.degree(),'adjugate resultant:',sp.factor(resultant),flush=True)
    factors=sp.factor_list(resultant,z)[1]
    candidates=[sp.Poly(factor,z) for factor,power in factors if sp.degree(factor,z)>1]
    divisor_poly=sp.Poly(sp.prod(p.as_expr() for p in candidates),z)
    assert divisor_poly.degree()==6,divisor_poly.degree()
    convert=lambda a:ctx.mpc(str(sp.N(sp.re(a),dps+8)),str(sp.N(sp.im(a),dps+8)))
    def coefficients(poly):return [convert(x) for x in sp.Poly(poly,z).all_coeffs()]
    av,bv,cv=map(coefficients,(aa,bb,cc));dv=coefficients(disc.as_expr())
    lv,sv=map(coefficients,(ls,cs))
    poly=lambda coeff,x:ctx.polyval(coeff,x)
    derivative=lambda coeff:[(len(coeff)-1-j)*a for j,a in enumerate(coeff[:-1])]
    X=ctx.matrix([[convert(a) for a in row] for row in XS.tolist()])
    initial=[convert(a) for a in initial_exact]
    bL,bM=convert(bl),convert(bm);bN=list(map(convert,bn))
    def adjugate(lam,n):
        M=lam*ctx.eye(4)-X-n*ctx.diag([1,-1,0,0])
        return ctx.matrix([[(-1)**(i+j)*ctx.det(ctx.matrix([[M[r,c] for c in range(4) if c!=i] for r in range(4) if r!=j])) for j in range(4)] for i in range(4)])
    branch=ctx.polyroots(dv,maxsteps=1000,extraprec=80)
    e=min(branch,key=lambda x:(float(abs(x)),float(ctx.re(x))))
    # x=1/(lambda-e), Y=y*x^3 sends this branch point to infinity.
    odd=[ctx.zero]*6
    deriv=dv[:]
    for k in range(1,7):
        deriv=derivative(deriv)
        odd[6-k]=poly(deriv,e)/ctx.factorial(k)
    f=lambda x:ctx.polyval(list(reversed(odd)),x)
    fp=lambda x:sum(j*odd[j]*x**(j-1) for j in range(1,6))
    def point(lam,y):
        x=1/(lam-e);Y=y*x**3
        # Use the point on the reconstructed odd model, retaining its sheet.
        root=ctx.sqrt(f(x));Y=min((root,-root),key=lambda a:abs(a-Y))
        return (x,Y)
    ordinary=[];dual=[];root_errors=[]
    for lam in ctx.polyroots(coefficients(divisor_poly.as_expr()),maxsteps=1000,extraprec=80):
        n=-poly(sv,lam)/poly(lv,lam);y=2*poly(av,lam)*n+poly(bv,lam);A=adjugate(lam,n)
        col=A*ctx.matrix([0,0,1,1]);row=ctx.matrix([[0,0,1,1]])*A
        col_size=max(abs(t) for t in col);row_size=max(abs(t) for t in row)
        root_errors.append(min(col_size,row_size)/(1+max(col_size,row_size)))
        (ordinary if col_size<row_size else dual).append(point(lam,y))
    assert len(ordinary)==len(dual)==3,(ordinary,dual)
    R=[];resnu=[];regnu=[]
    for b in bN:
        B=poly(bv,b);ap=poly(derivative(av),b);x=1/(b-e)
        R.append(point(b,-B))
        rl=-B/ap;resnu.append(-rl*x*x)
        yp=poly(derivative(dv),b)/(-2*B)
        c0=(-poly(derivative(bv),b)+yp)/(2*ap)+B*poly(derivative(derivative(av)),b)/(2*ap*ap)
        regnu.append(c0-rl*x)
    # At x=0, nu/lambda=+1 in the L group and -1 in the M group.
    qroot=ctx.sqrt(odd[0]);Q=[(ctx.zero,-qroot),(ctx.zero,qroot)]
    # Reusable ascending-order polynomial arithmetic on the odd model.
    polyval=lambda c,x:ctx.polyval(list(reversed(c)),x)
    def mul(a,b):
        result=[ctx.zero]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):result[i+j]+=x*y
        return result
    def roots_poly(points):
        out=[ctx.one]
        for x,y in points:out=mul(out,[-x,1])
        return out
    def interpolate(points):
        return list(ctx.lu_solve(ctx.matrix([[x**j for j in range(len(points))] for x,y in points]),ctx.matrix([y for x,y in points])))
    qD=roots_poly(ordinary);pD=interpolate(ordinary)
    qDual=roots_poly(dual);pDual=interpolate(dual)
    def basis(divisor,q,p):
        fn=lambda P:(P[1]+polyval(p,P[0]))/polyval(q,P[0])
        fR=[fn(P) for P in R]
        return lambda P,k:((fn(P)-fR[1])/(fR[0]-fR[1]) if k==0 else (fn(P)-fR[0])/(fR[1]-fR[0]))
    psi0=basis(ordinary,qD,pD);dual0=basis(dual,qDual,pDual)
    rows=[];targets=[]
    for k,(x,y) in enumerate(R):
        rows.append([x**j for j in range(6)]+[-y*x**j for j in range(3)]);targets.append(0)
        rows.append([0]+[j*x**(j-1) for j in range(1,6)]+[-(j*x**(j-1)*y if j else 0)-x**j*fp(x)/(2*y) for j in range(3)]);targets.append(0)
        other=R[1-k][0]
        rows.append([x**j for j in range(6)]+[y*x**j for j in range(3)]);targets.append(resnu[k]*y*(x-other)**2)
    for x,y in ordinary:rows.append([x**j for j in range(6)]+[y*x**j for j in range(3)]);targets.append(0)
    beta_coeff=list(ctx.lu_solve(ctx.matrix(rows),ctx.matrix(targets)))
    beta=lambda P:(polyval(beta_coeff[:6],P[0])+P[1]*polyval(beta_coeff[6:],P[0]))/(polyval(roots_poly(R),P[0])**2*P[1])
    def kernel0(P,T):
        p=interpolate(ordinary+[T])
        def fq(P):
            x,y=P
            if abs(x-T[0])<ctx.mpf('1e-25'):
                assert abs(y+T[1])<ctx.mpf('1e-25')
                pd=[j*p[j] for j in range(1,len(p))]
                return (fp(x)/(2*y)+polyval(pd,x))/polyval(qD,x)
            return (y+polyval(p,x))/(polyval(qD,x)*(x-T[0]))
        return (fq(P)-psi0(P,0)*fq(R[0])-psi0(P,1)*fq(R[1]))*polyval(qD,T[0])/(2*T[1]*beta(T))
    ci=sum(initial[6:8]);cj=sum(initial[2:4]);cbi=sum(initial[8:10]);cbj=sum(initial[4:6])
    ki=ci*cbi;kj=cj*cbj
    # Recover the original initial modes without any theta evaluations.
    reconstructed=[ctx.sqrt(2)*ki*kj/(ci*cj)*kernel0(Q[0],Q[1]),-ctx.sqrt(2)*ci*cj*kernel0(Q[1],Q[0])]
    reconstructed += [cj*psi0(Q[1],k) for k in range(2)]+[kj/cj*dual0(Q[1],k) for k in range(2)]
    reconstructed += [ci*dual0(Q[0],k) for k in range(2)]+[ki/ci*psi0(Q[0],k) for k in range(2)]
    initial_error=max(abs(a-b)/(1+abs(b)) for a,b in zip(reconstructed,initial))
    weight_error=max(abs(-beta(Q[0])-ki)/(1+abs(ki)),abs(beta(Q[1])-kj)/(1+abs(kj)))
    print('Algebraic initial recovery:',ctx.nstr(initial_error,8),'weight recovery:',ctx.nstr(weight_error,8),flush=True)
    assert initial_error<ctx.mpf('1e-25') and weight_error<ctx.mpf('1e-25')
    # Check the reconstructed sections and kernel at independent curve points.
    probes=[point(lam,ctx.sqrt(poly(dv,lam))) for lam in (ctx.mpc('.9','.3'),ctx.mpc('-.7','.5'))]
    def coordinate(P):
        x,Y=P;lam=e+1/x;y=Y/x**3;n=(-poly(bv,lam)+y)/(2*poly(av,lam))
        return lam,n
    eig=[];reference_errors=[]
    for P in probes:
        lam,n=coordinate(P);A=adjugate(lam,n);a=ctx.matrix([0,0,1,1]);S=(a.T*A*a)[0]
        right=A*a/S;left=a.T*A/S;eig.append((n,left,right))
        for k in range(2):
            reference_errors.extend([abs(psi0(P,k)-left[k+2])/(1+abs(left[k+2])),abs(dual0(P,k)-right[k+2])/(1+abs(right[k+2]))])
    n1,l1,_=eig[0];n2,_,r2=eig[1]
    expected=-(l1*r2)[0]/(n1-n2)
    reference_errors.append(abs(kernel0(probes[0],probes[1])-expected)/(1+abs(expected)))
    assert max(reference_errors)<ctx.mpf('1e-25'),max(reference_errors)
    print('Computing periods and BA evolution from the recovered curve.',flush=True)
    curve=Curve({**{(j,0):-c for j,c in enumerate(odd)},(0,2):ctx.one},ctx=ctx)
    periods=curve.periods_kind_1();inv=(2*periods.omega)**-1
    tau=(periods.tau+periods.tau.T)/2;char=curve.riemann_constant().characteristic
    @lru_cache(None)
    def abel(P):return tuple(inv*ctx.matrix(curve.abel_map_kind_1(P).value))
    ar=[abel(P) for P in R];ad=[abel(P) for P in ordinary]
    z0=tuple(sum(a[j] for a in ar)-sum(a[j] for a in ad) for j in range(2))
    ER=[tuple(inv*ctx.matrix([1/y,x/y])) for x,y in R]
    EdR=[tuple(inv*ctx.matrix([-fp(x)/(2*y**3),1/y-x*fp(x)/(2*y**3)])) for x,y in R]
    weights=[-a for a in resnu]
    direction=tuple(sum(weights[k]*ER[k][j] for k in range(2)) for j in range(2))
    @lru_cache(None)
    def jet(z):return rtheta_jet(list(z),tau,1,char,ctx=ctx)
    th=lambda z:jet(tuple(z))[(0,0)]
    grad=lambda z:[jet(tuple(z))[(1,0)],jet(tuple(z))[(0,1)]]
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
    add=lambda a,b:tuple(x+y for x,y in zip(a,b))
    grad0=grad((ctx.zero,ctx.zero))
    @lru_cache(None)
    def Xi(P):return sum(weights[k]*dot(ER[k],grad(sub(abel(P),ar[k])))/th(sub(abel(P),ar[k])) for k in range(2))
    finite=[]
    for k in range(2):
        other=1-k;delta=sub(ar[k],ar[other])
        finite.append(-weights[k]*dot(EdR[k],grad0)/(2*dot(ER[k],grad0))+weights[other]*dot(ER[other],grad(delta))/th(delta)+regnu[k])
    def evaluate(t,high_precision=False):
        t=ctx.mpc(t);z=tuple(z0[j]+ctx.j*t*direction[j] for j in range(2))
        def element(base,delta,rate,phase=0):
            zp=add(z,delta)
            value=base*ctx.exp(-ctx.j*t*(rate+phase))*th(zp)/th(z)*th(z0)/th(add(z0,delta))
            irate=rate+phase-dot(direction,grad(zp))/th(zp)+dot(direction,grad(z))/th(z)
            return value,-ctx.j*irate*value
        def psi(P,k,partner=False,phase=0):
            delta=sub(abel(P),ar[k]);rate=Xi(P)-finite[k]
            return element(dual0(P,k) if partner else psi0(P,k),tuple(-x for x in delta) if partner else delta,-rate if partner else rate,phase)
        def omega(P,T,phase=0):return element(kernel0(P,T),sub(abel(P),abel(T)),Xi(P)-Xi(T),phase)
        pairs=[]
        pairs.append(tuple(ctx.sqrt(2)*ki*kj/(ci*cj)*x for x in omega(Q[0],Q[1],bL+bM)))
        pairs.append(tuple(-ctx.sqrt(2)*ci*cj*x for x in omega(Q[1],Q[0],-bL-bM)))
        pairs.extend(tuple(cj*x for x in psi(Q[1],k,phase=-bM)) for k in range(2))
        pairs.extend(tuple(kj/cj*x for x in psi(Q[1],k,True,phase=bM)) for k in range(2))
        pairs.extend(tuple(ci*x for x in psi(Q[0],k,True,phase=-bL)) for k in range(2))
        pairs.extend(tuple(ki/ci*x for x in psi(Q[0],k,phase=bL)) for k in range(2))
        values,rates=map(list,zip(*pairs))
        if high_precision:return values,rates
        return np.array(list(map(complex,values))),np.array(list(map(complex,rates)))
    def rhs_mp(values):
        u,ub=values[:2];v=values[2:4];vb=values[4:6];w=values[6:8];wb=values[8:10];sq=ctx.sqrt(2)
        rates=[(bL+bM)*u+sq*sum(vb[k]*wb[k] for k in range(2)),-(bL+bM)*ub-sq*sum(v[k]*w[k] for k in range(2))]
        rates += [(bM+bN[k])*v[k]+sq*ub*wb[k] for k in range(2)]
        rates += [-(bM+bN[k])*vb[k]-sq*u*w[k] for k in range(2)]
        rates += [(bN[k]+bL)*w[k]+sq*ub*vb[k] for k in range(2)]
        rates += [-(bN[k]+bL)*wb[k]-sq*u*v[k] for k in range(2)]
        return [-ctx.j*x for x in rates]
    errors=[]
    for t in (0,.1,.4,.8,1.5):
        values,rates=evaluate(t,True);expected=rhs_mp(values)
        errors.append(max(abs(a-b)/(1+abs(b)) for a,b in zip(rates,expected)))
        print('Recovered BA time',t,'rate error',ctx.nstr(errors[-1],8),flush=True)
    assert max(errors)<ctx.mpf('1e-23'),max(errors)
    report={'dps':dps,'L_M_N':[1,1,2],'genus':2,'generic_network_genus':2,
            'launch_origin':'Direct independent canonical modes; no BA launch data',
            'source_modes':list(map(str,initial_exact)),'source_phases':list(map(str,[bl,bm]+bn)),
            'characteristic_polynomial':str(FS.as_expr()),'discriminant':str(disc.as_expr()),
            'exact_discriminant_squarefree':True,'algebraic_initial_recovery_error':str(initial_error),
            'residue_weight_recovery_error':str(weight_error),'adjugate_zero_residual':str(max(root_errors)),
            'odd_curve_coefficients':list(map(str,odd)),'branch_sent_to_infinity':str(e),
            'ordinary_divisor':[[str(a) for a in P] for P in ordinary],
            'dual_divisor':[[str(a) for a in P] for P in dual],
            'punctures':[[str(a) for a in P] for P in R],
            'generapy_version':importlib.metadata.version('generapy'),
            'reference_eigenline_and_kernel_error':str(max(reference_errors)),
            'all_ten_mode_rate_error':str(max(errors)),
            'normalized_Abel_direction':list(map(str,direction)),
            'puncture_local_parameters':'s_k=-1/nu',
            'puncture_principal_coefficients':list(map(str,weights)),
            'second_kind_regular_parts':list(map(str,finite)),
            'theta_reference':list(map(str,z0)),
            'period_matrix':[[str(tau[i,j]) for j in range(2)] for i in range(2)],
            'characteristic':[[str(x) for x in row] for row in char]}
    phases=np.array(list(map(complex,[bL,bM]+bN)))
    def rhs(values):
        u,ub=values[:2];v=values[2:4];vb=values[4:6];w=values[6:8];wb=values[8:10];bl,bm=phases[:2];bn=phases[2:];sq=np.sqrt(2)
        return -1j*np.r_[(bl+bm)*u+sq*np.dot(vb,wb),-(bl+bm)*ub-sq*np.dot(v,w),
                        (bm+bn)*v+sq*ub*wb,-(bm+bn)*vb-sq*u*w,
                        (bn+bl)*w+sq*ub*vb,-(bn+bl)*wb-sq*u*v]
    return report,evaluate,rhs


def numerics(report,evaluate,rhs,duration=2.,steps=1600):
    from ba_numeric_examples import rk4
    sparse=np.linspace(0,duration,65)
    direct=np.array([evaluate(float(t))[0] for t in sparse])
    def solve(n):
        h=duration/n;z=direct[0].copy();samples=[z.copy()]
        for j in range(1,n+1):
            z=rk4(rhs,z,h)
            if j%(n//64)==0:samples.append(z.copy())
        return np.array(samples)
    n=((steps+63)//64)*64;coarse=solve(n);fine=solve(2*n)
    error=lambda a,b:float(np.max(abs(a-b)/(1+np.maximum(abs(a),abs(b)))))
    e1,e2=error(coarse,direct),error(fine,direct)
    assert e2<1e-8,(e1,e2)
    phase=[complex(sp.sympify(x)) for x in report['source_phases']]
    bl,bm,b1,b2=phase;DM=np.diag([1,-1,0,0])
    def invariants(values):
        u,ub=values[:2];v=values[2:4];vb=values[4:6];w=values[6:8];wb=values[8:10]
        X=np.array([[-2*bl,np.sqrt(2)*u,*wb],[np.sqrt(2)*ub,-2*bm,*v],
                    [-w[0],-vb[0],b1,0],[-w[1],-vb[1],0,b2]],dtype=complex)
        return np.array([np.trace(np.linalg.matrix_power(X+n*DM,r)) for n in (0,.17) for r in (2,3,4)])
    iv=invariants(direct[0]);ba_drift=max(error(invariants(z),iv) for z in direct);ode_drift=max(error(invariants(z),iv) for z in fine)
    assert max(ba_drift,ode_drift)<1e-9,(ba_drift,ode_drift)
    report['numerics']={'duration':duration,'direct_samples':65,'coarse_steps':n,'fine_steps':2*n,
                        'coarse_RK4_error':e1,'fine_RK4_error':e2,'refinement_ratio':e1/e2,
                        'direct_spectral_drift':ba_drift,'RK4_spectral_drift':ode_drift}
    print('Converse RK4 comparison:',report['numerics'],flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dps',type=int,default=40)
    parser.add_argument('--numerics',action='store_true')
    args=parser.parse_args()
    report,evaluate,rhs=build(args.dps)
    if args.numerics:report=numerics(report,evaluate,rhs)
    path=ROOT/'results'/('ba_genus2_converse_numerics.json' if args.numerics else 'ba_genus2_converse_check.json')
    path.write_text(json.dumps(report,indent=2)+'\n')
