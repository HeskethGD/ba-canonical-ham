"""Two-puncture, genus-two L=M=N=2 BA construction.

Public Generapy Abel maps and complete theta jets supply the transcendental
part; rational functions on the curve supply reference sections and weights.
"""
from pathlib import Path
from functools import lru_cache
import argparse
import json
import importlib.metadata
import time
import numpy as np
from mpmath import mp
from generapy import Curve, rtheta_jet

ROOT=Path(__file__).resolve().parent


def build(dps=32):
    ctx=mp.clone();ctx.dps=dps
    f=lambda x:4*x*(x*x-1)*(x*x-4)
    fp=lambda x:20*x**4-60*x*x+16
    def point(x,y=None):
        x=ctx.mpc(x);return (x,2*ctx.sqrt(x*(x*x-1)*(x*x-4)) if y is None else ctx.mpc(y))
    R=list(map(point,['.25','2.4']))
    D=list(map(point,['.1+.15j','.7+.25j','1.4+.2j']))
    Q=list(map(point,['.4','2.6','.8','3.3']))
    poly=lambda a,x:ctx.polyval(list(reversed(a)),x)
    def deriv(a):return [k*a[k] for k in range(1,len(a))]
    def multiply(a,b):
        out=[ctx.zero]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):out[i+j]+=x*y
        return out
    def from_roots(points):
        p=[ctx.one]
        for x,_ in points:p=multiply(p,[-x,1])
        return p
    def interp(points):
        return list(ctx.lu_solve(ctx.matrix([[x**j for j in range(len(points))] for x,y in points]),ctx.matrix([y for x,y in points])))
    qD=from_roots(D);pD=interp(D)
    fun=lambda P:(P[1]+poly(pD,P[0]))/poly(qD,P[0])
    refR=[fun(P) for P in R]
    ordinary=lambda P,k:((fun(P)-refR[1])/(refR[0]-refR[1]) if k==0 else (fun(P)-refR[0])/(refR[1]-refR[0]))
    # beta dx has zeros D+Ddual and double-pole coefficient -1 at R.
    rows=[];rhs=[]
    for x,y in R:
        rows.append([x**j for j in range(6)]+[-y*x**j for j in range(3)]);rhs.append(0)
        rows.append([0]+[j*x**(j-1) for j in range(1,6)]+[-(j*x**(j-1)*y if j else 0)-x**j*fp(x)/(2*y) for j in range(3)]);rhs.append(0)
        other=R[1] if (x,y)==R[0] else R[0]
        rows.append([x**j for j in range(6)]+[y*x**j for j in range(3)]);rhs.append(-y*(x-other[0])**2)
    for x,y in D:
        rows.append([x**j for j in range(6)]+[y*x**j for j in range(3)]);rhs.append(0)
    coeff=list(ctx.lu_solve(ctx.matrix(rows),ctx.matrix(rhs)));A=coeff[:6];B=coeff[6:]
    qR=from_roots(R);den=multiply(qR,qR)
    beta=lambda P:(poly(A,P[0])+poly(B,P[0])*P[1])/(poly(den,P[0])*P[1])
    AA=multiply(A,A);BBf=multiply(multiply(B,B),[0,16,0,-20,0,4]);num=[AA[j]-(BBf[j] if j<len(BBf) else 0) for j in range(len(AA))]
    divisor=multiply(den,qD);rem=num[:];quot=[ctx.zero]*4
    for j in range(3,-1,-1):
        quot[j]=rem[j+len(divisor)-1]/divisor[-1]
        for k,a in enumerate(divisor):rem[j+k]-=quot[j]*a
    roots=ctx.polyroots(list(reversed(quot)),maxsteps=300,extraprec=60)
    Ddual=[point(x,-poly(A,x)/poly(B,x)) for x in roots]
    qDual=from_roots(Ddual);pDual=interp(Ddual)
    fund=lambda P:(P[1]+poly(pDual,P[0]))/poly(qDual,P[0])
    rd=[fund(P) for P in R]
    dual=lambda P,k:((fund(P)-rd[1])/(rd[0]-rd[1]) if k==0 else (fund(P)-rd[0])/(rd[1]-rd[0]))
    def kernel0(P,Q):
        q=Q[0];p=interp(D+[Q])
        fq=lambda T:(T[1]+poly(p,T[0]))/(poly(qD,T[0])*(T[0]-q))
        val=fq(P)-ordinary(P,0)*fq(R[0])-ordinary(P,1)*fq(R[1])
        return val*poly(qD,q)/(2*Q[1]*beta(Q))
    def hfunc(nodes):
        pts=R+nodes;q=from_roots(pts);qd=deriv(q)
        rows=[];rhs=[]
        for x,y in pts:
            rows.append([x**j for j in range(4)]+[-y,-x*y]);rhs.append(0)
        for x,y in R:
            rows.append([x**j for j in range(4)]+[y,x*y]);rhs.append(poly(qd,x))
        c=list(ctx.lu_solve(ctx.matrix(rows),ctx.matrix(rhs)));a,b=c[:4],c[4:]
        value=lambda P:(poly(a,P[0])+P[1]*poly(b,P[0]))/poly(q,P[0])
        residues=[2*y*poly(b,x)/poly(qd,x) for x,y in nodes]
        finite=[]
        qdd=deriv(qd);ad=deriv(a);bd=deriv(b)
        for x,y in R:
            nprime=poly(ad,x)+fp(x)/(2*y)*poly(b,x)+y*poly(bd,x)
            finite.append(nprime/poly(qd,x)-poly(qdd,x)/(2*poly(qd,x)))
        return value,residues,finite,(a,b,q)
    hL,rL,aL,packL=hfunc(Q[:2]);hM,rM,aM,packM=hfunc(Q[2:])
    kap=[beta(P)*r/2 for P,r in zip(Q,rL+rM)]
    bL=[-hM(P)/2 for P in Q[:2]];bM=[hL(P)/2 for P in Q[2:]];bN=[(aM[k]-aL[k])/2 for k in range(2)]
    ell=bL+[-x for x in bM]+[-(aM[k]+aL[k])/2 for k in range(2)]
    print('Computing genus-two periods and Abel maps.',flush=True)
    curve=Curve({(0,2):1,(5,0):-4,(3,0):20,(1,0):-16},ctx=ctx)
    periods=curve.periods_kind_1();inv=(2*periods.omega)**-1
    tau=(periods.tau+periods.tau.T)/2
    char=curve.riemann_constant().characteristic
    @lru_cache(None)
    def abel(P):return tuple(inv*ctx.matrix(curve.abel_map_kind_1(P).value))
    ar=[abel(P) for P in R];ad=[abel(P) for P in D]
    z0=tuple(sum(a[j] for a in ar)-sum(a[j] for a in ad) for j in range(2))
    eR=[tuple(inv*ctx.matrix([1/y,x/y])) for x,y in R]
    edR=[tuple(inv*ctx.matrix([-fp(x)/(2*y**3),1/y-x*fp(x)/(2*y**3)])) for x,y in R]
    direction=tuple(sum(e[j] for e in eR) for j in range(2))
    @lru_cache(None)
    def jet(z):return rtheta_jet(list(z),tau,1,char,ctx=ctx)
    th=lambda z:jet(tuple(z))[(0,0)]
    grad=lambda z:[jet(tuple(z))[(1,0)],jet(tuple(z))[(0,1)]]
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    def sub(a,b):return tuple(x-y for x,y in zip(a,b))
    def add(a,b):return tuple(x+y for x,y in zip(a,b))
    # The infinity Riemann characteristic is nonsingular odd in genus two.
    grad0=grad((ctx.zero,ctx.zero))
    Xi=lambda P:sum(dot(e,grad(sub(abel(P),ar[k])))/th(sub(abel(P),ar[k])) for k,e in enumerate(eR))
    finiteXi=[]
    for k in range(2):
        other=1-k;delta=sub(ar[k],ar[other])
        finiteXi.append(-dot(edR[k],grad0)/(2*dot(eR[k],grad0))+dot(eR[other],grad(delta))/th(delta))
    def local_rational(pack,P,order,extra=0):
        a,b,q=pack;x,y=P
        other=q[:]
        for _ in range(order):
            # Synthetic division by x-x0, in ascending order.
            nxt=[ctx.zero]*(len(other)-1);nxt[-1]=other[-1]
            for j in range(len(nxt)-2,-1,-1):nxt[j]=other[j+1]+x*nxt[j+1]
            other=nxt
        def regular(t):
            xx=x+t;yy=y*ctx.sqrt(f(xx)/f(x))
            return (poly(a,xx)+yy*poly(b,xx))/poly(other,xx)
        return list(ctx.taylor(regular,0,order+extra))
    def higher_flow():
        # Normalize q=lambda^3/(nu-shift) by a rational correction with
        # double poles at the four ordinary points and simple poles at R.
        shift=ctx.mpf('.17');principal=[];qfinite=[]
        for j,P in enumerate(Q):
            local=local_rational(packL if j<2 else packM,P,1,1)
            other=hM if j<2 else hL
            other0=other(P)
            other1=ctx.diff(lambda t:other((P[0]+t,P[1]*ctx.sqrt(f(P[0]+t)/f(P[0])))),0)
            if j<2:hl=local;hm=[0,other0,other1]
            else:hl=[0,other0,other1];hm=local
            lm=(hm[0]-hl[0])/2;nm=-(hl[0]+hm[0])/2
            l0=(hm[1]-hl[1])/2;n0=-(hl[1]+hm[1])/2-shift
            l1=(hm[2]-hl[2])/2;n1=-(hl[2]+hm[2])/2
            aa=lm**3/nm;bb=3*lm**2*l0/nm-lm**3*n0/nm**2
            cc=(3*lm**2*l1+3*lm*l0*l0)/nm-3*lm*lm*l0*n0/nm**2+lm**3*(n0*n0/nm**3-n1/nm**2)
            principal.append((-aa,-bb));qfinite.append(cc)
        denf=multiply(multiply(from_roots(Q),from_roots(Q)),qR)
        def row(P,sign=1,derivative=False):
            x,y=P;y*=sign
            if not derivative:return [x**j for j in range(10)]+[y*x**j for j in range(8)]
            return [0]+[j*x**(j-1) for j in range(1,10)]+[(j*y*x**(j-1) if j else 0)+fp(x)/(2*y)*x**j for j in range(8)]
        mat=[];targets=[]
        for P in Q:
            mat.extend([row(P,-1),row(P,-1,True)]);targets.extend([0,0])
        for P in R:mat.append(row(P,-1));targets.append(0)
        for P,(aa,bb) in zip(Q,principal):
            other=from_roots([T for T in Q if T!=P]);other=multiply(multiply(other,other),qR)
            mat.extend([row(P),row(P,derivative=True)])
            targets.extend([aa*poly(other,P[0]),bb*poly(other,P[0])+aa*poly(deriv(other),P[0])])
        co=list(ctx.lu_solve(ctx.matrix(mat),ctx.matrix(targets)));pack=(co[:10],co[10:],denf)
        rlocal=[local_rational(pack,P,1) for P in R]
        const=-rlocal[0][1]
        outer=[qfinite[j]+local_rational(pack,P,2)[2]+const for j,P in enumerate(Q)]
        central=[a[1]+const for a in rlocal]
        weights=[-a[0] for a in rlocal]
        return dict(weights=weights,phases=outer+central,shift=shift,
                    direction=tuple(sum(weights[k]*eR[k][j] for k in range(2)) for j in range(2)))
    flow4=higher_flow()
    print('Checking two-puncture theta normalization.',flush=True)
    points=Q+[point('.9+.9j'),point('2.1+.65j')]
    for P in points:abel(P);Xi(P)
    def evaluate(t,flow=None,high_precision=False):
        active=direction if flow is None else flow["direction"]
        phases=ell if flow is None else flow["phases"]
        weights=[1,1] if flow is None else flow["weights"]
        xi=lambda P:sum(weights[k]*dot(eR[k],grad(sub(abel(P),ar[k])))/th(sub(abel(P),ar[k])) for k in range(2))
        xi0=[]
        for k in range(2):
            other=1-k;delta=sub(ar[k],ar[other])
            xi0.append(-weights[k]*dot(edR[k],grad0)/(2*dot(eR[k],grad0))+weights[other]*dot(eR[other],grad(delta))/th(delta))
        t=ctx.mpc(t);z=tuple(z0[j]+ctx.j*t*active[j] for j in range(2))
        evaluate.last_theta_abs=abs(th(z))
        def element(base,delta,rate,phase=0):
            zd=add(z,delta)
            value=base*ctx.exp(-ctx.j*t*(rate+phase))*th(zd)/th(z)*th(z0)/th(add(z0,delta))
            lograte=rate+phase-dot(active,grad(zd))/th(zd)+dot(active,grad(z))/th(z)
            return value,lograte*value # second entry is i*d/dt
        def psi(P,k,conjugate=False):
            delta=sub(abel(P),ar[k]);rate=xi(P)-xi0[k]
            if conjugate:return element(dual(P,k),tuple(-x for x in delta),-rate)
            return element(ordinary(P,k),delta,rate)
        def omega(P,Q):return element(kernel0(P,Q),sub(abel(P),abel(Q)),xi(P)-xi(Q))
        def gauged(pair,phase):
            value,rate=pair;factor=ctx.exp(-ctx.j*t*phase)
            return value*factor,(rate+phase*value)*factor
        vals=[];rates=[]
        for name in ('u','ub','v','vb','w','wb'):
            vv=[];rr=[]
            for a in range(2):
                for b in range(2):
                    if name=='u':pair=omega(Q[a],Q[b+2]);scale=ctx.sqrt(2)*kap[a]*kap[b+2];phase=phases[a]-phases[b+2]
                    if name=='ub':pair=omega(Q[b+2],Q[a]);scale=-ctx.sqrt(2);phase=phases[b+2]-phases[a]
                    if name=='v':pair=psi(Q[a+2],b);scale=1;phase=phases[a+2]-phases[b+4]
                    if name=='vb':pair=psi(Q[a+2],b,True);scale=kap[a+2];phase=phases[b+4]-phases[a+2]
                    if name=='w':pair=psi(Q[b],a,True);scale=1;phase=phases[a+4]-phases[b]
                    if name=='wb':pair=psi(Q[b],a);scale=kap[b];phase=phases[b]-phases[a+4]
                    value,rate=gauged(pair,phase);vv.append(scale*value);rr.append(scale*rate)
            vals.extend(vv);rates.extend(rr)
        if high_precision:return vals,[-ctx.j*v for v in rates],(psi,omega)
        return np.array([complex(v) for v in vals]),np.array([complex(-ctx.j*v) for v in rates]),(psi,omega)
    phases=[np.array(list(map(complex,b))) for b in (bL,bM,bN)]
    def rhs(t,values):
        u,ub,v,vb,w,wb=[x.reshape(2,2) for x in np.split(values,6)]
        bl,bm,bn=phases;s=np.sqrt(2)
        return -1j*np.concatenate([a.ravel() for a in (
            (bl[:,None]+bm[None,:])*u+s*wb.T@vb.T,
            -(bl[:,None]+bm[None,:])*ub-s*w.T@v.T,
            (bm[:,None]+bn[None,:])*v+s*ub.T@wb.T,
            -(bm[:,None]+bn[None,:])*vb-s*u.T@w.T,
            (bn[:,None]+bl[None,:])*w+s*vb.T@ub.T,
            -(bn[:,None]+bl[None,:])*wb-s*v.T@u.T)])
    residuals=[];kernels=[]
    for t in (0,.003,.009):
        values,rates,(psi,omega)=evaluate(t);expected=rhs(t,values)
        residuals.append(float(np.max(abs(rates-expected)/(1+abs(expected)))))
        probeP,probeQ=points[-2:];value,rate=omega(probeP,probeQ)
        target=sum(psi(probeP,k)[0]*psi(probeQ,k,True)[0] for k in range(2))
        kernels.append(ctx.nstr(abs(rate-target)/(1+abs(target)),12))
        print('t',t,'rate residual',residuals[-1],'kernel residual',kernels[-1],flush=True)
    def matrices(values):
        return [ctx.matrix([[values[4*k+2*i+j] for j in range(2)] for i in range(2)]) for k in range(6)]
    def pencil(values,rate=False):
        u,ub,v,vb,w,wb=matrices(values);X=ctx.matrix(6)
        for i in range(2):
            if not rate:X[i,i]=-2*bL[i];X[i+2,i+2]=-2*bM[i];X[i+4,i+4]=bN[i]
            for j in range(2):
                X[i,j+2]=ctx.sqrt(2)*u[i,j];X[j+2,i]=ctx.sqrt(2)*ub[i,j]
                X[i+2,j+4]=v[i,j];X[j+4,i+2]=-vb[i,j]
                X[i+4,j]=-w[i,j];X[j,i+4]=wb[i,j]
        return X
    high_errors=[];spectral_errors=[]
    DM=ctx.diag([1,1,-1,-1,0,0]);Kph=ctx.diag([-3*x for x in bL]+[3*x for x in bM]+[0,0])
    for t in (0,.17,.63):
        values,rates,(psi,omega)=evaluate(t,high_precision=True)
        X=pencil(values);field=-(DM*X**2-X**2*DM)+Kph*X-X*Kph
        du=pencil(rates,True)
        high_errors.append(max(abs(ctx.j*du[i,j]-field[i,j])/(1+abs(field[i,j])) for i in range(6) for j in range(6)))
        for P in points[-2:]:
            nu=-(hL(P)+hM(P))/2;lam=(hM(P)-hL(P))/2
            right=ctx.matrix([-kap[i]*omega(Q[i],P)[0] for i in range(2)]+[omega(Q[j+2],P)[0] for j in range(2)]+[psi(P,k,True)[0] for k in range(2)])
            right=ctx.matrix([right[j]*ctx.exp(-ctx.j*t*ell[j]) for j in range(6)])
            left=ctx.matrix([[-omega(P,Q[i])[0] for i in range(2)]+[kap[j+2]*omega(P,Q[j+2])[0] for j in range(2)]+[psi(P,k)[0] for k in range(2)]])
            left=ctx.matrix([[left[j]*ctx.exp(ctx.j*t*ell[j]) for j in range(6)]])
            lx=X+nu*DM
            for residual,reference in ((lx*right-lam*right,lam*right),(left*lx-lam*left,lam*left)):
                spectral_errors.append(max(abs(x) for x in residual)/(1+max(abs(x) for x in reference)))
    higher_errors=[]
    for t in (0,.007,.021):
        values,rates,_=evaluate(t,flow4,True)
        X=pencil(values);gg=(X+flow4['shift']*DM)**3;field=DM*gg-gg*DM
        du=pencil(rates,True)
        higher_errors.append(max(abs(ctx.j*du[i,j]-field[i,j])/(1+abs(field[i,j])) for i in range(6) for j in range(6)))
    det=direction[0]*flow4['direction'][1]-direction[1]*flow4['direction'][0]
    print('Fourth trace rate errors:',list(map(str,higher_errors)),'direction determinant:',det,flush=True)
    metadata={'backend':'Generapy periods, Abel maps and complete first-order theta jets',
              'generapy_version':importlib.metadata.version('generapy'),'dps':dps,'genus':2,'L_M_N':[2,2,2],'curve':'y^2=4*x*(x^2-1)*(x^2-4)',
              'punctures':[[str(x),str(y)] for x,y in R],
              'ordinary_divisor':[[str(x),str(y)] for x,y in D],
              'dual_divisor':[[str(x),str(y)] for x,y in Ddual],
              'theta_reference':[str(x) for x in z0],'Abel_direction':[str(x) for x in direction],
              'phases':[[str(x) for x in b] for b in (bL,bM,bN)],
              'mode_rate_residuals':residuals,'kernel_rate_residuals':kernels,
              'high_precision_mode_rate_error':str(max(high_errors)),
              'spectral_eigenvector_error':str(max(spectral_errors)),
              'fourth_trace_rate_errors':list(map(str,higher_errors)),
              'fourth_trace_direction':list(map(str,flow4['direction'])),
              'Abel_direction_determinant':str(det)}
    metadata['evaluation_points']=[[str(x),str(y)] for x,y in Q]
    metadata['weights']=list(map(str,kap))
    metadata['fourth_trace_puncture_weights']=list(map(str,flow4['weights']))
    metadata['fourth_trace_vertex_phases']=list(map(str,flow4['phases']))
    metadata['beta_numerator_A']=list(map(str,A))
    metadata['beta_numerator_B']=list(map(str,B))
    metadata['generic_network_genus']=7
    metadata['spectral_defect']=5
    assert max(high_errors)<ctx.mpf('1e-20')
    assert max(spectral_errors)<ctx.mpf('1e-20')
    assert max(higher_errors)<ctx.mpf('1e-20')
    assert abs(det)>ctx.mpf('.001')
    evaluate.higher_flow=flow4
    evaluate.pencil=pencil
    evaluate.ctx=ctx
    evaluate.phases=phases
    evaluate.jet_cache_info=jet.cache_info
    return evaluate,rhs,metadata


def run_numerics(evaluate,rhs,data,duration=8.,higher_duration=1.5,steps=8000):
    """Sparse direct BA samples versus dense independent classical RK4."""
    from ba_numeric_examples import rk4
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    started=time.perf_counter()
    D=np.diag([1,1,-1,-1,0,0]);shift=complex(evaluate.higher_flow['shift'])
    def pencil_np(z):
        u,ub,v,vb,w,wb=[a.reshape(2,2) for a in np.split(z,6)]
        bl,bm,bn=evaluate.phases
        return np.block([[-2*np.diag(bl),np.sqrt(2)*u,wb.T],
                         [np.sqrt(2)*ub.T,-2*np.diag(bm),v],
                         [-w,-vb.T,np.diag(bn)]])
    def higher_rhs(t,z):
        X=pencil_np(z);G=np.linalg.matrix_power(X+shift*D,3);rate=-1j*(D@G-G@D)
        return np.concatenate([a.ravel() for a in (rate[:2,2:4]/np.sqrt(2),rate[2:4,:2].T/np.sqrt(2),
                               rate[2:4,4:],-rate[4:,2:4].T,-rate[4:,:2],rate[:2,4:].T)])
    def invariants(z):
        u,ub,v,vb,w,wb=[a.reshape(2,2) for a in np.split(z,6)]
        bl,bm,bn=evaluate.phases
        H=np.sum((bl[:,None]+bm[None,:])*u*ub)+np.sum((bm[:,None]+bn[None,:])*v*vb)+np.sum((bn[:,None]+bl[None,:])*w*wb)
        H+=np.sqrt(2)*(np.trace(u@v@w)+np.trace(ub@vb@wb))
        power_u=u*ub;power_v=v*vb;power_w=w*wb
        charges=np.r_[power_w.sum(axis=0)-power_u.sum(axis=1),power_v.sum(axis=1)-power_u.sum(axis=0),power_w.sum(axis=1)-power_v.sum(axis=0)]
        X=pencil_np(z)
        # Nonzero spectral shifts exercise invariants beyond tr(X^r).
        spectral=[np.trace(np.linalg.matrix_power(X+nu*D,r)) for nu in (.23,.11+.17j) for r in range(2,7)]
        return np.r_[H,charges,spectral]
    scaled=lambda a,b:float(np.max(abs(a-b)/(1+np.maximum(abs(a),abs(b)))))
    runs=[];reports=[]
    for name,end,flow,fun in [('cubic',duration,None,rhs),('fourth_trace',higher_duration,evaluate.higher_flow,higher_rhs)]:
        sparse=np.linspace(0,end,65);dense_t=np.linspace(0,end,801);direct=[];residuals=[];denoms=[]
        print(name,'evaluating 65 sparse direct BA samples.',flush=True)
        for t in sparse:
            z,dz,_=evaluate(float(t),flow);direct.append(z);residuals.append(scaled(dz,fun(t,z)));denoms.append(float(evaluate.last_theta_abs))
        direct=np.array(direct)
        def solve(n):
            h=end/n;z=direct[0].copy();samples=[z.copy()];dense=[z.copy()]
            for j in range(1,n+1):
                z=rk4(lambda z:fun(0,z),z,h)
                if not np.isfinite(z).all():raise ValueError('Nonfinite RK4 solution')
                if j%(n//64)==0:samples.append(z.copy())
                if j%(n//800)==0:dense.append(z.copy())
            return np.array(samples),np.array(dense)
        n=((steps+1599)//1600)*1600
        if flow is not None:n*=4
        coarse,_=solve(n);fine,dense=solve(2*n)
        err_coarse=scaled(coarse,direct);err_fine=scaled(fine,direct)
        initial=invariants(direct[0]);drift_BA=max(scaled(invariants(z),initial) for z in direct);drift_ODE=max(scaled(invariants(z),initial) for z in dense)
        report={'flow':name,'duration':end,'direct_samples':65,'dense_ODE_samples':801,'coarse_steps':n,'fine_steps':2*n,
                'all_24_rates_max_scaled_error':max(residuals),'coarse_RK4_error':err_coarse,'fine_RK4_error':err_fine,
                'coarse_to_fine_error_ratio':err_coarse/err_fine,'direct_integral_drift':drift_BA,'RK4_integral_drift':drift_ODE,
                'minimum_sampled_theta_denominator':min(denoms),'maximum_sampled_mode_size':float(np.max(abs(direct)))}
        assert max(residuals)<1e-10,report
        assert err_fine<1e-6,report
        assert drift_BA<1e-9 and drift_ODE<1e-6,report
        print(json.dumps(report,indent=2),flush=True)
        runs.append((sparse,direct,dense_t,dense));reports.append(report)
    plt.rcParams.update({'font.family':'serif','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(4,2,figsize=(7.1,9.2),layout='constrained')
    def obs(z,row):
        if row==0:return [(z[:,0].real,r'$\mathrm{Re}\,u_{11}$')]
        if row==1:return [(z[:,8].real,r'$\mathrm{Re}\,v_{11}$'),(z[:,9].real,r'$\mathrm{Re}\,v_{12}$')]
        if row==2:return [(z[:,16].real,r'$\mathrm{Re}\,w_{11}$'),(z[:,18].real,r'$\mathrm{Re}\,w_{21}$')]
        return [((z[:,0]*z[:,4]).real,r'$\mathrm{Re}(u_{11}\bar u_{11})$')]
    for col,(ts,z,td,dense) in enumerate(runs):
        for row in range(4):
            ax=axes[row,col]
            for k,((line,label),(markers,_)) in enumerate(zip(obs(dense,row),obs(z,row))):
                ax.plot(td,line,color='black',ls=('-','--')[k],lw=.9,label=label)
                ax.plot(ts,markers,ls='none',marker=('o','s')[k],ms=3,mfc='white',mec='black',mew=.65)
            ax.legend(loc='best',frameon=False,fontsize=8)
            ax.xaxis.set_major_locator(MaxNLocator(5));ax.tick_params(labelsize=8)
            if row==3:ax.set_xlabel(r'$t$' if col==0 else r'$t_4$')
        axes[0,col].set_title(('Cubic Hamiltonian $H$',r'Fourth-trace Hamiltonian $H_{4,\nu_0}$')[col],fontsize=10)
    (ROOT/'figures').mkdir(exist_ok=True)
    target=ROOT/'figures/multipoint_genus2_ba.pdf';fig.savefig(target);fig.savefig(target.with_suffix('.png'),dpi=150);plt.close(fig)
    np.savez_compressed(ROOT/'results/multipoint_genus2_trajectories.npz',**{f'{name}_{key}':value for name,run in zip(('cubic','fourth_trace'),runs) for key,value in zip(('direct_t','direct_modes','RK4_t','RK4_modes'),run)})
    data['numerical_runs']=reports;data['numerical_seconds']=time.perf_counter()-started
    data['figure']='figures/multipoint_genus2_ba.pdf'
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dps',type=int,default=32)
    parser.add_argument('--numerics',action='store_true')
    parser.add_argument('--duration',type=float,default=8.)
    parser.add_argument('--higher-duration',type=float,default=1.5)
    parser.add_argument('--steps',type=int,default=8000)
    parser.add_argument('--output',type=Path,help='Optional JSON result path')
    args=parser.parse_args()
    evaluate,rhs,data=build(args.dps)
    if args.numerics:data=run_numerics(evaluate,rhs,data,args.duration,args.higher_duration,args.steps)
    path=args.output or ROOT/'results'/('multipoint_BA_genus2_numerics.json' if args.numerics else 'multipoint_BA_genus2_check.json')
    path.write_text(json.dumps(data,indent=2)+'\n')
