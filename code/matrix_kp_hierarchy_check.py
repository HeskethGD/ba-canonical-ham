"""Compare the second matrix-KP time with canonical spectral flows, genus two.

Uses released Generapy and the existing two-puncture fixture. Paper sources
are unchanged. All derivatives of BA sections use complete theta jets.
"""
from pathlib import Path
from inspect import getclosurevars
import json
from generapy import rtheta_jet
import multipoint_BA_genus2_check as fixture

ROOT = Path(__file__).resolve().parent


def check_source_brackets():
    """Exact ambient-mode Poisson check; no KP field bracket is assumed."""
    import sympy as sp
    u=sp.symbols('u:4');ub=sp.symbols('ub:4')
    v=sp.symbols('v:4');vb=sp.symbols('vb:4')
    w=sp.symbols('w:4');wb=sp.symbols('wb:4')
    c=sp.symbols('c:4',nonzero=True);e=sp.symbols('e:2',nonzero=True)
    kap=sp.symbols('kap:4',nonzero=True)
    pairs=list(zip(u,ub))+list(zip(v,vb))+list(zip(w,wb))
    def bracket(f,g):
        return sp.simplify(-sp.I*sum(sp.diff(f,q)*sp.diff(g,p)-sp.diff(f,p)*sp.diff(g,q) for q,p in pairs))
    ordinary=[];dual=[]
    for a in range(4):
        for k in range(2):
            if a<2:
                ordinary.append(c[a]*e[k]*wb[2*k+a]/kap[a])
                dual.append(w[2*k+a]/(c[a]*e[k]))
            else:
                ordinary.append(e[k]*v[2*(a-2)+k]/c[a])
                dual.append(c[a]*vb[2*(a-2)+k]/(kap[a]*e[k]))
    checks=0
    for a in range(8):
        for b in range(8):
            expected=(sp.I if a//2<2 else -sp.I)/kap[a//2] if a==b else 0
            assert sp.simplify(bracket(ordinary[a],dual[b])-expected)==0
            assert bracket(ordinary[a],ordinary[b])==0
            assert bracket(dual[a],dual[b])==0
            checks+=3
    return {'exact_symbolic_bracket_checks':checks,
            'source_pairs':8,'additional_pair_composite_pairs':4,
            'scope':'Ambient finite canonical coordinates, not a KP field-Poisson reduction.'}


def run():
    source_brackets=check_source_brackets()
    evaluate, _, metadata = fixture.build(32)
    env = getclosurevars(evaluate).nonlocals
    ctx = evaluate.ctx
    Q, eR, edR = (env[k] for k in ('Q', 'eR', 'edR'))
    R=[tuple(map(ctx.mpc,p)) for p in metadata['punctures']]
    abel, ar, z0 = (env[k] for k in ('abel', 'ar', 'z0'))
    th, grad, dot = (env[k] for k in ('th', 'grad', 'dot'))
    jetenv = getclosurevars(getclosurevars(th).nonlocals['jet'].__wrapped__).nonlocals
    tau, char = jetenv['tau'], jetenv['char']
    direction, kap = env['direction'], env['kap']
    f = lambda x: 4*x*(x*x-1)*(x*x-4)
    fp = lambda x: 20*x**4-60*x*x+16
    poly = lambda a,x: ctx.polyval(list(reversed(a)),x)

    def rootpoly(points, orders):
        a = [ctx.one]
        for (x,_), order in zip(points, orders):
            for _ in range(order):
                b = [ctx.zero]*(len(a)+1)
                for j,v in enumerate(a): b[j]-=x*v; b[j+1]+=v
                a=b
        return a

    def basis(P, derivative=False):
        x,y=P
        if not derivative: return [1,x,x*x,x**3,y,x*y]
        return [0,1,2*x,3*x*x,fp(x)/(2*y),y+x*fp(x)/(2*y)]

    def rational(points,orders):
        den=rootpoly(points,orders);rows=[];targets=[]
        for (x,y),order in zip(points,orders):
            rows.append(basis((x,-y)));targets.append(0)
            if order==2: rows.append(basis((x,-y),True));targets.append(0)
        for P in R:
            others=[T for T in points if T!=P]
            otherorders=[n for T,n in zip(points,orders) if T!=P]
            rows.append(basis(P));targets.append(poly(rootpoly(others,otherorders),P[0]))
        assert len(rows)==6
        co=list(ctx.lu_solve(ctx.matrix(rows),ctx.matrix(targets)))
        value=lambda P: (poly(co[:4],P[0])+P[1]*poly(co[4:],P[0]))/poly(den,P[0])
        def local(P,order):
            others=[T for T in points if T!=P]
            otherorders=[n for T,n in zip(points,orders) if T!=P]
            other=rootpoly(others,otherorders);x,y=P
            def regular(s):
                xx=x+s;yy=y*ctx.sqrt(f(xx)/f(x))
                return (poly(co[:4],xx)+yy*poly(co[4:],xx))/poly(other,xx)
            return list(ctx.taylor(regular,0,order+1))
        return value,local

    hL,hlocal=rational(R+Q[:2],[1,1,1,1])
    correction,flocal=rational(R,[2,2])
    hcoeff=[hlocal(P,1) for P in R]
    fcoeff=[flocal(P,2) for P in R]
    a=[v[1] for v in hcoeff]
    constant2=[v[1]**2+2*v[2] for v in hcoeff]
    weights=[2*a[k]-fcoeff[k][1] for k in range(2)]
    dir2=[sum(edR[k][j]+2*a[k]*eR[k][j] for k in range(2)) for j in range(2)]
    reduced_direction=[sum(weights[k]*eR[k][j] for k in range(2)) for j in range(2)]
    assert max(abs(x-y) for x,y in zip(dir2,reduced_direction))<ctx.mpf('1e-27')
    AL=ctx.diag(a)
    grad0=grad((ctx.zero,ctx.zero))
    def xi(P):
        return sum(weights[k]*dot(eR[k],grad(tuple(x-y for x,y in zip(abel(P),ar[k]))))
                   /th(tuple(x-y for x,y in zip(abel(P),ar[k]))) for k in range(2))
    xi0=[]
    for k in range(2):
        other=1-k;delta=tuple(x-y for x,y in zip(ar[k],ar[other]))
        xi0.append(-weights[k]*dot(edR[k],grad0)/(2*dot(eR[k],grad0))
                   +weights[other]*dot(eR[other],grad(delta))/th(delta))

    def hessian(logjet):
        value=logjet[(0,0)];g=ctx.matrix([logjet[(1,0)],logjet[(0,1)]])/value
        H=ctx.matrix([[logjet[(2,0)],logjet[(1,1)]],
                      [logjet[(1,1)],logjet[(0,2)]]])/value-g*g.T
        return H,g

    norm=lambda v:max(abs(x) for x in v)
    error=lambda v,w:norm(v-w)/(1+norm(w))
    D=ctx.diag([1,1,-1,-1,0,0])
    higher=evaluate.higher_flow
    columns=ctx.matrix([[direction[0],higher['direction'][0]],
                        [direction[1],higher['direction'][1]]])
    combination=ctx.lu_solve(columns,ctx.matrix(dir2))
    source_errors=[];dual_errors=[];kernel_errors=[];flow_errors=[]
    phase_rates=None
    probes=[(x,ctx.sqrt(f(x))) for x in map(ctx.mpc,['.9+.9j','2.1+.65j'])]
    for t in map(ctx.mpf,['0','.017','.063']):
        values,base_rates,(psi,omega)=evaluate(t,high_precision=True)
        z=tuple(z0[j]+ctx.j*t*direction[j] for j in range(2))
        H0,g0=hessian(rtheta_jet(list(z),tau,2,char,ctx=ctx))
        E=ctx.diag([ctx.exp(-ctx.j*t*a[k]) for k in range(2)])
        def derivatives(P,dual=False):
            p=ctx.matrix(2,1);px=ctx.matrix(2,1);pxx=ctx.matrix(2,1);py=ctx.matrix(2,1)
            for k in range(2):
                delta=tuple(((-1 if dual else 1)*(x-y)) for x,y in zip(abel(P),ar[k]))
                zd=tuple(x+y for x,y in zip(z,delta))
                H,g=hessian(rtheta_jet(list(zd),tau,2,char,ctx=ctx))
                val,rate=psi(P,k,dual)
                sign=-1 if dual else 1
                shift=sign*a[k];factor=E[k,k]**sign
                p[k]=factor*val;px[k]=factor*(rate+shift*val)
                curvature=(ctx.matrix([list(direction)])*(H-H0)*ctx.matrix(direction))[0]
                rawxx=rate*rate/val+curvature*val
                pxx[k]=factor*(rawxx+2*shift*rate+shift**2*val)
                logy=sign*(xi(P)-xi0[k]+correction(P)-fcoeff[k][2]+constant2[k])
                logy-=dot(reduced_direction,list(g-g0))
                py[k]=p[k]*logy
            return p,px,pxx,py

        V=ctx.matrix(2)
        for i in range(2):
            p=derivatives(Q[i])[0];q=derivatives(Q[i],True)[0]
            V+=2*kap[i]*p*q.T
        for P in Q+probes:
            p,px,pxx,py=derivatives(P)
            q,qx,qxx,qy=derivatives(P,True)
            source_errors.append(error(py,pxx+2*V*p))
            dual_errors.append(error(qy,-qxx-2*V.T*q))
        def kernel_rate(P,Qpoint):
            delta=tuple(x-y for x,y in zip(abel(P),abel(Qpoint)))
            zd=tuple(x+y for x,y in zip(z,delta));_,g=hessian(rtheta_jet(list(zd),tau,2,char,ctx=ctx))
            return (xi(P)-xi(Qpoint)+correction(P)-correction(Qpoint)
                    -dot(reduced_direction,list(g-g0)))*omega(P,Qpoint)[0]
        for P,Qpoint in [(probes[0],probes[1]),(Q[0],Q[2])]:
            p,px,_,_=derivatives(P);q,qx,_,_=derivatives(Qpoint,True)
            expected=(q.T*px-qx.T*p)[0]
            kernel_errors.append(abs(kernel_rate(P,Qpoint)-expected)/(1+abs(expected)))
        yrates=[]
        for name in ('u','ub','v','vb','w','wb'):
            for i in range(2):
                for j in range(2):
                    if name=='u':lograte=kernel_rate(Q[i],Q[j+2])/omega(Q[i],Q[j+2])[0]
                    elif name=='ub':lograte=kernel_rate(Q[j+2],Q[i])/omega(Q[j+2],Q[i])[0]
                    else:
                        P,k,dual=((Q[i+2],j,False) if name=='v' else
                                  (Q[i+2],j,True) if name=='vb' else
                                  (Q[j],i,True) if name=='w' else (Q[j],i,False))
                        pair=derivatives(P,dual);lograte=pair[3][k]/pair[0][k]
                    yrates.append(lograte*values[len(yrates)])
        X=evaluate.pencil(values)
        G=(X+higher['shift']*D)**3;field=D*G-G*D
        def unpack(matrix):
            out=[]
            for name in ('u','ub','v','vb','w','wb'):
                for i in range(2):
                    for j in range(2):
                        if name=='u':v=matrix[i,j+2]/ctx.sqrt(2)
                        elif name=='ub':v=matrix[j+2,i]/ctx.sqrt(2)
                        elif name=='v':v=matrix[i+2,j+4]
                        elif name=='vb':v=-matrix[j+4,i+2]
                        elif name=='w':v=-matrix[i+4,j]
                        else:v=matrix[j,i+4]
                        out.append(v)
            return ctx.matrix(out)
        residual=ctx.matrix(yrates)-combination[0]*ctx.j*ctx.matrix(base_rates)-combination[1]*unpack(field)
        incidence=ctx.matrix(24,5)
        discrete=[]
        row=0
        for name in ('u','ub','v','vb','w','wb'):
            for i in range(2):
                for j in range(2):
                    if name in ('u','ub'):left,right=i,j+2
                    elif name in ('v','vb'):left,right=i+2,j+4
                    else:left,right=i+4,j
                    sign=-1 if name in ('ub','vb','wb') else 1
                    if left<5:incidence[row,left]=sign*values[row]
                    if right<5:incidence[row,right]=-sign*values[row]
                    discrete.append([sign*(int(left==k)-int(right==k)) for k in range(5)])
                    row+=1
        if phase_rates is None:
            # Pick five independent rows exactly by rational incidence structure.
            import sympy as sp
            independent=list(sp.Matrix(discrete).T.rref()[1])
            phase_rates=ctx.lu_solve(ctx.matrix([[incidence[i,j] for j in range(5)] for i in independent]),
                                     ctx.matrix([residual[i] for i in independent]))
        flow_errors.append(error(residual,incidence*phase_rates))
    output={'genus':2,'L_M_N':[2,2,2],'generapy_version':metadata['generapy_version'],
            'dps':32,'second_KP_source_error':str(max(source_errors)),
            'second_KP_dual_source_error':str(max(dual_errors)),
            'second_KP_kernel_error':str(max(kernel_errors)),
            'canonical_spectral_and_gauge_flow_error':str(max(flow_errors)),
            'cubic_and_fourth_trace_coefficients':list(map(str,combination)),
            'vertex_gauge_rates_last_fixed_zero':list(map(str,phase_rates))+['0'],
            'reduced_simple_pole_weights':list(map(str,weights)),
            'KP_local_time':'h_L(P)^2 at both punctures; geometric derivative i*d/dt2',
            'source_bracket_check':source_brackets,
            'Poisson_scope':'Canonical combination on this fixed spectral level; no KP field-bracket pullback asserted.'}
    print(json.dumps(output,indent=2),flush=True)
    assert max(source_errors+dual_errors+kernel_errors+flow_errors)<ctx.mpf('1e-20')
    (ROOT/'results/matrix_kp_hierarchy_check.json').write_text(json.dumps(output,indent=2)+'\n')
    return output


if __name__=='__main__':run()
