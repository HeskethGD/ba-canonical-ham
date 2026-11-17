"""Positive-genus, all-groups-nonsingleton check using Generapy theta jets.

L=M=N=2 on an elliptic curve. Checks all 24 canonical mode rates and
the prescribed spectral map. No AHH solution formula is used.
"""
from functools import lru_cache
from pathlib import Path
import json
import importlib.metadata
import mpmath as mp
from generapy import rtheta_jet


def check(higher_orders=()):
    mp.mp.dps=55
    tau=mp.mpc('.23','1.17'); period=[[tau]]
    half=mp.mpf('.5');char=([half],[half])
    @lru_cache(None)
    def jet(z):
        return rtheta_jet([mp.mpc(z)],period,3,char,ctx=mp.mp)
    th=lambda z:jet(z)[(0,)]
    e1=lambda z:jet(z)[(1,)]/jet(z)[(0,)]
    ep=lambda z:jet(z)[(2,)]/jet(z)[(0,)]-e1(z)**2
    phi=lambda z,a:jet(0)[(1,)]*th(z+a)/(th(z)*th(a))
    R=[mp.mpc('.13','.14'),mp.mpc('.41','.29')]
    beta=lambda x:sum(ep(x-r) for r in R)+mp.mpc('.23','-.19')
    def norm(u):
        b=mp.im(u)/mp.im(tau);a=mp.re(u)-mp.re(tau)*b
        return u-mp.floor(a)-mp.floor(b)*tau
    roots=[]
    for a in range(6):
        for b in range(6):
            seed=(mp.mpf(a)+mp.mpf('.31'))/6+tau*(mp.mpf(b)+mp.mpf('.47'))/6
            try:
                root=norm(mp.findroot(beta,(seed,seed+mp.mpf('.013')),maxsteps=50))
            except (ValueError,ZeroDivisionError):
                continue
            if abs(beta(root))>mp.mpf('1e-42'):
                continue
            for candidate in [root,norm(sum(R)-root)]:
                if all(abs(candidate-old)>mp.mpf('1e-35') for old in roots):
                    assert abs(beta(candidate))<mp.mpf('1e-40')
                    roots.append(candidate)
            if len(roots)==4:break
        if len(roots)==4:break
    assert len(roots)==4,len(roots)
    D=[roots[0],roots[2]]
    z0=sum(R)-sum(D)
    fixed=lambda P:mp.fprod(th(P-r) for r in R)/mp.fprod(th(P-d) for d in D)
    fixedp=[jet(0)[(1,)]*th(R[k]-R[1-k])/
            mp.fprod(th(R[k]-d) for d in D) for k in range(2)]
    dualfixed=lambda P:1/(beta(P)*fixed(P))
    chi=lambda P:sum(e1(P-r) for r in R)
    reg=[e1(R[k]-R[1-k]) for k in range(2)]
    nodes=[mp.mpc('.07','.03'),mp.mpc('.57','.08'),
           mp.mpc('.78','.21'),mp.mpc('.26','.63')]
    principal=list(map(mp.mpf,['-.7','-1.3','-.9','-1.1']))
    kap=[beta(nodes[a])*principal[a]/2 for a in range(4)]
    hL=lambda P:chi(P)+sum(principal[a]*e1(P-nodes[a]) for a in range(2))
    hM=lambda P:chi(P)+sum(principal[a]*e1(P-nodes[a]) for a in range(2,4))
    aL=[reg[k]+sum(principal[a]*e1(R[k]-nodes[a]) for a in range(2)) for k in range(2)]
    aM=[reg[k]+sum(principal[a]*e1(R[k]-nodes[a]) for a in range(2,4)) for k in range(2)]
    bL=[-hM(nodes[i])/2 for i in range(2)]
    bM=[hL(nodes[j+2])/2 for j in range(2)]
    bN=[(aM[k]-aL[k])/2 for k in range(2)]
    ellL=bL;ellM=[-x for x in bM]
    ellN=[-(aM[k]+aL[k])/2 for k in range(2)]
    root2=mp.sqrt(2)
    errors=[];spectral_errors=[];kernel_errors=[];bilinear_errors=[]
    initial_ordinary=lambda P:sum(fixed(P)/fixedp[k]*phi(z0,P-R[k]) for k in range(2))
    initial_dual=lambda P:sum(fixedp[k]*dualfixed(P)*phi(z0,R[k]-P) for k in range(2))
    partition_error=max(abs(fn(P)-1) for fn in (initial_ordinary,initial_dual) for P in nodes)
    assert partition_error<mp.mpf('1e-38'),partition_error
    # In genus one retain the original two-puncture time as the single
    # gap-complement direction. The higher-flow correction has principal
    # part -speed/s_k at both punctures.
    def normalize_higher(order,nu0):
        assert order in (3,4)
        theta_linear=jet(0)[(3,)]/(3*jet(0)[(1,)])
        A=[];B=[];qfinite=[]
        for a,node in enumerate(nodes):
            hl0=chi(node)+sum(principal[b]*e1(node-nodes[b]) for b in range(2) if b!=a)
            hm0=chi(node)+sum(principal[b]*e1(node-nodes[b]) for b in range(2,4) if b!=a)
            hl1=sum(ep(node-r) for r in R)+sum(principal[b]*ep(node-nodes[b]) for b in range(2) if b!=a)
            hm1=sum(ep(node-r) for r in R)+sum(principal[b]*ep(node-nodes[b]) for b in range(2,4) if b!=a)
            if a<2:hl1+=principal[a]*theta_linear
            else:hm1+=principal[a]*theta_linear
            lm=(-principal[a]/2 if a<2 else principal[a]/2)
            nm=-principal[a]/2
            l0=(hm0-hl0)/2;n0=-(hl0+hm0)/2-nu0
            l1=(hm1-hl1)/2;n1=-(hl1+hm1)/2
            if order==3:
                A.append(mp.mpc(0));B.append(lm**2/nm)
                qfinite.append(2*lm*l0/nm-lm**2*n0/nm**2)
            else:
                A.append(lm**3/nm)
                B.append(3*lm**2*l0/nm-lm**3*n0/nm**2)
                qfinite.append((3*lm**2*l1+3*lm*l0**2)/nm
                               -3*lm**2*l0*n0/nm**2
                               +lm**3*(n0**2/nm**3-n1/nm**2))
        speed=-sum(B)/2
        constant=-sum(A[a]*ep(R[0]-nodes[a])-B[a]*e1(R[0]-nodes[a]) for a in range(4))
        constant+=speed*e1(R[0]-R[1])
        outer=[]
        for a,node in enumerate(nodes):
            finite=A[a]*theta_linear+constant
            finite+=sum(A[b]*ep(node-nodes[b])-B[b]*e1(node-nodes[b]) for b in range(4) if b!=a)
            finite-=speed*sum(e1(node-r) for r in R)
            outer.append(qfinite[a]+finite)
        central=[sum(A[a]*ep(r-nodes[a])-B[a]*e1(r-nodes[a]) for a in range(4))
                 -speed*e1(r-R[1-k])+constant for k,r in enumerate(R)]
        assert abs(central[0])<mp.mpf('1e-45')
        return {'order':order,'nu0':nu0,'speed':speed,'outer':outer,
                'central':central,'errors':[]}
    higher=[normalize_higher(order,mp.mpc('.17','.09')) for order in higher_orders]
    for time in list(map(mp.mpf,['0','.031','.087'])):
        z=z0+2j*time
        logshift=lambda a:-2*(e1(z+a)-e1(z))
        def omega(P,Q):
            return fixed(P)*dualfixed(Q)*mp.exp(-1j*time*(chi(P)-chi(Q)))*phi(z,P-Q)
        def psi(P,k):
            return fixed(P)/fixedp[k]*mp.exp(-1j*time*(chi(P)-reg[k]))*phi(z,P-R[k])
        def psid(P,k):
            return fixedp[k]*dualfixed(P)*mp.exp(1j*time*(chi(P)-reg[k]))*phi(z,R[k]-P)
        def matrix(values):return mp.matrix(2,2,values)
        def make(value,rate):
            A=mp.matrix(2);dA=mp.matrix(2)
            for a in range(2):
                for b in range(2):
                    base,phase=value(a,b)
                    A[a,b]=mp.exp(-1j*time*phase)*base
                    dA[a,b]=(rate(a,b)+phase)*A[a,b]
            return A,dA
        u,du=make(lambda i,j:(root2*kap[i]*kap[j+2]*omega(nodes[i],nodes[j+2]),ellL[i]-ellM[j]),
                  lambda i,j:chi(nodes[i])-chi(nodes[j+2])+logshift(nodes[i]-nodes[j+2]))
        ub,dub=make(lambda i,j:(-root2*omega(nodes[j+2],nodes[i]),-(ellL[i]-ellM[j])),
                    lambda i,j:chi(nodes[j+2])-chi(nodes[i])+logshift(nodes[j+2]-nodes[i]))
        v,dv=make(lambda j,k:(psi(nodes[j+2],k),ellM[j]-ellN[k]),
                  lambda j,k:chi(nodes[j+2])-reg[k]+logshift(nodes[j+2]-R[k]))
        vb,dvb=make(lambda j,k:(kap[j+2]*psid(nodes[j+2],k),-(ellM[j]-ellN[k])),
                    lambda j,k:-(chi(nodes[j+2])-reg[k])+logshift(R[k]-nodes[j+2]))
        w,dw=make(lambda k,i:(psid(nodes[i],k),ellN[k]-ellL[i]),
                  lambda k,i:-(chi(nodes[i])-reg[k])+logshift(R[k]-nodes[i]))
        wb,dwb=make(lambda k,i:(kap[i]*psi(nodes[i],k),-(ellN[k]-ellL[i])),
                    lambda k,i:chi(nodes[i])-reg[k]+logshift(nodes[i]-R[k]))
        for i in range(2):
            for j in range(2):
                expected=(bL[i]+bM[j])*u[i,j]+root2*sum(vb[j,k]*wb[k,i] for k in range(2))
                expectedb=-(bL[i]+bM[j])*ub[i,j]-root2*sum(v[j,k]*w[k,i] for k in range(2))
                errors.extend([abs(du[i,j]-expected)/(1+abs(expected)),abs(dub[i,j]-expectedb)/(1+abs(expectedb))])
        for j in range(2):
            for k in range(2):
                expected=(bM[j]+bN[k])*v[j,k]+root2*sum(ub[i,j]*wb[k,i] for i in range(2))
                expectedb=-(bM[j]+bN[k])*vb[j,k]-root2*sum(u[i,j]*w[k,i] for i in range(2))
                errors.extend([abs(dv[j,k]-expected)/(1+abs(expected)),abs(dvb[j,k]-expectedb)/(1+abs(expectedb))])
        for k in range(2):
            for i in range(2):
                expected=(bN[k]+bL[i])*w[k,i]+root2*sum(ub[i,j]*vb[j,k] for j in range(2))
                expectedb=-(bN[k]+bL[i])*wb[k,i]-root2*sum(u[i,j]*v[j,k] for j in range(2))
                errors.extend([abs(dw[k,i]-expected)/(1+abs(expected)),abs(dwb[k,i]-expectedb)/(1+abs(expectedb))])
        for P,Q in [(nodes[0],nodes[2]),(nodes[1],nodes[3])]:
            lhs=(chi(P)-chi(Q)+logshift(P-Q))*omega(P,Q)
            rhs=sum(psi(P,k)*psid(Q,k) for k in range(2))
            kernel_errors.append(abs(lhs-rhs)/(1+abs(rhs)))
        X=mp.matrix(6)
        for i in range(2):
            X[i,i]=-2*bL[i]
            for j in range(2):X[i,j+2]=root2*u[i,j];X[j+2,i]=root2*ub[i,j]
        for j in range(2):
            X[j+2,j+2]=-2*bM[j]
            for k in range(2):X[j+2,k+4]=v[j,k];X[k+4,j+2]=-vb[j,k]
        for k in range(2):
            X[k+4,k+4]=bN[k]
            for i in range(2):X[i,k+4]=wb[k,i];X[k+4,i]=-w[k,i]
        Dmat=mp.diag([1,1,-1,-1,0,0])
        for flow in higher:
            gradient=(X+flow['nu0']*Dmat)**(flow['order']-1)
            field=Dmat*gradient-gradient*Dmat
            speed=flow['speed'];out=flow['outer'];central=flow['central']
            entries=[]
            for i in range(2):
                for j in range(2):
                    phase=ellL[i]-ellM[j];rho=out[i]-out[j+2]
                    entries.append((speed*(du[i,j]-phase*u[i,j])+rho*u[i,j],field[i,j+2]/root2))
                    entries.append((speed*(dub[i,j]+phase*ub[i,j])-rho*ub[i,j],field[j+2,i]/root2))
            for j in range(2):
                for k in range(2):
                    phase=ellM[j]-ellN[k];rho=out[j+2]-central[k]
                    entries.append((speed*(dv[j,k]-phase*v[j,k])+rho*v[j,k],field[j+2,k+4]))
                    entries.append((speed*(dvb[j,k]+phase*vb[j,k])-rho*vb[j,k],-field[k+4,j+2]))
            for k in range(2):
                for i in range(2):
                    phase=ellN[k]-ellL[i];rho=central[k]-out[i]
                    entries.append((speed*(dw[k,i]-phase*w[k,i])+rho*w[k,i],-field[k+4,i]))
                    entries.append((speed*(dwb[k,i]+phase*wb[k,i])-rho*wb[k,i],field[i,k+4]))
            flow['errors'].extend(abs(a-b)/(1+abs(b)) for a,b in entries)
        probes=[mp.mpc('.32','.07'),mp.mpc('.63','.42')]
        eigenvectors=[]
        for P in probes:
            nu=-(hL(P)+hM(P))/2;lam=(hM(P)-hL(P))/2
            r=mp.matrix([-kap[i]*omega(nodes[i],P) for i in range(2)]+
                        [omega(nodes[j+2],P) for j in range(2)]+
                        [psid(P,k) for k in range(2)])
            phases=ellL+ellM+ellN
            r=mp.matrix([r[a]*mp.exp(-1j*time*phases[a]) for a in range(6)])
            left=mp.matrix([[omega(P,nodes[i]) for i in range(2)]+
                            [-kap[j+2]*omega(P,nodes[j+2]) for j in range(2)]+
                            [-psi(P,k) for k in range(2)]])
            left=mp.matrix([[left[a]*mp.exp(1j*time*phases[a]) for a in range(6)]])
            pencil=X+mp.diag([nu,nu,-nu,-nu,0,0])
            residual=pencil*r-lam*r
            spectral_errors.append(max(abs(x) for x in residual)/(1+max(abs(x) for x in lam*r)))
            residual_left=left*pencil-lam*left
            spectral_errors.append(max(abs(x) for x in residual_left)/(1+max(abs(x) for x in lam*left)))
            eigenvectors.append((nu,-left/initial_ordinary(P),r/initial_dual(P)))
        for a,b in [(0,1),(1,0)]:
            nuP,left,_=eigenvectors[a];nuQ,_,right=eigenvectors[b]
            reconstructed=-(left*right)[0]/(nuP-nuQ)
            expected=omega(probes[a],probes[b])/(initial_ordinary(probes[a])*initial_dual(probes[b]))
            bilinear_errors.append(abs(reconstructed-expected)/(1+abs(expected)))
    assert max(errors)<mp.mpf('1e-38'),max(errors)
    assert max(kernel_errors)<mp.mpf('1e-38'),max(kernel_errors)
    assert max(spectral_errors)<mp.mpf('1e-38'),max(spectral_errors)
    assert max(bilinear_errors)<mp.mpf('1e-38'),max(bilinear_errors)
    for flow in higher:
        assert max(flow['errors'])<mp.mpf('1e-38'),(flow['order'],max(flow['errors']))
    return {'backend':'Generapy rtheta_jet, order 3; cached complete jets',
            'generapy_version':importlib.metadata.version('generapy'),
            'dps':mp.mp.dps,'genus':1,'L_M_N':[2,2,2],
            'curve_data':{'tau':str(tau),'punctures':[str(x) for x in R],
                          'ordinary_divisor':[str(x) for x in D],
                          'dual_divisor':[str(x) for x in roots if x not in D],
                          'reference':str(z0),'points':[str(x) for x in nodes],
                          'principal_coefficients':[str(x) for x in principal]},
            'all_24_rates_max_scaled_error':mp.nstr(max(errors),12),
            'rank_two_kernel_rate_max_scaled_error':mp.nstr(max(kernel_errors),12),
            'prescribed_curve_eigenvector_max_scaled_error':mp.nstr(max(spectral_errors),12),
            'normalized_eigenvector_kernel_max_scaled_error':mp.nstr(max(bilinear_errors),12),
            'reference_partition_of_unity_error':mp.nstr(partition_error,12),
            'times':['0','.031','.087'],'jet_cache_entries':jet.cache_info().currsize,
            'higher_flows':[{'trace_order':flow['order'],'spectral_shift':str(flow['nu0']),
                             'max_all_24_rate_error':mp.nstr(max(flow['errors']),12),
                             'BA_time_scale':str(flow['speed']),
                             'third_group_finite_parts':[str(x) for x in flow['central']]}
                            for flow in higher],
            'scope':'Positive-genus multipuncture realization and spectral map. '
                    'A genus-one special spectral level, not generic genus-seven initial data.'}


if __name__=='__main__':
    result=check()
    (Path(__file__).parent/'results'/'multipoint_BA_elliptic_check.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
