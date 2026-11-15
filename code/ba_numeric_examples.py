"""Generapy sigma BA validation on hyperelliptic and trigonal curves.

All sigma evaluations and derivatives call Generapy directly. Constants
are initialized once outside evaluate_modes; no alternate theta evaluator
or evaluation-cache layer is used. Sparse analytic samples are plotted as
open markers against dense independent Runge--Kutta curves.
"""
from pathlib import Path
from dataclasses import dataclass
import argparse
import importlib.metadata
import json
import math
from time import perf_counter
import numpy as np
from mpmath import mp
from generapy import Curve
from reference_sigma_ba import GenerapySigma,ReferenceSigmaBA

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results';FIGURES=ROOT/'figures'


def scalar(z):
    return [mp.nstr(z.real,40),mp.nstr(z.imag,40)]


def restore(z,ctx):
    return ctx.mpc(z[0],z[1])


def pack_matrix(A):
    return [[scalar(A[i,j]) for j in range(A.cols)] for i in range(A.rows)]


def unpack_matrix(A,ctx):
    return ctx.matrix([[restore(z,ctx) for z in row] for row in A])


def trigonal_local_points(ctx,kvalues,sheets):
    mu=[ctx.mpf(a) for a in (2,3,5,7)];points=[]
    for kval,sheet in zip(kvalues,sheets):
        ss=1/ctx.mpf(str(kval))
        def perturb(t):return sum(mu[j]*t**(3*(j+1)) for j in range(4))
        def w(t):return ctx.exp(ctx.log1p(perturb(t))/3)
        A=[ctx.quad(lambda t:-t**4/w(t)**2,[0,ss]),
           ctx.quad(lambda t:-t/w(t)**2,[0,ss]),
           ctx.quad(lambda t:-1/w(t),[0,ss])]
        def regular(t):
            if not t:return ctx.zero
            return -ctx.expm1(-ctx.mpf(2)/3*ctx.log1p(perturb(t)))/t**2
        rate=1/ss+ctx.quad(regular,[0,ss])
        root=ctx.exp(2*ctx.pi*ctx.j*sheet/3)
        A=[root*A[0],root*A[1],A[2]/root];rate*=root
        x=ss**-3;y=root*ss**-4*w(ss)
        points.append(dict(x=scalar(x),y=scalar(y),abel=[scalar(a) for a in A],
                           rate=scalar(rate),k=str(kval),sheet=sheet))
    return points


def trigonal_setup(ctx):
    mu=[ctx.mpf(a) for a in (2,3,5,7)]
    def fy(x,y):return 3*y*y
    first=(lambda x,y:1/fy(x,y),lambda x,y:x/fy(x,y),lambda x,y:y/fy(x,y))
    second=(lambda x,y:(5*x*x+3*mu[0]*x+mu[1])*y/fy(x,y),
            lambda x,y:2*x*y/fy(x,y),lambda x,y:x*x/fy(x,y))
    curve=Curve({(0,3):1,(4,0):-1,(3,0):-2,(2,0):-3,(1,0):-5,(0,0):-7},
                ctx=ctx,differentials_kind_1=first,differentials_kind_2=second)
    print('Trigonal: computing periods and infinity characteristic.',flush=True)
    one=curve.periods_kind_1();two=curve.periods_kind_2()
    chart=curve.monomial_chart(-3,-4)
    seed=min(curve.chart_fibre(chart,0),key=lambda z:abs(z-1))
    infinity=curve.chart_place(chart,seed,ctx.mpf('.12'))
    characteristic=curve.riemann_constant(base_place=infinity).characteristic
    tau=(one.tau+one.tau.T)/2;kappa=(two.kappa+two.kappa.T)/2
    points=trigonal_local_points(ctx,
        ('0.6','0.9','1.4','0.7','1.1','1.8'),(0,1,2,0,1,2))
    # Independent public Abel-map comparison, including the infinity tail.
    native=curve.abel_map_kind_1((restore(points[0]['x'],ctx),restore(points[0]['y'],ctx)),
                                base_place=infinity).value
    native_error=max(abs(native[i]-restore(points[0]['abel'][i],ctx)) for i in range(3))
    print('Trigonal: Abel-map comparison',ctx.nstr(native_error,6),flush=True)
    return dict(kind='trigonal',genus=3,L=3,M=3,
                curve='y^3=x^4+2*x^3+3*x^2+5*x+7',
                omega=pack_matrix(one.omega),tau=pack_matrix(tau),kappa=pack_matrix(kappa),
                characteristic=[[str(a) for a in row] for row in characteristic],
                points=points,period_symmetry=str(one.symmetry_residual),
                kappa_symmetry=str(two.kappa_symmetry_residual),
                native_abel_error=str(native_error),
                second_kind_rate='FP integral of x^2 dx/(3y^2), constant term zero in s=x^(-1/3)')


def hyperelliptic_setup(ctx):
    # 4*x*(x^2-1)*(x^2-4).
    spec={(0,2):1,(5,0):-4,(3,0):20,(1,0):-16}
    curve=Curve(spec,ctx=ctx)
    print('Hyperelliptic: computing periods and Abel integrals.',flush=True)
    one=curve.periods_kind_1();two=curve.periods_kind_2()
    characteristic=curve.riemann_constant().characteristic
    points=[]
    for j,aa in enumerate(('0.1','2.3','0.5','3.0')):
        x=ctx.mpf(aa);y=2*ctx.sqrt(x*(x*x-1)*(x*x-4))
        if j>=2:y=-y
        A=curve.abel_map_kind_1((x,y)).value
        R=curve.abel_map_kind_2((x,y)).value
        points.append(dict(x=scalar(x),y=scalar(y),abel=[scalar(a) for a in A],rate=scalar(R[1])))
    return dict(kind='hyperelliptic',genus=2,L=2,M=2,
                curve='y^2=4*x*(x^2-1)*(x^2-4)',
                omega=pack_matrix(one.omega),tau=pack_matrix((one.tau+one.tau.T)/2),
                kappa=pack_matrix((two.kappa+two.kappa.T)/2),
                characteristic=[[str(a) for a in row] for row in characteristic],points=points,
                period_symmetry=str(one.symmetry_residual),
                kappa_symmetry=str(two.kappa_symmetry_residual),
                second_kind_rate='Generapy automatic odd-degree second-kind Abel integral')


def setup(kind,dps):
    RESULTS.mkdir(exist_ok=True)
    path=RESULTS/f'{kind}_sigma_ba_geometry_{dps}.json'
    ctx=mp.clone();ctx.dps=dps
    # Geometry is computed once for this run; JSON is a reproducibility
    # record, not a separate cache used by the evaluator.
    data=(trigonal_setup if kind=='trigonal' else hyperelliptic_setup)(ctx)
    data['dps']=dps;data['generapy_version']=importlib.metadata.version('generapy')
    path.write_text(json.dumps(data,indent=2)+'\n')
    return ctx,data


def error(a,b):
    a=np.asarray(a);b=np.asarray(b)
    return float(np.max(abs(a-b)/(1+np.maximum(abs(a),abs(b)))))


@dataclass
class Network:
    data:dict
    sigma:object
    modes:list
    z0:tuple
    x:np.ndarray
    y:np.ndarray
    weights:np.ndarray
    gauges:np.ndarray
    phases:np.ndarray
    beta:np.ndarray

    @property
    def L(self):return self.data['L']
    @property
    def M(self):return self.data['M']

    def evaluate(self,t):
        L,M=self.L,self.M;K=L+M
        abel_time=-1j*t
        z=list(self.z0);z[-1]+=self.sigma.ctx.mpc(abel_time);d=self.sigma.jet(z,3)
        f,g=[np.array(a,dtype=complex) for a in zip(*(a.jets(abel_time,d) for a in self.modes))]
        p=complex(3*(d[2]/d[0]-(d[1]/d[0])**2))
        pd=complex(3*(d[3]/d[0]-3*d[2]*d[1]/d[0]**2+2*(d[1]/d[0])**3))
        def kernel(i,j):
            delta=self.x[i]-self.x[j]
            if self.data['kind']=='trigonal':
                value=(g[j,0]*f[i,2]-g[j,1]*f[i,1]+g[j,2]*f[i,0]+p*g[j,0]*f[i,0])/delta
                derivative=(g[j,0]*f[i,3]+g[j,3]*f[i,0]+pd*g[j,0]*f[i,0]
                            +p*(g[j,1]*f[i,0]+g[j,0]*f[i,1]))/delta
            else:
                value=(g[j,0]*f[i,1]-g[j,1]*f[i,0])/delta
                derivative=(g[j,0]*f[i,2]-g[j,2]*f[i,0])/delta
            return value,derivative
        u=np.concatenate((g[:L,0],f[L:,0]))*self.gauges
        ubar=np.concatenate((f[:L,0],g[L:,0]))*self.weights/self.gauges
        du=np.concatenate((g[:L,1],f[L:,1]))*self.gauges
        dubar=np.concatenate((f[:L,1],g[L:,1]))*self.weights/self.gauges
        C=np.empty((L,M),complex);Cbar=C.copy();dC=C.copy();dCbar=C.copy()
        kernel_error=0.
        for i in range(L):
            for jj in range(M):
                j=L+jj;O,Od=kernel(i,j);R,Rd=kernel(j,i)
                factor=math.sqrt(2)*self.weights[i]*self.weights[j]/(self.gauges[i]*self.gauges[j])
                reverse=-math.sqrt(2)*self.gauges[i]*self.gauges[j]
                C[i,jj]=factor*O;dC[i,jj]=factor*Od
                Cbar[i,jj]=reverse*R;dCbar[i,jj]=reverse*Rd
                kernel_error=max(kernel_error,error(Od,f[i,0]*g[j,0]),error(Rd,f[j,0]*g[i,0]))
        if self.data['kind']=='trigonal':
            beta=g[:,0]*f[:,2]-g[:,1]*f[:,1]+g[:,2]*f[:,0]+p*g[:,0]*f[:,0]
        else:beta=g[:,0]*f[:,1]-g[:,1]*f[:,0]
        values=np.concatenate((u,ubar,C.ravel(),Cbar.ravel()))
        rates=np.concatenate((du,dubar,dC.ravel(),dCbar.ravel()))
        return values,-1j*rates,dict(kernel_derivative=kernel_error,beta_drift=error(beta,self.beta),
                                 sigma_denominator=float(abs(d[0])))

    def geometric_rhs(self,z):
        L,M=self.L,self.M;K=L+M;P=L*M
        u,ubar=z[:K],z[K:2*K]
        C=z[2*K:2*K+P].reshape(L,M);Cbar=z[2*K+P:].reshape(L,M)
        du=self.phases*u;dubar=-self.phases*ubar
        du[:L]+=math.sqrt(2)*Cbar@ubar[L:];du[L:]+=math.sqrt(2)*Cbar.T@ubar[:L]
        dubar[:L]-=math.sqrt(2)*C@u[L:];dubar[L:]-=math.sqrt(2)*C.T@u[:L]
        return np.concatenate((du,dubar,(math.sqrt(2)*ubar[:L,None]*ubar[None,L:]).ravel(),
                               (-math.sqrt(2)*u[:L,None]*u[None,L:]).ravel()))

    def rhs(self,z):
        return -1j*self.geometric_rhs(z)

    def integrals(self,z):
        L,M=self.L,self.M;K=L+M;P=L*M
        u,ubar=z[:K],z[K:2*K]
        C=z[2*K:2*K+P].reshape(L,M);Cbar=z[2*K+P:].reshape(L,M)
        I=u*ubar;J=C*Cbar
        H=self.phases@I+math.sqrt(2)*(u[:L]@C@u[L:]+ubar[:L]@Cbar@ubar[L:])
        return np.concatenate(([H],I[:L]-J.sum(axis=1),I[L:]-J.sum(axis=0)))


def network(data,ctx,z0=None):
    L=data['L'];M=data['M'];K=L+M;G=data['genus']
    omega,tau,kappa=[unpack_matrix(data[n],ctx) for n in ('omega','tau','kappa')]
    original=[[ctx.mpf(a) for a in row] for row in data['characteristic']]
    characteristic=[[ctx.nint(2*a)/2 for a in row] for row in original]
    char_error=max(abs(a-b) for row,rounded in zip(original,characteristic) for a,b in zip(row,rounded))
    if char_error>ctx.mpf('1e-7'):raise ValueError('Infinity characteristic is not half-integral')
    if z0 is None:
        if data['kind']=='hyperelliptic':
            z0=tuple(omega[j,1]+ctx.mpc(0,'.8') for j in range(G))
        else:
            z0=tuple(ctx.mpc(a,b) for a,b in [('0.6','0.3'),('0.2','0.15'),('0.1','0.2')])
    else:z0=tuple(ctx.mpc(a) for a in z0)
    sigma=GenerapySigma(ctx,omega,tau,kappa,characteristic)
    points=data['points'];x=np.array([complex(restore(p['x'],ctx)) for p in points])
    y=np.array([complex(restore(p['y'],ctx)) for p in points])
    reference=sigma.jet(z0,0)[0]
    modes=[ReferenceSigmaBA(sigma,[restore(a,ctx) for a in p['abel']],restore(p['rate'],ctx),
                            z0,reference_denominator=reference) for p in points]
    d=sigma.jet(z0,3);p=complex(3*(d[2]/d[0]-(d[1]/d[0])**2))
    f,g=[np.array(a,dtype=complex) for a in zip(*(a.jets(0,d) for a in modes))]
    if data['kind']=='trigonal':
        beta=f[:,2]-g[:,1]*f[:,1]+g[:,2]+p
        residues=[]
        for idx in (slice(0,L),slice(L,K)):
            E=np.array([np.ones(len(x[idx])),x[idx],y[idx]])
            residues.extend(np.linalg.solve(E,np.array([0,0,1])))
        residues=np.array(residues);weights=residues*3*y*y/(2*beta)
        def h(group,a,b):return sum(residues[j]*(b*b+y[j]*b+y[j]*y[j])/(a-x[j]) for j in group)
    else:
        beta=f[:,1]-g[:,1]
        raw=np.array([1/(2*(x[0]-x[1])),1/(2*(x[1]-x[0])),
                      1/(2*(x[2]-x[3])),1/(2*(x[3]-x[2]))])
        weights=raw*y/beta
        def h(group,a,b):return sum(raw[j]*(b+y[j])/(a-x[j]) for j in group)
    phases=np.array([-h(range(L,K),x[i],y[i]) for i in range(L)]
                    +[h(range(L),x[j],y[j]) for j in range(L,K)])
    gauges=np.sqrt(weights)
    obj=Network(data,sigma,modes,z0,x,y,weights,gauges,phases,beta)
    return obj,float(char_error)


def rk4(fun,z,h):
    k1=fun(z);k2=fun(z+h*k1/2);k3=fun(z+h*k2/2);k4=fun(z+h*k3)
    return z+h*(k1+2*k2+2*k3+k4)/6


def prepare_example(kind,dps=24):
    ctx,data=setup(kind,dps)
    constants,char_error=network(data,ctx)
    constants.data['characteristic_rounding']=char_error
    return constants


def evaluate_modes(t,constants):
    """Physical time t: one call for modes, time rates, powers and checks."""
    values,rates,checks=constants.evaluate(t)
    K=constants.L+constants.M;P=constants.L*constants.M
    return dict(values=values,rates=rates,elementary_powers=values[:K]*values[K:2*K],
                composite_powers=values[2*K:2*K+P]*values[2*K+P:],checks=checks)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=['hyperelliptic','trigonal','both'],default='both')
    parser.add_argument('--dps',type=int,default=24)
    parser.add_argument('--duration',type=float,default=6.)
    parser.add_argument('--samples',type=int,default=41,help='Sparse direct Generapy BA samples')
    parser.add_argument('--dense-samples',type=int,default=601)
    parser.add_argument('--steps',type=int,default=6000)
    parser.add_argument('--check-only',action='store_true')
    args=parser.parse_args()
    for kind in (['hyperelliptic','trigonal'] if args.case=='both' else [args.case]):
        started=perf_counter();obj=prepare_example(kind,args.dps)
        setup_seconds=perf_counter()-started
        timings=[]
        for t in (0,.1,.1):
            started=perf_counter();sample=evaluate_modes(t,obj)
            timings.append(perf_counter()-started)
        obj.data['setup_seconds']=setup_seconds
        obj.data['evaluation_benchmark_seconds']=timings
        print(kind,'setup seconds',setup_seconds,'evaluation seconds',timings,flush=True)
        print(kind,'canonical residual',error(sample['rates'],obj.rhs(sample['values'])),
              sample['checks'],flush=True)
        if args.check_only:continue
        run(obj,args)


def run(obj,args):
    L,M=obj.L,obj.M;K=L+M
    mesh=args.samples-1;dense_mesh=args.dense_samples-1
    if min(mesh,dense_mesh)<2:raise ValueError('At least three samples are required')
    multiple=math.lcm(mesh,dense_mesh)
    steps=math.ceil(args.steps/multiple)*multiple
    grid=np.linspace(0,args.duration,args.samples)
    dense_grid=np.linspace(0,args.duration,args.dense_samples)
    print(obj.data['kind'],'evaluating direct Generapy BA modes on',args.samples,'points.',flush=True)
    analytic=[];residual=0.;kernel=0.;beta=0.;minimum=np.inf;sample_seconds=[]
    for plot_time in grid:
        started=perf_counter();data=evaluate_modes(plot_time,obj)
        sample_seconds.append(perf_counter()-started)
        values,rates,checks=data['values'],data['rates'],data['checks']
        analytic.append(values)
        residual=max(residual,error(rates,obj.rhs(values)))
        kernel=max(kernel,checks['kernel_derivative']);beta=max(beta,checks['beta_drift'])
        minimum=min(minimum,checks['sigma_denominator'])
    analytic=np.array(analytic)
    print('All-time canonical residual:',residual,'kernel:',kernel,'beta:',beta,flush=True)
    if not np.isfinite(analytic).all():raise ValueError('Nonfinite BA curve')
    def solve(n):
        h=args.duration/n;z=analytic[0].copy();sparse=[z.copy()];dense=[z.copy()]
        for step in range(1,n+1):
            z=rk4(obj.rhs,z,h)
            if not np.isfinite(z).all():raise ValueError('ODE left the finite chart')
            if step%(n//mesh)==0:sparse.append(z.copy())
            if step%(n//dense_mesh)==0:dense.append(z.copy())
        return np.array(sparse),np.array(dense)
    coarse,_=solve(steps);fine,dense=solve(2*steps)
    coarse_error=error(coarse,analytic);fine_error=error(fine,analytic)
    P=L*M
    def powers(z):
        return np.concatenate((z[:,:K]*z[:,K:2*K],
                               z[:,2*K:2*K+P]*z[:,2*K+P:]),axis=1)
    power_error=error(powers(fine),powers(analytic))
    initial=obj.integrals(analytic[0])
    def reality(z):
        z=np.asarray(z)
        return max(error(z[...,K:2*K],z[...,:K].conj()),
                   error(z[...,2*K+P:],z[...,2*K:2*K+P].conj()))
    BA_reality=reality(analytic);ODE_reality=reality(dense)
    physical=obj.data['kind']=='hyperelliptic'
    if physical:
        assert max(BA_reality,ODE_reality)<1e-8
        assert max(abs(obj.phases.imag))<1e-8 and min(obj.weights.real)>0
    BA_drift=max(error(obj.integrals(z),initial) for z in analytic)
    ODE_drift=max(error(obj.integrals(z),initial) for z in dense)
    print('Independent ODE errors:',coarse_error,fine_error,
          'conserved drift:',BA_drift,ODE_drift,flush=True)
    report=dict(case=obj.data['kind'],genus=obj.data['genus'],group_sizes=[L,M],
                canonical_pairs=K+L*M,curve=obj.data['curve'],dps=args.dps,
                generapy_version=obj.data['generapy_version'],duration=args.duration,
                analytic_samples=args.samples,dense_ode_samples=args.dense_samples,
                coarse_steps=steps,fine_steps=2*steps,
                time=('z(t)=z0-i*t*e_g; real t; brackets -i; '
                      + ('real Hamiltonian H0 on the all-positive real form'
                         if physical else 'complex Hamiltonian H0; independent partners')),
                z0=[scalar(z) for z in obj.z0],
                normalized_residue_weights=[[float(z.real),float(z.imag)] for z in obj.weights],
                phase_coefficients=[[float(z.real),float(z.imag)] for z in obj.phases],
                all_mode_canonical_residual=residual,kernel_derivative_residual=kernel,
                diagonal_pairing_drift=beta,coarse_ode_error=coarse_error,
                fine_ode_error=fine_error,step_halving_ratio=coarse_error/max(fine_error,1e-300),
                fine_modal_power_error=power_error,
                BA_conserved_quantity_drift=BA_drift,ODE_conserved_quantity_drift=ODE_drift,
                infinity_characteristic_rounding=obj.data['characteristic_rounding'],
                minimum_sigma_denominator=minimum,
                sigma_backend='Generapy kleinian_sigma_jet for every analytic sample',
                setup_seconds=obj.data['setup_seconds'],
                evaluation_benchmark_seconds=obj.data['evaluation_benchmark_seconds'],
                mean_sample_seconds=float(np.mean(sample_seconds)),
                plotted_modes=['Re xi_1',f'Re xi_{L+1}',f'Re xi_1,{L+1}','Re (xi_1 xi_1_bar)'],
                plot_style='black and white; sparse open BA markers, dense solid Runge--Kutta lines',
                all_modes_checked=True,all_powers_checked=True,
                physical_real_form=physical,
                BA_reality_error=BA_reality,ODE_reality_error=ODE_reality,
                minimum_elementary_power=float(np.min(powers(analytic)[:,:K].real)),
                maximum_imaginary_power=float(np.max(abs(powers(analytic).imag))))
    assert max(residual,kernel,beta,BA_drift)<1e-8,report
    assert max(fine_error,power_error,ODE_drift)<1e-6,report
    assert fine_error<coarse_error/5 or coarse_error<1e-10,report
    RESULTS.mkdir(exist_ok=True);FIGURES.mkdir(exist_ok=True)
    stem=obj.data['kind']+'_sigma_ba'
    (RESULTS/(stem+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(RESULTS/(stem+'.npz'),analytic_t=grid,ode_t=dense_grid,
                        BA=analytic,ODE=dense,weights=obj.weights,phases=obj.phases,
                        z0=np.asarray(obj.z0,dtype=complex))
    figure(obj,grid,analytic,dense_grid,dense,FIGURES/stem)
    print('PASS:',stem,flush=True)
    return report


def figure(obj,grid,analytic,dense_grid,direct,path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from cycler import cycler
    L,M=obj.L,obj.M;K=L+M
    selections=[(analytic[:,0].real,direct[:,0].real,r'$\operatorname{Re}\xi_1$'),
                (analytic[:,L].real,direct[:,L].real,rf'$\operatorname{{Re}}\xi_{{{L+1}}}$'),
                (analytic[:,2*K].real,direct[:,2*K].real,rf'$\operatorname{{Re}}\xi_{{1,{L+1}}}$'),
                ((analytic[:,0]*analytic[:,K]).real,(direct[:,0]*direct[:,K]).real,
                 (r'$|\xi_1|^2$' if obj.data.get('kind',obj.data.get('case'))=='hyperelliptic'
                  else r'$\operatorname{Re}(\xi_1\bar\xi_1)$'))]
    with plt.rc_context({'axes.prop_cycle':cycler(color=['black']),
                         'font.family':'serif','font.size':10,'savefig.facecolor':'white',
                         'axes.labelcolor':'black','text.color':'black',
                         'xtick.color':'black','ytick.color':'black'}):
        fig,axes=plt.subplots(4,1,figsize=(6.4,8.0),sharex=True,layout='constrained')
        kind=obj.data.get('kind') or obj.data['case']
        fig.suptitle(f"Genus {obj.data['genus']} {kind}: "
                     f"{K+L*M} mode pairs",fontsize=11)
        for i,(ax,(a,b,label)) in enumerate(zip(axes,selections)):
            ax.plot(dense_grid,b,color='black',linestyle='-',linewidth=1.0,
                    label="Runge--Kutta")
            ax.plot(grid,a,color='black',linestyle='None',marker=['o','s','^','D'][i],
                    markersize=3.5,markerfacecolor='white',markeredgecolor='black',
                    markeredgewidth=.7,label='Generapy BA formula')
            ax.set_ylabel(label);ax.text(.01,.91,f'({chr(97+i)})',transform=ax.transAxes)
            ax.spines[['top','right']].set_visible(False)
        axes[0].legend(loc='best',frameon=False,fontsize=8)
        axes[-1].set_xlabel(r'$t\quad (z=z_0-it e_g)$')
        fig.savefig(path.with_suffix('.pdf'))
        fig.savefig(path.with_suffix('.png'),dpi=180)
        plt.close(fig)


if __name__=='__main__':main()
