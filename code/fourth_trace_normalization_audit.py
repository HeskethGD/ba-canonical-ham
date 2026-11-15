"""Audit checks: coincident phases on a singular genus-two BA sector.

Uses only modal powers, the inversion polynomial and the canonical
matrix field. No theta evaluator, numerical differentiation or ODE solve.
"""
from pathlib import Path
import json
import mpmath as mp
import sympy as s


def branch_finite_part_identity():
    """Exact local check, symbolic in every parameter that enters the limit."""
    t,C,rho,K = s.symbols('t C rho K', nonzero=True)
    observed = []
    for HL,HM in [(1/t,C*t),(C*t,1/t)]:
        nu,lam = -(HL+HM)/2,rho+(HM-HL)/2
        correction = (HM**2-HL**2)/4+3*rho*(HL+HM)/2-K
        observed.append(s.limit(s.cancel(correction+lam**3/nu),t,0))
    assert s.expand(observed[0]-(3*rho**2-C-K)) == 0
    assert s.expand(observed[1]-(-3*rho**2+C-K)) == 0
    return {'symbolic_branch_finite_parts':list(map(str,observed)),
            'C_means':'x(Q_k)-a_fifth',
            'K_means':'(sum_second_nodes-sum_first_nodes)/2'}


def branch_point_check():
    ctx = mp.mp.clone()
    ctx.dps = 55
    # Five finite branch points; four supply the two elementary groups.
    branches = list(map(ctx.mpf, [-6, -2, 3, 8, 13]))
    nodes = branches[:4]
    s2 = ctx.sqrt(2)
    weights = [1/(2*(nodes[0]-nodes[1])), 1/(2*(nodes[1]-nodes[0])),
               1/(2*(nodes[2]-nodes[3])), 1/(2*(nodes[3]-nodes[2]))]
    def f(x):
        return 4*ctx.fprod(x-a for a in branches)
    records = []
    for roots in [(ctx.mpf(17),ctx.mpf(21)),(ctx.mpf(18),ctx.mpf(25))]:
        r1,r2 = roots
        y1,y2 = ctx.sqrt(f(r1)),ctx.sqrt(f(r2))
        def P(x):
            return (x-r1)*(x-r2)
        def J(x):
            return y1+(y2-y1)*(x-r1)/(r2-r1)
        u = list(map(ctx.mpf,[2,3,5,7]))
        v = [weights[k]*P(nodes[k])/u[k] for k in range(4)]
        rates = [J(a)/(2*P(a)) for a in nodes]
        C = [[s2*v[i]*v[j]*(rates[j]-rates[i])/(nodes[j]-nodes[i])
              for j in range(2,4)] for i in range(2)]
        E = [[-s2*u[i]*u[j]*(rates[j]-rates[i])/(nodes[j]-nodes[i])
              for j in range(2,4)] for i in range(2)]
        # Geometric H0 field, all b_k=0.
        du = [s2*sum(v[j]*E[i][j-2] for j in range(2,4)) for i in range(2)] \
             +[s2*sum(v[i]*E[i][j-2] for i in range(2)) for j in range(2,4)]
        dv = [-s2*sum(u[j]*C[i][j-2] for j in range(2,4)) for i in range(2)] \
             +[-s2*sum(u[i]*C[i][j-2] for i in range(2)) for j in range(2,4)]
        # Independent elementary Schrödinger/BA rate check at a branch.
        assert max(abs(du[k]-rates[k]*u[k])/(1+abs(du[k])) for k in range(4)) < ctx.mpf('1e-48')
        for rho in [ctx.mpf(0),ctx.mpf(3)/7]:
            X = ctx.eye(5)*rho
            dX0 = ctx.zeros(5)
            prediction = ctx.zeros(5)
            beta = [nodes[k]+sum(nodes[2:] if k < 2 else nodes[:2])
                    +sum(nodes)/2-sum(branches)-3*rho**2 for k in range(4)]
            for k in range(4):
                if k < 2:
                    X[k,4],X[4,k] = v[k],-u[k]
                    dX0[k,4],dX0[4,k] = dv[k],-du[k]
                    prediction[k,4] = -3*rho*dv[k]-beta[k]*v[k]
                    prediction[4,k] = 3*rho*du[k]-beta[k]*u[k]
                else:
                    X[k,4],X[4,k] = u[k],-v[k]
                    dX0[k,4],dX0[4,k] = du[k],-dv[k]
                    prediction[k,4] = -3*rho*du[k]+beta[k]*u[k]
                    prediction[4,k] = 3*rho*dv[k]+beta[k]*v[k]
            for i in range(2):
                for j in range(2,4):
                    X[i,j],X[j,i] = s2*C[i][j-2],s2*E[i][j-2]
                    dX0[i,j],dX0[j,i] = 2*v[i]*v[j],-2*u[i]*u[j]
                    prediction[i,j] = -3*rho*dX0[i,j]-(beta[i]+beta[j])*X[i,j]
                    prediction[j,i] = -3*rho*dX0[j,i]+(beta[i]+beta[j])*X[j,i]
            D = ctx.diag([1,1,-1,-1,0])
            actual = D*X**3-X**3*D
            error = max(abs(actual[i,j]-prediction[i,j])/(1+abs(actual[i,j]))
                        for i in range(5) for j in range(5))
            assert error < ctx.mpf('1e-47')
            records.append({'moving_divisor_x':list(map(str,roots)),
                            'central_shift':str(rho),
                            'max_scaled_matrix_field_error':str(error),
                            'all_16_mode_equations':True})
    return {'digits':ctx.dps,'curve':'4*prod(x-a), a in '+str(list(map(str,branches))),
            'all_phase_coefficients_zero':True,
            'two_equal_phase_blocks_of_size_2':True,
            'genus_two_BA_sector_of_size_five_pencil':True,
            'records':records}


if __name__ == '__main__':
    report = branch_point_check()
    report.update(branch_finite_part_identity())
    out = Path(__file__).parent/'results/fourth_trace_normalization_audit.json'
    out.write_text(json.dumps(report,indent=2)+'\n')
    print('Repeated-phase branch-point audit passed: four states, all 16 equations each.')
    print('Maximum scaled field error:',max(mp.mpf(r['max_scaled_matrix_field_error'])
                                           for r in report['records']))
