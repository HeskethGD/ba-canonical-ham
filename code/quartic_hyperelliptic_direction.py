"""Genus-two quartic-flow audit, keeping the prescribed BA data fixed.

The symbolic elementary-row reduction gives the Abel direction and phase
rates. Exact fixed-curve divisor tests check all sixteen equations, and an
independent Generapy calculation checks the actual sigma derivatives.
"""
from pathlib import Path
import json
import sympy as s

from support.bipartite_symbolic import build
from support.hyperelliptic_embedding import mode_embedding


def diagnostic():
    model = build(2, 2)
    u, v, z, b, c, H, pb = [model[k] for k in ('u','v','z','b','c','H','pb')]
    H = H.subs(dict.fromkeys(c, 0))
    F3 = s.Matrix([pb(q, H) for q in z])
    powers = s.Matrix([u[i]*v[i] for i in range(4)])
    power3 = powers.jacobian(z)*F3
    wronskians = s.Matrix([u[i]*pb(v[i], H)-v[i]*pb(u[i], H) for i in range(4)])
    x = s.Symbol('x')
    nodes = list(map(s.Integer, (6, 9, 8, 10)))
    ordinates = list(map(s.Integer, (1, 11, 1, 31)))
    weights = [-s.Rational(1,6), s.Rational(1,6), -s.Rational(1,4), s.Rational(1,4)]
    beta = [s.Rational(15,4), s.Integer(5), s.Rational(10,3), s.Rational(25,6)]
    p0,p1,j0,j1 = s.symbols('p0 p1 j0 j1')
    scales = s.symbols('scale1:5', nonzero=True)
    base = [p0,p1,j0,j1]+list(scales)
    ref = {p0:s.Rational(2813,64),p1:-s.Rational(213,16),
           j0:s.Rational(1023,128),j1:-s.Rational(35,32)}
    ref.update(dict.fromkeys(scales, 1))
    embedding = s.Matrix(mode_embedding(nodes, ordinates, p0,p1,j0,j1,scales))
    point = dict(zip(z, embedding.subs(ref)))
    point.update(dict(zip(b,beta)))
    point.update(dict.fromkeys(c, 0))
    X = model['X'].subs(point)
    DX = model['D']*X**3-X**3*model['D']
    field4 = []
    for variable in z:
        position = [(i,j,s.diff(model['X'][i,j],variable))
                    for i in range(5) for j in range(5)
                    if model['X'][i,j].has(variable)]
        i,j,factor = position[0]
        field4.append(s.simplify(DX[i,j]/factor))
    field4 = s.Matrix(field4)
    P, J = x*x+p1*x+p0, j1*x+j0
    f = 4*(x-7)**5-4*(x-7)+1
    rem = s.Poly(s.rem(f-J**2,P,x),x)
    constraints = s.Matrix([rem.nth(0),rem.nth(1)])
    E = embedding.jacobian(base).subs(ref)
    C = constraints.jacobian(base).subs(ref)
    tangent = E*s.Matrix.hstack(*C.nullspace())
    fixed_tangent = s.Matrix.hstack(tangent,field4).rank() == tangent.rank()
    print('Fixed-curve/marked-point tangent:', fixed_tangent, flush=True)
    w_rates = (wronskians.jacobian(z).subs(point)*field4).applyfunc(s.simplify)
    print('H4 derivative of the H0 Wronskians:', list(w_rates), flush=True)
    power4 = (powers.jacobian(z).subs(point)*field4).applyfunc(s.simplify)
    P4 = s.interpolate([(nodes[i],power4[i]/weights[i]) for i in range(2)],x)
    other = [s.simplify(P4.subs(x,nodes[i])-power4[i]/weights[i]) for i in range(4)]
    print('Same inversion polynomial at all marked nodes:', other, flush=True)
    P2 = J.subs(ref)
    P1 = (j0*x+p1*j0-p0*j1).subs(ref)
    A,B = s.symbols('A B')
    direction = s.solve(s.Poly(P4-A*P1-B*P2,x).all_coeffs(),[A,B])
    print('Candidate D4 = A partial_z1 + B partial_z2:', direction, flush=True)
    Js = s.interpolate([(nodes[i],power3[i]/weights[i]) for i in range(2)],x)
    J4 = (s.Matrix([s.expand(Js).coeff(x,k) for k in (0,1)]).jacobian(z)
          .subs(point)*field4).applyfunc(s.simplify)
    print('J4 coefficients:', list(J4), flush=True)
    if direction:
        rr = ref
        mode_phases=[]
        for i in range(4):
            Pi,Ji = P.subs(rr).subs(x,nodes[i]),J.subs(rr).subs(x,nodes[i])
            r2=(Ji-ordinates[i])/(2*Pi)
            r1=(nodes[i]+rr[p1])*r2-rr[j1]/2
            mode_phases.append(s.simplify(field4[i]/point[u[i]]-direction[A]*r1-direction[B]*r2))
        print('Candidate constant mode phase rates:', mode_phases, flush=True)
    return fixed_tangent


def symbolic_row_identity():
    """Reduce the quartic elementary row by degree-two interpolation moments.

    The only curve identity used before this polynomial calculation is
    P(a) r_minus(a) r_plus(a) = -W(a). This is f=J^2+4PW
    at a marked point. No characteristic-polynomial expansion is used.
    """
    x,a,own,opp,opp_product = s.symbols('x a own_sum opposite_sum opposite_product')
    h0,h1,h2,h3,p0,p1,j0,j1,lambda4,w0,w1,r = s.symbols(
        'h0 h1 h2 h3 p0 p1 j0 j1 lambda4 w0 w1 r_minus')
    P=x*x+p1*x+p0
    J=j1*x+j0
    h=h3*x**3+h2*x*x+h1*x+h0
    W=x**3+(lambda4/4-p1)*x*x+w1*x+w0
    moments={0:0,1:0,2:s.Rational(1,2),3:(a+opp)/2,
             4:(a*a+a*opp+opp*opp-opp_product)/2}
    def L(poly):
        poly=s.Poly(s.expand(poly),x)
        assert poly.degree()<=4
        return sum(coefficient*moments[degree[0]] for degree,coefficient in poly.terms())
    Yi=L(h)
    Yop=(h3*x+h2+h3*own)/2
    total_power=p1+(own+opp)/2
    observed=(-Yi*r-2*Yi**2+total_power+2*L(W)
              +r*L(J+h)+2*L(Yop*(J-h))-4*r*L(Yop*P))
    weighted_y=(h3*(opp*opp-opp_product)+h2*opp+h1)/2
    phase=(a+opp+(own+opp)/2+lambda4/4
           -Yi*(2*h3*a+2*h2+h3*(own+opp))-h3*weighted_y)
    A,B=-h3,-h2-h3*(own+opp)
    target=A*((a+p1)*r-j1/2)+B*r+phase
    assert s.expand(observed-target)==0
    print('General symbolic elementary-row identity passes.',flush=True)
    return {'all_curve_coefficients_and_node_symmetric_functions_symbolic':True,
            'Abel_direction':['-h3','-h2-h3*sum(a_k)'],
            'phase_formula':'a_k+s_opposite+sum(a)/2+lambda4/4'
                            '-Y_k*(2*h3*a_k+2*h2+h3*sum(a))'
                            '-h3*sum_opposite(kappa_j*y_j)'}


def quartic_field(model, point, central_shift=0):
    X=model['X'].subs(point)+central_shift*s.eye(5)
    DX=model['D']*X**3-X**3*model['D']
    rates=[]
    for variable in model['z']:
        position=[(i,j,s.diff(model['X'][i,j],variable))
                  for i in range(5) for j in range(5)
                  if model['X'][i,j].has(variable)]
        i,j,factor=position[0]
        rates.append(s.cancel(DX[i,j]/factor))
    return s.Matrix(rates)


def full_mode_test(nodes, P, J, h, f, rounds=3):
    """All sixteen equations on several divisors of one fixed curve."""
    x=s.Symbol('x')
    model=build(2,2)
    ys=[h.subs(x,a) for a in nodes]
    weights=[1/(2*(nodes[0]-nodes[1])),-1/(2*(nodes[0]-nodes[1])),
             1/(2*(nodes[2]-nodes[3])),-1/(2*(nodes[2]-nodes[3]))]
    sums=[sum(nodes[:2]),sum(nodes[2:])]
    Y=[s.cancel(sum(weights[j]*(ys[j]-ys[i])/(nodes[j]-nodes[i])
                    for j in (range(2,4) if i<2 else range(2)))) for i in range(4)]
    h3,h2=s.Poly(h,x).nth(3),s.Poly(h,x).nth(2)
    A,B=-h3,-h2-h3*sum(nodes)
    phase=[]
    for i in range(4):
        opp=range(2,4) if i<2 else range(2)
        phase.append(s.cancel(nodes[i]+sum(nodes[j] for j in opp)+sum(nodes)/2
                     +s.Poly(f,x).nth(4)/4-Y[i]*(2*h3*nodes[i]+2*h2+h3*sum(nodes))
                     -h3*sum(weights[j]*ys[j] for j in opp)))
    p0,p1,j0,j1=s.symbols('p0 p1 j0 j1')
    scales=s.symbols('c1:5',nonzero=True)
    base=[p0,p1,j0,j1]+list(scales)
    E=s.Matrix(mode_embedding(nodes,ys,p0,p1,j0,j1,scales))
    Eprime=E.jacobian(base)
    constraints=s.Matrix([s.Poly(s.rem(f-(j1*x+j0)**2,x*x+p1*x+p0,x),x).nth(k)
                          for k in (0,1)])
    for run in range(rounds):
        assert s.rem(f-J**2,P,x)==0
        assert all(P.subs(x,a)!=0 for a in nodes)
        ref={p0:P.subs(x,0),p1:P.coeff(x,1),j0:J.subs(x,0),j1:J.coeff(x,1)}
        ref.update(dict(zip(scales,map(s.Integer,(2,3,5,7)))))
        values=E.subs(ref)
        point=dict(zip(model['z'],values))
        point.update(dict(zip(model['b'],[2*y for y in Y])))
        point.update(dict.fromkeys(model['c'],0))
        central_shift=s.Rational(3,7) if run%2 else s.Integer(0)
        F4=quartic_field(model,point,central_shift)
        W=s.cancel((f-J**2)/(4*P)).expand()
        Q=s.Poly(f,x).nth(4)/4-2*ref[p1]
        K=s.expand(2*((x+Q)*P-W))
        assert s.degree(K,x)<=1
        k0,k1=K.subs(x,0),K.coeff(x,1)
        direction2=[ref[j0],ref[j1],k0,k1]
        direction1=[ref[p1]*ref[j0]-ref[p0]*ref[j1],ref[j0],
                    ref[p1]*k0-ref[p0]*k1,k0]
        for i,a in enumerate(nodes):
            r2=(J.subs(x,a)-ys[i])/(2*P.subs(x,a))
            direction2.append(ref[scales[i]]*r2)
            direction1.append(ref[scales[i]]*((a+ref[p1])*r2-ref[j1]/2))
        vector=A*s.Matrix(direction1)+(B-3*central_shift)*s.Matrix(direction2)
        for i in range(4):
            corrected_phase=phase[i]-6*central_shift*Y[i]-3*central_shift**2
            vector[4+i]+=corrected_phase*ref[scales[i]]
        predicted=Eprime.subs(ref)*vector
        assert all(s.cancel(q)==0 for q in F4-predicted)
        assert constraints.jacobian(base).subs(ref)*vector==s.zeros(2,1)
        # Change the moving divisor by adding one marked point on this same curve.
        a,y=nodes[run%4],ys[run%4]
        H=-J+s.cancel((y+J.subs(x,a))/P.subs(x,a))*P
        newP=s.cancel((f-H**2)/(4*P*(x-a))).expand()
        assert s.degree(newP,x)==2 and newP.coeff(x,2)==1
        J=s.rem(H,newP,x).expand()
        P=newP
    print(f'Fixed curve: all 16 mode equations pass at {rounds} divisor states.',flush=True)
    return {'curve':str(f),'marked_abscissae':list(map(str,nodes)),
            'direction':[str(A),str(B)],'phase_rates':list(map(str,phase)),
            'different_divisors_on_same_curve':rounds,'all_16_mode_equations':True,
            'fixed_curve_tangency':True,'nonzero_central_phase_checked':True}


def full_mode_certificates():
    x=s.Symbol('x')
    cases=[]
    for nodes,r1,r2,j1,j0,l in [
            ([-5,-2,3,7],s.Rational(1,3),s.Rational(5,2),s.Rational(2,3),-4,s.Rational(3,5)),
            ([-7,-3,2,9],s.Rational(-1,2),s.Rational(4,3),-2,s.Rational(7,5),s.Rational(-2,3))]:
        nodes=list(map(s.Integer,nodes))
        P=s.expand((x-r1)*(x-r2)); J=j1*x+j0
        R=s.prod(x-a for a in nodes)
        offset=(1-l*l*(P.coeff(x,1)+sum(nodes)))/(2*l)
        linear=l*x+offset
        h=s.expand(J+2*P*linear)
        f=s.expand(h*h-4*l*l*P*R)
        assert s.degree(f,x)==5 and f.coeff(x,5)==4
        assert s.gcd(f,s.diff(f,x))==1
        cases.append(full_mode_test(nodes,P,J,h,f))
    return cases


def generapy_direction_check():
    """Independently compare the elementary sigma derivatives with H4."""
    from mpmath import mp
    from generapy import Curve, kleinian_p
    from generapy.kleinian import kleinian_baker_akhiezer
    ctx=mp.clone();ctx.dps=32
    coefficients=list(map(ctx.mpf,['-3.04128','13.632','11.68','-18.4','-4','4']))
    curve=Curve({(0,2):1,**{(i,0):-c for i,c in enumerate(coefficients)}},ctx=ctx)
    first,second=curve.periods_kind_1(),curve.periods_kind_2()
    characteristic=curve.riemann_constant().characteristic
    z=[ctx.mpf('.13'),ctx.mpf('.27')]
    nodes=list(map(ctx.mpf,['3.1','3.7','4.3','4.9']))
    ordinates=[(-1 if i%2 else 1)*ctx.sqrt(ctx.polyval(coefficients[::-1],a))
               for i,a in enumerate(nodes)]
    weights=[1/(2*(nodes[0]-nodes[1])),-1/(2*(nodes[0]-nodes[1])),
             1/(2*(nodes[2]-nodes[3])),-1/(2*(nodes[2]-nodes[3]))]
    Y=[sum(weights[j]*(ordinates[j]-ordinates[i])/(nodes[j]-nodes[i])
           for j in (range(2,4) if i<2 else range(2))) for i in range(4)]
    vandermonde=ctx.matrix([[a**k for k in range(4)] for a in nodes])
    h=ctx.lu_solve(vandermonde,ctx.matrix(ordinates))
    central=ctx.mpf('.27')
    A,B=-h[3],-h[2]-h[3]*sum(nodes)-3*central
    phases=[]
    for i in range(4):
        opp=range(2,4) if i<2 else range(2)
        phases.append(nodes[i]+sum(nodes[j] for j in opp)+sum(nodes)/2
                      +coefficients[4]/4-Y[i]*(2*h[3]*nodes[i]+2*h[2]+h[3]*sum(nodes))
                      -h[3]*sum(weights[j]*ordinates[j] for j in opp)
                      -6*central*Y[i]-3*central**2)
    def kp(index):
        return kleinian_p(z,first.omega,first.tau,second.kappa,index,characteristic,ctx=ctx)
    p1,p0,j1,j0=-kp((1,1)),-kp((0,1)),-kp((1,1,1)),-kp((0,1,1))
    P=[a*a+p1*a+p0 for a in nodes]
    J=[j1*a+j0 for a in nodes]
    minus=[(J[i]-ordinates[i])/(2*P[i]) for i in range(4)]
    plus=[(J[i]+ordinates[i])/(2*P[i]) for i in range(4)]
    values=[];partners=[];predicted=[];predicted_star=[];direction_errors=[]
    for i in range(4):
        target=(nodes[i],ordinates[i])
        abel=curve.abel_map_kind_1(target).value
        integral=curve.abel_map_kind_2(target).value
        def psi(argument,sheet):
            return kleinian_baker_akhiezer(argument,sheet*abel,sheet*integral,
                       first.omega,first.tau,second.kappa,characteristic,ctx=ctx)
        q,qs=psi(z,1),weights[i]*psi(z,-1)
        d1=ctx.diff(lambda t:psi([t,z[1]],1),z[0])/q
        d2=ctx.diff(lambda t:psi([z[0],t],1),z[1])/q
        ds1=ctx.diff(lambda t:psi([t,z[1]],-1),z[0])/(qs/weights[i])
        ds2=ctx.diff(lambda t:psi([z[0],t],-1),z[1])/(qs/weights[i])
        direction_errors.extend([abs(d1-((nodes[i]+p1)*minus[i]-j1/2)),
                                 abs(ds1-((nodes[i]+p1)*plus[i]-j1/2))])
        values.append(q);partners.append(qs)
        predicted.append(A*d1+B*d2+phases[i])
        predicted_star.append(A*ds1+B*ds2-phases[i])
    X=ctx.matrix(5)
    for i in range(4):X[i,i]=Y[i]+central
    X[4,4]=central
    for i in range(2):X[i,4]=partners[i];X[4,i]=-values[i]
    for j in range(2,4):X[j,4]=values[j];X[4,j]=-partners[j]
    for i in range(2):
        for j in range(2,4):
            U=ctx.sqrt(2)*partners[i]*partners[j]*(plus[j]-plus[i])/(nodes[j]-nodes[i])
            V=-ctx.sqrt(2)*values[i]*values[j]*(minus[j]-minus[i])/(nodes[j]-nodes[i])
            X[i,j]=ctx.sqrt(2)*U;X[j,i]=ctx.sqrt(2)*V
    D=ctx.diag([1,1,-1,-1,0]);DX=D*X**3-X**3*D
    errors=[]
    for i in range(4):
        actual=(-DX[4,i] if i<2 else DX[i,4])/values[i]
        actual_star=(DX[i,4] if i<2 else -DX[4,i])/partners[i]
        errors.extend([abs(actual-predicted[i])/max(1,abs(actual),abs(predicted[i])),
                       abs(actual_star-predicted_star[i])/max(1,abs(actual_star),abs(predicted_star[i]))])
    assert max(errors)<ctx.mpf('1e-24')
    assert max(direction_errors)<ctx.mpf('1e-24')
    print('Generapy: all 8 elementary sigma derivatives agree with H4.',flush=True)
    return {'digits':ctx.dps,'all_eight_elementary_equations':True,
            'max_scaled_error':ctx.nstr(max(errors),10),
            'max_z1_identity_error':ctx.nstr(max(direction_errors),10),
            'central_phase':str(central)}


if __name__ == '__main__':
    report={'symbolic_row':symbolic_row_identity(),
            'fixed_data_certificate':diagnostic(),
            'all_mode_certificates':full_mode_certificates(),
            'generapy_check':generapy_direction_check()}
    destination=Path(__file__).with_name('results')/'quartic_hyperelliptic_direction.json'
    destination.write_text(json.dumps(report,indent=2)+'\n')
