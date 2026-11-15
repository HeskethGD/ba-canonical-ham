"""Check the eigenvector finite-part phase formula on exact genus-two data.

All curve data stay symbolic or rational; no sigma evaluations are needed.
The all-mode check uses the existing exact BA embedding, independently of
the meromorphic-function calculation performed here.
"""
from pathlib import Path
import json
import sympy as s

from quartic_hyperelliptic_direction import full_mode_test


def full_genus_one_field_check():
    """Independent symbolic check on all six canonical coordinates."""
    u1,u2,c,v1,v2,d = s.symbols('u1 u2 c v1 v2 d')
    b1,b2,rho = s.symbols('b1 b2 rho')
    e1,e2 = rho-b1,rho-b2
    X = s.Matrix([[e1,s.sqrt(2)*c,v1],
                  [s.sqrt(2)*d,e2,u2],[-u1,-v2,rho]])
    H4 = s.expand(s.trace(X**4)/4)
    H0 = b1*u1*v1+b2*u2*v2+s.sqrt(2)*(u1*u2*c+v1*v2*d)
    C1,C2 = u1*v1-c*d,u2*v2-c*d
    m0 = u1*v1-u2*v2
    correction_constant = m0-(e2-e1)*rho
    gamma1 = -3*e1**2+2*C1-(e2-e1)*e1-correction_constant
    gamma2 = -3*e2**2+2*C2+(e2-e1)*e2+correction_constant
    q,p = [u1,u2,c],[v1,v2,d]
    def field(f,H):
        return sum(s.diff(f,a)*s.diff(H,b)-s.diff(f,b)*s.diff(H,a)
                   for a,b in zip(q,p))
    for mode in q+p:
        assert s.expand(field(mode,H4)+(e1+e2+rho)*field(mode,H0)
                        -gamma1*field(mode,C1)-gamma2*field(mode,C2)) == 0
    print('Full genus-one field identity passes symbolically on all six coordinates.')
    return {'all_initial_modes_symbolic':True, 'six_field_components':True,
            'Abel_rate_relative_to_H0':str(s.expand(-(e1+e2+rho))),
            'gauge_rates':[str(s.expand(gamma1)),str(s.expand(gamma2))]}


def check_case(nodes, f, h, P, J):
    x, y, z, tau = s.symbols('x y local_x local_s')
    rho = s.Rational(3, 7)
    nodes = list(map(s.sympify, nodes))
    ordinates = [h.subs(x, a) for a in nodes]
    weights = [1/(2*(nodes[0]-nodes[1])), 1/(2*(nodes[1]-nodes[0])),
               1/(2*(nodes[2]-nodes[3])), 1/(2*(nodes[3]-nodes[2]))]
    phases = [-sum(weights[j]*(ordinates[j]-ordinates[i])/(nodes[j]-nodes[i])
                   for j in (range(2,4) if i < 2 else range(2)))
              for i in range(4)]
    RL = s.prod(x-a for a in nodes[:2])
    RM = s.prod(x-a for a in nodes[2:])
    lL = s.interpolate(list(zip(nodes[:2], ordinates[:2])), x)
    lM = s.interpolate(list(zip(nodes[2:], ordinates[2:])), x)
    hL, hM = (y+lL)/(2*RL), (y-lM)/(2*RM)
    nu, lam = -(hL+hM)/2, rho+(hM-hL)/2
    q = lam**3/nu
    A = s.interpolate(list(zip(nodes[:2], [-v for v in phases[:2]])), x)
    B = s.interpolate(list(zip(nodes[2:], [-v for v in phases[2:]])), x)
    # This function cancels q's principal parts at all four outer points.
    meromorphic = -hL**2/4+hM**2/4+(A+3*rho/2)*hL+(B+3*rho/2)*hM
    SL, SM = sum(nodes[:2]), sum(nodes[2:])
    h3, h2 = s.expand(h).coeff(x,3), s.expand(h).coeff(x,2)
    constant_at_infinity = (SM-SL)/2+h3*(lL.coeff(x,1)-lM.coeff(x,1))/4
    meromorphic -= constant_at_infinity

    # Principal part at the central point, in local parameter s=x^-1/2.
    lam4 = s.expand(f).coeff(x,4)
    aL, aM = lam4/8+SL, lam4/8+SM
    # Include the even s^2 terms, which determine the normalization constant.
    HL = 1/tau+aL*tau+lL.coeff(x,1)*tau**2/2
    HM = 1/tau+aM*tau-lM.coeff(x,1)*tau**2/2
    infinity = -HL**2/4+HM**2/4+(A.subs(x,tau**-2)+3*rho/2)*HL \
               +(B.subs(x,tau**-2)+3*rho/2)*HM-constant_at_infinity
    infinity = s.series(infinity,tau,0,1).removeO().expand()
    direction = [-h3, -h2-h3*sum(nodes)-3*rho]
    desired = -direction[0]*(tau**-3+lam4*tau**-1/8)-direction[1]/tau
    assert s.expand(infinity-desired) == 0

    predicted_phases = []
    for i, a in enumerate(nodes):
        # General BA points in the second group use the opposite sheet.
        Y = ordinates[i] if i < 2 else -ordinates[i]
        local_y = Y+s.diff(f,x).subs(x,a)*z/(2*Y) \
                  +(s.diff(f,x,2).subs(x,a)/(4*Y)
                     -s.diff(f,x).subs(x,a)**2/(8*Y**3))*z**2
        # Cancel the poles using the curve relation before taking the limit.
        combined = s.cancel(meromorphic+q)
        numerator, denominator = s.fraction(combined)
        numerator = s.rem(numerator, y*y-f, y)
        combined = s.cancel(numerator/denominator)
        local = combined.subs({x:a+z, y:local_y}, simultaneous=True)
        finite = s.cancel(local).subs(z,0)
        if finite.has(s.nan, s.zoo):
            finite = s.limit(local,z,0)
        assert finite.is_finite
        gamma = -finite if i < 2 else finite
        gamma = s.factor(gamma)
        opposite = range(2,4) if i < 2 else range(2)
        expected = (a+sum(nodes[j] for j in opposite)+sum(nodes)/2+lam4/4
                    +phases[i]*(2*h3*a+2*h2+h3*sum(nodes))
                    -h3*sum(weights[j]*ordinates[j] for j in opposite)
                    +6*rho*phases[i]-3*rho**2)
        assert s.cancel(gamma-expected) == 0, (i, gamma, expected)
        predicted_phases.append(str(gamma))

    # The previous all-mode certificate checks these same rates against
    # the Hamiltonian field at multiple divisors on this fixed curve.
    certificate = full_mode_test(nodes, P, J, h, f, rounds=3)
    return {'curve':str(f), 'central_shift':str(rho),
            'Abel_direction':list(map(str,direction)),
            'phase_rates_from_finite_parts':predicted_phases,
            'central_principal_part_and_constant':True,
            'outer_finite_parts_match_all_four_rates':True,
            'all_mode_certificate':certificate}


def main():
    x = s.Symbol('x')
    cases = []
    for nodes, r1, r2, j1, j0, ell in [
        ([-5,-2,3,7],s.Rational(1,3),s.Rational(5,2),s.Rational(2,3),-4,s.Rational(3,5)),
        ([-7,-3,2,9],s.Rational(-1,2),s.Rational(4,3),-2,s.Rational(7,5),s.Rational(-2,3)),
    ]:
        P = s.expand((x-r1)*(x-r2))
        J = j1*x+j0
        R = s.prod(x-a for a in nodes)
        offset = (1-ell**2*(P.coeff(x,1)+sum(nodes)))/(2*ell)
        h = s.expand(J+2*P*(ell*x+offset))
        f = s.expand(h*h-4*ell**2*P*R)
        assert s.degree(f,x)==5 and f.coeff(x,5)==4
        assert s.gcd(f,s.diff(f,x))==1
        cases.append(check_case(nodes, f, h, P, J))
    report = {'full_genus_one_field':full_genus_one_field_check(),
              'exact_genus_two_cases':cases,
              'total_fixed_curve_divisor_states':6,
              'equations_per_state':16}
    destination = Path(__file__).parent/'results/fourth_trace_phase_normalization.json'
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print('Finite-part phases and all-mode certificates passed for both fixed curves.')


if __name__ == '__main__':
    main()
