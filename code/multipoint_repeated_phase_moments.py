"""Check the conserved quadratic block moments with three repeated groups."""
from pathlib import Path
import json
import sympy as s


def check():
    n=2
    arr=lambda name:s.Matrix(n,n,lambda i,j:s.Symbol(f'{name}{i}{j}'))
    u,ub,v,vb,w,wb=[arr(name) for name in ('u','ub','v','vb','w','wb')]
    a,b,c=s.symbols('bL bM bN')
    eye=s.eye(n)
    X=s.BlockMatrix([[-2*a*eye,s.sqrt(2)*u,wb.T],
                     [s.sqrt(2)*ub.T,-2*b*eye,v],
                     [-w,-vb.T,c*eye]]).as_explicit()
    D=s.diag(1,1,-1,-1,0,0)
    K=s.diag(-3*a,-3*a,3*b,3*b,0,0)
    rate=-(D*X**2-X**2*D)+K*X-X*K
    slopes=[1,-1,0]
    indices=[list(range(2*r,2*r+2)) for r in range(3)]
    shift=s.Symbol('nu0')
    fields={'original_cubic':rate}
    for order in (3,4):
        gradient=(X+shift*D)**(order-1)
        fields[f'shifted_trace_{order}']=D*gradient-gradient*D
    checks={}
    for name,field in fields.items():
        checks[name]=[]
        for group in range(3):
            rows=indices[group]
            dQ=s.zeros(2)
            for other in range(3):
                if other==group:continue
                cols=indices[other]
                A=X.extract(rows,cols);B=X.extract(cols,rows)
                dA=field.extract(rows,cols);dB=field.extract(cols,rows)
                denominator=slopes[group]-slopes[other]
                dQ+=(dA*B+A*dB)/denominator
            assert dQ.applyfunc(s.expand)==s.zeros(2)
            checks[name].append(True)
    return {'L_M_N':[2,2,2],'all_phases_repeated_within_each_group':True,
            'all_three_quadratic_block_moments_conserved':checks,
            'scope':'Exact conservation check for the block moments used '
                    'to choose constant infinity-branch bases. General '
                    'regularity is established by the separate argument.'}


if __name__=='__main__':
    result=check()
    (Path(__file__).parent/'results'/'multipoint_repeated_phase_moments.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
