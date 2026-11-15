"""Exact checks of the all-size phase-cluster insertion proof.

Append a mode to one equal-phase group. At outgoing new coordinates u=0
and fixed incoming v, the new spectral line and principal-moment
differences vanish identically in all old variables. Their u-gradients
must span the outside-group coordinates. This is the induction block.
"""
from pathlib import Path
import json
import numpy as np
import sympy as sp
from jrw_repeated_phase_general import PRIME, mm, add_row


def inverse(A):
    n=len(A)
    B=[[int(A[i,j]) % PRIME for j in range(n)]
       +[int(i==j) for j in range(n)] for i in range(n)]
    determinant=1
    for j in range(n):
        pivot=next((i for i in range(j,n) if B[i][j]),None)
        if pivot is None:
            return None,0
        if pivot!=j:
            B[j],B[pivot]=B[pivot],B[j]
            determinant=-determinant
        t=B[j][j];determinant=determinant*t % PRIME
        inv=pow(t,-1,PRIME)
        B[j]=[x*inv % PRIME for x in B[j]]
        for i in range(n):
            if i!=j and B[i][j]:
                c=B[i][j]
                B[i]=[(x-c*y) % PRIME for x,y in zip(B[i],B[j])]
    return np.array([row[n:] for row in B],dtype=np.int64),determinant % PRIME


def check(sizes,which,phases,seed=1949,profile="scalar",insert_phase=None):
    group=[g for g,n in enumerate(sizes) for _ in range(n)]
    d=len(group)
    slopes=[1,-1,0]
    intercepts=[]
    counts=[0,0,0]
    for g in group:
        if profile=="scalar":
            intercepts.append(phases[g])
        elif profile=="partial":
            intercepts.append(phases[g]+counts[g]%2)
        elif profile=="source":
            intercepts.append(0 if g<2 else phases[g]+counts[g])
        elif profile=="distinct":
            intercepts.append(phases[g]+counts[g])
        else:
            raise ValueError(profile)
        counts[g]+=1
    e=phases[which] if insert_phase is None else insert_phase
    inside=[i for i,g in enumerate(group) if g==which and intercepts[i]==e]
    outside=[i for i,g in enumerate(group) if g!=which]
    r,n=len(inside),len(outside)
    rng=np.random.default_rng(seed)
    X=np.zeros((d,d),dtype=np.int64)
    for i in range(d):
        X[i,i]=intercepts[i] % PRIME
        for j in range(d):
            if group[i]!=group[j]:
                X[i,j]=rng.integers(1,100)
    A=X[np.ix_(inside,outside)]
    B=X[np.ix_(outside,inside)]
    W=np.diag([pow((slopes[which]-slopes[group[i]]) % PRIME,-1,PRIME)
               for i in outside]).astype(np.int64)
    v=rng.integers(1,100,size=(n,1),dtype=np.int64)
    J=mm(mm(A,W),B);p=mm(mm(A,W),v)
    power=np.eye(r,dtype=np.int64)
    moment=[]
    for k in range(2,min(r,n)+2):
        # d[tr(Jnew^k)/k-tr(Jold^k)/k]/du at u=0.
        moment.append(mm(mm(mm(W,B),power),p).reshape(-1))
        power=mm(power,J)
    moment_basis=[]
    for row in moment:
        add_row(moment_basis,row)
    assert len(moment_basis)==min(r,n)
    spectral=[];sample_nu=[]
    if r<n:
        # Cofactor formula: gradient of det(line(nu)I-Xnew-nu D)
        # in the outgoing new row is -det(K)*K^{-1}_{OO} v.
        for nu in range(d+5):
            lam=e+slopes[which]*nu
            K=np.diag([(lam-slopes[g]*nu) % PRIME for g in group])-X
            inv,det=inverse(K % PRIME)
            if inv is None:
                continue
            row=(-det*mm(inv[np.ix_(outside,outside)],v).reshape(-1)) % PRIME
            assert not np.any(mm(A,row.reshape(-1,1)))
            spectral.append(row);sample_nu.append(nu)
    spectral_basis=[]
    for row in spectral:
        add_row(spectral_basis,row)
    assert len(spectral_basis)==max(n-r,0)
    combined=[]
    for row in spectral+moment:
        add_row(combined,row)
    assert len(combined)==n
    return dict(old_sizes=list(sizes),insert_group=which,
                phases=intercepts,profile=profile,insert_phase=e,
                old_original_group_size=sizes[which],
                old_equal_phase_cluster_size=r,outside_size=n,
                spectral_new_rank=len(spectral_basis),
                moment_new_rank=len(moment_basis),combined_new_rank=n,
                all_spectral_gradients_in_kernel_of_A=True,
                prime=PRIME,seed=seed,spectral_arguments=sample_nu)


def symbolic_moment_check():
    r,n=2,3
    A=sp.Matrix(r,n,sp.symbols('a0:6'))
    B=sp.Matrix(n,r,sp.symbols('b0:6'))
    W=sp.diag(*sp.symbols('w0:3'))
    u=sp.Matrix(1,n,sp.symbols('u0:3'))
    v=sp.Matrix(n,1,sp.symbols('v0:3'))
    J=A*W*B;p=A*W*v
    Jnew=J.row_join(p).col_join((u*W*B).row_join(u*W*v))
    zero={z:0 for z in u}
    for k in range(2,r+2):
        delta=sp.trace(Jnew**k)/k-sp.trace(J**k)/k
        expected=W*B*J**(k-2)*p
        for i in range(n):
            assert sp.expand(sp.diff(delta,u[i]).subs(zero)-expected[i])==0
        assert sp.expand(delta.subs(zero))==0
    return True


def eliminated_phase_schur_check():
    """Polynomial-in-nu check including a same-slope different-phase vertex."""
    nu=sp.Symbol("nu")
    # Growing group: one e=0 vertex, one e=2 vertex. Outside slopes -1,-1,0.
    A=sp.Matrix([[1,2,3]]);B=sp.Matrix([4,5,6])
    AC=sp.Matrix([[2,3,1]]);BC=sp.Matrix([1,4,2])
    Xout=sp.Matrix([[0,0,7],[0,1,8],[9,10,0]])
    X=sp.zeros(5)
    X[1,1]=2
    X[0,2:]=A;X[1,2:]=AC
    X[2:,0]=B;X[2:,1]=BC;X[2:,2:]=Xout
    D=sp.diag(1,1,-1,-1,0)
    K=nu*sp.eye(5)-X-nu*D
    W=sp.diag(sp.Rational(1,2),sp.Rational(1,2),1)
    J=A*W*B
    Pi=sp.eye(3)-W*B*J.inv()*A
    Q=-Xout+BC*AC/2
    H=Pi*W*Q
    v=sp.Matrix([2,5,7]);b=Pi*W*v
    detK=sp.expand(K.det())
    adj=K.adjugate()
    y=adj[2:,2:]*v
    assert (A*y).applyfunc(sp.expand)==sp.zeros(1,1)
    assert ((nu*sp.eye(3)+H)*y-detK*b).applyfunc(sp.expand)==sp.zeros(3,1)
    coefficients=sp.Matrix.hstack(*[
        y.applyfunc(lambda f:sp.expand(f).coeff(nu,k)) for k in range(2)])
    assert coefficients.rank()==2
    assert all(sp.degree(sp.expand(f),nu)<=1 for f in y)
    return True


def diagonal_adjustment_check():
    """Different-phase vertices split an otherwise repeated zero spectrum."""
    # New phase e=0, three old same-slope phases e_C=1, outside sizes 4,1.
    W=sp.diag(sp.Rational(1,2),sp.Rational(1,2),
              sp.Rational(1,2),sp.Rational(1,2),1)
    AC=sp.zeros(3,5);BC=sp.zeros(5,3);Xout=sp.zeros(5)
    AC[0,1]=1;BC[1,0]=2
    AC[1,2]=1;BC[2,1]=6
    Xout[0,4]=-2;Xout[4,0]=-4
    H=W*(-Xout+BC*AC)
    assert sorted(H.eigenvals())==[-2,0,1,2,3]
    b=W*sp.Matrix([1,2,3,4,5])
    krylov=sp.Matrix.hstack(*[H**k*b for k in range(5)])
    determinant=krylov.det()
    assert determinant!=0
    return str(determinant)


if __name__=="__main__":
    checks=[]
    for sizes,which in [((1,1,1),0),((2,2,1),0),((3,2,2),0),
                        ((3,3,3),1),((4,3,3),2),((4,2,2),0),
                        ((5,2,2),0),((8,2,2),0),((2,7,2),1)]:
        for phases in [(0,0,0),(3,11,29)]:
            row=check(sizes,which,phases)
            checks.append(row)
            print(sizes,which,phases,":",row["spectral_new_rank"],"+",
                  row["moment_new_rank"],"=",row["combined_new_rank"],flush=True)
    for sizes,which,profile,e in [
            ((2,2,1),0,"partial",3),
            ((3,2,2),0,"partial",3),
            ((3,3,3),1,"partial",11),
            ((5,2,2),0,"partial",3),
            ((3,3,3),2,"source",31),
            ((3,3,3),2,"source",29),
            ((5,2,2),0,"partial",101),
            ((4,4,1),0,"distinct",3),
            ((4,4,1),0,"distinct",101)]:
        row=check(sizes,which,(3,11,29),profile=profile,insert_phase=e)
        checks.append(row)
        print(sizes,which,profile,e,":",row["spectral_new_rank"],"+",
              row["moment_new_rank"],"=",row["combined_new_rank"],flush=True)
    symbolic=symbolic_moment_check()
    schur=eliminated_phase_schur_check()
    diagonal=diagonal_adjustment_check()
    output=Path(__file__).parent/"results/jrw_phase_cluster_insertion.json"
    output.write_text(json.dumps(dict(checks=checks,
        symbolic_principal_moment_increment_identity=symbolic,
        different_phase_elimination_polynomial_identity=schur,
        different_phase_diagonal_adjustment_krylov_minor=diagonal),indent=2)+"\n")
