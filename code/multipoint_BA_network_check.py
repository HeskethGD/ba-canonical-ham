"""Exact genus-zero test of the prescribed-curve multipuncture construction.

Tests all 24 mode rates for L=M=N=2, plus kernel-rate and convolution
identities. The arbitrary-genus argument is in the accompanying note.
No AHH reconstruction is used.
"""
import json
from pathlib import Path
import sympy as s


def check():
    P,Q=s.symbols('P Q')
    R=[-s.Rational(2,5),s.Rational(3,7)]
    D=s.Rational(2,3);dual=sum(R)-D
    C=(R[0]-R[1])**2/((R[0]-D)*(R[0]-dual))
    beta=lambda x:-C*(x-D)*(x-dual)/s.prod((x-r)**2 for r in R)
    chi=lambda x:sum(1/(x-r) for r in R)
    reg=[1/(R[k]-R[1-k]) for k in range(2)]
    psi=lambda x,k:(R[k]-D)/(R[k]-R[1-k])*(x-R[1-k])/(x-D)
    psid=lambda x,k:(R[k]-dual)/(R[k]-R[1-k])*(x-R[1-k])/(x-dual)
    kernel=lambda p,q:-s.prod(p-r for r in R)*s.prod(q-r for r in R)/(C*(p-D)*(q-dual)*(p-q))
    assert s.factor((chi(P)-chi(Q))*kernel(P,Q)
                    -sum(psi(P,k)*psid(Q,k) for k in range(2)))==0
    nodes=[-s.Rational(7,6),-s.Rational(1,8),s.Rational(6,5),s.Rational(11,7)]
    weights=[s.Rational(2,7),-s.Rational(3,8),s.Rational(5,9),s.Rational(4,11)]
    kap=[beta(nodes[i])*weights[i]/2 for i in range(4)]
    hL=lambda x:chi(x)+sum(weights[i]/(x-nodes[i]) for i in range(2))
    hM=lambda x:chi(x)+sum(weights[j]/(x-nodes[j]) for j in range(2,4))
    aL=[reg[k]+sum(weights[i]/(R[k]-nodes[i]) for i in range(2)) for k in range(2)]
    aM=[reg[k]+sum(weights[j]/(R[k]-nodes[j]) for j in range(2,4)) for k in range(2)]
    for group,h in [(range(2),hL),(range(2,4),hM)]:
        residual=(h(P)-h(Q))*kernel(P,Q)
        residual-=2*sum(kap[i]*kernel(P,nodes[i])*kernel(nodes[i],Q) for i in group)
        residual-=sum(psi(P,k)*psid(Q,k) for k in range(2))
        assert s.factor(residual)==0
    bL=[-hM(nodes[i])/2 for i in range(2)]
    bM=[hL(nodes[j])/2 for j in range(2,4)]
    bN=[(aM[k]-aL[k])/2 for k in range(2)]
    cL=bL
    cM=[-x for x in bM]
    cN=[-(aM[k]+aL[k])/2 for k in range(2)]
    root=s.sqrt(2)
    data={};rates={};pairs=[]
    def put(name,rows,cols,value,rate):
        A=s.Matrix(rows,cols,lambda i,j:s.Symbol(f'{name}{i}{j}'))
        for i in range(rows):
            for j in range(cols):
                data[A[i,j]]=value(i,j)
                rates[A[i,j]]=rate(i,j)*value(i,j)
        return A
    u=put('u',2,2,lambda i,j:root*kap[i]*kap[j+2]*kernel(nodes[i],nodes[j+2]),
          lambda i,j:chi(nodes[i])-chi(nodes[j+2])+cL[i]-cM[j])
    ub=put('ub',2,2,lambda i,j:-root*kernel(nodes[j+2],nodes[i]),
           lambda i,j:-(chi(nodes[i])-chi(nodes[j+2])+cL[i]-cM[j]))
    v=put('v',2,2,lambda j,k:psi(nodes[j+2],k),
          lambda j,k:chi(nodes[j+2])-reg[k]+cM[j]-cN[k])
    vb=put('vb',2,2,lambda j,k:kap[j+2]*psid(nodes[j+2],k),
           lambda j,k:-(chi(nodes[j+2])-reg[k]+cM[j]-cN[k]))
    w=put('w',2,2,lambda k,i:psid(nodes[i],k),
          lambda k,i:-(chi(nodes[i])-reg[k])+cN[k]-cL[i])
    wb=put('wb',2,2,lambda k,i:kap[i]*psi(nodes[i],k),
           lambda k,i:-(-(chi(nodes[i])-reg[k])+cN[k]-cL[i]))
    H=sum((bL[i]+bM[j])*u[i,j]*ub[i,j] for i in range(2) for j in range(2))
    H+=sum((bM[j]+bN[k])*v[j,k]*vb[j,k] for j in range(2) for k in range(2))
    H+=sum((bN[k]+bL[i])*w[k,i]*wb[k,i] for k in range(2) for i in range(2))
    H+=root*sum(u[i,j]*v[j,k]*w[k,i]+ub[i,j]*vb[j,k]*wb[k,i]
                for i in range(2) for j in range(2) for k in range(2))
    # Each triad has zero exponential rate, so this exact coefficient check
    # verifies the equations for all t, rather than merely at one initial point.
    for i in range(2):
        for j in range(2):
            for k in range(2):
                assert s.simplify(rates[u[i,j]]/data[u[i,j]]+
                                  rates[v[j,k]]/data[v[j,k]]+
                                  rates[w[k,i]]/data[w[k,i]])==0
    for A,B in [(u,ub),(v,vb),(w,wb)]:
        for a,b in zip(A,B):
            assert s.simplify(s.diff(H,b).subs(data)-rates[a])==0
            assert s.simplify(-s.diff(H,a).subs(data)-rates[b])==0
    return {'genus':0,'L_M_N':[2,2,2],'kernel_rate_identity':True,
            'both_group_convolution_identities':True,'all_24_mode_rates':True,
            'triad_exponents_cancel_for_all_times':True,
            'scope':'Exact rational multipuncture check. Positive-genus '
                    'scope follows from the stated BA uniqueness/residue proof, '
                    'not from this genus-zero example.'}


if __name__=='__main__':
    result=check()
    (Path(__file__).parent/'results'/'multipoint_BA_network_check.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
