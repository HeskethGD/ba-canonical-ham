"""Exact genus-two mode embedding used by the higher-flow checks."""
import sympy as s

def mode_embedding(nodes, ordinates, p0, p1, j0, j1, phases):
    P = [a*a+p1*a+p0 for a in nodes]
    J = [j1*a+j0 for a in nodes]
    weights = [1/(2*(nodes[0]-nodes[1])), -1/(2*(nodes[0]-nodes[1])),
               1/(2*(nodes[2]-nodes[3])), -1/(2*(nodes[2]-nodes[3]))]
    U, V = [], []
    for i in range(2):
        for j in range(2, 4):
            plus = (P[i]*(ordinates[j]+J[j])-P[j]*(ordinates[i]+J[i]))/(2*(nodes[j]-nodes[i]))
            minus = (P[i]*(ordinates[j]-J[j])-P[j]*(ordinates[i]-J[i]))/(2*(nodes[j]-nodes[i]))
            U.append(s.sqrt(2)*weights[i]*weights[j]*plus/(phases[i]*phases[j]))
            V.append(s.sqrt(2)*minus*phases[i]*phases[j]/(P[i]*P[j]))
    return list(phases)+U+[weights[i]*P[i]/phases[i] for i in range(4)]+V
