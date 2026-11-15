"""Symbolic N=1 pencil and geometric unit-bracket reference model."""
import sympy as s

def build(m, n):
    count = m + n
    edges = [(i, j) for i in range(m) for j in range(m, count)]
    u = s.symbols(f'u1:{count+1}')
    v = s.symbols(f'v1:{count+1}')
    U = s.symbols(' '.join(f'u{i+1}_{j+1}' for i, j in edges), seq=True)
    V = s.symbols(' '.join(f'v{i+1}_{j+1}' for i, j in edges), seq=True)
    b = s.symbols(f'b1:{count+1}')
    c = s.symbols(f'c1:{count+1}')
    nu = s.Symbol('nu')
    z = u + U + v + V
    pairs = count + len(edges)

    def pb(f, g):
        return s.expand(sum(s.diff(f, z[k]) * s.diff(g, z[pairs+k])
                            - s.diff(f, z[pairs+k]) * s.diff(g, z[k])
                            for k in range(pairs)))

    H = -sum(b[i]*u[i]*v[i] for i in range(count))/2
    H += sum((c[i]+c[j])*U[k]*V[k] for k, (i, j) in enumerate(edges))
    H += s.sqrt(2)*sum(u[i]*u[j]*U[k]+v[i]*v[j]*V[k]
                       for k, (i, j) in enumerate(edges))
    charges = [u[i]*v[i]-sum(U[k]*V[k] for k, e in enumerate(edges) if i in e)
               for i in range(count)]
    d = count + 1
    X = s.zeros(d)
    D = s.diag(*([1]*m+[-1]*n+[0]))
    T = s.zeros(d)
    for i in range(count):
        X[i, i] = b[i]/2-c[i]
        T[i, i] = -c[i] if i < m else c[i]
    for i in range(m):
        X[i, count] = T[i, count] = v[i]
        X[count, i] = T[count, i] = -u[i]
    for j in range(m, count):
        X[j, count] = u[j]
        T[j, count] = -u[j]
        X[count, j] = -v[j]
        T[count, j] = v[j]
    for k, (i, j) in enumerate(edges):
        X[i, j] = s.sqrt(2)*U[k]
        X[j, i] = s.sqrt(2)*V[k]
    T[count, count] = -nu
    return locals()
