"""Spectral-evaluation independence for the symmetric three-group system.

Use generic symbolic diagonal intercepts for (2,2,2).  At X_upper=0,
X_lower=epsilon, the first nonzero differential coefficient is diagonal
in pairwise spectral-line evaluations.  Power charges give a weighted
incidence matrix.  A spanning tree selects charges plus G evaluations.
No characteristic polynomial in all modes is expanded.
"""
from pathlib import Path
import json
import sympy as s


def edges_for(sizes):
    group = [g for g, count in enumerate(sizes) for _ in range(count)]
    return group, [(a, b) for a in range(len(group))
                   for b in range(a + 1, len(group)) if group[a] != group[b]]


def tree_for(sizes, edges):
    L, M, N = sizes
    # Retain every A-B evaluation: connect all A/B vertices to C_1,
    # then connect the remaining C vertices to A_1.
    return ([(a, L + M) for a in range(L + M)]
            + [(0, L + M + n) for n in range(1, N)])


def incidence(group, edges):
    d = len(group)
    B = s.zeros(d, len(edges))
    # Row signs convert the charge gradients to oriented incidence.
    # A-B power = X_ab X_ba / 2; A-C, B-C power = -X_ab X_ba.
    for col, (a, b) in enumerate(edges):
        weight = {(0, 1): s.Rational(1, 2),
                  (0, 2): s.S.One, (1, 2): -s.S.One}[(group[a], group[b])]
        B[a, col], B[b, col] = -weight, weight
    return B


def check(sizes, symbolic=False):
    group, edges = edges_for(sizes)
    d, p = sum(sizes), len(edges)
    slopes = [1 if g == 0 else -1 if g == 1 else 0 for g in group]
    if symbolic:
        intercept = s.symbols(f'a0:{d}')
    else:
        # Distinct within blocks and no three concurrent cross-group lines.
        starts = [s.Integer(3), s.Integer(41), s.Integer(211)]
        intercept = tuple(starts[g] + 7 * i for g, count in enumerate(sizes)
                          for i in range(count))
    eps = s.Symbol('epsilon')
    X = s.diag(*intercept)
    for a, b in edges:
        X[b, a] = eps
    D = s.diag(*slopes)
    derivative = s.zeros(p)
    leading_diagonal = []
    intersections = []
    for row, (a, b) in enumerate(edges):
        nu = (intercept[b] - intercept[a]) / (slopes[a] - slopes[b])
        lam = intercept[a] + slopes[a] * nu
        intersections.append((lam, nu))
        A = lam * s.eye(d) - X - nu * D
        target = -s.prod(lam - slopes[k] * nu - intercept[k]
                         for k in range(d) if k not in (a, b))
        leading_diagonal.append(target)
        assert target != 0
        for col, (c, e) in enumerate(edges):
            # d det(A) / d X_ce = -cofactor_ce(A).
            value = -A.cofactor(c, e)
            derivative[row, col] = value
            assert s.expand(value).coeff(eps, 1).equals(target if row == col else 0)
            assert s.expand(value).coeff(eps, 0) == 0
    tree = tree_for(sizes, edges)
    chords = [e for e in edges if e not in tree]
    tree_cols = [edges.index(e) for e in tree]
    chord_rows = [edges.index(e) for e in chords]
    B = incidence(group, edges)[:-1, :]
    actual_charges = s.zeros(d, p)
    for col, (a, b) in enumerate(edges):
        ga, gb = group[a], group[b]
        if (ga, gb) == (0, 1):
            actual_charges[a, col] = actual_charges[b, col] = -s.Rational(1, 2)
        elif (ga, gb) == (0, 2):
            actual_charges[a, col] = actual_charges[b, col] = -1
        else:
            actual_charges[a, col], actual_charges[b, col] = -1, 1
    signs = s.diag(*([1] * sizes[0] + [-1] * (sizes[1] + sizes[2])))
    assert signs * actual_charges == incidence(group, edges)
    tree_minor = B[:, tree_cols].det()
    assert tree_minor != 0
    genus = p - d + 1
    assert len(chords) == genus
    # Factor each nonzero spectral diagonal from its row first. The
    # remaining minor is constant, so there is no large symbolic det.
    mixed_normalized = B.col_join(s.eye(p)[chord_rows, :])
    assert mixed_normalized.det() != 0
    numeric_check = None
    if not symbolic:
        actual = derivative.subs(eps, s.Rational(1, 97))
        mixed = (eps * B).col_join(derivative[chord_rows, :]).subs(eps, s.Rational(1, 97))
        assert actual.det() != 0
        assert mixed.det() != 0
        numeric_check = dict(all_evaluation_rank=actual.rank(),
                             charges_plus_chord_evaluation_rank=mixed.rank())
    return dict(sizes=sizes,canonical_pairs=p,independent_charges=d-1,
                required_spectral_integrals=genus,
                symbolic_intercepts=symbolic,
                all_evaluation_leading_matrix_diagonal=True,
                spanning_tree_edges=[[a+1,b+1] for a,b in tree],
                spectral_chord_edges=[[a+1,b+1] for a,b in chords],
                weighted_tree_minor=str(tree_minor),
                complete_leading_minor_nonzero=True,
                exact_rational_witness=numeric_check)


if __name__ == '__main__':
    report = dict(
        genericity='Distinct intercepts within each group; no cross-group triple concurrence.',
        checks=[check((2,2,2), symbolic=True)] +
               [check(shape) for shape in [(1,1,1),(2,2,1),(1,2,2),(2,2,2),(2,2,3),(3,2,2)]],
        scope='The all-size argument is in the manuscript; these are explicit spectral-evaluation and spanning-tree checks.')
    out = Path(__file__).with_name('results') / 'jrw_three_group_integral_count.json'
    out.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))
