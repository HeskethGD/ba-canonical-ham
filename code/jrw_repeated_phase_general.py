"""Exact finite-field rank certificates for spectral plus block-chain integrals.

Coordinates are the cross-block entries of X; their frozen bracket has
{X_hi,X_ih}=D_hh-D_ii. All matrices and derivatives are evaluated exactly
modulo a prime. No characteristic polynomial in symbolic modes is expanded.
"""
from pathlib import Path
import json
import numpy as np

PRIME = 1000003


def mm(a, b):
    return (a @ b) % PRIME


def add_row(basis, row):
    row = [int(x) % PRIME for x in row]
    for pivot, old in basis:
        c = row[pivot]
        if c:
            row = [(x-c*y) % PRIME for x, y in zip(row, old)]
    pivot = next((i for i, x in enumerate(row) if x), None)
    if pivot is None:
        return False
    inv = pow(row[pivot], -1, PRIME)
    row = [(x*inv) % PRIME for x in row]
    basis.append((pivot, row))
    basis.sort(key=lambda pair: pair[0])
    return True


def check(sizes, pattern="source", seed=913):
    group = [g for g, n in enumerate(sizes) for _ in range(n)]
    slopes = np.array([1 if g == 0 else -1 if g == 1 else 0 for g in group],
                      dtype=np.int64)
    d = len(group)
    pairs = [(i,j) for i in range(d) for j in range(i+1,d)
             if group[i] != group[j]]
    coords = pairs + [(j,i) for i,j in pairs]
    p = len(pairs)
    rng = np.random.default_rng(seed)
    X = np.zeros((d,d),dtype=np.int64)
    for i,j in coords:
        X[i,j] = int(rng.integers(1,100))
    if pattern == "source":
        phase = [0 if g < 2 else i+3 for i,g in enumerate(group)]
    elif pattern == "all_zero":
        phase = [0]*d
    elif pattern == "three_scalars":
        phase = [3 if g == 0 else 11 if g == 1 else 29 for g in group]
    elif pattern == "partial":
        phase = [7*g+(i%2) for i,g in enumerate(group)]
    else:
        raise ValueError(pattern)
    for i in range(d):
        X[i,i] = phase[i]
    clusters = []
    for g in range(3):
        for phase_value in sorted({phase[i] for i in range(d) if group[i]==g}):
            clusters.append([i for i in range(d)
                             if group[i]==g and phase[i]==phase_value])
    basis=[];selected=[];rows=[];names=[]
    # Traces at d+1 spectral arguments span all coefficient gradients.
    # The coefficient of nu^k in tr(X+nu D)^k is fixed, so d samples suffice.
    for nu in range(d+1):
        pencil=(X+nu*np.diag(slopes)) % PRIME
        power=np.eye(d,dtype=np.int64)
        for k in range(1,d+1):
            row=[(k*int(power[j,i])) % PRIME for i,j in coords]
            name=f"tr(L({nu})^{k})"
            rows.append(row);names.append(name)
            if add_row(basis,row):
                selected.append(name)
            power=mm(power,pencil)
    spectral_rank=len(basis)
    for cluster_index,cluster in enumerate(clusters):
        # J = sum outside X_block,j X_j,block /(slope_block-slope_j).
        r=len(cluster)
        J=np.zeros((r,r),dtype=np.int64)
        deriv=[np.zeros((r,r),dtype=np.int64) for _ in coords]
        for h,i in enumerate(cluster):
            for k,j in enumerate(cluster):
                for v in range(d):
                    if group[v]==group[i]:
                        continue
                    inv=pow(int(slopes[i]-slopes[v]) % PRIME,-1,PRIME)
                    J[h,k]=(int(J[h,k])+int(X[i,v])*int(X[v,j])*inv) % PRIME
                    for z,(c,e) in enumerate(coords):
                        value=((int(X[v,j]) if (c,e)==(i,v) else 0)
                               +(int(X[i,v]) if (c,e)==(v,j) else 0))
                        deriv[z][h,k]=(int(deriv[z][h,k])+value*inv) % PRIME
        # Gelfand--Zeitlin chain: leading r0-by-r0 blocks, powers 1..r0.
        for r0 in range(1,r+1):
            A=J[:r0,:r0]
            power=np.eye(r0,dtype=np.int64)
            for k in range(1,r0+1):
                row=[int(k*np.sum(power.T*B[:r0,:r0])) % PRIME for B in deriv]
                name=f"cluster{cluster_index}:tr(J_{r0}^{k})"
                rows.append(row);names.append(name)
                if add_row(basis,row):
                    selected.append(name)
                power=mm(power,A)
    # Check the full evaluated Poisson matrix as an implementation check.
    R=np.array(rows,dtype=np.int64)
    weights=np.array([int(slopes[i]-slopes[j]) % PRIME for i,j in pairs],
                     dtype=np.int64)
    poisson=(mm((R[:,:p]*weights) % PRIME,R[:,p:].T)
             -mm((R[:,p:]*weights) % PRIME,R[:,:p].T)) % PRIME
    assert not np.any(poisson), (sizes,pattern,"noncommuting candidate gradients")
    assert len(basis)<=p
    orbit_half_dimensions=[]
    for cluster in clusters:
        r=len(cluster)
        outside=d-sizes[group[cluster[0]]]
        k=min(r,outside)
        orbit_half_dimensions.append(r*k-k*(k+1)//2)
    predicted_spectral_rank=p-sum(orbit_half_dimensions)
    assert spectral_rank==predicted_spectral_rank
    return dict(sizes=list(sizes),pattern=pattern,prime=PRIME,seed=seed,
                diagonal_intercepts=phase,
                coordinate_order=[[i+1,j+1] for i,j in coords],
                coordinate_values=[int(X[i,j]) for i,j in coords],
                independent_coordinate_columns=[pivot for pivot,_ in basis],
                pairs=p,spectral_rank=spectral_rank,
                predicted_spectral_rank=predicted_spectral_rank,
                moment_orbit_half_dimensions=orbit_half_dimensions,
                spectral_plus_block_chain_rank=len(basis),
                full_count=len(basis)==p,
                evaluated_brackets_zero=True,clusters=clusters,
                selected_integrals=selected)


if __name__ == "__main__":
    shapes=[(1,1,1),(2,2,1),(3,2,1),(3,3,1),(4,2,1),(4,4,1),
            (5,2,1),(6,1,1),(2,2,2),(3,2,2),(3,3,2),(4,2,3)]
    results=[];second_seed=[]
    for sizes in shapes:
        for pattern in ["source","all_zero","three_scalars","partial"]:
            row=check(sizes,pattern)
            results.append(row)
            print(sizes,pattern,row["spectral_rank"],"->",
                  row["spectral_plus_block_chain_rank"],"/",row["pairs"],flush=True)
        for pattern in ["all_zero","three_scalars"]:
            second_seed.append(check(sizes,pattern,seed=1367))
    output=Path(__file__).parent/"results/jrw_repeated_phase_general.json"
    output.write_text(json.dumps(dict(checks=results,second_seed_checks=second_seed,
        scope="Finite-field certificates for listed cases, not an all-size independence proof."),
        indent=2)+"\n")
