"""Check MHDG per-rank decomposition files <mesh>_<rank>_<np>.h5 for consistency.

Usage: python tools/check_domain_decomposition.py "westCore_*_48.h5"
Checks: each element owned by one rank; each face assembled by one rank that holds
all its elements; ghostPro/ghostLoc and ghelsPro/ghelsLoc point to the owner.
"""
import sys, glob, h5py, numpy as np
from collections import defaultdict
pat = sys.argv[1]                      # e.g. run/westCoreP1_*_4.h5
files = sorted(glob.glob(pat), key=lambda s: int(s.rsplit('_', 2)[1]))
R = []
for fn in files:
    with h5py.File(fn) as f:
        d = {k: f[k][()] for k in f.keys() if isinstance(f[k], h5py.Dataset)}
    R.append(d)
P = len(R); errors = []
def err(msg):
    errors.append(msg)
    if len(errors) <= 10: print("  ERROR:", msg)
el_owner = defaultdict(list); fa_owner = defaultdict(list)
for p, d in enumerate(R):
    l2e, l2f = d['loc2glob_el'], d['loc2glob_fa']
    gE, gF = d['ghostElems'], d['ghostFaces']
    for e, g in zip(l2e, gE):
        if g == 0: el_owner[e].append(p)
    T = d['T']; nint = int(d['Nintfaces'])
    for i, (fa, g) in enumerate(zip(l2f, gF)):
        if g == 0: fa_owner[fa].append((p, i))
Nel = max(max(d['loc2glob_el']) for d in R); Nfa = max(max(d['loc2glob_fa']) for d in R)
print(f"{P} ranks, {Nel} elements, {Nfa} faces")
for e in range(1, Nel + 1):
    if len(el_owner[e]) != 1: err(f"element {e} owned by ranks {el_owner[e]}")
for fa in range(1, Nfa + 1):
    if len(fa_owner[fa]) != 1: err(f"face {fa} assembled by {[o[0] for o in fa_owner[fa]]}")
# owner holds all elements adjacent to the face: owner's local face is interior
# (two local elements) or a physical boundary face
for fa, owners in fa_owner.items():
    if len(owners) != 1: continue
    p, i = owners[0]; d = R[p]; nint = int(d['Nintfaces'])
    if i >= nint and d['boundaryFlag'][i - nint] == 0:
        err(f"face {fa} assembled on rank {p} but lies on the border of its local mesh")
# ghost pointers
for p, d in enumerate(R):
    gF = d['ghostFaces']; idx = np.nonzero(gF)[0]
    if len(idx):
        pro, loc = d['ghostPro'], d['ghostLoc']
        for j, i in enumerate(idx):
            q, l = int(pro[j]) - 1, int(loc[j]) - 1
            if q < 0 or q >= P: err(f"rank {p} ghost face {d['loc2glob_fa'][i]}: ghostPro={pro[j]}"); continue
            dq = R[q]
            if dq['loc2glob_fa'][l] != d['loc2glob_fa'][i] or dq['ghostFaces'][l] != 0:
                err(f"rank {p} ghost face {d['loc2glob_fa'][i]} -> rank {q} local {l+1} is wrong")
    gE = d['ghostElems']; idx = np.nonzero(gE)[0]
    if len(idx):
        pro, loc = d['ghelsPro'], d['ghelsLoc']
        for j, i in enumerate(idx):
            q, l = int(pro[j]) - 1, int(loc[j]) - 1
            if q < 0 or q >= P: err(f"rank {p} ghost element {d['loc2glob_el'][i]}: ghelsPro={pro[j]}"); continue
            dq = R[q]
            if dq['loc2glob_el'][l] != d['loc2glob_el'][i] or dq['ghostElems'][l] != 0:
                err(f"rank {p} ghost element {d['loc2glob_el'][i]} -> rank {q} local {l+1} is wrong")
own = [int((d['ghostElems'] == 0).sum()) for d in R]
print(f"owned elements per rank: min {min(own)} max {max(own)};  ghost faces total {sum(int(d['ghostFaces'].sum()) for d in R)}")
print("OK: decomposition consistent" if not errors else f"FAILED: {len(errors)} errors")
