#!/usr/bin/env python3
"""
phaseE_patch_mem.py — fix the two memory faults in the Phase E scripts.

1. q6_per_atom in phaseE_build.py built an N x N pairwise tensor. At 12,096
   atoms that is ~1.75 GB per intermediate and several exist at once, so it
   OOMs on a 44 GB card. Replaced with an ASE neighbour-list version
   (cell-list, linear scaling); phaseE_coex.py imports the same function, so
   both are fixed by this one change.
2. The calculators asked for default_dtype="float64" while the production
   models are float32, forcing an upconversion that doubles activation memory
   on a 12k-atom graph. Switched to float32 in the Phase E scripts (the gate
   scripts keep float64 — they run on small cells where it is free).

Run once, in ~/fe_melt:  python3.13 phaseE_patch_mem.py
"""
import re, sys

NEW_Q6 = '''def q6_per_atom(at, rcut=3.2, device="cpu"):
    """Mean-|q6m| local bond order per atom, neighbour-list based (linear
    scaling; the previous N x N formulation exhausted GPU memory at 12k atoms).
    Pure numpy/scipy on the host: called every few ps, not every step."""
    import numpy as np
    from ase.neighborlist import neighbor_list
    from scipy.special import sph_harm
    n = len(at)
    i, D, r = neighbor_list("iDd", at, rcut)
    if len(i) == 0:
        return np.zeros(n)
    order = np.argsort(i, kind="stable")
    i, D, r = i[order], D[order], r[order]
    starts = np.searchsorted(i, np.arange(n), side="left")
    ends = np.searchsorted(i, np.arange(n), side="right")
    out = np.zeros(n)
    for a in range(n):
        s, e = starts[a], ends[a]
        if e - s < 4:
            continue
        v, rr = D[s:e], r[s:e]
        th = np.arccos(np.clip(v[:, 2] / rr, -1.0, 1.0))
        ph = np.arctan2(v[:, 1], v[:, 0])
        qs = 0.0
        for m in range(-6, 7):
            y = sph_harm(m, 6, ph, th).mean()
            qs += abs(y) ** 2
        out[a] = np.sqrt(4 * np.pi / 13 * qs)
    return out
'''

src = open("phaseE_build.py").read()
start = src.index("def q6_per_atom")
end = src.index("def npt_lattice")
src = src[:start] + NEW_Q6 + "\n" + src[end:]
src = src.replace('default_dtype="float64"', 'default_dtype="float32"')
open("phaseE_build.py", "w").write(src)
print("phaseE_build.py: q6 replaced, calculator dtype -> float32")

for f in ("phaseE_coex.py", "phaseE_pin.py"):
    try:
        s = open(f).read()
        s2 = s.replace('default_dtype="float64"', 'default_dtype="float32"')
        open(f, "w").write(s2)
        print(f"{f}: calculator dtype -> float32" if s != s2 else f"{f}: no change needed")
    except FileNotFoundError:
        print(f"{f}: not found, skipped")

print("\nAlso export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True before running")
print("to limit fragmentation on the 44 GB card.")
