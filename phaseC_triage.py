#!/usr/bin/env python3
"""
phaseC_triage.py — salvage flagged hcp runs (fe_melt_v1)

For each run whose JSON carries FLAG_decay, recompute the per-frame mean-square
displacement from the first stored frame and split the trajectory at the first
frame where MSD > CUT:
  frames before the crossing -> traj/<tag>_pre.extxyz, new JSON, FLAG_decay
     false, phase unchanged (valid thermal hcp; enters selection normally)
  frames at/after            -> traj/<tag>_sf.extxyz, phase 'stackfault'
     (kept for optional labelling as a stacking-energetics class; NOT counted
     toward hcp coverage)

Run AFTER sampling, BEFORE phaseC_select.py. Original flagged JSONs are
renamed .json -> .json.flagged so selection ignores them.

Usage: python phaseC_triage.py
"""
import glob, json, os
import numpy as np
from ase.io import read, write

CUT = 0.45   # A^2 vs frame 0; below the 0.6 monitor limit for margin

n_pre = n_sf = 0
for fj in sorted(glob.glob("traj/*_*.json")):
    meta = json.load(open(fj))
    if not meta.get("FLAG_decay"):
        continue
    tag = meta["tag"]
    frames = read(fj.replace(".json", ".extxyz"), index=":")
    r0 = frames[0].get_positions()
    msd = [float(((f.get_positions() - r0) ** 2).sum(axis=1).mean()) for f in frames]
    cross = next((i for i, m in enumerate(msd) if m > CUT), len(frames))
    pre, post = frames[:cross], frames[cross:]
    if len(pre) >= 5:
        fx = f"traj/{tag}_pre.extxyz"
        if os.path.exists(fx): os.remove(fx)
        for f in pre: write(fx, f, append=True)
        m2 = dict(meta); m2.update(tag=f"{tag}_pre", frames=len(pre),
                                   msd_final=msd[max(cross - 1, 0)], FLAG_decay=False,
                                   note="pre-decay salvage")
        json.dump(m2, open(f"traj/{tag}_pre.json", "w"), indent=1)
        n_pre += len(pre)
    if len(post) >= 5:
        fx = f"traj/{tag}_sf.extxyz"
        if os.path.exists(fx): os.remove(fx)
        for f in post: write(fx, f, append=True)
        m3 = dict(meta); m3.update(tag=f"{tag}_sf", phase="stackfault",
                                   frames=len(post), FLAG_decay=False,
                                   note="post-decay stacking-fault class")
        json.dump(m3, open(f"traj/{tag}_sf.json", "w"), indent=1)
        n_sf += len(post)
    os.rename(fj, fj + ".flagged")
    print(f"{tag}: crossing at frame {cross}/{len(frames)} "
          f"-> {len(pre)} hcp frames, {len(post)} stackfault frames")
print(f"\nsalvaged {n_pre} hcp frames; {n_sf} stacking-fault frames set aside")
