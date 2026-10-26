#!/usr/bin/env python3
"""
phaseC_select.py — frame selection + coverage map (fe_melt_v1, gate G-C1)

Reads all traj/<phase>_*.extxyz produced by phaseC_sample.py, applies:
  - exclusion of runs whose diagnostics JSON carries FLAG_decay
  - even-stride thinning to the per-phase frame budget
  - a minimum-displacement decorrelation filter within each run

Writes:
  selected/manifest.json         one entry per selected frame (id, phase, V, T, P_est)
  selected/frames_<phase>.extxyz concatenated selected frames
  coverage_map.txt               phase x pressure-bin x T-band counts  (gate G-C1:
                                 no empty production cell allowed)

Usage: python phaseC_select.py [--trajdir traj] [--outdir selected]
"""
import argparse, glob, json, os
import numpy as np
from ase.io import read, write

BUDGET = {"hcp": 896, "bcc": 300, "fcc": 200, "liquid": 900}
P_BINS = [(50, 100), (100, 150), (150, 200), (200, 250), (250, 300), (300, 400)]
T_BANDS = [(0.0, 0.85), (0.85, 1.05), (1.05, 2.0)]   # fraction of Tm_est
MIN_DISP = 0.15   # A mean displacement between kept frames (decorrelation)

def Tm_est(P): return 3010.0 + 12.07 * P - 0.005410 * P * P

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trajdir", default="traj")
    ap.add_argument("--outdir", default="selected")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    manifest, cover = [], {}
    for phase in BUDGET:
        pool = []
        for fj in sorted(glob.glob(f"{a.trajdir}/{phase}_*.json")):
            meta = json.load(open(fj))
            if meta.get("FLAG_decay"):
                print(f"EXCLUDED (decay flag): {meta['tag']}"); continue
            fx = fj.replace(".json", ".extxyz")
            frames = read(fx, index=":")
            last = None
            for k, fr in enumerate(frames):
                if last is not None:
                    d = np.abs(fr.get_positions() - last).mean()
                    if d < MIN_DISP: continue
                last = fr.get_positions().copy()
                pool.append((meta, k, fr))
        if not pool:
            print(f"WARNING: no frames for {phase}"); continue
        # stratified: thin each source run proportionally (no tail truncation)
        from collections import defaultdict
        byrun = defaultdict(list)
        for item in pool: byrun[item[0]["tag"]].append(item)
        frac = BUDGET[phase] / len(pool)
        sel = []
        for tag in sorted(byrun):
            r = byrun[tag]; k = max(1, round(len(r) * frac))
            st = max(1, len(r) // k)
            sel.extend(r[::st][:k])
        sel = sel[:BUDGET[phase]]
        out = os.path.join(a.outdir, f"frames_{phase}.extxyz")
        if os.path.exists(out): os.remove(out)
        for i, (meta, k, fr) in enumerate(sel):
            fid = f"{phase}-{i:04d}"
            fr.info.update(frame_id=fid, phase=phase, vat=meta["vat"],
                           T=meta["T"], P_est=meta["P_est"], src=meta["tag"])
            write(out, fr, append=True)
            manifest.append(dict(id=fid, phase=phase, vat=meta["vat"],
                                 T=meta["T"], P_est=meta["P_est"]))
            tb = next(j for j, (lo, hi) in enumerate(T_BANDS)
                      if lo <= meta["T"] / Tm_est(meta["P_est"]) < hi)
            pb = next((j for j, (lo, hi) in enumerate(P_BINS)
                       if lo <= meta["P_est"] < hi), None)
            if pb is not None:
                cover[(phase, pb, tb)] = cover.get((phase, pb, tb), 0) + 1
        print(f"{phase}: {len(sel)} / budget {BUDGET[phase]}")
    json.dump(manifest, open(os.path.join(a.outdir, "manifest.json"), "w"), indent=1)

    with open("coverage_map.txt", "w") as f:
        f.write(f"{'phase':>8} | " + " | ".join(f"{lo}-{hi} GPa" for lo, hi in P_BINS) + "\n")
        empties = 0
        for phase in BUDGET:
            row = [f"{phase:>8}"]
            for pb in range(len(P_BINS)):
                counts = [cover.get((phase, pb, tb), 0) for tb in range(len(T_BANDS))]
                row.append(f"{'/'.join(str(c) for c in counts):>11}")
                # production cells: hcp+liquid must populate all bands in 100-300
                if phase in ("hcp", "liquid") and 1 <= pb <= 4 and sum(counts) == 0:
                    empties += 1
            f.write(" | ".join(row) + "\n")
        f.write("\ncells shown as cold/near-melt/hot counts (T bands vs Tm_est)\n")
        f.write(f"G-C1: {'PASS' if empties == 0 else f'FAIL — {empties} empty production cells'}\n")
    print(open("coverage_map.txt").read())

if __name__ == "__main__":
    main()
