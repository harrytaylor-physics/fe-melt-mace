"""cc_check.py - Clausius-Clapeyron check of the coexistence runs.

dV_m per cell: final configurations of the 12,096-atom hcp coexistence cells (runs_tel/*_ckpt.npz),
cut into 42 slabs along z (one two-layer hcp cell per slab); each atom classified solid or liquid
with q-bar-6 (rcut = 3.2 A, the radius the cell thresholds were calibrated with); a slab counts as
solid (liquid) when its solid fraction is > 0.85 (< 0.15) and the two slabs on each side agree.
V per atom = slab volume / atoms in slab. dV_m = <V_liq> - <V_sol>, standard error from both.
dS_m, T_m, P and their standard errors per cell: analysis_live_all_runs_tel.json.

Per pair of neighbouring cells:
  slope of the curve     (T_m2 - T_m1) / (P2 - P1), error hypot(sem1, sem2) / (P2 - P1)
  Clausius-Clapeyron     dV / dS with dV, dS the averages of the two cells,
                         error from the dV standard errors only (dS errors are not stored)
Output: out/cc_check.json {"cells": [...], "intervals": [...]}.
"""
import glob, json, os, sys
import numpy as np
from ase import Atoms
import paths

sys.path.insert(0, paths.Q6_MODULE_DIR)
import phaseE_q6 as q  # noqa: E402

KB = 1.380649e-23
CELLS = ["coex_hcp_%dGPa" % p for p in (100, 150, 200, 250, 300, 330, 350)]
NB, RCUT, F_SOL, F_LIQ = 42, 3.2, 0.85, 0.15


def cell_volumes(cell, thr):
    Vs, Vl = [], []
    files = sorted(glob.glob(os.path.join(paths.PHASEE, "runs_tel", f"{cell}_ti*_ckpt.npz")))
    for f in files:
        d = np.load(f)
        at = Atoms("Fe%d" % len(d["positions"]), positions=d["positions"], cell=d["cell"], pbc=True)
        at.wrap()
        sol = q.q6_per_atom(at, rcut=RCUT) > thr
        Lz = at.cell[2, 2]
        b = (at.positions[:, 2] % Lz / Lz * NB).astype(int) % NB
        cnt = np.bincount(b, minlength=NB)
        fs = np.bincount(b, weights=sol, minlength=NB) / np.maximum(cnt, 1)
        vslab = at.get_volume() / NB
        cls = np.where(fs > F_SOL, 1, np.where(fs < F_LIQ, -1, 0))
        for k in range(NB):
            nbh = [cls[(k + s) % NB] for s in (-2, -1, 0, 1, 2)]
            if all(c == 1 for c in nbh):
                Vs.append(vslab / cnt[k])
            elif all(c == -1 for c in nbh):
                Vl.append(vslab / cnt[k])
    return np.array(Vs), np.array(Vl), len(files)


def main():
    A = json.load(open(os.path.join(paths.PHASEE, "analysis_live_all_runs_tel.json")))["cells"]
    rows = []
    for cell in CELLS:
        thr = json.load(open(os.path.join(paths.PHASEE, "cells", f"{cell}.json")))["q6_threshold"]
        Vs, Vl, nrun = cell_volumes(cell, thr)
        vs, vl = Vs.mean(), Vl.mean()
        se = float(np.hypot(Vs.std() / np.sqrt(len(Vs)), Vl.std() / np.sqrt(len(Vl))))
        c = A[cell]
        r = dict(cell=cell, P=c["P"], Tm=c["Tm"], Tm_sem=c["sem"], dS=c["dS_kB"], Vs=float(vs), Vl=float(vl),
                 dV=float(vl - vs), dV_se=se, dV_pct=float(100 * (vl - vs) / vs), dV_pct_se=float(100 * se / vs),
                 ns=len(Vs), nl=len(Vl), runs=nrun)
        rows.append(r)
        print(f"{cell}: P {r['P']:.1f} GPa, T_m {r['Tm']:.0f} K, V_s {vs:.4f}, V_l {vl:.4f} A^3/atom, "
              f"dV {r['dV']:.4f} +- {se:.4f} ({r['dV_pct']:.2f} +- {r['dV_pct_se']:.2f}%), dS {r['dS']:.3f} kB "
              f"[{len(Vs)}/{len(Vl)} slabs, {nrun} runs]", flush=True)
    iv = []
    print("\nClausius-Clapeyron: slope of T_m(P) between neighbouring cells vs dV/dS")
    for a, b in zip(rows[:-1], rows[1:]):
        dP = b["P"] - a["P"]
        slope, slope_err = (b["Tm"] - a["Tm"]) / dP, float(np.hypot(a["Tm_sem"], b["Tm_sem"]) / dP)
        dV, dV_err = 0.5 * (a["dV"] + b["dV"]), 0.5 * float(np.hypot(a["dV_se"], b["dV_se"]))
        dS = 0.5 * (a["dS"] + b["dS"])
        cc = dV * 1e-30 / (dS * KB) * 1e9
        cc_err = cc * dV_err / dV
        z = (slope - cc) / np.hypot(slope_err, cc_err)
        iv.append(dict(P_lo=a["P"], P_hi=b["P"], curve=slope, curve_err=slope_err, cc=cc, cc_err=cc_err,
                       dV=dV, dS=dS, z=float(z), rel_diff=float((cc - slope) / slope)))
        print(f"  {a['P']:.0f}-{b['P']:.0f} GPa: curve {slope:5.2f} +- {slope_err:.2f}, dV/dS {cc:5.2f} +- {cc_err:.2f} K/GPa, "
              f"difference {z:+.1f} sigma ({100 * (cc - slope) / slope:+.0f}%)")
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(dict(cells=rows, intervals=iv), open(os.path.join(paths.OUT, "cc_check.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
