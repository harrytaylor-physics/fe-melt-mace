#!/usr/bin/env python3
"""
liquid_diag.py — where does the liquid force error live? (fe_melt_v1)

Groups held-out liquid test frames by sampling volume and temperature band and
reports force RMSE per group, plus the mean |F| in that group so absolute and
relative error can be compared. Decides whether a G-D2p liquid top-up should be
broad (uniform error -> more data everywhere) or narrow (error concentrated at
specific state points -> targeted resampling).

Usage (GPU node):
  python3.13 liquid_diag.py --model phaseD/models/fe_ft_Tel6000_f32.model --tel 6000
"""
import argparse
from collections import defaultdict
import numpy as np
from ase.io import read

ap = argparse.ArgumentParser()
ap.add_argument("--model", required=True)
ap.add_argument("--tel", type=int, default=6000)
ap.add_argument("--device", default="cuda")
a = ap.parse_args()

from mace.calculators import MACECalculator
calc = MACECalculator(model_paths=a.model, device=a.device, default_dtype="float64")

frames = [f for f in read(f"phaseD/test_Tel{a.tel}.xyz", index=":")
          if f.info.get("phase") == "liquid"]
if not frames:
    raise SystemExit("no liquid frames in the test set")

groups = defaultdict(list)
for fr in frames:
    vat = fr.info.get("vat", 0.0)
    P = fr.info.get("P_est", 0.0)
    T = fr.info.get("T", 0.0)
    fr2 = fr.copy(); fr2.calc = calc
    d = fr2.get_forces() - fr.arrays["REF_forces"]
    groups[(round(float(vat), 2), round(float(P)), round(float(T), -2))].append(
        (d, fr.arrays["REF_forces"]))

print(f"{'V':>6} {'P/GPa':>7} {'T/K':>7} {'n':>3} {'F_RMSE':>9} {'mean|F|':>9} {'rel%':>6}")
rows = []
for k in sorted(groups):
    ds = np.concatenate([g[0] for g in groups[k]])
    fs = np.concatenate([g[1] for g in groups[k]])
    rmse = 1000 * np.sqrt((ds ** 2).mean())
    meanf = 1000 * np.sqrt((fs ** 2).mean())
    rows.append((k, len(groups[k]), rmse, meanf))
    print(f"{k[0]:6.2f} {k[1]:7d} {k[2]:7.0f} {len(groups[k]):3d} "
          f"{rmse:9.1f} {meanf:9.1f} {100*rmse/meanf:6.1f}")

allrmse = 1000 * np.sqrt(np.concatenate([g[0] for gs in groups.values() for g in gs]) ** 2).mean() ** 0.5
worst = sorted(rows, key=lambda r: -r[2])[:3]
print("\nworst groups:", ", ".join(f"V={r[0][0]} P={r[0][1]} T={r[0][2]:.0f} ({r[2]:.0f})" for r in worst))
print("Interpretation: spread < ~25 meV/A between best and worst -> uniform,")
print("broad top-up; a clear tail -> targeted resampling at those state points.")
