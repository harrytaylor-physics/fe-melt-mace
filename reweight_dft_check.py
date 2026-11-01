#!/usr/bin/env python3
"""
reweight_dft_check.py — does the reweighting's first-order shift come from
physics or from differences between the three independently trained models?

The reweighting uses dU = F_model(Tel) - F_model(6000) on configurations from
bulk MD. The quantity that sets dTm is the liquid-minus-solid difference of
<dU>. DFT gives the same quantity directly on the 2,009 frames labelled at all
three Tels with identical geometry. This compares the two, pressure bin by
pressure bin (frames within +/-25 GPa of each measured melting pressure,
hcp and liquid only).

If model and DFT agree within errors, the first-order means are physical and
the large second-order term must be examined separately. If they disagree,
the models' own errors contaminate the reweighting.

Run (login node): python3.13 reweight_dft_check.py
"""
import glob, hashlib, json
import numpy as np
from ase.io import read

DS = 0.84 * 8.617333e-2      # meV/atom/K
CELLS = [100, 150, 200, 250, 300, 330, 350]


def key(a):
    h = hashlib.sha1()
    h.update(np.round(a.get_positions(), 4).tobytes())
    h.update(np.round(np.array(a.cell), 4).tobytes())
    return h.hexdigest()


def load(tel):
    out = {}
    for f in sorted(glob.glob(f"phaseD/*_Tel{tel}.xyz")):
        for a in read(f, index=":"):
            out[key(a)] = (a.info.get("phase"), float(a.info.get("P_est", np.nan)),
                           a.info["REF_energy"] / len(a))
    return out


d = {T: load(T) for T in (4000, 6000, 8000)}
common = set(d[4000]) & set(d[6000]) & set(d[8000])
cells = json.load(open("phaseE/analysis_live_all.json"))["cells"]

print(f"liquid-minus-hcp of F(Tel) - F(6000), meV/atom  (frames labelled at all Tels: {len(common)})")
print(f"{'cell':8s} {'P':>6s} {'n_hcp':>5s} {'n_liq':>5s} | {'DFT 4000':>12s} {'model 4000':>10s} | "
      f"{'DFT 8000':>12s} {'model 8000':>10s} | {'DFT dTm 4000/8000 (K)':>22s}")
for c in CELLS:
    name = f"coex_hcp_{c}GPa"
    P0 = cells[name]["P"]
    rows = {"hcp": {4000: [], 8000: []}, "liquid": {4000: [], 8000: []}}
    for k in common:
        ph, P, e6 = d[6000][k]
        if ph in rows and abs(P - P0) <= 25.0:
            for T in (4000, 8000):
                rows[ph][T].append(1000 * (d[T][k][2] - e6))
    nh, nl = len(rows["hcp"][4000]), len(rows["liquid"][4000])
    try:
        rw = json.load(open(f"phaseE/reweight/{name}.json"))
        model = {T: rw["liquid"][f"dU_{T}"] - rw["solid"][f"dU_{T}"] for T in (4000, 8000)}
    except FileNotFoundError:
        model = {4000: np.nan, 8000: np.nan}
    if nh < 3 or nl < 3:
        print(f"{name[9:]:8s} {P0:6.1f} {nh:5d} {nl:5d} | too few frames in this pressure bin")
        continue
    out, dtm = [], []
    for T in (4000, 8000):
        h, l = np.array(rows["hcp"][T]), np.array(rows["liquid"][T])
        diff = l.mean() - h.mean()
        se = np.hypot(l.std(ddof=1) / np.sqrt(len(l)), h.std(ddof=1) / np.sqrt(len(h)))
        out.append(f"{diff:+6.2f}+-{se:4.2f} {model[T]:+10.2f}")
        dtm.append(f"{diff / DS:+6.0f}+-{se / DS:3.0f}")
    print(f"{name[9:]:8s} {P0:6.1f} {nh:5d} {nl:5d} | {out[0]} | {out[1]} | {dtm[0]:>11s} {dtm[1]:>10s}")
print("\nDFT frames are training/test geometries near (not at) coexistence; the model column is"
      "\nthe reweighting value from bulk MD at the measured (Tm, P).")
