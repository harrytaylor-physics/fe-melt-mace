#!/usr/bin/env python3
"""
tel_diff_check.py — do the potentials reproduce how the DFT free energy
changes with electronic temperature, on IDENTICAL frames?

Uses the 2,009 frames labelled by DFT at Tel = 4000, 6000 and 8000 K with the
same geometry. For each frame and Tel in (4000, 8000):
    dU_DFT   = F_DFT(Tel)   - F_DFT(6000)
    dU_model = F_model(Tel) - F_model(6000)
Reports the per-frame error of dU_model, its per-phase means, and the
liquid-minus-hcp difference in pressure bins around each melting pressure,
for the model and for DFT side by side. Identical frames remove the
bin-composition caveat of reweight_dft_check.py: any disagreement here is the
model's error.

Two modes:
  three fixed-Tel specialists (production env):
      python3.13 tel_diff_check.py --specialists
  one Tel-conditioned model (tel_env):
      python tel_diff_check.py --tel-model phaseD/models/fe_tel_e3_stagetwo.model
"""
import argparse, glob, hashlib, json
import numpy as np
from ase.io import read

TELS = (4000, 6000, 8000)
CELLS = [100, 150, 200, 250, 300, 330, 350]
DS = 0.84 * 8.617333e-2      # meV/atom/K

ap = argparse.ArgumentParser()
ap.add_argument("--specialists", action="store_true")
ap.add_argument("--tel-model")
ap.add_argument("--device", default="cuda")
a = ap.parse_args()
from mace.calculators import MACECalculator


def key(at):
    h = hashlib.sha1()
    h.update(np.round(at.get_positions(), 4).tobytes())
    h.update(np.round(np.array(at.cell), 4).tobytes())
    return h.hexdigest()


frames = {T: {} for T in TELS}
for T in TELS:
    for f in sorted(glob.glob(f"phaseD/*_Tel{T}.xyz")):
        for at in read(f, index=":"):
            frames[T][key(at)] = at
common = sorted(set(frames[4000]) & set(frames[6000]) & set(frames[8000]))
print(f"frames labelled at all three Tels: {len(common)}")

if a.specialists:
    calcs = {T: MACECalculator(model_paths=f"phaseD/models/fe_ft_Tel{T}_f32.model",
                               device=a.device, default_dtype="float32") for T in TELS}
    label = "three fixed-Tel specialists"
else:
    c = MACECalculator(model_paths=a.tel_model, device=a.device, default_dtype="float32")
    calcs = {T: c for T in TELS}
    label = f"Tel model {a.tel_model}"

ph, P, dD, dM = [], [], {4000: [], 8000: []}, {4000: [], 8000: []}
for k in common:
    at6 = frames[6000][k]
    ph.append(at6.info.get("phase")); P.append(float(at6.info.get("P_est", np.nan)))
    E = {}
    for T in TELS:
        g = at6.copy(); g.info["elec_temp"] = float(T); g.calc = calcs[T]
        E[T] = g.get_potential_energy() / len(g)
    for T in (4000, 8000):
        dD[T].append(1000 * (frames[T][k].info["REF_energy"] - at6.info["REF_energy"]) / len(at6))
        dM[T].append(1000 * (E[T] - E[6000]))
ph, P = np.array(ph), np.array(P)
dD = {T: np.array(v) for T, v in dD.items()}
dM = {T: np.array(v) for T, v in dM.items()}

print(f"\n{label}")
print("per-frame error of dU (model - DFT), meV/atom:")
for T in (4000, 8000):
    e = dM[T] - dD[T]
    print(f"  Tel {T}: mean {e.mean():+6.2f}  rmse after mean {e.std():5.2f}   per phase: " +
          "  ".join(f"{p} {e[ph == p].mean():+6.2f}" for p in sorted(set(ph)) if p))

cells = json.load(open("phaseE/analysis_live_all.json"))["cells"]
print(f"\nliquid-minus-hcp of dU, same frames, +/-25 GPa bins (meV/atom; K at dS = 0.84 kB)")
print(f"{'cell':7s} {'n_h':>4s} {'n_l':>4s} | {'DFT 4000':>13s} {'model 4000':>13s} | {'DFT 8000':>13s} {'model 8000':>13s}")
for cc in CELLS:
    P0 = cells[f"coex_hcp_{cc}GPa"]["P"]
    m = np.abs(P - P0) <= 25.0
    h, l = m & (ph == "hcp"), m & (ph == "liquid")
    if h.sum() < 3 or l.sum() < 3:
        print(f"{cc:4d}GPa {h.sum():4d} {l.sum():4d} | too few frames of one phase")
        continue
    out = []
    for T in (4000, 8000):
        for src in (dD[T], dM[T]):
            d = src[l].mean() - src[h].mean()
            se = np.hypot(src[l].std(ddof=1) / np.sqrt(l.sum()), src[h].std(ddof=1) / np.sqrt(h.sum()))
            out.append(f"{d:+6.2f}+-{se:4.2f}")
    print(f"{cc:4d}GPa {h.sum():4d} {l.sum():4d} | {out[0]:>13s} {out[1]:>13s} | {out[2]:>13s} {out[3]:>13s}")
