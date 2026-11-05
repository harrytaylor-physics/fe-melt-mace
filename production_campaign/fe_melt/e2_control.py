#!/usr/bin/env python3
"""
e2_control.py — E2 comparator predictions (PRODUCTION env, python3.13, GPU node)

Evaluates the three fixed-Tel specialists on the held-out Tel-6000 test set:
  - fe_ft_Tel6000_f32 : the ceiling (specialist trained on this data's Tel)
  - fe_ft_Tel4000_f32, fe_ft_Tel8000_f32 : legs of the naive control
  - control = mean of the 4000/8000 predictions (energy and forces)
Saves per-frame predictions to phaseD/e2_control.npz for e2_eval.py (tel_env)
and prints per-phase force RMSE for each comparator.

Run: python3.13 e2_control.py --device cuda
"""
import argparse
import numpy as np
from ase.io import read

ap = argparse.ArgumentParser()
ap.add_argument("--device", default="cuda")
a = ap.parse_args()

from mace.calculators import MACECalculator

frames = read("phaseD/test_Tel6000.xyz", index=":")
phases = [f.info.get("phase", "?") for f in frames]
eref = np.array([f.info["REF_energy"] / len(f) for f in frames])
fref = [f.arrays["REF_forces"] for f in frames]

models = {T: MACECalculator(model_paths=f"phaseD/models/fe_ft_Tel{T}_f32.model",
                            device=a.device, default_dtype="float64")
          for T in (4000, 6000, 8000)}

E = {T: [] for T in models}
F = {T: [] for T in models}
for fr in frames:
    for T, c in models.items():
        g = fr.copy(); g.calc = c
        F[T].append(g.get_forces())
        E[T].append(g.get_potential_energy() / len(g))

def rmse_by_phase(fpred):
    out = {}
    for ph in sorted(set(phases)):
        ds = np.concatenate([fpred[i] - fref[i]
                             for i in range(len(frames)) if phases[i] == ph])
        out[ph] = 1000 * np.sqrt((ds ** 2).mean())
    ds = np.concatenate([fpred[i] - fref[i] for i in range(len(frames))])
    out["ALL"] = 1000 * np.sqrt((ds ** 2).mean())
    return out

Fmean = [(F[4000][i] + F[8000][i]) / 2 for i in range(len(frames))]
Emean = (np.array(E[4000]) + np.array(E[8000])) / 2

print(f"{'comparator':22s} " + " ".join(f"{p:>10s}" for p in sorted(set(phases))) + f" {'ALL':>10s}")
for name, fp in [("specialist Tel6000", F[6000]),
                 ("leg Tel4000", F[4000]),
                 ("leg Tel8000", F[8000]),
                 ("mean(4000,8000)", Fmean)]:
    r = rmse_by_phase(fp)
    print(f"{name:22s} " + " ".join(f"{r[p]:10.1f}" for p in sorted(set(phases))) + f" {r['ALL']:10.1f}")

e_spec = 1000 * np.sqrt(((np.array(E[6000]) - eref) ** 2).mean())
e_mean = 1000 * np.sqrt(((Emean - eref) ** 2).mean())
print(f"\nenergy RMSE meV/atom: specialist {e_spec:.2f}   mean-control {e_mean:.2f}")

np.savez("phaseD/e2_control.npz",
         phases=np.array(phases), eref=eref,
         e6000=np.array(E[6000]), e_mean=Emean,
         f6000=np.array(F[6000], dtype=object),
         f_mean=np.array(Fmean, dtype=object),
         fref=np.array(fref, dtype=object), allow_pickle=True)
print("saved phaseD/e2_control.npz")
