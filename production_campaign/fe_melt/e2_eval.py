#!/usr/bin/env python3
"""
e2_eval.py — E2 judgment (TEL_ENV: ~/tel_env/bin/python, GPU node)

Evaluates the Tel-conditioned model, trained on Tel 4000+8000 labels only, on
the held-out Tel-6000 test set it has never seen, and prints the comparison
table against the comparators saved by e2_control.py (run that first, in the
production environment).

A12 pass criterion: Tel-model within statistical distance of the Tel-6000
specialist (ceiling) and clearly better than the mean(4000,8000) control.

Run: ~/tel_env/bin/python e2_eval.py --model phaseD/models/fe_tel_e2.model --device cuda
"""
import argparse
import numpy as np
from ase.io import read

ap = argparse.ArgumentParser()
ap.add_argument("--model", required=True)
ap.add_argument("--device", default="cuda")
a = ap.parse_args()

from mace.calculators import MACECalculator

# stamped test file carries elec_temp=6000 in each frame's info; the
# calculator reads it from there (info key overrides any argument).
frames = read("phaseD/test_Tel6000_stamped.xyz", index=":")
phases = [f.info.get("phase", "?") for f in frames]
eref = np.array([f.info["REF_energy"] / len(f) for f in frames])
fref = [f.arrays["REF_forces"] for f in frames]

calc = MACECalculator(model_paths=a.model, device=a.device, default_dtype="float64")
Et, Ft = [], []
for fr in frames:
    g = fr.copy(); g.calc = calc
    Ft.append(g.get_forces())
    Et.append(g.get_potential_energy() / len(g))

z = np.load("phaseD/e2_control.npz", allow_pickle=True)

def rmse_by_phase(fpred):
    out = {}
    for ph in sorted(set(phases)):
        ds = np.concatenate([fpred[i] - fref[i]
                             for i in range(len(frames)) if phases[i] == ph])
        out[ph] = 1000 * np.sqrt((ds ** 2).mean())
    ds = np.concatenate([fpred[i] - fref[i] for i in range(len(frames))])
    out["ALL"] = 1000 * np.sqrt((ds ** 2).mean())
    return out

rt = rmse_by_phase(Ft)
rs = rmse_by_phase(list(z["f6000"]))
rm = rmse_by_phase(list(z["f_mean"]))

print(f"{'model':22s} " + " ".join(f"{p:>10s}" for p in sorted(set(phases))) + f" {'ALL':>10s}")
for name, r in [("Tel-model @6000 (E2)", rt),
                ("specialist Tel6000", rs),
                ("mean(4000,8000)", rm)]:
    print(f"{name:22s} " + " ".join(f"{r[p]:10.1f}" for p in sorted(set(phases))) + f" {r['ALL']:10.1f}")

KB_MEV = 0.08617333  # meV/K
DS_KB = 0.84          # entropy of fusion in kB, median of the Phase E transient estimates

def energy_metrics(e_pred):
    d = (np.asarray(e_pred, float) - eref) * 1000.0     # meV/atom
    ph = np.array(phases)
    per = {p: float(d[ph == p].mean()) for p in sorted(set(phases))}
    sem = {p: float(d[ph == p].std(ddof=1) / np.sqrt((ph == p).sum())) for p in sorted(set(phases))}
    lh = per.get("liquid", np.nan) - per.get("hcp", np.nan)
    lh_err = float(np.hypot(sem.get("liquid", np.nan), sem.get("hcp", np.nan)))
    return dict(rmse=float(np.sqrt((d ** 2).mean())), offset=float(d.mean()),
                rmse_off=float(d.std()), per=per, lh=float(lh), lh_err=lh_err)

mt = energy_metrics(Et)
ms = energy_metrics(z["e6000"])
mm = energy_metrics(z["e_mean"])

print("\nENERGY (meV/atom). liq-hcp = liquid mean error minus hcp mean error;")
print("first-order Tm bias = (liq-hcp) / dS_m with dS_m = %.2f kB" % DS_KB)
print(f"{'model':22s} {'rmse':>6s} {'offset':>7s} {'rmse-off':>8s} {'liq-hcp':>13s} {'Tm bias':>8s}")
for name, m in [("Tel-model @6000 (E2)", mt), ("specialist Tel6000", ms), ("mean(4000,8000)", mm)]:
    print(f"{name:22s} {m['rmse']:6.2f} {m['offset']:+7.2f} {m['rmse_off']:8.2f} "
          f"{m['lh']:+6.2f}+-{m['lh_err']:.2f} {m['lh'] / (DS_KB * KB_MEV):+7.0f} K")

# A17 criteria (registered 2026-09-29, before any Tel-model result existed):
# the Tel model must sit closer to the specialist than to the control on both
# energy measures, i.e. below the specialist/control midpoints, and must not
# lose more than 10% in force accuracy relative to the specialist.
thr_off = 0.5 * (ms["rmse_off"] + mm["rmse_off"])
thr_lh = 0.5 * (abs(ms["lh"]) + abs(mm["lh"]))
thr_f = 1.10 * rs["ALL"]
c1 = mt["rmse_off"] < thr_off
c2 = abs(mt["lh"]) < thr_lh
c3 = rt["ALL"] < thr_f
print(f"\nA17 criteria:")
print(f"  (1) offset-removed energy RMSE {mt['rmse_off']:.2f} < {thr_off:.2f} meV/atom : {'pass' if c1 else 'FAIL'}")
print(f"  (2) |liq-hcp| {abs(mt['lh']):.2f} < {thr_lh:.2f} meV/atom                : {'pass' if c2 else 'FAIL'}")
print(f"  (3) force RMSE {rt['ALL']:.1f} < {thr_f:.1f} meV/A                    : {'pass' if c3 else 'FAIL'}")
print(f"E2 VERDICT: {'PASS' if (c1 and c2 and c3) else 'FAIL'}")
