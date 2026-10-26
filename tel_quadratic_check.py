#!/usr/bin/env python3
"""
tel_quadratic_check.py — DFT-only test of how energies depend on electronic
temperature. No machine learning, no GPU; login node.

Every frame was labelled at Tel = 4000, 6000 and 8000 K with identical atomic
positions. This matches the three labels of each frame by geometry and asks:
if you know the DFT answer at 4000 K and 8000 K only, how well can you
predict the DFT answer at 6000 K?

Two interpolation rules are compared:
  linear in Tel     : F6 = (F4 + F8) / 2
  quadratic in Tel  : F6 = F4 + (F8 - F4) * (6000^2 - 4000^2) / (8000^2 - 4000^2)
The quadratic rule is the leading term of the Sommerfeld expansion of the
electronic free energy, F_el = -(1/2) gamma Tel^2, valid when k Tel is small
compared with the d-band width.

If the quadratic rule reproduces the 6000 K energies to a few meV/atom with no
phase dependence, a Tel model with that form built in would interpolate
reliably from the existing labels. If it does not, interpolation needs labels
at more electronic temperatures.

Run: python3 tel_quadratic_check.py
"""
import glob, hashlib
import numpy as np
from ase.io import read

W_QUAD = (6000.0**2 - 4000.0**2) / (8000.0**2 - 4000.0**2)   # = 0.4167


def key(at):
    h = hashlib.sha1()
    h.update(np.round(at.get_positions(), 4).tobytes())
    h.update(np.round(np.array(at.cell), 4).tobytes())
    return h.hexdigest()


def load(tel):
    out = {}
    for f in sorted(glob.glob(f"phaseD/*_Tel{tel}.xyz")):
        for at in read(f, index=":"):
            out[key(at)] = (at.info.get("phase", "?"),
                            at.info["REF_energy"] / len(at),
                            np.array(at.arrays["REF_forces"]))
    return out


d4, d6, d8 = load(4000), load(6000), load(8000)
common = sorted(set(d4) & set(d6) & set(d8))
print(f"frames labelled at all three Tels: {len(common)}  "
      f"(4000: {len(d4)}, 6000: {len(d6)}, 8000: {len(d8)})")

ph = np.array([d6[k][0] for k in common])
E4 = np.array([d4[k][1] for k in common])
E6 = np.array([d6[k][1] for k in common])
E8 = np.array([d8[k][1] for k in common])

print("\nDFT energy change between Tels (meV/atom):")
print(f"  F(4000) - F(8000): mean {1000*np.mean(E4-E8):.0f}, range {1000*np.min(E4-E8):.0f} to {1000*np.max(E4-E8):.0f}")

rules = {"linear in Tel": 0.5, "quadratic in Tel": W_QUAD}
print(f"\n{'rule':18s} {'offset':>8s} {'rmse-off':>9s} {'liq-hcp':>14s}   per-phase mean error (meV/atom)")
for name, w in rules.items():
    pred = E4 + (E8 - E4) * w
    d = 1000 * (pred - E6)
    per = {p: d[ph == p].mean() for p in sorted(set(ph))}
    sem = {p: d[ph == p].std(ddof=1) / np.sqrt((ph == p).sum()) for p in sorted(set(ph))}
    lh = per.get("liquid", np.nan) - per.get("hcp", np.nan)
    lhe = np.hypot(sem.get("liquid", np.nan), sem.get("hcp", np.nan))
    print(f"{name:18s} {d.mean():+8.2f} {d.std():9.2f} {lh:+7.2f}+-{lhe:.2f}   " +
          "  ".join(f"{p} {per[p]:+.2f}" for p in sorted(per)))

print("\nforces, same rules (meV/A, RMSE against DFT at 6000 K):")
for name, w in rules.items():
    dd = np.concatenate([d4[k][2] + (d8[k][2] - d4[k][2]) * w - d6[k][2] for k in common])
    print(f"  {name:18s} {1000*np.sqrt((dd**2).mean()):6.1f}")

# best-fit power: which exponent n in F = a + b Tel^n makes the interpolation exact on average?
pred_err = []
for n in (1.0, 1.5, 2.0, 2.5, 3.0):
    w = (6000**n - 4000**n) / (8000**n - 4000**n)
    d = 1000 * (E4 + (E8 - E4) * w - E6)
    pred_err.append((n, d.mean(), d.std()))
print("\nexponent n in F = a + b*Tel^n:  n, offset, rmse-off (meV/atom)")
for n, m, s in pred_err:
    print(f"  n = {n:.1f}   {m:+7.2f}   {s:6.2f}")
