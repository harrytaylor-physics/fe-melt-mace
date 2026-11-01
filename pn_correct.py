#!/usr/bin/env python3
"""
pn_correct.py — A29c rule applied: assign each coexistence point the pressure normal
to the interface (P_n, the liquid pressure) instead of the cell mean -tr(sigma)/3.

Input: phaseE/decomp_tel_tel.json (per-frame diagonal virial stress from decomp_tel.py,
stage tel, Tel model at the T_el used; kinetic part is isotropic and cancels in P_n - P_mean),
phaseE/analysis_live_all.json (Tel-6000) and phaseE/analysis_live_all_runs_tel.json (Tel model).
Per cell and per campaign: dP = mean over that campaign's final frames of (P_n - P_mean).
For the Tel-6000 frames the stress is evaluated with the Tel model (proxy; the deviatoric
stress is set by the cell geometry and elastic response, which both potentials share closely).
Each cell's pressure is moved by dP, both Simon curves are refitted, and the matched-pressure
difference and the literature comparison are reprinted. No torch; runs on a login node.
  python3.13 pn_correct.py
"""
import json
import numpy as np
from scipy.optimize import curve_fit

P_REF = 100.0
LIT = [("Sun2023BccHcp", 323, 6357, 45), ("Alfe2009Coexistence", 328, 6370, 100),
       ("Wu2024DeepPotential", 330, 6415, None), ("GonzalezCataldo2023", 330, 6523, 8),
       ("Rao2025DMFTMelting", 330, 6225, 42)]


def simon(P, T1, a, c):
    return T1 * (1.0 + (np.asarray(P, float) - P_REF) / a) ** (1.0 / c)


def fit(cells):
    h = sorted((v["P"], v["Tm"], max(v.get("sem", 10.0), 10.0)) for k, v in cells.items() if "_hcp_" in k)
    P, T, s = map(np.array, zip(*h))
    p, cov = curve_fit(simon, P, T, p0=[T[0], 120.0, 2.0], sigma=s, absolute_sigma=True, maxfev=20000)
    chi2 = float((((T - simon(P, *p)) / s) ** 2).sum()) / (len(P) - 3)
    return p, np.sqrt(np.diag(cov)), chi2


def main():
    fr = json.load(open("phaseE/decomp_tel_tel.json"))
    dP = {}
    for r in fr:
        if "p_axes" not in r:
            continue
        d = r["p_axes"][r["normal"]] - float(np.mean(r["p_axes"]))
        dP.setdefault((r["src"], r["cell"]), []).append(d)
    camp = {"runs": json.load(open("phaseE/analysis_live_all.json")),
            "runs_tel": json.load(open("phaseE/analysis_live_all_runs_tel.json"))}
    name = {"runs": "Tel-6000", "runs_tel": "Tel model"}

    print("P_n - P_mean per cell (GPa), mean +- sd over final frames:")
    print(f"{'cell':16s} {'Tel-6000':>14s} {'Tel model':>14s}")
    cells = sorted({c for (_, c) in dP}, key=lambda c: int(c.split("_")[2][:-3]))
    for c in cells:
        row = []
        for src in ("runs", "runs_tel"):
            v = np.array(dP.get((src, c), [np.nan]))
            row.append(f"{v.mean():+6.2f} +- {v.std():4.2f}")
        print(f"{c:16s} {row[0]:>14s} {row[1]:>14s}")

    new, fits = {}, {}
    for src, a in camp.items():
        cc = {}
        for c, v in a["cells"].items():
            if "_hcp_" not in c:
                continue
            d = float(np.mean(dP.get((src, c), [0.0])))
            cc[c] = dict(v, P_mean=v["P"], P=v["P"] + d, dP=d)
        new[src] = cc
        p, e, chi2 = fit(cc)
        fits[src] = p
        p0, e0, chi20 = fit({k: v for k, v in a["cells"].items() if "_hcp_" in k})
        print(f"\n{name[src]} Simon fit, cell-mean pressure: T1 = {p0[0]:.0f}+-{e0[0]:.0f} K, a = {p0[1]:.0f}+-{e0[1]:.0f} GPa, "
              f"c = {p0[2]:.2f}+-{e0[2]:.2f} (chi2/dof {chi20:.2f})")
        print(f"{name[src]} Simon fit, P_n:                T1 = {p[0]:.0f}+-{e[0]:.0f} K, a = {p[1]:.0f}+-{e[1]:.0f} GPa, "
              f"c = {p[2]:.2f}+-{e[2]:.2f} (chi2/dof {chi2:.2f})")
        print(f"  shift of the curve at fixed P: " + ", ".join(
            f"{P:.0f} GPa {simon(P, *p) - simon(P, *p0):+.0f} K" for P in (100, 200, 330, 360)))

    print("\nTel model minus Tel-6000 at matched pressure, P_n convention:")
    for c, v in new["runs_tel"].items():
        print(f"  {c:16s} P_n {v['P']:6.1f} GPa   Tm {v['Tm']:5.0f} K   dTm@P {v['Tm'] - simon(v['P'], *fits['runs']):+5.0f} K")

    print("\nLiterature (this work minus reference), P_n convention:")
    print(f"{'reference':22s} {'P':>4s} {'Tm':>10s} {'Tel-6000':>9s} {'Tel model':>9s}")
    for k, P, T, e in LIT:
        print(f"{k:22s} {P:4d} {T:6d}{('+-%d' % e) if e else '':4s} {simon(P, *fits['runs']) - T:+9.0f} "
              f"{simon(P, *fits['runs_tel']) - T:+9.0f}")
    json.dump({src: dict(cells=new[src], fit=dict(T1=float(fits[src][0]), a=float(fits[src][1]), c=float(fits[src][2])))
               for src in new}, open("phaseE/pn_correct.json", "w"), indent=1)
    print("\nwritten phaseE/pn_correct.json")


if __name__ == "__main__":
    main()
