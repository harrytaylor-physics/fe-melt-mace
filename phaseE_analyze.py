#!/usr/bin/env python3
"""
phaseE_analyze.py — aggregate Phase E runs into melting points and the curve
(fe_melt_v1). Reads phaseE/runs/*.json.

Per (branch, nominal pressure):
  NVE coexistence: COEXISTING runs give (P*, T*) points directly (mean over
    seeds; errors combined in quadrature with seed scatter).
  Pinning: linear fit of mean bias force vs T across the scanned temperatures;
    Tm = zero crossing, error from fit covariance; P from run average.
  G-E2: |Tm_NVE - Tm_pin| < 100 K after mapping both to a common pressure
    with the local Simon-fit slope.
Outputs: phaseE/melting_points.json, phaseE/summary.txt, Simon fit of the
NVE points, and anchor comparison at 323/330 GPa.
"""
import glob, json
import numpy as np

def simon(P, T0, P0, a, c): return T0 * ((P - P0) / a + 1) ** c

runs = [json.load(open(f)) for f in glob.glob("phaseE/runs/*.json")]
coex = {}; pins = {}
for r in runs:
    key = r["tag"]
    if r.get("method") == "pinning":
        pins.setdefault(key, []).append(r)
    elif r.get("verdict") == "COEXISTING":
        coex.setdefault(key, []).append(r)

points = []
lines = []
for key, rs in sorted(coex.items()):
    T = np.array([r["T"] for r in rs]); P = np.array([r["P"] for r in rs])
    Te = np.sqrt(np.mean([r["T_err"]**2 for r in rs]) + (T.std(ddof=1)**2 if len(rs)>1 else 0)/max(len(rs),1))
    pt = dict(tag=key, method="nve", P=float(P.mean()), T=float(T.mean()),
              T_err=float(Te), n_runs=len(rs))
    points.append(pt)
    lines.append(f"{key:24s} NVE     P={pt['P']:6.1f} GPa  Tm={pt['T']:6.0f} +/- {pt['T_err']:.0f} K  ({len(rs)} runs)")

for key, rs in sorted(pins.items()):
    rs = sorted(rs, key=lambda r: r["T"])
    T = np.array([r["T"] for r in rs]); F = np.array([r["force"] for r in rs])
    W = 1/np.array([max(r["force_err"],1e-12) for r in rs])**2
    if len(rs) < 3 or F.min() > 0 or F.max() < 0:
        lines.append(f"{key:24s} PIN     insufficient bracket ({len(rs)} temps)"); continue
    A = np.vstack([T, np.ones_like(T)]).T
    cov = np.linalg.inv(A.T @ (A * W[:,None]))
    m, b = cov @ A.T @ (W * F)
    Tm = -b/m
    Tm_err = float(np.sqrt(cov[0,0]*(b/m**2)**2 + cov[1,1]/m**2))
    Pm = float(np.mean([r["P"] for r in rs]))
    points.append(dict(tag=key, method="pinning", P=Pm, T=float(Tm), T_err=Tm_err, n_runs=len(rs)))
    lines.append(f"{key:24s} PIN     P={Pm:6.1f} GPa  Tm={Tm:6.0f} +/- {Tm_err:.0f} K  ({len(rs)} temps)")

# G-E2 cross-method agreement (map pin point to NVE pressure with dTm/dP ~ 8.8 K/GPa)
for key in sorted(set(p["tag"] for p in points)):
    got = {p["method"]: p for p in points if p["tag"] == key}
    if "nve" in got and "pinning" in got:
        dP = got["nve"]["P"] - got["pinning"]["P"]
        Tpin = got["pinning"]["T"] + 8.8 * dP
        diff = abs(got["nve"]["T"] - Tpin)
        ok = diff < 100
        lines.append(f"{key:24s} G-E2    |dTm|={diff:.0f} K -> {'pass' if ok else 'FAIL'}")

nve = [p for p in points if p["method"] == "nve"]
if len(nve) >= 4:
    from scipy.optimize import curve_fit
    P = np.array([p["P"] for p in nve]); T = np.array([p["T"] for p in nve])
    s = np.array([max(p["T_err"], 20) for p in nve])
    try:
        popt, pcov = curve_fit(simon, P, T, p0=[3700, 100, 30, 1.8], sigma=s, absolute_sigma=True)
        lines.append(f"\nSimon fit (ref P0={popt[1]:.0f}): T0={popt[0]:.0f} K a={popt[2]:.1f} GPa c={popt[3]:.2f}")
        for Pq, name, val, err in [(323, "Sun 2023", 6357, 45), (330, "consensus", 6370, 200)]:
            lines.append(f"  Tm({Pq} GPa) = {simon(Pq,*popt):.0f} K   [{name}: {val} +/- {err}]")
    except Exception as e:
        lines.append(f"Simon fit failed: {e}")

json.dump(points, open("phaseE/melting_points.json", "w"), indent=1)
open("phaseE/summary.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
