#!/usr/bin/env python3
"""
phaseE_analyze_live.py — running analysis of the Phase E coexistence campaign
(fe_melt_v1). Works on partial runs: reads every checkpoint in phaseE/runs/,
so it can be rerun at any time as segments land.

WHAT IT DOES
  1. Per run: classify as MELTED / FROZE (terminated brackets), EQUILIBRATING
     (too short), DRIFTING (solid fraction or T still trending), or CONVERGED
     (on the melting line). For converged runs, (T, P) are averaged over the
     second half of the available trajectory with a block-averaged error.
  2. Per cell: Tm and P from the converged runs, run-to-run scatter, standard
     error, and seed-pair agreement.
  3. Curve: Clapeyron slopes between adjacent pressures, a Simon fit
     Tm = T1 (1 + (P - 100)/a)^(1/c) once >= 3 hcp pressures have converged,
     and a comparison with published values at the inner-core boundary.

PHYSICAL REALISM CHECKS (converged is necessary, not sufficient)
  a. Temperature fluctuations. For N atoms the canonical kinetic-temperature
     fluctuation is sigma_T = T sqrt(2/(3N)) (~0.7% at N = 12,096); NVE values
     are somewhat smaller. Measured detrended sigma_T more than 3x this means
     the dynamics are unstable (integration or force errors), not thermal.
  b. Entropy of fusion from the transient. In NVE, total energy is fixed, so
     melting a fraction dx of the cell cools it by dT = (L / c) dx. With the
     classical heat capacity c = 3 kB per atom (the MD is classical; electronic
     heat capacity enters only through the Tel-fixed free-energy surface), the
     slope dT/d(solid fraction) gives dS_m / kB = 3 (dT/dx) / T. Metals melt
     with dS_m of order 1 kB (Richards' rule); values outside 0.5-2 kB are
     flagged. Caveat: the q6 solid fraction counts interface atoms imperfectly,
     so this is an order-of-magnitude check, not a measurement.
  c. Bracket ordering. Within a cell, runs started hotter (larger ti) should
     never FREEZE while cooler ones coexist, and vice versa for MELTED.
  d. Coexistence geometry. Mean solid fraction outside 0.15-0.85 means one
     slab is thin enough for the two interfaces to interact.
  e. Pressure. Coexistence P far (>15%) from the cell's nominal P means the
     cell was built at a lattice inconsistent with the potential.
  f. Curve shape. Tm must rise with P, and dTm/dP should fall with P (the
     Simon form's concavity). A violation points to a bad cell, not physics.
  g. Seed agreement. The two seeds at the same (cell, ti) should agree
     within ~80 K once converged.
  h. (--energy-check, GPU only) NVE energy conservation: from one converged
     checkpoint per cell, run a short NVE continuation with the production
     settings and fit the total-energy drift. Tolerance: in coexistence, drift
     changes the solid fraction rather than T. Over 200 ps a drift dE moves the
     fraction by dE/L with L ~ Tm dS_m ~ 0.34 eV/atom at 4,000 K; keeping that
     below 0.1 requires drift < ~0.17 meV/atom/ps.

USAGE (login node, no torch needed):
  python3 phaseE_analyze_live.py            # report + phaseE/analysis_live.json
  python3 phaseE_analyze_live.py --plot     # also phaseE/melt_curve_live.png
Energy check (GPU node, production env):
  python3.13 phaseE_analyze_live.py --energy-check --device cuda
"""
import argparse, glob, json, os, re
import numpy as np

TAG = re.compile(r"coex_(hcp|fcc)_(\d+)GPa_ti([\d.]+)_s(\d+)$")
P_REF = 100.0  # GPa, Simon reference pressure
SUFFIX = ""

# Published melting temperatures near the inner-core boundary. Values taken
# from sources retrieved during this project; verify each against the paper
# before publication.
LIT = [
    ("Alfe 2009, DFT coexistence", 328.0, 6370.0, 100.0),
    ("Wu et al. 2024, DeepMD", 330.0, 6415.0, None),
    ("DFT+DMFT MLIP, arXiv 2512.25061", 330.0, 6225.0, None),
    ("Sun et al. 2023", 323.0, 6357.0, 45.0),
]

# thresholds
MIN_PS = 40          # ps before a run can be judged
MIN_WIN = 20         # ps minimum averaging window
FR_TOL = 0.04        # max solid-fraction change across the window
T_TOL = 40.0         # K, max linear T change across the window
SCATTER_WARN = 100.0 # K, run-to-run scatter within a cell
SEED_WARN = 80.0     # K, seed-pair disagreement
FLUCT_WARN = 3.0     # measured / canonical sigma_T
DS_RANGE = (0.5, 2.0)
FR_RANGE = (0.15, 0.85)
P_DEV_WARN = 0.15


EV_A3_GPA = 160.21766208      # 1 eV/A^3 in GPa
KB_EV = 8.617333262e-5


def pressure_total(z):
    """Total pressure series. The driver stored -tr(sigma)/3 from ASE's
    get_stress(), whose default excludes the kinetic (ideal-gas) term, so the
    stored series is the virial pressure only. The kinetic term is
    N kB T / V with the same kinetic temperature the driver recorded each ps,
    and the NVE cell volume is fixed, so the correction is exact."""
    P = np.asarray(z["P_t"], float)
    if "cell" not in z.files or "positions" not in z.files:
        return P
    V = abs(np.linalg.det(np.asarray(z["cell"], float)))
    N = len(z["positions"])
    T = np.asarray(z["T_t"], float)
    return P + N * KB_EV * T / V * EV_A3_GPA


def load_runs(rundir):
    runs = []
    for f in sorted(glob.glob(os.path.join(rundir, "coex_*_ckpt.npz"))):
        tag = os.path.basename(f)[:-9]
        m = TAG.match(tag)
        if not m:
            continue
        z = np.load(f, allow_pickle=True)
        fr = np.asarray(z["fr_t"], float).reshape(-1, 2)
        st = verdict = ""
        js = f[:-9] + ".json"
        if os.path.exists(js):
            d = json.load(open(js))
            st, verdict = d.get("status", ""), d.get("verdict", "") or ""
        runs.append(dict(
            tag=tag, ckpt=f, branch=m.group(1), Pnom=int(m.group(2)),
            ti=float(m.group(3)), seed=int(m.group(4)),
            cell=f"coex_{m.group(1)}_{m.group(2)}GPa",
            T=np.asarray(z["T_t"], float), P=pressure_total(z),
            fr=fr, ps=int(z["blk"]), status=st, verdict=verdict))
    return runs


def cell_meta(celldir, cell):
    meta, N = {}, None
    p = os.path.join(celldir, cell + ".json")
    if os.path.exists(p):
        meta = json.load(open(p))
    x = os.path.join(celldir, cell + ".extxyz")
    if os.path.exists(x):
        with open(x) as fh:
            N = int(fh.readline().split()[0])
    return meta, N


def block_sem(x, nb=4):
    b = [c.mean() for c in np.array_split(x, nb) if len(c)]
    return float(np.std(b, ddof=1) / np.sqrt(len(b))) if len(b) > 1 else float("nan")


def entropy_of_fusion(r):
    """dS_m/kB from the NVE transient: 3 (dT/dx) / T, c = 3 kB classical."""
    if len(r["fr"]) < 4:
        return None
    b = r["fr"][:, 0].astype(int)
    f = r["fr"][:, 1]
    ok = b < len(r["T"])
    b, f = b[ok], f[ok]
    if len(f) < 4 or f.max() - f.min() < 0.10:
        return None
    Tb = r["T"][b]
    slope = np.polyfit(f, Tb, 1)[0]
    return float(3.0 * slope / Tb.mean())


def analyse_run(r, N):
    if r["verdict"] in ("MELTED", "FROZE"):
        r["state"] = r["verdict"]
        return
    n = len(r["T"])
    if n < MIN_PS or len(r["fr"]) < 3:
        r["state"] = "EQUILIBRATING"
        return
    w0 = min(n // 2, n - MIN_WIN)
    Tw, Pw = r["T"][w0:], r["P"][w0:]
    frw = r["fr"][r["fr"][:, 0] >= w0]
    if len(frw) < 3:
        r["state"] = "EQUILIBRATING"
        return
    dfr = np.polyfit(frw[:, 0], frw[:, 1], 1)[0] * (frw[-1, 0] - frw[0, 0])
    tt = np.arange(len(Tw))
    pT = np.polyfit(tt, Tw, 1)
    dT = pT[0] * len(Tw)
    resid = Tw - np.polyval(pT, tt)
    r.update(win=(w0, n), Tm=float(Tw.mean()), Tsem=block_sem(Tw),
             Pm=float(Pw.mean()), Psem=block_sem(Pw),
             frm=float(frw[:, 1].mean()), dfr=float(dfr), dT=float(dT),
             sigT=float(resid.std(ddof=1)))
    r["sigT_can"] = r["Tm"] * np.sqrt(2.0 / (3.0 * N)) if N else float("nan")
    r["state"] = "CONVERGED" if (abs(dfr) <= FR_TOL and abs(dT) <= T_TOL) else "DRIFTING"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rundir", default="phaseE/runs")
    ap.add_argument("--celldir", default="phaseE/cells")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--include-drifting", action="store_true",
                    help="average DRIFTING runs with CONVERGED ones (justified once the "
                         "energy check shows the slow fraction change is not integration drift)")
    ap.add_argument("--energy-check", action="store_true")
    ap.add_argument("--echeck-cells", default="", help="comma-separated cells for --energy-check (default all)")
    ap.add_argument("--model", default="phaseD/models/fe_ft_Tel6000_f32.model")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--steps", type=int, default=1000, help="energy-check NVE steps (2 fs); 1000 steps = 2 ps, ~20 min per cell at 12k atoms")
    a = ap.parse_args()

    global SUFFIX
    SUFFIX = "_all" if a.include_drifting else ""
    if a.rundir.rstrip("/") != "phaseE/runs":
        SUFFIX += "_" + os.path.basename(a.rundir.rstrip("/"))
    runs = load_runs(a.rundir)
    if not runs:
        raise SystemExit("no checkpoints found")
    cells = sorted({r["cell"] for r in runs},
                   key=lambda c: (c.split("_")[1] != "hcp", int(c.split("_")[2][:-3])))
    meta = {c: cell_meta(a.celldir, c) for c in cells}
    for r in runs:
        analyse_run(r, meta[r["cell"]][1])
        r["dS"] = entropy_of_fusion(r)

    # ---------- progress ----------
    states = {}
    for r in runs:
        states[r["state"]] = states.get(r["state"], 0) + 1
    print(f"runs: {len(runs)}   ps banked: {sum(r['ps'] for r in runs)}   " +
          "  ".join(f"{k}: {v}" for k, v in sorted(states.items())))

    # ---------- per run ----------
    print(f"\n{'run':32s} {'ps':>4s} {'state':>13s} {'T/K':>7s} {'+-':>4s} {'P/GPa':>7s} "
          f"{'solid':>6s} {'sigT/can':>8s} {'dS/kB':>6s}")
    for c in cells:
        for r in sorted([x for x in runs if x["cell"] == c], key=lambda x: (x["ti"], x["seed"])):
            if r["state"] in ("CONVERGED", "DRIFTING"):
                ratio = r["sigT"] / r["sigT_can"] if r["sigT_can"] == r["sigT_can"] else float("nan")
                print(f"{r['tag']:32s} {r['ps']:4d} {r['state']:>13s} {r['Tm']:7.0f} {r['Tsem']:4.0f} "
                      f"{r['Pm']:7.1f} {r['frm']:6.2f} {ratio:8.2f} "
                      f"{(r['dS'] if r['dS'] is not None else float('nan')):6.2f}")
            else:
                print(f"{r['tag']:32s} {r['ps']:4d} {r['state']:>13s}")

    # ---------- per cell + checks ----------
    summary, flags = {}, []
    print(f"\n{'cell':18s} {'n':>2s} {'Tm/K':>7s} {'scat':>5s} {'sem':>4s} {'P/GPa':>7s} "
          f"{'Tm_est':>7s} {'dev%':>5s}")
    for c in cells:
        cr = [r for r in runs if r["cell"] == c]
        conv = [r for r in cr if r["state"] == "CONVERGED"
                or (a.include_drifting and r["state"] == "DRIFTING")]
        Tm_est = meta[c][0].get("Tm_est")
        Pnom = cr[0]["Pnom"]
        if not conv:
            print(f"{c:18s}  0   (no converged runs yet)")
            continue
        Ts = np.array([r["Tm"] for r in conv])
        Ps = np.array([r["Pm"] for r in conv])
        scat = float(Ts.std(ddof=1)) if len(Ts) > 1 else float("nan")
        rsem = float(np.sqrt(np.nanmean([r["Tsem"] ** 2 for r in conv])))
        sem = float(max(scat / np.sqrt(len(Ts)) if len(Ts) > 1 else 0.0, rsem / np.sqrt(len(Ts))))
        Tm, Pm = float(Ts.mean()), float(Ps.mean())
        dev = 100 * (Tm - Tm_est) / Tm_est if Tm_est else float("nan")
        summary[c] = dict(branch=cr[0]["branch"], Pnom=Pnom, n=len(conv), Tm=Tm, scatter=scat,
                          sem=sem, P=Pm, Psem=float(Ps.std(ddof=1) / np.sqrt(len(Ps))) if len(Ps) > 1 else float("nan"),
                          Tm_est=Tm_est)
        print(f"{c:18s} {len(conv):2d} {Tm:7.0f} {scat:5.0f} {sem:4.0f} {Pm:7.1f} "
              f"{(Tm_est or float('nan')):7.0f} {dev:5.1f}")

        # checks
        if len(Ts) > 1 and scat > SCATTER_WARN:
            flags.append(f"{c}: run-to-run scatter {scat:.0f} K > {SCATTER_WARN:.0f} K")
        if abs(Pm - Pnom) / Pnom > P_DEV_WARN:
            flags.append(f"{c}: coexistence P {Pm:.1f} GPa is >{100*P_DEV_WARN:.0f}% from nominal {Pnom}")
        for r in conv:
            ratio = r["sigT"] / r["sigT_can"] if r["sigT_can"] == r["sigT_can"] else 0
            if ratio > FLUCT_WARN:
                flags.append(f"{r['tag']}: T fluctuation {ratio:.1f}x canonical (instability?)")
            if not (FR_RANGE[0] <= r["frm"] <= FR_RANGE[1]):
                flags.append(f"{r['tag']}: mean solid fraction {r['frm']:.2f} outside {FR_RANGE} (interfaces close)")
        for ti in sorted({r["ti"] for r in conv}):
            pair = [r for r in conv if r["ti"] == ti]
            if len(pair) == 2 and abs(pair[0]["Tm"] - pair[1]["Tm"]) > SEED_WARN:
                flags.append(f"{c} ti{ti:.2f}: seeds disagree by {abs(pair[0]['Tm']-pair[1]['Tm']):.0f} K")
        conv_ti = [r["ti"] for r in conv]
        for r in cr:
            if r["state"] == "FROZE" and r["ti"] > min(conv_ti):
                flags.append(f"{r['tag']}: FROZE although a cooler bracket coexists (ordering violated)")
            if r["state"] == "MELTED" and r["ti"] < max(conv_ti):
                flags.append(f"{r['tag']}: MELTED although a hotter bracket coexists (ordering violated)")
        dS = [r["dS"] for r in cr if r["dS"] is not None]
        if dS:
            med = float(np.median(dS))
            summary[c]["dS_kB"] = med
            if not (DS_RANGE[0] <= med <= DS_RANGE[1]):
                flags.append(f"{c}: transient-derived dS_m = {med:.2f} kB outside {DS_RANGE} (check)")

    # ---------- curve ----------
    hcp = sorted([(v["P"], v["Tm"], v["sem"], k) for k, v in summary.items() if v["branch"] == "hcp"])
    fit = None
    if len(hcp) >= 2:
        print("\nClapeyron slopes (hcp):")
        slopes = []
        for (P1, T1, _, _), (P2, T2, _, _) in zip(hcp, hcp[1:]):
            s = (T2 - T1) / (P2 - P1)
            slopes.append(s)
            print(f"  {P1:6.1f} -> {P2:6.1f} GPa   dTm/dP = {s:6.2f} K/GPa")
            if s <= 0:
                flags.append(f"Tm does not increase between {P1:.0f} and {P2:.0f} GPa")
        for s1, s2, (_, _, _, k) in zip(slopes, slopes[1:], hcp[1:]):
            if s2 > s1 * 1.25:
                flags.append(f"slope rises through {k} ({s1:.1f} -> {s2:.1f} K/GPa): curve not concave")
    if len(hcp) >= 3:
        try:
            from scipy.optimize import curve_fit
            P = np.array([h[0] for h in hcp]); T = np.array([h[1] for h in hcp])
            sig = np.array([max(h[2], 10.0) for h in hcp])

            def simon(P, T1, a, c):
                return T1 * (1.0 + (P - P_REF) / a) ** (1.0 / c)
            popt, pcov = curve_fit(simon, P, T, p0=[T[0], 50.0, 3.0], sigma=sig,
                                   absolute_sigma=True, bounds=([1000, 1, 0.5], [10000, 3000, 20]),
                                   maxfev=20000)
            perr = np.sqrt(np.diag(pcov))
            dof = len(P) - 3
            chi2 = float((((T - simon(P, *popt)) / sig) ** 2).sum())
            fit = dict(T1=float(popt[0]), a=float(popt[1]), c=float(popt[2]),
                       err=[float(x) for x in perr], dof=dof, chi2=chi2)
            print(f"\nSimon fit  Tm = T1 (1 + (P-{P_REF:.0f})/a)^(1/c):  "
                  f"T1 = {popt[0]:.0f}+-{perr[0]:.0f} K, a = {popt[1]:.0f}+-{perr[1]:.0f} GPa, "
                  f"c = {popt[2]:.2f}+-{perr[2]:.2f}   (dof {dof}"
                  + (f", chi2/dof {chi2/dof:.2f})" if dof > 0 else ", exact fit)"))
            print("Comparison near the inner-core boundary:")
            for name, Pl, Tl, el in LIT:
                Tf = simon(Pl, *popt)
                print(f"  {name:34s} {Pl:5.0f} GPa  {Tl:6.0f}{('+-%.0f' % el) if el else '':6s} K   "
                      f"this work {Tf:6.0f} K   diff {Tf - Tl:+5.0f} K")
        except Exception as e:
            print(f"\nSimon fit failed: {e}")

    # ---------- flags ----------
    print("\nPHYSICAL CHECKS:")
    print("  none flagged" if not flags else "\n".join("  WARN " + f for f in flags))

    json.dump(dict(cells=summary, fit=fit, flags=flags,
                   runs={r["tag"]: {k: (v if not isinstance(v, np.ndarray) else None)
                                    for k, v in r.items() if k not in ("T", "P", "fr", "ckpt")}
                         for r in runs}),
              open(f"phaseE/analysis_live{SUFFIX}.json", "w"), indent=1, default=float)
    print(f"\nwritten phaseE/analysis_live{SUFFIX}.json")

    if a.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(6.5, 4.5))
            for k, v in summary.items():
                mk, col = ("o", "k") if v["branch"] == "hcp" else ("s", "tab:green")
                ax.errorbar(v["P"], v["Tm"], yerr=v["sem"], fmt=mk, color=col, capsize=3,
                            label=None)
            if fit:
                Pp = np.linspace(min(h[0] for h in hcp) - 5, 360, 200)
                ax.plot(Pp, fit["T1"] * (1 + (Pp - P_REF) / fit["a"]) ** (1 / fit["c"]),
                        "k-", lw=1, label="this work, Simon fit")
            for name, Pl, Tl, el in LIT:
                ax.errorbar(Pl, Tl, yerr=el or 0, fmt="D", ms=5, capsize=3, label=name)
            ax.plot([], [], "ko", label="this work, hcp coexistence")
            if any(v["branch"] == "fcc" for v in summary.values()):
                ax.plot([], [], "s", color="tab:green", label="this work, fcc coexistence")
            ax.set_xlabel("P (GPa)"); ax.set_ylabel("T$_m$ (K)")
            ax.legend(fontsize=7, loc="lower right")
            fig.tight_layout(); fig.savefig(f"phaseE/melt_curve_live{SUFFIX}.png", dpi=200)
            print(f"written phaseE/melt_curve_live{SUFFIX}.png")
        except Exception as e:
            print(f"plot skipped: {e}")

    if a.energy_check:
        energy_check(runs, a)


def energy_check(runs, a):
    """Short NVE continuation from one converged checkpoint per cell; total
    energy drift with production settings (float32, cuEquivariance, skin 0.3)."""
    from phaseE_nlcache import enable_nl_cache
    enable_nl_cache(skin=0.3)
    from ase import units
    from ase.io import read
    from ase.md.verlet import VelocityVerlet
    from mace.calculators import MACECalculator
    calc = MACECalculator(model_paths=a.model, device=a.device,
                          default_dtype="float32", enable_cueq=True)
    print(f"\nENERGY CONSERVATION ({a.steps} steps x 2 fs; tolerance 0.17 meV/atom/ps)")
    only = [x for x in a.echeck_cells.split(",") if x]
    for c in sorted({r["cell"] for r in runs}):
        if only and c not in only:
            continue
        cand = [r for r in runs if r["cell"] == c and r["state"] == "CONVERGED"]
        if not cand:
            continue
        r = max(cand, key=lambda x: x["ps"])
        at = read(os.path.join(a.celldir, c + ".extxyz"))
        z = np.load(r["ckpt"], allow_pickle=True)
        at.set_cell(z["cell"]); at.set_positions(z["positions"]); at.set_momenta(z["momenta"])
        at.pbc = True
        # Tel-conditioned model: use the electronic temperature the run used
        jf = os.path.join(a.rundir, r["tag"] + ".json")
        et = json.load(open(jf)).get("elec_temp") if os.path.exists(jf) else None
        if et is None and a.rundir.rstrip("/") != "phaseE/runs" and os.path.exists("phaseE/analysis_live_all.json"):
            et = round(json.load(open("phaseE/analysis_live_all.json"))["cells"][c]["Tm"])
        if et is not None:
            at.info["elec_temp"] = float(et)
        at.calc = calc
        dyn = VelocityVerlet(at, 2.0 * units.fs)
        E = []
        every = 50
        for _ in range(a.steps // every + 1):
            E.append((at.get_potential_energy() + at.get_kinetic_energy()) / len(at))
            dyn.run(every)
        E = 1000.0 * np.array(E)
        t = np.arange(len(E)) * every * 0.002
        p = np.polyfit(t, E, 1)
        fl = float((E - np.polyval(p, t)).std())
        verdict = "ok" if abs(p[0]) < 0.17 else "DRIFT"
        print(f"  {r['tag']:32s} T_el {at.info.get('elec_temp', float('nan')):5.0f}  drift {p[0]:+7.3f} meV/atom/ps   fluct {fl:6.3f} meV/atom   {verdict}")


if __name__ == "__main__":
    main()
