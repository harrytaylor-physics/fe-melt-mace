import sys, json, glob
sys.path.insert(0, '.')
from figstyle import *
import numpy as np
KB, EVA3 = 8.617333262e-5, 160.21766208
pn = json.load(open(DATA + "pn_correct.json"))
at = json.load(open(DATA + "analysis_live_all_runs_tel.json"))
fr = json.load(open(DATA + "decomp_tel_tel.json"))
CELLS = ["coex_hcp_100GPa", "coex_hcp_150GPa", "coex_hcp_200GPa", "coex_hcp_250GPa", "coex_hcp_300GPa", "coex_hcp_330GPa", "coex_hcp_350GPa"]
nom = [100, 150, 200, 250, 300, 330, 350]
f = pn["runs_tel"]["fit"]
def simon(P, T1=f["T1"], a=f["a"], c=f["c"]): return T1 * (1 + (np.asarray(P) - 100) / a) ** (1 / c)
def dsimon(P, T1=f["T1"], a=f["a"], c=f["c"]): return T1 / (c * a) * (1 + (np.asarray(P) - 100) / a) ** (1 / c - 1)
out = {}

# kinetic term and anisotropy per cell
kin, dpn, dpn_sd = [], [], []
for c in CELLS:
    ks = []
    for ff in sorted(glob.glob(DATA + f"runs_tel/{c}_ti*_ckpt.npz")):
        z = np.load(ff, allow_pickle=True); V = abs(np.linalg.det(z["cell"])); N = len(z["positions"])
        T = z["T_t"][100:]; ks.append((N * KB * T / V * EVA3).mean())
    kin.append(np.mean(ks))
    d = [r["p_axes"][r["normal"]] - np.mean(r["p_axes"]) for r in fr if r["cell"] == c and r["src"] == "runs_tel"]
    dpn.append(np.mean(d)); dpn_sd.append(np.std(d))
kin, dpn, dpn_sd = map(np.array, (kin, dpn, dpn_sd))
Pn = np.array([pn["runs_tel"]["cells"][c]["P"] for c in CELLS])
shift_kin = -dsimon(Pn) * kin; shift_pn = -dsimon(Pn) * dpn
out["pressure"] = dict(kin=kin.round(2).tolist(), dpn=dpn.round(2).tolist(), dpn_sd=dpn_sd.round(2).tolist(),
                       shift_kin=shift_kin.round(0).tolist(), shift_pn=shift_pn.round(0).tolist())

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.3), gridspec_kw=dict(wspace=0.32))
x = np.arange(len(CELLS)); lab = [f"{n}" for n in nom]
a1.bar(x, kin, width=0.62, color=THEORY)
for i, v in enumerate(kin): a1.text(i, v + 0.25, f"{v:.1f}", ha="center", fontsize=11, color=INK)
a1.set_xticks(x, lab); a1.set_xlabel("Cell (nominal pressure, GPa)"); a1.set_ylabel("N k_B T / V  (GPa)")
a1.set_title("Kinetic pressure added to the virial", fontsize=13); a1.set_ylim(0, 15.5); a1.grid(axis="x", visible=False)
a2.bar(x, dpn, width=0.62, color=OURS, yerr=dpn_sd, ecolor=GREY, capsize=3)
for i, v in enumerate(dpn): a2.text(i, v + dpn_sd[i] + 0.15, f"{v:.1f}", ha="center", fontsize=11, color=INK)
a2.set_xticks(x, lab); a2.set_xlabel("Cell (nominal pressure, GPa)"); a2.set_ylabel("P_n − P_mean  (GPa)")
a2.set_title("Normal minus cell-mean pressure", fontsize=13); a2.set_ylim(0, 6.6); a2.grid(axis="x", visible=False)
save(fig, "f_pressure")

# melting curve
Tm = np.array([pn["runs_tel"]["cells"][c]["Tm"] for c in CELLS]); sem = np.array([pn["runs_tel"]["cells"][c]["sem"] for c in CELLS])
fcc = at["cells"]["coex_fcc_100GPa"]
fig, ax = plt.subplots(figsize=(8.6, 5.6))
PP = np.linspace(80, 380, 300)
ax.plot(PP, simon(PP), color=OURS, lw=2.2, label="Simon fit, hcp")
ax.errorbar(Pn, Tm, yerr=sem, fmt="o", ms=8, color=OURS, mec="white", mew=1.2, capsize=0, label="hcp–liquid coexistence (6 runs each)")
ax.plot([fcc["P"]], [fcc["Tm"]], "s", ms=8, mfc="white", mec=EXPT, mew=2, label="fcc–liquid (2 runs, cell-mean P)")
for P_, T_, n_ in zip(Pn, Tm, nom):
    ax.annotate(f"{T_:,.0f} K", (P_, T_), xytext=(8, -16), textcoords="offset points", fontsize=10.5, color=MUTED)
ax.axvspan(329, 364, color="#F1F2F4", zorder=0); ax.text(346.5, 6950, "inner core\n329–364 GPa", ha="center", va="top", fontsize=10.5, color=MUTED)
ax.set_xlabel("Pressure P_n (GPa)"); ax.set_ylabel("Melting temperature (K)"); ax.set_xlim(80, 380); ax.set_ylim(3500, 7000)
ax.legend(loc="upper left", fontsize=11)
ax.text(372, 3700, f"T_m = {f['T1']:,.0f} K (1 + (P − 100 GPa)/{f['a']:.0f} GPa)^(1/{f['c']:.2f})", ha="right", fontsize=11, color=INK)
save(fig, "f_curve")
out["curve"] = dict(Pn=Pn.round(1).tolist(), Tm=Tm.round(0).tolist(), sem=sem.round(1).tolist(),
                    T330=float(simon(330)), T323=float(simon(323)), T328=float(simon(328)), T364=float(simon(364)), T100=float(simon(100)),
                    slope100=float(dsimon(100)), slope200=float(dsimon(200)), slope330=float(dsimon(330)),
                    fcc=[fcc["P"], fcc["Tm"], fcc["sem"]], fcc_minus_hcp=float(fcc["Tm"] - simon(fcc["P"])),
                    dS=[at["cells"][c].get("dS_kB") for c in CELLS])

# comparison with literature
from refs import r
fig, (ax, axi) = plt.subplots(1, 2, figsize=(12.8, 5.6), gridspec_kw=dict(width_ratios=[1.3, 1], wspace=0.22))
PP = np.linspace(90, 375, 300)
LIT = {x["id"]: x for x in json.load(open(layout.LITERATURE))["melting_curve_figure"]}   # melt_curve/validation/refs/literature.json
pt = lambda k, i=0: LIT[k]["points"][i]
wc = LIT["Wu2024"]["fit_coeffs"]; rc = LIT["Rao2025"]["fit_coeffs"]
wu = lambda P: wc[0] + wc[1] * P + wc[2] * P ** 2
rz = lambda P: rc[0] * (1 + (P - rc[1]) / rc[2]) ** (1 / rc[3])
def draw(A, dx=0.0):
    A.plot(PP, simon(PP), color=OURS, lw=2.6, label="This work: PBE, MACE, T_el = T", zorder=5)
    A.plot(PP, wu(PP), color=THEORY, lw=1.6, ls="--", label="Wu et al. 2024, fit (PBE, DeepMD) " + r("Wu2024"))
    Pr = np.linspace(*LIT["Rao2025"]["fit_P_range"], 50); A.plot(Pr, rz(Pr), color=BEYOND, lw=1.8, ls="-.", label="Rao & Zhu 2025, fit (DFT+DMFT) " + r("Rao2025"))
    s = LIT["Sun2023"]["points"]
    A.errorbar([s[0][0], s[1][0]], [s[0][1], s[1][1]], yerr=[s[0][2], s[1][2]], fmt="D", ms=7, color=THEORY, capsize=3, label="Sun et al. 2023 (PBE, TI) " + r("Sun2023"))
    A.errorbar([pt("Alfe2009")[0] - dx], [pt("Alfe2009")[1]], yerr=[pt("Alfe2009")[2]], fmt="o", ms=7, color=THEORY, mfc="white", mew=1.8, capsize=3, label="Alfè 2009 (PW91, coexistence) " + r("Alfe2009"))
    A.errorbar([pt("Wu2024")[0] + dx], [pt("Wu2024")[1]], yerr=[pt("Wu2024")[2]], fmt="^", ms=8, color=THEORY, mfc="white", mew=1.8, capsize=3, label="Wu et al. 2024 (coexistence; ±200 K quoted) " + r("Wu2024"))
    A.errorbar([pt("Sola2009")[0] - 2 * dx], [pt("Sola2009")[1]], yerr=[pt("Sola2009")[2]], fmt="v", ms=8, color=BEYOND, capsize=3, label="Sola & Alfè 2009 (DMC) " + r("Sola2009"))
    A.errorbar([pt("Anzellini2013")[0] + 2 * dx], [pt("Anzellini2013")[1]], yerr=[pt("Anzellini2013")[2]], fmt="s", ms=7, color=EXPT, mfc="white", mew=1.8, capsize=3, label="Anzellini et al. 2013 (DAC, extrapolated) " + r("Anzellini2013"))
    A.errorbar([pt("Sinmyo2019")[0] + 3 * dx], [pt("Sinmyo2019")[1]], yerr=[pt("Sinmyo2019")[2]], fmt="P", ms=8, color=EXPT, capsize=3, label="Sinmyo et al. 2019 (DAC, extrapolated) " + r("Sinmyo2019"))
    A.errorbar([pt("Balugani2024")[0] - 3 * dx], [pt("Balugani2024")[1]], yerr=[pt("Balugani2024")[2]], fmt="X", ms=8, color=EXPT, capsize=3, label="Balugani et al. 2024 (shock + XAS, extrapolated) " + r("Balugani2024"))
draw(ax); draw(axi, dx=1.4)
ax.set_xlim(90, 375); ax.set_ylim(3400, 7500); ax.set_xlabel("Pressure (GPa)"); ax.set_ylabel("Melting temperature (K)")
h, l = ax.get_legend_handles_labels()
ax.legend(h, l, loc="upper left", fontsize=9.3)
axi.set_xlim(318, 365); axi.set_ylim(5200, 7400); axi.set_xlabel("Pressure (GPa)"); axi.set_title("Inner-core range (symbols at 330 GPa offset for clarity)", fontsize=11)
save(fig, "f_compare")
NAMES = [("Sun et al. 2023, PBE TI", "Sun2023", 0), ("Sun et al. 2023, PBE TI", "Sun2023", 1),
         ("Alfè 2009, PW91 coexistence", "Alfe2009", 0), ("Wu et al. 2024, PBE DeepMD coexistence", "Wu2024", 0),
         ("Rao & Zhu 2025, DFT+DMFT", "Rao2025", 0), ("Sola & Alfè 2009, DMC", "Sola2009", 0),
         ("Anzellini et al. 2013, DAC", "Anzellini2013", 0), ("Sinmyo et al. 2019, DAC", "Sinmyo2019", 0),
         ("Balugani et al. 2024, shock + XAS", "Balugani2024", 0)]
lit = [(n, *pt(k, i), r(k)) for n, k, i in NAMES]
out["lit"] = [dict(name=n, P=P, T=T, err=e, ref=rr, ours=float(simon(P)), diff=float(simon(P) - T), nsig=(float((simon(P) - T) / e) if e else None)) for n, P, T, e, rr in lit]
json.dump(out, open(FIG + "results.json", "w"), indent=1)
print(json.dumps(out["curve"], indent=0)); print(out["pressure"])
for r in out["lit"]: print(r)
