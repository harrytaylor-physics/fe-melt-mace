"""Phases of iron in P-T: our hcp melting curve, fcc point, Komabayashi fcc-hcp boundary, ambient transitions."""
import sys, json
sys.path.insert(0, '.')
from figstyle import *
from refs import r
import numpy as np
pn = json.load(open(DATA + "pn_correct.json"))["runs_tel"]["fit"]
simon = lambda P: pn["T1"] * (1 + (np.asarray(P) - 100) / pn["a"]) ** (1 / pn["c"])
res = json.load(open(FIG + "results.json"))["curve"]
fig, ax = plt.subplots(figsize=(10.6, 5.6))
# inner core
ax.axvspan(329, 364, color=GREY, alpha=0.12, lw=0)
ax.text(346.5, 1150, "inner\ncore", ha="center", fontsize=10.5, color=MUTED)
ax.axvline(136, color=GREY, lw=0.8, ls=":")
ax.text(138, 1150, "core–mantle\nboundary", fontsize=10, color=MUTED)
# our hcp melting curve
Pm = np.linspace(104.3, 371, 200); Pe = np.linspace(50, 104.3, 60)
ax.plot(Pm, simon(Pm), color=OURS, lw=2.4, label="hcp–liquid melting, this work (PBE)")
ax.plot(Pe, simon(Pe), color=OURS, lw=1.8, ls="--", label="this work, extrapolated below 104 GPa")
ax.errorbar(res["Pn"], res["Tm"], yerr=res["sem"], fmt="o", color=OURS, ms=5, capsize=2)
P_f, T_f, e_f = res["fcc"]
ax.plot([P_f], [T_f], "s", ms=8, mfc="white", mec=EXPT, mew=2, label="fcc–liquid, this work (91 GPa)")
# Komabayashi fcc-hcp boundary: measured 21-69 GPa, dP/dT = 0.0394 GPa/K, triple point 88 GPa / 2,800 K (extrapolated)
kb = lambda P: 2800 - (88 - np.asarray(P)) / 0.0394
Pk = np.linspace(21, 69, 50); Pk2 = np.linspace(69, 88, 20)
ax.plot(Pk, kb(Pk), color=THEORY, lw=2.0, label=f"fcc–hcp boundary, measured 21–69 GPa {r('Komabayashi2009')}")
ax.plot(Pk2, kb(Pk2), color=THEORY, lw=1.6, ls=":")
ax.plot([88], [2800], "D", color=THEORY, ms=7, label="fcc–hcp–liquid triple point, extrapolated (88 GPa, 2,800 K)")
# ambient
ax.plot([0, 0, 0], [1185, 1667, 1811], "k_", ms=14, mew=2, clip_on=False)
ax.text(6, 470, "Ticks at 0 GPa (1 bar):\nα (bcc) → γ (fcc) 1,185 K; γ → δ (bcc) 1,667 K;\nmelting 1,811 K", fontsize=10, va="center")
# field labels
ax.text(230, 6400, "liquid", fontsize=15, color=OURS, weight="bold")
ax.text(230, 3200, "ε (hcp)", fontsize=15, color=INK, weight="bold")
ax.text(22, 2650, "γ (fcc)", fontsize=13, color=EXPT, weight="bold")
ax.set_xlim(0, 380); ax.set_ylim(0, 7400)
ax.set_xlabel("Pressure (GPa)"); ax.set_ylabel("Temperature (K)")
ax.legend(loc="upper left", fontsize=9.5, framealpha=1, edgecolor="#D0D4D9")
fig.savefig(FIG + "f_phases.png"); plt.close(fig)
print(kb(21), kb(69))
