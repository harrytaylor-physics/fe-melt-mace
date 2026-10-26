"""Triple-point slide: hcp melting curve (fit, extrapolated below 104 GPa), the fcc point, planned state points."""
import sys, json
sys.path.insert(0, '.')
from figstyle import *
import numpy as np
pn = json.load(open(DATA + "pn_correct.json"))["runs_tel"]["fit"]
simon = lambda P: pn["T1"] * (1 + (np.asarray(P) - 100) / pn["a"]) ** (1 / pn["c"])
res = json.load(open(FIG + "results.json"))["curve"]
fig, ax = plt.subplots(figsize=(7.4, 4.9))
Pm = np.linspace(104.3, 160, 100); Pe = np.linspace(45, 104.3, 100)
ax.plot(Pm, simon(Pm), color=OURS, lw=2.2, label="hcp–liquid, Simon fit (this work)")
ax.plot(Pe, simon(Pe), color=OURS, lw=2.0, ls="--", label="hcp–liquid, fit extrapolated below 104 GPa")
ax.errorbar(res["Pn"][:2], res["Tm"][:2], yerr=res["sem"][:2], fmt="o", color=OURS, ms=7, capsize=3, label="hcp–liquid, measured")
P_f, T_f, e_f = res["fcc"]
ax.errorbar([P_f], [T_f], yerr=[e_f], fmt="s", color=EXPT, ms=8, capsize=3, mfc="white", mew=2, label="fcc–liquid, measured (cell-mean P)")
for P in (60, 80):
    ax.plot([P], [simon(P)], "o", ms=12, mfc="none", mec=INK, mew=1.4)
ax.plot([], [], "o", ms=12, mfc="none", mec=INK, mew=1.4, label="planned: hcp and fcc coexistence")
ax.annotate(f"fcc {T_f:,.0f} K at {P_f:.0f} GPa\n{simon(P_f) - T_f:.0f} K below the hcp curve", (P_f, T_f), xytext=(P_f + 12, T_f - 80), fontsize=10.5, va="top",
            arrowprops=dict(arrowstyle="-", color=GREY, lw=0.8))
ax.axvspan(45, P_f, color=GREY, alpha=0.08, lw=0)
ax.text(78, 3020, "triple point expected\nbelow ~91 GPa", fontsize=10.5, color=MUTED, va="bottom", ha="left")
ax.set_xlim(45, 160); ax.set_ylim(2900, 4900)
ax.set_xlabel("Pressure (GPa)"); ax.set_ylabel("Temperature (K)")
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.8), fontsize=9.8, framealpha=1, edgecolor="#D0D4D9")
fig.savefig(FIG + "f_triple.png"); plt.close(fig)
print({P: round(float(simon(P))) for P in (60, 80, 91, 104.3)}, P_f, T_f)
