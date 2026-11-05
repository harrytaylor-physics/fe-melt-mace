import sys, json
sys.path.insert(0, '.')
from figstyle import *
from refs import r
import numpy as np
res = json.load(open(FIG + "results.json"))
# ---- timeline of Tm at the ICB ----
# year, T, err, kind, two-line label, ref, label anchor (x, y), ha
pts = [
 (1993.0, 4850, 200, "E", "Boehler 1993\nDAC, speckle", r("Boehler1993"), (1993.6, 5060), "left"),
 (1999.0, 6700, 600, "T", "Alfè et al. 1999\nDFT free energy", r("Alfe1999"), (1999.6, 6960), "left"),
 (2009.0, 6370, 100, "T", "Alfè 2009\nDFT coexistence, 328 GPa", r("Alfe2009"), (2008.4, 5990), "right"),
 (2009.3, 6900, 400, "B", "Sola & Alfè 2009\nquantum Monte Carlo", r("Sola2009"), (2008.7, 7300), "right"),
 (2013.0, 6230, 500, "E", "Anzellini et al. 2013\nlaser-heated DAC + XRD", r("Anzellini2013"), (2013.6, 7300), "left"),
 (2019.0, 5500, 220, "E", "Sinmyo et al. 2019\nresistance-heated DAC", r("Sinmyo2019"), (2018.4, 5500), "right"),
 (2022.0, 6230, 540, "E", "Kraus et al. 2022\nlaser ramp + XRD (interpolated)", r("Kraus2022"), (2022.6, 4980), "left"),
 (2023.0, 6357, 45, "T", "Sun et al. 2023\nDFT TI, 323 GPa", r("Sun2023"), (2021.6, 7050), "right"),
 (2024.0, 6420, 200, "T", "Wu et al. 2024\nMLIP coexistence", r("Wu2024"), (2024.5, 6950), "left"),
 (2024.55, 6202, 514, "E", "Balugani et al. 2024\nlaser shock + XAS", r("Balugani2024"), (2024.9, 5460), "left"),
 (2025.9, 6225, 42, "B", "Rao & Zhu 2025\nDFT+DMFT MLIP", r("Rao2025"), (2026.4, 5990), "left"),
]
col = {"E": EXPT, "T": THEORY, "B": BEYOND}; mk = {"E": "s", "T": "o", "B": "D"}
fig, ax = plt.subplots(figsize=(11.2, 5.6))
for y, T, e, k, lab, rf, (lx, ly), ha in pts:
    ax.errorbar([y], [T], yerr=[e], fmt=mk[k], ms=8, color=col[k], capsize=3, mew=1.5)
    ax.annotate(f"{lab} {rf}", (y, T), xytext=(lx, ly), textcoords="data", fontsize=9.8, color=INK, ha=ha, va="center",
                linespacing=1.15, arrowprops=dict(arrowstyle="-", color="#9AA0A6", lw=0.8, shrinkA=2, shrinkB=5))
T330 = res["curve"]["T330"]
T330_se = json.load(open(os.path.join(layout.VAL_OUT, "melting_point.json")))["at"]["330"]["Tm_fit_se"]  # Simon-fit s.e.
ax.errorbar([2026.9], [T330], yerr=[T330_se], fmt="*", ms=16, color=OURS, capsize=3)
ax.annotate(f"This work 2026\n{T330:,.0f} ± {T330_se:.0f} K", (2026.9, T330), xytext=(2028.0, 6560), textcoords="data", fontsize=11,
            color=OURS, ha="left", va="center", weight="bold", arrowprops=dict(arrowstyle="-", color="#9AA0A6", lw=0.8, shrinkA=2, shrinkB=7))
for k, lab in (("E", "experiment (at or extrapolated to 330 GPa)"), ("T", "DFT (PBE or PW91)"), ("B", "beyond DFT (QMC, DMFT)")):
    ax.plot([], [], mk[k], color=col[k], label=lab, ms=8)
ax.plot([], [], "*", color=OURS, ms=13, label="this work (PBE, MACE)")
ax.legend(loc="lower left", bbox_to_anchor=(0.25, 0.0), fontsize=10.5, framealpha=1, edgecolor="#D0D4D9")
ax.set_xlim(1990, 2037); ax.set_ylim(4400, 7600); ax.set_xticks(range(1990, 2030, 5))
ax.set_xlabel("Year of publication"); ax.set_ylabel("T_m of iron at the ICB (K)")
save(fig, "f_timeline")

# ---- uncertainty budget at 330 GPa ----
items = [("Statistics: 6 runs, 100 ps windows", 13, "±"), ("Potential vs DFT, liquid − hcp\n(40 unseen test structures)", 13, "±"),
         ("Electronic temperature term\n(first order, same frames)", 6, "±"), ("DFT pressure scale\n(≤ 2 GPa × 8.4 K/GPa)", 17, "±"),
         ("Elastic strain of solid slab\n(estimate, one-sided)", 20, "+10 to +")]
fig, ax = plt.subplots(figsize=(8.8, 4.6))
y = np.arange(len(items))[::-1]
for yi, (lab, v, sgn) in zip(y, items):
    if sgn == "±":
        ax.barh(yi, 2 * v, left=-v, height=0.55, color="#B9C7DA")
        ax.text(v + 1.5, yi, f"± {v} K", va="center", fontsize=12)
    else:
        ax.barh(yi, 10, left=10, height=0.55, color="#E7B8AC")
        ax.text(21.5, yi, "+10 to +20 K", va="center", fontsize=12)
ax.axvline(0, color=INK, lw=1)
ax.set_yticks(y, [i[0] for i in items], fontsize=11.5); ax.set_xlim(-32, 40); ax.set_xlabel("Effect on T_m at 330 GPa (K)")
ax.grid(axis="y", visible=False)
q = np.sqrt(13**2 + 13**2 + 6**2 + 17**2)
ax.set_title(f"Quadrature sum of the ± terms: ± {q:.0f} K.  Finite size: not yet measured.", fontsize=12)
save(fig, "f_budget")
print("quad", q)
