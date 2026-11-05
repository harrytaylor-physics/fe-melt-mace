"""figs_validation.py - figures of the validation producers' outputs (melt_curve/validation/out/*.json).

f_heldout_errors   energy and force RMSE of the held-out test set by phase and T_el   (test_stats.json)
f_eos_error        MACE minus PBE pressure along the hcp EOS, and the energy error     (static_compare.json)
f_cc_check         slope of T_m(P) against dV_m/dS_m between neighbouring cells          (cc_check.json)
f_pbe_vs_exp       static PBE + vibrational pressure minus three 300 K experimental EOS  (pbe_vs_exp.json)
Reads only the committed producer outputs; writes nothing but the figures.
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *  # noqa: F401,F403
import numpy as np

J = lambda n: json.load(open(os.path.join(layout.VAL_OUT, n)))
TEL_COL = {4000: THEORY, 6000: EXPT, 8000: OURS}
PHASES = [("hcp", "hcp"), ("liquid", "liquid"), ("stackfault", "stacking fault"), ("fcc", "fcc"), ("bcc", "bcc"), ("all", "all")]

# ---- held-out test set: errors by phase ----
t = J("test_stats.json")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.2), gridspec_kw=dict(wspace=0.25))
x = np.arange(len(PHASES)); w = 0.26
for k, T in enumerate((4000, 6000, 8000)):
    e = [t[f"{T}_{p}"]["rmse"] for p, _ in PHASES]; f = [t[f"{T}_{p}"]["frmse"] for p, _ in PHASES]
    a1.bar(x + (k - 1) * w, e, w * 0.92, color=TEL_COL[T], label=f"T_el = {T:,} K")
    a2.bar(x + (k - 1) * w, f, w * 0.92, color=TEL_COL[T])
for a, yl in ((a1, "Energy RMSE (meV/atom)\nafter one offset per T_el"), (a2, "Force RMSE (meV/Å)")):
    a.set_xticks(x, [n for _, n in PHASES], rotation=20); a.set_ylabel(yl); a.grid(axis="x", visible=False)
a1.legend(fontsize=10, loc="upper left")
n = sum(t[f"{T}_all"]["n"] for T in (4000, 6000, 8000))
a1.set_title(f"Held-out test set: {n} frames", fontsize=12)
save(fig, "f_heldout_errors")

# ---- hcp EOS: pressure and energy error ----
sc = J("static_compare.json")
eos = np.array(sc["eos"])                      # V, P_dft, P_mace, dP, dE_dft, dE_mace
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 4.0), gridspec_kw=dict(wspace=0.3))
a1.axhspan(-0.5, 0.5, color="#F1F2F4", zorder=0, label="±0.5 GPa (PBE basis convergence)")
a1.plot(eos[:, 1], eos[:, 3], "o-", color=OURS, ms=6, label="MACE − PBE")
a1.axhline(0, color=INK, lw=0.8)
a1.set_xlabel("PBE pressure (GPa)"); a1.set_ylabel("P(MACE) − P(PBE)  (GPa)"); a1.set_ylim(-1, 1)
a1.legend(fontsize=10, loc="lower left")
a2.plot(eos[:, 0], eos[:, 5] - eos[:, 4], "o-", color=THEORY, ms=6)
a2.axhline(0, color=INK, lw=0.8)
a2.set_xlabel("Volume (Å³/atom)"); a2.set_ylabel("ΔE(MACE) − ΔE(PBE)  (meV/atom)\nrelative to 7.0 Å³/atom")
fig.suptitle("hcp, static lattice, c/a 1.60, T_el = 6,000 K", fontsize=12)
save(fig, "f_eos_error")

# ---- Clausius-Clapeyron ----
cc = J("cc_check.json")
iv = cc["intervals"]
Pm = np.array([0.5 * (r["P_lo"] + r["P_hi"]) for r in iv])
fig, ax = plt.subplots(figsize=(7.6, 4.6))
ax.errorbar(Pm - 1.5, [r["curve"] for r in iv], yerr=[r["curve_err"] for r in iv], fmt="o", ms=7, color=OURS, capsize=3,
            label="slope of T_m(P) between neighbouring cells")
ax.errorbar(Pm + 1.5, [r["cc"] for r in iv], yerr=[r["cc_err"] for r in iv], fmt="s", ms=7, color=THEORY, mfc="white", mew=1.8,
            capsize=3, label="ΔV_m / ΔS_m (Clausius-Clapeyron)")
for r, p in zip(iv, Pm):
    ax.text(p, max(r["curve"], r["cc"]) + max(r["curve_err"], r["cc_err"]) + 0.4, f"{r['z']:+.1f}σ", ha="center", fontsize=10, color=MUTED)
top = max(max(r["curve"] + r["curve_err"], r["cc"] + r["cc_err"]) for r in iv)
ax.set_ylim(top=top + 1.2)
ax.set_xlabel("Pressure, midpoint of the interval (GPa)"); ax.set_ylabel("dT_m/dP  (K/GPa)")
ax.legend(fontsize=10, loc="upper right")
save(fig, "f_cc_check")

# ---- PBE against experiment ----
pe = J("pbe_vs_exp.json")
rows = sorted(pe["rows"], key=lambda r: -r["V"])
fig, ax = plt.subplots(figsize=(7.6, 4.6))
for (k, lab), c, m in zip((("Dewaele2006", "Dewaele et al. 2006 (Vinet)"), ("Fei2016", "Fei et al. 2016 (BM3)"),
                           ("Hirao2022", "Hirao et al. 2022 (Vinet)")), (THEORY, EXPT, BEYOND), ("o", "s", "^")):
    ax.plot([r["P_exp"][k] for r in rows], [r["diff"][k] for r in rows], marker=m, ms=6, color=c, label=lab)
ax.axhline(0, color=INK, lw=0.8)
core = [r for r in rows if abs(r["V"] - pe["V_core"]) < 1e-6][0]
ax.axvspan(min(core["P_exp"].values()), max(core["P_exp"].values()), color="#F1F2F4", zorder=0)
ax.text(min(core["P_exp"].values()) - 3, ax.get_ylim()[0] + 1.0, f"{pe['V_core']:.1f} Å³/atom", ha="right", fontsize=10, color=MUTED)
ax.set_xlabel("Experimental pressure at the same volume (GPa)")
ax.set_ylabel("P(PBE static + vib) − P(exp)  (GPa)")
ax.legend(fontsize=10, loc="upper left")
save(fig, "f_pbe_vs_exp")
print("validation figures written to", OUT)
