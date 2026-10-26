import sys, json, re
sys.path.insert(0, '.')
from figstyle import *
import numpy as np
M = json.load(open(FIG + "train_meta.json"))
pn = json.load(open(DATA + "pn_correct.json"))["runs_tel"]["fit"]
simon = lambda P: pn["T1"] * (1 + (np.asarray(P) - 100) / pn["a"]) ** (1 / pn["c"])
PH = [("hcp", THEORY, "o"), ("liquid", OURS, "o"), ("stackfault", BEYOND, "s"), ("fcc", EXPT, "^"), ("bcc", "#B8860B", "D")]
NAME = {"hcp": "hcp", "liquid": "liquid", "stackfault": "stacking fault", "fcc": "fcc", "bcc": "bcc"}

# ---- coverage map (T_el = 6000 K set; the 4000/8000 K sets use the same geometries) ----
fig, ax = plt.subplots(figsize=(8.6, 5.4))
rows = M["train_Tel6000"] + M["test_Tel6000"]
from collections import Counter
off = {"hcp": (-2.5, 0), "liquid": (2.5, 0), "fcc": (0, 0), "bcc": (0, 0)}
for ph, c, m in PH:
    r = [x for x in rows if x["phase"] == ph and x["P"] == x["P"]]
    if not r: continue
    cnt = Counter((round(x["P"], 1), round(x["T"], 0)) for x in r)
    P = np.array([k[0] for k in cnt]) + off[ph][0]; T = np.array([k[1] for k in cnt]); n = np.array(list(cnt.values()))
    ax.scatter(P, T, s=8 + 3.0 * n, marker=m, color=c, alpha=0.6, lw=0, label=f"{NAME[ph]} ({len(r)} frames, {len(cnt)} runs)")
sp = M["spot_frames"]
ax.scatter([x["P"] for x in sp], [x["T"] for x in sp], s=60, facecolor="none", edgecolor=INK, lw=1.4, label="unseen test structures (20)")
PP = np.linspace(95, 375, 200)
ax.plot(PP, simon(PP), color=INK, lw=2.2, label="melting curve (this work)")
ax.axvspan(329, 364, color=GREY, alpha=0.10, lw=0)
ax.text(346.5, 2350, "inner\ncore", ha="center", fontsize=10.5, color=MUTED)
ax.set_xlabel("Pressure of sampling run (GPa)"); ax.set_ylabel("Ionic temperature of sampling run (K)")
ax.set_xlim(40, 470); ax.set_ylim(1800, 9800)
ax.legend(loc="upper left", fontsize=9.5, markerscale=0.8, framealpha=1, edgecolor="#D0D4D9")
fig.savefig(FIG + "f_coverage.png"); plt.close(fig)

# ---- force magnitude distribution by phase ----
edges = np.linspace(0, 12, 61); mid = 0.5 * (edges[1:] + edges[:-1])
fig, ax = plt.subplots(figsize=(6.0, 4.0))
for ph, c, m in PH:
    h = np.sum([x["fhist"] for x in rows if x["phase"] == ph], axis=0).astype(float)
    ax.plot(mid, h / h.sum() / (edges[1] - edges[0]), color=c, lw=2, label=NAME[ph])
ax.set_xlabel("|F| per atom, DFT (eV/Å)"); ax.set_ylabel("Probability density (Å/eV)")
ax.set_xlim(0, 10); ax.legend(fontsize=10)
fig.savefig(FIG + "f_forces.png"); plt.close(fig)
for ph, c, m in PH:
    r = [x for x in rows if x["phase"] == ph]
    print(ph, len(r), "mean|F| %.2f" % np.mean([x["fmean"] for x in r]), "P %.0f-%.0f" % (min(x["P"] for x in r), max(x["P"] for x in r)),
          "T %.0f-%.0f" % (min(x["T"] for x in r), max(x["T"] for x in r)))
allF = np.sum([x["fhist"] for x in rows], axis=0); print("overall mean |F| approx", (allF * mid).sum() / allF.sum())

# ---- training curve from the MACE log ----
ep, E, F = [], [], []
for l in open("arc/logs/fe_tel_e3_run-7.log"):
    m = re.search(r"Epoch (\d+): .*RMSE_E_per_atom=\s*([\d.]+) meV, RMSE_F=\s*([\d.]+)", l)
    if m:
        ep.append(int(m.group(1))); E.append(float(m.group(2))); F.append(float(m.group(3)))
d = {}
for e, a, b in zip(ep, E, F): d[e] = (a, b)
ep = sorted(d); E = [d[e][0] for e in ep]; F = [d[e][1] for e in ep]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.0, 3.8))
a1.plot(ep, E, color=THEORY); a1.set_yscale("log"); a1.set_xlabel("Epoch"); a1.set_ylabel("Validation energy RMSE (meV/atom)")
a2.plot(ep, F, color=OURS); a2.set_xlabel("Epoch"); a2.set_ylabel("Validation force RMSE (meV/Å)"); a2.set_ylim(85, 130)
for a in (a1, a2):
    a.axvline(110, color=GREY, ls="--", lw=1)
a1.text(112, 200, "stage two:\nenergy weight\n1 → 1,000", fontsize=10, color=MUTED, va="top")
a1.annotate(f"{E[-1]:.1f} meV/atom", (ep[-1], E[-1]), xytext=(-30, 55), textcoords="offset points", ha="right", fontsize=10.5,
            arrowprops=dict(arrowstyle="-", color=GREY, lw=0.8))
a2.annotate(f"{F[-1]:.1f} meV/Å", (ep[-1], F[-1]), xytext=(-30, 45), textcoords="offset points", ha="right", fontsize=10.5,
            arrowprops=dict(arrowstyle="-", color=GREY, lw=0.8))
fig.tight_layout(); fig.savefig(FIG + "f_train_curve.png"); plt.close(fig)
print("epochs", ep[0], ep[-1], "E", E[0], E[-1], "F", F[0], F[-1], "E@109", d[109])

# ---- same-frame electronic-temperature check (e3_test_final.out) ----
txt = open("arc/logs/e3_test_final.out").read()
rows = re.findall(r"\s(\d+)GPa\s+(\d+)\s+(\d+)\s+\|\s+([+-][\d.]+)\+-([\d.]+)\s+([+-][\d.]+)\+-([\d.]+)\s+\|\s+([+-][\d.]+)\+-([\d.]+)\s+([+-][\d.]+)\+-([\d.]+)", txt)
lab = [f"{r[0]} GPa" for r in rows]; x = np.arange(len(rows)); w = 0.2
fig, ax = plt.subplots(figsize=(7.6, 4.2))
sets = [(3, 4, "DFT, 4,000 K", "#9DB5D8", None), (5, 6, "model, 4,000 K", THEORY, None), (7, 8, "DFT, 8,000 K", "#E7B8AC", None), (9, 10, "model, 8,000 K", OURS, None)]
for k, (iv, ie, name, c, _) in enumerate(sets):
    v = [float(r[iv]) for r in rows]; e = [float(r[ie]) for r in rows]
    ax.bar(x + (k - 1.5) * w, v, w * 0.92, yerr=e, color=c, label=name, capsize=2, error_kw=dict(lw=0.8, ecolor=INK))
ax.axhline(0, color=INK, lw=0.8)
ax.set_xticks(x, lab); ax.set_ylabel("Liquid minus hcp of\nF(T_el) − F(6,000 K) (meV/atom)")
ax.legend(ncol=2, fontsize=10, loc="lower left", framealpha=1, edgecolor="#D0D4D9"); ax.grid(axis="x", visible=False)
fig.savefig(FIG + "f_tel_check.png"); plt.close(fig)
json.dump(rows, open(FIG + "tel_check.json", "w"))
print(rows)
