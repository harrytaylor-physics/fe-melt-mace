import sys, json, glob, os
sys.path.insert(0, '.')
from figstyle import *
import numpy as np
from ase.io import read
KB, EVA3 = 8.617333262e-5, 160.21766208
R = DATA + "runs_tel/"

def load(tag):
    z = np.load(R + tag + "_ckpt.npz", allow_pickle=True)
    V = abs(np.linalg.det(z["cell"])); N = len(z["positions"])
    T = z["T_t"]; P = z["P_t"] + N * KB * T / V * EVA3
    return dict(T=T, P=P, fr=z["fr_t"], N=N, V=V, z=z)

# ---------- F1: 330 GPa, six runs ----------
fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.6, 6.4), sharex=True, gridspec_kw=dict(hspace=0.12))
for ti in (0.92, 1.00, 1.08):
    for s in (0, 1):
        d = load(f"coex_hcp_330GPa_ti{ti:.2f}_s{s}")
        t = np.arange(len(d["T"])) + 1
        ls = "-" if s == 0 else (0, (4, 2))
        lab = f"start {ti:.2f} × 6,404 K" if s == 0 else None
        # 5-ps running mean for readability
        Ts = np.convolve(d["T"], np.ones(5) / 5, mode="valid")
        a1.plot(t[2:-2], Ts, color=TI_COL[ti], ls=ls, lw=1.6, label=lab)
        a2.plot(d["fr"][:, 0], d["fr"][:, 1], color=TI_COL[ti], ls=ls, lw=1.6, marker="o", ms=3.5)
a1.axvspan(100, 200, color="#F1F2F4", zorder=0)
a1.text(150, 6745, "averaging window\n100–200 ps", ha="center", va="top", color=MUTED, fontsize=11.5)
a1.set_ylabel("Kinetic temperature (K)\n5-ps running mean")
a1.set_ylim(6100, 6760)
h, l = a1.get_legend_handles_labels()
a2.legend(h, l, loc="lower center", ncol=3, fontsize=11, bbox_to_anchor=(0.5, -0.01))
a1.text(198, 6130, "solid line: seed 0;  dashed line: seed 1", ha="right", fontsize=11, color=MUTED)
a2.set_ylabel("Solid fraction\n(q̄6 > 0.142)")
a2.set_xlabel("Time (ps)"); a2.set_xlim(0, 200); a2.set_ylim(0.4, 1.0)
a2.axvspan(100, 200, color="#F1F2F4", zorder=0)
save(fig, "f_330_runs")

# ---------- F2: temperature extraction, one run ----------
d = load("coex_hcp_330GPa_ti1.00_s0"); T = d["T"]; N = d["N"]
t = np.arange(1, len(T) + 1)
w = T[100:]; Tm = w.mean()
blocks = np.array_split(w, 4); bm = [b.mean() for b in blocks]
sem = np.std(bm, ddof=1) / 2.0
sig_can = Tm * np.sqrt(2.0 / (3 * N))
fig = plt.figure(figsize=(10.4, 4.6))
ax = fig.add_axes([0.07, 0.14, 0.62, 0.78]); axh = fig.add_axes([0.72, 0.14, 0.2, 0.78], sharey=ax)
ax.plot(t, T, color=GREY, lw=0.9, marker="o", ms=2.6, label="sample every 1 ps (500 steps)")
for k, b in enumerate(blocks):
    x0 = 101 + 25 * k
    ax.hlines(b.mean(), x0, x0 + 24, color=OURS, lw=3, label="25-ps block means" if k == 0 else None)
ax.hlines(Tm, 101, 200, color=INK, lw=1.2, ls="--", label=f"mean {Tm:,.0f} ± {sem:.0f} K")
ax.axvspan(100, 200, color="#F1F2F4", zorder=0)
ax.set_xlabel("Time (ps)"); ax.set_ylabel("T = 2 E_kin / (3 N k_B)  (K)"); ax.set_xlim(0, 201)
ax.legend(loc="upper right", fontsize=10.5)
cnt, edges, _ = axh.hist(w, bins=18, orientation="horizontal", color="#D9DCE0", edgecolor="white")
bw = edges[1] - edges[0]
yy = np.linspace(w.min() - 60, w.max() + 60, 200)
axh.plot(len(w) * bw * np.exp(-(yy - Tm) ** 2 / (2 * sig_can ** 2)) / (sig_can * np.sqrt(2 * np.pi)), yy, color=OURS, lw=1.8)
axh.set_xlabel("count"); plt.setp(axh.get_yticklabels(), visible=False); axh.grid(False)
axh.text(0.98, 0.02, f"σ measured {w.std(ddof=1):.0f} K\nσ canonical {sig_can:.0f} K", transform=axh.transAxes, ha="right", va="bottom", fontsize=10.5, color=MUTED)
save(fig, "f_T_extract")
json.dump(dict(Tm=float(Tm), sem=float(sem), sig_meas=float(w.std(ddof=1)), sig_can=float(sig_can), bm=[float(x) for x in bm], N=int(N)),
          open(FIG + "f_T_extract.json", "w"), indent=1)

# ---------- F3: cell snapshot coloured by qbar6 + profile ----------
sys.path.insert(0, layout.Q6_MODULE_DIR)
from phaseE_q6 import q6_per_atom
at = read(DATA + "cells/coex_hcp_330GPa.extxyz")
z = np.load(R + "coex_hcp_330GPa_ti1.00_s0_ckpt.npz", allow_pickle=True)
at.set_cell(z["cell"]); at.set_positions(z["positions"]); at.pbc = True
at.wrap()
Q6 = FIG + "q6_330.npy"                # cache, not committed
q = np.load(Q6) if os.path.exists(Q6) else q6_per_atom(at)
np.save(Q6, q)
pos = at.get_positions(); L = at.cell.lengths()
fig = plt.figure(figsize=(12.0, 4.6))
gs = fig.add_gridspec(2, 2, height_ratios=[3.2, 0.9], width_ratios=[40, 0.8], hspace=0.08, wspace=0.03)
ax = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[1, 0], sharex=ax); cb = fig.add_subplot(gs[0, 1])
sel = (at.get_scaled_positions()[:, 1] < 0.5)
sc = ax.scatter(pos[sel, 2], pos[sel, 0], c=q[sel], s=5, cmap="cividis", vmin=0.04, vmax=0.32, linewidths=0)
ax.set_xlim(0, L[2]); ax.set_ylabel("x (Å)"); ax.grid(False)
plt.setp(ax.get_xticklabels(), visible=False)
c = fig.colorbar(sc, cax=cb); c.set_label("q̄6 (Lechner–Dellago, r_c = 3.2 Å)", fontsize=11)
zb = np.linspace(0, L[2], 61); zc = 0.5 * (zb[1:] + zb[:-1])
fs = [np.mean(q[(pos[:, 2] >= zb[k]) & (pos[:, 2] < zb[k + 1])] > 0.1421) for k in range(60)]
ax2.plot(zc, fs, color=OURS, lw=1.8); ax2.set_ylim(-0.1, 1.1); ax2.set_ylabel("solid\nfraction", fontsize=10.5)
ax2.set_yticks([0, 1]); ax2.set_xlabel("z, normal to the interfaces (Å)")
save(fig, "f_cell_q6")
json.dump(dict(L=[float(x) for x in L], frac=float((q > 0.1421).mean()), q_solid=float(np.median(q[q > 0.1421])),
               q_liq=float(np.median(q[q <= 0.1421]))), open(FIG + "f_cell_q6.json", "w"), indent=1)

# q6 distribution figure (bimodal)
fig, ax = plt.subplots(figsize=(5.4, 3.8))
ax.hist(q, bins=70, color="#B9BEC5", edgecolor="white")
ax.axvline(0.1421, color=OURS, lw=2); ax.text(0.146, ax.get_ylim()[1] * 0.9, "threshold 0.142", color=OURS, fontsize=11)
ax.set_xlabel("q̄6 per atom"); ax.set_ylabel("atoms"); ax.set_title("330 GPa cell, 12,096 atoms, t = 200 ps", fontsize=12)
save(fig, "f_q6_hist")
print("done", Tm, sem, sig_can, w.std(ddof=1))
