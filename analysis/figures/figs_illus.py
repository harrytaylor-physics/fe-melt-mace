"""Illustrations for the group-meeting deck (schematics, no simulation data except f_scale)."""
import sys
sys.path.insert(0, '.')
from figstyle import *
from refs import r
import numpy as np
from matplotlib.patches import Wedge, Polygon, Rectangle, FancyArrowPatch, Circle, FancyBboxPatch


def clean(ax):
    ax.set_axis_off(); ax.set_aspect("equal")


# ---------------- 1. Earth cross-section (PREM) ----------------
fig, ax = plt.subplots(figsize=(7.6, 6.4))
R = {"surf": 6371, "cmb": 3480, "icb": 1221}
ax.add_patch(Wedge((0, 0), R["surf"], 90, 270, facecolor="#E8E2D8", edgecolor="#8A8F96", lw=1))
ax.add_patch(Wedge((0, 0), R["cmb"], 90, 270, facecolor="#F3CFAE", edgecolor="#8A8F96", lw=1))
ax.add_patch(Wedge((0, 0), R["icb"], 90, 270, facecolor="#D0785A", edgecolor=OURS, lw=2.2))
ax.plot([0, 0], [-R["surf"], R["surf"]], color="#8A8F96", lw=1)
ax.text(-4300, 3300, "Mantle\n(solid silicate)", ha="center", va="center", fontsize=13)
ax.text(-2150, -1700, "Outer core\n(liquid Fe alloy)", ha="center", va="center", fontsize=13)
ax.annotate("Inner core\n(solid Fe alloy)", xy=(-700, 500), xytext=(-2300, 1500), ha="center", va="center", fontsize=13,
            arrowprops=dict(arrowstyle="-", color=INK, lw=1))
rows = [(R["surf"], "Surface", "depth 0 km", "0 GPa"),
        (R["cmb"], "Core–mantle boundary", "depth 2,891 km", "136 GPa"),
        (R["icb"], "Inner-core boundary (ICB)", "depth 5,150 km", "329 GPa"),
        (0, "Centre", "depth 6,371 km", "364 GPa")]
for y, name, d, p in rows:
    hi = name.startswith("Inner")
    ax.plot([0, 900], [y, y], color=OURS if hi else "#8A8F96", lw=1.6 if hi else 1)
    ax.text(1050, y, f"{name}\n{d},  {p}", va="center", fontsize=12.5, color=OURS if hi else INK,
            weight="bold" if hi else "normal", linespacing=1.25)
ax.text(1050, -1500, "At the ICB solid and liquid iron alloy\nare in contact, so T(ICB) = T_m of the\nalloy at 329 GPa.",
        va="top", fontsize=12, color=INK, linespacing=1.3)
ax.text(1050, -6100, f"Radii and pressures: PREM {r('PREM')}", fontsize=10.5, color=MUTED)
ax.set_xlim(-6600, 6400); ax.set_ylim(-6500, 6600); clean(ax)
fig.savefig(FIG + "f_earth.png"); plt.close(fig)

# ---------------- 2. Fermi-Dirac occupation at the three T_el ----------------
kB = 8.617333e-5
E = np.linspace(-3, 3, 600)
fig, ax = plt.subplots(figsize=(6.6, 4.2))
ax.plot(E, (E < 0).astype(float), color=GREY, lw=1.4, ls=":", label="T_el = 0 K")
for T, c in ((4000, THEORY), (6000, EXPT), (8000, OURS)):
    ax.plot(E, 1 / (np.exp(E / (kB * T)) + 1), color=c, lw=2.2, label=f"T_el = {T:,} K   (k_B T_el = {kB * T:.2f} eV)")
ax.set_xlabel("E − μ (eV)"); ax.set_ylabel("Occupation f(E)")
ax.set_xlim(-3, 3); ax.set_ylim(-0.03, 1.05)
ax.legend(loc="upper right", fontsize=10, framealpha=1, edgecolor="#D0D4D9")
ax.text(-2.9, 0.41, r"$f(E) = 1\,/\,[\,e^{(E-\mu)/k_B T_{\rm el}} + 1\,]$", fontsize=12.5, color=INK, va="center")
ax.text(-2.9, 0.22, "probability that a Kohn–Sham state\nat energy E is occupied;\nμ = chemical potential (Fermi level)", fontsize=10, color=MUTED, va="center", linespacing=1.25)
fig.savefig(FIG + "f_fermi.png"); plt.close(fig)

# ---------------- 3. Free energy vs T: superheating / undercooling ----------------
fig, ax = plt.subplots(figsize=(7.0, 4.4))
T = np.linspace(0.6, 1.5, 300)
Gs = -1.0 * (T - 1); Gl = -1.9 * (T - 1)
ax.plot(T[T <= 1], Gs[T <= 1], color=THEORY, lw=2.6); ax.plot(T[T >= 1], Gs[T >= 1], color=THEORY, lw=1.8, ls="--")
ax.plot(T[T >= 1], Gl[T >= 1], color=OURS, lw=2.6); ax.plot(T[T <= 1], Gl[T <= 1], color=OURS, lw=1.8, ls="--")
ax.axvspan(1.06, 1.43, color=THEORY, alpha=0.10, lw=0)
ax.text(1.245, 0.44, f"heated perfect crystal melts\nat 1.06–1.43 T_m {r('Luo2003')}", ha="center", fontsize=10.5, color=THEORY)
ax.axvline(1, color=GREY, lw=1)
ax.text(1.01, 0.88, "T_m: G_solid = G_liquid", fontsize=10.5, color=INK)
ax.text(0.62, 0.0, "solid (stable)", color=THEORY, fontsize=11.5, weight="bold")
ax.text(1.03, -0.66, "liquid (stable)", color=OURS, fontsize=11.5, weight="bold")
ax.text(1.33, -0.24, "solid, superheated\n(metastable)", color=THEORY, fontsize=10.5)
ax.text(0.70, 0.74, "liquid, undercooled\n(metastable)", color=OURS, fontsize=10.5)
ax.set_xlabel("T / T_m"); ax.set_ylabel("Gibbs free energy (arb. units)")
ax.set_xlim(0.6, 1.5); ax.set_ylim(-1.0, 1.0); ax.set_yticks([]); ax.grid(False)
fig.savefig(FIG + "f_hysteresis.png"); plt.close(fig)

# ---------------- 4. Stress in the slab: P_n vs P_mean ----------------
fig, ax = plt.subplots(figsize=(9.4, 3.9))
L, H = 10.0, 3.0
ax.add_patch(Rectangle((0, 0), 3, H, facecolor="#F3CFAE", edgecolor="none"))
ax.add_patch(Rectangle((3, 0), 4, H, facecolor="#C9D6EA", edgecolor="none"))
ax.add_patch(Rectangle((7, 0), 3, H, facecolor="#F3CFAE", edgecolor="none"))
ax.add_patch(Rectangle((0, 0), L, H, facecolor="none", edgecolor="#8A8F96", lw=1.2))
for x in (3, 7):
    ax.plot([x, x], [0, H], color=INK, lw=1.2, ls="--")
ax.text(1.5, 2.55, "liquid", ha="center", fontsize=13, color=OURS, weight="bold")
ax.text(8.5, 2.55, "liquid", ha="center", fontsize=13, color=OURS, weight="bold")
ax.text(5.0, 2.55, "solid (hcp)", ha="center", fontsize=13, color=THEORY, weight="bold")


def arrow(x0, y0, x1, y1, c, lw=2.2):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=15, color=c, lw=lw))


# normal (z) stress: equal on both sides of each interface
for x in (3, 7):
    arrow(x - 1.05, 1.15, x - 0.08, 1.15, INK); arrow(x + 1.05, 1.15, x + 0.08, 1.15, INK)
ax.text(5.0, 1.30, "σ_zz is continuous across\nthe interface: −σ_zz = P_liquid", ha="center", fontsize=10.5, va="bottom")
# in-plane (x) stress in the solid: different
arrow(6.4, 0.08, 6.4, 0.62, THEORY); arrow(6.4, H - 0.08, 6.4, H - 0.62, THEORY)
ax.text(3.15, 0.18, "in-plane σ_xx = σ_yy\nset by the box width", fontsize=10.5, color=THEORY, va="bottom")
ax.annotate("", xy=(10.6, H), xytext=(10.6, 0), arrowprops=dict(arrowstyle="-|>", color=GREY))
ax.annotate("", xy=(L, -0.45), xytext=(0, -0.45), arrowprops=dict(arrowstyle="-|>", color=GREY))
ax.text(5.0, -0.85, "z (interface normal, 145 Å)", ha="center", fontsize=11, color=MUTED)
ax.text(10.75, 1.5, "x, y\n(25.9 Å)", va="center", fontsize=11, color=MUTED)
ax.set_xlim(-0.2, 12.0); ax.set_ylim(-1.1, 3.3); clean(ax)
fig.savefig(FIG + "f_slab.png"); plt.close(fig)

# ---------------- 5. MACE neighbourhood ----------------
rng = np.random.default_rng(3)
pts = []
for i in range(-6, 7):
    for j in range(-6, 7):
        x = i + 0.5 * j; y = j * np.sqrt(3) / 2
        if abs(x) < 3.1 and abs(y) < 2.9:
            pts.append((x, y))
pts = np.array(pts) + rng.normal(0, 0.07, (len(pts), 2))
pts[np.argmin(np.hypot(*pts.T))] = (0, 0)
rc = 2.15
fig, ax = plt.subplots(figsize=(7.4, 4.8))
ax.add_patch(Circle((0, 0), rc, facecolor="#EEF3FA", edgecolor=THEORY, lw=1.6, ls="--"))
d = np.hypot(*pts.T)
for (x, y), dd in zip(pts, d):
    if 0 < dd < rc:
        ax.plot([0, x], [0, y], color="#9DB5D8", lw=1.0, zorder=1)
for (x, y), dd in zip(pts, d):
    ax.add_patch(Circle((x, y), 0.22 if dd == 0 else 0.18, color=OURS if dd == 0 else (THEORY if dd < rc else "#C4C8CE"), zorder=3))
j = int(np.argmin(np.where((pts[:, 0] > 0.3) & (pts[:, 1] > 0.3) & (d < rc), np.abs(d - 1.9), 9)))
ax.add_patch(FancyArrowPatch((0, 0), tuple(pts[j]), arrowstyle="-|>", mutation_scale=14, color=INK, lw=2.0, zorder=4))
lab = dict(fontsize=12, va="center", ha="left")
ln = dict(arrowstyle="-", color="#8A8F96", lw=0.9)
ax.annotate("atom i (energy E_i)", xy=(0.2, -0.1), xytext=(3.6, -0.6), color=OURS, arrowprops=ln, **lab)
ax.annotate("vector r_ij to neighbour j", xy=tuple(pts[j] * 0.6), xytext=(3.6, 1.9), color=INK, arrowprops=ln, **lab)
ax.annotate("neighbours inside r_c", xy=(1.0, -1.75), xytext=(3.6, -1.6), color=THEORY, arrowprops=ln, **lab)
ax.annotate("cutoff sphere, radius r_c", xy=(rc * np.cos(0.55), rc * np.sin(0.55)), xytext=(3.6, 0.9), color=THEORY, arrowprops=ln, **lab)
ax.annotate("atoms outside r_c\n(reached via message passing)", xy=(2.75, -2.55), xytext=(3.6, -2.6), color=MUTED, arrowprops=ln, **lab)
ax.set_xlim(-3.3, 8.6); ax.set_ylim(-3.1, 3.1); clean(ax)
fig.savefig(FIG + "f_mace_nbhd.png"); plt.close(fig)

# ---------------- 6. Diamond-anvil cell ----------------
fig, ax = plt.subplots(figsize=(5.6, 4.6))
up = np.array([(-1.6, 2.4), (1.6, 2.4), (1.9, 1.75), (0.32, 0.16), (-0.32, 0.16), (-1.9, 1.75)])
ax.add_patch(Polygon(up, closed=True, facecolor="#E9F1F8", edgecolor="#7E93AE", lw=1.4))
ax.add_patch(Polygon(up * [1, -1], closed=True, facecolor="#E9F1F8", edgecolor="#7E93AE", lw=1.4))
ax.add_patch(Rectangle((-2.6, -0.16), 2.38, 0.32, facecolor="#9AA0A6")); ax.add_patch(Rectangle((0.22, -0.16), 2.38, 0.32, facecolor="#9AA0A6"))
ax.add_patch(Rectangle((-0.16, -0.07), 0.32, 0.14, facecolor=OURS))
for s in (1, -1):
    ax.annotate("", xy=(-0.05, 0.08 * s), xytext=(-1.25, 2.9 * s), arrowprops=dict(arrowstyle="-|>", color="#D64F2A", lw=2.2))
ax.annotate("", xy=(0, 0), xytext=(0, -3.0), arrowprops=dict(arrowstyle="-", color=BEYOND, lw=2.6))
for ang in (12, -12):
    t = np.deg2rad(ang); ax.plot([0, 3.0 * np.sin(t)], [0, 3.0 * np.cos(t)], color=BEYOND, lw=1.4, ls="--")
ax.annotate("", xy=(1.55, 3.0), xytext=(0.1, 0.1), arrowprops=dict(arrowstyle="-|>", color="#B8860B", lw=1.8))
ax.text(-1.35, 3.0, "IR heating laser\n(both sides)", ha="right", va="center", fontsize=10.5, color="#D64F2A")
ax.text(0.15, -2.95, "X-rays in", fontsize=10.5, color=BEYOND, va="bottom")
ax.text(0.0, 3.08, "diffraction", fontsize=10.5, color=BEYOND, va="bottom", ha="center")
ax.text(1.7, 3.05, "thermal emission\n→ T (spectroradiometry)", fontsize=10.5, color="#8A6508", va="bottom")
ax.text(2.0, 1.2, "diamond", fontsize=11, color="#4E627C")
ax.text(2.7, -0.0, "gasket", fontsize=11, color=MUTED, va="center")
ax.text(0.45, -0.55, "Fe sample", fontsize=11, color=OURS)
ax.set_xlim(-3.3, 4.6); ax.set_ylim(-3.2, 4.0); clean(ax)
fig.savefig(FIG + "f_ap_dac.png"); plt.close(fig)

# ---------------- 7. Dynamic (shock / ramp) compression ----------------
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.add_patch(Rectangle((0.0, 0), 0.7, 3, facecolor="#D9D3C7", edgecolor="#8A8F96"))
ax.add_patch(Rectangle((0.7, 0), 2.4, 3, facecolor="#C9D6EA", edgecolor="#8A8F96"))
ax.add_patch(Rectangle((1.7, 0), 1.4, 3, facecolor="#9FB6D6", edgecolor="none"))
ax.add_patch(Rectangle((3.1, 0), 1.2, 3, facecolor="#F4F6F8", edgecolor="#8A8F96"))
ax.plot([1.7, 1.7], [0, 3], color=INK, lw=2.0)
ax.annotate("", xy=(2.25, 1.5), xytext=(1.75, 1.5), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6))
for y in (0.7, 1.5, 2.3):
    ax.annotate("", xy=(-0.02, y), xytext=(-1.4, y), arrowprops=dict(arrowstyle="-|>", color="#D64F2A", lw=2.2))
ax.annotate("", xy=(4.32, 2.2), xytext=(5.6, 2.2), arrowprops=dict(arrowstyle="<|-|>", color=EXPT, lw=1.8))
ax.plot([0.9, 3.9], [-0.9, 3.9], color=BEYOND, lw=2.0, ls="-")
ax.text(-1.45, 3.15, "drive laser", fontsize=10.5, color="#D64F2A")
ax.text(0.35, -0.25, "ablator", fontsize=10.5, ha="center", va="top", color=MUTED)
ax.text(1.2, 3.12, "Fe", fontsize=12, ha="center", color=INK, weight="bold")
ax.text(2.4, 3.12, "compressed Fe", fontsize=10.5, ha="center", color=INK)
ax.text(3.7, -0.25, "window", fontsize=10.5, ha="center", va="top", color=MUTED)
ax.text(1.95, 1.62, "shock\nfront", fontsize=10, va="bottom", color=INK)
ax.text(4.5, 2.45, "velocimetry (VISAR)\n→ pressure", fontsize=10.5, color=EXPT, va="bottom")
ax.text(3.95, 3.95, "X-ray probe\n→ XRD / XAS", fontsize=10.5, color=BEYOND, va="top", ha="left")
ax.set_xlim(-1.6, 6.6); ax.set_ylim(-1.0, 4.2); clean(ax)
fig.savefig(FIG + "f_ap_shock.png"); plt.close(fig)

# ---------------- 8. Cell sizes at the same scale ----------------
from ase.io import read
at = read(DATA + "cells/coex_hcp_330GPa.extxyz")
q = np.load("q6_330.npy"); assert len(q) == len(at)
at.wrap()
P = at.get_positions(); F = at.get_scaled_positions(); Lz = at.cell[2, 2]; Hy = at.cell[1, 1]
vat = at.get_volume() / len(at)
fig, ax = plt.subplots(figsize=(11.0, 3.0))
sel = F[:, 0] < 0.25
ax.scatter(P[sel, 2], P[sel, 1], s=1.6, c=np.where(q[sel] > 0.142, THEORY, OURS), lw=0)
ax.add_patch(Rectangle((0, 0), Lz, Hy, fill=False, edgecolor=GREY, lw=1))
x = -6.0
for N, lab in ((980, f"N = 980\nAlfè 2009 {r('Alfe2009')}"), (250, f"N = 250\nSun 2023 {r('Sun2023')}")):
    Lb = (N * vat) ** (1 / 3); x -= Lb
    ax.add_patch(Rectangle((x, 0), Lb, Lb, facecolor="#E9EDF3", edgecolor=GREY, lw=1))
    ax.text(x + Lb / 2, -3, lab, ha="center", va="top", fontsize=10.5)
    x -= 6.0
ax.text(Lz / 2, -3, f"this work: N = {len(at):,}; hexagonal cell, a = {at.cell.lengths()[0]:.1f} Å, c = {Lz:.0f} Å  (blue: solid, red: liquid, by q̄6)",
        ha="center", va="top", fontsize=10.5)
ax.plot([Lz - 20, Lz], [Hy + 3, Hy + 3], color=INK, lw=2)
ax.text(Lz - 10, Hy + 4.2, "20 Å", ha="center", va="bottom", fontsize=10.5)
Lc = [0, 0, Lz]; nax = 2; oth = 1
ax.set_xlim(x + 2, Lz + 2); ax.set_ylim(-11, Hy + 10); clean(ax)
fig.savefig(FIG + "f_scale.png"); plt.close(fig)
print("done; vat", vat)
