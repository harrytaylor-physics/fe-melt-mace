"""Parity plot: MACE (fe_tel_e3_stagetwo) vs DFT forces on the held-out test labels (parity_tel.npz from ARC)."""
import sys, json
sys.path.insert(0, '.')
from figstyle import *
import numpy as np
from matplotlib.colors import LogNorm
z = np.load(sys.argv[1], allow_pickle=True)
Fr, Fm, fph, ftel = z["F_ref"], z["F_mod"], z["f_phase"], z["f_tel"]
Er, Em, eph, etel = z["E_ref"], z["E_mod"], z["e_phase"], z["e_tel"]
x = Fr.ravel(); y = Fm.ravel()
rm = lambda a, b: 1000 * np.sqrt(np.mean((a - b) ** 2))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.8, 4.9), gridspec_kw=dict(width_ratios=[1.15, 1]))
lim = np.percentile(np.abs(x), 99.95) * 1.05
hb = a1.hexbin(x, y, gridsize=170, cmap="Greys", norm=LogNorm(), mincnt=1, extent=(-lim, lim, -lim, lim), linewidths=0, zorder=2)
a1.plot([-lim, lim], [-lim, lim], color=OURS, lw=0.8, ls="--", zorder=0)
a1.set_xlim(-lim, lim); a1.set_ylim(-lim, lim); a1.set_aspect("equal")
a1.set_xlabel("DFT force component (eV/Å)"); a1.set_ylabel("MACE force component (eV/Å)")
a1.text(0.04, 0.96, f"{len(Er)} structures, {len(x):,} components\nRMSE {rm(x, y):.1f} meV/Å", transform=a1.transAxes, va="top", fontsize=11.5)
cb = fig.colorbar(hb, ax=a1, fraction=0.046, pad=0.03); cb.set_label("components per bin")
# right: error distribution per T_el
bins = np.linspace(-0.6, 0.6, 121)
for tel, c in ((4000, THEORY), (6000, EXPT), (8000, OURS)):
    s = ftel == tel
    d = (Fm[s] - Fr[s]).ravel()
    a2.hist(d, bins=bins, histtype="step", lw=2, color=c, density=True, label=f"T_el = {tel:,} K\nRMSE {rm(Fm[s], Fr[s]):.0f} meV/Å")
a2.set_xlabel("MACE − DFT, force component (eV/Å)"); a2.set_ylabel("probability density (Å/eV)")
a2.legend(fontsize=10, loc="upper left", framealpha=1, edgecolor="#D0D4D9"); a2.set_xlim(-0.6, 0.6); a2.set_ylim(0, 10)
fig.tight_layout(); fig.savefig(FIG + "f_parity.png"); plt.close(fig)
# per-phase table
rows = []
for ph in ["hcp", "liquid", "stackfault", "fcc", "bcc", "all"]:
    r = dict(phase=ph)
    for tel in (4000, 6000, 8000):
        fs = (ftel == tel) & ((fph == ph) | (ph == "all")); es = (etel == tel) & ((eph == ph) | (ph == "all"))
        off = np.mean(Em[etel == tel] - Er[etel == tel])
        r[f"F{tel}"] = round(rm(Fm[fs], Fr[fs]), 1) if fs.any() else None
        r[f"E{tel}"] = round(1000 * np.sqrt(np.mean((Em[es] - Er[es] - off) ** 2)), 2) if es.any() else None
        r[f"n{tel}"] = int(es.sum())
    rows.append(r)
json.dump(dict(rows=rows, n=int(len(Er)), ncomp=int(len(x)), F_all=round(rm(x, y), 1)), open(FIG + "parity.json", "w"), indent=1)
for r in rows: print(r)
