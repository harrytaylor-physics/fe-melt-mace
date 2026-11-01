"""Schematic for the DFT primer: all-electron vs pseudo (PAW) wavefunction and potential."""
import sys
sys.path.insert(0, '.')
from figstyle import *
import numpy as np
r = np.linspace(0.02, 4.0, 800); rc = 1.2
# all-electron-like radial function: oscillates in the core, decays outside
ae = np.sin(5.2 * r) * np.exp(-0.15 * r) * (r < rc) + 0
tail = lambda x: 0.9 * np.exp(-0.55 * (x - rc)) * np.cos(0.9 * (x - rc) + 0.25)
ae = np.where(r < rc, 0.9 * np.cos(0.25) * np.sin(7.0 * r) / np.sin(7.0 * rc) * (r / rc) ** 0.3, tail(r))
ps = np.where(r < rc, 0.9 * np.cos(0.25) * (1.6 * (r / rc) ** 2 - 0.6 * (r / rc) ** 4), tail(r))
fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 3.8))
a1.plot(r, ae, color=INK, lw=2, label="all-electron (nodes in the core)")
a1.plot(r, ps, color=OURS, lw=2.2, ls="--", label="smooth pseudo function")
a1.axvline(rc, color=GREY, lw=1); a1.text(rc + 0.05, -1.05, "r_c", color=MUTED, fontsize=12)
a1.set_xlabel("r (schematic)"); a1.set_ylabel("radial wavefunction"); a1.set_yticks([]); a1.set_xticks([])
a1.legend(fontsize=10.5, loc="upper right"); a1.set_ylim(-1.2, 1.4)
v = -2.0 / r
vp = np.where(r < rc, -2.0 / rc * (1.5 - 0.5 * (r / rc) ** 2), v)
a2.plot(r, v, color=INK, lw=2, label="−Z_eff / r (bare ion)")
a2.plot(r, vp, color=OURS, lw=2.2, ls="--", label="pseudopotential (finite at r = 0)")
a2.axvline(rc, color=GREY, lw=1); a2.text(rc + 0.05, -9.3, "r_c", color=MUTED, fontsize=12)
a2.set_ylim(-10, 0.5); a2.set_xlabel("r (schematic)"); a2.set_ylabel("potential"); a2.set_yticks([]); a2.set_xticks([])
a2.legend(fontsize=10.5, loc="lower right")
fig.tight_layout(); save(fig, "f_paw")
print("ok")
