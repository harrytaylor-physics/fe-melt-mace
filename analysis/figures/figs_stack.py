"""Stacking schematic: hcp ABAB, a stacking fault inside hcp, fcc ABCABC (close-packed layers seen edge-on)."""
import sys
sys.path.insert(0, '.')
from figstyle import *
import numpy as np
from matplotlib.patches import Rectangle
COL = {"A": THEORY, "B": EXPT, "C": OURS}
SHIFT = {"A": 0.0, "B": 1 / 3, "C": 2 / 3}
panels = [("hcp", "ABABABAB", None), ("stacking fault in hcp", "ABABCACA", (2.5, 4.5)), ("fcc", "ABCABCAB", None)]
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.9))
for ax, (title, seq, band) in zip(axes, panels):
    for k, L in enumerate(seq):
        y = k * 0.62
        xs = np.arange(-1, 5) + SHIFT[L]
        xs = xs[(xs > -0.6) & (xs < 4.4)]
        ax.scatter(xs, [y] * len(xs), s=260, color=COL[L], edgecolor="white", lw=1, zorder=3)
        ax.text(-1.0, y, L, ha="center", va="center", fontsize=13, weight="bold", color=COL[L])
    if band:
        ax.add_patch(Rectangle((-1.4, band[0] * 0.62 - 0.31 + 0.31), 6.0, (band[1] - band[0]) * 0.62, color=GREY, alpha=0.13, lw=0, zorder=1))
        ax.text(4.75, (band[0] + band[1]) / 2 * 0.62, "fcc-like\nlayers (B, C)", fontsize=10.5, va="center", color=INK)
    ax.set_title(title, fontsize=13)
    ax.set_xlim(-1.5, 5.9); ax.set_ylim(-0.5, 7 * 0.62 + 0.5); ax.set_axis_off()
axes[0].text(4.75, 3.5 * 0.62, "ABAB…", fontsize=11, va="center")
axes[2].text(4.75, 3.5 * 0.62, "ABCABC…", fontsize=11, va="center")
fig.text(0.5, 0.01, "Close-packed layers viewed edge-on; A, B, C are the three possible lateral positions of a layer.", ha="center", fontsize=10.5, color=MUTED)
fig.savefig(FIG + "f_stacking.png"); plt.close(fig)
print("ok")
