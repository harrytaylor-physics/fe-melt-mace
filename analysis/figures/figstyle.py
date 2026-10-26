import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OURS, THEORY, EXPT, BEYOND, GREY, INK, MUTED = "#C0442A", "#2F6DB5", "#1B8A6B", "#7A4FB0", "#6B7280", "#1F2328", "#555B63"
TI_COL = {0.92: "#2F6DB5", 1.00: "#1B8A6B", 1.08: "#C0442A"}
plt.rcParams.update({
    "font.family": ["Liberation Sans", "Arial", "DejaVu Sans"], "font.size": 13, "axes.labelsize": 14, "axes.titlesize": 14,
    "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 11.5, "legend.frameon": False,
    "axes.edgecolor": "#8A8F96", "axes.linewidth": 0.8, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.labelcolor": INK, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#E6E8EB", "grid.linewidth": 0.7, "lines.linewidth": 2.0,
    "savefig.dpi": 220, "savefig.bbox": "tight", "savefig.facecolor": "white", "figure.facecolor": "white",
})
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import layout  # noqa: E402

FIG = layout.FIG_INPUTS + "/"       # intermediate json (analysis/figures/inputs)
DATA = layout.PHASEE + "/"          # coexistence analysis outputs (phaseE)
OUT = layout.FIGURES + "/"          # figures (docs/figures)
plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none"})


def save(fig, name):
    """Write docs/figures/<name>.pdf (vector) and <name>.png (200 dpi)."""
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(OUT + name + ".pdf", metadata={"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None})
    fig.savefig(OUT + name + ".png", dpi=200, metadata={"Software": None})
    plt.close(fig)
