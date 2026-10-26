import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OURS, THEORY, EXPT, BEYOND, GREY, INK, MUTED = "#C0442A", "#2F6DB5", "#1B8A6B", "#7A4FB0", "#6B7280", "#1F2328", "#555B63"
TI_COL = {0.92: "#2F6DB5", 1.00: "#1B8A6B", 1.08: "#C0442A"}
plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 13, "axes.labelsize": 14, "axes.titlesize": 14,
    "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 11.5, "legend.frameon": False,
    "axes.edgecolor": "#8A8F96", "axes.linewidth": 0.8, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.labelcolor": INK, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": "#E6E8EB", "grid.linewidth": 0.7, "lines.linewidth": 2.0,
    "savefig.dpi": 220, "savefig.bbox": "tight", "savefig.facecolor": "white", "figure.facecolor": "white",
})
import os
FIG = os.environ.get("FE_FIG_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "inputs")) + "/"   # intermediate json and figure output
DATA = os.environ.get("FE_PHASEE", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "phaseE")) + "/"   # coexistence analysis outputs
