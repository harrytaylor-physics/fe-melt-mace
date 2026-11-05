"""coex_cells.py - the coexistence-cell table of the T_el-conditioned campaign.

Input: data/coexistence/analysis_live_all_runs_tel.json (phaseE_analyze_live.py --rundir phaseE/runs_tel
--include-drifting), cells and runs. P is the cell-mean pressure (kinetic term included), before
the interface-normal correction of pn_correct.py; T_m and its standard error over the runs of a cell;
dS_m in k_B (cell value; per-run values where the analysis could form them).
Output: out/coex_cells.json - per cell: P, P_sem, T_m, T_m_sem, scatter, dS_kB, runs, verdicts,
per-run dS range; campaign totals (runs, ps, verdict counts) and the analysis's own Simon fit.
"""
import json, os, sys
from collections import Counter
import paths

sys.path.insert(0, os.path.dirname(paths.HERE))
import layout  # noqa: E402


def main():
    A = json.load(open(layout.ANALYSIS_TEL))
    runs = A["runs"]
    cells = {}
    for c, v in A["cells"].items():
        rs = [r for r in runs.values() if r["cell"] == c]
        ds = [r["dS"] for r in rs if r.get("dS") is not None]
        cells[c] = dict(branch=v["branch"], P_nominal=v["Pnom"], P=v["P"], P_sem=v["Psem"], Tm=v["Tm"], Tm_sem=v["sem"],
                        scatter=v["scatter"], dS_kB=v.get("dS_kB"), Tm_est=v["Tm_est"], runs=v["n"],
                        verdicts=dict(Counter(r["verdict"] for r in rs)),
                        dS_runs_range=[min(ds), max(ds)] if ds else None, n_dS_runs=len(ds))
    out = dict(cells=cells, n_runs=len(runs), ps_total=sum(r["ps"] for r in runs.values()),
               verdicts=dict(Counter(r["verdict"] for r in runs.values())),
               hcp_runs_per_cell=sorted({v["runs"] for v in cells.values() if v["branch"] == "hcp"}),
               fit_cell_mean_P=A["fit"], flags=A["flags"])
    for c, v in cells.items():
        print(f"{c:18s} P {v['P']:6.1f}  T_m {v['Tm']:6.0f} +- {v['Tm_sem']:4.1f}  dS {v['dS_kB']:.3f}  runs {v['runs']}  {v['verdicts']}")
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "coex_cells.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
