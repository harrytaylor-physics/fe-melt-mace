"""tel_check.py - the same-frame electronic-temperature check of the production potential.

Input: production_campaign/arc_scripts/logs/e3_test_final.out (layout.TEL_CHECK_LOG), part 2, written on ARC by tel_diff_check.py
--tel-model phaseD/models/fe_tel_e3_stagetwo.model. Frames labelled at all three T_el; for each
pressure bin (+-25 GPa) the liquid-minus-hcp mean of dU = F(T_el) - F(6,000 K), DFT and model,
at T_el = 4,000 and 8,000 K.
Output: out/tel_check.json - per cell the four values with standard errors; the range of the DFT
values; the model-minus-DFT difference per cell and T_el, its smallest and largest magnitude, and
the largest converted to K with the dS_m the log states.
"""
import json, os, re, sys
import paths

sys.path.insert(0, os.path.dirname(paths.HERE))
import layout  # noqa: E402

KB_MEV = 8.617333262e-2          # k_B in meV/K


def main():
    txt = open(layout.TEL_CHECK_LOG).read()
    part = txt.split("PART 2: SAME-FRAME CHECK")[1]
    out = dict(n_frames_all_three_tel=int(re.search(r"frames labelled at all three Tels: (\d+)", part).group(1)),
               model=re.search(r"Tel model (\S+)", part).group(1),
               dS_kB=float(re.search(r"K at dS = ([\d.]+) kB", part).group(1)))
    for T, s in (("4000", re.search(r"Tel 4000: mean\s+([+-][\d.]+)\s+rmse after mean\s+([\d.]+)", part)),
                 ("8000", re.search(r"Tel 8000: mean\s+([+-][\d.]+)\s+rmse after mean\s+([\d.]+)", part))):
        out[f"dU_error_{T}"] = dict(mean=float(s.group(1)), rmse_after_mean=float(s.group(2)))
    cells = {}
    row = re.compile(r"^\s*(\d+)GPa\s+(\d+)\s+(\d+)\s+\|\s+([+-][\d.]+)\+-([\d.]+)\s+([+-][\d.]+)\+-([\d.]+)\s+\|"
                     r"\s+([+-][\d.]+)\+-([\d.]+)\s+([+-][\d.]+)\+-([\d.]+)", re.M)
    for m in row.finditer(part):
        g = m.groups()
        v = list(map(float, g[3:]))
        cells[f"{g[0]}GPa"] = dict(n_hcp=int(g[1]), n_liquid=int(g[2]),
                                   dft_4000=v[0], dft_4000_se=v[1], model_4000=v[2], model_4000_se=v[3],
                                   dft_8000=v[4], dft_8000_se=v[5], model_8000=v[6], model_8000_se=v[7],
                                   diff_4000=v[2] - v[0], diff_8000=v[6] - v[4])
    out["cells"] = cells
    dft = [c[k] for c in cells.values() for k in ("dft_4000", "dft_8000")]
    mod = [c[k] for c in cells.values() for k in ("model_4000", "model_8000")]
    dif = [abs(c[k]) for c in cells.values() for k in ("diff_4000", "diff_8000")]
    out["dft_range_meV"] = [min(dft), max(dft)]
    out["model_range_meV"] = [min(mod), max(mod)]
    out["abs_diff_range_meV"] = [min(dif), max(dif)]
    out["max_abs_diff_K"] = max(dif) / (out["dS_kB"] * KB_MEV)
    print(json.dumps(out, indent=1))
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "tel_check.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
