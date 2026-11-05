#!/usr/bin/env python3
"""
compare_curves.py — Tel-6000 curve against the T_el-conditioned curve (A26).

Inputs (written by phaseE_analyze_live.py --include-drifting, which already adds
the kinetic pressure term, A27):
  phaseE/analysis_live_all.json           Tel-6000 campaign (phaseE/runs)
  phaseE/analysis_live_all_runs_tel.json  Tel-model campaign (phaseE/runs_tel)

Reports, per cell:
  1. both melting points, and the difference at matching pressure (the Tel-6000
     Simon fit evaluated at the Tel-model pressure);
  2. the electronic temperature used (the Tel-6000 Tm, as in fe_phaseE_tel.sbatch)
     and the mismatch Tm(Tel model) - T_el used;
  3. the sensitivity dTm/dT_el measured from the two campaigns where T_el moved
     by more than 200 K, and the shift a second self-consistency round would
     produce: sensitivity x mismatch. A26: a second round is needed if any cell
     exceeds 15 K;
  4. the DFT-bin expectation of the shift (A23, liquid-minus-hcp energy change
     on training frames near coexistence), where available.
Then compares both fits with published values and writes phaseE/compare_curves.png.

Run (login node): python3.13 compare_curves.py
"""
import json, os
import numpy as np

A6, AT = "phaseE/analysis_live_all.json", "phaseE/analysis_live_all_runs_tel.json"
SECOND_ROUND_K = 15.0
P_REF = 100.0

# A23: DFT-bin first-order dTm (K) at T_el = 4000 and 8000 K relative to 6000 K,
# keyed by cell; frames near (not at) coexistence, liquid sampled hotter than hcp.
DFT_BINS = {"coex_hcp_100GPa": (-70, 185), "coex_hcp_150GPa": (-58, 161),
            "coex_hcp_250GPa": (183, -132), "coex_hcp_330GPa": (249, -218)}

# Published values (Mendeley keys); P in GPa, T in K, error or None.
LIT = [("Wu2024DeepPotential", 330, 6415, None), ("Sun2023BccHcp", 323, 6357, 45),
       ("Alfe2009Coexistence", 328, 6370, 100), ("Rao2025DMFTMelting", 330, 6225, 42),
       ("Sola2009DMC", 330, 6900, 400), ("Anzellini2013Melting", 330, 6230, 500)]


def simon(P, f):
    return f["T1"] * (1.0 + (np.asarray(P, float) - P_REF) / f["a"]) ** (1.0 / f["c"])


# Cells without a DFT bin use the nearest binned pressure (sensitivity bound only).
BIN_FOR = {"coex_hcp_100GPa": "coex_hcp_100GPa", "coex_fcc_100GPa": "coex_hcp_100GPa",
           "coex_hcp_150GPa": "coex_hcp_150GPa", "coex_hcp_200GPa": "coex_hcp_150GPa",
           "coex_hcp_250GPa": "coex_hcp_250GPa", "coex_hcp_300GPa": "coex_hcp_330GPa",
           "coex_hcp_330GPa": "coex_hcp_330GPa", "coex_hcp_350GPa": "coex_hcp_330GPa"}


def dft_expect(cell, tel, nearest=False):
    key = BIN_FOR.get(cell) if nearest else cell
    if key not in DFT_BINS:
        return None
    d4, d8 = DFT_BINS[key]
    p = np.polyfit([4000.0, 6000.0, 8000.0], [d4, 0.0, d8], 2)
    return float(np.polyval(p, tel))


def main():
    for f in (A6, AT):
        if not os.path.exists(f):
            raise SystemExit(f"missing {f}: run phaseE_analyze_live.py first "
                             f"({'--rundir phaseE/runs_tel ' if f == AT else ''}--include-drifting)")
    a6, at = json.load(open(A6)), json.load(open(AT))
    c6, ct = a6["cells"], at["cells"]
    if not a6.get("fit"):
        raise SystemExit("Tel-6000 analysis has no Simon fit")

    print(f"{'cell':16s} {'P6000':>6s} {'Tm6000':>11s} {'Ptel':>6s} {'Tm_tel':>11s} "
          f"{'dTm@P':>6s} {'Tel used':>8s} {'mismatch':>8s} {'dTm/dTel':>8s} {'2nd rnd':>7s} {'DFT exp':>7s} {'ratio':>7s}")
    flags, rows = [], []
    for cell in sorted(c6, key=lambda c: (c.split("_")[1] != "hcp", int(c.split("_")[2][:-3]))):
        r6 = c6[cell]
        tel = round(r6["Tm"])                       # what fe_phaseE_tel.sbatch used
        if cell not in ct:
            print(f"{cell:16s} {r6['P']:6.1f} {r6['Tm']:6.0f}+-{r6['sem']:3.0f}   (no converged Tel-model runs yet)")
            continue
        rt = ct[cell]
        if r6["branch"] == "hcp":
            d_at_p = rt["Tm"] - float(simon(rt["P"], a6["fit"]))
        else:   # fcc: shift the Tel-6000 fcc point to the Tel-model pressure with the hcp Clapeyron slope
            f = a6["fit"]; Pm = 0.5 * (rt["P"] + r6["P"])
            slope = f["T1"] / (f["c"] * f["a"]) * (1 + (Pm - P_REF) / f["a"]) ** (1 / f["c"] - 1)
            d_at_p = rt["Tm"] - (r6["Tm"] + slope * (rt["P"] - r6["P"]))
        mism = rt["Tm"] - tel
        dT = tel - 6000.0
        # Campaign ratio: includes the model difference (Tel model vs Tel-6000 potential),
        # so it is NOT a T_el sensitivity; printed for information only.
        ratio = d_at_p / dT if abs(dT) > 200 else None
        # Second round = rerun with T_el moved from `tel` to the new Tm, same model, same P:
        # estimated from the DFT-bin T_el dependence (A23), nearest binned pressure.
        e0, e1 = dft_expect(cell, tel, True), dft_expect(cell, rt["Tm"], True)
        second = (e1 - e0) if e0 is not None else None
        sens = second / mism if (second is not None and abs(mism) > 1) else None
        exp = dft_expect(cell, tel)
        if second is not None and abs(second) > SECOND_ROUND_K:
            flags.append(f"{cell}: second-round shift {second:+.0f} K exceeds {SECOND_ROUND_K:.0f} K")
        rows.append((rt["P"], rt["Tm"], rt["sem"], r6["branch"]))
        print(f"{cell:16s} {r6['P']:6.1f} {r6['Tm']:6.0f}+-{r6['sem']:3.0f} {rt['P']:6.1f} "
              f"{rt['Tm']:6.0f}+-{rt['sem']:3.0f} {d_at_p:+6.0f} {tel:8d} {mism:+8.0f} "
              f"{(f'{sens:+.3f}' if sens is not None else '-'):>8s} "
              f"{(f'{second:+.0f}' if second is not None else '-'):>7s} "
              f"{(f'{exp:+.0f}' if exp is not None else '-'):>7s} "
              f"{(f'{ratio:+.3f}' if ratio is not None else '-'):>7s}")

    print("\nColumns: dTm@P = Tm(Tel model) - Tel-6000 fit at the same pressure (fcc: Tel-6000 point moved to the Tel-model pressure with the hcp fit slope);")
    print("mismatch = Tm(Tel model) - T_el used; dTm/dTel = DFT-bin (A23) sensitivity at fixed P, nearest binned pressure;")
    print("2nd rnd = DFT-bin shift for moving T_el from 'Tel used' to Tm(Tel model);")
    print("DFT exp = A23 DFT-bin expectation of the T_el part of dTm@P (binned pressures only);")
    print("ratio = dTm@P / (T_el used - 6000): includes the model difference, not a sensitivity.")
    print("\nSELF-CONSISTENCY:", "no second round needed" if not flags else "")
    for f in flags:
        print("  WARN", f)

    print("\nLITERATURE (K):")
    print(f"{'reference':24s} {'P':>5s} {'Tm':>11s} {'Tel-6000':>9s} {'Tel model':>9s}")
    ft = at.get("fit")
    for key, P, T, e in LIT:
        s6 = f"{float(simon(P, a6['fit'])) - T:+.0f}"
        st = f"{float(simon(P, ft)) - T:+.0f}" if ft else "-"
        print(f"{key:24s} {P:5d} {T:6d}{('+-%d' % e) if e else '':5s} {s6:>9s} {st:>9s}")
    print("(this work minus reference; DMC and DMFT include physics beyond PBE)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        Pp = np.linspace(95, 365, 200)
        ax.plot(Pp, simon(Pp, a6["fit"]), "-", color="0.4", lw=1.2, label="Tel-6000, Simon fit")
        h6 = [(v["P"], v["Tm"], v["sem"]) for v in c6.values() if v["branch"] == "hcp"]
        ax.errorbar(*zip(*h6), fmt="o", color="0.4", mfc="white", capsize=3, label="Tel-6000, hcp")
        if ft:
            ax.plot(Pp, simon(Pp, ft), "-", color="tab:red", lw=1.2, label="T_el-conditioned, Simon fit")
        ht = [(p, t, e) for p, t, e, b in rows if b == "hcp"]
        if ht:
            ax.errorbar(*zip(*ht), fmt="o", color="tab:red", capsize=3, label="T_el-conditioned, hcp")
        for key, P, T, e in LIT:
            ax.errorbar(P, T, yerr=e or 0, fmt="D", ms=4, capsize=2, label=key)
        ax.set_xlabel("P (GPa)"); ax.set_ylabel("T$_m$ (K)")
        ax.legend(fontsize=6, loc="lower right")
        fig.tight_layout(); fig.savefig("phaseE/compare_curves.png", dpi=200)
        print("\nwritten phaseE/compare_curves.png")
    except Exception as e:
        print(f"plot skipped: {e}")


if __name__ == "__main__":
    main()
