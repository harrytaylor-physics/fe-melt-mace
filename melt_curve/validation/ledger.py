"""ledger.py - every number of the potential-validation report as a named key.

Reads the producer outputs in out/ and writes ledger.json:
  {"keys": {key: value}, "summary": {key: value}, "inputs": {file: {sha256, producer}}}
"keys" holds every leaf of every out/*.json, named <file stem>.<path>. "summary" holds the
headline numbers derived from them (largest errors, ranges). A number quoted in a report must
be one of these keys; a number that is not here needs a producer first.

Producers (run from this directory):
  stats_test.py   -> out/test_stats.json     needs parity_tel.npz
  val_dft.py      -> refs/dft/*.json         needs pw.x and the PBE PAW file (about 3 core-hours per structure)
  val_mace.py     -> out/static_compare.json needs torch, mace-torch, the model file
  cc_check.py     -> out/cc_check.json       needs data/coexistence/ (analysis json, cells, runs_tel checkpoints)
  pbe_vs_exp.py   -> out/pbe_vs_exp.json     needs refs/eos_pbe_cold/*.json (and data/coexistence/pn_correct.json)
  train_table.py  -> out/training.json       needs the e3 training logs, data/training/dataset_report.txt, data/training/train_e3.xyz
  tel_check.py    -> out/tel_check.json      needs production_campaign/arc_scripts/logs/e3_test_final.out
  coex_cells.py   -> out/coex_cells.json     needs data/coexistence/analysis_live_all_runs_tel.json
  melting_point.py -> out/melting_point.json needs data/coexistence/pn_correct.json
  repro_check.py  -> out/repro_check.json    needs torch, mace-torch, the model, data/training/test_Tel*.xyz, parity_tel.npz
  settings.py     -> out/settings.json       needs refs/*/*.in, the producers' sources, production_campaign/fe_melt/pseudo/, data/coexistence/cells, production_campaign/fe_melt/labels/ldaX_*.in
  conv_check.py   -> out/convergence.json    needs refs/convergence/*.json (QE 6.6 runs)
  lda_setup_check.py -> out/lda_setup.json   needs refs/lda_setup/*.json (QE 6.6 runs)
  ldashift_summary.py -> out/ldashift.json   needs studies/lda_shift/ldashift_result.json, md_meta.json (ldashift.py eval on ARC)
Committed constants (no producer): refs/literature.json -> keys "literature.*".
"""
import hashlib, json, math, os
import paths

PRODUCER = {"test_stats": "stats_test.py", "static_compare": "val_mace.py",
            "cc_check": "cc_check.py", "pbe_vs_exp": "pbe_vs_exp.py",
            "training": "train_table.py", "tel_check": "tel_check.py", "coex_cells": "coex_cells.py",
            "melting_point": "melting_point.py", "repro_check": "repro_check.py", "settings": "settings.py",
            "convergence": "conv_check.py", "lda_setup": "lda_setup_check.py", "ldashift": "ldashift_summary.py"}
CONSTANTS = {"literature": paths.LITERATURE}
KB_MEV = 8.617333262e-2                                   # k_B in meV/K


def flatten(x, pre, out):
    if isinstance(x, dict):
        for k, v in x.items():
            flatten(v, f"{pre}.{k}", out)
    elif isinstance(x, list) and x and isinstance(x[0], (dict, list)):
        for i, v in enumerate(x):
            name = v.get("cell") if isinstance(v, dict) and "cell" in v else v["id"] if isinstance(v, dict) and "id" in v else v["tag"] if isinstance(v, dict) and "tag" in v else (
                f"{v['P_lo']:.0f}-{v['P_hi']:.0f}GPa" if isinstance(v, dict) and "P_lo" in v else (
                    f"V{v[0]:.2f}" if isinstance(v, list) else (f"V{v['V']:.2f}" if isinstance(v, dict) and "V" in v else str(i))))
            flatten(v, f"{pre}.{name}", out)
    else:
        out[pre] = x


def summary(J):
    s = {}
    if "static_compare" in J:
        sc = J["static_compare"]
        eos = sc["eos"]                                   # [V, P_dft, P_mace, dP, dE_dft, dE_mace]
        s["eos_max_abs_dP_GPa"] = max(abs(r[3]) for r in eos)
        s["eos_max_abs_dE_meV"] = max(abs(r[5] - r[4]) for r in eos)
        s["eos_P_range_GPa"] = [min(r[1] for r in eos), max(r[1] for r in eos)]
        s["eos_E_span_eV"] = (max(r[4] for r in eos) - min(r[4] for r in eos)) / 1e3
        ph = {k: v for k, v in sc.items() if k.endswith(("_7.0", "_8.0")) and "-hcp" in k}
        s["phase_max_abs_diff_meV"] = max(abs(v[1] - v[0]) for v in ph.values())
        s["phase_min_abs_diff_meV"] = min(abs(v[1] - v[0]) for v in ph.values())
        if "elastic" in sc:
            e = sc["elastic"]
            rel = {k: (e["mace"][k] - e["dft"][k]) / e["dft"][k] for k in e["dft"]}
            s["elastic_rel_err"] = rel
            s["elastic_max_abs_rel_err_excl_C66"] = max(abs(v) for k, v in rel.items() if k != "C66")
            if "E2g_DFT" in sc:
                s["E2g_rel_err"] = (sc["E2g_MACE"] - sc["E2g_DFT"]) / sc["E2g_DFT"]
        if "bulk_modulus" in sc:
            b = sc["bulk_modulus"]
            s["bulk_modulus_max_abs_rel_err"] = max(abs(b["mace"]["K_GPa"][v] / b["dft"]["K_GPa"][v] - 1) for v in b["dft"]["K_GPa"])
    if "cc_check" in J:
        iv = J["cc_check"]["intervals"]
        hi = [r for r in iv if r["P_lo"] >= 250]
        s["cc_max_abs_z_above_258GPa"] = max(abs(r["z"]) for r in hi)
        lo = [r for r in iv if r["P_hi"] <= 260]
        s["cc_rel_diff_range_below_258GPa"] = [min(r["rel_diff"] for r in lo), max(r["rel_diff"] for r in lo)]
    if "test_stats" in J:
        t = J["test_stats"]
        s["test_E_rmse_range"] = [min(t[f"{T}_all"]["rmse"] for T in (4000, 6000, 8000)),
                                  max(t[f"{T}_all"]["rmse"] for T in (4000, 6000, 8000))]
        s["test_F_rmse_range_meV_A"] = [min(t[f"{T}_all"]["frmse"] for T in (4000, 6000, 8000)),
                                        max(t[f"{T}_all"]["frmse"] for T in (4000, 6000, 8000))]
        lh = [t[f"{T}_liquid-hcp"][0] for T in (4000, 6000, 8000)]
        s["test_liquid_minus_hcp_range_meV"] = [min(lh), max(lh)]
        s["test_liquid_minus_hcp_mean_meV"] = sum(lh) / 3
    if "pbe_vs_exp" in J:
        p = J["pbe_vs_exp"]
        s["pbe_minus_exp_core_GPa"] = p["core_range"]
        s["Tm_true_minus_pbe_K"] = p["dTm_true_minus_pbe_range"]
        s["dTm_dP_330_K_GPa"] = p["dtdp_330"]
    extra(J, s)
    followup(J, s)
    return s


def extra(J, s):
    """Derived numbers added on 26 Oct 2026 for the quantities the report states as ranges or ratios."""
    s["K_per_meV_at_dS_0.80"] = 1.0 / (0.80 * KB_MEV)
    if "test_stats" in J:
        t = J["test_stats"]
        s["test_F_rel_range"] = [min(t[f"{T}_all"]["frel"] for T in (4000, 6000, 8000)),
                                 max(t[f"{T}_all"]["frel"] for T in (4000, 6000, 8000))]
        s["test_liquid_minus_hcp_mean_K"] = s["test_liquid_minus_hcp_mean_meV"] * s["K_per_meV_at_dS_0.80"]
        for T in (4000, 6000, 8000):
            o = [t[f"{T}_{p}"] for p in ("fcc", "bcc", "stackfault")]
            s[f"test_{T}_fcc_bcc_sf"] = dict(n=sum(x["n"] for x in o), rmse=[min(x["rmse"] for x in o), max(x["rmse"] for x in o)],
                                           frmse=[min(x["frmse"] for x in o), max(x["frmse"] for x in o)],
                                           frel=[min(x["frel"] for x in o), max(x["frel"] for x in o)])
        s["test_n_frames"] = sum(t[f"{T}_all"]["n"] for T in (4000, 6000, 8000))
    if "static_compare" in J:
        sc = J["static_compare"]
        if "tel_response" in sc:
            tr = sc["tel_response"]
            s["tel_response_max_abs_dP_GPa"] = max(abs(v["P_mace"] - v["P_dft"]) for v in tr.values())
            s["tel_response_max_abs_dF_meV"] = max(abs(v["dF_mace_meV"] - v["dF_dft_meV"]) for v in tr.values())
        s["phase_E_diff_meV"] = {k: v[1] - v[0] for k, v in sc.items() if "-hcp_" in k}
        s["phase_P_diff_GPa"] = {k: v[1] - v[0] for k, v in sc.items() if "_P_" in k}
        if "coa_min" in sc:
            s["coa_rel_err"] = sc["coa_min"]["mace"] / sc["coa_min"]["dft"] - 1
        allp = [abs(r[3]) for r in sc["eos"]] + [abs(v) for v in s["phase_P_diff_GPa"].values()]
        if "tel_response" in sc:
            allp.append(s["tel_response_max_abs_dP_GPa"])
        s["static_max_abs_dP_all_GPa"] = max(allp)
        v = sc["voigt"]
        s["voigt_rel_err"] = {k: v["mace"][k] / v["dft"][k] - 1 for k in ("K_V", "G_V", "vp_km_s", "vs_km_s")}
        b = sc["bulk_modulus"]["dft"]["K_GPa"]["7.00"], sc["bulk_modulus"]["mace"]["K_GPa"]["7.00"]
        s["bulk_modulus_7.00_rel_err"] = b[1] / b[0] - 1
    if "cc_check" in J:
        iv = J["cc_check"]["intervals"]
        s["cc_P_range_consistent_GPa"] = [min(r["P_lo"] for r in iv if r["P_lo"] >= 250), max(r["P_hi"] for r in iv)]
    if "pbe_vs_exp" in J:
        rows = J["pbe_vs_exp"]["rows"]
        lo = max(rows, key=lambda r: r["V"])                     # largest volume, lowest pressure
        s["pbe_minus_exp_lowest_P"] = dict(V=lo["V"], P_exp=lo["P_exp"], diff=lo["diff"])
        spread = lambda r: max(r["P_exp"].values()) - min(r["P_exp"].values())
        s["exp_spread_below_55GPa_max_GPa"] = max(spread(r) for r in rows if max(r["P_exp"].values()) < 55)
        s["exp_spread_at_V_core_GPa"] = spread([r for r in rows if abs(r["V"] - J["pbe_vs_exp"]["V_core"]) < 1e-6][0])
        s["P_vib_range_GPa"] = [min(r["P_vib"] for r in rows), max(r["P_vib"] for r in rows)]
        # pressure where PBE + vib - experiment changes sign (linear between computed volumes), per scale
        z = {}
        for k in rows[0]["diff"]:
            pts = sorted((r["P_exp"][k], r["diff"][k]) for r in rows)
            for (p1, d1), (p2, d2) in zip(pts[:-1], pts[1:]):
                if d1 * d2 <= 0 and d1 != d2:
                    z[k] = p1 + (0 - d1) * (p2 - p1) / (d2 - d1)
                    break
        s["pbe_minus_exp_zero_crossing_P_GPa"] = z
    if "melting_point" in J and "test_stats" in J and "static_compare" in J:
        m = J["melting_point"]
        s["Tm_330_K"] = m["Tm_330"]
        s["Tm_error_budget_K"] = (s["eos_max_abs_dP_GPa"] * m["dTdP_330"]
                                  + s["test_liquid_minus_hcp_mean_meV"] * s["K_per_meV_at_dS_0.80"])
    if "coex_cells" in J and "cc_check" in J:
        lit = J.get("literature", {}).get("melting_330GPa", [])
        dS_lit = [x["dS"] for x in lit if x.get("dS")]
        if dS_lit:
            dS330 = J["coex_cells"]["cells"]["coex_hcp_330GPa"]["dS_kB"]
            s["dS_330_rel_to_literature"] = [dS330 / max(dS_lit) - 1, dS330 / min(dS_lit) - 1]


def followup(J, s):
    """Derived numbers added on 29 Oct 2026: the 330 GPa snapshot test and the statistical error of T_m(330 GPa)."""
    if "melting_point" in J:
        m = J["melting_point"]
        s["Tm_330_stat_err_K"] = m["at"]["330"]["Tm_fit_se"]          # standard error of the Simon fit at 330 GPa
    if "test_stats" in J:
        t = J["test_stats"]
        s["test_liquid_minus_hcp_mean_se_meV"] = sum(t[f"{T}_liquid-hcp"][1] ** 2 for T in (4000, 6000, 8000)) ** 0.5 / 3
    if "ldashift" in J:
        x = J["ldashift"]
        s["direct_test_dE_rel_hcp_meV"] = dict(hcp=0.0, liquid=x["mace_dE_liq_hcp_meV"])
        s["lda_dP_mean_GPa"] = 0.5 * (x["lda_dP_hcp"] + x["lda_dP_liquid"])
        # largest |dT_m| within one standard error, over the two bulk-modulus treatments
        s["lda_dTm_max_abs_within_1se_K"] = max(abs(x["lda_dTm_K"]), abs(x["lda_dTm_equalB_K"])) + x["lda_dTm_se_K"]
        if "melting_point" in J:
            # PBE melting temperature: the potential's T_m minus its error found by the direct test
            s["Tm_PBE_330_K"] = J["melting_point"]["Tm_330"] - x["mace_Tm_total_K"]
            s["Tm_PBE_330_err_K"] = math.hypot(s["Tm_330_stat_err_K"], x["mace_Tm_se_K"])
            rz = {r["id"]: r for r in J.get("literature", {}).get("melting_330GPa", [])}.get("Rao2025")
            if rz:
                s["RaoZhu_minus_Tm_PBE_330_K"] = rz["Tm"] - s["Tm_PBE_330_K"]
                s["RaoZhu_minus_Tm_PBE_330_err_K"] = math.hypot(rz["Tm_err"], s["Tm_PBE_330_err_K"])


def main():
    J, keys, inputs = {}, {}, {}
    for stem in PRODUCER:
        f = os.path.join(paths.OUT, stem + ".json")
        if not os.path.exists(f):
            continue
        J[stem] = json.load(open(f))
        flatten(J[stem], stem, keys)
        inputs[os.path.relpath(f, os.path.dirname(paths.LEDGER))] = dict(
            sha256=hashlib.sha256(open(f, "rb").read()).hexdigest(), producer=PRODUCER[stem])
    for stem, f in CONSTANTS.items():
        J[stem] = json.load(open(f))
        flatten({k: v for k, v in J[stem].items() if not k.startswith("_")}, stem, keys)
        inputs[os.path.relpath(f, os.path.dirname(paths.LEDGER))] = dict(
            sha256=hashlib.sha256(open(f, "rb").read()).hexdigest(), producer="committed constants")
    led = dict(keys=keys, summary=summary(J), inputs=inputs)
    json.dump(led, open(paths.LEDGER, "w"), indent=1, sort_keys=True)
    print(f"{len(keys)} keys, {len(led['summary'])} summary entries from {len(inputs)} producer outputs -> {paths.LEDGER}")
    for k, v in led["summary"].items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
