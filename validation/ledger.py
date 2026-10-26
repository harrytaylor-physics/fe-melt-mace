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
  cc_check.py     -> out/cc_check.json       needs phaseE/ (analysis json, cells, runs_tel checkpoints)
  pbe_vs_exp.py   -> out/pbe_vs_exp.json     needs refs/eos_pbe_cold/*.json (and phaseE/pn_correct.json)
"""
import hashlib, json, os
import paths

PRODUCER = {"test_stats": "stats_test.py", "static_compare": "val_mace.py",
            "cc_check": "cc_check.py", "pbe_vs_exp": "pbe_vs_exp.py"}


def flatten(x, pre, out):
    if isinstance(x, dict):
        for k, v in x.items():
            flatten(v, f"{pre}.{k}", out)
    elif isinstance(x, list) and x and isinstance(x[0], (dict, list)):
        for i, v in enumerate(x):
            name = v.get("cell") if isinstance(v, dict) and "cell" in v else (
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
    return s


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
    led = dict(keys=keys, summary=summary(J), inputs=inputs)
    json.dump(led, open(paths.LEDGER, "w"), indent=1, sort_keys=True)
    print(f"{len(keys)} keys, {len(led['summary'])} summary entries from {len(inputs)} producer outputs -> {paths.LEDGER}")
    for k, v in led["summary"].items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
