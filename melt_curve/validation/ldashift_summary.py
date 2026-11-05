"""ldashift_summary.py - report numbers from the 330 GPa snapshot test (ldashift.py eval on ARC).

Input: studies/lda_shift/ldashift_result.json and md_meta.json, written on ARC by
`python ldashift.py eval` from 80 QE 7.5 runs (40 snapshots x PBE, LDA; raw outputs not committed,
sha256 in data/). Energies in the input are eV/atom, pressures GPa, volumes A^3/atom.
Output: out/ldashift.json, in meV/atom and K.

Two results:
- potential vs PBE at 330 GPa, 6,329 K: pressure error, force RMSE, liquid - hcp energy bias and the
  T_m shift it implies (first order, split into pressure and equal-volume parts, plus second order);
- LDA - PBE shift of T_m by free-energy perturbation (a1 + a2 + a3 per phase). a3 = -V dP^2 / (2B)
  depends on the two bulk moduli of the potential (NVT at V(1 +- 0.01)); the equal-modulus variant
  uses their mean for both phases and shows how much of the result rests on their difference.
Conversion to K at dS_m = 0.80 k_B, as in ldashift.py.
"""
import json, os
import numpy as np
import paths

KB = 8.617333262e-5                  # eV/K
DS = 0.80 * KB                       # eV/K per atom
EV_A3_PER_GPA = 1.0 / 160.21766208   # eV/A^3 per GPa
SRC = getattr(paths, "LDASHIFT", os.path.join(paths.ROOT, "studies", "lda_shift"))


def a3(V, B, dP):
    return -(V / (2.0 * B * EV_A3_PER_GPA)) * (dP * EV_A3_PER_GPA) ** 2


def main():
    r = json.load(open(os.path.join(SRC, "ldashift_result.json")))
    meta = json.load(open(os.path.join(SRC, "md_meta.json")))
    h, l, m = r["results"]["hcp"], r["results"]["liquid"], r["mace_validation"]
    K = lambda e: e / DS
    d1, d2, d3 = l["a1"] - h["a1"], l["a2"] - h["a2"], l["a3"] - h["a3"]
    se1 = float(np.hypot(l["a1_se"], h["a1_se"]))
    Bc = 0.5 * (h["B"] + l["B"])
    d3c = a3(l["V"], Bc, l["dP"]) - a3(h["V"], Bc, h["dP"])
    out = dict(
        n_hcp=h["n"], n_liquid=l["n"], V_hcp=h["V"], V_liquid=l["V"], dV_m=l["V"] - h["V"],
        B_hcp=h["B"], B_liquid=l["B"], P_static_pbe_hcp=h["P_pbe"], P_static_pbe_liquid=l["P_pbe"],
        # potential vs PBE
        mace_dP_hcp=m["hcp"]["dP"], mace_dP_hcp_se=m["hcp"]["dP_se"],
        mace_dP_liquid=m["liquid"]["dP"], mace_dP_liquid_se=m["liquid"]["dP_se"],
        mace_F_rmse_hcp_meV=1e3 * m["hcp"]["F_rmse"], mace_F_rmse_liquid_meV=1e3 * m["liquid"]["F_rmse"],
        mace_F_rms_pbe_hcp=m["hcp"]["F_rms_pbe"], mace_F_rms_pbe_liquid=m["liquid"]["F_rms_pbe"],
        mace_F_rel_hcp_pct=100 * m["hcp"]["F_rmse"] / m["hcp"]["F_rms_pbe"],
        mace_F_rel_liquid_pct=100 * m["liquid"]["F_rmse"] / m["liquid"]["F_rms_pbe"],
        mace_dE_liq_hcp_meV=1e3 * m["dE_bias"], mace_dE_liq_hcp_se_meV=1e3 * m["dE_bias_se"],
        mace_Tm_first_order_K=m["T_m_shift_K"], mace_Tm_se_K=m["T_m_shift_se_K"],
        mace_Tm_pressure_part_K=m["T_m_pressure_part_K"], mace_Tm_equal_volume_part_K=m["T_m_structure_part_K"],
        mace_Tm_second_order_K=-m["second_order"] / DS, mace_Tm_total_K=m["T_m_shift_total_K"],
        # LDA - PBE
        lda_dP_hcp=h["dP"], lda_dP_liquid=l["dP"],
        lda_d_a1_meV=1e3 * d1, lda_d_a1_se_meV=1e3 * se1, lda_d_a2_meV=1e3 * d2, lda_d_a3_meV=1e3 * d3,
        lda_dG_meV=1e3 * (d1 + d2 + d3),
        lda_dTm_first_order_K=K(d1), lda_dTm_a2_K=K(d2), lda_dTm_a3_K=K(d3),
        lda_dTm_K=K(d1 + d2 + d3), lda_dTm_se_K=K(se1),
        lda_pressure_part_K=r["dT_pressure_part"], lda_equal_volume_part_K=r["dT_structure_part"],
        lda_B_common=Bc, lda_d_a3_equalB_meV=1e3 * d3c, lda_dTm_equalB_K=K(d1 + d2 + d3c))
    assert abs(out["lda_dTm_K"] - r["dT"]) < 1e-6 and abs(out["lda_dTm_se_K"] - r["dT_se"]) < 1e-6
    assert abs(meta["hcp"]["B_GPa"] - h["B"]) < 1e-9 and abs(meta["liquid"]["B_GPa"] - l["B"]) < 1e-9
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "ldashift.json"), "w"), indent=1)
    print(f"potential - PBE: T_m {out['mace_Tm_total_K']:+.0f} +- {out['mace_Tm_se_K']:.0f} K "
          f"(first order {out['mace_Tm_first_order_K']:+.0f}: pressure {out['mace_Tm_pressure_part_K']:+.0f}, "
          f"equal volume {out['mace_Tm_equal_volume_part_K']:+.0f}; second order {out['mace_Tm_second_order_K']:+.0f})")
    print(f"LDA - PBE: dT_m {out['lda_dTm_K']:+.0f} +- {out['lda_dTm_se_K']:.0f} K (a1 {out['lda_dTm_first_order_K']:+.0f}, "
          f"a2 {out['lda_dTm_a2_K']:+.0f}, a3 {out['lda_dTm_a3_K']:+.0f}); equal B = {Bc:.0f} GPa: "
          f"a3 {out['lda_d_a3_equalB_meV']:+.2f} meV/atom, dT_m {out['lda_dTm_equalB_K']:+.0f} K")


if __name__ == "__main__":
    main()
