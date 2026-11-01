"""lda_setup_check.py - cost of running LDA (PZ) with the PBE PAW dataset.

The PZ PAW dataset generated with ld1.x from the pslibrary 1.0.0 'spn' input gives a divergent
SCF in displaced and liquid configurations at 330 GPa, so the LDA snapshot runs (ldashift.py)
use the official PBE dataset with input_dft = 'PZ'. These perfect-lattice runs at 7.0 A^3/atom
(refs/lda_setup/*.json, inputs alongside; QE 6.6, ecutwfc 100 Ry, ecutrho 1,000 Ry, Fermi-Dirac
at T_el = 6,329 K; hcp k 6x6x4, fcc k 8x8x8) measure what the mixed setup changes.
Output: out/lda_setup.json.
"""
import glob, json, os
import paths

RY_EV = 13.605693122994
DS_M_KB = 0.80
KB_MEV = 8.617333262e-2          # meV/K


def main():
    d = os.path.join(os.path.dirname(paths.DFT_REFS), "lda_setup")
    R = {r["tag"]: r for r in (json.load(open(f)) for f in glob.glob(os.path.join(d, "*.json")))}
    e = lambda t: R[t]["E_Ry"] * RY_EV * 1e3 / R[t]["nat"]       # meV/atom
    dfh_pz = e("fcc_pz") - e("hcp_pz")
    dfh_mix = e("fcc_mix") - e("hcp_mix")
    out = dict(
        hcp_P_pz=R["hcp_pz"]["P_GPa"], hcp_P_mix=R["hcp_mix"]["P_GPa"], hcp_P_pbe=R["hcp_pbe"]["P_GPa"],
        hcp_P_mix_minus_pz=R["hcp_mix"]["P_GPa"] - R["hcp_pz"]["P_GPa"],
        fcc_P_mix_minus_pz=R["fcc_mix"]["P_GPa"] - R["fcc_pz"]["P_GPa"],
        hcp_P_lda_minus_pbe=R["hcp_pz"]["P_GPa"] - R["hcp_pbe"]["P_GPa"],
        fcc_minus_hcp_pz_meV=dfh_pz, fcc_minus_hcp_mix_meV=dfh_mix,
        fcc_minus_hcp_mix_minus_pz_meV=dfh_mix - dfh_pz,
        fcc_minus_hcp_mix_minus_pz_K=abs(dfh_mix - dfh_pz) / (DS_M_KB * KB_MEV))
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "lda_setup.json"), "w"), indent=1)
    print(f"hcp P: PZ file {out['hcp_P_pz']:.2f}, PBE file + PZ {out['hcp_P_mix']:.2f}, PBE {out['hcp_P_pbe']:.2f} GPa "
          f"(mixed - consistent {out['hcp_P_mix_minus_pz']:+.2f} GPa; LDA - PBE {out['hcp_P_lda_minus_pbe']:+.1f} GPa)")
    print(f"fcc P: mixed - consistent {out['fcc_P_mix_minus_pz']:+.2f} GPa")
    print(f"fcc - hcp: PZ file {dfh_pz:.2f}, mixed {dfh_mix:.2f} meV/atom; difference "
          f"{dfh_mix - dfh_pz:+.2f} meV/atom = {out['fcc_minus_hcp_mix_minus_pz_K']:.1f} K at dS_m = 0.80 k_B")


if __name__ == "__main__":
    main()
