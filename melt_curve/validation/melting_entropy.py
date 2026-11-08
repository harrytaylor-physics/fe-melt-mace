"""melting_entropy.py - entropy and volume of melting at 330 GPa, 6,329 K from the direct-test snapshots.

Input: studies/lda_shift/pbe_snapshot_energies_330.json (PBE energies and static pressures of the potential
on 20 hcp and 20 liquid snapshots, with V and B of each phase from the NVT runs).
Output: out/melting_entropy.json. Energies meV/atom, pressures GPa, entropies k_B/atom.

dH = <E_int>_liq - <E_int>_hcp + P (V_liq - V_hcp) at P = 330 GPa; dS = dH / (k_B T). The error is the standard
error of the difference of the two means. P_total = <P_mace_static> + k_B T / V. Volume at 330 GPa:
V' = V (1 + (P_total - 330) / B). Volume ratios are relative to hcp.

The input's expected_results block also lists dS_corrected_kB (0.956) and clapeyron_dTdP_K_per_GPa (8.43). Neither
is written here: the formula behind them is not recorded, and V' in dH gives 1.008 and 8.00 instead.
"""
import json, os
import numpy as np
import paths

KB = 8.617333262e-5                  # eV/K
EV_A3_PER_GPA = 1.0 / 160.21766208   # eV/A^3 per GPa
SRC = getattr(paths, "LDASHIFT", os.path.join(paths.ROOT, "studies", "lda_shift"))
INPUT = os.path.join(SRC, "pbe_snapshot_energies_330.json")


def mean_se(x):
    x = np.asarray(x, float)
    return float(x.mean()), float(x.std(ddof=1) / np.sqrt(len(x)))


def at_target(V, P_total, B, P):
    """Volume at pressure P from the volume, total pressure and bulk modulus of a run."""
    return V * (1.0 + (P_total - P) / B)


def compute():
    d = json.load(open(INPUT))
    T, P = d["T_K"], d["P_target_GPa"]
    h, l = d["md_meta"]["hcp"], d["md_meta"]["liquid"]
    E_h, se_h = mean_se([s["E_int_eV_per_atom"] for s in d["snapshots"]["hcp"]])
    E_l, se_l = mean_se([s["E_int_eV_per_atom"] for s in d["snapshots"]["liquid"]])
    dE = E_l - E_h
    dH = dE + P * EV_A3_PER_GPA * (l["V_at"] - h["V_at"])
    P_h = mean_se([s["P_mace_static_GPa"] for s in d["snapshots"]["hcp"]])[0] + KB * T / h["V_at"] / EV_A3_PER_GPA
    P_l = mean_se([s["P_mace_static_GPa"] for s in d["snapshots"]["liquid"]])[0] + KB * T / l["V_at"] / EV_A3_PER_GPA
    V_h, V_l = at_target(h["V_at"], P_h, h["B_GPa"], P), at_target(l["V_at"], P_l, l["B_GPa"], P)
    return dict(
        dE_meV=1e3 * dE, dH_meV=1e3 * dH, dS_kB=dH / (KB * T), dS_err_kB=float(np.hypot(se_h, se_l)) / (KB * T),
        dV_over_V_pct=100 * (l["V_at"] - h["V_at"]) / h["V_at"], P_total_hcp_GPa=P_h, P_total_liq_GPa=P_l,
        dV_over_V_330_pct=100 * (V_l - V_h) / V_h)


def main():
    out = compute()
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "melting_entropy.json"), "w"), indent=1)
    print(f"dH = {out['dH_meV']:.0f} meV/atom, dS_m = {out['dS_kB']:.3f} +- {out['dS_err_kB']:.3f} k_B; "
          f"dV/V = {out['dV_over_V_pct']:.2f}% ({out['dV_over_V_330_pct']:.2f}% at 330 GPa)")


if __name__ == "__main__":
    main()
