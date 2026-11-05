"""melting_point.py - T_m(330 GPa) and dT_m/dP from the Simon fit of the T_el-conditioned campaign.

Input: data/coexistence/pn_correct.json, runs_tel (pn_correct.py): per-cell P_n (interface-normal pressure),
T_m, standard error, and the Simon fit T_m(P) = T1 (1 + (P - 100 GPa)/a)^(1/c).
The fit is repeated here (same weights as pn_correct.py: sigma = max(sem, 10 K), absolute) to check
it reproduces the stored parameters and to give its covariance, which pn_correct.json does not store.
Output: out/melting_point.json - T_m and dT_m/dP at 330 GPa (and 323, 328, 360, 364 GPa), the fit,
the fit's standard error of T_m(330 GPa), and the per-cell P_n and T_m used.
"""
import json, os, sys
import numpy as np
from scipy.optimize import curve_fit
import paths

sys.path.insert(0, os.path.dirname(paths.HERE))
import layout  # noqa: E402

P_REF = 100.0


def simon(P, T1, a, c):
    return T1 * (1.0 + (np.asarray(P, float) - P_REF) / a) ** (1.0 / c)


def main():
    pn = json.load(open(layout.PN_CORRECT))["runs_tel"]
    f = pn["fit"]
    p0 = np.array([f["T1"], f["a"], f["c"]])
    h = sorted((v["P"], v["Tm"], max(v.get("sem", 10.0), 10.0)) for v in pn["cells"].values())
    P, T, s = map(np.array, zip(*h))
    p, cov = curve_fit(simon, P, T, p0=p0, sigma=s, absolute_sigma=True, maxfev=20000)
    assert np.allclose(p, p0, rtol=1e-6), (p, p0)
    chi2 = float((((T - simon(P, *p0)) / s) ** 2).sum() / (len(P) - 3))

    def at(Px):
        g = np.array([(simon(Px, *(p0 + d)) - simon(Px, *(p0 - d))) / (2 * d[i])
                      for i, d in enumerate(np.diag(1e-6 * np.abs(p0)))])
        return dict(Tm=float(simon(Px, *p0)), dTdP=float((simon(Px + 0.01, *p0) - simon(Px - 0.01, *p0)) / 0.02),
                    Tm_fit_se=float(np.sqrt(g @ cov @ g)))
    out = dict(fit=f, fit_param_se=dict(zip(("T1", "a", "c"), np.sqrt(np.diag(cov)).tolist())), fit_chi2_dof=chi2,
               at={f"{Px:.0f}": at(Px) for Px in (100.0, 200.0, 323.0, 328.0, 330.0, 360.0, 364.0)},
               cells={c: dict(P_n=v["P"], P_mean=v["P_mean"], dP=v["dP"], Tm=v["Tm"], sem=v["sem"])
                      for c, v in pn["cells"].items()})
    out["Tm_330"] = out["at"]["330"]["Tm"]
    out["dTdP_330"] = out["at"]["330"]["dTdP"]
    print(f"T_m(330 GPa) = {out['Tm_330']:.1f} K (fit s.e. {out['at']['330']['Tm_fit_se']:.1f} K), "
          f"dT_m/dP = {out['dTdP_330']:.3f} K/GPa; chi2/dof {chi2:.2f}")
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "melting_point.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
