"""pbe_vs_exp.py - PBE pressure of hcp Fe against three 300 K experimental equations of state.

PBE: static, non-magnetic hcp, c/a 1.60, QE 6.6, ecutwfc 100 Ry, ecutrho 1,000 Ry,
Marzari-Vanderbilt smearing 0.01 Ry, k 16x16x10 (refs/eos_pbe_cold/*.json, QE 'P=' line).
Vibrational pressure (zero point + 300 K) added to PBE: P_vib(V) = 2.5 + (10.4 - V)/3.4 GPa,
a Debye-Grueneisen estimate (theta ~ 850 K, gamma ~ 1.5), uncertain by about 1 GPa.

Experiment (300 K, hcp):
  Dewaele et al. 2006, PRL 97, 215504: Vinet, V0 = 75.7 bohr^3 per atom, K0 = 163.4 GPa, K0' = 5.38
  Fei et al. 2016, GRL 43, 6837:       BM3, K0 = 172.7 GPa, K0' = 4.79, V0 = 11.214 A^3 (assumed)
  Hirao et al. 2022, MRE 7, 038403:    Vinet, V0 = 11.215 A^3, K0 = 159.27 GPa, K0' = 5.86 (Pt scale)

The T_m consequence uses dT_m/dP at 330 GPa from the Simon fit (pn_correct.json, if present;
8.41 K/GPa when this was written). A functional that is too stiff by dP = P_PBE - P_exp at the
core-density volume runs a coexistence cell at nominal 330 GPa whose true pressure is 330 - dP,
so the true T_m(330 GPa) is higher than the PBE value by (dT_m/dP) * dP:
    T_m(true) - T_m(PBE) = +(dT_m/dP) * (P_PBE - P_exp).
Output: out/pbe_vs_exp.json.
"""
import glob, json, os
import numpy as np
import paths

B3 = 0.148184711                     # A^3 per bohr^3
V_CORE = 7.0                         # A^3/atom, closest computed volume to the 330 GPa coexistence volumes


def p_vib(V):
    return 2.5 + (10.4 - V) / 3.4


def vinet(V, V0, K0, Kp):
    x = (V / V0) ** (1 / 3)
    return 3 * K0 * (1 - x) / x ** 2 * np.exp(1.5 * (Kp - 1) * (1 - x))


def bm3(V, V0, K0, Kp):
    f = (V0 / V) ** (2 / 3)
    return 1.5 * K0 * (f ** 3.5 - f ** 2.5) * (1 + 0.75 * (Kp - 4) * (f - 1))


EXP = {"Dewaele2006": lambda V: vinet(V, 75.7 * B3, 163.4, 5.38),
       "Fei2016": lambda V: bm3(V, 11.214, 172.7, 4.79),
       "Hirao2022": lambda V: vinet(V, 11.215, 159.27, 5.86)}


def dtdp_330():
    f = os.path.join(paths.PHASEE, "pn_correct.json")
    if not os.path.exists(f):
        return 8.41, "default (pn_correct.json not found)"
    s = json.load(open(f))["runs_tel"]["fit"]
    T = lambda P: s["T1"] * ((P - 100.0) / s["a"] + 1) ** (1 / s["c"])
    return float((T(330.01) - T(329.99)) / 0.02), "Simon fit, pn_correct.json runs_tel"


def main():
    rows = sorted((json.load(open(f)) for f in glob.glob(os.path.join(paths.EOS_COLD, "PBE_cold_V*.json"))),
                  key=lambda r: -r["V"])
    dtdp, src = dtdp_330()
    out = dict(rows=[], dtdp_330=dtdp, dtdp_source=src, V_core=V_CORE)
    print(f"{'V':>6} {'P_PBE':>7} {'P_vib':>5} " + " ".join(f"{k:>12}" for k in EXP) + "   PBE+vib - experiment (GPa)")
    for r in rows:
        V, P = r["V"], r["P"]
        pe = {k: float(f(V)) for k, f in EXP.items()}
        d = {k: P + p_vib(V) - v for k, v in pe.items()}
        out["rows"].append(dict(V=V, P_pbe_static=P, P_vib=p_vib(V), P_exp=pe, diff=d))
        print(f"{V:6.2f} {P:7.1f} {p_vib(V):5.2f} " + " ".join(f"{pe[k]:12.1f}" for k in EXP) + "   "
              + " ".join(f"{d[k]:+6.1f}" for k in EXP))
    core = [x for x in out["rows"] if abs(x["V"] - V_CORE) < 1e-6][0]["diff"]
    lo, hi = min(core.values()), max(core.values())
    out.update(core_diff=core, core_range=[lo, hi], dTm_true_minus_pbe_range=[dtdp * lo, dtdp * hi])
    print(f"\nat {V_CORE} A^3/atom: PBE - experiment {lo:+.1f} to {hi:+.1f} GPa; "
          f"T_m(true) - T_m(PBE) = +(dT_m/dP) dP with dT_m/dP = {dtdp:.2f} K/GPa ({src}): "
          f"{dtdp * lo:+.0f} to {dtdp * hi:+.0f} K")
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "pbe_vs_exp.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
