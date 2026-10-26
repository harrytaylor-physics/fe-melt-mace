"""val_mace.py - the production potential against the static PBE references in refs/dft/*.json.

float32 as in production; elec_temp input = T_el / 1e4 (upstream mace-torch, see load_model.py).
Sections: 1 hcp EOS, 2 T_el response, 3 fcc/bcc vs hcp, 4 c/a, 5 clamped-ion stress-strain
coefficients plus EOS-fit bulk modulus, Voigt moduli and wave speeds, 6 E2g phonon.
Writes out/static_compare.json. Needs torch, mace-torch and the model file (paths.MODEL).
One 2-atom evaluation takes about 10 s on 2 CPU threads; the 34 structures take about 6 min.
"""
import glob, json, os, sys
import numpy as np
from ase import Atoms
from scipy.optimize import curve_fit
import paths

EVA3_GPA = 160.21766
AMU = 1.66053907e-27
C_CM = 2.99792458e10


TEL_SCALE = 1e4


def vinet(V, V0, K0, Kp):
    x = (V / V0) ** (1 / 3)
    return 3 * K0 * (1 - x) / x ** 2 * np.exp(1.5 * (Kp - 1) * (1 - x))


def load():
    R = {}
    for f in glob.glob(os.path.join(paths.DFT_REFS, "*.json")):
        r = json.load(open(f))
        R[r["tag"]] = r
    return R


def mace_eval(R, calc):
    for tag, r in R.items():
        at = Atoms("Fe%d" % len(r["positions"]), positions=r["positions"], cell=r["cell"], pbc=True)
        at.info["elec_temp"] = r["tel"] / TEL_SCALE
        at.calc = calc
        e = at.get_potential_energy() / len(at)
        s = at.get_stress(voigt=False) * EVA3_GPA            # physical stress, GPa (tensile > 0)
        f = at.get_forces()
        r["mace"] = dict(E=e, P=-np.trace(s) / 3, sigma=s.tolist(), forces=f.tolist())
        # DFT physical stress = - QE printed stress
        r["sigma_phys"] = (-np.array(r["sigma_qe"])).tolist()
    return R


def report(R):
    out = {}
    print("1. hcp EOS, static lattice, T_el = 6000 K (c/a 1.60)")
    print(f"{'V':>6} {'P_DFT':>8} {'P_MACE':>8} {'dP':>6} | {'dE_DFT':>8} {'dE_MACE':>8} {'diff':>6}  (energies relative to 7.00 A^3, meV/atom)")
    e0d, e0m = R["eos_V7.00"]["F_per_atom"], R["eos_V7.00"]["mace"]["E"]
    rows = []
    for tag in sorted([t for t in R if t.startswith("eos_")], key=lambda t: R[t]["V"]):
        r = R[tag]
        dP = r["mace"]["P"] - r["P"]
        de_d, de_m = 1e3 * (r["F_per_atom"] - e0d), 1e3 * (r["mace"]["E"] - e0m)
        rows.append((r["V"], r["P"], r["mace"]["P"], dP, de_d, de_m))
        print(f"{r['V']:6.2f} {r['P']:8.2f} {r['mace']['P']:8.2f} {dP:+6.2f} | {de_d:8.1f} {de_m:8.1f} {de_m - de_d:+6.1f}")
    out["eos"] = rows
    # bulk modulus K = -V dP/dV from a Vinet fit to each set of eight pressures
    out["bulk_modulus"] = {}
    V = np.array([x[0] for x in rows])
    for lab, col in (("dft", 1), ("mace", 2)):
        P = np.array([x[col] for x in rows])
        p, _ = curve_fit(vinet, V, P, p0=(10.0, 250.0, 4.5), maxfev=20000)
        K = {f"{v:.2f}": float(-v * (vinet(v + 1e-4, *p) - vinet(v - 1e-4, *p)) / 2e-4) for v in (8.0, 7.0, 6.8)}
        out["bulk_modulus"][lab] = dict(K_GPa=K, fit_rms_GPa=float(np.sqrt(np.mean((vinet(V, *p) - P) ** 2))))
    print("   bulk modulus from Vinet fit (GPa): " + ", ".join(
        f"V {v}: DFT {out['bulk_modulus']['dft']['K_GPa'][v]:.1f}, MACE {out['bulk_modulus']['mace']['K_GPa'][v]:.1f}"
        for v in ("8.00", "7.00", "6.80")))
    if all(t in R for t in ("tel_T4000", "eos_V7.00", "tel_T8000")):
        print("\n2. Electronic-temperature dependence, hcp 7.00 A^3")
        for tag in ("tel_T4000", "eos_V7.00", "tel_T8000"):
            r = R[tag]
            print(f"   T_el {r['tel']:.0f} K: P_DFT {r['P']:.2f}, P_MACE {r['mace']['P']:.2f} GPa; "
                  f"F - F(6000): DFT {1e3 * (r['F_per_atom'] - e0d):+.1f}, MACE {1e3 * (r['mace']['E'] - e0m):+.1f} meV/atom")
    print("\n3. Phase energy differences at fixed volume, T_el = 6000 K (meV/atom)")
    for V in (7.0, 8.0):
        h = f"eos_V{V:.2f}"
        for ph in ("fcc", "bcc"):
            t = f"{ph}_V{V:.2f}"
            if t in R and h in R:
                dd = 1e3 * (R[t]["F_per_atom"] - R[h]["F_per_atom"]); dm = 1e3 * (R[t]["mace"]["E"] - R[h]["mace"]["E"])
                out[f"{ph}-hcp_{V}"] = (dd, dm)
                print(f"   V {V:.2f}: {ph}-hcp DFT {dd:+.1f}, MACE {dm:+.1f}, diff {dm - dd:+.1f};  "
                      f"P {ph}: DFT {R[t]['P']:.1f}, MACE {R[t]['mace']['P']:.1f} GPa")
    if all(t in R for t in ("coa_1.56", "eos_V6.80", "coa_1.64")):
        print("\n4. c/a at 6.80 A^3, T_el = 6000 K")
        x = np.array([1.56, 1.60, 1.64])
        for lab, key in (("DFT", None), ("MACE", "mace")):
            y = [R[t]["F_per_atom"] if key is None else R[t]["mace"]["E"] for t in ("coa_1.56", "eos_V6.80", "coa_1.64")]
            c = np.polyfit(x, y, 2)
            print(f"   {lab}: c/a minimum {-c[1] / (2 * c[0]):.4f}")
    el = [t for t in R if t.startswith("el_")]
    if len(el) == 8:
        print("\n5. Clamped-ion stress-strain coefficients, hcp 7.00 A^3, T_el 6000 K (GPa)")
        eps = 0.01
        def ds(tag, src, i, j):
            a = np.array(R[tag + "p"]["sigma_phys"] if src == "dft" else R[tag + "p"]["mace"]["sigma"])
            b = np.array(R[tag + "m"]["sigma_phys"] if src == "dft" else R[tag + "m"]["mace"]["sigma"])
            return (a[i, j] - b[i, j]) / (2 * eps)
        res = {}
        for src in ("dft", "mace"):
            res[src] = dict(C11=ds("el_e1", src, 0, 0), C12=ds("el_e1", src, 1, 1), C13=ds("el_e1", src, 2, 2),
                            C33=ds("el_e3", src, 2, 2), C44=ds("el_e4", src, 1, 2), C66=ds("el_e6", src, 0, 1))
        for k in res["dft"]:
            d, m = res["dft"][k], res["mace"][k]
            print(f"   {k}: DFT {d:7.1f}  MACE {m:7.1f}  ({100 * (m - d) / d:+.1f}%)")
        out["elastic"] = res
        rho = 55.845 * AMU / 7.0e-30                                            # kg/m^3 at 7.0 A^3/atom
        out["voigt"] = {}
        for src in ("dft", "mace"):
            c = res[src]
            K = (2 * c["C11"] + c["C33"] + 2 * c["C12"] + 4 * c["C13"]) / 9
            M = c["C11"] + c["C12"] + 2 * c["C33"] - 4 * c["C13"]
            G = (M + 3 * c["C11"] - 3 * c["C12"] + 12 * c["C44"] + 12 * c["C66"]) / 30
            vp = np.sqrt((K + 4 * G / 3) * 1e9 / rho) / 1e3
            vs = np.sqrt(G * 1e9 / rho) / 1e3
            out["voigt"][src] = dict(K_V=K, G_V=G, vp_km_s=vp, vs_km_s=vs, half_C11_minus_C12=(c["C11"] - c["C12"]) / 2,
                                     density_g_cm3=rho / 1e3)
            print(f"   {src.upper()}: K_V {K:.1f}, G_V {G:.1f} GPa, v_p {vp:.2f}, v_s {vs:.2f} km/s, "
                  f"(C11-C12)/2 {(c['C11'] - c['C12']) / 2:.1f} GPa (C66 {c['C66']:.1f})")
    if "ph_p" in R and "ph_m" in R:
        print("\n6. E2g optical phonon at Gamma (frozen phonon, u = 0.03 A), hcp 7.00 A^3, T_el 6000 K")
        u = 0.03
        for lab, key in (("DFT", "forces"), ("MACE", None)):
            fp = np.array(R["ph_p"][key] if key else R["ph_p"]["mace"]["forces"])[0, 0]
            fm = np.array(R["ph_m"][key] if key else R["ph_m"]["mace"]["forces"])[0, 0]
            k = -(fp - fm) / (2 * u)                       # eV/A^2
            w = np.sqrt(k * 16.0217663 / (55.845 * AMU)) / (2 * np.pi)
            print(f"   {lab}: k {k:.2f} eV/A^2, frequency {w / 1e12:.2f} THz = {w / C_CM:.0f} cm^-1")
            out[f"E2g_{lab}"] = w / C_CM
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "static_compare.json"), "w"), indent=1, default=float)


if __name__ == "__main__":
    import load_model as lm
    TEL_SCALE = lm.TEL_SCALE
    R = mace_eval(load(), lm.calc("float32"))
    report(R)
