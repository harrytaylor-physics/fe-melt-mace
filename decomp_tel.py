#!/usr/bin/env python3
"""
decomp_tel.py — split the Tel-model minus Tel-6000 melting-point difference (A29)
into an electronic-temperature part and a model part, first order, on the
coexistence configurations themselves.

Physics. For a coexistence cell with solid fraction x, the total energy per atom
under a change of potential (or of T_el) changes by
    d(E/N) = de_l + x (de_s - de_l) + (interface term, same for every run of a cell)
so a straight-line fit of d(E/N) against x over the runs of one cell gives the slope
de_s - de_l. First-order Gibbs-energy perturbation then moves the melting point by
    dTm = (de_l - de_s) / dS_m = -slope / dS_m ,   dS_m = 0.84 kB = 0.0724 meV/K.
Two differences, on the same frames:
    A (T_el part):  Tel model at T_el = T_el used  minus  Tel model at T_el = 6000 K
    B (model part): Tel model at T_el = 6000 K     minus  Tel-6000 potential
A + B is the first-order prediction of the observed matched-pressure difference.
Frames: the final checkpoint of every hcp run of both campaigns (12 per cell),
solid fraction from the run's own q6 classification.

Stages (each GPU stage in its own environment):
  python  decomp_tel.py --stage tel --device cuda   # in ~/tel_env  (Tel model, two T_el)
  python3.13 decomp_tel.py --stage ref --device cuda  # production env (Tel-6000 potential)
  python3.13 decomp_tel.py --stage combine           # any node, no torch
"""
import argparse, glob, json, os
import numpy as np

CELLS = ["coex_hcp_100GPa", "coex_hcp_150GPa", "coex_hcp_200GPa", "coex_hcp_250GPa",
         "coex_hcp_300GPa", "coex_hcp_330GPa", "coex_hcp_350GPa"]
RUNDIRS = ["phaseE/runs_tel", "phaseE/runs"]
TEL_MODEL = "phaseD/models/fe_tel_e3_stagetwo.model"
REF_MODEL = "phaseD/models/fe_ft_Tel6000_f32.model"
DS_MEV_PER_K = 0.84 * 0.08617333            # 0.0724 meV/K
GPA = 1.0 / 160.21766208                     # eV/A^3 per GPa
OUT = "phaseE/decomp_tel_{}.json"


def frames():
    """(cell, source, tag, x, ckpt path) for every hcp run with 0.05 < x < 0.95."""
    out = []
    for c in CELLS:
        for rd in RUNDIRS:
            for f in sorted(glob.glob(os.path.join(rd, c + "_ti*_ckpt.npz"))):
                z = np.load(f, allow_pickle=True)
                fr = z["fr_t"]
                if len(fr) == 0:
                    continue
                x = float(np.mean([r[1] for r in fr[-4:]]))
                if 0.05 < x < 0.95:
                    out.append((c, os.path.basename(rd), os.path.basename(f)[:-9], x, f))
    return out


def load_atoms(cell, ckpt):
    from ase.io import read
    at = read(os.path.join("phaseE/cells", cell + ".extxyz"))
    z = np.load(ckpt, allow_pickle=True)
    at.set_cell(z["cell"]); at.set_positions(z["positions"]); at.pbc = True
    return at


def energy(at, calc, tel=None):
    if tel is not None:
        at.info["elec_temp"] = float(tel)
    calc.reset()                               # info is not part of ASE's change check
    at.calc = calc
    return at.get_potential_energy() / len(at)  # eV/atom


def stage_tel(a):
    from mace.calculators import MACECalculator
    calc = MACECalculator(model_paths=TEL_MODEL, device=a.device, default_dtype="float32", enable_cueq=True)
    tel_used = {c: round(v["Tm"]) for c, v in json.load(open("phaseE/analysis_live_all.json"))["cells"].items()}
    res = []
    for c, src, tag, x, f in frames():
        at = load_atoms(c, f)
        e_used = energy(at, calc, tel_used[c])
        sv = at.get_stress(voigt=False) / GPA          # virial stress (eV/A^3 -> GPa), T_el used
        n = int(np.argmax(at.cell.lengths()))          # interface normal = long axis of the slab
        p_axes = [-float(sv[i, i]) for i in range(3)]
        e_6000 = energy(at, calc, 6000.0)
        if abs(tel_used[c] - 6000) > 200 and abs(e_used - e_6000) < 1e-7:
            raise SystemExit(f"{tag}: energy does not change with elec_temp; the model is not reading atoms.info")
        res.append(dict(cell=c, src=src, tag=tag, x=x, tel_used=tel_used[c], e_used=e_used, e_6000=e_6000,
                        p_axes=p_axes, normal=n))
        print(f"{tag:32s} {src:9s} x={x:.3f}  E(Tel {tel_used[c]})-E(6000) = {1000*(e_used-e_6000):+8.3f} meV/atom", flush=True)
    json.dump(res, open(OUT.format("tel"), "w"), indent=1)
    print("written", OUT.format("tel"))


def stage_ref(a):
    from mace.calculators import MACECalculator
    calc = MACECalculator(model_paths=REF_MODEL, device=a.device, default_dtype="float32", enable_cueq=True)
    res = []
    for c, src, tag, x, f in frames():
        at = load_atoms(c, f)
        res.append(dict(cell=c, src=src, tag=tag, x=x, e_ref=energy(at, calc)))
        print(f"{tag:32s} {src:9s} x={x:.3f}  E_ref = {res[-1]['e_ref']:.6f} eV/atom", flush=True)
    json.dump(res, open(OUT.format("ref"), "w"), indent=1)
    print("written", OUT.format("ref"))


def slope_dtm(x, d):
    """d in meV/atom; returns dTm (K) and its standard error from the slope."""
    if len(x) < 3 or np.ptp(x) < 0.1:
        return float("nan"), float("nan")
    p, cov = np.polyfit(x, d, 1, cov=True)
    return -p[0] / DS_MEV_PER_K, np.sqrt(cov[0, 0]) / DS_MEV_PER_K


def stage_combine(a):
    T = {(r["src"], r["tag"]): r for r in json.load(open(OUT.format("tel")))}
    R = {(r["src"], r["tag"]): r for r in json.load(open(OUT.format("ref")))}
    a6 = json.load(open("phaseE/analysis_live_all.json"))
    at = json.load(open("phaseE/analysis_live_all_runs_tel.json"))
    f6 = a6["fit"]
    simon = lambda P: f6["T1"] * (1 + (P - 100.0) / f6["a"]) ** (1 / f6["c"])
    print(f"First-order decomposition, dS_m = 0.84 kB; frames = final checkpoints of both campaigns\n")
    print(f"{'cell':16s} {'n':>3s} {'x range':>11s} {'A: T_el':>12s} {'B: model':>12s} {'A+B':>8s} {'observed':>9s} {'obs-(A+B)':>10s}")
    rows = []
    for c in CELLS:
        keys = [k for k in T if T[k]["cell"] == c and k in R]
        if not keys:
            continue
        x = np.array([T[k]["x"] for k in keys])
        dA = np.array([1000 * (T[k]["e_used"] - T[k]["e_6000"]) for k in keys])
        dB = np.array([1000 * (T[k]["e_6000"] - R[k]["e_ref"]) for k in keys])
        tA, sA = slope_dtm(x, dA)
        tB, sB = slope_dtm(x, dB)
        obs = at["cells"][c]["Tm"] - simon(at["cells"][c]["P"]) if c in at["cells"] else float("nan")
        rows.append(dict(cell=c, n=len(keys), dTm_tel=tA, se_tel=sA, dTm_model=tB, se_model=sB, observed=obs))
        print(f"{c:16s} {len(keys):3d} {x.min():5.2f}-{x.max():4.2f} {tA:+6.0f}+-{sA:3.0f} {tB:+6.0f}+-{sB:3.0f} "
              f"{tA + tB:+8.0f} {obs:+9.0f} {obs - tA - tB:+10.0f}")
    print("\nA: Tel model at T_el used minus at 6000 K (electronic temperature).")
    print("B: Tel model at 6000 K minus Tel-6000 potential (model change at equal T_el).")
    print("observed: Tel-model Tm minus Tel-6000 Simon fit at the Tel-model pressure (A29).")
    print("A23 DFT-bin expectation of A for comparison: 100 GPa -70 K, 150 GPa -51 K, 250 GPa +19 K, 330 GPa -60 K.")
    # Pressure convention check: the liquid carries only hydrostatic stress, so the pressure normal
    # to the interface (P_n) is the liquid pressure. The analysis assigns the cell-mean pressure
    # -tr(sigma)/3. If the solid slab carries in-plane stress, the two differ and the point moves.
    ft = at.get("fit")
    print("\nPressure convention (Tel-model frames, virial part; the kinetic part is isotropic):")
    print(f"{'cell':16s} {'n':>3s} {'P_n - P_mean (GPa)':>19s} {'curve shift at fixed P (K)':>27s}")
    for c in CELLS:
        ks = [k for k in T if T[k]["cell"] == c and "p_axes" in T[k] and T[k]["src"] == "runs_tel"]
        if not ks:
            continue
        d = np.array([T[k]["p_axes"][T[k]["normal"]] - np.mean(T[k]["p_axes"]) for k in ks])
        if ft:
            P = at["cells"][c]["P"]
            slope = ft["T1"] / (ft["c"] * ft["a"]) * (1 + (P - 100.0) / ft["a"]) ** (1 / ft["c"] - 1)
            shift = -slope * d.mean()
        else:
            shift = float("nan")
        rows.append(dict(cell=c, dP_normal_minus_mean=float(d.mean()), dP_sd=float(d.std()), curve_shift_K=float(shift)))
        print(f"{c:16s} {len(ks):3d} {d.mean():+10.2f} +- {d.std():4.2f} {shift:+20.0f}")
    print("A positive P_n - P_mean means each point belongs at a higher pressure, which lowers the curve at fixed P.")
    json.dump(rows, open(OUT.format("result"), "w"), indent=1)
    print("written", OUT.format("result"))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["tel", "ref", "combine"], required=True)
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()
    {"tel": stage_tel, "ref": stage_ref, "combine": stage_combine}[a.stage](a)
