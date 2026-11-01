#!/usr/bin/env python3
"""
phaseC_fill.py — fill the two empty coverage cells found by G-C1 (fe_melt_v1)

  hcp 150-200 GPa: no hcp sampling volume mapped into the bin (7.9 A^3/atom
    lands at 149.7 GPa). Add V = 7.65 (P ~ 176 GPa), four T fractions,
    TETHERED (both neighbouring volumes glided when untethered).
  liquid 250-300 GPa: liquid volumes 7.2 and 6.8 land at 238 and 308 GPa.
    Add V = 7.0 (P ~ 270 GPa), three T fractions, standard liquid protocol.

Usage: python phaseC_fill.py --model <path> --job N --device cuda   (N = 0..6)
"""
import argparse, json, os
import numpy as np
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as MBD
from ase.io import write
from ase.calculators.mixing import SumCalculator

from phaseC_sample import hcp_cell, P_of_V, Tm_est, GRID, run_one
from phaseC_tether import Tether, K_TETHER, DT_FS, T_EQ_PS, T_PROD_PS, FRICTION

JOBS = [("hcp_t", 7.65, 0.60), ("hcp_t", 7.65, 0.80),
        ("hcp_t", 7.65, 0.95), ("hcp_t", 7.65, 1.05),
        ("liquid", 7.0, 1.00), ("liquid", 7.0, 1.15), ("liquid", 7.0, 1.35)]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--job", type=int, required=True)
    ap.add_argument("--device", default="cpu")
    a = ap.parse_args()
    kind, vat, tfrac = JOBS[a.job]
    from mace.calculators import mace_mp
    calc = mace_mp(model=os.path.expanduser(a.model),
                   default_dtype="float64", device=a.device)
    os.makedirs("traj", exist_ok=True)
    if kind == "liquid":
        run_one(calc, "liquid", vat, tfrac, "traj")
        return
    # tethered hcp, same protocol as phaseC_tether
    P = P_of_V(vat); T = tfrac * Tm_est(P)
    tag = f"hcp_v{vat:.2f}_T{int(T)}"
    fx, fj = f"traj/{tag}.extxyz", f"traj/{tag}.json"
    for f in (fx, fj):
        if os.path.exists(f): os.remove(f)
    at = hcp_cell(vat)
    at.calc = SumCalculator([calc, Tether(at.get_positions(), K_TETHER)])
    MBD(at, temperature_K=T)
    dyn = Langevin(at, DT_FS * units.fs, temperature_K=T, friction=FRICTION)
    dyn.run(int(T_EQ_PS * 1000 / DT_FS))
    r0 = at.get_positions().copy()
    nsteps = int(T_PROD_PS * 1000 / DT_FS); stride = max(1, nsteps // 32)
    msd, kept = [], 0
    for i in range(nsteps):
        dyn.run(1)
        if i % stride == 0:
            d = at.get_positions() - r0
            msd.append(float((d * d).sum(axis=1).mean()))
            write(fx, at.copy(), append=True); kept += 1
    m = float(np.mean(msd[-5:]))
    json.dump(dict(tag=tag, phase="hcp", vat=vat, P_est=P, T=T, frames=kept,
                   msd_final=m, FLAG_decay=bool(m > 0.6),
                   note=f"fill run, tethered k={K_TETHER}"), open(fj, "w"), indent=1)
    print(f"{tag}: {kept} frames, MSD_end={m:.3f}")

if __name__ == "__main__":
    main()
