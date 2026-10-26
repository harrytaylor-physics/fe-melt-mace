#!/usr/bin/env python3
"""
phaseC_tether.py — rerun the 13 decayed hcp grid points with harmonic tethers
(fe_melt_v1, Phase C amendment)

Each atom is tethered to its ideal lattice site by a harmonic spring,
k = 2.0 eV/A^2. Effect on the two relevant motions:
  - thermal vibration (u ~ 0.2-0.3 A): tether force ~0.5 eV/A against thermal
    forces of ~3 eV/A; the sampled vibrational ensemble is mildly narrowed
    (quote as a sampling note, not a physics error - labels are exact DFT
    for whatever geometry is visited).
  - basal glide (u ~ 1.4 A per glided atom): tether energy ~2 eV/atom, far
    above kT; the stacking transformation is suppressed by construction.

The 13 (volume, T-fraction) points are hard-coded from the flagged runs of
the first sampling campaign. Frames are written with the same tags (fresh
files), so phaseC_select.py picks them up directly.

Usage (one point per task via --job 0..12):
  python phaseC_tether.py --model <path> --job N --device cuda
"""
import argparse, json, os
import numpy as np
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as MBD
from ase.io import write
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.mixing import SumCalculator

from phaseC_sample import hcp_cell, P_of_V, Tm_est   # reuse builders/maps

JOBS = [  # (vat, tfrac) of the 13 flagged runs
 (6.6,0.80),(6.6,0.95),(6.6,1.05),
 (7.0,0.60),(7.0,1.05),
 (7.4,0.60),(7.4,0.80),(7.4,0.95),(7.4,1.05),
 (7.9,0.95),(7.9,1.05),
 (8.5,0.80),
 (9.2,0.95),
]
K_TETHER = 2.0     # eV/A^2
DT_FS, T_EQ_PS, T_PROD_PS, FRICTION, NFRAMES = 2.0, 2.0, 8.0, 0.02, 32

class Tether(Calculator):
    implemented_properties = ["energy", "forces"]
    def __init__(self, r_ref, k):
        super().__init__(); self.r_ref = r_ref.copy(); self.k = k
    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        d = atoms.get_positions() - self.r_ref
        self.results["energy"] = 0.5 * self.k * float((d * d).sum())
        self.results["forces"] = -self.k * d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--job", type=int, required=True)
    ap.add_argument("--device", default="cpu")
    a = ap.parse_args()
    vat, tfrac = JOBS[a.job]
    P = P_of_V(vat); T = tfrac * Tm_est(P)
    tag = f"hcp_v{vat:.1f}_T{int(T)}"
    os.makedirs("traj", exist_ok=True)
    fx, fj = f"traj/{tag}.extxyz", f"traj/{tag}.json"
    for f in (fx, fj):
        if os.path.exists(f): os.remove(f)
    from mace.calculators import mace_mp
    mace = mace_mp(model=os.path.expanduser(a.model),
                   default_dtype="float64", device=a.device)
    at = hcp_cell(vat)
    at.calc = SumCalculator([mace, Tether(at.get_positions(), K_TETHER)])
    MBD(at, temperature_K=T)
    dyn = Langevin(at, DT_FS * units.fs, temperature_K=T, friction=FRICTION)
    dyn.run(int(T_EQ_PS * 1000 / DT_FS))
    r0 = at.get_positions().copy()
    nsteps = int(T_PROD_PS * 1000 / DT_FS)
    stride = max(1, nsteps // NFRAMES)
    msd, kept = [], 0
    for i in range(nsteps):
        dyn.run(1)
        if i % stride == 0:
            d = at.get_positions() - r0
            msd.append(float((d * d).sum(axis=1).mean()))
            fr = at.copy()           # strip the SumCalculator before writing
            write(fx, fr, append=True); kept += 1
    m_end = float(np.mean(msd[-5:]))
    json.dump(dict(tag=tag, phase="hcp", vat=vat, P_est=P, T=T, frames=kept,
                   msd_final=m_end, FLAG_decay=bool(m_end > 0.6),
                   note=f"tethered k={K_TETHER} eV/A^2"),
              open(fj, "w"), indent=1)
    print(f"{tag}: {kept} frames, MSD_end={m_end:.3f} (tethered)")

if __name__ == "__main__":
    main()
