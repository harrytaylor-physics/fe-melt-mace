#!/usr/bin/env python3
"""
phaseC_sample.py — training-set sampling with MACE-MatPES-PBE-0 (fe_melt_v1, Phase C)

Generates MD trajectories whose frames become QE label targets. The sampler's
only job is to visit realistic configurations; ALL physics comes from the QE
labels later. Design rules implemented here:
  - SOLIDS run at FIXED CELL (NVT): the foundation model's documented bcc bias
    (Phase B: bcc 254 meV/atom below hcp at 300 GPa) must not be allowed to
    decay hcp trajectories; fixed lattice vectors block the martensitic path.
    Internal decay is still monitored via MSD (flagged, run excluded if hit).
  - LIQUID runs free NVT at target volumes, prepared by overheat-scramble.
  - Temperatures are set per-volume as fractions of Tm_est(P) using the
    Sun/Wu Simon fit and this work's Tel=6000 K EOS for the P(V) map.

Budget (design C3): hcp 896, bcc 300, fcc 200, liquid 900 frames.

Usage (run per phase so jobs are restartable / schedulable):
  python phaseC_sample.py --model ~/fe_melt/MACE-matpes-pbe-omat-ft.model --phase hcp
  ... --phase bcc | fcc | liquid
  Optional: --device cuda   (ARC GPU nodes; cpu default)
            --run-index N   (run only the Nth (vol,T) job of that phase)

Outputs: traj/<phase>_vXX_TYYYY.extxyz (frames) + diagnostics JSON per run.
"""
import argparse, json, os
import numpy as np
from ase import Atoms, units
from ase.md.langevin import Langevin
from ase.io import write

# ---------- this work's P(V) map (Tel = 6000 K, static ions), campaigns B/C2 ----
V_MAP = np.array([6.0, 6.4, 6.8, 7.2, 7.6, 8.0, 8.5, 9.0, 9.5, 10.0])
P_MAP = np.array([517.0, 399.0, 308.2, 237.8, 182.3, 138.9, 96.8, 65.1, 41.1, 22.9])

def P_of_V(v):   return float(np.interp(v, V_MAP, P_MAP))
def Tm_est(P):   return 3010.0 + 12.07 * P - 0.005410 * P * P   # Sun/Wu Simon fit (K, GPa)

# ---------- structures (volume per atom in A^3) ---------------------------------
def hcp_cell(vat, nx=4, ny=4, nz=4, coa=1.60):
    a = (4.0 * vat / (np.sqrt(3.0) * coa)) ** (1 / 3.0); c = coa * a
    cell = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]])
    pos = [[0, 0, 0], [1 / 3., 2 / 3., 0.5]]
    at = Atoms("Fe2", cell=cell, scaled_positions=pos, pbc=True)
    return at.repeat((nx, ny, nz))                      # 128 atoms

def bcc_cell(vat, n=4):
    a = (2.0 * vat) ** (1 / 3.0)
    at = Atoms("Fe2", cell=np.eye(3) * a, pbc=True,
               scaled_positions=[[0, 0, 0], [.5, .5, .5]])
    return at.repeat((n, n, n))                         # 128 atoms

def fcc_cell(vat, n=3):
    a = (4.0 * vat) ** (1 / 3.0)
    at = Atoms("Fe4", cell=np.eye(3) * a, pbc=True,
               scaled_positions=[[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    return at.repeat((n, n, n))                         # 108 atoms

# ---------- sampling grids (design C1/C3) ---------------------------------------
GRID = {
  "hcp":    dict(vols=[6.2, 6.6, 7.0, 7.4, 7.9, 8.5, 9.2],
                 tfrac=[0.60, 0.80, 0.95, 1.05], frames=32, builder=hcp_cell),
  "bcc":    dict(vols=[6.4, 7.0, 7.9, 8.8],
                 tfrac=[0.60, 0.85, 1.00],       frames=25, builder=bcc_cell),
  "fcc":    dict(vols=[7.9, 8.5, 9.2],           # 80-130 GPa weighting (A2)
                 tfrac=[0.80, 0.95, 1.00],       frames=23, builder=fcc_cell),
  "liquid": dict(vols=[6.4, 6.8, 7.2, 7.7, 8.3, 9.0],
                 tfrac=[1.00, 1.15, 1.35],       frames=50, builder=hcp_cell),
}
DT_FS, T_EQ_PS, T_PROD_PS, FRICTION = 2.0, 2.0, 8.0, 0.02
MSD_SOLID_LIMIT = 0.6  # A^2 sustained -> flag possible internal decay / premelt

def run_one(calc, phase, vat, tfrac, outdir):
    P = P_of_V(vat); T = tfrac * Tm_est(P)
    tag = f"{phase}_v{vat:.1f}_T{int(T)}"
    fx = os.path.join(outdir, f"{tag}.extxyz")
    fj = os.path.join(outdir, f"{tag}.json")
    if os.path.exists(fj):
        print(f"skip {tag} (done)"); return
    at = GRID[phase]["builder"](vat)
    if phase == "liquid":                    # overheat-scramble preparation
        at.calc = calc
        from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as MBD
        MBD(at, temperature_K=3.0 * T)
        Langevin(at, DT_FS * units.fs, temperature_K=2.5 * T,
                 friction=FRICTION).run(1500)
    at.calc = calc
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as MBD
    MBD(at, temperature_K=T)
    dyn = Langevin(at, DT_FS * units.fs, temperature_K=T, friction=FRICTION)
    dyn.run(int(T_EQ_PS * 1000 / DT_FS))                       # discard equil
    r0 = at.get_positions().copy()
    nsteps = int(T_PROD_PS * 1000 / DT_FS)
    stride = max(1, nsteps // GRID[phase]["frames"])
    msd_series, kept = [], 0
    for i in range(nsteps):
        dyn.run(1)
        if i % stride == 0:
            d = at.get_positions() - r0
            msd_series.append(float((d * d).sum(axis=1).mean()))
            write(fx, at, append=True)
            kept += 1
    msd_end = float(np.mean(msd_series[-5:]))
    flagged = (phase != "liquid") and (msd_end > MSD_SOLID_LIMIT)
    json.dump(dict(tag=tag, phase=phase, vat=vat, P_est=P, T=T, frames=kept,
                   msd_final=msd_end, FLAG_decay=flagged), open(fj, "w"), indent=1)
    print(f"{tag}: {kept} frames, MSD_end={msd_end:.3f} A^2"
          + ("  ** FLAGGED (possible solid decay) — exclude & inspect **" if flagged else ""))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--phase", required=True, choices=list(GRID))
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--run-index", type=int, default=-1)
    ap.add_argument("--outdir", default="traj")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    from mace.calculators import mace_mp
    calc = mace_mp(model=os.path.expanduser(a.model),
                   default_dtype="float64", device=a.device)
    jobs = [(v, t) for v in GRID[a.phase]["vols"] for t in GRID[a.phase]["tfrac"]]
    if a.run_index >= 0:
        jobs = [jobs[a.run_index]]
    print(f"{a.phase}: {len(jobs)} runs")
    for v, t in jobs:
        run_one(calc, a.phase, v, t, a.outdir)

if __name__ == "__main__":
    main()
