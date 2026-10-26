#!/usr/bin/env python3
"""
phaseE_build.py — construct solid-liquid coexistence cells (fe_melt_v1, Phase E)

Protocol (per pressure, per solid branch):
 1. Equilibrate bulk solid in NPT at (P, 0.95 Tm_est) to get the lattice
    parameter and c/a at temperature; build an elongated supercell
    (~12x12x42 A cross-section, z-long, >=12,000 atoms).
 2. Freeze the lower half; overheat and melt the upper half at 2 Tm_est
    with fixed cell; re-equilibrate the liquid at Tm_est.
 3. Release, brief full-cell NVT anneal at Tm_est to heal the interfaces.
 4. Calibrate the per-atom order parameter (mean |q6| via Steinhardt-like
    local bond order) on bulk solid and bulk liquid at the same (V,T);
    store thresholds with the cell.
Outputs: phaseE/cells/coex_{branch}_{P}GPa.extxyz + .json (thresholds, box info)

Usage:  python phaseE_build.py --model <model> --pressure 300 --branch hcp
        [--natoms 12096] [--device cuda]
Branches: hcp, fcc (fcc only needed at 100 GPa), bcc available for checks.
"""
from phaseE_nlcache import enable_nl_cache
enable_nl_cache(skin=0.0)
import argparse, json, os
import numpy as np
import torch
from ase import Atoms, units
from ase.io import write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as MBD
from ase.units import GPa

def Tm_est(P): return 3010.0 + 12.07 * P - 0.005410 * P * P

def hcp_cell(a, coa, n):
    cell = np.array([[a,0,0],[-a/2,a*np.sqrt(3)/2,0],[0,0,coa*a]])
    at = Atoms("Fe2", cell=cell, pbc=True,
               scaled_positions=[[0,0,0],[1/3.,2/3.,.5]])
    return at.repeat(n)

def fcc_cell(a, n):
    at = Atoms("Fe4", cell=np.eye(3)*a, pbc=True,
               scaled_positions=[[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])
    return at.repeat(n)

def q6_per_atom(at, rcut=3.2, device="cpu"):
    """Mean-|q6m| local bond order per atom, neighbour-list based (linear
    scaling; the previous N x N formulation exhausted GPU memory at 12k atoms).
    Pure numpy/scipy on the host: called every few ps, not every step."""
    import numpy as np
    from ase.neighborlist import neighbor_list
    from scipy.special import sph_harm
    n = len(at)
    i, D, r = neighbor_list("iDd", at, rcut)
    if len(i) == 0:
        return np.zeros(n)
    order = np.argsort(i, kind="stable")
    i, D, r = i[order], D[order], r[order]
    starts = np.searchsorted(i, np.arange(n), side="left")
    ends = np.searchsorted(i, np.arange(n), side="right")
    out = np.zeros(n)
    for a in range(n):
        s, e = starts[a], ends[a]
        if e - s < 4:
            continue
        v, rr = D[s:e], r[s:e]
        th = np.arccos(np.clip(v[:, 2] / rr, -1.0, 1.0))
        ph = np.arctan2(v[:, 1], v[:, 0])
        qs = 0.0
        for m in range(-6, 7):
            y = sph_harm(m, 6, ph, th).mean()
            qs += abs(y) ** 2
        out[a] = np.sqrt(4 * np.pi / 13 * qs)
    return out

def npt_lattice(model_calc, builder, P, T, steps=3000):
    from ase.md.nptberendsen import NPTBerendsen
    at = builder()
    at.calc = model_calc
    MBD(at, temperature_K=T)
    dyn = NPTBerendsen(at, 2*units.fs, temperature_K=T, taut=100*units.fs,
                       pressure_au=P*GPa, taup=1000*units.fs, compressibility_au=1/(250*GPa))
    dyn.run(steps)
    return at

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--pressure", type=float, required=True)
    ap.add_argument("--branch", default="hcp", choices=["hcp","fcc"])
    ap.add_argument("--natoms", type=int, default=12096)
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()
    os.makedirs("phaseE/cells", exist_ok=True)
    from mace.calculators import MACECalculator
    calc = MACECalculator(enable_cueq=True, model_paths=os.path.expanduser(a.model),
                          device=a.device, default_dtype="float32")
    P, T0 = a.pressure, Tm_est(a.pressure)

    # 1: lattice at temperature (small cell NPT)
    if a.branch == "hcp":
        probe = npt_lattice(calc, lambda: hcp_cell(2.35, 1.60, (6,6,4)), P, 0.95*T0)
        L = probe.cell.lengths(); a_lat = L[0]/6; coa = (L[2]/4)/a_lat
        nz = int(round(a.natoms / (2*12*12*1.0)))            # 12x12xnz primitive
        solid = hcp_cell(a_lat, coa, (12,12,nz))
    else:
        probe = npt_lattice(calc, lambda: fcc_cell(3.7, (5,5,5)), P, 0.95*T0)
        a_lat = probe.cell.lengths()[0]/5
        nz = int(round(a.natoms / (4*8*8)))
        solid = fcc_cell(a_lat, (8,8,nz))
    print(f"[build] lattice a={a_lat:.4f}  N={len(solid)}")

    # 2: melt upper half with lower half frozen
    solid.calc = calc
    z = solid.get_positions()[:,2]; zmid = z.mean()
    lower = z < zmid
    from ase.constraints import FixAtoms
    solid.set_constraint(FixAtoms(mask=lower))
    MBD(solid, temperature_K=2.2*T0)
    Langevin(solid, 2*units.fs, temperature_K=2.2*T0, friction=0.05).run(4000)
    Langevin(solid, 2*units.fs, temperature_K=T0, friction=0.02).run(3000)
    solid.set_constraint()
    # 3: heal interfaces
    MBD(solid, temperature_K=T0)
    Langevin(solid, 2*units.fs, temperature_K=T0, friction=0.02).run(2000)

    # 4: order-parameter thresholds from small bulk references at same a,T
    ref_s = hcp_cell(a_lat, coa, (5,5,4)) if a.branch=="hcp" else fcc_cell(a_lat,(4,4,4))
    ref_s.calc = calc; MBD(ref_s, temperature_K=T0)
    Langevin(ref_s, 2*units.fs, temperature_K=T0, friction=0.02).run(1500)
    q_s = q6_per_atom(ref_s)
    ref_l = ref_s.copy(); ref_l.calc = calc
    MBD(ref_l, temperature_K=2.5*T0)
    Langevin(ref_l, 2*units.fs, temperature_K=2.5*T0, friction=0.05).run(3000)
    Langevin(ref_l, 2*units.fs, temperature_K=T0, friction=0.02).run(1500)
    q_l = q6_per_atom(ref_l)
    thr = 0.5*(np.median(q_s)+np.median(q_l))
    tag = f"coex_{a.branch}_{int(P)}GPa"
    _o = solid.copy(); _o.calc = None; write(f"phaseE/cells/{tag}.extxyz", _o)
    json.dump(dict(tag=tag, P=P, branch=a.branch, N=len(solid), a=a_lat,
                   Tm_est=T0, q6_solid=float(np.median(q_s)),
                   q6_liquid=float(np.median(q_l)), q6_threshold=float(thr)),
              open(f"phaseE/cells/{tag}.json","w"), indent=1)
    print(f"[build] {tag}: q6 solid {np.median(q_s):.3f} / liquid {np.median(q_l):.3f} -> threshold {thr:.3f}")

from phaseE_q6 import q6_per_atom
if __name__ == "__main__":
    main()
