#!/usr/bin/env python3
"""
phaseE_reweight.py — self-consistent electronic-temperature correction to the
Tel-6000 melting curve, using the three fixed-Tel potentials (A21, step 3).

PHYSICS
  Tm(6000) is measured. Changing the electronic temperature to Tel changes the
  Gibbs energies of solid and liquid by dG_sol and dG_liq. To first order in
  the perturbation, dG = <F_Tel(R) - F_6000(R)> over configurations sampled
  with the Tel-6000 potential at (Tm, P). The melting point moves by
      dTm(Tel) = (dG_liq - dG_sol) / dS_m
  (a liquid made relatively more stable melts lower). dS_m per atom is the
  entropy of fusion from the Phase E transients. Doing this at Tel = 4000 and
  8000 K gives Tm at three electronic temperatures per pressure; the
  self-consistent point solves Tm(Tel) = Tel.

PER CELL (array task):
  1. bulk hcp, orthorhombic 4-atom cell repeated (10, 6, 6) = 1,440 atoms,
     lattice a from the coexistence cell, c/a = 1.60; NPT at the measured
     (Tm, P): 5 ps equilibration, 10 ps sampling.
  2. bulk liquid: same cell melted at 1.6 Tm (3 ps NVT), then the same NPT
     protocol at (Tm, P).
  3. every 0.25 ps, energies of the snapshot with the Tel-4000, Tel-6000 and
     Tel-8000 potentials. Reports <dU> per phase, its standard error, the
     second-order (variance) correction, and the solid fraction at the end
     (checks the solid stayed solid and the liquid stayed liquid).

USAGE
  task:    python3.13 phaseE_reweight.py --cell coex_hcp_300GPa
  collect: python3.13 phaseE_reweight.py --collect
"""
import argparse, glob, json, os
import numpy as np

KB = 8.617333262e-5          # eV/K
MODELS = {T: f"phaseD/models/fe_ft_Tel{T}_f32.model" for T in (4000, 6000, 8000)}


def ortho_hcp(a, coa, reps):
    from ase import Atoms
    c = coa * a
    b = a * np.sqrt(3.0)
    base = np.array([[0, 0, 0], [a / 2, b / 2, 0], [a / 2, b / 6, c / 2], [0, 2 * b / 3, c / 2]])
    at = Atoms("Fe4", positions=base, cell=[a, b, c], pbc=True)
    return at.repeat(reps)


def run_cell(cellname, device):
    from ase import units
    from ase.io import read
    from ase.md.npt import NPT
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as MBD
    from mace.calculators import MACECalculator
    from phaseE_q6 import q6_per_atom

    res = json.load(open("phaseE/analysis_live_all.json"))["cells"][cellname]
    Tm, P = res["Tm"], res["P"]
    meta = json.load(open(f"phaseE/cells/{cellname}.json"))
    thr = meta["q6_threshold"]
    coex = read(f"phaseE/cells/{cellname}.extxyz")
    a0 = coex.cell.lengths()[0] / 12.0

    calcs = {T: MACECalculator(model_paths=p, device=device, default_dtype="float32",
                               enable_cueq=True) for T, p in MODELS.items()}
    out = dict(cell=cellname, Tm=Tm, P=P, a0=a0)
    pfactor = (750 * units.fs) ** 2 * 1000 * units.GPa

    for phase in ("solid", "liquid"):
        at = ortho_hcp(a0, 1.60, (10, 6, 6))
        at.calc = calcs[6000]
        rng = np.random.RandomState(11 if phase == "solid" else 12)
        if phase == "liquid":
            MBD(at, temperature_K=1.6 * Tm, rng=rng)
            Langevin(at, 2 * units.fs, temperature_K=1.6 * Tm, friction=0.02).run(1500)
        else:
            MBD(at, temperature_K=Tm, rng=rng)
        dyn = NPT(at, 2 * units.fs, temperature_K=Tm, externalstress=P * units.GPa,
                  ttime=25 * units.fs, pfactor=pfactor)
        dyn.run(2500)                                   # 5 ps equilibration
        dU = {4000: [], 8000: []}
        vols, Ts = [], []
        for _ in range(40):                             # 40 x 0.25 ps = 10 ps
            dyn.run(125)
            n = len(at)
            e6 = at.get_potential_energy() / n
            for T in (4000, 8000):
                b = at.copy(); b.calc = calcs[T]
                dU[T].append(b.get_potential_energy() / n - e6)
            vols.append(at.get_volume() / n); Ts.append(at.get_temperature())
        fs = float((q6_per_atom(at) > thr).mean())
        ph = dict(V=float(np.mean(vols)), T=float(np.mean(Ts)), solid_frac=fs)
        for T in (4000, 8000):
            x = np.array(dU[T])
            blocks = [b.mean() for b in np.array_split(x, 4)]
            ph[f"dU_{T}"] = float(1000 * x.mean())                       # meV/atom
            ph[f"dU_{T}_sem"] = float(1000 * np.std(blocks, ddof=1) / 2.0)
            ph[f"dU_{T}_var2"] = float(1000 * len(at) * x.var() / (2 * KB * Tm))  # meV/atom
        out[phase] = ph
        print(f"[{cellname}] {phase}: V={ph['V']:.3f} A^3/atom T={ph['T']:.0f} K "
              f"solid_frac={fs:.2f} dU4000={ph['dU_4000']:+.2f} dU8000={ph['dU_8000']:+.2f} meV/atom",
              flush=True)
    os.makedirs("phaseE/reweight", exist_ok=True)
    json.dump(out, open(f"phaseE/reweight/{cellname}.json", "w"), indent=1)


def collect():
    cells = json.load(open("phaseE/analysis_live_all.json"))["cells"]
    rows = []
    print(f"{'cell':16s} {'Tm6000':>7s} {'dTm@4000':>9s} {'dTm@8000':>9s} {'Tm_sc':>7s}  checks")
    for f in sorted(glob.glob("phaseE/reweight/coex_hcp_*.json"),
                    key=lambda p: int(p.split("_")[-1][:-8])):
        d = json.load(open(f))
        c = d["cell"]
        dS = cells[c].get("dS_kB", 0.84) * KB * 1000          # meV/atom/K
        Tm6 = d["Tm"]
        pts = {6000: Tm6}
        errs = {}
        for T in (4000, 8000):
            dd = d["liquid"][f"dU_{T}"] - d["solid"][f"dU_{T}"]
            se = np.hypot(d["liquid"][f"dU_{T}_sem"], d["solid"][f"dU_{T}_sem"])
            sec = d["liquid"][f"dU_{T}_var2"] - d["solid"][f"dU_{T}_var2"]
            pts[T] = Tm6 + dd / dS
            # second-order cumulant: dG2 = -(beta/2) N Var(du) per atom, so a larger
            # variance in the liquid lowers its G and lowers Tm: the sign is negative
            errs[T] = (se / dS, -sec / dS)
        Tels = np.array([4000.0, 6000.0, 8000.0])
        Tms = np.array([pts[4000], pts[6000], pts[8000]])
        p = np.polyfit(Tels, Tms, 2)
        roots = np.roots(p - np.array([0, 1, 0]))
        roots = [r.real for r in roots if abs(r.imag) < 1e-6 and 3000 < r.real < 9000]
        Tsc = float(roots[0]) if roots else float("nan")
        chk = []
        if d["solid"]["solid_frac"] < 0.8: chk.append("solid partly melted")
        if d["liquid"]["solid_frac"] > 0.2: chk.append("liquid partly crystallised")
        if max(abs(errs[T][1]) for T in errs) > 20: chk.append("2nd-order term >20 K")
        rows.append((cells[c]["P"], Tsc))
        print(f"{c:16s} {Tm6:7.0f} {pts[4000]-Tm6:+9.0f} {pts[8000]-Tm6:+9.0f} {Tsc:7.0f}  "
              + ("; ".join(chk) if chk else "ok")
              + f"   (stat +-{errs[4000][0]:.0f}/{errs[8000][0]:.0f} K, 2nd-order {errs[4000][1]:+.0f}/{errs[8000][1]:+.0f} K)")
    json.dump([dict(P=P, Tm_selfconsistent=T) for P, T in rows],
              open("phaseE/reweight/selfconsistent.json", "w"), indent=1)
    print("\nwritten phaseE/reweight/selfconsistent.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell")
    ap.add_argument("--collect", action="store_true")
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()
    collect() if a.collect else run_cell(a.cell, a.device)
