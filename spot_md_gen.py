#!/usr/bin/env python3
"""
spot_md_gen.py — fresh spot-check frames for the Tel-model test (A21/A22).

The held-out Tel-6000 test frames all appear in some training set at another
electronic temperature (independent splits per Tel), so none are unseen
geometries. This script makes new ones, at exactly the conditions that matter:
128-atom hcp and liquid cells run with the Tel-6000 potential at the measured
melting point (Tm, P) of five pressures, two snapshots per phase per pressure,
2 ps apart. Each snapshot's phase is confirmed with the cell's q-bar-6
threshold (solid fraction > 0.8 for hcp, < 0.2 for liquid); failing snapshots
are skipped and the run continues until two pass.

Then writes pw.x inputs at Tel = 4600 and 6500 K by copying an existing label
input of the same phase (identical settings; only cell, positions, stem and
degauss change), labels/spot_list.txt, and phaseD/spot_frames.xyz.

Run on a GPU node, production env:  python3.13 spot_md_gen.py
"""
import glob, json, os, re
import numpy as np
from ase import Atoms, units
from ase.io import read, write
from ase.md.npt import NPT
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as MBD
from mace.calculators import MACECalculator
from phaseE_q6 import q6_per_atom

CELLS = ["coex_hcp_100GPa", "coex_hcp_150GPa", "coex_hcp_250GPa",
         "coex_hcp_300GPa", "coex_hcp_350GPa"]
TELS = (4600, 6500)
KB_RY = 8.617333262e-5 / 13.605693123


def ortho_hcp(a, coa, reps):
    c, b = coa * a, a * np.sqrt(3.0)
    base = np.array([[0, 0, 0], [a / 2, b / 2, 0], [a / 2, b / 6, c / 2], [0, 2 * b / 3, c / 2]])
    return Atoms("Fe4", positions=base, cell=[a, b, c], pbc=True).repeat(reps)


def make_input(tpl_path, at, stem, tel):
    src = open(tpl_path).read()
    nat_old = int(re.search(r"nat\s*=\s*(\d+)", src).group(1))
    txt = src.replace(os.path.basename(tpl_path)[:-3], stem)
    txt, n = re.subn(r"degauss\s*=\s*[0-9.eEdD+-]+", f"degauss = {KB_RY * tel:.8f}", txt)
    if n != 1:
        raise SystemExit(f"{tpl_path}: expected one degauss line, found {n}")
    txt = re.sub(r"nat\s*=\s*\d+", f"nat = {len(at)}", txt)
    lines, out, i = txt.splitlines(), [], 0
    while i < len(lines):
        u = lines[i].strip().upper()
        if u.startswith("CELL_PARAMETERS"):
            out.append("CELL_PARAMETERS angstrom")
            out += ["  " + "  ".join(f"{x:.10f}" for x in row) for row in np.array(at.cell)]
            i += 4
            continue
        if u.startswith("ATOMIC_POSITIONS"):
            out.append("ATOMIC_POSITIONS angstrom")
            out += [f"  {s}  " + "  ".join(f"{x:.10f}" for x in p)
                    for s, p in zip(at.get_chemical_symbols(), at.get_positions())]
            i += 1 + nat_old
            continue
        out.append(lines[i]); i += 1
    if not any(l.startswith("CELL_PARAMETERS") for l in out):
        raise SystemExit(f"{tpl_path}: no CELL_PARAMETERS block (ibrav != 0?)")
    return "\n".join(out) + "\n"


res = json.load(open("phaseE/analysis_live_all.json"))["cells"]
calc = MACECalculator(model_paths="phaseD/models/fe_ft_Tel6000_f32.model", device="cuda",
                      default_dtype="float32", enable_cueq=True)
pfactor = (750 * units.fs) ** 2 * 1000 * units.GPa
frames = []
count = {"hcp": 0, "liquid": 0}

for cname in CELLS:
    Tm, P = res[cname]["Tm"], res[cname]["P"]
    thr = json.load(open(f"phaseE/cells/{cname}.json"))["q6_threshold"]
    a0 = read(f"phaseE/cells/{cname}.extxyz").cell.lengths()[0] / 12.0
    for phase in ("hcp", "liquid"):
        at = ortho_hcp(a0, 1.60, (4, 2, 4))          # 128 atoms
        at.calc = calc
        rng = np.random.RandomState(100 + len(frames))
        if phase == "liquid":
            MBD(at, temperature_K=1.6 * Tm, rng=rng)
            Langevin(at, 2 * units.fs, temperature_K=1.6 * Tm, friction=0.02).run(1000)
        else:
            MBD(at, temperature_K=Tm, rng=rng)
        dyn = NPT(at, 2 * units.fs, temperature_K=Tm, externalstress=P * units.GPa,
                  ttime=25 * units.fs, pfactor=pfactor)
        dyn.run(1500)                                  # 3 ps at (Tm, P)
        got, tries = 0, 0
        while got < 2 and tries < 6:
            dyn.run(1000)                              # 2 ps between snapshots
            tries += 1
            fs = float((q6_per_atom(at) > thr).mean())
            ok = fs > 0.8 if phase == "hcp" else fs < 0.2
            print(f"{cname} {phase:6s} snapshot {tries}: solid_frac {fs:.2f} "
                  f"V {at.get_volume()/len(at):.3f} A^3/atom -> {'kept' if ok else 'skipped'}", flush=True)
            if not ok:
                continue
            f = Atoms(at.get_chemical_symbols(), positions=at.get_positions(),
                      cell=np.array(at.cell), pbc=True)
            f.info.update(phase=phase, P_est=float(P), T=float(Tm), source=cname,
                          spot_name=f"spot_{phase[0]}{count[phase]:02d}")
            count[phase] += 1
            frames.append(f)
            got += 1

print(f"\nframes kept: hcp {count['hcp']}, liquid {count['liquid']}")
write("phaseD/spot_frames.xyz", frames)

tpl = {}
for phase in ("hcp", "liquid"):
    t = sorted(glob.glob(f"labels/{phase}-*_Tel6000.in"))
    if not t:
        raise SystemExit(f"no template labels/{phase}-*_Tel6000.in")
    tpl[phase] = t[0]
names = []
for f in frames:
    for tel in TELS:
        stem = f"{f.info['spot_name']}_T{tel}"
        open(f"labels/{stem}.in", "w").write(make_input(tpl[f.info["phase"]], f, stem, tel))
        names.append(f"{stem}.in")
open("labels/spot_list.txt", "w").write("\n".join(names) + "\n")
print(f"wrote {len(names)} inputs and labels/spot_list.txt")
print(f"templates: hcp <- {os.path.basename(tpl['hcp'])}, liquid <- {os.path.basename(tpl['liquid'])}")
print(f"degauss: 4600 K = {KB_RY*4600:.6f} Ry, 6500 K = {KB_RY*6500:.6f} Ry")
