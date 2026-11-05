#!/usr/bin/env python3
"""
spot_dft_gen.py — inputs for the Tel-model spot check (A21).

Picks 10 hcp and 10 liquid frames from the held-out Tel-6000 test set, at
pressures across the melting band (95-360 GPa), that appear in NO training set
at ANY electronic temperature (geometry match against all train_*.xyz files).
hcp: the hottest frames (closest to melting). liquid: the coolest frames.

Writes pw.x inputs at Tel = 4600 K and 6500 K (Tm at ~150 and ~330 GPa). Each
input is a copy of an existing label input of the same phase, so k-points,
cutoffs, pseudopotential, mixing and convergence settings are identical to the
campaign; only the cell, positions, file stem and degauss change.

Outputs (in labels/):
  spot_{h|l}{NN}_T{4600|6500}.in      40 inputs
  spot_list.txt                       the 40 names, one per line
and phaseD/spot_frames.xyz (the 20 geometries with phase, pressure, name).

Run (login node, production env): python3.13 spot_dft_gen.py
"""
import glob, hashlib, os, re
import numpy as np
from ase.io import read, write

TELS = (4600, 6500)
KB_RY = 8.617333262e-5 / 13.605693123   # Ry/K
P_BAND = (95.0, 360.0)
N_PER_PHASE = 10


def key(at):
    h = hashlib.sha1()
    h.update(np.round(at.get_positions(), 4).tobytes())
    h.update(np.round(np.array(at.cell), 4).tobytes())
    return h.hexdigest()


def pressure(at):
    for k in ("P_est", "P", "pressure"):
        if k in at.info:
            return float(at.info[k])
    return float("nan")


trained = set()
for f in glob.glob("phaseD/train_Tel*.xyz"):
    if "_stamped" in f:
        continue
    trained.update(key(a) for a in read(f, index=":"))
print(f"training geometries (all Tels): {len(trained)}")

test = [a for a in read("phaseD/test_Tel6000.xyz", index=":") if key(a) not in trained]
print(f"held-out test frames not in any training set: {len(test)}")

picked = []
for phase, hottest in (("hcp", True), ("liquid", False)):
    cand = [a for a in test if a.info.get("phase") == phase
            and P_BAND[0] <= pressure(a) <= P_BAND[1]]
    if len(cand) < N_PER_PHASE:
        raise SystemExit(f"only {len(cand)} {phase} candidates in {P_BAND} GPa; "
                         "check the pressure key in the frame info")
    # keep the hottest (hcp) / coolest (liquid) half, then spread across pressure
    cand.sort(key=lambda a: float(a.info.get("T", 0)), reverse=hottest)
    cand = cand[:max(N_PER_PHASE, len(cand) // 2)]
    cand.sort(key=pressure)
    idx = np.unique(np.linspace(0, len(cand) - 1, N_PER_PHASE).round().astype(int))
    for n, i in enumerate(idx):
        a = cand[i].copy()
        a.info["spot_name"] = f"spot_{phase[0]}{n:02d}"
        picked.append(a)
        print(f"  {a.info['spot_name']}  {phase:6s} P = {pressure(a):6.1f} GPa  "
              f"T = {float(a.info.get('T', 0)):6.0f} K")


def template_for(phase):
    t = sorted(glob.glob(f"labels/{phase}-*_Tel6000.in"))
    if not t:
        raise SystemExit(f"no template labels/{phase}-*_Tel6000.in")
    return t[0]


def make_input(tpl_path, at, stem, tel):
    txt = open(tpl_path).read()
    old_stem = os.path.basename(tpl_path)[:-3]
    txt = txt.replace(old_stem, stem)
    txt, n = re.subn(r"degauss\s*=\s*[0-9.eEdD+-]+", f"degauss = {KB_RY * tel:.8f}", txt)
    if n != 1:
        raise SystemExit(f"{tpl_path}: expected one degauss line, found {n}")
    txt = re.sub(r"nat\s*=\s*\d+", f"nat = {len(at)}", txt)
    lines = txt.splitlines()
    out, i = [], 0
    nat_old = None
    m = re.search(r"nat\s*=\s*(\d+)", open(tpl_path).read())
    nat_old = int(m.group(1))
    while i < len(lines):
        L = lines[i]
        u = L.strip().upper()
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
        out.append(L)
        i += 1
    if not any(l.startswith("CELL_PARAMETERS") for l in out):
        raise SystemExit(f"{tpl_path}: no CELL_PARAMETERS block (ibrav != 0?)")
    return "\n".join(out) + "\n"


names = []
for a in picked:
    tpl = template_for(a.info["phase"])
    for tel in TELS:
        stem = f"{a.info['spot_name']}_T{tel}"
        open(f"labels/{stem}.in", "w").write(make_input(tpl, a, stem, tel))
        names.append(f"{stem}.in")
open("labels/spot_list.txt", "w").write("\n".join(names) + "\n")
write("phaseD/spot_frames.xyz", picked)
print(f"\nwrote {len(names)} inputs, labels/spot_list.txt, phaseD/spot_frames.xyz")
print(f"templates: hcp <- {os.path.basename(template_for('hcp'))}, "
      f"liquid <- {os.path.basename(template_for('liquid'))}")
print(f"degauss: 4600 K = {KB_RY*4600:.6f} Ry, 6500 K = {KB_RY*6500:.6f} Ry")
