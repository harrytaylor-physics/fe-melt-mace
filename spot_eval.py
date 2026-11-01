#!/usr/bin/env python3
"""
spot_eval.py — A21 step 4: the 3-temperature Tel model against new DFT at
Tel = 4600 and 6500 K, on 20 structures generated at the measured melting
conditions (10 hcp, 10 liquid; spot_md_gen.py).

DFT energy convention. With Fermi-Dirac smearing, the QE line
"!    total energy" is the Mermin free energy F = E - TS (QE also prints
"internal energy E=F+TS"). The campaign labels are checked against this
before anything is judged: a few campaign outputs are matched to training
frames by geometry and their REF_energy compared with the parsed value.

A21 criteria, applied separately at each Tel; PASS needs all three at both:
  (1) |liquid - hcp mean energy error| < 2.2 meV/atom  (~30 K of Tm)
  (2) energy RMSE after removing the uniform offset < 5 meV/atom
  (3) force RMSE < 110 meV/A

Run (tel_env, GPU): python spot_eval.py --model phaseD/models/fe_tel_e3_stagetwo.model
"""
import argparse, glob, hashlib, re
import numpy as np
from ase.io import read

RY = 13.605693123                 # eV
RY_BOHR = RY / 0.529177210903     # eV/A
TELS = (4600, 6500)

ap = argparse.ArgumentParser()
ap.add_argument("--model", required=True)
ap.add_argument("--device", default="cuda")
a = ap.parse_args()


def parse_out(path, nat):
    txt = open(path, errors="ignore").read()
    if "JOB DONE" not in txt:
        return None
    F = float(re.findall(r"^!\s+total energy\s+=\s+(\S+)\s+Ry", txt, re.M)[-1]) * RY
    blk = txt.split("Forces acting on atoms")[-1]
    f = np.array([[float(x) for x in m] for m in
                  re.findall(r"atom\s+\d+\s+type\s+\d+\s+force =\s+(\S+)\s+(\S+)\s+(\S+)", blk)[:nat]])
    return F, f * RY_BOHR


def parse_in_geometry(path):
    L = open(path).read().splitlines()
    nat = int(re.search(r"nat\s*=\s*(\d+)", "\n".join(L)).group(1))
    for i, l in enumerate(L):
        u = l.strip().upper()
        if u.startswith("CELL_PARAMETERS"):
            unit = u.lower(); cell = np.array([[float(x) for x in L[i + k].split()[:3]] for k in (1, 2, 3)])
            if "bohr" in unit: cell *= 0.529177210903
        if u.startswith("ATOMIC_POSITIONS"):
            punit = u.lower(); pos = np.array([[float(x) for x in L[i + 1 + k].split()[1:4]] for k in range(nat)])
    if "crystal" in punit: pos = pos @ cell
    elif "bohr" in punit: pos *= 0.529177210903
    return nat, cell, pos


def key(pos, cell):
    h = hashlib.sha1()
    h.update(np.round(pos, 4).tobytes()); h.update(np.round(cell, 4).tobytes())
    return h.hexdigest()


# ---- convention check against the campaign ----
ref = {}
for f in glob.glob("phaseD/*_Tel6000.xyz"):
    if "_stamped" in f:
        continue
    for at in read(f, index=":"):
        ref[key(at.get_positions(), np.array(at.cell))] = at.info["REF_energy"]
diffs = []
for p in sorted(glob.glob("labels/hcp-*_Tel6000.in"))[:300]:
    nat, cell, pos = parse_in_geometry(p)
    k = key(pos, cell)
    o = p[:-3] + ".out"
    if k in ref and glob.glob(o):
        r = parse_out(o, nat)
        if r:
            diffs.append(r[0] - ref[k])
    if len(diffs) >= 5:
        break
if diffs:
    print(f"convention check: parsed F minus campaign REF_energy over {len(diffs)} frames: "
          f"{np.mean(diffs):+.6f} eV (spread {np.ptp(diffs):.2e})")
    if np.ptp(diffs) > 1e-3:
        raise SystemExit("energy convention does not match the campaign labels: stopping")
else:
    print("convention check: no campaign frame matched by geometry; using the '!' free energy")

# ---- evaluate ----
from mace.calculators import MACECalculator
calc = MACECalculator(model_paths=a.model, device=a.device, default_dtype="float32")
frames = read("phaseD/spot_frames.xyz", index=":")
allpass = True
for tel in TELS:
    ph, eE, fE = [], [], []
    for fr in frames:
        r = parse_out(f"labels/{fr.info['spot_name']}_T{tel}.out", len(fr))
        if r is None:
            continue
        g = fr.copy(); g.info["elec_temp"] = float(tel); g.calc = calc
        eE.append(1000 * (g.get_potential_energy() - r[0]) / len(g))
        fE.append(g.get_forces() - r[1])
        ph.append(fr.info["phase"])
    ph, eE = np.array(ph), np.array(eE)
    if len(eE) == 0:
        print(f"\nTel {tel}: no finished DFT outputs"); allpass = False; continue
    nh, nl = (ph == "hcp").sum(), (ph == "liquid").sum()
    lh = eE[ph == "liquid"].mean() - eE[ph == "hcp"].mean()
    lhe = np.hypot(eE[ph == "liquid"].std(ddof=1) / np.sqrt(nl), eE[ph == "hcp"].std(ddof=1) / np.sqrt(nh))
    roff = eE.std()
    frmse = 1000 * np.sqrt((np.concatenate(fE) ** 2).mean())
    c1, c2, c3 = abs(lh) < 2.2, roff < 5.0, frmse < 110.0
    allpass &= c1 and c2 and c3
    print(f"\nTel {tel} K  (hcp {nh}, liquid {nl} frames)")
    print(f"  uniform offset {eE.mean():+.2f} meV/atom (does not affect Tm)")
    print(f"  (1) |liquid - hcp| {abs(lh):.2f} +- {lhe:.2f} < 2.2 meV/atom "
          f"(= {lh / (0.84 * 0.08617333):+.0f} K)        : {'pass' if c1 else 'FAIL'}")
    print(f"  (2) energy RMSE after offset {roff:.2f} < 5.0 meV/atom : {'pass' if c2 else 'FAIL'}")
    print(f"  (3) force RMSE {frmse:.1f} < 110 meV/A                 : {'pass' if c3 else 'FAIL'}")
    if abs(abs(lh) - 2.2) < lhe:
        print("  note: criterion (1) is within one standard error of its threshold: inconclusive")
print(f"\nSPOT CHECK VERDICT: {'PASS' if allpass else 'FAIL'}")
