#!/usr/bin/env python3
"""
phaseC_qe_labels.py — QE label inputs at three electronic temperatures (fe_melt_v1, C4)

Reads selected/manifest.json + frames_<phase>.extxyz and writes, per frame and
per Tel in {4000, 6000, 8000} K, a QE scf input using the Phase A locked recipe:
  PBE, Fe.pbe-spn-kjpaw_psl.1.0.0 (16-valence), ecutwfc 100 Ry, ecutrho 1000 Ry,
  Fermi-Dirac smearing with degauss = k_B * Tel, NM, tprnfor+tstress.
k-mesh: 2x2x2 for >=100-atom cells, 3x3x3 below (design C3 / Phase A density rule).

degauss (Ry): 4000 K -> 0.025334 ; 6000 K -> 0.038002 ; 8000 K -> 0.050669
Same geometry is written for all three Tel — only SCF is repeated, per design.

Also writes:
  labels/run_array.sbatch   SLURM array over all inputs (ARC), with the FIRST
                            task timed and reported — that timing IS gate G-A4
                            and feeds the budget kill-switch (>1.5 M core-hours
                            projected -> drop to 2 Tel / 2,000 frames, log C-amend).
  labels/collect.py         parses outputs -> labels.jsonl (F_mermin eV, forces
                            eV/A, stress GPa, per frame per Tel) for MACE training.

Usage: python phaseC_qe_labels.py [--seldir selected] [--outdir labels]
"""
import argparse, glob, json, os
import numpy as np
from ase.io import read

DEGAUSS = {4000: 0.025334, 6000: 0.038002, 8000: 0.050669}
PSEUDO = "~/scw/fe_melt/pseudo"

TEMPLATE = """&CONTROL
  calculation = 'scf'
  prefix = '{tag}'
  pseudo_dir = '{pseudo}'
  outdir = './tmp_{tag}'
  tprnfor = .true.
  tstress = .true.
  disk_io = 'none'
/
&SYSTEM
  ibrav = 0
  nat = {nat}
  ntyp = 1
  ecutwfc = 100
  ecutrho = 1000
  nbnd = {nbnd}
  occupations = 'smearing'
  smearing = 'fermi-dirac'
  degauss = {degauss}
  nspin = 1
/
&ELECTRONS
  conv_thr = 1.0d-8
  mixing_beta = 0.3
  electron_maxstep = 250
/
ATOMIC_SPECIES
  Fe 55.845 Fe.pbe-spn-kjpaw_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
{cell}
ATOMIC_POSITIONS angstrom
{pos}
K_POINTS automatic
  {k} {k} {k} 0 0 0
"""

SBATCH = """#!/bin/bash
#SBATCH --job-name=fe_labels
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --array=0-{nmax}%40
#SBATCH --output=lab_%A_%a.out

module purge
module load QuantumESPRESSO/7.5-foss-2025a
cd $SLURM_SUBMIT_DIR
INP=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" input_list.txt)
T0=$SECONDS
srun pw.x -nk 4 -i "$INP" > "${{INP%.in}}.out"
EL=$((SECONDS-T0))
echo "TIMING $INP ${{EL}} s on 48 cores = $((EL*48/3600)) core-h"
# G-A4: the first completed task's TIMING line, x {ntot} inputs, is the
# projected label cost. Kill-switch: if projection > 1.5 M core-hours,
# STOP the array (scancel on the job id) and invoke the 2-Tel amendment.
"""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seldir", default="selected")
    ap.add_argument("--outdir", default="labels")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    inputs = []
    for fx in sorted(glob.glob(f"{a.seldir}/frames_*.extxyz")):
        for fr in read(fx, index=":"):
            fid = fr.info["frame_id"]; nat = len(fr)
            nbnd = int(np.ceil(nat * 16 / 2 * 1.30))          # 30% virtual headroom
            k = 2 if nat >= 100 else 3
            cell = "\n".join("  " + " ".join(f"{x:14.8f}" for x in v) for v in fr.cell[:])
            pos = "\n".join("  Fe " + " ".join(f"{x:14.8f}" for x in p)
                            for p in fr.get_positions())
            for tel, dg in DEGAUSS.items():
                tag = f"{fid}_Tel{tel}"
                path = os.path.join(a.outdir, f"{tag}.in")
                open(path, "w").write(TEMPLATE.format(
                    tag=tag, pseudo=PSEUDO, nat=nat, nbnd=nbnd,
                    degauss=dg, cell=cell, pos=pos, k=k))
                inputs.append(f"{tag}.in")
    open(os.path.join(a.outdir, "input_list.txt"), "w").write("\n".join(inputs) + "\n")
    open(os.path.join(a.outdir, "run_array.sbatch"), "w").write(
        SBATCH.format(nmax=len(inputs) - 1, ntot=len(inputs)))
    collect = '''#!/usr/bin/env python3
"""collect.py — parse QE outputs -> labels.jsonl (one line per frame per Tel)."""
import glob, json, re
RY, RYA, RYB3 = 13.6056980659, 25.71104309541616, 14710.507848  # eV, eV/A, GPa
with open("labels.jsonl", "w") as out:
    for f in sorted(glob.glob("*.out")):
        t = open(f).read()
        if "JOB DONE" not in t:
            print("INCOMPLETE:", f); continue
        # '!' line total energy is the Mermin free energy F = E - TS with smearing
        F = float(re.findall(r"^!.*?=\\s*(-?[\\d.]+)\\s*Ry", t, re.M)[-1]) * RY
        forces = [[float(x) * RYA for x in m]
                  for m in re.findall(r"atom\\s+\\d+ type\\s+\\d+\\s+force =\\s*"
                                      r"(-?[\\d.]+)\\s+(-?[\\d.]+)\\s+(-?[\\d.]+)", t)]
        P = float(re.findall(r"P=\\s*(-?[\\d.]+)", t)[-1]) / 10.0
        tag = f[:-4]
        fid, tel = tag.rsplit("_Tel", 1)
        out.write(json.dumps(dict(id=fid, Tel=int(tel), F_eV=F,
                                  forces_eVA=forces, P_GPa=P)) + "\\n")
print("wrote labels.jsonl")
'''
    open(os.path.join(a.outdir, "collect.py"), "w").write(collect)
    n = len(inputs)
    print(f"{n} inputs written ({n//3} frames x 3 Tel) -> {a.outdir}/")
    print("Order of operations on ARC:")
    print("  1) sbatch run_array.sbatch  — but FIRST submit with --array=0-0 alone;")
    print("     the TIMING line of that single task is gate G-A4.")
    print(f"  2) projected cost = timing x {n}; kill-switch at 1.5 M core-hours.")
    print("  3) full array, then: python collect.py")

if __name__ == "__main__":
    main()
