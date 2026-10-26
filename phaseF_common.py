#!/usr/bin/env python3
"""
phaseF_common.py — shared pieces for Phase F of fe_melt: extending the Tel
potential to 50–90 GPa so the fcc–hcp–liquid triple point can be located.

Why Phase F exists
  On the melting line the training data stop at 9.20 A^3/atom for the solids
  and 9.00 A^3/atom for the liquid (about 80–85 GPa and 100–105 GPa), and at
  T_el = 4,000 K. The triple-point cells at 60 and 80 GPa need 9.2–9.8
  A^3/atom and T_el = 3,200–3,550 K.

Melting-temperature estimates (this work only, no literature fit)
  simon_T(P)      Simon fit to the seven hcp cells (pressure P_n):
                  T = 3846 K (1 + (P - 100 GPa)/120 GPa)^(1/2.15),
                  extrapolated below 100 GPa.
  tcell_of_v(V)   ln T_m quadratic in ln V, fitted to the same seven cells
                  (backup table B1, V = cell-mean volume per atom).
  vcell_of_p(P)   ln V quadratic in ln P_mean, same cells.
  A solid sits ~1% below the cell-mean volume and the liquid ~1% above
  (Delta V_m / V ~ 2%), so
     tm_solid(V)  = tcell_of_v(1.01 V)
     tm_liquid(V) = tcell_of_v(0.99 V)
  Check: tcell_of_v(vcell_of_p(60)) = 3,213 K vs simon_T(60) = 3,185 K;
         at 80 GPa 3,560 K vs 3,533 K.
"""
import glob, os, re
import numpy as np
from ase import Atoms, units

KB_RY = 8.617333262e-5 / 13.605693123          # Ry per K
RY = 13.605693123                              # eV
RY_BOHR = RY / 0.529177210903                  # eV/A

# seven hcp coexistence cells (backup table B1)
_V = np.array([8.94, 8.29, 7.82, 7.44, 7.13, 6.97, 6.87])           # A^3/atom, cell mean
_T = np.array([3911., 4602., 5196., 5734., 6211., 6488., 6656.])    # K
_P = np.array([102.5, 153.6, 206.0, 258.2, 311.7, 343.7, 365.6])   # GPa, cell mean incl. kinetic
_CT = np.polyfit(np.log(_V), np.log(_T), 2)
_CV = np.polyfit(np.log(_P), np.log(_V), 2)

# in-range reference errors of fe_tel_e3_stagetwo at T_el = 4,000 K (346 held-out labels)
REF_FRMSE_4000 = {"hcp": 89.0, "fcc": 75.0, "liquid": 156.0}   # meV/A
REF_ERMSE_4000 = {"hcp": 1.4, "fcc": 1.6, "liquid": 1.8}       # meV/atom, offset removed
DS_M_KB = 0.84                                                 # entropy of melting, k_B/atom
K_PER_MEV = 1.0 / (DS_M_KB * 0.08617333262)                    # 13.8 K per meV/atom


def simon_T(P):
    return 3846.0 * (1.0 + (P - 100.0) / 120.0) ** (1.0 / 2.15)


def tcell_of_v(V):
    return float(np.exp(np.polyval(_CT, np.log(V))))


def vcell_of_p(P):
    return float(np.exp(np.polyval(_CV, np.log(P))))


def tm_solid(V):
    return tcell_of_v(1.01 * V)


def tm_liquid(V):
    return tcell_of_v(0.99 * V)


# ---------------------------------------------------------------- structures
def hcp_hex(vat, reps=(4, 4, 4), coa=1.60):
    """Hexagonal hcp cell as in Phase C (128 atoms). Not upper-triangular: NVT only."""
    a = (4.0 * vat / (np.sqrt(3.0) * coa)) ** (1.0 / 3.0)
    cell = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, coa * a]])
    at = Atoms("Fe2", cell=cell, pbc=True, scaled_positions=[[0, 0, 0], [1 / 3., 2 / 3., 0.5]])
    return at.repeat(reps)


def hcp_ortho(vat, reps=(4, 2, 4), coa=1.60):
    """Orthorhombic hcp cell (128 atoms), usable with ASE NPT (upper-triangular)."""
    a = (2.0 * vat / (np.sqrt(3.0) / 2.0 * coa)) ** (1.0 / 3.0)
    c, b = coa * a, a * np.sqrt(3.0)
    base = np.array([[0, 0, 0], [a / 2, b / 2, 0], [a / 2, b / 6, c / 2], [0, 2 * b / 3, c / 2]])
    return Atoms("Fe4", positions=base, cell=[a, b, c], pbc=True).repeat(reps)


def fcc_cubic(vat, n=3):
    """Cubic fcc cell (108 atoms), as in Phase C."""
    a = (4.0 * vat) ** (1.0 / 3.0)
    at = Atoms("Fe4", cell=np.eye(3) * a, pbc=True,
               scaled_positions=[[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    return at.repeat((n, n, n))


def frac_msd(at, s0):
    """Mean-square displacement (A^2) from reference fractional positions s0,
    minimum image, centre-of-mass drift removed. Homogeneous strain of the cell
    (NPT) does not count as displacement."""
    ds = at.get_scaled_positions(wrap=False) - s0
    ds -= np.round(ds)
    d = ds @ np.array(at.cell)
    d -= d.mean(axis=0)
    return float((d * d).sum(axis=1).mean())


# ---------------------------------------------------------------- order parameter
def q6_reference(cells_dir="phaseE/cells"):
    """q-bar-6 threshold and cell-mean volume of the 100 GPa hcp cell. The cut-off
    radius is scaled with (V / V_ref)^(1/3) so the same neighbour shells are used."""
    import json
    from ase.io import read
    thr, vref = 0.212, 8.94
    js = os.path.join(cells_dir, "coex_hcp_100GPa.json")
    xs = os.path.join(cells_dir, "coex_hcp_100GPa.extxyz")
    if os.path.exists(js):
        thr = float(json.load(open(js))["q6_threshold"])
    if os.path.exists(xs):
        c = read(xs)
        vref = c.get_volume() / len(c)
    return thr, vref


def solid_fraction(at, thr, vref):
    from phaseE_q6 import q6_per_atom
    rcut = 3.2 * (at.get_volume() / len(at) / vref) ** (1.0 / 3.0)
    return float((q6_per_atom(at, rcut=rcut) > thr).mean())


# ---------------------------------------------------------------- calculator
def load_calc(model, device="cuda", dry_run=False):
    if dry_run:
        from ase.calculators.lj import LennardJones
        return LennardJones(sigma=2.10, epsilon=0.60, rc=5.0, smooth=True)
    from mace.calculators import MACECalculator
    return MACECalculator(model_paths=os.path.expanduser(model), device=device,
                          default_dtype="float32")


# ---------------------------------------------------------------- QE inputs
def template_for(phase, nat, labels_dir="labels"):
    """An existing campaign input of the same phase and atom count. Copying it keeps
    k-points, cutoffs, bands, pseudopotential, mixing and convergence identical."""
    for pat in (f"{phase}-*_Tel4000.in", f"{phase}-*_Tel6000.in", f"{phase}-*_Tel*.in"):
        for t in sorted(glob.glob(os.path.join(labels_dir, pat))):
            m = re.search(r"nat\s*=\s*(\d+)", open(t).read())
            if m and int(m.group(1)) == nat:
                return t
    raise SystemExit(f"no campaign input for phase {phase} with nat = {nat} in {labels_dir}/")


def make_input(tpl_path, at, stem, tel):
    src = open(tpl_path).read()
    nat_old = int(re.search(r"nat\s*=\s*(\d+)", src).group(1))
    if nat_old != len(at):
        raise SystemExit(f"{tpl_path}: nat {nat_old} != {len(at)}")
    txt = src.replace(os.path.basename(tpl_path)[:-3], stem)
    txt, n = re.subn(r"degauss\s*=\s*[0-9.eEdD+-]+", f"degauss = {KB_RY * tel:.8f}", txt)
    if n != 1:
        raise SystemExit(f"{tpl_path}: expected one degauss line, found {n}")
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
        out.append(lines[i])
        i += 1
    if not any(l.startswith("CELL_PARAMETERS") for l in out):
        raise SystemExit(f"{tpl_path}: no CELL_PARAMETERS block (ibrav != 0?)")
    return "\n".join(out) + "\n"


_FALLBACK_SBATCH = """#!/bin/bash
#SBATCH --job-name={job}
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --array=0-{nmax}%{throttle}
#SBATCH --output={job}_%A_%a.out

module purge
module load QuantumESPRESSO/7.5-foss-2025a
cd $SLURM_SUBMIT_DIR
INP=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" {listname})
T0=$SECONDS
srun pw.x -nk 4 -i "$INP" > "${{INP%.in}}.out"
EL=$((SECONDS-T0))
echo "TIMING $INP ${{EL}} s on 48 cores = $((EL*48/3600)) core-h"
"""


def write_array_sbatch(labels_dir, listname, sbname, job, n, throttle=40):
    """Copy the campaign's labels/run_array.sbatch (so module, partition and pw.x
    flags match the labels already computed) with a new list, array range and
    job name. Falls back to the Phase C template if that file is missing."""
    src = os.path.join(labels_dir, "run_array.sbatch")
    txt = open(src).read() if os.path.exists(src) else ""
    if "input_list.txt" in txt and "--array" in txt:
        txt = txt.replace("input_list.txt", listname)
        txt = re.sub(r"#SBATCH --array=\S+", f"#SBATCH --array=0-{n - 1}%{throttle}", txt)
        txt = re.sub(r"#SBATCH --job-name=\S+", f"#SBATCH --job-name={job}", txt)
        txt = re.sub(r"#SBATCH --output=\S+", f"#SBATCH --output={job}_%A_%a.out", txt)
        origin = "copied from labels/run_array.sbatch"
    else:
        txt = _FALLBACK_SBATCH.format(job=job, nmax=n - 1, throttle=throttle, listname=listname)
        origin = "Phase C template (labels/run_array.sbatch not found)"
    path = os.path.join(labels_dir, sbname)
    open(path, "w").write(txt)
    return path, origin


# ---------------------------------------------------------------- QE output
def parse_out(path, nat):
    """Mermin free energy F (eV, the '!' line with smearing), forces (eV/A) and
    static pressure (GPa) from a pw.x output; None if unfinished."""
    if not os.path.exists(path):
        return None
    txt = open(path, errors="ignore").read()
    if "JOB DONE" not in txt:
        return None
    F = float(re.findall(r"^!\s+total energy\s+=\s+(\S+)\s+Ry", txt, re.M)[-1]) * RY
    blk = txt.split("Forces acting on atoms")[-1]
    f = np.array([[float(x) for x in m] for m in
                  re.findall(r"atom\s+\d+\s+type\s+\d+\s+force =\s+(\S+)\s+(\S+)\s+(\S+)", blk)[:nat]])
    pm = re.findall(r"P=\s*(-?[\d.]+)", txt)
    P = float(pm[-1]) / 10.0 if pm else float("nan")
    return F, f * RY_BOHR, P
