"""
val_dft.py - static DFT references for validating the MACE potential (QE 6.6, pslibrary 1.0.0 PBE PAW,
ecutwfc 100 / ecutrho 1000 Ry, Fermi-Dirac smearing at the electronic temperature).

Structures (all perfect or homogeneously strained lattices, static ions):
  eos_V*      hcp, c/a 1.60, V = 6.6 ... 8.9 A^3/atom, T_el 6000 K
  tel_T*      hcp, V 7.0, T_el 4000 / 8000 K   (6000 K is eos_V7.00)
  fcc_V*, bcc_V*   1-atom cells at V 7.0 and 8.0, T_el 6000 K
  coa_*       hcp V 6.8, c/a 1.56 / 1.64 (1.60 from eos)
  el_*        hcp V 7.0 strained by +-1% (e1 = xx, e3 = zz, e4 = yz shear, e6 = xy shear), clamped ions
  ph_*        hcp V 7.0, E2g displacement +-0.03 A along x (atoms moved in opposite directions)
Each run writes <tag>.json: Mermin free energy F_per_atom (eV/atom), pressure P (GPa, QE 'P=' line),
sigma_qe (GPa, the stress tensor as QE prints it; physical stress = -sigma_qe), forces (eV/A),
cell, positions, tel, V. Raw .out files are not kept in the repository.
Usage: python val_dft.py [group ...]   (groups: eos tel fcc bcc coa el ph)
The committed references were made with QE 6.6 and a pslibrary 1.0.0 PBE PAW file regenerated
with ld1.x 6.6 (identical to the official file to 2e-6 relative; see pseudo_lda/).
"""
import os, sys, json, re, subprocess, time
import numpy as np
from ase import Atoms
from ase.build import bulk

import paths
QEBIN = paths.QE_PW
PSEUDO = paths.PSEUDO_DIR
UPF = "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF"
ROOT = paths.DFT_REFS
RY = 13.605693123
KB_RY = 8.617333262e-5 / RY
BOHR = 0.529177210903
RYB3_GPA = 14710.507848260711   # Ry/bohr^3 -> GPa


def hcp(V, coa=1.60):
    a = (4.0 * V / (np.sqrt(3) * coa)) ** (1 / 3)
    return bulk("Fe", "hcp", a=a, c=coa * a)


def structures():
    S = {}
    for V in (6.6, 6.8, 7.0, 7.2, 7.5, 8.0, 8.6, 8.9):
        S[f"eos_V{V:.2f}"] = (hcp(V), 6000.0, (12, 12, 8))
    for T in (4000.0, 8000.0):
        S[f"tel_T{T:.0f}"] = (hcp(7.0), T, (12, 12, 8))
    for V in (7.0, 8.0):
        S[f"fcc_V{V:.2f}"] = (bulk("Fe", "fcc", a=(4 * V) ** (1 / 3)), 6000.0, (20, 20, 20))
        S[f"bcc_V{V:.2f}"] = (bulk("Fe", "bcc", a=(2 * V) ** (1 / 3)), 6000.0, (20, 20, 20))
    for coa in (1.56, 1.64):
        S[f"coa_{coa:.2f}"] = (hcp(6.8, coa), 6000.0, (12, 12, 8))
    eps = 0.01
    for mode in ("e1", "e3", "e4", "e6"):
        for sgn in (+1, -1):
            at = hcp(7.0)
            E = np.zeros((3, 3))
            if mode == "e1":
                E[0, 0] = sgn * eps
            elif mode == "e3":
                E[2, 2] = sgn * eps
            elif mode == "e4":
                E[1, 2] = E[2, 1] = sgn * eps / 2
            else:
                E[0, 1] = E[1, 0] = sgn * eps / 2
            at.set_cell(np.array(at.cell) @ (np.eye(3) + E), scale_atoms=True)
            S[f"el_{mode}{'p' if sgn > 0 else 'm'}"] = (at, 6000.0, (12, 12, 8))
    for sgn in (+1, -1):
        at = hcp(7.0)
        p = at.get_positions()
        p[0, 0] += sgn * 0.03
        p[1, 0] -= sgn * 0.03
        at.set_positions(p)
        S[f"ph_{'p' if sgn > 0 else 'm'}"] = (at, 6000.0, (12, 12, 8))
    return S


def qe_input(tag, at, tel, k):
    nb = 20 * len(at)
    cell = "\n".join("  " + "  ".join(f"{x:.10f}" for x in row) for row in np.array(at.cell))
    pos = "\n".join(f"  Fe  " + "  ".join(f"{x:.10f}" for x in p) for p in at.get_positions())
    return f"""&CONTROL
  calculation = 'scf', prefix = '{tag}', outdir = './tmp_{tag}', pseudo_dir = '{PSEUDO}',
  tstress = .true., tprnfor = .true., disk_io = 'none'
/
&SYSTEM
  ibrav = 0, nat = {len(at)}, ntyp = 1, ecutwfc = 100, ecutrho = 1000,
  occupations = 'smearing', smearing = 'fd', degauss = {KB_RY * tel:.8f}, nbnd = {nb}
/
&ELECTRONS
  conv_thr = 1.0d-10, mixing_beta = 0.3
/
ATOMIC_SPECIES
  Fe 55.845 {UPF}
CELL_PARAMETERS angstrom
{cell}
ATOMIC_POSITIONS angstrom
{pos}
K_POINTS automatic
  {k[0]} {k[1]} {k[2]} 0 0 0
"""


def parse(path, nat):
    t = open(path, errors="ignore").read()
    if "JOB DONE" not in t:
        return None
    F = float(re.findall(r"^!\s+total energy\s+=\s+(\S+)\s+Ry", t, re.M)[-1]) * RY
    blk = t.split("total   stress")[-1].splitlines()[1:4]
    sig = np.array([[float(x) for x in l.split()[:3]] for l in blk]) * RYB3_GPA   # as printed by QE: -physical stress
    fb = t.split("Forces acting on atoms")[-1]
    f = np.array([[float(x) for x in m] for m in
                  re.findall(r"atom\s+\d+\s+type\s+\d+\s+force =\s+(\S+)\s+(\S+)\s+(\S+)", fb)[:nat]]) * RY / BOHR
    P = float(re.findall(r"P=\s*(-?[\d.]+)", t)[-1]) / 10.0
    return dict(F_per_atom=F / nat, P=P, sigma_qe=sig.tolist(), forces=f.tolist())


def main(groups):
    os.makedirs(ROOT, exist_ok=True)
    S = structures()
    for tag, (at, tel, k) in S.items():
        if groups and not any(tag.startswith(g) for g in groups):
            continue
        js = os.path.join(ROOT, tag + ".json")
        if os.path.exists(js):
            continue
        open(os.path.join(ROOT, tag + ".in"), "w").write(qe_input(tag, at, tel, k))
        t = time.time()
        with open(os.path.join(ROOT, tag + ".out"), "w") as fo:
            subprocess.run([QEBIN, "-i", tag + ".in"], cwd=ROOT, stdout=fo, stderr=subprocess.STDOUT,
                           env=dict(os.environ, OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "2")))
        r = parse(os.path.join(ROOT, tag + ".out"), len(at))
        if r is None:
            print(f"FAILED {tag}", flush=True)
            continue
        r.update(tag=tag, tel=tel, V=at.get_volume() / len(at), cell=np.array(at.cell).tolist(),
                 positions=at.get_positions().tolist(), t=time.time() - t)
        json.dump(r, open(js, "w"))
        subprocess.run(["rm", "-rf", os.path.join(ROOT, f"tmp_{tag}")])
        print(f"{tag}: F {r['F_per_atom']:.5f} eV/atom, P {r['P']:.2f} GPa ({r['t']:.0f} s)", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
