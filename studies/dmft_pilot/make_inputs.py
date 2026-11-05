#!/usr/bin/env python3
"""
make_inputs.py — QE + Wannier90 inputs for the DFT+DMFT pilot (fe_melt, DMFT step D0).

System: hcp Fe, 2 atoms, 7.05 A^3/atom (inner-core density used by Pourovskii's
DFT+DMFT work), c/a = 1.60, T = 5,800 K (Fermi-Dirac smearing = k_B T, as in the
campaign labels). Pseudopotential, cutoffs and pseudo_dir are copied from a
campaign label input, so the DFT part matches the labels.

Stages
  scf  writes scf.in (16x16x10 k, 40 bands).
  w90  reads E_F from scf.out and writes
         nscf.in        full 10x10x6 grid, no symmetry, verbosity high (output
                        name fe.nscf.out, as TRIQS/DFTTools expects later)
         pw2wan.in      amn/mmn/eig for the Fe 3d projections
         fe.win         base Wannier90 input (used for the -pp step)
         fe_w1.win, fe_w2.win   projected-orbital runs (no localisation:
                        dis_num_iter = num_iter = 0) in two outer windows,
                        w1 = [E_F - 10, E_F + 4] eV, w2 = [E_F - 12, E_F + 8] eV.
                        Z is computed for both; the spread between them is the
                        projector-choice uncertainty.
  The 8 lowest bands (3s, 3p semicore of 2 atoms) are excluded from Wannier90.

Usage (run by pilot_qe_w90.sbatch):
  python3 make_inputs.py scf --labels ~/fe_melt/labels
  python3 make_inputs.py w90
"""
import argparse, glob, json, math, os, re, sys

V_AT, COA, T_K = 7.05, 1.60, 5800.0
KB_RY = 8.617333262e-5 / 13.605693123
NK_SCF = (16, 16, 10)
NK_W90 = (10, 10, 6)
NBND, NEXCL, NWANN = 40, 8, 10
WINDOWS = {"w1": (-10.0, 4.0), "w2": (-12.0, 8.0)}   # eV relative to E_F


def cell():
    a = (4.0 * V_AT / (math.sqrt(3.0) * COA)) ** (1.0 / 3.0)
    c = COA * a
    return a, c, [[a, 0.0, 0.0], [-a / 2, a * math.sqrt(3.0) / 2, 0.0], [0.0, 0.0, c]]


FRAC = [(1 / 3., 2 / 3., 0.25), (2 / 3., 1 / 3., 0.75)]


def campaign_settings(labels):
    tpls = sorted(glob.glob(os.path.join(os.path.expanduser(labels), "hcp-*_Tel4000.in")))
    if not tpls:
        sys.exit(f"no campaign input labels/hcp-*_Tel4000.in under {labels}")
    t = open(tpls[0]).read()
    pd = re.search(r"pseudo_dir\s*=\s*'([^']+)'", t).group(1)
    pd = os.path.expanduser(pd)
    if not os.path.isabs(pd):
        pd = os.path.abspath(os.path.join(os.path.expanduser(labels), pd))
    spec = re.search(r"ATOMIC_SPECIES\s*\n\s*(\S+)\s+(\S+)\s+(\S+)", t)
    ecw = re.search(r"ecutwfc\s*=\s*([0-9.dDeE+-]+)", t).group(1)
    ecr = re.search(r"ecutrho\s*=\s*([0-9.dDeE+-]+)", t).group(1)
    s = dict(pseudo_dir=pd, upf=spec.group(3), mass=spec.group(2), ecutwfc=ecw, ecutrho=ecr,
             template=os.path.basename(tpls[0]))
    if not os.path.exists(os.path.join(pd, s["upf"])):
        sys.exit(f"pseudopotential not found: {os.path.join(pd, s['upf'])}")
    return s


def pw_input(calc, s, kblock, extra_sys="", extra_ctl=""):
    a, c, L = cell()
    cellblk = "\n".join("  " + "  ".join(f"{x:.10f}" for x in v) for v in L)
    pos = "\n".join(f"  Fe  {x:.10f}  {y:.10f}  {z:.10f}" for x, y, z in FRAC)
    return f"""&CONTROL
  calculation = '{calc}'
  prefix = 'fe'
  outdir = './tmp'
  pseudo_dir = '{s["pseudo_dir"]}'
  verbosity = 'high'
{extra_ctl}/
&SYSTEM
  ibrav = 0
  nat = 2
  ntyp = 1
  ecutwfc = {s["ecutwfc"]}
  ecutrho = {s["ecutrho"]}
  nbnd = {NBND}
  occupations = 'smearing'
  smearing = 'fermi-dirac'
  degauss = {KB_RY * T_K:.8f}
  nspin = 1
{extra_sys}/
&ELECTRONS
  conv_thr = 1.0d-10
  mixing_beta = 0.3
/
ATOMIC_SPECIES
  Fe {s["mass"]} {s["upf"]}
CELL_PARAMETERS angstrom
{cellblk}
ATOMIC_POSITIONS crystal
{pos}
{kblock}
"""


def kgrid(nk):
    pts = []
    for i in range(nk[0]):
        for j in range(nk[1]):
            for k in range(nk[2]):
                pts.append((i / nk[0], j / nk[1], k / nk[2]))
    return pts


def win_text(ef, window=None):
    a, c, L = cell()
    ks = kgrid(NK_W90)
    lines = [f"num_wann = {NWANN}", f"num_bands = {NBND - NEXCL}", f"exclude_bands = 1-{NEXCL}",
             f"mp_grid = {NK_W90[0]} {NK_W90[1]} {NK_W90[2]}", "",
             "begin unit_cell_cart", "ang"] + \
            ["  " + "  ".join(f"{x:.10f}" for x in v) for v in L] + ["end unit_cell_cart", "",
             "begin atoms_frac"] + [f"  Fe  {x:.10f}  {y:.10f}  {z:.10f}" for x, y, z in FRAC] + \
            ["end atoms_frac", "", "begin projections", "  Fe:d", "end projections", ""]
    if window is not None:
        lo, hi = window
        lines += [f"dis_win_min = {ef + lo:.4f}", f"dis_win_max = {ef + hi:.4f}",
                  "dis_num_iter = 0", "num_iter = 0", "write_hr = .true.",
                  "write_u_matrices = .true.", "write_xyz = .true.", ""]
    lines += ["begin kpoints"] + [f"  {x:.10f}  {y:.10f}  {z:.10f}" for x, y, z in ks] + ["end kpoints"]
    return "\n".join(lines) + "\n"


def fermi(out="scf.out"):
    t = open(out).read()
    if "JOB DONE" not in t:
        sys.exit("scf.out is not finished")
    return float(re.findall(r"the Fermi energy is\s+(-?[\d.]+)\s+ev", t)[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["scf", "w90"])
    ap.add_argument("--labels", default="~/fe_melt/labels")
    a = ap.parse_args()
    if a.stage == "scf":
        s = campaign_settings(a.labels)
        json.dump(s, open("settings.json", "w"), indent=1)
        k = f"K_POINTS automatic\n  {NK_SCF[0]} {NK_SCF[1]} {NK_SCF[2]} 0 0 0"
        open("scf.in", "w").write(pw_input("scf", s, k))
        ac = cell()
        print(f"hcp Fe: a = {ac[0]:.4f} A, c = {ac[1]:.4f} A, V = {V_AT} A^3/atom, T = {T_K:.0f} K, "
              f"degauss = {KB_RY * T_K:.6f} Ry; settings from {s['template']}: {s['upf']}, "
              f"ecutwfc {s['ecutwfc']}, ecutrho {s['ecutrho']}")
        return
    s = json.load(open("settings.json"))
    ef = fermi()
    ks = kgrid(NK_W90)
    w = 1.0 / len(ks)
    k = "K_POINTS crystal\n" + f"  {len(ks)}\n" + "\n".join(f"  {x:.10f} {y:.10f} {z:.10f} {w:.10e}" for x, y, z in ks)
    open("nscf.in", "w").write(pw_input("nscf", s, k, extra_sys="  nosym = .true.\n  noinv = .true.\n"))
    open("pw2wan.in", "w").write("&inputpp\n  outdir = './tmp'\n  prefix = 'fe'\n  seedname = 'fe'\n"
                                 "  write_mmn = .true.\n  write_amn = .true.\n  write_unk = .false.\n/\n")
    open("fe.win", "w").write(win_text(ef))
    for name, win in WINDOWS.items():
        open(f"fe_{name}.win", "w").write(win_text(ef, win))
    json.dump(dict(E_F_eV=ef, windows={k_: [ef + v[0], ef + v[1]] for k_, v in WINDOWS.items()}),
              open("fermi.json", "w"), indent=1)
    print(f"E_F = {ef:.4f} eV; windows " + ", ".join(f"{k_}: [{ef + v[0]:.2f}, {ef + v[1]:.2f}] eV"
                                                     for k_, v in WINDOWS.items()))


if __name__ == "__main__":
    main()
