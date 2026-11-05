#!/usr/bin/env python3
"""
wan_check.py — checks on the Fe 3d Wannier Hamiltonians before any DMFT
(fe_melt, DMFT step D0). Needs numpy (run inside the TRIQS environment).

For each window (fe_w1, fe_w2) it builds H(k) from <seed>_hr.dat on a
20x20x12 grid and reports, at the DFT Fermi level and T = 5,800 K:
  - total 3d occupation per atom and per orbital (5 d orbitals x 2 atoms)
  - the 3d band centre relative to E_F and the 3d band extent
  - the final Wannier spreads from <seed>.wout
What to look for: the two atoms identical (they are equivalent by symmetry);
d occupation per atom between ~6.5 and ~7.5; band centre a few eV below E_F.
A large difference between w1 and w2 means the projector choice matters and
both must be carried into DMFT.

Usage: python wan_check.py
"""
import json, re
import numpy as np

T_K, KB = 5800.0, 8.617333262e-5
NK = (20, 20, 12)


def read_hr(path):
    L = open(path).read().split("\n")
    nw, nr = int(L[1]), int(L[2])
    deg, i = [], 3
    while len(deg) < nr:
        deg += [int(x) for x in L[i].split()]
        i += 1
    H, R = {}, []
    for line in L[i:]:
        p = line.split()
        if len(p) < 7:
            continue
        r = (int(p[0]), int(p[1]), int(p[2]))
        m, n = int(p[3]) - 1, int(p[4]) - 1
        if r not in H:
            H[r] = np.zeros((nw, nw), complex)
            R.append(r)
        H[r][m, n] = float(p[5]) + 1j * float(p[6])
    Rs = np.array(R, float)
    Hs = np.array([H[r] for r in R]) / np.array(deg, float)[:, None, None]
    return nw, Rs, Hs


def bands(Rs, Hs):
    ks = np.array([(i / NK[0], j / NK[1], k / NK[2]) for i in range(NK[0])
                   for j in range(NK[1]) for k in range(NK[2])])
    ph = np.exp(2j * np.pi * ks @ Rs.T)                       # (nk, nR)
    Hk = np.einsum("kr,rmn->kmn", ph, Hs)
    Hk = 0.5 * (Hk + np.conj(np.transpose(Hk, (0, 2, 1))))
    return np.linalg.eigh(Hk)


ef = json.load(open("fermi.json"))["E_F_eV"]
for seed in ("fe_w1", "fe_w2"):
    try:
        nw, Rs, Hs = read_hr(f"{seed}_hr.dat")
    except FileNotFoundError:
        print(f"{seed}: no {seed}_hr.dat")
        continue
    e, U = bands(Rs, Hs)
    f = 1.0 / (np.exp(np.clip((e - ef) / (KB * T_K), -200, 200)) + 1.0)
    occ_orb = 2.0 * np.einsum("kmb,kb->m", np.abs(U) ** 2, f) / len(e)   # spin x 2
    dens = np.abs(U) ** 2
    centre = float((e[:, None, :] * dens).sum() / dens.sum()) - ef
    spreads = re.findall(r"WF centre and spread\s+\d+\s+\(.*?\)\s+([\d.]+)", open(f"{seed}.wout").read())
    sp = [float(x) for x in spreads[-nw:]] if spreads else []
    print(f"\n{seed}: {nw} Wannier orbitals")
    print(f"  3d occupation per atom: {occ_orb[:5].sum():.3f} (atom 1)  {occ_orb[5:].sum():.3f} (atom 2)")
    print("  per orbital, atom 1:   " + "  ".join(f"{x:.3f}" for x in occ_orb[:5]))
    print("  per orbital, atom 2:   " + "  ".join(f"{x:.3f}" for x in occ_orb[5:]))
    print(f"  3d band centre {centre:+.2f} eV vs E_F; 3d bands span {e.min() - ef:+.2f} to {e.max() - ef:+.2f} eV")
    if sp:
        print(f"  final spreads {min(sp):.3f}-{max(sp):.3f} A^2")
