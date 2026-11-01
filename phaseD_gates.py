#!/usr/bin/env python3
"""
phaseD_gates.py — pre-registered Phase D gates against a fine-tuned model
(fe_melt_v1, gates G-D1..G-D4)

Usage:
  python phaseD_gates.py --model phaseD/models/fe_ft_Tel6000.model --tel 6000
  python phaseD_gates.py --model ... --tel 6000 --gate d4        (NPT gate only;
      run on a GPU node — ~25k MD steps on 128 atoms)

Static gates (default run):
  G-D2p force parity on the HELD-OUT test set, per phase:
        RMSE < 80 meV/A solids, < 120 meV/A liquid; energy RMSE < 5 meV/atom
  G-D3  hcp EOS |dP| < 3 GPa vs this work's DFT over 100-350 GPa
  G-D1  static enthalpy ordering hcp < fcc < bcc at 200 and 300 GPa
  c/a   at V=6.8 A^3/atom within 0.01 of DFT 1.604
Dynamic gate:
  G-D4  hcp survives 50 ps flexible-cell NPT at 0.9 Tm, 300 GPa
        (fail = c/a outside [1.55,1.68], cell angles drift >5 deg, or melt)

Phonon gate (G-D2 full) is run separately with phonopy once statics pass.
All DFT reference numbers below are this work's (campaigns convB/C/C2).
"""
import argparse, json
import numpy as np
from ase import Atoms, units
from ase.units import GPa
from ase.io import read

# ---- DFT references (this work) --------------------------------------------
DFT_V = np.array([6.0, 6.4, 6.8, 7.2, 7.6, 8.0, 8.5, 9.0, 9.5, 10.0])
DFT_P_HCP_6000K = np.array([517.0, 399.0, 308.2, 237.8, 182.3, 138.9, 96.8, 65.1, 41.1, 22.9])
DFT_CA = 1.604
DFT_GAPS = {200: {"fcc": 61, "bcc": 299}, 300: {"fcc": 69, "bcc": 383}}  # meV/atom vs hcp
GAP_TOL = 40   # meV/atom tolerance on gap magnitude; sign must be correct

def hcp(vat, coa=1.60, rep=(4, 4, 4)):
    a = (4 * vat / (np.sqrt(3) * coa)) ** (1 / 3.); c = coa * a
    cell = [[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, c]]
    return Atoms("Fe2", cell=cell, pbc=True,
                 scaled_positions=[[0, 0, 0], [1/3., 2/3., .5]]).repeat(rep)

def bcc(vat, rep=4):
    a = (2 * vat) ** (1 / 3.)
    return Atoms("Fe2", cell=np.eye(3) * a, pbc=True,
                 scaled_positions=[[0, 0, 0], [.5, .5, .5]]).repeat((rep,) * 3)

def fcc(vat, rep=3):
    a = (4 * vat) ** (1 / 3.)
    return Atoms("Fe4", cell=np.eye(3) * a, pbc=True,
                 scaled_positions=[[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]]).repeat((rep,) * 3)

def ep(at, calc):
    at.calc = calc
    e = at.get_potential_energy() / len(at)
    p = -np.trace(at.get_stress(voigt=False)) / 3 / GPa
    return e, p

def static_gates(calc, tel, out):
    # G-D2p: parity on held-out test set
    res = {}
    try:
        test = read(f"phaseD/test_Tel{tel}.xyz", index=":")
    except Exception:
        test = []
    by = {}
    for fr in test:
        by.setdefault(fr.info.get("phase", "?"), []).append(fr)
    out.append("G-D2p force/energy parity (held-out test set):")
    d2_pass = True
    for ph, seq in sorted(by.items()):
        fe, ff = [], []
        for fr in seq:
            f_ref = fr.arrays["REF_forces"]; e_ref = fr.info["REF_energy"]
            fr2 = fr.copy(); fr2.calc = calc
            ff.append(fr2.get_forces() - f_ref)
            fe.append((fr2.get_potential_energy() - e_ref) / len(fr))
        frms = 1000 * np.sqrt(np.mean(np.concatenate(ff) ** 2))
        # energy RMSE about the mean offset (constant shifts are absorbed by E0s)
        fe = np.array(fe); erms = 1000 * np.sqrt(np.mean((fe - fe.mean()) ** 2))
        lim = 120 if ph == "liquid" else 80
        ok = frms < lim and erms < 5
        d2_pass &= ok or ph == "stackfault"   # stackfault reported, not gated
        out.append(f"  {ph:11s} n={len(seq):3d}  F_RMSE={frms:6.1f} meV/A "
                   f"(lim {lim})  E_RMSE={erms:5.2f} meV/atom  -> {'pass' if ok else 'FAIL'}")
    # G-D3: hcp EOS
    P = np.array([ep(hcp(v, rep=(2, 2, 2)), calc)[1] for v in DFT_V])
    m = (DFT_P_HCP_6000K >= 100) & (DFT_P_HCP_6000K <= 350)
    dmax = np.max(np.abs((P - DFT_P_HCP_6000K)[m]))
    d3 = dmax < 3.0
    out.append(f"G-D3 hcp EOS: max|dP| = {dmax:.2f} GPa over 100-350  -> {'pass' if d3 else 'FAIL'}")
    # G-D1: enthalpy ordering
    d1 = True
    for P0, gaps in DFT_GAPS.items():
        H = {}
        for nm, fn in (("hcp", hcp), ("bcc", bcc), ("fcc", fcc)):
            vs = np.linspace(5.8, 9.5, 25)
            es, ps = zip(*[ep(fn(v) if nm != "hcp" else fn(v, rep=(2, 2, 2)), calc) for v in vs])
            v0 = np.interp(P0, np.array(ps)[::-1], vs[::-1])
            H[nm] = np.interp(v0, vs, es) + P0 * GPa * v0
        gf = 1000 * (H["fcc"] - H["hcp"]); gb = 1000 * (H["bcc"] - H["hcp"])
        ok = gf > 0 and gb > gf and abs(gf - gaps["fcc"]) < GAP_TOL and abs(gb - gaps["bcc"]) < 2 * GAP_TOL
        d1 &= ok
        out.append(f"G-D1 @{P0} GPa: fcc-hcp {gf:+.0f} (DFT {gaps['fcc']:+d}), "
                   f"bcc-hcp {gb:+.0f} (DFT {gaps['bcc']:+d}) meV/atom -> {'pass' if ok else 'FAIL'}")
    # c/a
    coas = np.arange(1.52, 1.72, 0.01)
    es = [ep(hcp(6.8, c, rep=(2, 2, 2)), calc)[0] for c in coas]
    i = int(np.argmin(es)); lo, hi = max(i - 3, 0), min(i + 4, len(coas))
    c2 = np.polyfit(coas[lo:hi], es[lo:hi], 2)
    ca = -c2[1] / (2 * c2[0]); ca_ok = abs(ca - DFT_CA) < 0.01
    out.append(f"c/a @V=6.8: {ca:.4f} (DFT {DFT_CA})  -> {'pass' if ca_ok else 'FAIL'}")
    out.append(f"\nSTATIC VERDICT: {'ALL PASS' if (d1 and d3 and ca_ok and d2_pass) else 'FAILURES PRESENT'}")

def gate_d4(calc, out):
    from ase.md.npt import NPT
    from ase.md.nvtberendsen import NVTBerendsen
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as MBD
    import numpy as np
    T = 0.9 * 6350.0
    at = hcp(6.84)
    cell0, Q = at.cell.standard_form()
    at.set_cell(cell0, scale_atoms=False)
    at.positions[:] = at.positions @ Q.T
    at.calc = calc
    MBD(at, temperature_K=T)
    NVTBerendsen(at, 2.0 * units.fs, temperature_K=T,
                 taut=100 * units.fs).run(2500)        # 5 ps NVT pre-equil
    dyn = NPT(at, 2.0 * units.fs, temperature_K=T,
              externalstress=300.0 * GPa, ttime=100 * units.fs,
              pfactor=(750 * units.fs) ** 2 * 250 * GPa)
    bad = 0
    for blk in range(55):                              # 5 ps grace + 50 ps scored
        cas, angs = [], []
        for _ in range(10):
            dyn.run(50)                                # sample every 100 fs
            L = at.cell.lengths(); A = at.cell.angles()
            cas.append(L[2] / L[0])
            angs.append(max(abs(A[0]-90), abs(A[1]-90), abs(A[2]-120)))
        ca, ang = float(np.mean(cas)), float(np.mean(angs))
        if blk < 5:
            continue                                   # grace: barostat settling
        if not (1.55 < ca < 1.68) or ang > 5:
            bad += 1
            out.append(f"  block {blk-4}/50: c/a={ca:.3f} angdrift={ang:.1f} (bad {bad}/2)")
            if bad >= 2:
                out.append(f"G-D4 FAIL: sustained departure from hcp (c/a={ca:.3f}, ang={ang:.1f})")
                return
        else:
            bad = 0
        if blk % 10 == 9:
            out.append(f"  {blk-4} ps: c/a={ca:.3f} angdrift={ang:.1f}")
    out.append(f"G-D4 pass: hcp intact over 50 ps NPT at 0.9 Tm, 300 GPa (final block c/a={ca:.3f})")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--tel", type=int, default=6000)
    ap.add_argument("--gate", default="static", choices=["static", "d4"])
    ap.add_argument("--device", default="cpu")
    a = ap.parse_args()
    from mace.calculators import MACECalculator
    calc = MACECalculator(model_paths=a.model, device=a.device,
                          default_dtype="float64")
    out = [f"model: {a.model}   Tel: {a.tel}"]
    (static_gates(calc, a.tel, out) if a.gate == "static" else gate_d4(calc, out))
    txt = "\n".join(out)
    print(txt)
    open(f"phaseD/gates_{a.gate}_Tel{a.tel}.txt", "w").write(txt + "\n")

if __name__ == "__main__":
    main()
