#!/usr/bin/env python3
"""
parity_tel.py - MACE (fe_tel_e3_stagetwo) against DFT on the held-out test labels.

Test sets: phaseD/test_Tel4000.xyz, test_Tel6000.xyz, test_Tel8000.xyz (346 DFT calculations,
128-atom cells, every phase). Each frame is evaluated at its own electronic temperature
(atoms.info["elec_temp"]), the same way the production MD sets it.
Writes phaseD/parity_tel.npz (per-atom force components, per-frame energies, phase, T_el, P_est)
and prints per-phase force and energy RMSE for each T_el.
Run in ~/tel_env:   python -u parity_tel.py [--device cpu|cuda]
"""
import argparse
import numpy as np
from ase.io import read

ap = argparse.ArgumentParser()
ap.add_argument("--device", default="cpu")
ap.add_argument("--model", default="phaseD/models/fe_tel_e3_stagetwo.model")
a = ap.parse_args()

from mace.calculators import MACECalculator
calc = MACECalculator(model_paths=a.model, device=a.device, default_dtype="float32")

F_ref, F_mod, f_phase, f_tel = [], [], [], []
E_ref, E_mod, e_phase, e_tel, e_P, e_n = [], [], [], [], [], []
for tel in (4000, 6000, 8000):
    frames = read(f"phaseD/test_Tel{tel}.xyz", ":")
    print(f"T_el {tel} K: {len(frames)} frames", flush=True)
    for k, at in enumerate(frames):
        fr = np.array(at.arrays["REF_forces"]); er = float(at.info["REF_energy"]) / len(at)
        at.info["elec_temp"] = float(tel)
        calc.reset()                          # info is not part of ASE's change check
        at.calc = calc
        em = at.get_potential_energy() / len(at); fm = at.get_forces()
        ph = str(at.info.get("phase", "?"))
        F_ref.append(fr); F_mod.append(fm)
        f_phase += [ph] * len(at); f_tel += [tel] * len(at)
        E_ref.append(er); E_mod.append(em); e_phase.append(ph); e_tel.append(tel)
        e_P.append(float(at.info.get("P_est", np.nan))); e_n.append(len(at))
        if k % 20 == 0:
            print(f"  {k:4d}  {ph:10s}  F_rmse {1000 * np.sqrt(np.mean((fm - fr) ** 2)):6.1f} meV/A", flush=True)

F_ref = np.concatenate(F_ref); F_mod = np.concatenate(F_mod)
f_phase = np.array(f_phase); f_tel = np.array(f_tel)
E_ref = np.array(E_ref); E_mod = np.array(E_mod); e_phase = np.array(e_phase); e_tel = np.array(e_tel)
np.savez_compressed("phaseD/parity_tel.npz", F_ref=F_ref, F_mod=F_mod, f_phase=f_phase, f_tel=f_tel,
                    E_ref=E_ref, E_mod=E_mod, e_phase=e_phase, e_tel=e_tel, e_P=np.array(e_P), e_n=np.array(e_n))

print("\nRMSE against DFT (forces: meV/A per component; energy: meV/atom, offset removed per T_el)")
print(f"{'T_el':>6s} {'phase':12s} {'frames':>6s} {'F_RMSE':>8s} {'E_RMSE':>8s} {'E_offset':>9s}")
for tel in (4000, 6000, 8000, None):
    fm = np.ones(len(f_tel), bool) if tel is None else f_tel == tel
    em = np.ones(len(e_tel), bool) if tel is None else e_tel == tel
    off = np.mean(E_mod[em] - E_ref[em])
    for ph in sorted(set(e_phase)) + ["all"]:
        f_sel = fm & ((f_phase == ph) | (ph == "all")); e_sel = em & ((e_phase == ph) | (ph == "all"))
        if e_sel.sum() == 0:
            continue
        frmse = 1000 * np.sqrt(np.mean((F_mod[f_sel] - F_ref[f_sel]) ** 2))
        dE = E_mod[e_sel] - E_ref[e_sel]
        ermse = 1000 * np.sqrt(np.mean((dE - off) ** 2))
        print(f"{('all' if tel is None else tel):>6} {ph:12s} {e_sel.sum():6d} {frmse:8.1f} {ermse:8.2f} {1000 * off:+9.2f}")
print("\nwritten phaseD/parity_tel.npz")
