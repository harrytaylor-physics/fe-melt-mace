"""stats_test.py - held-out test-set statistics of the production potential.

Input: parity_tel.npz (paths.PARITY), written by parity_tel.py on ARC: per-frame reference and
model energies (eV/atom), forces (eV/A), phase, T_el and pressure of the 346 test frames.
Output: out/test_stats.json - per T_el and phase: energy bias and RMSE after one uniform offset
per T_el, force RMSE, relative force error, max force error; phase-minus-hcp energy errors with
standard errors; energy RMSE in pressure bins.
"""
import json, os
import numpy as np
import paths
d=np.load(paths.PARITY)
dE=(d['E_mod']-d['E_ref'])*1e3          # meV/atom
ph, tel, P = d['e_phase'], d['e_tel'], d['e_P']
dF=d['F_mod']-d['F_ref']; fph, ftel = d['f_phase'], d['f_tel']
out={}
print("Held-out test set: 346 frames x 128 atoms (DFT never seen in training)")
print(f"{'T_el':>5} {'phase':>10} {'n':>3} {'E bias':>8} {'E RMSE*':>8} {'F RMSE':>8} {'F rel':>6} {'F max':>6}")
for T in (4000,6000,8000):
    m=tel==T; off=dE[m].mean()
    for p in ('hcp','liquid','fcc','bcc','stackfault','all'):
        mm = m & ((ph==p) if p!='all' else True)
        fm = (ftel==T) & ((fph==p) if p!='all' else True)
        if mm.sum()==0: continue
        e=dE[mm]-off
        frm=np.sqrt((dF[fm]**2).mean())*1e3; rel=frm/1e3/np.sqrt((d['F_ref'][fm]**2).mean())
        fmax=np.abs(dF[fm]).max()
        out[f'{T}_{p}']=dict(n=int(mm.sum()), bias=float(e.mean()), rmse=float(np.sqrt((e**2).mean())), frmse=float(frm), frel=float(rel), fmax=float(fmax))
        print(f"{T:5d} {p:>10} {mm.sum():3d} {e.mean():+8.2f} {np.sqrt((e**2).mean()):8.2f} {frm:8.1f} {100*rel:5.1f}% {fmax:6.2f}")
print("* energies after removing one uniform offset per T_el (a uniform shift does not affect any free-energy difference)")
print("\nPhase-difference error vs hcp (meV/atom, mean +- s.e.; 1 meV/atom = 14.5 K at dS_m = 0.80 kB):")
for T in (4000,6000,8000):
    m=tel==T; h=dE[m&(ph=='hcp')]
    row=[]
    for p in ('liquid','fcc','bcc','stackfault'):
        x=dE[m&(ph==p)]
        if len(x)<2: continue
        d0=x.mean()-h.mean(); se=np.hypot(x.std(ddof=1)/np.sqrt(len(x)), h.std(ddof=1)/np.sqrt(len(h)))
        out[f'{T}_{p}-hcp']=(float(d0),float(se)); row.append(f"{p}-hcp {d0:+.2f}+-{se:.2f}")
    print(f"  T_el {T}: "+"; ".join(row))
print("\nBy pressure (all T_el, all phases; offset removed per T_el):")
dEc=dE.copy()
for T in (4000,6000,8000): dEc[tel==T]-=dE[tel==T].mean()
for lo,hi in ((50,150),(150,250),(250,320),(320,400)):
    m=(P>=lo)&(P<hi)
    if m.sum():
        rm=float(np.sqrt((dEc[m]**2).mean())); lh=float(dEc[m&(ph=='liquid')].mean()-dEc[m&(ph=='hcp')].mean())
        out[f'Pbin_{lo}-{hi}']=dict(n=int(m.sum()), rmse=rm, liquid_minus_hcp=lh)
        print(f"  {lo:3d}-{hi:3d} GPa: n {m.sum():3d}, E RMSE {rm:.2f} meV/atom, liquid-hcp {lh:+.2f} meV/atom")
print(f"\nP range of test frames: {np.nanmin(P):.0f}-{np.nanmax(P):.0f} GPa; {int(np.isnan(P).sum())} of {len(P)} frames have no stored pressure")
out['P_range_GPa']=[float(np.nanmin(P)), float(np.nanmax(P))]; out['n_without_P']=int(np.isnan(P).sum())
os.makedirs(paths.OUT, exist_ok=True)
json.dump(out, open(os.path.join(paths.OUT,'test_stats.json'),'w'), indent=1)
