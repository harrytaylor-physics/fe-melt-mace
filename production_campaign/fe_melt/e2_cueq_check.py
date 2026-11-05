#!/usr/bin/env python3
"""
e2_cueq_check.py — can the Tel-conditioned model run production-size MD?
(tel_env, GPU node)

The 12,096-atom coexistence cells fit on ARC's GPUs only with cuEquivariance
(~1.85 MB/atom with it, ~8.7 MB/atom without). The Tel embedding has never been
run with cuEquivariance. Three checks:

  1. Identity: one 128-atom test frame with cuEquivariance off and on.
     Energies and forces must agree to float32 precision.
  2. Tel response: the same frame at Tel = 4000, 6000, 8000 K. The Mermin free
     energy must fall as Tel rises (F = E - Tel*S_el, S_el > 0).
  3. Memory: the 12,096-atom 300 GPa coexistence cell at Tel = 6000 K with
     cuEquivariance, one force + stress evaluation, peak GPU memory.
     The fixed-Tel production model measured 21.6 GiB here (skin 0).

Run: ~/tel_env/bin/python e2_cueq_check.py
"""
import numpy as np
import torch
from ase.io import read
from mace.calculators import MACECalculator

M = "phaseD/models/fe_tel_e2.model"
frame = read("phaseD/test_Tel6000_stamped.xyz", index=0)

print("1. cuEquivariance identity (128-atom frame, Tel = 6000 K)")
res = {}
for cq in (False, True):
    try:
        c = MACECalculator(model_paths=M, device="cuda", default_dtype="float32", enable_cueq=cq)
        a = frame.copy(); a.info["elec_temp"] = 6000.0; a.calc = c
        res[cq] = (a.get_potential_energy() / len(a), a.get_forces())
        print(f"   cueq={cq!s:5s}: E = {res[cq][0]:.6f} eV/atom")
    except Exception as e:
        print(f"   cueq={cq!s:5s}: FAILED  {type(e).__name__}: {e}")
if len(res) == 2:
    dE = 1000 * abs(res[True][0] - res[False][0])
    dF = 1000 * np.abs(res[True][1] - res[False][1]).max()
    ok = dE < 0.05 and dF < 5.0
    print(f"   difference: dE = {dE:.4f} meV/atom, max dF = {dF:.3f} meV/A  -> {'ok' if ok else 'MISMATCH'}")

print("\n2. Tel response (same frame)")
c = MACECalculator(model_paths=M, device="cuda", default_dtype="float32",
                   enable_cueq=(True in res))
E = {}
for T in (4000.0, 6000.0, 8000.0):
    a = frame.copy(); a.info["elec_temp"] = T; a.calc = c
    E[T] = a.get_potential_energy() / len(a)
    print(f"   Tel {T:6.0f} K: F = {E[T]:.5f} eV/atom")
mono = E[4000.0] > E[6000.0] > E[8000.0]
curv = (E[4000.0] + E[8000.0]) / 2 - E[6000.0]
print(f"   falls with Tel: {'yes' if mono else 'NO'};  "
      f"midpoint curvature (mean of ends minus middle) = {1000*curv:+.1f} meV/atom "
      f"(negative expected: F is concave in Tel)")

print("\n3. memory, 12,096-atom 300 GPa cell, Tel = 6000 K")
if True not in res:
    print("   skipped: cuEquivariance unavailable for this model")
else:
    at = read("phaseE/cells/coex_hcp_300GPa.extxyz")
    at.info["elec_temp"] = 6000.0
    at.calc = c
    torch.cuda.reset_peak_memory_stats()
    try:
        at.get_forces(); at.get_stress()
        print(f"   peak {torch.cuda.max_memory_allocated()/2**30:.1f} GiB on "
              f"{torch.cuda.get_device_name(0)} (fixed-Tel model: 21.6 GiB)")
    except torch.OutOfMemoryError:
        print("   OUT OF MEMORY")
