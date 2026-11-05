#!/usr/bin/env python3
"""
phaseE_pin.py — interface-pinning cross-check (fe_melt_v1, Phase E method 2)

Interface pinning (Pedersen et al., PRB 88, 094101): run a two-phase cell in
NVT at trial temperature T with a harmonic bias U = (k/2)(Q - a0)^2 on a
global crystallinity order parameter Q. The bias holds the interface in
place; the mean bias force <k(Q - a0)> equals the thermodynamic driving
force mu_l - mu_s per unit Q. Scanning T and finding where the mean force
crosses zero gives Tm at the run's pressure, WITHOUT waiting for slow
interface migration. Complements the NVE method with different systematics
(thermostatted, biased, smaller cell).

Q here: sum of per-atom Steinhardt-like q6 computed in torch with autograd,
so bias forces are exact gradients. Cell: ~4,000 atoms (half the linear
size of the NVE cells).

One run = one (T, a0) point -> mean force. Drive several T per pressure:
  python phaseE_pin.py --model <m> --cell phaseE/cells/coex_hcp_300GPa \
      --temp 6200 [--k 5.0] [--ps 60] [--device cuda]
Aggregation to Tm happens in phaseE_analyze.py.
"""
from phaseE_nlcache import enable_nl_cache
enable_nl_cache(skin=0.0)
import argparse, json, os
import numpy as np
import torch
from ase import units
from ase.io import read
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as MBD
from ase.calculators.calculator import Calculator, all_changes
from ase.calculators.mixing import SumCalculator
from ase.units import GPa

def torch_q6_sum(pos, cell, rcut=3.2):
    """Differentiable global Q = sum_i |q6_i| (simplified: rotationally
    averaged pair-angle kernel P6(cos theta_ijk) contracted; keeps autograd
    tractable at 4k atoms via neighbor sparsity)."""
    inv = torch.linalg.inv(cell)
    d = pos[None, :, :] - pos[:, None, :]
    frac = d @ inv; frac = frac - torch.round(frac); d = frac @ cell
    r = torch.linalg.norm(d, dim=-1)
    n = pos.shape[0]
    mask = (r < rcut) & (r > 0.1)
    # smooth cutoff
    w = torch.where(mask, 0.5 * (1 + torch.cos(np.pi * r.clamp(0.1, rcut) / rcut)), torch.zeros_like(r))
    dn = torch.nn.functional.normalize(torch.where(mask[..., None], d, torch.ones_like(d)), dim=-1)
    # bond-order via Legendre P6 of pairwise bond-vector dot products, per atom
    Q_atoms = []
    idx = torch.nonzero(mask)
    for i in range(n):
        js = idx[idx[:, 0] == i, 1]
        if js.numel() < 4:
            Q_atoms.append(pos.sum() * 0.0)
            continue
        v = dn[i, js]                        # (m,3) unit bond vectors
        wi = w[i, js]
        c = v @ v.T                          # cos(theta) between bonds
        x = c.clamp(-1, 1)
        P6 = (231 * x**6 - 315 * x**4 + 105 * x**2 - 5) / 16.0
        ww = wi[:, None] * wi[None, :]
        q2 = (P6 * ww).sum() / ww.sum()
        Q_atoms.append(torch.sqrt(torch.clamp(q2, min=1e-8)))
    return torch.stack(Q_atoms).sum()

class Q6Bias(Calculator):
    implemented_properties = ["energy", "forces"]
    def __init__(self, k, a0, device="cuda"):
        super().__init__(); self.k, self.a0, self.dev = k, a0, device
        self.last_Q = None
    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        pos = torch.tensor(atoms.get_positions(), dtype=torch.float32,
                           device=self.dev, requires_grad=True)
        cell = torch.tensor(atoms.cell[:].astype(np.float32), device=self.dev)
        Q = torch_q6_sum(pos, cell)
        E = 0.5 * self.k * (Q - self.a0) ** 2
        E.backward()
        self.last_Q = float(Q.detach())
        self.results["energy"] = float(E.detach())
        self.results["forces"] = -pos.grad.detach().cpu().numpy()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--cell", required=True)
    ap.add_argument("--temp", type=float, required=True)
    ap.add_argument("--k", type=float, default=5.0e-4)   # eV per Q^2 (Q ~ N-scaled)
    ap.add_argument("--ps", type=float, default=60.0)
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()
    meta = json.load(open(a.cell + ".json"))
    at = read(a.cell + ".extxyz")
    # reduce to ~4k atoms: keep central z-slab spanning the interface
    z = at.get_positions()[:, 2]; zc = z.mean(); Lz = at.cell.lengths()[2]
    keep = np.abs(z - zc) < Lz * (2000.0 / len(at))       # slab -> ~4000 atoms
    small = at[keep]
    small.cell[2, 2] = small.get_positions()[:, 2].ptp() + 2.2
    small.pbc = True
    from mace.calculators import MACECalculator
    mace = MACECalculator(enable_cueq=True, model_paths=os.path.expanduser(a.model),
                          device=a.device, default_dtype="float32")
    bias = Q6Bias(a.k, a0=None, device=a.device)
    small.calc = mace
    # set a0 to the current Q (interface at half-and-half)
    posT = torch.tensor(small.get_positions(), dtype=torch.float32, device=a.device)
    cellT = torch.tensor(small.cell[:].astype(np.float32), device=a.device)
    with torch.no_grad():
        bias.a0 = float(torch_q6_sum(posT, cellT))
    small.calc = SumCalculator([mace, bias])
    MBD(small, temperature_K=a.temp)
    dyn = Langevin(small, 2.0 * units.fs, temperature_K=a.temp, friction=0.02)
    dyn.run(2500)                                          # 5 ps equil under bias
    forces_Q, Ps = [], []
    nblk = int(a.ps)
    for b in range(nblk):
        dyn.run(500)
        forces_Q.append(a.k * (bias.last_Q - bias.a0))
        Ps.append(-np.trace(small.get_stress(voigt=False)) / 3 / GPa)
        if b % 10 == 0:
            print(f"  {b} ps  <kdQ>={np.mean(forces_Q):+.3e}  P={np.mean(Ps):.1f}", flush=True)
    blocks = np.array_split(np.array(forces_Q[len(forces_Q)//4:]), 6)
    f_mean = float(np.mean([b.mean() for b in blocks]))
    f_err = float(np.std([b.mean() for b in blocks]) / np.sqrt(6))
    out = dict(tag=meta["tag"], method="pinning", T=a.temp,
               P=float(np.mean(Ps[len(Ps)//4:])), k=a.k, a0=bias.a0,
               force=f_mean, force_err=f_err, N=len(small))
    os.makedirs("phaseE/runs", exist_ok=True)
    fn = f"phaseE/runs/pin_{meta['tag']}_T{int(a.temp)}.json"
    json.dump(out, open(fn, "w"), indent=1)
    print(f"[pin] T={a.temp:.0f}: <force>={f_mean:+.3e} +/- {f_err:.1e} -> {fn}")

if __name__ == "__main__":
    main()
