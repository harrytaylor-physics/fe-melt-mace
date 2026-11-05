# DFT+DMFT pilot (step D0): hcp Fe at inner-core density, TRIQS + Quantum ESPRESSO

The pilot reproduces one published result before the method is used on liquid and
solid snapshots.

| Quantity | Published setup (Pourovskii 2019 review) | Pilot |
|---|---|---|
| Structure | hcp Fe, fixed lattice | hcp Fe, 2 atoms, c/a = 1.60 |
| Volume | 7.05 Å³/atom | 7.05 Å³/atom (a = 2.167 Å, c = 3.467 Å) |
| Temperature | up to 5,800 K | 5,800 K (DFT smearing = k_B T; β = 2.00 eV⁻¹ in DMFT) |
| U, J | 3.4 eV, 0.94 eV (constrained RPA) | same |
| Interaction | density-density | density-density |
| Double counting | around mean-field | around mean-field |
| Solver | CT-QMC | TRIQS/CTHYB |
| Result | Z = 0.7 | pass if 0.65–0.75 |

The DFT part copies the pseudopotential and cutoffs from a campaign label input
(PAW, 16 valence electrons, 100/1,000 Ry), so it matches the melting-curve labels.

## What each step does

- **D0a, environment.** TRIQS 3.3, DFTTools, CTHYB and solid_dmft, installed in a
  conda environment. ARC has no TRIQS module.
- **D0b, DFT and Wannier orbitals.**
  1. QE SCF, then NSCF on a full 10×10×6 grid.
  2. Wannier90 builds 10 Fe 3d orbitals (5 per atom), projected from the Bloch states
     in an energy window, with no further localisation. This matches the projector
     approach of the published work.
  3. Two windows are run: w1 = [E_F − 10, E_F + 4] eV and w2 = [E_F − 12, E_F + 8] eV.
     The difference in Z between them is the projector-choice uncertainty.
  4. The 3s and 3p semicore bands are left out.
- **D0c, checks before DMFT.** The script reports:
  - the 3d occupation per atom and per orbital (the two atoms must be identical);
  - the 3d band centre and band extent;
  - the Wannier spreads.
- **D1 (next).** One-shot DMFT, then Z from the slope of Im Σ(iω_n) at the lowest
  Matsubara frequencies. The input is written once the installed solid_dmft version
  is known.

## Commands

**D0a.** On the login node, once:

```bash
mkdir -p ~/fe_melt/dmft/pilot && cd ~/fe_melt/dmft/pilot && tar xzf ~/fe_melt/dmft_pilot.tgz && \
bash env_install.sh 2>&1 | tee env_install.log
```

**D0b.** Submit the DFT and Wannier90 job (about 30–60 min):

```bash
cd ~/fe_melt/dmft/pilot && sbatch pilot_qe_w90.sbatch
```

**D0c.** When the job has finished:

```bash
cd ~/fe_melt/dmft/pilot && tail -30 D0_*.out && \
module purge && module load Anaconda3 && source "$(conda info --base)/etc/profile.d/conda.sh" && \
conda activate ${DATA:-$HOME}/envs/triqs33 && python wan_check.py
```

Send back `env_install.log`, `~/solid_dmft_params.tgz` and the output of D0c.

## If a Wannier90 window fails

"Energy window contains fewer states than number of target WFs" means that window is
too narrow at some k-point. The other window still runs. Report which one failed.
