# validation/ — accuracy of the production potential against PBE and against experiment

Checks of `fe_tel_e3_stagetwo.model` (T_el-conditioned MACE) made in October 2026, outside
the training pipeline. Every number is written to `ledger.json` by `ledger.py`.

## Producers

| Script | Output | Inputs | Cost |
|---|---|---|---|
| `stats_test.py` | `out/test_stats.json` | `parity_tel.npz` (346 held-out frames, from `parity_tel.py`) | seconds |
| `val_dft.py` | `refs/dft/<tag>.json` | `pw.x`, PBE PAW file | 26 runs, 1.5-12 min each on 2 cores |
| `val_mace.py` | `out/static_compare.json` | `refs/dft/*.json`, model file, torch, mace-torch | about 6 min on 2 CPU threads |
| `cc_check.py` | `out/cc_check.json` | `phaseE/analysis_live_all_runs_tel.json`, `phaseE/cells/*.json`, `phaseE/runs_tel/*_ckpt.npz`, `phaseE_q6.py` | about 1 min |
| `pbe_vs_exp.py` | `out/pbe_vs_exp.json` | `refs/eos_pbe_cold/*.json`, `phaseE/pn_correct.json` (optional) | seconds |
| `ledger.py` | `ledger.json` | `out/*.json` | seconds |

Paths: `paths.py` (relative to the repository root, laid out as `~/fe_melt`; each overridable
by an environment variable of the same name).

## What each check is

- **Held-out test set** (`stats_test.py`): energy and force errors per T_el and phase after one
  uniform energy offset per T_el; phase-minus-hcp mean energy errors with standard errors;
  energy RMSE in pressure bins. 1 meV/atom in the liquid-hcp difference moves T_m by 14.5 K at
  dS_m = 0.80 k_B.
- **Static PBE references** (`val_dft.py`, `val_mace.py`): QE 6.6, pslibrary 1.0.0 PBE PAW
  (16 valence electrons), ecutwfc 100 Ry, ecutrho 1,000 Ry, Fermi-Dirac smearing at T_el, nbnd = 20 per atom.
  - hcp EOS, c/a 1.60, V = 6.6-8.9 A^3/atom, T_el 6,000 K, k 12x12x8
  - T_el = 4,000 and 8,000 K at 7.0 A^3/atom
  - fcc and bcc at 7.0 and 8.0 A^3/atom, k 20x20x20
  - c/a = 1.56, 1.60, 1.64 at 6.8 A^3/atom (quadratic minimum)
  - clamped-ion stress-strain coefficients at 7.0 A^3/atom from +-1% strains (e1 xx, e3 zz,
    e4 yz, e6 xy); C_ij = (sigma+ - sigma-) / 0.02 with the physical stress (= -QE printed stress)
  - E2g phonon at Gamma: the two atoms displaced +-0.03 A along x in opposite directions,
    k = -(F+ - F-) / 0.06, omega = sqrt(k / m_Fe)
  - bulk modulus -V dP/dV from a Vinet fit to the eight EOS pressures; Voigt K, G and
    wave speeds at 13.25 g/cm^3 from the C_ij
- **Clausius-Clapeyron** (`cc_check.py`): dV_m from q-bar-6 slabs of the final coexistence
  configurations, dS_m from the coexistence analysis; slope of T_m(P) between neighbouring
  cells against dV/dS. The dV/dS error includes the dV standard errors only.
- **PBE against experiment** (`pbe_vs_exp.py`): static non-magnetic PBE (k 16x16x10,
  Marzari-Vanderbilt 0.01 Ry) plus a 2.5-3.5 GPa vibrational pressure (Debye-Grueneisen estimate,
  +-1 GPa), minus the 300 K EOS of Dewaele et al. 2006, Fei et al. 2016 and Hirao et al. 2022.
  T_m(true) - T_m(PBE) = +(dT_m/dP)(P_PBE - P_exp) at 7.0 A^3/atom.

## Conventions

- Model input: upstream mace-torch 0.3.16 needs `atoms.info["elec_temp"] = T_el / 1e4` and
  `info_keys={"elec_temp": "elec_temp", ...}`; float32 reproduces the stored ARC predictions
  exactly. See `load_model.py`.
- Energies are Mermin free energies at T_el (the training target).
- QE stress: QE prints the stress with compression positive; physical stress = -sigma_qe;
  pressure P = trace(sigma_qe) / 3.

The committed `refs/*/*.in` files are the inputs as run, with `pseudo_dir` and `outdir`
rewritten to `./pseudo` and `./tmp` (they pointed at directories of the machine that ran them).

## Tests

`python -m pytest validation/tests -q`. Producers whose inputs are absent skip. The model test
runs `val_mace.py` and compares with `out/static_compare.json` to 2e-4 relative.
