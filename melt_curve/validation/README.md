# melt_curve/validation/ — accuracy of the production potential against PBE and against experiment

Checks of `fe_tel_e3_stagetwo.model` (T_el-conditioned MACE) made in October 2026, outside
the training pipeline. Every number is written to `ledger.json` by `ledger.py`.

## Producers

| Script | Output | Inputs | Cost |
|---|---|---|---|
| `stats_test.py` | `out/test_stats.json` | `parity_tel.npz` (346 held-out frames, from `parity_tel.py`) | seconds |
| `val_dft.py` | `refs/dft/<tag>.json` | `pw.x`, PBE PAW file | 26 runs, 1.5-12 min each on 2 cores |
| `val_mace.py` | `out/static_compare.json` | `refs/dft/*.json`, model file, torch, mace-torch | about 6 min on 2 CPU threads |
| `cc_check.py` | `out/cc_check.json` | `data/coexistence/analysis_live_all_runs_tel.json`, `data/coexistence/cells/*.json`, `data/coexistence/runs_tel/*_ckpt.npz`, `phaseE_q6.py` | about 1 min |
| `pbe_vs_exp.py` | `out/pbe_vs_exp.json` | `refs/eos_pbe_cold/*.json`, `data/coexistence/pn_correct.json` (optional) | seconds |
| `train_table.py` | `out/training.json` | `production_campaign/arc_scripts/logs/fe_tel_e3_run-7.log`, `production_campaign/arc_scripts/logs/fe_tel_e3_run-7_debug.log`, `data/training/dataset_report.txt`, `data/training/train_e3.xyz`, the model file (class name) | 10 s |
| `tel_check.py` | `out/tel_check.json` | `production_campaign/arc_scripts/logs/e3_test_final.out` (same-frame T_el check) | seconds |
| `coex_cells.py` | `out/coex_cells.json` | `data/coexistence/analysis_live_all_runs_tel.json` | seconds |
| `melting_point.py` | `out/melting_point.json` | `data/coexistence/pn_correct.json` (Simon fit, refitted for its covariance) | seconds |
| `repro_check.py` | `out/repro_check.json` | the model, `data/training/test_Tel*.xyz`, `data/training/parity_tel.npz`; torch, mace-torch | about 8 min on 4 CPU threads |
| `conv_check.py` | `out/convergence.json` | `refs/convergence/*.json` (QE 6.6, hcp at 7.0 A^3/atom) | seconds |
| `lda_setup_check.py` | `out/lda_setup.json` | `refs/lda_setup/*.json` (QE 6.6, hcp and fcc at 7.0 A^3/atom) | seconds |
| `ldashift_summary.py` | `out/ldashift.json` | `studies/lda_shift/ldashift_result.json`, `md_meta.json` (`ldashift.py eval` on ARC) | seconds |
| `settings.py` | `out/settings.json` | `refs/*/*.in`, the producers' sources, `production_campaign/fe_melt/pseudo/*.UPF`, `data/coexistence/cells/*.json`, `data/coexistence/runs_tel/` file names, `production_campaign/fe_melt/labels/ldaX_*`, `production_campaign/fe_melt/ldashift.py`, `refs/lda_setup/*.json` | seconds |
| `ledger.py` | `ledger.json` | `out/*.json`, `refs/literature.json` (committed constants) | seconds |
| `record_inputs.py` | `ledger_inputs.tsv` | the committed inputs that no producer here rebuilds | seconds |

Paths: `paths.py`, a thin wrapper over `../layout.py` (paths relative to the repository root; each
overridable by an environment variable).

`report_numbers.tsv` lists every number of the October 2026 validation report (text at rev 47,
29 Oct 2026) with the ledger key it comes from (columns: number as printed, unit, key, ledger value; `keys:` for producer outputs,
`summary:` for derived numbers; an empty key means no producer). `tests/test_repository.py` checks it
against `ledger.json`.

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
- **Basis and k-point convergence** (`conv_check.py`): static PBE pressure of hcp at 7.0 A^3/atom
  (c/a 1.60, Marzari-Vanderbilt 0.01 Ry). 100/1,000 -> 130/1,300 Ry changes it by -0.53 GPa;
  k 12x12x8 -> 16x16x10 by +0.24 GPa.
- **Mixed LDA setup** (`lda_setup_check.py`): LDA (PZ) on the PBE PAW file, as in the snapshot test,
  against the generated PZ file in perfect lattices at 7.0 A^3/atom: hcp pressure +0.22 GPa,
  fcc - hcp energy -0.14 meV/atom (2 K at dS_m = 0.80 k_B).
- **330 GPa snapshot test** (`ldashift_summary.py`, from `ldashift.py eval`): 20 hcp and 20 liquid
  snapshots of 128 atoms at 330 GPa and 6,329 K, labelled with PBE and LDA (QE 7.5). The potential
  against PBE gives the pressure error, the force RMSE and the liquid - hcp energy error as a T_m error
  (first and second order). LDA - PBE gives the T_m shift by free-energy perturbation (a1 + a2 + a3),
  with an equal-bulk-modulus variant of a3. The ledger derives T_m(PBE, 330 GPa) = T_m(330 GPa) minus
  the potential's error.
- **PBE against experiment** (`pbe_vs_exp.py`): static non-magnetic PBE (k 16x16x10,
  Marzari-Vanderbilt 0.01 Ry) plus a 2.5-3.5 GPa vibrational pressure (Debye-Grueneisen estimate,
  +-1 GPa), minus the 300 K EOS of Dewaele et al. 2006, Fei et al. 2016 and Hirao et al. 2022.
  T_m(true) - T_m(PBE) = +(dT_m/dP)(P_PBE - P_exp) at 7.0 A^3/atom.

## Conventions

- Model input: upstream mace-torch 0.3.16 needs `atoms.info["elec_temp"] = T_el / 1e4` and
  `info_keys={"elec_temp": "elec_temp", ...}`. The checkpoint stores the scale (`elec_temp` embedding
  spec, `scale = 1e-4`); the tel-mace fork applies it, upstream does not. float32 reproduces the stored
  ARC predictions to 0.0 meV/atom and a few 1e-5 eV/A (`out/repro_check.json`). See `load_model.py`.
- Energies are Mermin free energies at T_el (the training target).
- QE stress: QE prints the stress with compression positive; physical stress = -sigma_qe;
  pressure P = trace(sigma_qe) / 3.

The committed `refs/*/*.in` files are the inputs as run, with `pseudo_dir` and `outdir`
rewritten to `./pseudo` and `./tmp` (they pointed at directories of the machine that ran them).

## Tests

`python -m pytest melt_curve/validation/tests -q`, or the whole repository from its root. Producers whose
inputs are absent skip. The model test runs `val_mace.py` on CPU (CUDA hidden) and compares with
`out/static_compare.json` to 2e-4 relative. On a CUDA machine a second test runs it on the GPU and
allows max(2e-4 |b|, 2e-3), with smaller absolute tolerances (`GPU_TIGHT`) for the values whose report
numbers are printed finer than that; `test_gpu_tolerance_below_printed_precision` checks that no report
number taken or derived from `out/static_compare.json` can move by half its last printed digit within
those tolerances. `out/static_compare.json` sections `tel_response`, `fcc_P_*`, `bcc_P_*` and
`coa_min` were added on 26 Oct 2026 from a run on a Mac (CPU, float32); the other values are from the
20-23 Oct run, which that run reproduces to 9e-5 relative. `stats_test.py` matches pytest's
`*_test.py` pattern; the root `conftest.py` keeps pytest from importing (and so running) it.
