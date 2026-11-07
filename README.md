# fe-melt-mace

The melting curve of iron to 370 GPa from two-phase coexistence simulations with a MACE potential
trained on PBE (Mermin free energies) whose input includes the electronic temperature T_el.

**Result** (keys of [`melt_curve/validation/ledger.json`](melt_curve/validation/ledger.json)):
- **T_m(330 GPa) = 6,329 ± 6 K from the potential.** The value is `keys.melting_point.Tm_330`, from the
  Simon fit in `data/coexistence/pn_correct.json`. The error is statistical, the fit's standard error at
  330 GPa: `summary.Tm_330_stat_err_K`. dT_m/dP = 8.41 K/GPa (`keys.melting_point.dTdP_330`).
- **T_m(330 GPa) = 6,306 ± 8 K for PBE.** A direct PBE test on 40 snapshots at 330 GPa and 6,329 K puts the
  potential's T_m 23 ± 6 K above PBE (`keys.ldashift.mace_Tm_total_K`). The PBE value is
  `summary.Tm_PBE_330_K`. Its error, `summary.Tm_PBE_330_err_K`, is the Simon-fit error and the
  direct test's error added in quadrature.

Every number quoted from this work has a ledger key; the list is
[`melt_curve/validation/report_numbers.tsv`](melt_curve/validation/report_numbers.tsv).

## What you can do with it

- Evaluate the T_el-conditioned iron potential at a chosen T_el (`melt_curve/validation/load_model.py`)
- Reproduce the Simon fit from the coexistence analysis
- Rerun the checks of the potential against PBE and of PBE against experiment
- Rebuild the ledger and the figures from the committed inputs

## Quick start

    git lfs install                  # required once; models and datasets are LFS objects
    git clone git@github.com:harrytaylor-physics/fe-melt-mace.git
    cd fe-melt-mace
    pip install -r requirements.txt
    python -m pytest -q -p no:cacheprovider

Without the LFS files the binaries are 132-byte pointers and the tests that need them skip. To fetch them:
`git lfs install --local && git lfs pull`. Conditions of use are in [TESTING.md](TESTING.md).

## Repository layout

| Where | What |
|---|---|
| [models/](models/README.md) | MACE potentials with T_el input, and their provenance. |
| [melt_curve/](melt_curve/README.md) | Validation against PBE and experiment, ledger and figures. |
| [production_campaign/](production_campaign/README.md) | ARC pipeline as run: QE labels, training, coexistence runs. |
| [data/](data/README.md) | Training and test sets with QE labels, and coexistence runs. |
| [studies/](studies/README.md) | Side studies: the 330 GPa LDA - PBE test, a DFT+DMFT pilot. |
| [tests/](tests/README.md) | Test suite for the producers, the ledger and the layout. |
| [docs/](docs/README.md) | Methods reference and figures of the iron melting curve. |
| [MIGRATION.md](MIGRATION.md) · [CHANGELOG.md](CHANGELOG.md) | how `~/fe_melt` came in and what moved where; what was done when |
| [TESTING.md](TESTING.md) | how to run the tests, what they cover, conditions of use |

## Status

| Component | Status |
|---|---|
| `fe_tel_e3_stagetwo.model` (T_el-conditioned potential) | production; SWA stage two; trained at T_el = 4,000, 6,000 and 8,000 K |
| `fe_ft_Tel{4000,6000,8000}_f32.model` (fixed-T_el potentials) | v1; the 6,000 K model ran the campaign in `data/coexistence/runs/`, which set each production cell's T_el |
| Melting curve, 100-350 GPa | complete; 12,096-atom cells, 6 NVE runs of 200 ps per cell |
| Validation against PBE and experiment | complete; every number in the ledger |
| 330 GPa snapshot test | complete ([studies/lda_shift/](studies/lda_shift/README.md)) |
| EOS anchor (`eos_anchor.py`) | not used; its 28 runs were cancelled |
| DFT+DMFT pilot | exploratory; not used in any result ([studies/dmft_pilot/](studies/dmft_pilot/README.md)) |

## Methods in one paragraph

MD sampling of bcc, fcc, hcp, liquid and stacking-faulted iron with a MACE foundation model; PBE labels in
Quantum ESPRESSO at T_el = 4,000, 6,000 and 8,000 K (Fermi-Dirac smearing, pslibrary 1.0.0 PAW with 16
valence electrons); fine-tuning of the foundation model with T_el as a learned per-structure embedding;
12,096-atom hcp-liquid coexistence at 100-350 GPa, each cell at the T_el of its melting temperature;
interface-normal pressure correction and a Simon fit. Full details: [docs/METHODS.md](docs/METHODS.md).

## Dependencies

- **Analysis, validation, figures and tests:** Python ≥ 3.11 with numpy, scipy, ase ≥ 3.23,
  matplotlib and pytest (`pip install -r requirements.txt`). For the model checks, also torch and
  upstream `mace-torch` 0.3.16. The validation evaluated the production model with upstream 0.3.16:
  `atoms.info["elec_temp"] = T_el / 1e4` (see `melt_curve/validation/load_model.py`).
- **The ARC pipeline (training, MD with the T_el-conditioned model):** the T_el-conditioned MACE fork
  **tel-mace**, a separate private repository pinned by commit:

  | | |
  |---|---|
  | remote | `git@github.com:harrytaylor-physics/tel-mace.git` |
  | branch | `tel-conditioning` |
  | commit | `5168bc0d08bfb0d85ac4504710ee00792d2da04c` (22 Sep 2026, "Tel conditioning: elec_temp stamping, scale-in-checkpoint, zero-init, tests"; `v0.3.16-114-g5168bc0`) |
  | uncommitted changes | none (ARC, 25 Oct 2026) |
  | installed version | `mace 0.3.17`, editable install in `~/tel_env` (Python 3.13.5, torch 2.14.0+cu130, ase 3.29.0, cuequivariance 0.12.0) |

      pip install git+ssh://git@github.com/harrytaylor-physics/tel-mace.git@5168bc0d08bfb0d85ac4504710ee00792d2da04c

  The checkpoint stores the T_el input scale: the `elec_temp` embedding spec has `scale = 1e-4`. The
  fork applies it inside the model, so the ARC scripts pass T_el in kelvin. Upstream mace-torch does not
  apply the stored scale, so callers pass T_el / 1e4. With the fork installed, set `FE_TEL_SCALE=1` for
  `melt_curve/validation/` and `tests/`.
- **QE:** 7.5 on ARC for the labels (`provenance/qe_version.txt`); 6.6 for the validation references.
- **Foundation model:** `MACE-matpes-pbe-omat-ft.model`. It is not in this repository; its source and
  sha256 are in `models/README.md`.
- **LFS:** models, coexistence checkpoints, the training and test sets, `inputs_as_run.tgz`.

## Reproducing

These run on any machine, from the repository root. The ARC stages are in
[production_campaign/README.md](production_campaign/README.md).

| Stage | Command |
|---|---|
| Simon fit | see [melt_curve/README.md](melt_curve/README.md) |
| Validation producers and ledger | `cd melt_curve/validation && python stats_test.py && python cc_check.py && python pbe_vs_exp.py && python train_table.py && python tel_check.py && python coex_cells.py && python melting_point.py && python conv_check.py && python lda_setup_check.py && python ldashift_summary.py && python settings.py && python ledger.py` (`val_mace.py` and `repro_check.py` need mace-torch, `val_dft.py` needs `pw.x`) |
| Figures | `python melt_curve/figures/build_figures.py` (about 2 min) |
| Tests | `python -m pytest -q -p no:cacheprovider` |

---
Maintainer: Harry Taylor (Oxford).
