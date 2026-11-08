# Testing: how to run, what is covered, conditions of use

How to run the tests and what they check.

## Running the tests

    python -m pytest -q -p no:cacheprovider          # from the repository root

The root `conftest.py` limits collection to `melt_curve/validation/tests/` and `tests/`. Tests skip when
their inputs are absent:
- tests needing the LFS files skip until `git lfs pull` has run;
- the model tests skip without torch and mace-torch;
- `tests/test_runs_tel_model.py` and the GPU run of `val_mace.py` need a CUDA GPU.

A full run on a GPU node with the LFS files (ARC H100): 69 passed. On a Mac with the LFS files: 68 passed,
2 skipped (the two GPU tests).

## What the suites cover

- **`melt_curve/validation/tests/`:** the producers reproduce their committed outputs, and the ledger
  matches them; `val_mace.py` also on a GPU, to tolerances below the printed precision of its report numbers.
- **`tests/`:** see [`tests/README.md`](tests/README.md).

## Conditions of use

- **Training range of the potential:** bcc, fcc, hcp, liquid and stacking-faulted iron, at T_el = 4,000,
  6,000 and 8,000 K. The coexistence cells cover 100-350 GPa.
- **T_el input:** upstream mace-torch needs `atoms.info["elec_temp"] = T_el / 1e4`; the tel-mace fork takes
  T_el in kelvin. See `melt_curve/validation/load_model.py` and Dependencies below.
- **Precision:** the checkpoints are float32. See `melt_curve/validation/README.md`, Conventions.
- **Energies** are Mermin free energies at T_el, the training target.

## Dependencies

- **Analysis, validation, figures and tests:** Python ≥ 3.11 with numpy, scipy, ase ≥ 3.23,
  matplotlib and pytest (`pip install -r requirements.txt`). For the model checks, also torch and
  upstream `mace-torch` 0.3.16. The validation evaluated the production model with upstream 0.3.16:
  `atoms.info["elec_temp"] = T_el / 1e4` (see `melt_curve/validation/load_model.py`).
- **The ARC pipeline (training, MD with the T_el-conditioned model):** the T_el-conditioned MACE fork
  **tel-mace**, a separate repository pinned by commit:

  | | |
  |---|---|
  | remote | `git@github.com:harrytaylor-physics/tel-mace.git` |
  | branch | `tel-conditioning` |
  | commit | `69047bebe27b690b6d52572f5953dde264e4fd11` (22 Sep 2026, "Tel conditioning: elec_temp stamping, scale-in-checkpoint, zero-init, tests") |
  | commit as run | `5168bc0d08bfb0d85ac4504710ee00792d2da04c` (`v0.3.16-114-g5168bc0`); same files, the hash changed when the history was rewritten on 8 Nov 2026 |
  | uncommitted changes | none (ARC, 25 Oct 2026) |
  | installed version | `mace 0.3.17`, editable install in `~/tel_env` (Python 3.13.5, torch 2.14.0+cu130, ase 3.29.0, cuequivariance 0.12.0) |

      pip install git+ssh://git@github.com/harrytaylor-physics/tel-mace.git@69047bebe27b690b6d52572f5953dde264e4fd11

  `5168bc0`, cited in `models/README.md`, `melt_curve/validation/load_model.py`, `MIGRATION.md` and the ARC
  scripts, is the commit as run; use `69047be` to install it.

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
