# Testing: how to run, what is covered, conditions of use

## Running the tests

    python -m pytest -q -p no:cacheprovider          # from the repository root

The root `conftest.py` limits collection to `melt_curve/validation/tests/` and `tests/`. Tests skip when
their inputs are absent:
- tests needing the LFS files skip until `git lfs pull` has run;
- the model tests skip without torch and mace-torch;
- `tests/test_runs_tel_model.py` and the GPU run of `val_mace.py` need a CUDA GPU.

A full run on a GPU node with the LFS files (ARC H100): 69 passed. On a Mac with the LFS files: 67 passed,
2 skipped (the two GPU tests).

## What the suites cover

- **`melt_curve/validation/tests/`:** the producers reproduce their committed outputs, and the ledger
  matches them; `val_mace.py` also on a GPU, to tolerances below the printed precision of its report numbers.
- **`tests/`:** see [`tests/README.md`](tests/README.md).

## Conditions of use

- **Training range of the potential:** bcc, fcc, hcp, liquid and stacking-faulted iron, at T_el = 4,000,
  6,000 and 8,000 K. The coexistence cells cover 100-350 GPa.
- **T_el input:** upstream mace-torch needs `atoms.info["elec_temp"] = T_el / 1e4`; the tel-mace fork takes
  T_el in kelvin. See `melt_curve/validation/load_model.py` and [`README.md`](README.md).
- **Precision:** the checkpoints are float32. See `melt_curve/validation/README.md`, Conventions.
- **Energies** are Mermin free energies at T_el, the training target.
