# fe-melt-mace

The melting curve of iron to 370 GPa from two-phase coexistence simulations with a MACE potential
trained on PBE (Mermin free energies) whose input includes the electronic temperature T_el.

**Result:** T_m(330 GPa) = 6,329 K, with dT_m/dP = 8.41 K/GPa at 330 GPa. These are ledger keys
`melting_point.Tm_330` and `melting_point.dTdP_330` in [`validation/ledger.json`](validation/ledger.json),
from the Simon fit in `phaseE/pn_correct.json`. Every number quoted from this work has a ledger key;
the list is [`validation/report_numbers.tsv`](validation/report_numbers.tsv).

## What the pipeline does

1. **QE labels (Phase C).** MD sampling of bcc, fcc, hcp, liquid and stacking-faulted iron with a MACE
   foundation model. The selected 128-atom frames (108 for fcc) are labelled with PBE in Quantum
   ESPRESSO at T_el = 4,000, 6,000 and 8,000 K, using Fermi-Dirac smearing and pslibrary 1.0.0 PAW
   with 16 valence electrons. Scripts: `phaseC_*.py`, `labels/`.
2. **Training (Phase D).** Datasets (`phaseD_dataset.py`), then fine-tuning of the foundation model
   (`phaseD_train*.sbatch`) and gates (`phaseD_gates.py`). The production potential
   `phaseD/models/fe_tel_e3_stagetwo.model` takes T_el as a learned per-structure embedding. It was
   trained on all three T_el (`phaseD/train_e3.xyz`) and saved after the SWA stage two. Provenance:
   [`models/README.md`](models/README.md).
3. **Coexistence (Phase E).** 12,096-atom hcp-liquid cells at 100-350 GPa, 6 NVE runs of 200 ps per
   cell, each cell at the T_el of its melting temperature (`phaseE_build.py`, `phaseE_coex_tel.py`,
   `fe_phaseE_tel.sbatch`, output in `phaseE/runs_tel/`).
4. **Analysis.**
   - Melting temperature and entropy per cell: `phaseE_analyze_live.py` →
     `phaseE/analysis_live_all_runs_tel.json`.
   - Interface-normal pressure and the Simon fit: `pn_correct.py` → `phaseE/pn_correct.json`.
   - T_el decomposition: `decomp_tel.py`.
   - Validation against PBE and experiment, and the ledger: [`validation/`](validation/README.md).
   - Figures: [`docs/figures/`](docs/figures/README.md).

The LDA − PBE test at 330 GPa (`ldashift.py`, `labels/ldaX_*`) and the EOS anchor (`eos_anchor.py`)
were running on ARC when this repository was made; their results come in a later change.

## Dependencies

- **Analysis, validation, figures and tests:** Python ≥ 3.11 with numpy, scipy, ase ≥ 3.23,
  matplotlib and pytest (`pip install -r requirements.txt`). For the model checks, also torch and
  upstream `mace-torch` 0.3.16. The validation evaluated the production model with upstream 0.3.16:
  `atoms.info["elec_temp"] = T_el / 1e4` (see `validation/load_model.py`).
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
  `validation/` and `tests/`.
- **QE:** 7.5 on ARC for the labels (`provenance/qe_version.txt`); 6.6 for the validation references.
- **Foundation model:** `MACE-matpes-pbe-omat-ft.model`. It is not in this repository; its source and
  sha256 are in `models/README.md`.
- **LFS:** models, coexistence checkpoints, the training and test sets, `labels/inputs_as_run.tgz`.
  After cloning, run `git lfs install --local && git lfs pull`; otherwise those files are 132-byte
  pointers.

## How to run each stage

The ARC scripts run from the repository root on ARC, as they did in `~/fe_melt`. Their Slurm headers
use the `phys-plasma-dft` account, `cd ~/fe_melt`, and `~/tel_env` for the fork.

| Stage | Command (ARC unless noted) |
|---|---|
| Sampling | `sbatch fe_sampleC.sbatch` (array over phases) |
| QE labels | `phaseC_qe_labels.py` writes `labels/*.in`; `labels/run_*.sbatch` run them; `labels/collect.py` collects them |
| Datasets | `python phaseD_dataset.py` (needs `selected/` and `labels/labels.jsonl`, which stay on ARC) |
| Training | `sbatch phaseD_train_e3.sbatch` (production); `phaseD_train_f32.sbatch` with `TEL=4000/6000/8000` (v1) |
| Gates | `python phaseD_gates.py --model phaseD/models/<m>.model --tel <T>` |
| Coexistence | `python phaseE_build.py ...` (`build_cells.sbatch`); `sbatch fe_phaseE_tel.sbatch` (48-task array, resubmits itself) |
| Analysis | `python phaseE_analyze_live.py --rundir phaseE/runs_tel --include-drifting`; `python pn_correct.py` (any machine) |
| Held-out predictions | `python parity_tel.py --device cuda` (`~/tel_env`) |
| Validation producers | `cd validation && python stats_test.py && python cc_check.py && python pbe_vs_exp.py && python train_table.py && python tel_check.py && python coex_cells.py && python melting_point.py && python settings.py && python ledger.py` (any machine; `val_mace.py` and `repro_check.py` need mace-torch, `val_dft.py` needs `pw.x`) |
| Figures | `python analysis/figures/build_figures.py` (any machine, about 2 min) |

## Tests

    python -m pytest -q -p no:cacheprovider          # from the repository root

The root `conftest.py` limits collection to `validation/tests/` and `tests/`. Tests skip when their
inputs are absent:
- tests needing the LFS files skip until `git lfs pull` has run;
- the model tests skip without torch and mace-torch;
- `tests/test_runs_tel_model.py` needs a CUDA GPU.

What the suites cover:
- **`validation/tests/`:** the producers reproduce their committed outputs, and the ledger matches them.
- **`tests/`:**
  - the new producers and `pn_correct.py` reproduce their outputs;
  - the literature constants match the producers that use them;
  - the report-number list matches the ledger;
  - no committed code holds an absolute path;
  - the staleness guard (`validation/ledger_inputs.tsv`);
  - on a GPU, which model the T_el-conditioned runs loaded.

## Layout

| Where | What |
|---|---|
| `*.py`, `*.sbatch` (root) | the ARC pipeline, verbatim |
| [`labels/`](labels/) | QE input templates, the running jobs' inputs, queues; every input as run in `inputs_as_run.tgz` |
| [`pseudo/`](pseudo/) | PAW files; [`lda_generation/`](pseudo/lda_generation/README.md) |
| `phaseD/` | training and test sets, gates, predictions; `models/` (LFS) |
| `phaseE/` | coexistence cells, runs, final configurations (LFS), analysis, Simon fit |
| `phaseX/`, `traj/`, `results/`, `logs/` | LDA-shift frames, sampling metadata, e3 training metrics, every log |
| [`validation/`](validation/README.md) | checks against PBE and experiment, producers, `ledger.json`, `report_numbers.tsv` |
| [`analysis/figures/`](analysis/figures/README.md) | figure producers and their intermediate json |
| [`docs/figures/`](docs/figures/README.md) | figures (PDF + 200-dpi PNG) |
| [`models/README.md`](models/README.md) | which model each run loaded, and how it was made |
| [`data/`](data/README.md) | SHA-256 manifest of every file left out |
| [`exploratory/`](exploratory/README.md) | the DFT+DMFT pilot, not used in any result |
| `tests/`, `conftest.py` | repository tests |
| `layout.py` | every path the analysis, validation and figure code uses |
| `tools/stage_fe_melt.sh` | the read-only ARC snapshot script |
| [MIGRATION.md](MIGRATION.md) · [CHANGELOG.md](CHANGELOG.md) | how `~/fe_melt` came in; what was done when |
