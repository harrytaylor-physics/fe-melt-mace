# production_campaign/

ARC pipeline as run: QE labels, training, coexistence runs.

The pipeline that produced the iron melting curve, as it ran on ARC: sampling, QE labels, training,
coexistence and the snapshot test. **Kept unmodified, as provenance.** Paths inside the scripts are those of
`~/fe_melt` on ARC and are left as they ran. Slurm headers use the `phys-plasma-dft` account,
`cd ~/fe_melt`, and `~/tel_env` for the tel-mace fork.

    fe_melt/        the ARC working directory ~/fe_melt: Python drivers and analysis scripts, plus
                    labels/, pseudo/, traj/ and results/ at their ARC relative paths
    arc_scripts/    the Slurm submission scripts, stage_fe_melt.sh (the read-only ARC snapshot script),
                    and logs/, every log of the campaign

## ARC paths and this repository

The scripts run from `~/fe_melt`, not from a clone of this repository. To rerun one from a clone,
recreate the tree on the left.

| In the scripts and logs | Here |
|---|---|
| `*.py`, `labels/`, `pseudo/`, `traj/`, `results/` | `production_campaign/fe_melt/` |
| `*.sbatch`, `logs/` | `production_campaign/arc_scripts/` |
| `phaseD/models/` | `models/` |
| `phaseD/` (other files) | `data/training/` |
| `phaseE/` | `data/coexistence/` |
| `phaseX/ldashift/` | `studies/lda_shift/` |
| `selected/`, `labels/labels.jsonl`, QE outputs, trajectories | not in the repository; see `data/README.md` |

## How each stage ran

| Stage | Command (ARC, in `~/fe_melt`, unless noted) |
|---|---|
| Sampling | `sbatch fe_sampleC.sbatch` (array over phases) |
| QE labels | `phaseC_qe_labels.py` writes `labels/*.in`; `labels/run_*.sbatch` run them; `labels/collect.py` collects them |
| Datasets | `python phaseD_dataset.py` (needs `selected/` and `labels/labels.jsonl`, which stay on ARC) |
| Training | `sbatch phaseD_train_e3.sbatch` (production); `phaseD_train_f32.sbatch` with `TEL=4000/6000/8000` (v1) |
| Gates | `python phaseD_gates.py --model phaseD/models/<m>.model --tel <T>` |
| Coexistence | `python phaseE_build.py ...` (`build_cells.sbatch`); `sbatch fe_phaseE_tel.sbatch` (48-task array, resubmits itself) |
| Analysis | `python phaseE_analyze_live.py --rundir phaseE/runs_tel --include-drifting`; `python pn_correct.py` (any machine) |
| Held-out predictions | `python parity_tel.py --device cuda` (`~/tel_env`) |
| 330 GPa snapshot test | `sbatch ldashift_md.sbatch` (GPU; writes `phaseX/ldashift/frames.xyz` and the inputs); `cd labels && sbatch ldaXc.sbatch` (one input per task, any number of tasks); `python ldashift.py eval` |

Which model each run loaded is in [`../models/README.md`](../models/README.md).
