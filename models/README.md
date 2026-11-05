# models/ — the potentials and where they came from

The model files are here (git LFS). On ARC they are in `~/fe_melt/phaseD/models/`, the path the scripts,
logs and gate files name. This file records which file each run loaded and how that file was made. Superseded models are not in the repository;
their sha256 are in [`../data/excluded_sha256.tsv`](../data/excluded_sha256.tsv).

## How MACE names its output, and how files were matched

`mace_run_train` writes `<name>_run-<seed>.model` and `<name>.model` at the end of stage one (the best
validation checkpoint) and `<name>_run-<seed>_stagetwo.model` and `<name>_stagetwo.model` after the SWA
stage. The two files of a pair hold identical weights but different sha256, because the archive stores
the name it was saved under. Files were therefore matched by per-tensor comparison (all tensors equal),
by the saved-as name inside the archive, by file time against the "Saving model" lines of the training
log, and by the model paths in the scripts and logs. Source: `provenance/` of the 25 Oct 2026 ARC
snapshot (`model_files.tsv`, `model_tensor_compare.txt`, `training_log_lines.txt`,
`model_paths_in_scripts.txt`, `model_paths_in_run_records.txt`), made by
[`production_campaign/arc_scripts/stage_fe_melt.sh`](../production_campaign/arc_scripts/stage_fe_melt.sh).

## Production model provenance

| Potential | File | sha256 | Saved as | Stage | Training run and log | Loaded by |
|---|---|---|---|---|---|---|
| **Production**: T_el-conditioned e3 | `models/fe_tel_e3_stagetwo.model` (5,883,607 B, 4 Oct 2026 01:07:00) | `6b327a52f9d075510024091125d3f5b63e16120da891903e9479d9b9b7cf3df5` | `fe_tel_e3_stagetwo` | SWA stage two, epoch 149 | `phaseD_train_e3.sbatch`, jobs 8959136-8959138 (tel-mace 5168bc0, float32, `data/training/train_e3.xyz`); `production_campaign/arc_scripts/logs/train_e3_8959138.out` line 215, `production_campaign/arc_scripts/logs/fe_tel_e3_run-7.log` line 366: "Loaded Stage two model from epoch 149 ... Saving model to phaseD/models/fe_tel_e3_run-7_stagetwo.model" at 01:07:00.321 | `fe_phaseE_tel.sbatch` → `data/coexistence/runs_tel/` (48 runs; see below); `parity_tel.py` → `data/training/parity_tel.npz`; `tel_diff_check.py` → `production_campaign/arc_scripts/logs/e3_test_final.out`; `decomp_tel.py`, `decomp_tel.sbatch` energy check; `spot_eval.py`; `ldashift.py` (running); `melt_curve/validation/` (`paths.MODEL`) |
| e3 stage one | `models/fe_tel_e3.model` (5,882,635 B, 4 Oct 00:53:23) | `cd73202d427e729d5d79fc11e8eb97defc0c8220c1652c498d4fbdb5b58563a1` | `fe_tel_e3` | stage one, best checkpoint epoch 109 | same run; `production_campaign/arc_scripts/logs/fe_tel_e3_run-7.log` line 341 (00:53:23.008) | no run; kept for `tests/test_runs_tel_model.py` |
| v1, T_el = 4,000 K | `models/fe_ft_Tel4000_f32.model` (17 Sep 11:20:58) | `bfdbb165990a1807aae3646b60fd08407225a3e342eb0382aa1d09b09c78cc96` | `fe_ft_Tel4000_f32` | stage one (no SWA stage was run) | `phaseD_train_f32.sbatch`, `phaseD/train_Tel4000.xyz` on ARC; `production_campaign/arc_scripts/logs/train_Tel4000_f32_resume.out` line 111, `production_campaign/arc_scripts/logs/fe_ft_Tel4000_f32_run-7.log` | `data/training/gates_static_Tel4000.txt`, `gates_v1_Tel4000_static.txt`; `phaseE_reweight.py`; `e2_control.py`; `tel_diff_check.py` (reference) |
| v1, T_el = 6,000 K | `models/fe_ft_Tel6000_f32.model` (16 Sep 10:23:47) | `931313c45e03fdfebc4ee719985d4f20520eaf86a774f5fd25c34b4d8a208ad7` | `fe_ft_Tel6000_f32` | stage one | `phaseD_train_f32.sbatch`, `phaseD/train_Tel6000.xyz` on ARC; `production_campaign/arc_scripts/logs/fe_ft_Tel6000_f32_run-7.log` line 267 | `fe_phaseE.sbatch` → `data/coexistence/runs/` (T_el = 6,000 K campaign, 22-28 Sep), whose `analysis_live_all.json` sets each `runs_tel` cell's T_el; `build_cells*.sbatch` → `data/coexistence/cells/`; gates `gates_static_Tel6000.txt`, `gates_d4_Tel6000.txt`, `gates_v1_*.txt`; `phaseE_reweight.py`; `decomp_tel.py` (`REF_MODEL`); `spot_md_gen.py`; `probe_skin.py`; `liquid_diag.py`; `phaseE_analyze_live.py` default |
| v1, T_el = 8,000 K | `models/fe_ft_Tel8000_f32.model` (21 Sep 21:45:10) | `6b657e97c8dcdd70646a24fe729e4698fbe6812c0c5b50df3b50b0a4349a6c49` | `fe_ft_Tel8000_f32` | stage one | `phaseD_train_f32.sbatch`, `phaseD/train_Tel8000.xyz` on ARC; `production_campaign/arc_scripts/logs/train_Tel8000_f32_resume.out` line 147 | gates `gates_static_Tel8000.txt`, `gates_v1_Tel8000_static.txt`; `phaseE_reweight.py`; `e2_control.py`; `tel_diff_check.py` |
| v0, T_el = 6,000 K (float64) | not in the repository: `fe_ft_Tel6000.model` = `fe_ft_Tel6000_v0.model` (same file, 11,443,236 B, 2 Sep 02:42:59; the `_v0` copy made 13 Sep) | `764f3c0499f293dba8e7a7e7154d7cd9fc64c989500296315425666c010a2951` | `fe_ft_Tel6000` | stage one | `phaseD_train.sbatch` (float64), `production_campaign/arc_scripts/logs/train_Tel6000_interactive.out` line 245, `production_campaign/arc_scripts/logs/fe_ft_Tel6000_run-7.log` line 210 | `data/training/gates_v0_static.txt`, `gates_v0_d4.txt`, `gates_v0_d4b.txt` only |

The v1 SWA stage failed. `phaseD_train_f32swa.sbatch` (run name `fe_ft_Tel6000_f32swa`, 16 Sep) resumed
from `fe_ft_Tel6000_f32swa_run-7_epoch-149.pt` and stopped with `KeyError: 'swa_lr'`
(`production_campaign/arc_scripts/logs/train_Tel6000_f32swa.out` line 126; the "training done: phaseD/models/fe_ft_Tel6000.model" line
after it is a fixed `echo` in the script). It wrote no `.model` file. No SWA attempt was made at
4,000 or 8,000 K. Every run that used v1 loaded the stage-one files above, by name in its script
(`provenance/model_paths_in_scripts.txt`) and in its gate-file header (`model: phaseD/models/...`).

The production file is the SWA stage-two output of the e3 training:
- `fe_tel_e3_stagetwo.model` and `fe_tel_e3_run-7_stagetwo.model` agree in all 102 tensors.
- It agrees with stage one (`fe_tel_e3.model`, `fe_tel_e3_run-7.model`) in only 66 of 102.
- Its time, 01:07:00, is the time of the log's stage-two save.
- The compiled copy (`fe_tel_e3_stagetwo_compiled.model`) was not written. The log says
  "Can't redefine method: forward on class ... GenericJointEmbedding".

`fe_tel_e3_run-7_stagetwo.model` (`19f1e842…`, identical weights) is left out as a duplicate.

### Electronic-temperature input

The checkpoint stores the scale of the T_el input. The `elec_temp` embedding spec saved with the model is
`{'type': 'continuous', 'per': 'graph', 'in_dim': 1, 'emb_dim': 32, 'scale': 0.0001, 'zero_init': True}`
(`production_campaign/arc_scripts/logs/fe_tel_e3_run-7.log`, "Embedding specifications"; `melt_curve/validation/out/training.json`, key
`args.elec_temp_embedding`). In the tel-mace fork the T_el in kelvin is multiplied by that scale inside
the model (`mace/calculators/elec_temp.py`: "(Tel/10000) lives in that spec's `scale` and therefore inside
the checkpoint"; `provenance/tel_mace_elec_temp_lines.txt`). The ARC scripts therefore set
`atoms.info["elec_temp"]` in kelvin. Upstream mace-torch 0.3.16 loads the same file but does not apply
the stored scale, so `melt_curve/validation/load_model.py` passes `T_el / 1e4`. On 346 test frames this reproduces
the predictions stored on ARC: 0.0 meV/atom in energy and at most 4.6e-5 eV/Å per force component on a
Mac (`melt_curve/validation/out/repro_check.json`; the force figure varies by a few 1e-6 between runs with CPU threading).

### Which model the T_el-conditioned coexistence runs loaded

`data/coexistence/runs_tel/` gives T_m(330 GPa) = 6,329 K (`melt_curve/validation/out/melting_point.json`).

- **The script.** `fe_phaseE_tel.sbatch` sets `M=phaseD/models/fe_tel_e3_stagetwo.model`. The script file
  was last changed on 12 Oct, during the campaign (4-13 Oct). The model file is from 4 Oct 01:07, before
  the first test run at 4 Oct 19:57 (`production_campaign/arc_scripts/logs/coxT_test.out`).
- **The records.** No run record names the model: the run json files have no model field, the
  `production_campaign/arc_scripts/logs/coxT_*` logs do not print it, and `provenance/model_paths_in_run_records.txt` is empty.
- **No logged energy.** No log or run file holds a potential energy of a saved configuration, so no
  logged energy can be matched to a checkpoint.
- **The logged pressure.** Each checkpoint does hold one model-dependent number for its final
  configuration: the last entry of `P_t`. This is the virial pressure that `phaseE_coex_tel.py` computed
  with the loaded model for exactly the positions and cell it then saved.
- **The test.** [`tests/test_runs_tel_model.py`](../tests/test_runs_tel_model.py) evaluates that pressure
  for all 48 checkpoints with `fe_tel_e3_stagetwo.model` and with `fe_tel_e3.model`, at each run's T_el.
  T_el is checked against the script's rule: round(T_m) of the cell in `data/coexistence/analysis_live_all.json`.
  The test requires stage two to reproduce the logged pressure to 0.02 GPa and stage one not to.
- **Status.** The test needs a GPU: a 12,096-atom evaluation on the Mac's CPU ran out of memory. It had
  not been run when this file was written; its result goes to `tests/out/runs_tel_model.json`.

## Third-party foundation model (not in the repository)

| File | sha256 | Size | Source |
|---|---|---|---|
| `MACE-matpes-pbe-omat-ft.model` | `e618ad582b84239905b9c3b77ce6e9ce111b0ecd1533223a1a6aac7a696b8aa0` | 79,471,284 B | MACE-MATPES-PBE foundation model fine-tuned from OMAT (MACE foundation-model release; Batatia et al., J. Chem. Phys. 163, 184110 (2025); MatPES: Kaplan et al., arXiv:2503.04070). All Phase C sampling and every fine-tuning run started from it (`--foundation_model ~/fe_melt/MACE-matpes-pbe-omat-ft.model`). Check a download against the sha256 above. |

## Training and test data of the production model

- **Training:** `data/training/train_e3.xyz`, 6,609 frames. This is the concatenation of the T_el-stamped
  4,000, 6,000 and 8,000 K training sets (2,173 + 2,218 + 2,218 frames), checked byte for byte. MACE held
  out a random 5% for validation (330 frames, indices in `fe_tel_e3_valid_indices_7.txt`).
- **Held-out test:** `data/training/test_Tel{4000,6000,8000}.xyz` (114 + 116 + 116 frames). These are the files
  `parity_tel.py` read; it sets `elec_temp` per frame itself.
- **Frame size:** 128 atoms, except fcc frames, which have 108.
- **Numbers:** `melt_curve/validation/out/training.json` (`melt_curve/validation/train_table.py`).
