# MIGRATION — ~/fe_melt on ARC into fe-melt-mace

On 26 Oct 2026 the iron melting-curve project was brought into this repository on the branch
`import/fe_melt`. Its working directory on ARC, `~/fe_melt`, **was not under version control**: it
contains no `.git`, and no git repository sits inside it (`nested_git_repos.txt` is empty). The
repository therefore starts from a **snapshot without history**, staged read-only on ARC at
2026-10-25T17:46:37Z by [`tools/stage_fe_melt.sh`](tools/stage_fe_melt.sh) (`SUMMARY.txt` of the
stage). `~/fe_melt` itself was not moved or changed. Jobs were still running in `~/fe_melt/labels`
(80 LDA-shift snapshot calculations, 28 EOS-anchor calculations); their results come in a later change.

**Dates.** The only record of when each file last changed is the file time stored in the snapshot
archive `worktree.tgz`, which extraction preserves. `worktree_manifest.tsv` lists every file with size,
sha256 and included/excluded status, but no time. File dates quoted in this repository (CHANGELOG,
`validation/ledger_inputs.tsv`, `models/README.md`) are those archive times, in ARC local time.

## Sources

| Source | What it is |
|---|---|
| `~/fe_melt_stage/worktree.tgz` | every file of `~/fe_melt` except QE raw outputs under `labels/` and bytecode (the stage script's exclusions; no file exceeded its 200 MB limit and none matched its credential scan) |
| `~/fe_melt_stage/worktree_manifest.tsv` | path, size, sha256 and status of all 24,849 files |
| `~/fe_melt_stage/provenance/` | model files and per-tensor comparison, training-log lines, model paths in scripts and run records, tel-mace commit and environment, QE version, pseudopotential hashes |
| `~/fe_melt_stage/chat_handoff/` | files written in the iron chat outside `~/fe_melt` (23 Oct 2026): `validation/`, `pseudo_lda/`, `figures_src/`, `arc_scripts/`, `stage_fe_melt.sh` |

`~/fe_melt_backup/` (an older partial copy) was not needed: every file this repository keeps was in the snapshot.

## Layout after the migration

The paths of `~/fe_melt` are kept at the repository root, so the ARC scripts' relative paths
(`phaseD/models/...`, `phaseE/runs_tel/...`) and the provenance records stay valid.

    *.py, *.sbatch           the ARC pipeline (Phase C sampling and labels, D dataset/training/gates,
                             E coexistence, F LDA shift and EOS anchor), verbatim
    labels/                  QE inputs (templates, running jobs), queue lists, inputs_as_run.tgz (LFS)
    pseudo/                  PAW files; pseudo/lda_generation/ the LDA PAW generation (new)
    phaseD/                  training/test sets, gates, parity predictions; models/ (LFS)
    phaseE/                  coexistence cells, runs, final configurations (LFS), analysis
    phaseX/, traj/, results/ LDA-shift frames, sampling metadata, the e3 training metrics
    logs/                    every log, verbatim
    exploratory/dmft_pilot/  the DFT+DMFT pilot (moved from dmft/pilot), not used in any result
    validation/              potential-vs-PBE and PBE-vs-experiment checks, producers, ledger.json (new)
    analysis/figures/        figure producers and their intermediate json (new)
    docs/figures/            the figures, PDF + 200-dpi PNG (new)
    models/README.md         model provenance (new)
    data/                    SHA-256 manifest of everything left out (new)
    tests/                   repository tests (new)
    tools/stage_fe_melt.sh   the snapshot script (new)
    layout.py, conftest.py   the one layout module; pytest collection (new)

Every location that the analysis, validation and figure code reads or writes is named once, in
[`layout.py`](layout.py). The ARC scripts at the root do not import it.

## Mapping

Source paths are relative to `~/fe_melt` unless marked *handoff* (`~/fe_melt_stage/chat_handoff/`).

| Source | Destination | Reason |
|---|---|---|
| `*.py`, `*.sbatch`, `RUN_LDASHIFT.md` (54 files) | same paths | the pipeline as it ran; verbatim (sha256 checked against the manifest) |
| `*_valid_indices_7.txt`, `coverage_map.txt`, `env_candidates.txt` | same | the validation split of each training run; sampling configuration |
| `labels/{bcc,fcc,hcp,liquid,stackfault}-0000_Tel{4000,6000,8000}.in` (15) | same | one QE input template per phase and T_el. `labels/*.in` total 52.6 MB, over the 50 MB limit for keeping all of them as plain files |
| `labels/{spot_,ldaX_,eosA_,mixtest}*.in` (121) | same | inputs of the spot-check, LDA-shift, EOS-anchor and mixing-test jobs |
| all 7,445 `labels/*.in` | `labels/inputs_as_run.tgz` (LFS, 6.4 MB) | every input unmodified, file times kept; members verified against the manifest sha256 |
| `labels/*.sbatch`, `labels/*.txt`, `labels/collect.py`, `labels/{q,mq}_*`, `labels/old_queue/` | same | label queues and job scripts |
| `pseudo/*.UPF` | same | the PBE (pslibrary 1.0.0) and LDA PAW files the labels used |
| `Fe.pbe-spn-kjpaw_psl.1.0.0.regen.UPF` | `pseudo/lda_generation/` | the PBE file regenerated with ld1.x, beside the LDA generation that checked it |
| `phaseD/gates_*.txt`, `phaseD/dataset_report.txt` | same | gate results; dataset sizes |
| `phaseD/train_e3.xyz` (LFS) | same | the production model's training file (`train_file` in `logs/fe_tel_e3_run-7_debug.log`) |
| `phaseD/test_Tel{4000,6000,8000}.xyz` (LFS) | same | the held-out test sets; `parity_tel.py` reads these (unstamped) files |
| `phaseD/parity_tel.npz`, `phaseD/e2_control.npz` (LFS), `phaseD/spot_frames.xyz` | same | model predictions on the test frames; spot-check frames |
| `phaseD/models/fe_tel_e3_stagetwo.model`, `fe_tel_e3.model`, `fe_ft_Tel{4000,6000,8000}_f32.model` (LFS) | same | the production model, e3 stage one (identity test) and v1; see [`models/README.md`](models/README.md) |
| `results/fe_tel_e3_run-7_train.txt` (LFS, 20.6 MB) | same | per-epoch metrics of the production training run |
| `phaseE/` (228 files; `*_ckpt.npz` in LFS) | same | coexistence cells, both campaigns' runs and final configurations, reweighting, decomposition, analysis, Simon fit |
| `phaseX/ldashift/` | same | frames of the running LDA-shift jobs |
| `traj/*.json`, `traj/*.flagged` | same | metadata of the sampling runs (the trajectories are excluded) |
| `logs/` (1,110 files), `parity_tel.out` | same | every log, verbatim |
| `dmft/pilot/` (9 of 11 files) | `exploratory/dmft_pilot/` | exploratory; feeds no result. Its scripts' relative paths are unchanged and now point one level deeper |
| *handoff* `validation/` | `validation/` | the validation suite; see changes below |
| *handoff* `pseudo_lda/` (without its `Fe.pz-spn-kjpaw_psl.1.0.0.UPF`, identical to `pseudo/`) | `pseudo/lda_generation/` | the repository's pseudopotential directory |
| *handoff* `figures_src/` | `analysis/figures/` | with the analysis code |
| *handoff* `stage_fe_melt.sh` | `tools/stage_fe_melt.sh` | the script that made the snapshot |
| *handoff* `arc_scripts/` | — (compared only) | `ldashift.py`, `phaseF_common.py`, `ldashift_md.sbatch`, `RUN_LDASHIFT.md`, `eos_anchor.py` match the snapshot copies byte for byte (sha256); the snapshot copies are committed |
| *handoff* `NOT_FOR_REPO/` | — | report text; never committed. Only the numbers' ledger keys are (`validation/report_numbers.tsv`) |

## Changes made to imported files

- `labels/*.in` kept as plain files: `pseudo_dir = '/home/shug7909/fe_melt/pseudo'` became
  `pseudo_dir = './pseudo'`. Nothing else changed; the originals are in `labels/inputs_as_run.tgz`.
- `validation/paths.py`: a thin wrapper over `layout.py`, with the same names and environment
  variables. `PARITY` now points at `phaseD/parity_tel.npz`, where `parity_tel.py` writes it; the
  handoff default was the root.
- `validation/val_mace.py`: also writes the T_el response, the fcc/bcc pressures and the c/a minima
  (sections 2-4), which it printed but did not save. These keys of `validation/out/static_compare.json`
  come from a run on 26 Oct 2026 (Mac, CPU, float32, mace-torch 0.3.16). The other keys are the handoff
  values from the 20-23 Oct run, which the 26 Oct run reproduces to 9e-5 relative.
- `validation/load_model.py`: `calc()` takes a model path. The note on the T_el scale now cites the
  checkpoint's stored scale and the fork's code instead of inferring it.
- `validation/ledger.py`: reads six new producers' outputs and `refs/literature.json`, and adds derived
  summary numbers. All 418 earlier keys and 21 summary values are unchanged.
- New producers: `validation/train_table.py`, `tel_check.py`, `coex_cells.py`, `melting_point.py`,
  `repro_check.py` and `settings.py`. New supporting files: `validation/refs/literature.json`,
  `validation/record_inputs.py`, `validation/ledger_inputs.tsv` and `validation/report_numbers.tsv`.
- `analysis/figures/`:
  - `figstyle.py` takes its paths from `layout.py` and saves PDF plus 200-dpi PNG to `docs/figures/`.
  - `figs_train.py` reads `logs/` through the layout module instead of the chat's `arc/logs`.
  - `parse_train.py` reads `phaseD/train_e3.xyz` split by `elec_temp`. It previously read the per-T_el
    files, which are left out; the rebuilt `inputs/train_meta.json` is byte-identical.
  - `figs_results.py` draws the literature from `validation/refs/literature.json`. The rebuilt
    `inputs/results.json` is identical.
  - `figs_coex.py` and `figs_illus.py` cache `q6_330.npy` in `inputs/` (not committed).
  - `compare_curves_v1.py` reads and writes through the layout module.
  - New: `figs_validation.py` and `build_figures.py`.
- `tools/stage_fe_melt.sh`: `DATA` must now be set in the environment (ARC sets it) instead of
  defaulting to an absolute path.

## Snapshot scan and exclusions

The extracted snapshot (9,225 files) was scanned for each kind of file this repository must not hold:

| Kind of file | Found | Handling |
|---|---|---|
| Report or manuscript files: `*.tex`, `*.bib`, `*.pptx`, `*.docx`, `*.pdf`, drafts, proofs, cover letters | none | — |
| Meeting-deck material | `slides_manifest.txt`, `slides_logtypes.txt` (deck inventories) | excluded |
| QE raw outputs | the 15,617 `labels/*.out` (already left out of the archive by the stage script) and `dmft/pilot/scf.out` with its `CRASH` report | excluded |
| QE scratch (`*.save/`, `tmp*/`) and wavefunctions | none | — |
| Trajectories | `traj/*.extxyz` (98 sampling trajectories) | excluded |
| Files over 50 MB | `MACE-matpes-pbe-omat-ft.model` (79 MB), `phaseD/train_e3.xyz` (132 MB), `phaseD/train_e2.xyz` (88 MB), `labels/labels.jsonl` (55 MB) | all excluded except `train_e3.xyz`, the production training set, which is in LFS |
| Credential patterns in file contents (tokens, keys, passwords) | none; the stage's own `credential_scan.txt` is empty too | — |

`data/excluded_sha256.tsv` lists every excluded file with its sha256 and reason (15,788 files);
[`data/README.md`](data/README.md) summarises it. Further exclusions:

| Excluded | Reason |
|---|---|
| `phaseD/train_Tel{4000,6000,8000}.xyz`, `train_Tel*_stamped.xyz`, `train_e2.xyz`, `test_Tel6000_stamped.xyz` | superseded. `train_e3.xyz` is the concatenation of the three stamped files; each stamped file is its unstamped one plus `elec_temp` (checked byte for byte) |
| `selected/*.extxyz`, `labels/labels.jsonl` | the frame pool and collected labels; every labelled frame (6,955) is in `train_e3.xyz` + `test_Tel*.xyz` |
| 14 superseded models (v0, e2, e2b, the `_run-7` duplicates, `fe_tel_e3_run-7_stagetwo.model`), 4 `*_compiled.model`, 11 `*.pt` checkpoints | not loaded by any production, coexistence or validation run, or duplicates of a kept file; sha256 in the manifest; provenance in `models/README.md` |
| `MACE-matpes-pbe-omat-ft.model` | third-party foundation model; source and sha256 in `models/README.md` |
| `results/` except `fe_tel_e3_run-7_train.txt` | metrics and plots of the v0, v1, e2 and e2b runs, and the e3 run's plots |
| `dmft_pilot.tgz`, `ldashift.tgz` | archives of files kept unpacked (contents identical, checked with `cmp`) |
| `__pycache__/` | bytecode |

## Large files (git LFS)

`.gitattributes` tracks `*.model`, `*.npz`, `phaseD/train_e3.xyz`, `phaseD/test_Tel*.xyz`,
`labels/inputs_as_run.tgz` and `results/fe_tel_e3_run-7_train.txt`. That is 110 objects, 255 MB
(`git lfs ls-files -s`):

| Files | Size |
|---|---|
| `phaseD/train_e3.xyz` | 132 MB |
| the 5 models | 29 MB |
| 97 coexistence checkpoints (`runs_tel/`, `runs/`, one in `runs_teltest/`) | 54 MB |
| `results/fe_tel_e3_run-7_train.txt` | 21 MB |
| `labels/inputs_as_run.tgz` | 6 MB |
| 3 test sets, 2 prediction files | 10 MB |

`water-scw-mace` holds 135 MB in LFS (57 objects). Together they come to about 390 MB.

## Separate repositories

The T_el-conditioned MACE fork (tel-mace) is not copied in. It is pinned by commit in the README
("Dependencies"): `git@github.com:harrytaylor-physics/tel-mace.git`, branch `tel-conditioning`, commit
`5168bc0d08bfb0d85ac4504710ee00792d2da04c`.
