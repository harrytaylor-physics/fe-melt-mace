# MIGRATION — ~/fe_melt on ARC into fe-melt-mace

How ~/fe_melt on ARC was imported here.

## 5 Nov 2026: layout matched to water-scw-mace

The repository first kept the paths of `~/fe_melt` at its root (73 files), so that the ARC scripts could run
from a clone. On 5 Nov 2026 it moved to the layout of `water-scw-mace`, on the branch `water-layout`. The ARC
scripts are now provenance, as `production_campaign/` is there. `~/fe_melt` on ARC was not changed.

    production_campaign/     the ARC pipeline, verbatim (fe_melt/, arc_scripts/)
    models/                  the potentials (LFS) and their provenance
    data/                    training and test sets, coexistence output, the SHA-256 manifest
    melt_curve/              layout.py, validation/, figures/ (the role of subneptune/ in water-scw-mace)
    studies/                 lda_shift/, dmft_pilot/
    docs/                    METHODS.md, figures/
    tests/                   repository tests

| Source | Destination | Reason |
|---|---|---|
| root `*.py`, `*.txt`, `parity_tel.out` | `production_campaign/fe_melt/` | the ARC working directory, verbatim |
| `labels/`, `pseudo/`, `traj/`, `results/` | `production_campaign/fe_melt/` (same names) | keeps their relative paths to the scripts |
| root `*.sbatch`, `tools/stage_fe_melt.sh` | `production_campaign/arc_scripts/` | where water-scw-mace keeps Slurm scripts |
| `logs/` | `production_campaign/arc_scripts/logs/` | the same |
| `phaseD/models/*.model` | `models/` | models sit with their provenance |
| `phaseD/` (other files) | `data/training/` | derived data products |
| `phaseE/` | `data/coexistence/` | derived data products |
| `layout.py`, `validation/`, `analysis/figures/` | `melt_curve/layout.py`, `melt_curve/validation/`, `melt_curve/figures/` | one self-contained component that runs on any machine |
| `phaseX/ldashift/`, `RUN_LDASHIFT.md` | `studies/lda_shift/` (`RUN_LDASHIFT.md` becomes its `README.md`) | satellite study |
| `exploratory/dmft_pilot/` | `studies/dmft_pilot/` | satellite study |
| `exploratory/README.md` | removed; its text is in `studies/README.md` | |

What changed besides the paths:
- **Unchanged content.** Everything in `production_campaign/` except its new README, the models, datasets, coexistence output, `ledger.json` and
  `validation/out/*.json` are byte-identical to `main` before the move. No number changed.
- **`ledger_inputs.tsv`.** Paths rewritten; every sha256 and date is unchanged.
- **Edited.** `layout.py` (the new paths, plus `ARC` and `LABELS`), `conftest.py`, `.gitattributes`,
  `.gitignore`, the import paths in `tests/` and `melt_curve/`, and the path references in the READMEs and
  docstrings. Text quoting an ARC command, log or gate file keeps the ARC path.
- **Moved text.** The README's method description is now `docs/METHODS.md`, its stage table is in
  `production_campaign/README.md`, and its test section is in `TESTING.md`.
- **History.** No commit was rewritten. `git log --follow` traces a file across the move.

Water-scw-mace has `benchmarks/` and `convergence/`. Here their contents are in
`melt_curve/validation/` (`refs/dft/`, `refs/convergence/`), as `subneptune/validation/` holds that
component's gates.

The rest of this file describes the import of 26 Oct 2026 and uses the paths of that time.

On 26 Oct 2026 the iron melting-curve project was brought into this repository on the branch
`import/fe_melt`. Its working directory on ARC, `~/fe_melt`, **was not under version control**: it
contains no `.git`, and no git repository sits inside it (`nested_git_repos.txt` is empty). The
repository therefore starts from a **snapshot without history**, staged read-only on ARC at
2026-10-25T17:46:37Z by [`tools/stage_fe_melt.sh`](production_campaign/arc_scripts/stage_fe_melt.sh) (`SUMMARY.txt` of the
stage). `~/fe_melt` itself was not moved or changed. Jobs were still running in `~/fe_melt/labels`
(80 LDA-shift snapshot calculations, 28 EOS-anchor calculations). The snapshot calculations finished on
29 Oct and came in with the follow-up below. The EOS-anchor calculations were cancelled: `eos_anchor.py`
was not run to completion and nothing uses it.

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
[`layout.py`](melt_curve/layout.py). The ARC scripts at the root do not import it.

## Mapping

Source paths are relative to `~/fe_melt` unless marked *handoff* (`~/fe_melt_stage/chat_handoff/`).

| Source | Destination | Reason |
|---|---|---|
| `*.py`, `*.sbatch`, `RUN_LDASHIFT.md` (54 files) | same paths | the pipeline as it ran; verbatim (sha256 checked against the manifest) |
| `*_valid_indices_7.txt`, `coverage_map.txt`, `env_candidates.txt` | same | the validation split of each training run; sampling configuration |
| `labels/{bcc,fcc,hcp,liquid,stackfault}-0000_Tel{4000,6000,8000}.in` (15) | same | one QE input template per phase and T_el. `labels/*.in` total 52.6 MB, over the 50 MB limit for keeping all of them as plain files |
| `labels/{spot_,ldaX_,eosA_,mixtest}*.in` (121) | same | inputs of the spot-check, LDA-shift, EOS-anchor and mixing-test jobs (EOS anchor: not run to completion, unused) |
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

## Follow-up of 29 Oct 2026: the 330 GPa snapshot test

Added on 29 Oct 2026 to the same branch, from two sources:

| Source | What it is |
|---|---|
| `~/fe_melt_stage/followup/arc/` | files copied from `~/fe_melt` on ARC on 29 Oct 2026, at their `~/fe_melt` paths; `fe_melt_followup_arc.tgz` beside it has the same content (checked with `diff -r`) |
| `~/fe_melt_stage/followup/fe_followup/` | files from the iron chat (29 Oct): `validation/` producers, references, outputs and test; `tools/ldaXc.sbatch` |

| Source | Destination | Reason |
|---|---|---|
| *arc* `phaseX/ldashift/ldashift_result.json` | same | `python ldashift.py eval` on ARC, 29 Oct |
| *arc* `ldashift_eval.txt` | `phaseX/ldashift/ldashift_eval.txt` | stdout of that eval |
| *arc* `phaseX/ldashift/frames.xyz`, `md_meta.json` | same (already committed) | identical to the 25 Oct snapshot copies |
| *arc* `ldashift.py` | same | ARC copy of 27 Oct; differs from the snapshot copy only where T_m(330 GPa)'s error was 26 K (now 6 K) |
| *arc* `labels/ldaX_*_LDA_Tel6329.in` (40) | same, `pseudo_dir` rewritten to `./pseudo` | the 27 Oct LDA inputs (PBE PAW file, `input_dft = 'PZ'`). The originals are in `labels/ldaX_LDA_inputs_as_run.tgz`. `labels/inputs_as_run.tgz` holds the superseded 22 Oct versions, which used the generated PZ file |
| *arc* `labels/ldaX_*_PBE_Tel6329.in` (40), `ldaX_list.txt`, `ldaX_array.sbatch` | same (already committed) | unchanged since the snapshot (apart from `pseudo_dir`) |
| *arc* `labels/ldaX_list2.txt`, `ldaX_array2.sbatch`, `ldaXc.sbatch` | same | job scripts of 27-29 Oct. *chat* `tools/ldaXc.sbatch` is byte-identical to `labels/ldaXc.sbatch` and is not committed twice |
| *arc* `ldaX_out_sha256.txt` | rows of `data/excluded_sha256.tsv` | sha256 of the 80 QE outputs, which are not committed (no sizes were recorded) |
| *chat* `validation/conv_check.py`, `lda_setup_check.py`, `ldashift_summary.py` | `validation/` | new producers |
| *chat* `validation/refs/convergence/`, `refs/lda_setup/` | same | QE 6.6 runs (json and inputs) |
| *chat* `validation/out/convergence.json`, `lda_setup.json` | same | rebuilt here byte-identical; `out/ldashift.json` was written here by `ldashift_summary.py` |
| *chat* `validation/tests/test_followup.py` | same | the three producers against their outputs |
| *chat* `NOT_FOR_REPO/awe_validation_report_text_rev46.md` | — | report text; never committed. `validation/report_numbers.tsv` was rebuilt against it, then updated for the two passages rev 47 changed (from Harry's description, 29 Oct) |

Changes made to these and to earlier files:
- `ldashift.py`: `lda_input()` writes the LDA setup that ran (PBE PAW file, `input_dft = 'PZ'` as the
  first line of `&SYSTEM`); the module docstring explains the mixed setup; the check for the generated
  PZ file is gone. `RUN_LDASHIFT.md` says the same and quotes the Simon-fit error.
- `validation/tests/test_followup.py`: puts `validation/` on `sys.path`, so it runs from the repository root.
- `layout.py`: `LDASHIFT` (`phaseX/ldashift`, `FE_LDASHIFT`), also in `validation/paths.py`.
- `validation/settings.py`: direct-test, mixed-setup and coexistence-run settings (new keys only).
- `validation/ledger.py`: reads the three new outputs; derived keys for the PBE T_m(330 GPa), the
  Rao & Zhu gap and the Simon-fit error. All earlier keys and summary values are unchanged.
- `validation/record_inputs.py`, `ledger_inputs.tsv`: the new inputs are under the staleness guard.
- `analysis/figures/figs_lit.py`: draws the Simon-fit error instead of 26 K.
- `pseudo/lda_generation/README.md`: the generated PZ dataset diverges in displaced configurations and
  is not used by any result.

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
