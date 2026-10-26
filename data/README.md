# data/ — what was left out, and how to check a copy

The iron campaign's raw data stay on ARC in `~/fe_melt` (and its backups); this repository holds the
inputs, the processed datasets, the models and the results. `excluded_sha256.tsv` lists every file of
the 25 Oct 2026 snapshot of `~/fe_melt` that is **not** in the repository, with its size, sha256 and the
reason. Use it to check that a file found elsewhere is the one the campaign used:
`shasum -a 256 <file>` and look the hash up in the third column.

Columns: `path` (relative to `~/fe_melt`), `bytes`, `sha256`, `reason`.

| Reason | Files | Size |
|---|---|---|
| QE raw output (`labels/*.out`: the PBE labels, the spot checks, the LDA-shift and EOS-anchor jobs) | 15,617 | 978 MB |
| sampling trajectory (`traj/*.extxyz`) | 98 | 62 MB |
| frame pool (`selected/`; the labelled frames are in `phaseD/train_e3.xyz`, `phaseD/test_Tel*.xyz`) | 6 | 52 MB |
| `labels/labels.jsonl` (over 50 MB; its 6,955 labels are the train and test frames above) | 1 | 55 MB |
| superseded datasets (`phaseD/train_Tel*.xyz`, `train_Tel*_stamped.xyz`, `train_e2.xyz`, `test_Tel6000_stamped.xyz`) | 8 | 354 MB |
| superseded models, compiled copies and training checkpoints (`phaseD/models/`) | 29 | 290 MB |
| training metrics and plots other than the e3 run's `_train.txt` (`results/`) | 16 | 55 MB |
| third-party foundation model `MACE-matpes-pbe-omat-ft.model` (source in `models/README.md`) | 1 | 79 MB |
| QE output and crash report of the exploratory DMFT pilot | 2 | 0.1 MB |
| archives of files kept unpacked (`dmft_pilot.tgz`, `ldashift.tgz`), deck inventories, bytecode | 10 | 0.2 MB |

The superseded datasets were checked byte for byte before they were left out:
- `train_e3.xyz` is the concatenation of `train_Tel{4000,6000,8000}_stamped.xyz`.
- `train_e2.xyz` is the concatenation of the 4,000 K and 8,000 K stamped files.
- Each `*_stamped.xyz` is its unstamped file with ` elec_temp=<T>` added to every header line.

`labels/inputs_as_run.tgz` holds all 7,445 QE inputs of `labels/` unmodified. The plain `labels/*.in` in
the repository are 15 templates and the 121 inputs of the spot, LDA-shift and EOS-anchor jobs, with
`pseudo_dir` rewritten to `./pseudo`.

The QE outputs of the 80 LDA-shift (`ldaX_*`) and 28 EOS-anchor (`eosA_*`) jobs were still running when
the snapshot was taken, and the EOS-anchor runs have since failed and will be rerun. Their results come
in a later change.
