# melt_curve/figures/ — figure producers (iron chat, October 2026; wired to the repository 26 Oct 2026)

These made the figures of the 19 Oct 2026 group-meeting deck. The deck itself is not part of the
repository. Every path now comes from `../layout.py` through `figstyle.py`: `FIG` is
`melt_curve/figures/inputs/` (intermediate json), `DATA` is `data/coexistence/`, and `save(fig, name)` writes
`docs/figures/<name>.pdf` and `.png` (200 dpi). `build_figures.py` rebuilds `docs/figures/` (see its
README there); `figs_validation.py` draws the validation outputs.

| Script | Writes | Reads |
|---|---|---|
| `parse_train.py` | `inputs/train_meta.json` | `data/training/train_e3.xyz` (split by `elec_temp`), `data/training/test_Tel6000.xyz`, `data/training/spot_frames.xyz` |
| `figs_train.py` | training figures, `inputs/tel_check.json` | `inputs/train_meta.json`, `data/coexistence/pn_correct.json`, `production_campaign/arc_scripts/logs/fe_tel_e3_run-7.log`, `production_campaign/arc_scripts/logs/e3_test_final.out` |
| `parity_tel.py` (ARC) then `figs_parity.py [parity_tel.npz]` | parity figure, `inputs/parity.json` | `data/training/parity_tel.npz` (default) |
| `figs_coex.py` | coexistence-cell figures, `inputs/f_T_extract.json`, `inputs/f_cell_q6.json` | `data/coexistence/runs_tel/*_ckpt.npz`, `inputs/q6_330.npy` (per-atom q-bar-6 of one 330 GPa frame; a cache, not committed; recomputed if absent) |
| `figs_results.py` | melting-curve figures, `inputs/results.json` | `melt_curve/validation/refs/literature.json`, `data/coexistence/pn_correct.json`, `data/coexistence/analysis_live_all_runs_tel.json`, `data/coexistence/decomp_tel_tel.json`, run checkpoints |
| `figs_phase.py`, `figs_triple.py`, `figs_lit.py` | phase-diagram and literature figures | `data/coexistence/pn_correct.json`, `inputs/results.json`, `refs.py` |
| `figs_dft.py`, `figs_stack.py`, `figs_illus.py` | method illustrations | none, or `q6_330.npy` |
| `compare_curves_v1.py` | Tel-6000 vs T_el-conditioned melting curves | `data/coexistence/analysis_live_all.json`, `data/coexistence/analysis_live_all_runs_tel.json` |
| `figs_validation.py` | `f_heldout_errors`, `f_eos_error`, `f_cc_check`, `f_pbe_vs_exp` | `melt_curve/validation/out/*.json` |
| `build_figures.py` | `docs/figures/` | runs the producers above in order; checks the ledger is unchanged |
| `figstyle.py` | shared style; `FIG`, `DATA`, `OUT` and `save()` from `layout.py` | |
| `refs.py` | literature reference list | |

`inputs/*.json` were the intermediate files as they stood on 19-20 Oct 2026. `build_figures.py`
rebuilds `train_meta.json`, `parity.json`, `results.json` and `tel_check.json`; on 26 Oct 2026 the
rebuilt files were identical to the chat's. `f_T_extract.json` and `f_cell_q6.json` come from
`figs_coex.py`, which is not in the build; the staleness guard (`melt_curve/validation/ledger_inputs.tsv`) tracks
them.
