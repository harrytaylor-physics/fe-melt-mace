# figures_src/ — figure producers written in the iron chat (October 2026)

These made the figures of the 19 Oct 2026 group-meeting deck. The deck itself is not part of the
repository. The scripts still use relative paths from the chat's working directory; wire them to
the repository layout module before use.

| Script | Writes | Reads |
|---|---|---|
| `parse_train.py` | `inputs/train_meta.json` | `phaseD/{train,test}_Tel*.xyz`, `phaseD/spot_frames.xyz` (written as `arc/phaseD/...`) |
| `figs_train.py` | training figures, `inputs/tel_check.json` | `inputs/train_meta.json`, `phaseE/pn_correct.json`, `logs/fe_tel_e3_run-7.log`, `logs/e3_test_final.out` (written as `arc/logs/...`) |
| `parity_tel.py` (ARC) then `figs_parity.py <parity_tel.npz>` | parity figure, `inputs/parity.json` | `parity_tel.npz` |
| `figs_coex.py` | coexistence-cell figures, `inputs/f_T_extract.json`, `inputs/f_cell_q6.json` | `phaseE/runs_tel/*_ckpt.npz`, `q6_330.npy` (per-atom q-bar-6 of one 330 GPa frame; recomputed if absent) |
| `figs_results.py` | melting-curve figures, `inputs/results.json` | `phaseE/pn_correct.json`, `phaseE/analysis_live_all_runs_tel.json`, `phaseE/decomp_tel_tel.json`, run checkpoints |
| `figs_phase.py`, `figs_triple.py`, `figs_lit.py` | phase-diagram and literature figures | `phaseE/pn_correct.json`, `inputs/results.json`, `refs.py` |
| `figs_dft.py`, `figs_stack.py`, `figs_illus.py` | method illustrations | none, or `q6_330.npy` |
| `compare_curves_v1.py` | Tel-6000 vs T_el-conditioned melting curves | `phaseE/analysis_live_all.json`, `phaseE/analysis_live_all_runs_tel.json` |
| `figstyle.py` | shared style; `FIG` (default `inputs/`, env `FE_FIG_DIR`) and `DATA` (default `../phaseE/`, env `FE_PHASEE`) | |
| `refs.py` | literature reference list | |

`inputs/*.json` are the intermediate files as they stood on 19-20 Oct 2026. Each has a producer
above; rebuild them rather than trusting the committed copies, and put them under the staleness guard.
