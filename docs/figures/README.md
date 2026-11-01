# docs/figures — figures rebuilt from the committed producers

Each figure is here as a vector PDF and a 200-dpi PNG. To rebuild them, run
`python analysis/figures/build_figures.py` from the repository root (about 2 min; needs the LFS files).
The build also rewrites the intermediate json in `analysis/figures/inputs/` and refuses to finish if
`validation/ledger.json` or any `validation/out/*.json` changed. Key names are those of
`validation/ledger.json`: `keys.<…>` for producer outputs, `summary.<…>` for derived numbers.

| Figure | Shows | Producer | Ledger keys |
|---|---|---|---|
| `f_compare` | melting curve of this work with published results to 375 GPa, and the inner-core range | `analysis/figures/figs_results.py` | `melting_point.fit.*`, `melting_point.Tm_330`, `literature.melting_curve_figure.*` |
| `f_curve` | the coexistence cells (P_n, T_m, s.e.), the Simon fit and the fcc-liquid point | `figs_results.py` | `melting_point.cells.*.{P_n,Tm,sem}`, `melting_point.fit.*`, `coex_cells.cells.coex_fcc_100GPa.{P,Tm}` |
| `f_pressure` | kinetic pressure per cell, and interface-normal minus cell-mean pressure | `figs_results.py` | `melting_point.cells.*.dP`; the kinetic term is only in `analysis/figures/inputs/results.json` |
| `f_parity` | MACE vs DFT force components on the held-out test set, and the error distribution per T_el | `figs_parity.py` (reads `phaseD/parity_tel.npz`) | `test_stats.{4000,6000,8000}_all.frmse`, `summary.test_n_frames` |
| `f_heldout_errors` | energy and force RMSE of the held-out test set by phase and T_el | `figs_validation.py` | `test_stats.{T}_{phase}.rmse`, `test_stats.{T}_{phase}.frmse` |
| `f_eos_error` | MACE − PBE pressure along the hcp EOS (with the ±0.5 GPa basis band) and the energy error | `figs_validation.py` | `static_compare.eos.*`, `summary.eos_max_abs_dP_GPa`, `summary.eos_max_abs_dE_meV` |
| `f_cc_check` | Clausius-Clapeyron check: slope of T_m(P) against ΔV_m/ΔS_m per pair of cells | `figs_validation.py` | `cc_check.intervals.*.{curve,curve_err,cc,cc_err,z}`, `summary.cc_max_abs_z_above_258GPa` |
| `f_pbe_vs_exp` | static PBE + vibrational pressure minus three 300 K experimental EOS | `figs_validation.py` | `pbe_vs_exp.rows.*.{P_exp,diff}`, `pbe_vs_exp.V_core`, `summary.pbe_minus_exp_core_GPa` |
| `f_train_curve` | validation energy and force RMSE per epoch of the e3 training | `figs_train.py` (reads `logs/fe_tel_e3_run-7.log`) | end points: `training.model.final_errors.stage_two.valid_*`; the per-epoch curve is not in the ledger |
| `f_tel_check` | same-frame T_el check: liquid − hcp of F(T_el) − F(6,000 K), DFT and model | `figs_train.py` (reads `logs/e3_test_final.out`) | `tel_check.cells.*` |
| `f_coverage` | (P, T) of the sampling runs behind the training frames, by phase | `figs_train.py` | `training.train_file.per_phase`; per-frame P and T only in `inputs/train_meta.json` |
| `f_forces` | DFT force-magnitude distribution by phase | `figs_train.py` | none (from `inputs/train_meta.json`) |

Other producers in `analysis/figures/` are not part of this build. `figs_coex.py`, `figs_phase.py`,
`figs_triple.py`, `figs_lit.py`, `figs_dft.py`, `figs_stack.py`, `figs_illus.py` and
`compare_curves_v1.py` made method illustrations and talk figures. They use the layout module and write
here if run. `figs_coex.py` also writes `inputs/f_T_extract.json` and `inputs/f_cell_q6.json`; those two
committed files are not rebuilt, and the staleness guard tracks them (`validation/ledger_inputs.tsv`).
