# CHANGELOG

Oldest first. `~/fe_melt` had no version control, so this history is rebuilt from three sources:
- the file times in the 25 Oct 2026 snapshot archive (see [MIGRATION.md](MIGRATION.md));
- the dated logs and result files;
- the amendment numbers (A*n*) that the script headers cite.

No amendment register exists in `~/fe_melt`; each A*n* below is quoted from the header that cites it.
Entries from 20 Oct onwards follow `HANDOFF.md` of the iron chat.

## 2026-08-04 to 2026-08-27 — Phase C: sampling and PBE labels

- **Sampling.** Started from the MACE-MATPES-PBE foundation model fine-tuned on OMAT
  (`MACE-matpes-pbe-omat-ft.model`, 4 Aug), which A4 chose as the base. Sampled bcc, fcc, hcp and liquid
  iron (`phaseC_sample.py`, `fe_sampleC.sbatch`, `logs/sampleC_*`, 4 Aug). The fcc volumes weight
  80-130 GPa (A2).
- **Tethered resampling and fill runs.** `phaseC_tether.py`, `phaseC_fill.py` (4-5 Aug).
- **Frame selection.** `selected/` (5 Aug); stacking-fault frames added 27 Aug.
- **QE inputs.** PBE Mermin free-energy labels at T_el = 4,000, 6,000 and 8,000 K: 128-atom cells
  (108 for fcc), pslibrary 1.0.0 PAW with 16 valence electrons, ecutwfc 100 Ry
  (`phaseC_qe_labels.py`, `labels/*.in`, 18 Aug).

## 2026-08-18 to 2026-09-13 — Phase D: datasets and the first potentials

- **Datasets.** `phaseD_dataset.py` (18 Aug) matched the labels to the frames and split 95/5 into
  `phaseD/train_Tel*.xyz` and `test_Tel*.xyz`, written 13 Sep. Per T_el: 2,173-2,218 training and
  114-116 test frames (`phaseD/dataset_report.txt`).
- **v0** (float64, T_el = 6,000 K). `fe_ft_Tel6000.model`, 2 Sep. Gates `phaseD/gates_v0_*.txt`
  (`phaseD_gates.py`, 2 Sep).

## 2026-09-14 to 2026-09-22 — v1: one float32 potential per T_el

- **Training.** `fe_ft_Tel6000_f32.model` (16 Sep), `fe_ft_Tel4000_f32.model` (17 Sep) and
  `fe_ft_Tel8000_f32.model` (21 Sep), all stage one (`phaseD_train_f32.sbatch`).
- **Gates.** `phaseD/gates_static_Tel*.txt`, `gates_v1_*.txt`, `gates_d4_Tel6000.txt` (16-22 Sep).
- **SWA failure.** An SWA continuation at 6,000 K (`phaseD_train_f32swa.sbatch`, 16 Sep) stopped with
  `KeyError: 'swa_lr'` and wrote no model.

## 2026-09-17 to 2026-10-05 — Phase E: coexistence with the T_el = 6,000 K potential

- **Cells.** Eight two-phase cells: hcp at 100-350 GPa and fcc at 100 GPa, 12,096 and 12,032 atoms
  (`phaseE_build.py`, `phaseE/cells/`, 20-25 Sep).
- **Runs.** 48 NVE runs of 200 ps, segmented and self-resubmitting (`phaseE_coex.py`, `fe_phaseE.sbatch`,
  `phaseE/runs/`, `logs/coex_*`, 22-28 Sep).
- **Analysis.** `phaseE/analysis_live_all.json` (5 Oct).
- **Reweighting.** Reweighting to T_el = 4,000 and 8,000 K with the three v1 potentials (A21, step 3):
  `phaseE_reweight.py`, `phaseE/reweight/` (2-3 Oct).

## 2026-09-28 to 2026-10-04 — the T_el-conditioned potential

- **Fork.** tel-mace `tel-conditioning`, commit 5168bc0 (22 Sep): T_el enters as a learned
  per-structure embedding, kelvin in, with the 1e-4 scale stored in the checkpoint and zero
  initialisation. It runs only in its own environment `~/tel_env` (A12 rule).
- **E2.** Trained on 4,000 and 8,000 K only, judged on the held-out 6,000 K set (A12, A14; criteria
  registered 29 Sep, A17). Runs `fe_tel_e2`, `fe_tel_e2b` (with SWA stage two) on 28 Sep-1 Oct;
  `logs/e2_result.out`, `logs/e2b_result.out`, `logs/e2b_stagetwo_result.out`.
- **e3, the production potential.** Trained on all three T_el (`phaseD/train_e3.xyz`, 6,609 frames;
  `phaseD_train_e3.sbatch`, 2-4 Oct). The SWA stage two output `fe_tel_e3_stagetwo.model` was saved
  4 Oct 01:07. Final errors: 1.5 / 1.5 meV/atom and 89.1 / 92.4 meV/Å (train / validation).
- **Checks (A21, A22).**
  - Spot check against new DFT frames: `spot_md_gen.py`, `spot_dft_gen.py`, `spot_eval.py`.
  - Same-frame T_el check on the 2,009 frames labelled at all three T_el: `tel_diff_check.py`,
    `logs/e3_test_final.out` (4 Oct).

## 2026-10-04 to 2026-10-14 — coexistence with the T_el-conditioned potential

- **Runs.** The same 48 runs with `fe_tel_e3_stagetwo.model`, each cell at the T_el of its 6,000 K
  melting temperature (A26). `fe_phaseE_tel.sbatch`, `phaseE_coex_tel.py`, `phaseE/runs_tel/`,
  `logs/coxT_*`, 4-13 Oct. On 12 Oct the segment wall time was tied to the job's remaining time.
- **Analysis.** Kinetic pressure term included (A27): `phaseE/analysis_live_all_runs_tel.json`, 13 Oct.
- **Decomposition (A29).** First-order T_el and model decomposition, and the 350 GPa energy-conservation
  check: `decomp_tel.py`, `decomp_tel.sbatch`, `phaseE/decomp_tel_*.json`, 13 Oct.
- **Melting point.** Interface-normal pressure P_n (A29c) and the Simon fit: `pn_correct.py`,
  `phaseE/pn_correct.json`, 13 Oct. T_m(330 GPa) = 6,329 K, dT_m/dP = 8.41 K/GPa.
- **Held-out test set.** Predictions of the production potential: `parity_tel.py`,
  `phaseD/parity_tel.npz`, 14 Oct.

## 2026-10-17 — DFT+DMFT pilot (exploratory)

- hcp Fe at 7.05 Å³/atom, 5,800 K, TRIQS + QE (`dmft/pilot/`, now `exploratory/dmft_pilot/`). The QE
  step stopped with an MPI pool error. It feeds no result.

## 2026-10-20 — validation suite (iron chat)

- `validation/` compares the potential with PBE and PBE with experiment:
  - 26 new static PBE references (QE 6.6, `val_dft.py`) and 8 cold-EOS references;
  - the held-out test statistics (`stats_test.py`), the static comparison (`val_mace.py`), the
    Clausius-Clapeyron check (`cc_check.py`) and PBE against three 300 K EOS (`pbe_vs_exp.py`);
  - every number written to `ledger.json` (`ledger.py`), with 6 tests.

## 2026-10-21 to 2026-10-22 — LDA PAW, the LDA − PBE test, EOS anchor (iron chat, then ARC)

- **LDA PAW.** An LDA PAW file in the pslibrary 1.0.0 construction (`pseudo/Fe.pz-spn-kjpaw_psl.1.0.0.UPF`,
  ld1.x inputs and logs in `pseudo/lda_generation/`), checked against the official PBE file
  (`upf_compare.py`).
- **LDA − PBE test.** `ldashift.py`, `phaseF_common.py`, `ldashift_md.sbatch`, `RUN_LDASHIFT.md`
  (21 Oct): 20 hcp and 20 liquid snapshots at 330 GPa and T_el = 6,329 K (`phaseX/ldashift/`), labelled
  with PBE and with LDA (`labels/ldaX_*`, 80 calculations).
- **EOS anchor.** The 28 runs (`eos_anchor.py`, `labels/eosA_*`) were submitted on 22 Oct. They were
  later cancelled and dropped (see 29 Oct); nothing uses them.

## 2026-10-23 — three report numbers corrected (iron chat)

- Bulk modulus "within 0.1%" became 0.2%: at 8.0 Å³/atom it is 774.2 vs 775.6 GPa.
- E2g "+0.9%" became +1.0%: 333.8 vs 337.0 cm⁻¹.
- The Clausius-Clapeyron shortfall below 258 GPa "30-50%" became 28-48%.

## 2026-10-26 — this repository

- **Import.** `~/fe_melt` imported from the 25 Oct snapshot without history (MIGRATION.md). The chat's
  validation, LDA generation and figure producers were added and wired to one layout module (`layout.py`).
- **New producers** for the numbers that had none:
  - the training table (`train_table.py`);
  - the same-frame T_el check (`tel_check.py`);
  - the coexistence cells (`coex_cells.py`);
  - T_m(330 GPa) and dT_m/dP (`melting_point.py`);
  - the reproduction of the stored predictions and the float64 offset (`repro_check.py`);
  - the method settings (`settings.py`);
  - the literature values (`refs/literature.json`).
- **Report numbers.** Every number in the validation report is listed with its ledger key
  (`validation/report_numbers.tsv`).
- **Guards and tests.**
  - A staleness guard over the inputs no producer here rebuilds (`validation/ledger_inputs.tsv`).
  - A GPU test of which model the T_el-conditioned runs loaded (`tests/test_runs_tel_model.py`).
- **Figures.** Rebuilt as PDF plus 200-dpi PNG in `docs/figures/`.

## 2026-10-26 — the LDA snapshot runs diverge (ARC)

- All 40 LDA runs of the 330 GPa snapshot test timed out with a diverging SCF: the first step's estimated
  SCF accuracy was about 5,500 Ry, with negative charge.
- A local cross test on a displaced 16-atom hcp cell (σ = 0.15 Å, 7.0 Å³/atom) located the cause in the
  PZ PAW dataset generated with ld1.x. First-step SCF accuracy, by PAW file and functional:

  | PAW file | Functional | Ry |
  |---|---|---|
  | PBE | PBE | 4.49 |
  | PBE | PZ | 4.47 |
  | PZ | PZ | 302.8, negative charge |
  | PZ | PBE | 8.63, negative charge |

  The perfect-crystal tests of 22 Oct had converged, which hid it (`pseudo/lda_generation/README.md`).

## 2026-10-27 — LDA runs on the PBE PAW file (ARC)

- The 40 LDA inputs were switched to `Fe.pbe-spn-kjpaw_psl.1.0.0.UPF` with `input_dft = 'PZ'`
  (`labels/ldaX_*_LDA_Tel6329.in`). `labels/ldaX_list2.txt` (59 of the 80 inputs) and
  `ldaX_array2.sbatch` date from 27 Oct.
- The PBE runs had taken at most 58 min, so the time limit went from 6 h to 2 h. The limits that applied
  were set on the sbatch command line and override the script headers (`--time=01:30:00` in both array
  scripts, none in `ldaXc.sbatch`): 2 h for the HTC arrays 9194880, 9196535 and 9205402, and 3 h for the
  ARC array 13569880 (29 Oct).
- `ldashift.py` on ARC quotes T_m(330 GPa) with the Simon-fit error, ±6 K, instead of ±26 K.

## 2026-10-29 — the 330 GPa snapshot test finishes (ARC)

- The HTC queue estimate for the remaining 61 runs was 6 Nov, so they went to the ARC cluster with
  `labels/ldaXc.sbatch`. Each task claims the next unclaimed input with an atomic `mkdir` in
  `labels/claims/`, so one array per cluster can run without duplicates.
- On ARC, `srun` started 48 single-rank copies of pw.x (`Number of MPI processes: 1`, then
  `invalid number of pools`). The script therefore uses `mpirun -np $SLURM_NTASKS` on ARC and `srun` on
  HTC. It stops a cluster's tasks after a run that fails in under 5 minutes.
- All 80 runs finished on 29 Oct. `python ldashift.py eval` wrote `phaseX/ldashift/ldashift_result.json`
  (stdout in `phaseX/ldashift/ldashift_eval.txt`). The 80 QE outputs stay on ARC; their sha256 are in
  `data/excluded_sha256.tsv`.
- The 28 EOS-anchor runs had been cancelled and dropped. Nothing uses them.

## 2026-10-29 — snapshot test, convergence and report rev 46-47 in this repository

- **Snapshot test.** The ARC inputs, job scripts and result are committed. `ldashift.py` is the ARC copy,
  with `lda_input()` changed to write the setup that ran; a test checks it against the 40 committed
  input pairs.
- **New producers** (iron chat):
  - `validation/conv_check.py`: basis and k-point convergence of the static PBE pressure. Going from
    100/1,000 to 130/1,300 Ry changes it by -0.53 GPa; going from k 12×12×8 to 16×16×10, by +0.24 GPa.
  - `validation/lda_setup_check.py`: what LDA with the PBE PAW file changes. The hcp pressure moves by
    +0.22 GPa and the fcc − hcp energy by -0.14 meV/atom (2 K).
  - `validation/ldashift_summary.py`: report numbers of the snapshot test.
- **Ledger.**
  - The PBE melting temperature: T_m(PBE, 330 GPa) = 6,329 − 23 = 6,306 ± 8 K (`summary.Tm_PBE_330_K`).
  - The Rao & Zhu gap: 6,225 − 6,306 = −81 ± 43 K.
  - The statistical error of T_m(330 GPa) is now the Simon fit's standard error, 6.25 K
    (`summary.Tm_330_stat_err_K`). The earlier ±26 K had no producer, and every use of it was replaced.
- **Report numbers** (`validation/report_numbers.tsv`) rebuilt against the report text at rev 46:
  - Removed: every ΔS_m value, the Clausius-Clapeyron table, chart and bullet, open items 1 and 2 (the
    direct test and the functional sensitivity, now done), the "about 20 K" budget sentence and ±26 K.
  - Changed: T_m(330 GPa) = 6,329 ± 6 K (Simon-fit error). The rev 23 numbers that did not match their
    keys now do, at the printed precision: 36.1/37.2/10.3 %, 108-128 atoms, +18.0 meV/atom, 1,974 GPa,
    0.26 %, about 15 K, +0.08/−0.15 GPa, 175-180 GPa, 5×10⁻⁵ eV/Å, +2.6 and 0.4 meV/atom.
  - New: the direct-test section and summary bullet; the LDA sensitivity paragraph and table; the
    mixed-setup and snapshot numbers in Methods; the convergence numbers (0.53 and 0.24 GPa); the
    literature row 6,306 ± 8 (potential: 6,329 ± 6); the Rao & Zhu sentence (81 ± 43 K).
  - Constant: 14.5 K per meV/atom at ΔS_m = 0.80 k_B.
- **Report rev 47.**
  - The equal-modulus volume term is printed as +0.53 meV/atom, matching `ldashift.lda_d_a3_equalB_meV`.
  - The Rao & Zhu sentence adds their DMFT U = 5.0 eV and J = 0.93 eV, their potential's mean absolute
    energy error against its labels (69.2 meV/atom) and their 4-8-atom label cells. These are new
    constants in `validation/refs/literature.json` (`Rao2025`), from arXiv:2512.25061 v4 (Appendix A,
    Sec. IV A, Appendix B).
  - "Within about 40 K of zero" is the larger |δT_m(LDA − PBE)| plus its standard error,
    `summary.lda_dTm_max_abs_within_1se_K` = 41.7 K.

## 2026-11-01 — tests fixed after the ARC GPU run

The first full run on an ARC H100 (80 GB; torch 2.12.0+cu126, mace-torch 0.3.16) gave 65 passed and
2 failed; on ARC CPU both pass. The tests changed, not the numbers:
- **`tests/test_runs_tel_model.py`** ran out of GPU memory. It now evaluates as `phaseE_coex_tel.py` did:
  float32 with the cuEquivariance kernels (e3nn if cuequivariance is missing; the output records which),
  `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`, and one model on the GPU at a time (all
  checkpoints with stage two, then the model is freed, then stage one). The pressure still comes from
  the stress, which needs the autograd pass, so the evaluation is not energies-only. The assertion is
  unchanged, and the files are opened with `with open`.
- **`test_static_compare_reproduces`** gave `eos[0][3]` = 0.24710 GPa on CUDA against 0.24738 GPa (CPU,
  float32) in `out/static_compare.json`. It now always runs on CPU, at 2e-4 relative as before. A new
  test, `test_static_compare_reproduces_on_gpu` (CUDA only), allows max(2e-4 |b|, 2e-3). That bound
  is not below half the last printed digit for 31 report numbers taken or derived from
  `out/static_compare.json`. The 2e-3 floor fails for c/a (printed to 1e-4); the 2e-4 relative term
  fails for pressures near 270 GPa printed to 0.01 GPa, for phase energies and for the MACE − DFT
  differences. Those values get absolute tolerances (`GPU_TIGHT`): 2e-3 GPa for the `eos`,
  `tel_response` and `*_P_*` pressures (the check fails from about 2.5e-3), 0.02 meV/atom for the `eos`,
  `tel_response` and `*-hcp_*` energies (stored float32 energies step by 1/512 meV/atom; the check fails
  from 0.025), 1e-3 km/s for v_s and 2e-5 for `coa_min`. The printed c/a 1.5978 is 2.8e-5 from a rounding
  edge, so a GPU failure there would mean the printed digit is not reproducible on a GPU. A third test,
  `test_gpu_tolerance_below_printed_precision` (any machine), checks the bound on every one of the 94
  report numbers concerned, one value at a time, propagated through `ledger.summary`.
- **ARC GPU run passed** at `19a52f7` (1 Nov, H100, torch 2.12.0+cu126, mace-torch 0.3.16): 69 passed, 0 skipped; the runs_tel runs loaded stage two (`coex_hcp_350GPa_ti1.08_s1`: stage two reproduces the logged P to 9.2e-5 GPa, stage one is off by −1.04 GPa).

## 2026-11-05: layout matched to water-scw-mace

- **Layout.** The ARC pipeline moved to `production_campaign/`, the models to `models/`, the datasets and
  coexistence output to `data/`, `layout.py`, `validation/` and `analysis/figures/` to `melt_curve/`, and
  the LDA test and the DMFT pilot to `studies/`. Mapping and reasons: [MIGRATION.md](MIGRATION.md).
- **No results changed.** The ARC scripts, models, datasets, outputs and the ledger are byte-identical.
  The paths in `ledger_inputs.tsv` were rewritten; its hashes and dates are unchanged.
- **Docs.** New: `docs/METHODS.md`, `TESTING.md`, and READMEs for `production_campaign/`, `melt_curve/`,
  `studies/`, `tests/` and `docs/`.
- **ARC GPU run passed** at `75c2974` after the layout change (5 Nov, job 9272282, torch 2.12.0+cu126, mace-torch 0.3.16): 69 passed, 0 skipped.
