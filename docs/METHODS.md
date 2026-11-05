# Methods

What the pipeline does, stage by stage. The numerical settings of each check are in
[`melt_curve/validation/README.md`](../melt_curve/validation/README.md). The scripts are in
[`production_campaign/`](../production_campaign/README.md).

1. **QE labels (Phase C).** MD sampling of bcc, fcc, hcp, liquid and stacking-faulted iron with a MACE
   foundation model. The selected 128-atom frames (108 for fcc) are labelled with PBE in Quantum
   ESPRESSO at T_el = 4,000, 6,000 and 8,000 K, using Fermi-Dirac smearing and pslibrary 1.0.0 PAW
   with 16 valence electrons. Scripts: `phaseC_*.py`, `labels/`.
2. **Training (Phase D).** Datasets (`phaseD_dataset.py`), then fine-tuning of the foundation model
   (`phaseD_train*.sbatch`) and gates (`phaseD_gates.py`). The production potential
   `models/fe_tel_e3_stagetwo.model` takes T_el as a learned per-structure embedding. It was
   trained on all three T_el (`data/training/train_e3.xyz`) and saved after the SWA stage two. Provenance:
   [`models/README.md`](../models/README.md).
3. **Coexistence (Phase E).** 12,096-atom hcp-liquid cells at 100-350 GPa, 6 NVE runs of 200 ps per
   cell, each cell at the T_el of its melting temperature (`phaseE_build.py`, `phaseE_coex_tel.py`,
   `fe_phaseE_tel.sbatch`, output in `data/coexistence/runs_tel/`).
4. **Analysis.**
   - Melting temperature and entropy per cell: `phaseE_analyze_live.py` →
     `data/coexistence/analysis_live_all_runs_tel.json`.
   - Interface-normal pressure and the Simon fit: `pn_correct.py` → `data/coexistence/pn_correct.json`.
   - T_el decomposition: `decomp_tel.py`.
   - Validation against PBE and experiment, and the ledger:
     [`melt_curve/validation/`](../melt_curve/validation/README.md).
   - Figures: [`docs/figures/`](figures/README.md).
5. **330 GPa snapshot test** (`ldashift.py`, `labels/ldaX_*`, [`studies/lda_shift/`](../studies/lda_shift/README.md)).
   20 hcp and 20 liquid snapshots from the potential at 330 GPa and 6,329 K are labelled with PBE and with
   LDA in QE 7.5. The test gives two results:
   - the potential's error against PBE at the state point of the result;
   - the LDA − PBE shift of T_m by free-energy perturbation.

   The LDA runs use the PBE PAW file with `input_dft = 'PZ'`, because the generated PZ dataset gives a
   divergent SCF (`production_campaign/fe_melt/pseudo/lda_generation/README.md`). Report numbers:
   `melt_curve/validation/ldashift_summary.py`.

The EOS anchor (`eos_anchor.py`, `labels/eosA_*`) was not run to completion and is not used: its 28
runs were cancelled.
