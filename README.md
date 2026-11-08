# fe-melt-mace

The melting curve of iron to 370 GPa from two-phase coexistence simulations with a MACE potential
trained on PBE (Mermin free energies) whose input includes the electronic temperature T_el.

**Result** (keys of [`melt_curve/validation/ledger.json`](melt_curve/validation/ledger.json)):
- **T_m(330 GPa) = 6,329 ± 6 K from the potential.** The value is `keys.melting_point.Tm_330`, from the
  Simon fit in `data/coexistence/pn_correct.json`. The error is statistical, the fit's standard error at
  330 GPa: `summary.Tm_330_stat_err_K`. dT_m/dP = 8.41 K/GPa (`keys.melting_point.dTdP_330`).
- **T_m(330 GPa) = 6,306 ± 8 K for PBE.** A direct PBE test on 40 snapshots at 330 GPa and 6,329 K puts the
  potential's T_m 23 ± 6 K above PBE (`keys.ldashift.mace_Tm_total_K`). The PBE value is
  `summary.Tm_PBE_330_K`. Its error, `summary.Tm_PBE_330_err_K`, is the Simon-fit error and the
  direct test's error added in quadrature.

Every number quoted from this work has a ledger key; the list is
[`melt_curve/validation/report_numbers.tsv`](melt_curve/validation/report_numbers.tsv).

## What you can do with it

- Evaluate the T_el-conditioned iron potential at a chosen T_el (`melt_curve/validation/load_model.py`)
- Reproduce the Simon fit from the coexistence analysis
- Rerun the checks of the potential against PBE and of PBE against experiment
- Rebuild the ledger and the figures from the committed inputs

## Quick start

    git lfs install                  # required once; models and datasets are LFS objects
    git clone git@github.com:harrytaylor-physics/fe-melt-mace.git
    cd fe-melt-mace
    pip install -r requirements.txt
    python -m pytest -q -p no:cacheprovider

Without the LFS files the binaries are 132-byte pointers and the tests that need them skip. To fetch them:
`git lfs install --local && git lfs pull`. Conditions of use are in [TESTING.md](TESTING.md).

## Map

| Where | What |
|---|---|
| [models/](models/README.md) | MACE potentials with T_el input, and their provenance. |
| [benchmarks/](benchmarks/README.md) | Rerunnable checks of the potential against held-out data, PBE and experiment. |
| [tests/](tests/README.md) | Test suite for the producers, the ledger and the layout. |
| [convergence/](convergence/README.md) | Convergence checks: what was converged, how, and the result. |
| [studies/](studies/README.md) | Side studies: the 330 GPa LDA - PBE test, a DFT+DMFT pilot. |
| [docs/](docs/README.md) | Methods reference and figures of the iron melting curve. |
| [melt_curve/](melt_curve/README.md) | Validation against PBE and experiment, ledger and figures. |
| [production_campaign/](production_campaign/README.md) | ARC pipeline as run: QE labels, training, coexistence runs. |
| [data/](data/README.md) | Training and test sets with QE labels, and coexistence runs. |
| [MIGRATION.md](MIGRATION.md) | how `~/fe_melt` came in and what moved where |
| [TODO.md](TODO.md) | technical task list |
| [TESTING.md](TESTING.md) | how to run the tests, what they cover, conditions of use, dependencies, reproducing |

## Status

| Component | Status |
|---|---|
| `fe_tel_e3_stagetwo.model` (T_el-conditioned potential) | production; SWA stage two; trained at T_el = 4,000, 6,000 and 8,000 K |
| `fe_ft_Tel{4000,6000,8000}_f32.model` (fixed-T_el potentials) | v1; the 6,000 K model ran the campaign in `data/coexistence/runs/`, which set each production cell's T_el |
| Melting curve, 100-350 GPa | complete; 12,096-atom cells, 6 NVE runs of 200 ps per cell |
| Validation against PBE and experiment | complete; every number in the ledger |
| 330 GPa snapshot test | complete ([studies/lda_shift/](studies/lda_shift/README.md)) |
| EOS anchor (`eos_anchor.py`) | not used; its 28 runs were cancelled |
| DFT+DMFT pilot | exploratory; not used in any result ([studies/dmft_pilot/](studies/dmft_pilot/README.md)) |

## Methods in one paragraph

MD sampling of bcc, fcc, hcp, liquid and stacking-faulted iron with a MACE foundation model; PBE labels in
Quantum ESPRESSO at T_el = 4,000, 6,000 and 8,000 K (Fermi-Dirac smearing, pslibrary 1.0.0 PAW with 16
valence electrons); fine-tuning of the foundation model with T_el as a learned per-structure embedding;
12,096-atom hcp-liquid coexistence at 100-350 GPa, each cell at the T_el of its melting temperature;
interface-normal pressure correction and a Simon fit. Full details: [docs/METHODS.md](docs/METHODS.md).

## References

**MACE (architecture and foundation models)**
- I. Batatia, D. P. Kovács, G. N. C. Simm, C. Ortner and G. Csányi, Adv. Neural Inf. Process. Syst. 35, 11423 (2022)
- I. Batatia et al., J. Chem. Phys. 163, 184110 (2025)
- A. D. Kaplan et al., arXiv:2503.04070 (2025) - MatPES

**Electronic structure**
- N. D. Mermin, Phys. Rev. 137, A1441 (1965) - finite-temperature DFT
- P. Giannozzi et al., J. Phys.: Condens. Matter 21, 395502 (2009); 29, 465901 (2017) - Quantum ESPRESSO
- J. P. Perdew, K. Burke and M. Ernzerhof, Phys. Rev. Lett. 77, 3865 (1996) - PBE
- P. E. Blöchl, Phys. Rev. B 50, 17953 (1994) - PAW
- A. Dal Corso, Comput. Mater. Sci. 95, 337 (2014) - pslibrary
- H. J. Monkhorst and J. D. Pack, Phys. Rev. B 13, 5188 (1976) - k-point meshes

**Coexistence and analysis**
- J. R. Morris, C. Z. Wang, K. M. Ho and C. T. Chan, Phys. Rev. B 49, 3109 (1994) - two-phase coexistence
- P. J. Steinhardt, D. R. Nelson and M. Ronchetti, Phys. Rev. B 28, 784 (1983) - bond-order parameters
- W. Lechner and C. Dellago, J. Chem. Phys. 129, 114707 (2008) - averaged bond-order parameters
- F. Simon and G. Glatzel, Z. Anorg. Allg. Chem. 178, 309 (1929) - Simon fit

**Iron melting: comparisons**
- S. Anzellini, A. Dewaele, M. Mezouar, P. Loubeyre and G. Morard, Science 340, 464 (2013)
- D. Alfè, Phys. Rev. B 79, 060101(R) (2009)
- E. Sola and D. Alfè, Phys. Rev. Lett. 103, 078501 (2009)
- Y. Sun et al., Geophys. Res. Lett. 50, e2022GL102447 (2023)
- R. Rao and L. Zhu, arXiv:2512.25061 (v1 Dec 2025; v4 Aug 2026)

The full list used in the figures and the report is in
[`melt_curve/figures/refs.py`](melt_curve/figures/refs.py).

---
Maintainer: Harry Taylor (Oxford).
