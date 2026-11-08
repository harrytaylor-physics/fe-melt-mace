# fe-melt-mace

The melting curve of hcp iron from 100 to 370 GPa, from 12,096-atom solid–liquid coexistence
simulations with a MACE potential fine-tuned to PBE. The potential takes the electronic temperature
T_el as an input, so one model covers the Mermin free-energy surface from 4,000 to 8,000 K.

## Result

**T_m(330 GPa) = 6,329 ± 6 K** from the potential, with dT_m/dP = 8.41 K/GPa.

Value: `keys.melting_point.Tm_330`, from the Simon fit in `data/coexistence/pn_correct.json`.
Error: statistical, the standard error of the fit at 330 GPa (`summary.Tm_330_stat_err_K`).
Slope: `keys.melting_point.dTdP_330`.

How far this value can move, and why:

| Test | Shift of T_m(330 GPa) | Ledger key |
|---|---|---|
| Potential against PBE: free-energy perturbation on 40 snapshots (20 hcp, 20 liquid, 128 atoms) relabelled with PBE | potential is 19.4 ± 4.8 K above PBE (PBE: 6,310 ± 8 K) | `keys.ldashift.mace_Tm_total_K` |
| PBE pressure against three 300 K experimental equations of state at 7.0 Å³/atom | −110 to +60 K | `summary.Tm_true_minus_pbe_K` |
| LDA instead of PBE, same 40 snapshots | +11 ± 14 K | `keys.ldashift.lda_dTm_K` |

Melting entropy and volume at 330 GPa: the PBE snapshot enthalpies give ΔS_m = 0.95 ± 0.04 k_B
(`keys.melting_entropy.dS_kB`, `keys.melting_entropy.dS_err_kB`) and the bulk volumes give ΔV_m/V = 1.6%
(`keys.melting_entropy.dV_over_V_330_pct`).

Melting curve: Simon fit T_m = T_1 [1 + (P − 100 GPa)/a]^(1/c) to seven hcp cells between 104 and
371 GPa, with T_1 = 3,846 ± 11 K, a = 120 ± 6 GPa, c = 2.15 ± 0.06 (`keys.melting_point.fit`, errors
`keys.melting_point.fit_param_se`). fcc at 91 GPa melts at 3,658 K
(`keys.coex_cells.cells.coex_fcc_100GPa.Tm`).

Every number in this section has a ledger key in `melt_curve/validation/ledger.json`; the list
is `melt_curve/validation/report_numbers.tsv`.

## What is in the repository

- **The potential** (`models/`). MACE, fine-tuned from the foundation model `MACE-matpes-pbe-omat-ft`.
  Trained on hcp, fcc, bcc, liquid and stacking-faulted iron at 55–458 GPa and 2,200–9,400 K, with T_el from 4,000 to 8,000 K.
  Validation errors 1.5 meV/atom (energy) and 92 meV/Å (forces).
- **The training data** (`data/`). 6,955 PBE calculations (Mermin free energy and forces) on 2,422
  structures of 108–128 atoms, at T_el = 4,000, 6,000 and 8,000 K, with the held-out test split.
- **The coexistence results** (`data/coexistence/`). 48 NVE runs of 12,032–12,096 atoms, 9.1 ns in
  total: temperature, pressure and solid fraction of every run, the interface-normal pressure of each
  cell, and the Simon fit.
- **The pipeline as run** (`production_campaign/`). Foundation-model sampling, Quantum ESPRESSO inputs,
  training, coexistence build and analysis, as run on the ARC cluster.
- **The checks** (`melt_curve/`, `studies/`). Free-energy perturbation of the potential and of LDA
  against PBE, PBE against experimental equations of state, the response to T_el, the ledger and the
  figures.

## What you can do with it

- Run molecular dynamics of solid or liquid iron at core conditions, at any T_el in the trained range,
  including T_el = T (`melt_curve/validation/load_model.py`).
- Train another potential on the PBE labels, or extend them to other electronic temperatures.
- Test another potential or functional against PBE at 330 GPa with the same snapshots and
  free-energy-perturbation script.
- Extend the melting curve to other pressures with the coexistence build and analysis scripts.
- Check the ledger against its committed inputs (`pytest`, see [`TESTING.md`](TESTING.md)) and rebuild the
  figures (`python melt_curve/figures/build_figures.py`).

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
