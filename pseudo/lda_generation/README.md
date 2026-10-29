# pseudo_lda/ — LDA (Perdew-Zunger) PAW for Fe, pslibrary 1.0.0 construction

The QE pseudopotential site lists Fe LDA (pz) PAW files only for pslibrary 0.2.1; version 1.0.0
has PBE and PBEsol only. This LDA file was generated for the LDA - PBE test (`ldashift.py`) and
`eos_anchor.py`, built exactly like the PBE file used for all training labels, so that the LDA - PBE
difference would contain the functional only.

**The generated PZ dataset is not used by any result.** It gives a divergent SCF in displaced and
liquid configurations (below), so the LDA snapshot runs of the 330 GPa test use the official PBE file
with `input_dft = 'PZ'`; `validation/lda_setup_check.py` measures what that mixed setup changes.

| File | What it is |
|---|---|
| `Fe.pz-spn-kjpaw_psl.1.0.0.in` | ld1.x input: the pslibrary 1.0.0 Fe input with `dft='PZ'` |
| `Fe.pbe-spn-kjpaw_psl.1.0.0.in` | the same input with `dft='PBE'` (the official recipe) |
| `Fe.pz-spn-kjpaw_psl.1.0.0.UPF` | generated with ld1.x from QE 6.6, 22 Oct 2026; installed in `~/fe_melt/pseudo/` on ARC (here `pseudo/`). Not used by any result |
| `Fe.pz.ld1.out`, `Fe.pbe.ld1.out` | ld1.x logs |
| `upf_compare.py` | compares two UPF files array by array |
| `upf_compare_ARC_2026-10-22.txt` | official PBE file on ARC vs PBE regenerated with the same ld1.x |

Check (22 Oct 2026, on ARC): the regenerated PBE file reproduces the official one (ld1.x 6.3):
total pseudo-energy -328.60578004 vs -328.60578006 Ry; all 67 numeric arrays used by pw.x
(local potential, projectors, D_ij, augmentation, core charge, PAW data) agree to <= 2.2e-6
relative. The script reported MISMATCH for two reasons that do not affect pw.x results:
PP_CHI.2 has the opposite overall sign (the sign of a pseudo-atomic orbital is arbitrary; pw.x
uses it only for starting wavefunctions and projections), and the official file carries three
GIPAW core-orbital arrays used only by gipaw.x.

## Divergent SCF with the generated PZ dataset (26 Oct 2026)

All 40 LDA snapshot runs of the 330 GPa test (128 atoms, hcp and liquid, 6,329 K) with this file
timed out with a diverging SCF: the first step's estimated SCF accuracy was about 5,500 Ry, with
negative charge. A local cross test on a displaced 16-atom hcp cell (7.0 Å³/atom, displacements
σ = 0.15 Å) located the cause in the dataset, not the functional:

| PAW file | Functional | First-step SCF accuracy (Ry) |
|---|---|---|
| `Fe.pbe-spn-kjpaw_psl.1.0.0.UPF` (official) | PBE | 4.49 |
| `Fe.pbe-spn-kjpaw_psl.1.0.0.UPF` (official) | PZ | 4.47 |
| `Fe.pz-spn-kjpaw_psl.1.0.0.UPF` (generated) | PZ | 302.8, negative charge |
| `Fe.pz-spn-kjpaw_psl.1.0.0.UPF` (generated) | PBE | 8.63, negative charge |

The perfect-crystal tests of 22 Oct had converged, which hid the problem. On 27 Oct the 40 LDA inputs
were switched to the official PBE file with `input_dft = 'PZ'` (`labels/ldaX_*_LDA_Tel6329.in`).
