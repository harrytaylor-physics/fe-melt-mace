# pseudo_lda/ — LDA (Perdew-Zunger) PAW for Fe, pslibrary 1.0.0 construction

The QE pseudopotential site lists Fe LDA (pz) PAW files only for pslibrary 0.2.1; version 1.0.0
has PBE and PBEsol only. The LDA - PBE test (`ldashift.py`) and `eos_anchor.py` need an LDA file
built exactly like the PBE file used for all training labels, so the LDA - PBE difference
contains the functional only.

| File | What it is |
|---|---|
| `Fe.pz-spn-kjpaw_psl.1.0.0.in` | ld1.x input: the pslibrary 1.0.0 Fe input with `dft='PZ'` |
| `Fe.pbe-spn-kjpaw_psl.1.0.0.in` | the same input with `dft='PBE'` (the official recipe) |
| `Fe.pz-spn-kjpaw_psl.1.0.0.UPF` | generated with ld1.x from QE 6.6, 22 Oct 2026; installed in `~/fe_melt/pseudo/` on ARC |
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
