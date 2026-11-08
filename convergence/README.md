# convergence/

The convergence checks behind the iron melting curve: what was converged, how, and the result.

The scripts and their inputs live in [`melt_curve/validation/`](../melt_curve/validation/README.md), because
`layout.py` and the tests use those paths. Run them from that directory:

    cd melt_curve/validation
    python conv_check.py          # writes out/convergence.json
    python lda_setup_check.py     # writes out/lda_setup.json

| Stage | Converge | Script | Inputs | Result |
|---|---|---|---|---|
| Electronic structure | plane-wave cutoff, 100/1,000 -> 130/1,300 Ry | `conv_check.py` | `refs/convergence/*.json` | static pressure changes by -0.53 GPa |
| Electronic structure | k mesh, 12x12x8 -> 16x16x10 | `conv_check.py` | `refs/convergence/*.json` | static pressure changes by +0.24 GPa |
| Electronic structure | LDA (PZ) on the PBE PAW file against the generated PZ file | `lda_setup_check.py` | `refs/lda_setup/*.json` | hcp pressure +0.22 GPa; fcc - hcp energy -0.14 meV/atom |

The cutoff and k-mesh checks are static PBE on hcp iron at 7.0 A^3/atom (c/a 1.60, Marzari-Vanderbilt 0.01 Ry,
QE 6.6, pslibrary 1.0.0 PBE PAW). The LDA setup check uses perfect hcp and fcc lattices at 7.0 A^3/atom.

There is no dedicated convergence script here for training-set size or coexistence cell size.
`melt_curve/figures/figs_lit.py` labels finite size "not yet measured".
