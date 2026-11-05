# melt_curve/

The code that turns the campaign's output into the reported numbers and figures. It runs on any machine.

    layout.py     every location this code reads or writes, in one place
    validation/   checks against PBE and experiment, producers, ledger.json, report_numbers.tsv
    figures/      figure producers and their intermediate json; they write to docs/figures/

Inputs come from `../data/`, `../models/`, `../studies/lda_shift/` and the logs in
`../production_campaign/arc_scripts/logs/`. Every path in `layout.py` is relative to the repository root and
can be overridden by the environment variable named beside it.

**Simon fit.** `pn_correct.py` is an ARC script, kept verbatim, that reads and writes `phaseE/` relative to
its working directory. To rerun it, make a directory with a `phaseE/` folder holding `decomp_tel_tel.json`,
`analysis_live_all.json` and `analysis_live_all_runs_tel.json` from `../data/coexistence/`, and run
`production_campaign/fe_melt/pn_correct.py` from there. `tests/test_producers.py` does this.

More: [`validation/README.md`](validation/README.md), [`figures/README.md`](figures/README.md).
