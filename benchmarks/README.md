# benchmarks/

Rerunnable checks of the production potential, `fe_tel_e3_stagetwo.model`, against held-out data, PBE and
experiment. Anyone should be able to clone, fetch the LFS files, and reproduce these numbers.

The scripts live in [`melt_curve/validation/`](../melt_curve/validation/README.md), because `layout.py` and the
tests use those paths. Run them from that directory:

    cd melt_curve/validation
    python stats_test.py

| Check | Script | Compared against | Output (`out/`) | Cost |
|---|---|---|---|---|
| Held-out test set | `stats_test.py` | 346 held-out frames (`parity_tel.npz`) | `test_stats.json` | seconds |
| Static PBE references | `val_mace.py` | `refs/dft/*.json` (QE 6.6) | `static_compare.json` | about 6 min on 2 CPU threads |
| Re-evaluation of the test sets | `repro_check.py` | `data/training/test_Tel*.xyz`, `data/training/parity_tel.npz` | `repro_check.json` | about 8 min on 4 CPU threads |
| Same-frame T_el check | `tel_check.py` | `production_campaign/arc_scripts/logs/e3_test_final.out` | `tel_check.json` | seconds |
| 330 GPa snapshot test | `ldashift_summary.py` | `studies/lda_shift/ldashift_result.json` | `ldashift.json` | seconds |
| PBE against experiment | `pbe_vs_exp.py` | `refs/eos_pbe_cold/*.json` and published 300 K EOS | `pbe_vs_exp.json` | seconds |

`val_mace.py` and `repro_check.py` need torch and mace-torch.

**Reference values.** Every result is written to
[`ledger.json`](../melt_curve/validation/ledger.json) by `ledger.py`. The tests in
`melt_curve/validation/tests/` check that the producers reproduce their committed outputs and that the ledger
matches them. See [`TESTING.md`](../TESTING.md) and
[`melt_curve/validation/README.md`](../melt_curve/validation/README.md) for what each check does.
