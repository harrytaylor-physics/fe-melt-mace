# tests/

Test suite for the producers, the ledger and the layout.

Repository test suite. How to run it and when tests skip: [`../TESTING.md`](../TESTING.md).

| File | Covers |
|---|---|
| `test_producers.py` | the new producers and `pn_correct.py` reproduce their outputs |
| `test_ldashift_inputs.py` | `ldashift.lda_input()` reproduces the committed LDA inputs from the PBE ones |
| `test_repository.py` | the report-number list matches the ledger; no committed code holds an absolute path; the paths in `melt_curve/layout.py` exist |
| `test_staleness.py` | the staleness guard (`melt_curve/validation/ledger_inputs.tsv`) |
| `test_runs_tel_model.py` | on a GPU, which model the T_el-conditioned runs loaded |
