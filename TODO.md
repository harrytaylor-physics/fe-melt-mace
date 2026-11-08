# TODO — technical (not physics)

Priority order.

1. **Add a licence.** The repository has no licence file, so others have no stated right to reuse the code or
   the data.
2. **Move ARC checkouts to the rewritten history.** `~/tel_env` is an editable install of tel-mace at `5168bc0`;
   the same files are now at `69047be`. Any clone of this repository needs resetting to the new `main`.
3. **One command for the benchmarks.** [`benchmarks/README.md`](benchmarks/README.md) lists the checks one by
   one. Add a `run_benchmarks.py` that runs them and compares each result with `ledger.json`.
4. **CI.** Add a GitHub Action running pytest. Tests that need the LFS files or torch already skip without them.
5. **Pin versions.** `requirements.txt` leaves numpy, scipy, matplotlib and pytest unpinned and sets only a
   minimum for ase. The versions used on ARC are listed under Dependencies in [`TESTING.md`](TESTING.md).
6. **Move the convergence checks into `convergence/`.** `conv_check.py`, `lda_setup_check.py` and their inputs
   are still in `melt_curve/validation/`. Moving them means updating `layout.py` and `paths.py`.
