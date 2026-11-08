"""record_inputs.py - record the committed inputs of the ledger that no producer here rebuilds.

Writes melt_curve/validation/ledger_inputs.tsv: one row per (input file, output it feeds) with the input's sha256
and date, its producer and where that producer runs, and the output's sha256 and date.
tests/test_staleness.py fails when an input's content differs from the record or its date is later
than the output's, i.e. when an output was not rebuilt after its input changed.

Dates: an unchanged file keeps its recorded date; a new or changed file gets today's date, or with
--dates-from SNAPSHOT HANDOFF the modification time of the same file in the ~/fe_melt snapshot (its tar
timestamps are the only record of when ARC files changed) or in the chat handoff (validation/ files).
--dates-from looks a file up by its repository path, and the snapshot keeps the ~/fe_melt paths, so it
applies only to the paths from before the 5 Nov 2026 layout change (MIGRATION.md).
After rebuilding an output from changed inputs, run this again and commit the new record.
"""
import argparse, csv, datetime, glob, hashlib, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import layout  # noqa: E402

R = os.path.relpath
ARC_GPU = "ARC, GPU, tel-mace"
# output (repository path) -> [(input glob, producer of the input, where it runs)]
FEEDS = {
    "melt_curve/validation/out/test_stats.json": [("data/training/parity_tel.npz", "parity_tel.py", ARC_GPU)],
    "melt_curve/validation/out/static_compare.json": [("melt_curve/validation/refs/dft/*.json", "melt_curve/validation/val_dft.py", "pw.x (QE 6.6), 1.5-12 min per structure"),
                                           ("models/fe_tel_e3_stagetwo.model", "phaseD_train_e3.sbatch", ARC_GPU)],
    "melt_curve/validation/out/cc_check.json": [("data/coexistence/analysis_live_all_runs_tel.json", "phaseE_analyze_live.py", ARC_GPU),
                                     ("data/coexistence/cells/coex_hcp_*.json", "phaseE_build.py", "ARC, GPU"),
                                     ("data/coexistence/runs_tel/*_ckpt.npz", "phaseE_coex_tel.py (fe_phaseE_tel.sbatch)", ARC_GPU)],
    "melt_curve/validation/out/pbe_vs_exp.json": [("melt_curve/validation/refs/eos_pbe_cold/*.json", "none in this repository (iron chat, QE 6.6)", "pw.x"),
                                       ("data/coexistence/pn_correct.json", "pn_correct.py (rebuilt by tests/test_producers.py)", "any")],
    "melt_curve/validation/out/training.json": [("production_campaign/arc_scripts/logs/fe_tel_e3_run-7.log", "phaseD_train_e3.sbatch", ARC_GPU),
                                     ("production_campaign/arc_scripts/logs/fe_tel_e3_run-7_debug.log", "phaseD_train_e3.sbatch", ARC_GPU),
                                     ("data/training/dataset_report.txt", "phaseD_dataset.py", "needs selected/ and labels/labels.jsonl (not in the repository)"),
                                     ("data/training/train_e3.xyz", "phaseD_dataset.py + mace_stamp_elec_temp", "as above")],
    "melt_curve/validation/out/tel_check.json": [("production_campaign/arc_scripts/logs/e3_test_final.out", "tel_diff_check.py", ARC_GPU)],
    "melt_curve/validation/out/coex_cells.json": [("data/coexistence/analysis_live_all_runs_tel.json", "phaseE_analyze_live.py", ARC_GPU)],
    "melt_curve/validation/out/melting_point.json": [("data/coexistence/pn_correct.json", "pn_correct.py (rebuilt by tests/test_producers.py)", "any")],
    "melt_curve/validation/out/repro_check.json": [("data/training/parity_tel.npz", "parity_tel.py", ARC_GPU),
                                        ("data/training/test_Tel*.xyz", "phaseD_dataset.py", "needs selected/ and labels/labels.jsonl")],
    "data/coexistence/pn_correct.json": [("data/coexistence/decomp_tel_tel.json", "decomp_tel.py (decomp_tel.sbatch)", ARC_GPU),
                               ("data/coexistence/analysis_live_all.json", "phaseE_analyze_live.py", "ARC, GPU"),
                               ("data/coexistence/analysis_live_all_runs_tel.json", "phaseE_analyze_live.py", ARC_GPU)],
    "melt_curve/validation/out/settings.json": [("melt_curve/validation/refs/dft/*.in", "melt_curve/validation/val_dft.py", "pw.x (QE 6.6)"),
                                     ("melt_curve/validation/refs/eos_pbe_cold/*.in", "none in this repository (iron chat)", "pw.x"),
                                     ("production_campaign/fe_melt/pseudo/*.UPF", "pslibrary 1.0.0 / ld1.x (production_campaign/fe_melt/pseudo/lda_generation)", "ld1.x"),
                                     ("data/coexistence/cells/*.json", "phaseE_build.py", "ARC, GPU"),
                                     ("production_campaign/fe_melt/labels/ldaX_*_Tel6329.in", "ldashift.py inputs (LDA inputs switched to the PBE PAW file on ARC, 27 Oct)", "ARC"),
                                     ("melt_curve/validation/refs/lda_setup/*.json", "none in this repository (iron chat, QE 6.6)", "pw.x")],
    "melt_curve/validation/out/convergence.json": [("melt_curve/validation/refs/convergence/*.json", "none in this repository (iron chat, QE 6.6)", "pw.x")],
    "melt_curve/validation/out/lda_setup.json": [("melt_curve/validation/refs/lda_setup/*.json", "none in this repository (iron chat, QE 6.6)", "pw.x")],
    "melt_curve/validation/out/ldashift.json": [("studies/lda_shift/ldashift_result.json", "ldashift.py eval", "ARC, from the 80 QE 7.5 outputs (not committed)"),
                                     ("studies/lda_shift/md_meta.json", "ldashift.py md (ldashift_md.sbatch)", ARC_GPU),
                                     ("studies/lda_shift/pbe_snapshot_energies_330.json", "ldashift.py eval outputs (80 QE 7.5 runs)", "ARC")],
    "melt_curve/validation/out/melting_entropy.json": [("studies/lda_shift/pbe_snapshot_energies_330.json", "ldashift.py eval outputs (80 QE 7.5 runs)", "ARC")],
    "melt_curve/validation/ledger.json": [("melt_curve/validation/refs/literature.json", "none (committed constants)", "edited by hand")],
    "melt_curve/figures/inputs/f_T_extract.json": [("data/coexistence/runs_tel/coex_hcp_330GPa_ti1.00_s0_ckpt.npz", "phaseE_coex_tel.py", ARC_GPU)],
    "melt_curve/figures/inputs/f_cell_q6.json": [("data/coexistence/cells/coex_hcp_330GPa.extxyz", "phaseE_build.py", "ARC, GPU")],
}
NOT_REBUILT_OUTPUTS = {"melt_curve/figures/inputs/f_T_extract.json": "figs_coex.py (not in build_figures.py)",
                       "melt_curve/figures/inputs/f_cell_q6.json": "figs_coex.py (not in build_figures.py)"}
TSV = os.path.join(HERE, "ledger_inputs.tsv")
COLS = ["input", "input_sha256", "input_date", "input_producer", "runs_on", "output", "output_sha256", "output_date"]


def sha(p):
    return hashlib.sha256(open(os.path.join(layout.ROOT, p), "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dates-from", nargs=2, metavar=("SNAPSHOT", "HANDOFF"))
    a = ap.parse_args()
    old = {}
    if os.path.exists(TSV):
        for r in csv.DictReader(open(TSV), delimiter="\t"):
            old[r["input"]] = (r["input_sha256"], r["input_date"])
            old[r["output"]] = (r["output_sha256"], r["output_date"])
    today = datetime.datetime.now().isoformat(timespec="seconds")

    def date(p, h):
        if p in old and old[p][0] == h:
            return old[p][1]
        if a.dates_from:
            snap, hand = a.dates_from
            for base, rel in ((snap, p), (hand, p), (hand, p.replace("melt_curve/figures/", "figures_src/"))):
                f = os.path.join(base, rel)
                if os.path.exists(f) and hashlib.sha256(open(f, "rb").read()).hexdigest() == h:
                    return datetime.datetime.fromtimestamp(os.path.getmtime(f)).isoformat(timespec="seconds")
        return today

    rows = []
    for out, ins in FEEDS.items():
        ho = sha(out)
        do = date(out, ho)
        for pat, prod, where in ins:
            files = sorted(R(f, layout.ROOT) for f in glob.glob(os.path.join(layout.ROOT, pat)))
            assert files, pat
            for f in files:
                h = sha(f)
                rows.append([f, h, date(f, h), prod, where, out, ho, do])
    with open(TSV, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"{len(rows)} rows -> {R(TSV, layout.ROOT)}")


if __name__ == "__main__":
    main()
