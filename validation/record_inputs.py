"""record_inputs.py - record the committed inputs of the ledger that no producer here rebuilds.

Writes validation/ledger_inputs.tsv: one row per (input file, output it feeds) with the input's sha256
and date, its producer and where that producer runs, and the output's sha256 and date.
tests/test_staleness.py fails when an input's content differs from the record or its date is later
than the output's, i.e. when an output was not rebuilt after its input changed.

Dates: an unchanged file keeps its recorded date; a new or changed file gets today's date, or with
--dates-from SNAPSHOT HANDOFF the modification time of the same file in the ~/fe_melt snapshot (its tar
timestamps are the only record of when ARC files changed) or in the chat handoff (validation/ files).
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
    "validation/out/test_stats.json": [("phaseD/parity_tel.npz", "parity_tel.py", ARC_GPU)],
    "validation/out/static_compare.json": [("validation/refs/dft/*.json", "validation/val_dft.py", "pw.x (QE 6.6), 1.5-12 min per structure"),
                                           ("phaseD/models/fe_tel_e3_stagetwo.model", "phaseD_train_e3.sbatch", ARC_GPU)],
    "validation/out/cc_check.json": [("phaseE/analysis_live_all_runs_tel.json", "phaseE_analyze_live.py", ARC_GPU),
                                     ("phaseE/cells/coex_hcp_*.json", "phaseE_build.py", "ARC, GPU"),
                                     ("phaseE/runs_tel/*_ckpt.npz", "phaseE_coex_tel.py (fe_phaseE_tel.sbatch)", ARC_GPU)],
    "validation/out/pbe_vs_exp.json": [("validation/refs/eos_pbe_cold/*.json", "none in this repository (iron chat, QE 6.6)", "pw.x"),
                                       ("phaseE/pn_correct.json", "pn_correct.py (rebuilt by tests/test_producers.py)", "any")],
    "validation/out/training.json": [("logs/fe_tel_e3_run-7.log", "phaseD_train_e3.sbatch", ARC_GPU),
                                     ("logs/fe_tel_e3_run-7_debug.log", "phaseD_train_e3.sbatch", ARC_GPU),
                                     ("phaseD/dataset_report.txt", "phaseD_dataset.py", "needs selected/ and labels/labels.jsonl (not in the repository)"),
                                     ("phaseD/train_e3.xyz", "phaseD_dataset.py + mace_stamp_elec_temp", "as above")],
    "validation/out/tel_check.json": [("logs/e3_test_final.out", "tel_diff_check.py", ARC_GPU)],
    "validation/out/coex_cells.json": [("phaseE/analysis_live_all_runs_tel.json", "phaseE_analyze_live.py", ARC_GPU)],
    "validation/out/melting_point.json": [("phaseE/pn_correct.json", "pn_correct.py (rebuilt by tests/test_producers.py)", "any")],
    "validation/out/repro_check.json": [("phaseD/parity_tel.npz", "parity_tel.py", ARC_GPU),
                                        ("phaseD/test_Tel*.xyz", "phaseD_dataset.py", "needs selected/ and labels/labels.jsonl")],
    "phaseE/pn_correct.json": [("phaseE/decomp_tel_tel.json", "decomp_tel.py (decomp_tel.sbatch)", ARC_GPU),
                               ("phaseE/analysis_live_all.json", "phaseE_analyze_live.py", "ARC, GPU"),
                               ("phaseE/analysis_live_all_runs_tel.json", "phaseE_analyze_live.py", ARC_GPU)],
    "validation/out/settings.json": [("validation/refs/dft/*.in", "validation/val_dft.py", "pw.x (QE 6.6)"),
                                     ("validation/refs/eos_pbe_cold/*.in", "none in this repository (iron chat)", "pw.x"),
                                     ("pseudo/*.UPF", "pslibrary 1.0.0 / ld1.x (pseudo/lda_generation)", "ld1.x"),
                                     ("phaseE/cells/*.json", "phaseE_build.py", "ARC, GPU")],
    "validation/ledger.json": [("validation/refs/literature.json", "none (committed constants)", "edited by hand")],
    "analysis/figures/inputs/f_T_extract.json": [("phaseE/runs_tel/coex_hcp_330GPa_ti1.00_s0_ckpt.npz", "phaseE_coex_tel.py", ARC_GPU)],
    "analysis/figures/inputs/f_cell_q6.json": [("phaseE/cells/coex_hcp_330GPa.extxyz", "phaseE_build.py", "ARC, GPU")],
}
NOT_REBUILT_OUTPUTS = {"analysis/figures/inputs/f_T_extract.json": "figs_coex.py (not in build_figures.py)",
                       "analysis/figures/inputs/f_cell_q6.json": "figs_coex.py (not in build_figures.py)"}
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
            for base, rel in ((snap, p), (hand, p), (hand, p.replace("analysis/figures/", "figures_src/"))):
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
