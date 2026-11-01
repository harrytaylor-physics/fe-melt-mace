"""build_figures.py - rebuild docs/figures/ (PDF + 200-dpi PNG) and the intermediate json in inputs/.

Runs the producers below in order and checks that validation/ledger.json and validation/out/*.json
are byte-for-byte unchanged afterwards (the figures read the ledger's sources; they never write them).
  python analysis/figures/build_figures.py            # from the repository root, about 2 min
Needs the LFS files (phaseD/train_e3.xyz, test sets, parity_tel.npz, phaseE/runs_tel/*_ckpt.npz).
figs_coex.py, figs_phase.py, figs_triple.py, figs_lit.py, figs_dft.py, figs_stack.py, figs_illus.py and
compare_curves_v1.py made illustrations for a talk; they are wired to the layout but not part of this build.
"""
import glob, hashlib, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
import layout  # noqa: E402

STEPS = ["parse_train.py", "figs_parity.py", "figs_results.py", "figs_train.py", "figs_validation.py"]


def digest():
    files = [layout.LEDGER] + sorted(glob.glob(os.path.join(layout.VAL_OUT, "*.json")))
    return {f: hashlib.sha256(open(f, "rb").read()).hexdigest() for f in files}


def main():
    before = digest()
    for s in STEPS:
        print(f"== {s}", flush=True)
        r = subprocess.run([sys.executable, s], cwd=HERE, capture_output=True, text=True)
        if r.returncode:
            sys.exit(r.stdout[-2000:] + r.stderr[-2000:])
    after = digest()
    changed = [f for f in before if before[f] != after.get(f)]
    if changed:
        sys.exit(f"the figure build changed ledger sources: {changed}")
    print(f"docs/figures rebuilt; ledger.json and {len(before) - 1} validation outputs unchanged")


if __name__ == "__main__":
    main()
