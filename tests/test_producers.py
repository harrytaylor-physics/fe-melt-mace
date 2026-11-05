"""The producers added for the ledger (melt_curve/validation/train_table.py, tel_check.py, coex_cells.py,
melting_point.py, settings.py) and pn_correct.py reproduce their committed outputs.
Producers whose LFS inputs are absent (git lfs pull not run) skip."""
import json, os, shutil, subprocess, sys
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VAL = os.path.join(ROOT, "melt_curve", "validation")
sys.path.insert(0, os.path.join(ROOT, "melt_curve"))
import layout  # noqa: E402


def lfs(path):
    return os.path.exists(path) and open(path, "rb").read(7) != b"version"


def close(a, b, tol, path=""):
    if isinstance(a, dict):
        assert set(a) == set(b), f"{path}: keys {sorted(set(a) ^ set(b))}"
        for k in a:
            close(a[k], b[k], tol, f"{path}.{k}")
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            close(x, y, tol, f"{path}[{i}]")
    elif isinstance(a, float) or isinstance(b, float):
        assert abs(a - b) <= tol * max(1.0, abs(b)), f"{path}: {a} vs {b}"
    else:
        assert a == b, f"{path}: {a} vs {b}"


def run(script, tmp):
    e = dict(os.environ, FE_VAL_OUT=str(tmp), FE_VAL_LEDGER=str(tmp / "ledger.json"))
    r = subprocess.run([sys.executable, script], cwd=VAL, env=e, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]


@pytest.mark.parametrize("script,out", [("tel_check.py", "tel_check.json"), ("coex_cells.py", "coex_cells.json"),
                                        ("melting_point.py", "melting_point.json"), ("settings.py", "settings.json")])
def test_producer_reproduces(script, out, tmp_path):
    run(script, tmp_path)
    close(json.load(open(tmp_path / out)), json.load(open(os.path.join(layout.VAL_OUT, out))), 1e-9)


@pytest.mark.skipif(not (lfs(layout.TRAIN_SET) and lfs(layout.MODEL)), reason="data/training/train_e3.xyz or the model not present (LFS)")
def test_train_table_reproduces(tmp_path):
    run("train_table.py", tmp_path)
    close(json.load(open(tmp_path / "training.json")), json.load(open(os.path.join(layout.VAL_OUT, "training.json"))), 1e-12)


def test_pn_correct_reproduces(tmp_path):
    """pn_correct.py (ARC, verbatim) writes phaseE/pn_correct.json relative to its working directory."""
    (tmp_path / "phaseE").mkdir()
    for f in ("decomp_tel_tel.json", "analysis_live_all.json", "analysis_live_all_runs_tel.json"):
        shutil.copy(os.path.join(layout.PHASEE, f), tmp_path / "phaseE" / f)
    r = subprocess.run([sys.executable, os.path.join(layout.ARC, "pn_correct.py")], cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    close(json.load(open(tmp_path / "phaseE" / "pn_correct.json")), json.load(open(layout.PN_CORRECT)), 1e-6)


def test_literature_matches_pbe_vs_exp():
    """The experimental EOS parameters in refs/literature.json are the ones pbe_vs_exp.py uses."""
    import numpy as np
    sys.path.insert(0, VAL)
    import pbe_vs_exp as p
    L = json.load(open(layout.LITERATURE))["pressure_scales_300K"]
    V = np.linspace(6.8, 10.6, 12)
    for k, c in L.items():
        V0 = c["V0_bohr3"] * p.B3 if "V0_bohr3" in c else c["V0_A3"]
        f = p.vinet if c["form"] == "vinet" else p.bm3
        assert np.allclose(p.EXP[k](V), f(V, V0, c["K0_GPa"], c["Kp"]), rtol=1e-12), k


def test_literature_matches_figure_producer():
    """figs_results.py draws literature.json; its intermediate results.json holds the same points."""
    L = {x["id"]: x for x in json.load(open(layout.LITERATURE))["melting_curve_figure"]}
    res = json.load(open(os.path.join(layout.FIG_INPUTS, "results.json")))["lit"]
    pts = sorted(tuple(p) for x in L.values() for p in x["points"])
    assert sorted((r["P"], r["T"], r["err"]) for r in res) == pts
