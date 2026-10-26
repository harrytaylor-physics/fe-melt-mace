"""Tests for validation/: committed outputs reproduce from their producers, and the ledger
matches the committed outputs. Producers whose inputs are not present skip."""
import glob, json, os, subprocess, sys
import numpy as np
import pytest

VAL = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, VAL)
import paths  # noqa: E402


def present(path):
    """True if the file exists and is not a git-LFS pointer (git lfs pull not run)."""
    return os.path.exists(path) and open(path, "rb").read(7) != b"version"


def run(script, tmp, **env):
    e = dict(os.environ, FE_VAL_OUT=str(tmp), FE_VAL_LEDGER=str(tmp / "ledger.json"), **env)
    r = subprocess.run([sys.executable, script], cwd=VAL, env=e, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    return r.stdout


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


def committed(name):
    return json.load(open(os.path.join(VAL, "out", name)))


def test_dft_refs_complete_and_sign_convention():
    fs = sorted(glob.glob(os.path.join(paths.DFT_REFS, "*.json")))
    assert len(fs) == 26
    for f in fs:
        r = json.load(open(f))
        for k in ("F_per_atom", "P", "sigma_qe", "forces", "cell", "positions", "tel", "V"):
            assert k in r, (f, k)
        # QE prints the stress with compression positive; its trace / 3 is the 'P=' value
        assert abs(np.trace(np.array(r["sigma_qe"])) / 3 - r["P"]) < 0.02, f


def test_ledger_matches_committed_outputs(tmp_path):
    for f in glob.glob(os.path.join(VAL, "out", "*.json")):
        (tmp_path / os.path.basename(f)).write_text(open(f).read())
    run("ledger.py", tmp_path)
    new = json.load(open(tmp_path / "ledger.json"))
    old = json.load(open(os.path.join(VAL, "ledger.json")))
    close(new["keys"], old["keys"], 1e-12)
    close(new["summary"], old["summary"], 1e-12)


def test_pbe_vs_exp_reproduces(tmp_path):
    run("pbe_vs_exp.py", tmp_path)
    new, old = json.load(open(tmp_path / "pbe_vs_exp.json")), committed("pbe_vs_exp.json")
    close(new["rows"], old["rows"], 1e-9)
    close(new["core_range"], old["core_range"], 1e-9)
    assert abs(new["dtdp_330"] - old["dtdp_330"]) < 0.01      # 8.41 default when pn_correct.json is absent


@pytest.mark.skipif(not present(paths.PARITY), reason="parity_tel.npz not present (LFS)")
def test_stats_reproduce(tmp_path):
    run("stats_test.py", tmp_path)
    close(json.load(open(tmp_path / "test_stats.json")), committed("test_stats.json"), 1e-9)


@pytest.mark.skipif(not all(present(f) for f in glob.glob(os.path.join(paths.PHASEE, "runs_tel", "*_ckpt.npz")) or [""]),
                    reason="coexistence checkpoints not present (LFS)")
def test_cc_check_reproduces(tmp_path):
    run("cc_check.py", tmp_path)
    close(json.load(open(tmp_path / "cc_check.json")), committed("cc_check.json"), 1e-6)


def _mace_available():
    try:
        import torch, mace  # noqa: F401
        return present(paths.MODEL)
    except ImportError:
        return False


@pytest.mark.skipif(not _mace_available(), reason="torch, mace-torch or the model file not present")
def test_static_compare_reproduces(tmp_path):
    run("val_mace.py", tmp_path)
    new, old = json.load(open(tmp_path / "static_compare.json")), committed("static_compare.json")
    # CPU and GPU float32 kernels differ at the 1e-4 relative level
    close(new, old, 2e-4)
