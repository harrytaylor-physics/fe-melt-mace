"""Tests for validation/: committed outputs reproduce from their producers, and the ledger
matches the committed outputs. Producers whose inputs are not present skip."""
import copy, csv, glob, json, math, os, re, subprocess, sys
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


def close(a, b, tol, path="", tolf=None):
    """Compare nested outputs: |a - b| <= tol * max(1, |b|), or <= tolf(path, b) if tolf is given."""
    if isinstance(a, dict):
        assert set(a) == set(b), f"{path}: keys {sorted(set(a) ^ set(b))}"
        for k in a:
            close(a[k], b[k], tol, f"{path}.{k}", tolf)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            close(x, y, tol, f"{path}[{i}]", tolf)
    elif isinstance(a, float) or isinstance(b, float):
        t = tolf(path, b) if tolf else tol * max(1.0, abs(b))
        assert abs(a - b) <= t, f"{path}: {a} vs {b} (tol {t:.3g})"
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
    # always on CPU (CUDA hidden): out/static_compare.json was computed on CPU in float32
    run("val_mace.py", tmp_path, CUDA_VISIBLE_DEVICES="")
    new, old = json.load(open(tmp_path / "static_compare.json")), committed("static_compare.json")
    close(new, old, 2e-4)


# GPU float32 kernels against the CPU reference: |a - b| <= max(2e-4 |b|, 2e-3), except where that is
# not below half the last printed digit of a report number taken from (or derived from) the value;
# there the absolute tolerance below applies (test_gpu_tolerance_below_printed_precision checks this).
GPU_REL, GPU_ABS = 2e-4, 2e-3
GPU_TIGHT = {
    ".eos": 2e-3,                 # P to 0.01 GPa at 350 GPa; summary dP and dE maxima, P range
    ".tel_response": 2e-3,        # P to 0.01 GPa at 263-281 GPa, dF to 0.1 meV/atom at 325 meV/atom
    ".fcc-hcp_": 2e-3, ".bcc-hcp_": 2e-3,     # 0.1 meV/atom at 250-355 meV/atom; MACE - DFT differences
    ".fcc_P_": 2e-3, ".bcc_P_": 2e-3,         # 0.1 GPa at 290 GPa; MACE - DFT differences to 0.01 GPa
    ".voigt.dft.vs_km_s": 1e-3, ".voigt.mace.vs_km_s": 1e-3,   # v_s relative error to 0.1 %
    ".coa_min": 2e-5,             # c/a to 1e-4
}


def gpu_tol(path, b):
    for prefix, t in GPU_TIGHT.items():
        if path.startswith(prefix):
            return t
    return max(GPU_REL * abs(b), GPU_ABS)


def _cuda():
    try:
        import torch
        return torch.cuda.is_available()
    except ImportError:
        return False


@pytest.mark.skipif(not (_mace_available() and _cuda()), reason="needs a CUDA GPU, torch, mace-torch and the model file")
def test_static_compare_reproduces_on_gpu(tmp_path):
    env = {k: v for k, v in os.environ.items() if k != "CUDA_VISIBLE_DEVICES" or v}
    e = dict(env, FE_VAL_OUT=str(tmp_path), FE_VAL_LEDGER=str(tmp_path / "ledger.json"))
    r = subprocess.run([sys.executable, "val_mace.py"], cwd=VAL, env=e, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    close(json.load(open(tmp_path / "static_compare.json")), committed("static_compare.json"), None, tolf=gpu_tol)


def _floats(x, p="", where=()):
    """(path as close() names it, key sequence, value) of every float leaf."""
    if isinstance(x, dict):
        for k, v in x.items():
            yield from _floats(v, f"{p}.{k}", where + (k,))
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from _floats(v, f"{p}[{i}]", where + (i,))
    elif isinstance(x, float):
        yield p, where, x


def _set(x, where, v):
    for k in where[:-1]:
        x = x[k]
    x[where[-1]] = v


def _flat_summary(s, pre, out):
    if isinstance(s, dict):
        for k, v in s.items():
            _flat_summary(v, f"{pre}/{k}" if pre else k, out)
    elif isinstance(s, list):
        for i, v in enumerate(s):
            _flat_summary(v, f"{pre}[{i}]", out)
    elif isinstance(s, (int, float)):
        out["summary:" + pre] = s
    return out


def _ledger_values(ledger, J):
    keys = {}
    ledger.flatten(J["static_compare"], "static_compare", keys)
    out = {"keys:" + k: v for k, v in keys.items()}
    for k, v in list(out.items()):
        if isinstance(v, list):
            out.update({f"{k}[{i}]": x for i, x in enumerate(v)})
    return _flat_summary(ledger.summary(J), "", out)


def test_gpu_tolerance_below_printed_precision():
    """Every report number taken from or derived from out/static_compare.json moves by less than half its
    last printed digit when each static_compare value moves by gpu_tol (linear worst case, one value at a
    time through ledger.summary)."""
    import ledger
    J = {stem: committed(stem + ".json") for stem in ledger.PRODUCER if os.path.exists(os.path.join(VAL, "out", stem + ".json"))}
    J.update({stem: json.load(open(f)) for stem, f in ledger.CONSTANTS.items()})
    base = _ledger_values(ledger, J)
    bound = dict.fromkeys(base, 0.0)
    for p, where, b in _floats(J["static_compare"]):
        worst = {}
        for sign in (1, -1):
            Jp = dict(J, static_compare=copy.deepcopy(J["static_compare"]))
            _set(Jp["static_compare"], where, b + sign * gpu_tol(p, b))
            new = _ledger_values(ledger, Jp)
            for k in base:
                if k in new and isinstance(base[k], (int, float)):
                    worst[k] = max(worst.get(k, 0.0), abs(new[k] - base[k]))
        for k, w in worst.items():
            bound[k] += w
    n, bad = 0, []
    with open(os.path.join(VAL, "report_numbers.tsv")) as fh:
        for num, unit, key, value in list(csv.reader(fh, delimiter="\t"))[1:]:
            if not bound.get(key):
                continue
            L = float(value)
            parts = re.findall(r"\d[\d,]*(?:\.\d+)?", num)          # "1-2.5" prints two values: take the matching one
            if not parts:
                continue
            txt = min(parts, key=lambda t: abs(abs(float(t.replace(",", ""))) - abs(L) * (100 if unit == "%" else 1)))
            v = float(txt.replace(",", ""))
            half = 0.5 * 10.0 ** -len(txt.partition(".")[2])
            scale = 1.0 if L == 0 or v == 0 else 10.0 ** round(math.log10(v / abs(L)))   # e.g. a fraction printed in %
            n += 1
            if bound[key] * scale >= half:
                bad.append(f"{num} {unit} {key}: moves by up to {bound[key] * scale:.3g}, half digit {half:.3g}")
    assert n >= 80, n
    assert not bad, "\n".join(bad)
