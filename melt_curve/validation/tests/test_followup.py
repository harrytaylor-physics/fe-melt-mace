"""Producers added 29 Oct 2026: convergence, mixed LDA setup, 330 GPa snapshot test."""
import json, os, subprocess, sys
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)


def run(script, tmp_path):
    env = dict(os.environ, FE_VAL_OUT=str(tmp_path))
    subprocess.run([sys.executable, os.path.join(HERE, script)], cwd=HERE, env=env, check=True,
                   capture_output=True)


def same(a, b):
    assert a.keys() == b.keys()
    for k in a:
        if isinstance(a[k], (int, float)):
            assert a[k] == pytest.approx(b[k], rel=1e-9, abs=1e-12), k


@pytest.mark.parametrize("script,out", [("conv_check.py", "convergence.json"),
                                        ("lda_setup_check.py", "lda_setup.json")])
def test_reproduces_committed_output(script, out, tmp_path):
    run(script, tmp_path)
    new = json.load(open(tmp_path / out))
    old = json.load(open(os.path.join(HERE, "out", out)))
    if out == "convergence.json":
        new.pop("rows"); old.pop("rows")
    same(new, old)


def test_ldashift_summary(tmp_path):
    import paths
    src = getattr(paths, "LDASHIFT", os.path.join(paths.ROOT, "studies", "lda_shift"))
    if not os.path.exists(os.path.join(src, "ldashift_result.json")):
        pytest.skip("studies/lda_shift/ldashift_result.json absent")
    run("ldashift_summary.py", tmp_path)
    same(json.load(open(tmp_path / "ldashift.json")), json.load(open(os.path.join(HERE, "out", "ldashift.json"))))
