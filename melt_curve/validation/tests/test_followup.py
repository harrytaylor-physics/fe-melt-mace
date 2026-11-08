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
    out = json.load(open(tmp_path / "ldashift.json"))
    assert out["mace_Tm_total_K"] == pytest.approx(19.42, abs=0.01) and out["mace_Tm_se_K"] == pytest.approx(4.80, abs=0.01)
    assert out["lda_dTm_K"] == pytest.approx(11.28, abs=0.01) and out["lda_dTm_se_K"] == pytest.approx(14.49, abs=0.01)


def test_melting_entropy(tmp_path):
    import melting_entropy as me
    if not os.path.exists(me.INPUT):
        pytest.skip("studies/lda_shift/pbe_snapshot_energies_330.json absent")
    run("melting_entropy.py", tmp_path)
    out = json.load(open(tmp_path / "melting_entropy.json"))
    same(out, json.load(open(os.path.join(HERE, "out", "melting_entropy.json"))))
    d = json.load(open(me.INPUT))
    x, liq = d["expected_results"], d["md_meta"]["liquid"]
    V_liq = me.at_target(liq["V_at"], out["P_total_liq_GPa"], liq["B_GPa"], d["P_target_GPa"])
    for got, want, digits in [(out["dE_meV"], x["dE_meV"], 0), (out["dH_meV"], x["dH_uncorrected_meV"], 0),
                              (out["dS_kB"], x["dS_kB"], 3), (out["dS_err_kB"], x["dS_err_kB"], 3),
                              (out["dV_over_V_pct"], x["dV_over_V_uncorrected_pct"], 2),
                              (out["P_total_hcp_GPa"], x["P_total_GPa"]["hcp"], 1),
                              (out["P_total_liq_GPa"], x["P_total_GPa"]["liquid"], 2),
                              (V_liq, x["liquid_V_at_330GPa"], 4),
                              (out["dV_over_V_330_pct"], x["dV_over_V_corrected_pct"], 2)]:
        assert round(got, digits) == pytest.approx(want), (got, want)
