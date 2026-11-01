"""ldashift.lda_input() writes the LDA setup that ran on ARC: each committed LDA input of the 330 GPa
snapshot test is its committed PBE input with the job name changed and input_dft = 'PZ' added as the
first line of &SYSTEM, on the same PBE PAW file."""
import glob, os, re, sys
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
PBE = sorted(glob.glob(os.path.join(ROOT, "labels", "ldaX_*_PBE_Tel*.in")))


def test_forty_pairs():
    assert len(PBE) == 40
    for p in PBE:
        assert os.path.exists(p.replace("_PBE_", "_LDA_")), p


@pytest.mark.parametrize("pbe", PBE, ids=os.path.basename)
def test_lda_input_reproduces_committed(pbe):
    ldashift = pytest.importorskip("ldashift")
    stem_p = os.path.basename(pbe)[:-3]
    stem_l = stem_p.replace("_PBE_", "_LDA_")
    lda = open(os.path.join(ROOT, "labels", stem_l + ".in")).read()
    assert ldashift.lda_input(open(pbe).read().replace(stem_p, stem_l)) == lda
    assert re.search(r"&SYSTEM\ninput_dft = 'PZ'\n", lda)
    assert "Fe.pbe-spn-kjpaw_psl.1.0.0.UPF" in lda and "Fe.pz-spn" not in lda
