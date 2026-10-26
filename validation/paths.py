"""paths.py - every file the validation producers read or write.

Defaults are relative to the repository root, laid out as ~/fe_melt on ARC.
Any path can be overridden with the environment variable of the same name;
the repository's own layout module can replace this file.
"""
import os

ROOT = os.environ.get("FE_ROOT", os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


def _p(env, rel):
    return os.environ.get(env, os.path.join(ROOT, rel))


# inputs produced elsewhere in the pipeline
MODEL = _p("FE_MODEL", "phaseD/models/fe_tel_e3_stagetwo.model")      # production T_el-conditioned potential
PARITY = _p("FE_PARITY", "parity_tel.npz")                            # held-out test-set predictions (parity_tel.py)
PHASEE = _p("FE_PHASEE", "phaseE")                                    # analysis_live_all_runs_tel.json, cells/, runs_tel/
Q6_MODULE_DIR = _p("FE_Q6_DIR", ".")                                  # directory holding phaseE_q6.py

# inputs committed with this directory
HERE = os.path.dirname(os.path.abspath(__file__))
DFT_REFS = os.environ.get("FE_VAL_DFT", os.path.join(HERE, "refs", "dft"))           # val_dft.py output
EOS_COLD = os.environ.get("FE_VAL_EOS", os.path.join(HERE, "refs", "eos_pbe_cold"))  # static PBE, 300 K comparison

# outputs
OUT = os.environ.get("FE_VAL_OUT", os.path.join(HERE, "out"))
LEDGER = os.environ.get("FE_VAL_LEDGER", os.path.join(HERE, "ledger.json"))

# DFT engine for val_dft.py
QE_PW = os.environ.get("FE_QE_PW", "pw.x")
PSEUDO_DIR = _p("FE_PSEUDO", "pseudo")
