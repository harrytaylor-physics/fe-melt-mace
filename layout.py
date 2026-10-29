"""layout.py - every location the analysis, validation and figure code reads or writes.

Paths are relative to the repository root, which keeps the layout of ~/fe_melt on ARC (the
ARC scripts at the root run from there and use these same relative paths). Each path can be
overridden with the environment variable named in the comment; FE_ROOT moves the root.

The ARC pipeline scripts at the root (phase*.py, *.sbatch) are kept verbatim as provenance and
do not import this module.
"""
import os

ROOT = os.environ.get("FE_ROOT", os.path.dirname(os.path.abspath(__file__)))


def _p(env, rel):
    return os.environ.get(env, os.path.join(ROOT, rel))


# models (phaseD/models; provenance in models/README.md)
MODEL = _p("FE_MODEL", "phaseD/models/fe_tel_e3_stagetwo.model")              # production potential
MODEL_STAGE_ONE = _p("FE_MODEL_STAGE_ONE", "phaseD/models/fe_tel_e3.model")   # e3 stage one (identity test only)
MODEL_V1 = {T: _p(f"FE_MODEL_V1_{T}", f"phaseD/models/fe_ft_Tel{T}_f32.model") for T in (4000, 6000, 8000)}

# training and test data (phaseD)
TRAIN_SET = _p("FE_TRAIN_SET", "phaseD/train_e3.xyz")
TEST_SETS = {T: _p(f"FE_TEST_SET_{T}", f"phaseD/test_Tel{T}.xyz") for T in (4000, 6000, 8000)}
SPOT_FRAMES = _p("FE_SPOT_FRAMES", "phaseD/spot_frames.xyz")
DATASET_REPORT = _p("FE_DATASET_REPORT", "phaseD/dataset_report.txt")
VALID_INDICES = _p("FE_VALID_INDICES", "fe_tel_e3_valid_indices_7.txt")
PARITY = _p("FE_PARITY", "phaseD/parity_tel.npz")                             # parity_tel.py output

# training records (logs, results)
TRAIN_LOG = _p("FE_TRAIN_LOG", "logs/fe_tel_e3_run-7.log")
TRAIN_DEBUG_LOG = _p("FE_TRAIN_DEBUG_LOG", "logs/fe_tel_e3_run-7_debug.log")
TEL_CHECK_LOG = _p("FE_TEL_CHECK_LOG", "logs/e3_test_final.out")             # same-frame T_el check
LOGS = _p("FE_LOGS", "logs")

# coexistence (phaseE)
PHASEE = _p("FE_PHASEE", "phaseE")
CELLS = os.path.join(PHASEE, "cells")
RUNS_TEL = os.path.join(PHASEE, "runs_tel")
RUNS_T6000 = os.path.join(PHASEE, "runs")
ANALYSIS_TEL = os.path.join(PHASEE, "analysis_live_all_runs_tel.json")
ANALYSIS_T6000 = os.path.join(PHASEE, "analysis_live_all.json")               # sets each runs_tel cell's T_el
PN_CORRECT = os.path.join(PHASEE, "pn_correct.json")                           # Simon fit (runs_tel.fit)
DECOMP_TEL = os.path.join(PHASEE, "decomp_tel_tel.json")
Q6_MODULE_DIR = _p("FE_Q6_DIR", ".")                                          # phaseE_q6.py

# 330 GPa snapshot test: potential vs PBE, LDA - PBE (ldashift.py; phaseX)
LDASHIFT = _p("FE_LDASHIFT", "phaseX/ldashift")                              # ldashift_result.json, md_meta.json, frames.xyz

# DFT engine and pseudopotentials
PSEUDO = _p("FE_PSEUDO", "pseudo")
QE_PW = os.environ.get("FE_QE_PW", "pw.x")

# validation (inputs committed there, outputs, ledger)
VALIDATION = os.path.join(ROOT, "validation")
DFT_REFS = os.environ.get("FE_VAL_DFT", os.path.join(VALIDATION, "refs", "dft"))
EOS_COLD = os.environ.get("FE_VAL_EOS", os.path.join(VALIDATION, "refs", "eos_pbe_cold"))
LITERATURE = os.environ.get("FE_LITERATURE", os.path.join(VALIDATION, "refs", "literature.json"))
METHOD_SETTINGS = os.environ.get("FE_METHOD_SETTINGS", os.path.join(VALIDATION, "refs", "method_settings.json"))
VAL_OUT = os.environ.get("FE_VAL_OUT", os.path.join(VALIDATION, "out"))
LEDGER = os.environ.get("FE_VAL_LEDGER", os.path.join(VALIDATION, "ledger.json"))
REPORT_NUMBERS = os.path.join(VALIDATION, "report_numbers.tsv")

# figures
FIG_SRC = os.path.join(ROOT, "analysis", "figures")
FIG_INPUTS = os.environ.get("FE_FIG_INPUTS", os.path.join(FIG_SRC, "inputs"))  # intermediate json
FIGURES = os.environ.get("FE_FIG_DIR", os.path.join(ROOT, "docs", "figures"))  # PDF + 200-dpi PNG
