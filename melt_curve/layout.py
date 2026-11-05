"""layout.py - every location the analysis, validation and figure code reads or writes.

Paths are relative to the repository root, one level above this file. Each path can be overridden
with the environment variable named in the comment; FE_ROOT moves the root.

The ARC pipeline in production_campaign/ is kept verbatim as provenance and does not import this
module. Its relative paths are those of ~/fe_melt on ARC; production_campaign/README.md maps them to
this repository.
"""
import os

ROOT = os.environ.get("FE_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _p(env, rel):
    return os.environ.get(env, os.path.join(ROOT, rel))


# ARC working directory, verbatim (production_campaign/fe_melt)
ARC = _p("FE_ARC", "production_campaign/fe_melt")
LABELS = os.path.join(ARC, "labels")
Q6_MODULE_DIR = _p("FE_Q6_DIR", "production_campaign/fe_melt")                # phaseE_q6.py

# models (provenance in models/README.md)
MODEL = _p("FE_MODEL", "models/fe_tel_e3_stagetwo.model")                     # production potential
MODEL_STAGE_ONE = _p("FE_MODEL_STAGE_ONE", "models/fe_tel_e3.model")          # e3 stage one (identity test only)
MODEL_V1 = {T: _p(f"FE_MODEL_V1_{T}", f"models/fe_ft_Tel{T}_f32.model") for T in (4000, 6000, 8000)}

# training and test data (data/training)
TRAIN_SET = _p("FE_TRAIN_SET", "data/training/train_e3.xyz")
TEST_SETS = {T: _p(f"FE_TEST_SET_{T}", f"data/training/test_Tel{T}.xyz") for T in (4000, 6000, 8000)}
SPOT_FRAMES = _p("FE_SPOT_FRAMES", "data/training/spot_frames.xyz")
DATASET_REPORT = _p("FE_DATASET_REPORT", "data/training/dataset_report.txt")
VALID_INDICES = _p("FE_VALID_INDICES", "production_campaign/fe_melt/fe_tel_e3_valid_indices_7.txt")
PARITY = _p("FE_PARITY", "data/training/parity_tel.npz")                      # parity_tel.py output

# training records (production_campaign/arc_scripts/logs)
LOGS = _p("FE_LOGS", "production_campaign/arc_scripts/logs")
TRAIN_LOG = _p("FE_TRAIN_LOG", os.path.join(LOGS, "fe_tel_e3_run-7.log"))
TRAIN_DEBUG_LOG = _p("FE_TRAIN_DEBUG_LOG", os.path.join(LOGS, "fe_tel_e3_run-7_debug.log"))
TEL_CHECK_LOG = _p("FE_TEL_CHECK_LOG", os.path.join(LOGS, "e3_test_final.out"))  # same-frame T_el check

# coexistence (data/coexistence; phaseE/ on ARC)
PHASEE = _p("FE_PHASEE", "data/coexistence")
CELLS = os.path.join(PHASEE, "cells")
RUNS_TEL = os.path.join(PHASEE, "runs_tel")
RUNS_T6000 = os.path.join(PHASEE, "runs")
ANALYSIS_TEL = os.path.join(PHASEE, "analysis_live_all_runs_tel.json")
ANALYSIS_T6000 = os.path.join(PHASEE, "analysis_live_all.json")               # sets each runs_tel cell's T_el
PN_CORRECT = os.path.join(PHASEE, "pn_correct.json")                           # Simon fit (runs_tel.fit)
DECOMP_TEL = os.path.join(PHASEE, "decomp_tel_tel.json")

# 330 GPa snapshot test: potential vs PBE, LDA - PBE (ldashift.py; studies/lda_shift)
LDASHIFT = _p("FE_LDASHIFT", "studies/lda_shift")                              # ldashift_result.json, md_meta.json, frames.xyz

# DFT engine and pseudopotentials
PSEUDO = _p("FE_PSEUDO", "production_campaign/fe_melt/pseudo")
QE_PW = os.environ.get("FE_QE_PW", "pw.x")

# validation (inputs committed there, outputs, ledger)
VALIDATION = os.path.join(ROOT, "melt_curve", "validation")
DFT_REFS = os.environ.get("FE_VAL_DFT", os.path.join(VALIDATION, "refs", "dft"))
EOS_COLD = os.environ.get("FE_VAL_EOS", os.path.join(VALIDATION, "refs", "eos_pbe_cold"))
LITERATURE = os.environ.get("FE_LITERATURE", os.path.join(VALIDATION, "refs", "literature.json"))
METHOD_SETTINGS = os.environ.get("FE_METHOD_SETTINGS", os.path.join(VALIDATION, "refs", "method_settings.json"))
VAL_OUT = os.environ.get("FE_VAL_OUT", os.path.join(VALIDATION, "out"))
LEDGER = os.environ.get("FE_VAL_LEDGER", os.path.join(VALIDATION, "ledger.json"))
REPORT_NUMBERS = os.path.join(VALIDATION, "report_numbers.tsv")

# figures
FIG_SRC = os.path.join(ROOT, "melt_curve", "figures")
FIG_INPUTS = os.environ.get("FE_FIG_INPUTS", os.path.join(FIG_SRC, "inputs"))  # intermediate json
FIGURES = os.environ.get("FE_FIG_DIR", os.path.join(ROOT, "docs", "figures"))  # PDF + 200-dpi PNG
