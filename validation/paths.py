"""paths.py - the names the validation producers use, taken from the repository's layout module.

Every path is defined once, in ../layout.py, and can be overridden there by the environment
variable of the same name (FE_MODEL, FE_PARITY, FE_PHASEE, FE_VAL_OUT, FE_VAL_LEDGER, ...).
"""
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import layout  # noqa: E402

ROOT = layout.ROOT

# inputs produced elsewhere in the pipeline
MODEL = layout.MODEL                  # production T_el-conditioned potential
MODEL_STAGE_ONE = layout.MODEL_STAGE_ONE
PARITY = layout.PARITY                # held-out test-set predictions (parity_tel.py)
PHASEE = layout.PHASEE                # analysis_live_all_runs_tel.json, cells/, runs_tel/
Q6_MODULE_DIR = layout.Q6_MODULE_DIR  # directory holding phaseE_q6.py

# inputs committed with this directory
DFT_REFS = layout.DFT_REFS            # val_dft.py output
EOS_COLD = layout.EOS_COLD            # static PBE, 300 K comparison
LITERATURE = layout.LITERATURE
METHOD_SETTINGS = layout.METHOD_SETTINGS

# outputs
OUT = layout.VAL_OUT
LEDGER = layout.LEDGER

# DFT engine for val_dft.py
QE_PW = layout.QE_PW
PSEUDO_DIR = layout.PSEUDO
