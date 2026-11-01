"""Which model did the T_el-conditioned coexistence runs (phaseE/runs_tel) load?

No run record names the model, and no log holds a potential energy of a saved configuration.
Each checkpoint does hold one model-dependent number for its final configuration: the last entry
of P_t, the virial pressure -tr(sigma)/3 that phaseE_coex_tel.py computed with the loaded model
for exactly the positions and cell it then saved (P_t is appended after the 1 ps block, the
checkpoint is written before the next step).

For every runs_tel checkpoint this test re-evaluates that pressure with fe_tel_e3_stagetwo.model
(production, SWA stage two) and with fe_tel_e3.model (stage one) at the run's T_el and checks that
stage two reproduces the logged value and stage one does not. T_el is taken from the run json and
checked against the rule in fe_phaseE_tel.sbatch (round of the cell's T_m in analysis_live_all.json).

GPU only (12,096-atom cells; on CPU the evaluation exhausts a 24 GB laptop). Results are written to
tests/out/runs_tel_model.json. With the tel-mace build (scales T_el internally) set FE_TEL_SCALE=1.

Memory: the pressure needs the stress, so each evaluation runs the autograd pass (energies-only
evaluation under torch.inference_mode would drop the quantity tested). As in phaseE_coex_tel.py the
calculator is float32 with the cuEquivariance kernels (enable_cueq=True; e3nn if cuequivariance is not
installed, recorded as "cueq" in the output), the allocator uses expandable segments (as
phaseE_patch_mem.py advises), and only one model is on the GPU at a time: all checkpoints with stage
two, then the model is freed, then all with stage one. With e3nn kernels and both models loaded the
evaluation ran out of memory on an 80 GB H100 (ARC, 1 Nov 2026).
"""
import gc, glob, json, os, sys

os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")    # read at the first CUDA allocation
import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import layout  # noqa: E402

TOL_GPA = 0.02          # GPU float32 / kernel differences at 12,096 atoms are far below this


def _lfs(path):
    if not os.path.exists(path):
        return False
    with open(path, "rb") as fh:
        return fh.read(7) != b"version"


def _json(path):
    with open(path) as fh:
        return json.load(fh)


def _gpu():
    try:
        import torch
        import mace  # noqa: F401
        return torch.cuda.is_available()
    except ImportError:
        return False


CKPTS = sorted(glob.glob(os.path.join(layout.RUNS_TEL, "*_ckpt.npz")))
need = pytest.mark.skipif(
    not (_gpu() and CKPTS and all(_lfs(f) for f in CKPTS[:1]) and _lfs(layout.MODEL) and _lfs(layout.MODEL_STAGE_ONE)),
    reason="needs a CUDA GPU, torch, mace-torch and the LFS files (runs_tel checkpoints, both e3 models)")


def run_elec_temp(tag):
    """T_el of a runs_tel run: fe_phaseE_tel.sbatch sets round(T_m) of the cell in analysis_live_all.json.
    The run json records it too, except for runs left 'partial' (stopped before 200 ps)."""
    cell = tag.split("_ti")[0]
    rule = round(_json(layout.ANALYSIS_T6000)["cells"][cell]["Tm"])
    rec = _json(os.path.join(layout.RUNS_TEL, tag + ".json")).get("elec_temp")
    return rule, rec


def test_runs_tel_elec_temp_follows_sbatch_rule():
    runs = sorted(glob.glob(os.path.join(layout.RUNS_TEL, "coex_*_s[0-9].json")))
    assert len(runs) == 48
    recorded = 0
    for f in runs:
        rule, rec = run_elec_temp(os.path.basename(f)[:-5])
        if rec is not None:
            assert rec == rule, f
            recorded += 1
    assert recorded >= 40


@need
def test_runs_tel_loaded_e3_stagetwo():
    import torch
    from ase import Atoms
    from ase.units import GPa
    from mace.calculators import MACECalculator
    scale = float(os.environ.get("FE_TEL_SCALE", "1e4"))
    keys = {"elec_temp": "elec_temp", "total_spin": "spin", "total_charge": "charge"}

    def calculator(path):
        try:
            return MACECalculator(model_paths=path, device="cuda", default_dtype="float32", info_keys=keys,
                                  enable_cueq=True), True
        except ImportError:                     # cuequivariance not installed: e3nn kernels
            return MACECalculator(model_paths=path, device="cuda", default_dtype="float32", info_keys=keys), False

    rows = {}
    for name, path in (("stagetwo", layout.MODEL), ("stageone", layout.MODEL_STAGE_ONE)):
        calc, cueq = calculator(path)
        for f in CKPTS:
            tag = os.path.basename(f)[:-len("_ckpt.npz")]
            with np.load(f, allow_pickle=True) as z:
                pos, cell, P_logged, blk = z["positions"], z["cell"], float(z["P_t"][-1]), int(z["blk"])
            tel = run_elec_temp(tag)[0]
            at = Atoms("Fe%d" % len(pos), positions=pos, cell=cell, pbc=True)
            at.info["elec_temp"] = tel / scale
            row = rows.setdefault(tag, dict(run=tag, elec_temp=tel, P_logged=P_logged, ps=blk))
            calc.reset()
            at.calc = calc
            row[f"P_{name}"] = float(-np.trace(at.get_stress(voigt=False)) / 3 / GPa)
            row[f"dP_{name}"] = row[f"P_{name}"] - row["P_logged"]
            at.calc = None
            torch.cuda.empty_cache()
        del calc
        gc.collect()
        torch.cuda.empty_cache()
    rows = list(rows.values())
    os.makedirs(os.path.join(ROOT, "tests", "out"), exist_ok=True)
    s = dict(n=len(rows), max_abs_dP_stagetwo=max(abs(r["dP_stagetwo"]) for r in rows),
             min_abs_dP_stageone=min(abs(r["dP_stageone"]) for r in rows), tel_scale=scale, cueq=cueq,
             peak_gpu_GiB=torch.cuda.max_memory_allocated() / 2**30, rows=rows)
    with open(os.path.join(ROOT, "tests", "out", "runs_tel_model.json"), "w") as fh:
        json.dump(s, fh, indent=1)
    print(f"\n{s['n']} runs: max |P(stage two) - logged| {s['max_abs_dP_stagetwo']:.4f} GPa, "
          f"min |P(stage one) - logged| {s['min_abs_dP_stageone']:.4f} GPa")
    assert len(rows) == 48
    assert s["max_abs_dP_stagetwo"] < TOL_GPA, "stage two does not reproduce the logged pressures"
    assert s["min_abs_dP_stageone"] > TOL_GPA, "stage one is not distinguishable from stage two on some run"
