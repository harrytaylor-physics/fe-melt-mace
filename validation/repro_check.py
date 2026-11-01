"""repro_check.py - re-evaluate the held-out test frames with the production potential and compare
with the predictions stored on ARC (parity_tel.npz, from parity_tel.py with the tel-mace build).

Same frames and order as parity_tel.py: phaseD/test_Tel{4000,6000,8000}.xyz, each frame at its own
T_el (load_model.py: upstream mace-torch takes elec_temp = T_el / 1e4), float32.
Output: out/repro_check.json - largest |dE| (meV/atom) and |dF| (eV/A, per component) between this
evaluation and the stored one, overall and per T_el; and the float64 - float32 energy difference
(mean, spread). Needs torch, mace-torch and the model file; about 8 min on 4 CPU threads.
"""
import json, os, sys
import numpy as np
from ase.io import read
import paths
import load_model as lm
import mace  # noqa: E402

sys.path.insert(0, os.path.dirname(paths.HERE))
import layout  # noqa: E402


def evaluate(dtype):
    calc = lm.calc(dtype, threads=int(os.environ.get("FE_THREADS", "4")))
    E, F = [], []
    for tel in (4000, 6000, 8000):
        for at in read(layout.TEST_SETS[tel], ":"):
            at.info["elec_temp"] = tel / lm.TEL_SCALE
            calc.reset()
            at.calc = calc
            E.append(at.get_potential_energy() / len(at))
            F.append(at.get_forces())
    return np.array(E), np.concatenate(F)


def main():
    z = np.load(paths.PARITY, allow_pickle=True)
    E, F = evaluate("float32")
    E64, _ = evaluate("float64")
    d64 = 1e3 * (E64 - E)
    dE = 1e3 * np.abs(E - z["E_mod"])
    dF = np.abs(F - z["F_mod"])
    out = dict(n_frames=int(len(E)), n_force_components=int(dF.size), max_abs_dE_meV_atom=float(dE.max()),
               max_abs_dF_eV_A=float(dF.max()), rms_dF_eV_A=float(np.sqrt((dF ** 2).mean())),
               per_tel={str(T): dict(max_abs_dE_meV_atom=float(dE[z["e_tel"] == T].max()),
                                     max_abs_dF_eV_A=float(dF[z["f_tel"] == T].max())) for T in (4000, 6000, 8000)},
               float64_minus_float32_meV_atom=dict(mean=float(d64.mean()), min=float(d64.min()), max=float(d64.max()),
                                                    std=float(d64.std())),
               tel_scale=lm.TEL_SCALE, mace_version=mace.__version__)
    print(json.dumps(out, indent=1))
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "repro_check.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
