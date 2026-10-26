"""parse_train.py - per-frame metadata of the training, test and spot-check frames -> inputs/train_meta.json.

The training frames come from phaseD/train_e3.xyz (the production training file, every frame
stamped with its elec_temp), split by T_el into train_Tel4000/6000/8000; this replaces the
per-T_el files train_Tel*.xyz the script first read (same frames, not kept in the repository).
"""
import json, os, sys
import numpy as np
from ase.io import iread

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import layout  # noqa: E402


def row(at):
    i = at.info
    F = at.arrays.get("REF_forces")
    fn = np.linalg.norm(F, axis=1) if F is not None else None
    return dict(phase=i.get("phase"), T=float(i.get("T", np.nan)), P=float(i.get("P_est", np.nan)),
                vat=float(i.get("vat", at.get_volume() / len(at))), N=len(at), id=i.get("frame_id", i.get("spot_name")),
                fmean=float(fn.mean()) if fn is not None else None, fmax=float(fn.max()) if fn is not None else None,
                E=float(i["REF_energy"]) / len(at) if "REF_energy" in i else None,
                fhist=np.histogram(fn, bins=np.linspace(0, 12, 61))[0].tolist() if fn is not None else None)


out = {f"train_Tel{T}": [] for T in (6000, 4000, 8000)}
for at in iread(layout.TRAIN_SET, index=":"):
    out[f"train_Tel{int(at.info['elec_temp'])}"].append(row(at))
out["test_Tel6000"] = [row(at) for at in iread(layout.TEST_SETS[6000], index=":")]
out["spot_frames"] = [row(at) for at in iread(layout.SPOT_FRAMES, index=":")]
out = {k: out[k] for k in ("train_Tel6000", "test_Tel6000", "train_Tel4000", "train_Tel8000", "spot_frames")}
for tag, rows in out.items():
    ph = {}
    for r in rows:
        ph[r["phase"]] = ph.get(r["phase"], 0) + 1
    P = np.array([r["P"] for r in rows]); T = np.array([r["T"] for r in rows]); N = np.array([r["N"] for r in rows])
    print(tag, len(rows), ph, "P", np.nanmin(P).round(1), np.nanmax(P).round(1), "T", np.nanmin(T).round(0), np.nanmax(T).round(0), "N", sorted(set(N.tolist())))
json.dump(out, open(os.path.join(layout.FIG_INPUTS, "train_meta.json"), "w"))
