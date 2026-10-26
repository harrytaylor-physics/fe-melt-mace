import sys, json, re, numpy as np
from ase.io import iread
out = {}
for tag in ["train_Tel6000", "test_Tel6000", "train_Tel4000", "train_Tel8000", "spot_frames"]:
    rows = []
    for at in iread(f"arc/phaseD/{tag}.xyz", index=":"):
        i = at.info
        F = at.arrays.get("REF_forces")
        fn = np.linalg.norm(F, axis=1) if F is not None else None
        rows.append(dict(phase=i.get("phase"), T=float(i.get("T", np.nan)), P=float(i.get("P_est", np.nan)),
                         vat=float(i.get("vat", at.get_volume() / len(at))), N=len(at), id=i.get("frame_id", i.get("spot_name")),
                         fmean=float(fn.mean()) if fn is not None else None, fmax=float(fn.max()) if fn is not None else None,
                         E=float(i["REF_energy"]) / len(at) if "REF_energy" in i else None,
                         fhist=np.histogram(fn, bins=np.linspace(0, 12, 61))[0].tolist() if fn is not None else None))
    out[tag] = rows
    ph = {}
    for r in rows: ph[r["phase"]] = ph.get(r["phase"], 0) + 1
    P = np.array([r["P"] for r in rows]); T = np.array([r["T"] for r in rows]); N = np.array([r["N"] for r in rows])
    print(tag, len(rows), ph, "P", np.nanmin(P).round(1), np.nanmax(P).round(1), "T", np.nanmin(T).round(0), np.nanmax(T).round(0), "N", sorted(set(N.tolist())))
json.dump(out, open("figs/train_meta.json", "w"))
