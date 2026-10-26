#!/usr/bin/env python3
"""
phaseD_dataset.py — build MACE training sets from completed QE labels (fe_melt_v1)

Reads:  labels/labels.jsonl   (run labels/collect.py first)
        selected/frames_*.extxyz
Writes: phaseD/train_Tel{4000,6000,8000}.xyz   (95% of matched frames)
        phaseD/test_Tel{...}.xyz               (5% held out — used ONLY by the
                                                gate suite, never by training)
        phaseD/dataset_report.txt

Energy written is the Mermin free energy F (info key REF_energy, eV);
forces to arrays key REF_forces (eV/A). No stress (not collected).
Safe to rerun daily as labels accumulate — output is rebuilt from scratch.

Usage (ARC login node):  python phaseD_dataset.py
"""
import glob, json, os, random
from ase.io import read, write

random.seed(11)
os.makedirs("phaseD", exist_ok=True)

labels = {}
with open("labels/labels.jsonl") as f:
    for line in f:
        d = json.loads(line)
        labels[(d["id"], d["Tel"])] = d

frames = {}
for fx in sorted(glob.glob("selected/frames_*.extxyz")):
    for fr in read(fx, index=":"):
        frames[fr.info["frame_id"]] = fr

report = []
for tel in (4000, 6000, 8000):
    matched = []
    for (fid, t), d in labels.items():
        if t != tel or fid not in frames:
            continue
        fr = frames[fid].copy()
        if len(d["forces_eVA"]) != len(fr):
            report.append(f"SKIP {fid} Tel{tel}: force count mismatch")
            continue
        fr.info = {k: fr.info.get(k) for k in
                   ("frame_id", "phase", "vat", "T", "P_est") if k in fr.info}
        fr.info["REF_energy"] = d["F_eV"]
        fr.new_array("REF_forces",
                     __import__("numpy").array(d["forces_eVA"]))
        matched.append(fr)
    random.shuffle(matched)
    ntest = max(1, len(matched) // 20) if matched else 0
    test, train = matched[:ntest], matched[ntest:]
    for name, seq in (("train", train), ("test", test)):
        path = f"phaseD/{name}_Tel{tel}.xyz"
        if os.path.exists(path): os.remove(path)
        for fr in seq: write(path, fr, append=True)
    per_phase = {}
    for fr in matched:
        per_phase[fr.info.get("phase", "?")] = per_phase.get(fr.info.get("phase", "?"), 0) + 1
    report.append(f"Tel {tel}: {len(train)} train / {len(test)} test   by phase: {per_phase}")

open("phaseD/dataset_report.txt", "w").write("\n".join(report) + "\n")
print("\n".join(report))
print("\nGuidance: train v0 (Tel6000) once its train count exceeds ~800 and")
print("every phase is represented; retrain on the full set when labels finish.")
