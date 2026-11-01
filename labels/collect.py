#!/usr/bin/env python3
"""collect.py — parse QE label outputs -> labels.jsonl (one line per frame per Tel)."""
import glob, json, re
RY, RYA = 13.6056980659, 25.71104309541616
PAT = re.compile(r"^(.+)_Tel(4000|6000|8000)\.out$")
n_ok = n_inc = 0
with open("labels.jsonl", "w") as out:
    for f in sorted(glob.glob("*.out")):
        m = PAT.match(f)
        if not m:
            continue                      # slurm logs, mixtest, crosscheck, etc.
        t = open(f, errors="ignore").read()
        if "JOB DONE" not in t:
            n_inc += 1
            continue
        F = float(re.findall(r"^!.*?=\s*(-?[\d.]+)\s*Ry", t, re.M)[-1]) * RY
        forces = [[float(x) * RYA for x in g]
                  for g in re.findall(r"atom\s+\d+ type\s+\d+\s+force =\s*"
                                      r"(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", t)]
        P = float(re.findall(r"P=\s*(-?[\d.]+)", t)[-1]) / 10.0
        out.write(json.dumps(dict(id=m.group(1), Tel=int(m.group(2)), F_eV=F,
                                  forces_eVA=forces, P_GPa=P)) + "\n")
        n_ok += 1
print(f"parsed {n_ok} labels; {n_inc} label runs incomplete/failed")
