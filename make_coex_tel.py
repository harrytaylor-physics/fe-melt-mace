#!/usr/bin/env python3
"""
make_coex_tel.py — derive phaseE_coex_tel.py from the validated Tel-6000
coexistence driver on ARC (phaseE_coex.py, which carries the ARC-side edits:
float32, cuEquivariance, neighbour-list skin 0.3, velocities at T0). The
original file is left untouched.

Changes, each required to match exactly once:
  - new arguments --elec-temp (K, required) and --rundir (default phaseE/runs_tel)
  - checkpoints and run JSONs written under --rundir
  - atoms.info["elec_temp"] set when the cell is read (the Tel-conditioned
    calculator reads it on every force call; it survives checkpoint resumes
    because positions/momenta are loaded into the same Atoms object)
  - elec_temp printed at start and stored in the final JSON

Run once in ~/fe_melt:  python3 make_coex_tel.py
"""
src = open("phaseE_coex.py").read()
edits = [
    ('ap.add_argument("--device", default="cuda")',
     'ap.add_argument("--device", default="cuda")\n'
     '    ap.add_argument("--elec-temp", type=float, required=True)\n'
     '    ap.add_argument("--rundir", default="phaseE/runs_tel")'),
    ('os.makedirs("phaseE/runs", exist_ok=True)', 'os.makedirs(a.rundir, exist_ok=True)'),
    ('ckpt = f"phaseE/runs/{tag}_ckpt.npz"', 'ckpt = f"{a.rundir}/{tag}_ckpt.npz"'),
    ('out_json = f"phaseE/runs/{tag}.json"', 'out_json = f"{a.rundir}/{tag}.json"'),
    ('at = read(a.cell + ".extxyz")',
     'at = read(a.cell + ".extxyz")\n    at.info["elec_temp"] = float(a.elec_temp)\n'
     '    print(f"[coex] elec_temp = {a.elec_temp:.0f} K", flush=True)'),
    ('json.dump(dict(tag=meta["tag"], run=tag,', 'json.dump(dict(tag=meta["tag"], run=tag, elec_temp=a.elec_temp,'),
]
for old, new in edits:
    n = src.count(old)
    if n != 1:
        raise SystemExit(f"expected exactly one match, found {n}: {old!r}")
    src = src.replace(old, new)
open("phaseE_coex_tel.py", "w").write(src)
print("written phaseE_coex_tel.py; phaseE_coex.py unchanged")
