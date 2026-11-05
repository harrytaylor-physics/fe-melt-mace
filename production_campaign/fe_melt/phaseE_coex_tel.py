#!/usr/bin/env python3
"""
phaseE_coex.py — two-phase NVE coexistence run, segmented (fe_melt_v1, Phase E)

Method unchanged (Morris-type): solid+liquid cell integrated in NVE; latent
heat drives the system onto the melting line, and the time-averaged (P*, T*)
after the transient is a point ON the melting curve.

New in v2: the run is segmented so it survives short scheduler windows.
State (positions, momenta, cell, accumulated T/P/solid-fraction series) is
checkpointed every few blocks; the process exits cleanly when --wall-hours is
nearly spent and the next invocation resumes from the checkpoint. The sbatch
resubmits itself until the JSON reports status "complete", so a 200 ps run can
be assembled from ~4-hour segments that backfill readily.

Usage:
  python phaseE_coex.py --model <m> --cell phaseE/cells/coex_hcp_300GPa \
      [--tinit-frac 1.0] [--ps 200] [--seed 0] [--wall-hours 3.5] [--device cuda]
Exit code 0 = complete, 64 = partial (resubmit), other = error.
"""
from phaseE_nlcache import enable_nl_cache
enable_nl_cache(skin=0.0)
import argparse, json, os, sys, time
import numpy as np
from ase import units
from ase.io import read
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as MBD
from ase.units import GPa
from phaseE_build import q6_per_atom

T0_WALL = time.time()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--cell", required=True, help="path prefix without extension")
    ap.add_argument("--tinit-frac", type=float, default=1.0)
    ap.add_argument("--ps", type=float, default=200.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--wall-hours", type=float, default=3.5,
                    help="exit cleanly this many hours after start")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--elec-temp", type=float, required=True)
    ap.add_argument("--rundir", default="phaseE/runs_tel")
    a = ap.parse_args()

    meta = json.load(open(a.cell + ".json"))
    tag = f"{meta['tag']}_ti{a.tinit_frac:.2f}_s{a.seed}"
    os.makedirs(a.rundir, exist_ok=True)
    ckpt = f"{a.rundir}/{tag}_ckpt.npz"
    out_json = f"{a.rundir}/{tag}.json"

    at = read(a.cell + ".extxyz")
    at.info["elec_temp"] = float(a.elec_temp)
    print(f"[coex] elec_temp = {a.elec_temp:.0f} K", flush=True)
    from mace.calculators import MACECalculator
    at.calc = MACECalculator(enable_cueq=True, model_paths=os.path.expanduser(a.model),
                             device=a.device, default_dtype="float32")

    if os.path.exists(ckpt):
        z = np.load(ckpt, allow_pickle=True)
        at.set_positions(z["positions"]); at.set_momenta(z["momenta"])
        at.set_cell(z["cell"]); at.pbc = True
        T_t = list(z["T_t"]); P_t = list(z["P_t"])
        fr_t = [tuple(x) for x in z["fr_t"]]
        blk0 = int(z["blk"])
        print(f"[coex] resumed {tag} at {blk0} ps", flush=True)
    else:
        T0 = a.tinit_frac * meta["Tm_est"]
        MBD(at, temperature_K=T0, rng=np.random.RandomState(a.seed))
        T_t, P_t, fr_t, blk0 = [], [], [], 0
        print(f"[coex] started {tag}, Tinit={T0:.0f} K", flush=True)

    dyn = VelocityVerlet(at, 2.0 * units.fs)
    nblocks = int(a.ps)
    thr = meta["q6_threshold"]
    verdict = None

    def save(blk):
        np.savez(ckpt, positions=at.get_positions(), momenta=at.get_momenta(),
                 cell=np.array(at.cell), T_t=np.array(T_t), P_t=np.array(P_t),
                 fr_t=np.array(fr_t, dtype=float) if fr_t else np.zeros((0, 2)),
                 blk=blk)

    for b in range(blk0, nblocks):
        dyn.run(500)                                   # 1 ps
        T_t.append(at.get_temperature())
        P_t.append(-np.trace(at.get_stress(voigt=False)) / 3 / GPa)
        if b % 5 == 0:
            fr = float((q6_per_atom(at) > thr).mean())
            fr_t.append((b, fr))
            print(f"  {b:4d} ps  T={T_t[-1]:6.0f} K  P={P_t[-1]:6.1f} GPa  "
                  f"solid_frac={fr:.2f}", flush=True)
            if len(fr_t) >= 5 and all(f < 0.05 for _, f in fr_t[-5:]):
                verdict = "MELTED"; save(b + 1); break
            if len(fr_t) >= 5 and all(f > 0.95 for _, f in fr_t[-5:]):
                verdict = "FROZE"; save(b + 1); break
        if b % 10 == 9:
            save(b + 1)
        if (time.time() - T0_WALL) / 3600.0 > a.wall_hours:
            save(b + 1)
            json.dump(dict(tag=tag, status="partial", ps_done=b + 1,
                           ps_target=nblocks), open(out_json, "w"), indent=1)
            print(f"[coex] wall reached at {b+1}/{nblocks} ps — checkpointed", flush=True)
            sys.exit(64)
    else:
        fr_late = [f for b_, f in fr_t if b_ > nblocks // 2]
        verdict = ("COEXISTING" if fr_late and 0.25 < np.mean(fr_late) < 0.75
                   and (max(fr_late) - min(fr_late)) < 0.25 else "DRIFTING")
        save(nblocks)

    half = len(T_t) // 2
    Tm = float(np.mean(T_t[half:])); Pm = float(np.mean(P_t[half:]))
    nb = 8
    Terr = float(np.std([x.mean() for x in np.array_split(np.array(T_t[half:]), nb)]) / np.sqrt(nb))
    Perr = float(np.std([x.mean() for x in np.array_split(np.array(P_t[half:]), nb)]) / np.sqrt(nb))
    json.dump(dict(tag=meta["tag"], run=tag, elec_temp=a.elec_temp, status="complete", seed=a.seed,
                   tinit_frac=a.tinit_frac, verdict=verdict, T=Tm, T_err=Terr,
                   P=Pm, P_err=Perr, solid_frac=fr_t, ps=len(T_t)),
              open(out_json, "w"), indent=1)
    print(f"[coex] {verdict}: T={Tm:.0f}({Terr:.0f}) K  P={Pm:.1f}({Perr:.1f}) GPa", flush=True)

if __name__ == "__main__":
    main()
