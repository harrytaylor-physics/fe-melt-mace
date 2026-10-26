#!/usr/bin/env python3
"""
ldashift.py — shift of the melting temperature at 330 GPa from replacing PBE by
LDA, by free-energy perturbation on PBE-potential snapshots (fe_melt, step X1).

Why
  Rao & Zhu (LDA+DMFT) find T_m(330 GPa) = 6,225 +- 42 K, 104 K below this work
  (PBE, 6,329 +- 26 K). Their DFT part is LDA and they report no LDA-only
  baseline, so the 104 K mixes the LDA - PBE difference with the correlation
  (DMFT) correction. This calculation measures the LDA - PBE part directly.

Physics
  Snapshots of bulk hcp and bulk liquid at the PBE melting point
  (P = 330 GPa, T = 6,329 K, T_el = T) are taken from the Tel MACE potential
  in NVT at each phase's NPT mean volume V_p. On each snapshot both functionals
  give the Mermin free energy F at T_el (same PAW construction, pslibrary 1.0.0
  'spn', same cutoffs, k-points and smearing; only the functional differs).
  Per phase p, with dU = F_LDA - F_PBE per cell and N atoms:
    first order      a1 = < dU > / N
    configurational  a2 = -(beta / 2) var(dU) / N
    volume term      a3 = -(V_p / (2 B_p)) dP^2,  dP = < P_LDA - P_PBE >
  a3 converts the fixed-volume result to fixed pressure: at V_p the LDA
  pressure differs from 330 GPa by dP, and G_LDA(P) = F_LDA(V_p) + P V_p
  - V_p dP^2 / (2 B) to second order. B_p is the isothermal bulk modulus of
  phase p from the potential (NVT at V_p x (1 +- 0.01)).
    dG_p = a1 + a2 + a3,   dT_m = +(dG_liq - dG_hcp) / dS_m,
  (a correction that raises the liquid's free energy relative to the solid
  raises T_m; from G_liq - G_hcp = -dS_m (T - T_m) near T_m),
  with dS_m = 0.80 k_B (this work's 330 GPa cell).

Subcommands
  md      (GPU)   2 phases x 2 runs x 10 snapshots (1 ps apart) = 40
                  snapshots of 128 atoms, plus B_p; then writes the inputs.
  inputs  (login) pw.x inputs for each snapshot: PBE (campaign template) and
                  LDA (same template, LDA pseudopotential, any input_dft line
                  removed). labels/ldaX_*, labels/ldaX_list.txt,
                  labels/ldaX_array.sbatch. 80 calculations, ~3,600 core-hours.
  eval    (login) checks every output used the intended functional, then
                  computes a1, a2, a3, dT_m and its standard error; also the
                  potential vs PBE at this state point: pressure error and force
                  RMSE per phase, and the liquid - hcp energy bias as a T_m error.

Run from ~/fe_melt (tel_env for md):
  sbatch ldashift_md.sbatch
  cd labels && sbatch ldaX_array.sbatch && cd ..
  python ldashift.py eval
"""
import argparse, glob, json, os, re, sys
import numpy as np
from ase import Atoms, units
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as MBD

from phaseF_common import (hcp_ortho, frac_msd, q6_reference, solid_fraction, load_calc,
                           template_for, make_input, write_array_sbatch, parse_out)

MODEL = "phaseD/models/fe_tel_e3_stagetwo.model"
P_GPA, T_K = 330.0, 6329.0
V_GUESS = 7.05                       # A^3/atom, cell mean near 330 GPa on the melting line
RUNS, SNAPS = 2, 10
DS_M_KB = 0.80
KB = 8.617333262e-5                  # eV/K
EV_A3_PER_GPA = 1.0 / 160.21766      # eV/A^3 per GPa
LDA_UPF = "Fe.pz-spn-kjpaw_psl.1.0.0.UPF"
OUT = "phaseX/ldashift"
FRAMES = f"{OUT}/frames.xyz"


def npt_class():
    try:
        from ase.md.npt import NPT
    except ImportError:
        from ase.md.melchionna import MelchionnaNPT as NPT
    return NPT


def mean_pressure(at, nsteps, every=50):
    dyn = Langevin(at, 2 * units.fs, temperature_K=T_K, friction=0.02)
    ps = []
    for _ in range(max(1, nsteps // every)):
        dyn.run(every)
        ps.append(-np.trace(at.get_stress(voigt=False, include_ideal_gas=True)) / 3 / units.GPa)
    return float(np.mean(ps))


def cmd_md(a):
    os.makedirs(OUT, exist_ok=True)
    NPT = npt_class()
    calc = load_calc(a.model, a.device, a.dry_run)
    thr, vref = q6_reference()
    sc = 0.01 if a.dry_run else 1.0
    n_melt, n_npt, n_avg, n_eq, n_gap, n_b = (int(x * sc) for x in (1000, 2000, 1500, 1000, 500, 1000))
    pfactor = (750 * units.fs) ** 2 * 1300 * units.GPa
    frames, meta = [], {}
    for phase in ("hcp", "liquid"):
        at = hcp_ortho((1.01 if phase == "liquid" else 0.99) * V_GUESS)
        at.info["elec_temp"] = T_K
        at.calc = calc
        rng = np.random.RandomState(7 if phase == "hcp" else 8)
        if phase == "liquid":
            MBD(at, temperature_K=2.0 * T_K, rng=rng)
            hot = Langevin(at, 2 * units.fs, temperature_K=2.0 * T_K, friction=0.02)
            for _ in range(4):
                hot.run(n_melt)
                if solid_fraction(at, thr, vref) < 0.2:
                    break
        else:
            MBD(at, temperature_K=T_K, rng=rng)
        dyn = NPT(at, 2 * units.fs, temperature_K=T_K, externalstress=P_GPA * units.GPa,
                  ttime=25 * units.fs, pfactor=pfactor)
        dyn.run(n_npt)                                         # equilibrate
        cells = []
        for _ in range(max(1, n_avg // 50)):
            dyn.run(50)
            cells.append(np.array(at.cell))
        cell = np.mean(cells, axis=0)
        at.set_cell(cell, scale_atoms=True)
        v = at.get_volume() / len(at)
        # isothermal bulk modulus of this phase from the potential: NVT at V (1 +- 0.01)
        pm = {}
        for s in (+1, -1):
            b = at.copy(); b.info["elec_temp"] = T_K; b.calc = calc
            b.set_cell(cell * (1 + s * 0.01) ** (1 / 3), scale_atoms=True)
            MBD(b, temperature_K=T_K, rng=np.random.RandomState(50 + s))
            Langevin(b, 2 * units.fs, temperature_K=T_K, friction=0.02).run(n_b // 2)
            pm[s] = mean_pressure(b, n_b)
        B = -(pm[+1] - pm[-1]) / 0.02                          # GPa, B = -V dP/dV
        s0 = at.get_scaled_positions(wrap=False).copy()
        start = at.copy()
        kept = 0
        for run in range(RUNS):
            r = start.copy(); r.info["elec_temp"] = T_K; r.calc = calc
            MBD(r, temperature_K=T_K, rng=np.random.RandomState(100 * (run + 1) + len(phase)))
            nvt = Langevin(r, 2 * units.fs, temperature_K=T_K, friction=0.02)
            nvt.run(n_eq)
            got, tries = 0, 0
            while got < SNAPS and tries < SNAPS + 6:
                nvt.run(n_gap)                                 # 1 ps between snapshots
                tries += 1
                fs = solid_fraction(r, thr, vref)
                msd = frac_msd(r, s0) if phase == "hcp" else -1.0
                ok = (fs > 0.8 and msd < 0.6) if phase == "hcp" else fs < 0.2
                if not (ok or a.dry_run):
                    print(f"  {phase} run {run} snap {tries}: solid_frac {fs:.2f} MSD {msd:.2f} -> skipped")
                    continue
                r.calc.reset() if hasattr(r.calc, "reset") else None
                e_mace = float(r.get_potential_energy())
                # static (virial) pressure and forces of the potential, for the MACE - PBE check in eval;
                # QE's P= line is also static (no ideal-gas term)
                p_mace = float(-np.trace(r.get_stress(voigt=False)) / 3 / units.GPa)
                f_mace = np.array(r.get_forces())
                f = Atoms(r.get_chemical_symbols(), positions=r.get_positions(), cell=np.array(r.cell), pbc=True)
                f.info.update(phase=phase, name=f"ldaX_{phase[0]}{kept:02d}", T=T_K, tel=T_K,
                              vat=v, solid_frac=fs, msd=msd, E_mace=e_mace, P_mace=p_mace)
                f.arrays["F_mace"] = f_mace
                frames.append(f)
                kept += 1; got += 1
        meta[phase] = dict(V_at=v, B_GPa=B, P_plus=pm[+1], P_minus=pm[-1], n=kept)
        print(f"{phase}: V = {v:.4f} A^3/atom, B_T = {B:.0f} GPa, {kept} snapshots", flush=True)
    write(FRAMES, frames)
    json.dump(meta, open(f"{OUT}/md_meta.json", "w"), indent=1)
    print(f"written {FRAMES} ({len(frames)} snapshots) and {OUT}/md_meta.json")
    if not a.no_inputs:
        cmd_inputs(a)


def lda_input(pbe_text, pbe_upf):
    t = pbe_text.replace(pbe_upf, LDA_UPF)
    t = re.sub(r"^\s*input_dft\s*=.*\n", "", t, flags=re.M | re.I)
    return t


def cmd_inputs(a):
    frames = read(FRAMES, index=":")
    names, tpls = [], {}
    for f in frames:
        ph = f.info["phase"]
        if ph not in tpls:
            tpls[ph] = template_for(ph, len(f), a.labels_dir)
        tel = float(f.info["tel"])
        stem_p = f"{f.info['name']}_PBE_Tel{int(tel)}"
        txt = make_input(tpls[ph], f, stem_p, tel)
        upf = re.search(r"ATOMIC_SPECIES\s*\n\s*\S+\s+\S+\s+(\S+)", txt).group(1)
        pd = os.path.expanduser(re.search(r"pseudo_dir\s*=\s*'([^']+)'", txt).group(1))
        if not os.path.isabs(pd):
            pd = os.path.join(a.labels_dir, pd)
        if not os.path.exists(os.path.join(pd, LDA_UPF)):
            sys.exit(f"LDA pseudopotential missing: {os.path.join(pd, LDA_UPF)} (see RUN_LDASHIFT.md step 0)")
        stem_l = f"{f.info['name']}_LDA_Tel{int(tel)}"
        open(os.path.join(a.labels_dir, stem_p + ".in"), "w").write(txt)
        open(os.path.join(a.labels_dir, stem_l + ".in"), "w").write(
            lda_input(txt.replace(stem_p, stem_l), upf))
        names += [stem_p + ".in", stem_l + ".in"]
    open(os.path.join(a.labels_dir, "ldaX_list.txt"), "w").write("\n".join(names) + "\n")
    path, origin = write_array_sbatch(a.labels_dir, "ldaX_list.txt", "ldaX_array.sbatch",
                                      "ldaX", len(names), throttle=40)
    print(f"wrote {len(names)} inputs ({len(frames)} snapshots x PBE, LDA), {path} ({origin})")
    print(f"  estimated cost ~{45 * len(names):,} core-hours")


def xc_of(path):
    t = open(path, errors="ignore").read()
    m = re.search(r"Exchange-correlation\s*=\s*(.+)", t)
    return m.group(1).strip() if m else "?"


def mace_check(rows, res, mb, mbe, dS):
    """Potential vs PBE on the same snapshots, at the reported state point (330 GPa, 6,329 K).
    Per phase: static pressure error and force RMSE. Then the liquid - hcp energy bias
    mb = <E_MACE - F_PBE>_liq - <E_MACE - F_PBE>_hcp (per atom) as a T_m error,
    T_m(potential) - T_m(PBE) = +mb / dS_m, split into the part from the pressure error,
    -<P_MACE - P_PBE> (V_liq - V_hcp), and the equal-volume rest. The pressure part equals
    -(dT_m/dP) <dP>: a potential that is too stiff by dP runs at a PBE pressure dP lower."""
    beta = 1.0 / (KB * T_K)
    out = dict(dE_bias=mb, dE_bias_se=mbe, T_m_shift_K=mb / dS, T_m_shift_se_K=mbe / dS)
    print(f"\n  Potential (MACE) vs PBE on the same snapshots, {P_GPA:.0f} GPa, {T_K:.0f} K")
    if not all("mace_dP" in r for r in rows):
        print("    frames.xyz has no stored pressures/forces of the potential (md step from an older "
              "ldashift.py); only the energy bias is available")
        print(f"    liquid - hcp energy bias -> T_m(potential) - T_m(PBE) = {mb / dS:+.0f} +- {mbe / dS:.0f} K")
        return out
    for ph in ("hcp", "liquid"):
        R = [r for r in rows if r["phase"] == ph]
        dp = np.array([r["mace_dP"] for r in R])
        rmse = float(np.sqrt(np.mean([r["mace_f_mse"] for r in R])))
        fref = float(np.sqrt(np.mean([r["pbe_f_ms"] for r in R])))
        cell_bias = np.array([r["mace"] for r in R]) * R[0]["n"]                 # eV per cell
        a2 = float(-0.5 * beta * cell_bias.var(ddof=1) / R[0]["n"])              # MACE -> PBE, eV/atom
        out[ph] = dict(n=len(R), dP=float(dp.mean()), dP_se=float(dp.std(ddof=1) / np.sqrt(len(dp))),
                       dP_sd=float(dp.std(ddof=1)), F_rmse=rmse, F_rms_pbe=fref, a2=a2)
        pct = f"{100 * rmse / fref:.1f}%" if fref > 0 else "n/a"
        print(f"    {ph:6s} n={len(R):2d}  P_MACE - P_PBE = {dp.mean():+.2f} +- {out[ph]['dP_se']:.2f} GPa "
              f"(sd {out[ph]['dP_sd']:.2f})   force RMSE {1000 * rmse:.0f} meV/A "
              f"({pct} of the PBE rms force, {fref:.2f} eV/A)   "
              f"second-order term {1000 * a2:+.2f} meV/atom")
    dPm = 0.5 * (out["hcp"]["dP"] + out["liquid"]["dP"])
    dVm = res["liquid"]["V"] - res["hcp"]["V"]
    e_p = -dPm * dVm * EV_A3_PER_GPA
    e_s = mb - e_p
    d2 = out["liquid"]["a2"] - out["hcp"]["a2"]
    print(f"    liquid - hcp energy bias {1000 * mb:+.2f} +- {1000 * mbe:.2f} meV/atom "
          f"-> T_m(potential) - T_m(PBE) = {mb / dS:+.0f} +- {mbe / dS:.0f} K")
    print(f"      from the pressure error, -<dP> dV_m = {1000 * e_p:+.2f} meV/atom ({e_p / dS:+.0f} K)"
          f"  [mean dP {dPm:+.2f} GPa x dV_m {dVm:.3f} A^3/atom]")
    print(f"      equal-volume part                   {1000 * e_s:+.2f} meV/atom ({e_s / dS:+.0f} K)")
    print(f"      second-order term, liquid - hcp     {1000 * d2:+.2f} meV/atom ({-d2 / dS:+.0f} K)")
    sh = (mb - d2) / dS                                    # T_m(potential) - T_m(PBE), both orders
    tp = 6329 - sh
    print(f"    total, first + second order: T_m(potential) - T_m(PBE) = {sh:+.0f} +- {mbe / dS:.0f} K")
    print(f"    PBE T_m(330 GPa) with the potential's error removed: 6,329 {-sh:+.0f} = {tp:,.0f} "
          f"+- {np.hypot(26, mbe / dS):.0f} K")
    out.update(dP_mean=dPm, dV_m=dVm, pressure_part=e_p, structure_part=e_s, second_order=d2,
               T_m_pressure_part_K=e_p / dS, T_m_structure_part_K=e_s / dS, T_m_shift_total_K=sh,
               T_m_PBE_corrected=tp)
    return out


def cmd_eval(a):
    frames = read(FRAMES, index=":")
    meta = json.load(open(f"{OUT}/md_meta.json"))
    rows = []
    for f in frames:
        n = len(f)
        tel = int(f.info["tel"])
        pp = os.path.join(a.labels_dir, f"{f.info['name']}_PBE_Tel{tel}.out")
        pl = os.path.join(a.labels_dir, f"{f.info['name']}_LDA_Tel{tel}.out")
        rp, rl = parse_out(pp, n), parse_out(pl, n)
        if rp is None or rl is None:
            continue
        xp, xl = xc_of(pp).upper(), xc_of(pl).upper()
        is_pbe = lambda s: "PBE" in s or "PBX" in s
        if not is_pbe(xp) or "PZ" not in xl or is_pbe(xl):
            sys.exit(f"functional check failed for {f.info['name']}: PBE run '{xp}', LDA run '{xl}'")
        row = dict(phase=f.info["phase"], n=n, dU=rl[0] - rp[0], dP=rl[2] - rp[2],
                   P_pbe=rp[2], mace=(float(f.info["E_mace"]) - rp[0]) / n)
        if "P_mace" in f.info and "F_mace" in f.arrays and np.shape(rp[1]) == f.arrays["F_mace"].shape:
            df = np.array(f.arrays["F_mace"]) - rp[1]
            row.update(mace_dP=float(f.info["P_mace"]) - rp[2],
                       mace_f_mse=float(np.mean(df ** 2)), pbe_f_ms=float(np.mean(rp[1] ** 2)))
        rows.append(row)
    if not rows:
        sys.exit("no finished PBE/LDA output pairs")
    beta = 1.0 / (KB * T_K)
    res = {}
    print(f"LDA - PBE free-energy perturbation at {P_GPA:.0f} GPa, {T_K:.0f} K (T_el = T)")
    for ph in ("hcp", "liquid"):
        R = [r for r in rows if r["phase"] == ph]
        if len(R) < 2:
            sys.exit(f"need at least 2 finished pairs for {ph}")
        N = R[0]["n"]
        dU = np.array([r["dU"] for r in R])
        dP = np.array([r["dP"] for r in R])
        a1 = dU.mean() / N
        a1_se = dU.std(ddof=1) / np.sqrt(len(dU)) / N
        a2 = -0.5 * beta * dU.var(ddof=1) / N
        V, B = meta[ph]["V_at"], meta[ph]["B_GPa"]
        a3 = -(V / (2 * B * EV_A3_PER_GPA)) * (dP.mean() * EV_A3_PER_GPA) ** 2
        res[ph] = dict(n=len(R), a1=a1, a1_se=a1_se, a2=a2, a3=a3, dG=a1 + a2 + a3,
                       dP=float(dP.mean()), P_pbe=float(np.mean([r["P_pbe"] for r in R])),
                       mace_bias=float(np.mean([r["mace"] for r in R])),
                       mace_bias_se=float(np.std([r["mace"] for r in R], ddof=1) / np.sqrt(len(R))), V=V, B=B)
        print(f"  {ph:6s} n={len(R):2d}  a1 = {1000 * a1:+9.2f} +- {1000 * a1_se:.2f} meV/atom   "
              f"a2 = {1000 * a2:+.2f}   a3 = {1000 * a3:+.2f} meV/atom   "
              f"P_LDA - P_PBE = {dP.mean():+.1f} GPa (PBE static P {res[ph]['P_pbe']:.1f} GPa)   "
              f"V = {V:.3f} A^3/atom, B = {B:.0f} GPa")
    dS = DS_M_KB * KB
    d1 = res["liquid"]["a1"] - res["hcp"]["a1"]
    se = np.hypot(res["liquid"]["a1_se"], res["hcp"]["a1_se"])
    dg = res["liquid"]["dG"] - res["hcp"]["dG"]
    dT1, dT, dT_se = d1 / dS, dg / dS, se / dS        # raising the liquid's G raises T_m
    print(f"\n  liquid - hcp, first order:      {1000 * d1:+.2f} +- {1000 * se:.2f} meV/atom -> dT_m = {dT1:+.0f} +- {dT_se:.0f} K")
    print(f"  liquid - hcp, with a2 and a3:   {1000 * dg:+.2f} meV/atom -> dT_m(LDA - PBE) = {dT:+.0f} +- {dT_se:.0f} K")
    # Clausius-Clapeyron split of the first-order term:
    #   pressure part  = -<dP> (V_liq - V_hcp)   (the correction changes P at fixed V)
    #   structure part = the rest (liquid - hcp difference at equal volume)
    dPm = 0.5 * (res["liquid"]["dP"] + res["hcp"]["dP"])
    dVm = res["liquid"]["V"] - res["hcp"]["V"]
    e_p = -dPm * dVm * EV_A3_PER_GPA
    e_s = d1 - e_p
    print(f"  split of the first-order term:  pressure part {1000 * e_p:+.2f} meV/atom ({e_p / dS:+.0f} K)"
          f"  [mean dP {dPm:+.1f} GPa x dV {dVm:.3f} A^3/atom]")
    print(f"                                  equal-volume structure part {1000 * e_s:+.2f} meV/atom ({e_s / dS:+.0f} K)")
    print(f"  T_m(LDA, 330 GPa) ~ 6,329 {dT:+.0f} = {6329 + dT:,.0f} K")
    gap = 6225 - 6329
    print(f"\n  Rao & Zhu (LDA+DMFT) - this work (PBE) = {gap:+d} K (+- {np.hypot(42, 26):.0f} K)")
    print(f"  of which LDA - PBE (this calculation):  {dT:+.0f} +- {dT_se:.0f} K")
    print(f"  remainder (DMFT + cell size + potential + protocol): {gap - dT:+.0f} +- "
          f"{np.hypot(np.hypot(42, 26), dT_se):.0f} K")
    mb = res["liquid"]["mace_bias"] - res["hcp"]["mace_bias"]
    mbe = np.hypot(res["liquid"]["mace_bias_se"], res["hcp"]["mace_bias_se"])
    print(f"\n  check: MACE - PBE liquid - hcp energy bias at this state point "
          f"{1000 * mb:+.2f} +- {1000 * mbe:.2f} meV/atom (held-out labels gave +1.0 +- 0.2)")
    val = mace_check(rows, res, mb, mbe, dS)
    json.dump(dict(results=res, dT_first_order=dT1, dT=dT, dT_se=dT_se, rao_zhu_gap=gap,
                   remainder=gap - dT, dT_pressure_part=e_p / dS, dT_structure_part=e_s / dS,
                   mace_validation=val),
              open(f"{OUT}/ldashift_result.json", "w"), indent=1)
    print(f"\nwritten {OUT}/ldashift_result.json")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["md", "inputs", "eval"])
    ap.add_argument("--model", default=MODEL)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--labels-dir", default="labels")
    ap.add_argument("--no-inputs", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="Lennard-Jones stand-in, 1/100 steps")
    a = ap.parse_args()
    {"md": cmd_md, "inputs": cmd_inputs, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    main()
