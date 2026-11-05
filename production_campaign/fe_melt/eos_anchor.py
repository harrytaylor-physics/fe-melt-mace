#!/usr/bin/env python3
"""
eos_anchor.py - cold and hot equation of state of hcp Fe near the 330 GPa coexistence
volume, PBE and LDA, 2-atom cell. Places our PBE pressure against the experimental
300 K equations of state and gives the crystal-only LDA - PBE pressure difference.

Why
  A correction to the exchange-correlation treatment changes T_m at fixed P in two
  ways: through the pressure it adds at fixed volume (dP), and through a liquid - hcp
  difference at equal volume. The first part is
        dT_m(pressure) = -(dT_m/dP) * dP,   dT_m/dP = 8.4 K/GPa at 330 GPa (our Simon fit)
  and needs only the crystal. This script measures dP(PBE -> experiment) and
  dP(PBE -> LDA) at the coexistence volume.

Cases (all hcp, 2 atoms, c/a = 1.60)
  volumes 6.60 6.75 6.90 7.00 7.10 7.25 7.40 A^3/atom
  cold : Marzari-Vanderbilt, degauss 0.01 Ry, k 20x20x13   (compare with 300 K experiment)
  hot  : Fermi-Dirac at T_el = 6,329 K,      k 14x14x9    (static lattice, hot electrons)
  PBE and LDA (pslibrary 1.0.0 spn PAW, the same pair as the LDA - PBE test)
  -> 28 pw.x runs, each a few core-hours.

Experimental 300 K EOS used (hcp Fe)
  Dewaele et al. 2006 PRL 97 215504   Vinet V0 = 75.7 bohr^3, K0 = 163.4 GPa, K0' = 5.38 (data to 205 GPa)
  Fei et al. 2016 GRL                 BM3 K0 = 172.7 GPa, K0' = 4.79, V0 = 11.214 A^3 (assumed 22.428/2)
  Hirao et al. 2022 MRE 7 038403      Vinet V0 = 11.215 A^3, K0 = 159.27 GPa, K0' = 5.86 (data to 354 GPa, Pt scale)
  Vibrational pressure (zero point + 300 K) added to the static PBE value: 3.5 +- 1 GPa
  (Debye-Grueneisen estimate, theta ~ 850 K, gamma ~ 1.5 at 7 A^3/atom).

Usage (~/fe_melt, login node)
  python eos_anchor.py inputs        writes labels/eosA_*.in, labels/eosA_list.txt, labels/eosA_array.sbatch
  cd labels && sbatch eosA_array.sbatch
  python eos_anchor.py eval          after the 28 runs finish
"""
import argparse, glob, json, os, re, sys
import numpy as np

KB_RY = 8.617333262e-5 / 13.605693123
VOLS = [6.60, 6.75, 6.90, 7.00, 7.10, 7.25, 7.40]
COA = 1.60
T_HOT = 6329.0
CASES = {"cold": dict(smear="mv", degauss=0.01, k=(20, 20, 13)),
         "hot": dict(smear="fd", degauss=KB_RY * T_HOT, k=(14, 14, 9))}
UPF = {"PBE": None, "LDA": "Fe.pz-spn-kjpaw_psl.1.0.0.UPF"}
DTDP = 8.4            # K/GPa at 330 GPa, Simon fit to the seven hcp cells
P_VIB = 3.5           # GPa, zero point + 300 K, estimate
V_REPORT = [6.86, 6.97, 7.15]   # hcp at 350, ~330, 300 GPa on our curve (A^3/atom)
B3 = 0.148184711


def vinet(V, V0, K0, Kp):
    x = (V / V0) ** (1 / 3)
    return 3 * K0 * (1 - x) / x ** 2 * np.exp(1.5 * (Kp - 1) * (1 - x))


def bm3(V, V0, K0, Kp):
    f = (V0 / V) ** (2 / 3)
    return 1.5 * K0 * (f ** 3.5 - f ** 2.5) * (1 + 0.75 * (Kp - 4) * (f - 1))


EXP = {"Dewaele 2006": lambda V: vinet(V, 75.7 * B3, 163.4, 5.38),
       "Fei 2016": lambda V: bm3(V, 11.214, 172.7, 4.79),
       "Hirao 2022": lambda V: vinet(V, 11.215, 159.27, 5.86)}


def template_settings(labels_dir):
    tpls = sorted(glob.glob(os.path.join(labels_dir, "hcp-*_Tel*.in")))
    if not tpls:
        sys.exit(f"no campaign hcp input in {labels_dir}/")
    t = open(tpls[0]).read()
    get = lambda k: re.search(rf"{k}\s*=\s*([^,\n]+)", t).group(1).strip()
    pseudo_dir = get("pseudo_dir").strip("'\"")
    upf = re.search(r"^\s*Fe\s+[\d.eEdD+-]+\s+(\S+)", t, re.M).group(1)
    return dict(src=tpls[0], ecutwfc=get("ecutwfc"), ecutrho=get("ecutrho"), pseudo_dir=pseudo_dir, upf=upf)


def write_input(stem, s, upf, V, case):
    a = (4.0 * V / (np.sqrt(3) * COA)) ** (1 / 3)
    c = CASES[case]
    return f"""&CONTROL
  calculation = 'scf', prefix = '{stem}', outdir = './tmp_eosA/{stem}',
  pseudo_dir = '{s['pseudo_dir']}', tstress = .true., tprnfor = .true., disk_io = 'nowf'
/
&SYSTEM
  ibrav = 4, celldm(1) = {a / 0.529177210903:.10f}, celldm(3) = {COA:.6f},
  nat = 2, ntyp = 1, ecutwfc = {s['ecutwfc']}, ecutrho = {s['ecutrho']},
  occupations = 'smearing', smearing = '{c['smear']}', degauss = {c['degauss']:.8f}, nbnd = 32
/
&ELECTRONS
  conv_thr = 1.0d-10, mixing_beta = 0.3
/
ATOMIC_SPECIES
  Fe 55.845 {upf}
ATOMIC_POSITIONS crystal
  Fe 0.3333333333 0.6666666667 0.25
  Fe 0.6666666667 0.3333333333 0.75
K_POINTS automatic
  {c['k'][0]} {c['k'][1]} {c['k'][2]} 0 0 0
"""


def cmd_inputs(a):
    from phaseF_common import write_array_sbatch
    s = template_settings(a.labels_dir)
    upf = dict(UPF, PBE=s["upf"])
    pd = os.path.expanduser(s["pseudo_dir"])
    for xc, u in upf.items():
        if not os.path.exists(os.path.join(pd, u)):
            sys.exit(f"missing {u} in {pd} (step 0 of RUN_LDASHIFT.md downloads the LDA file)")
    names = []
    for xc in upf:
        for case in CASES:
            for V in VOLS:
                stem = f"eosA_{xc}_{case}_V{V:.2f}"
                open(os.path.join(a.labels_dir, stem + ".in"), "w").write(write_input(stem, s, upf[xc], V, case))
                names.append(stem + ".in")
    open(os.path.join(a.labels_dir, "eosA_list.txt"), "w").write("\n".join(names) + "\n")
    path, origin = write_array_sbatch(a.labels_dir, "eosA_list.txt", "eosA_array.sbatch", "eosA", len(names), throttle=28)
    print(f"settings from {os.path.basename(s['src'])}: ecutwfc {s['ecutwfc']}, ecutrho {s['ecutrho']}, PBE {upf['PBE']}, LDA {upf['LDA']}")
    print(f"wrote {len(names)} inputs, {path} ({origin})")


def parse(path):
    if not os.path.exists(path):
        return None
    t = open(path, errors="ignore").read()
    if "JOB DONE" not in t:
        return None
    P = float(re.findall(r"P=\s*(-?[\d.]+)", t)[-1]) / 10.0
    xc = re.search(r"Exchange-correlation\s*=\s*(.+)", t)
    return P, (xc.group(1).upper() if xc else "?")


def p_at(Vs, Ps, V):
    c = np.polyfit(np.log(Vs), np.log(Ps), 2)
    return float(np.exp(np.polyval(c, np.log(V))))


def cmd_eval(a):
    P = {}
    for xc in UPF:
        for case in CASES:
            Vs, Ps = [], []
            for V in VOLS:
                r = parse(os.path.join(a.labels_dir, f"eosA_{xc}_{case}_V{V:.2f}.out"))
                if r is None:
                    continue
                isp = "PBE" in r[1] or "PBX" in r[1]
                if (xc == "PBE") != isp:
                    sys.exit(f"functional check failed: eosA_{xc}_{case}_V{V:.2f} reports '{r[1]}'")
                Vs.append(V); Ps.append(r[0])
            if len(Vs) < 4:
                sys.exit(f"need at least 4 finished volumes for {xc} {case}, have {len(Vs)}")
            P[(xc, case)] = (np.array(Vs), np.array(Ps))
    out = {}
    print("Static pressure (GPa) at the hcp coexistence volumes")
    print(f"{'':28s}" + "".join(f"{V:>10.2f}" for V in V_REPORT))
    for k, (Vs, Ps) in P.items():
        print(f"{k[0] + ' ' + k[1]:28s}" + "".join(f"{p_at(Vs, Ps, V):10.1f}" for V in V_REPORT))
    print(f"\nPBE cold + {P_VIB} GPa (zero point + 300 K) vs experiment at 300 K:")
    for name, f in EXP.items():
        row = []
        for V in V_REPORT:
            dP = f(V) - (p_at(*P[("PBE", "cold")], V) + P_VIB)
            row.append(dP)
        out[name] = row
        print(f"  {name:13s} P_exp - P_PBE = " + ", ".join(f"{d:+.1f}" for d in row)
              + f" GPa  -> dT_m(pressure part) = " + ", ".join(f"{-DTDP * d:+.0f}" for d in row) + " K")
    print("\nLDA - PBE at fixed volume (crystal, static lattice):")
    for case in CASES:
        row = [p_at(*P[("LDA", case)], V) - p_at(*P[("PBE", case)], V) for V in V_REPORT]
        out[f"LDA-PBE {case}"] = row
        print(f"  {case:4s}: dP = " + ", ".join(f"{d:+.1f}" for d in row)
              + " GPa  -> dT_m(pressure part) = " + ", ".join(f"{-DTDP * d:+.0f}" for d in row) + " K")
    row = [p_at(*P[("PBE", "hot")], V) - p_at(*P[("PBE", "cold")], V) for V in V_REPORT]
    print(f"\nPBE electronic thermal pressure at {T_HOT:.0f} K (static lattice): " + ", ".join(f"{d:+.1f}" for d in row) + " GPa")
    print("Compare the LDA - PBE 'hot' row with the pressure part printed by ldashift.py eval.")
    json.dump(dict(V=V_REPORT, results={k: list(map(float, v)) for k, v in out.items()},
                   P={f"{k[0]}_{k[1]}": [list(map(float, x)) for x in v] for k, v in P.items()}),
              open("eosA_result.json", "w"), indent=1)
    print("written eosA_result.json")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["inputs", "eval"])
    ap.add_argument("--labels-dir", default="labels")
    a = ap.parse_args()
    {"inputs": cmd_inputs, "eval": cmd_eval}[a.cmd](a)


if __name__ == "__main__":
    main()
