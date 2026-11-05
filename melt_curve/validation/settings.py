"""settings.py - method settings quoted with the validation results, read from where they are set.

Each value comes from the input or code that used it, not from a hand-written list:
  refs/dft/*.in, refs/eos_pbe_cold/*.in   QE cut-offs, smearing, k-meshes, volumes, c/a, T_el
  val_dft.py docstring                     QE version of the references
  production_campaign/fe_melt/pseudo/*.UPF                             valence electrons of the PAW datasets
  val_mace.py, cc_check.py, pbe_vs_exp.py  strain, phonon displacement, slab and q-bar-6 settings,
                                           vibrational pressure and its uncertainty
  data/coexistence/cells/*.json                      atoms per coexistence cell
  data/coexistence/runs_tel/ file names              starting temperatures and seeds of the coexistence runs
  production_campaign/fe_melt/labels/ldaX_*.in                         the 330 GPa direct-test snapshots (phase, functional, T_el,
                                           atoms, PAW file, input_dft of the LDA runs)
  production_campaign/fe_melt/labels/ldaX*.sbatch                      QE version of the direct-test runs
  production_campaign/fe_melt/ldashift.py                           state point, runs and snapshots per phase, bulk-modulus strain
  refs/lda_setup/*.json                    volume of the mixed-setup (LDA with the PBE PAW file) check
Output: out/settings.json.
"""
import glob, json, os, re, sys
import numpy as np
import paths

sys.path.insert(0, os.path.dirname(paths.HERE))
import layout  # noqa: E402

RY_K = 157887.51                      # 1 Ry / k_B in K
BOHR_A = 0.529177210903


def qe_in(path):
    t = open(path).read()
    g = lambda k: re.search(rf"\b{k}\s*=\s*'?([^,'\s/]+)", t).group(1)
    k = re.search(r"K_POINTS automatic\s*\n\s*(\d+) (\d+) (\d+)", t).groups()
    nat = int(g("nat"))
    m = re.search(r"CELL_PARAMETERS angstrom\s*\n((?:.*\n){3})", t)
    if m:
        cell = np.array([[float(x) for x in l.split()] for l in m.group(1).strip().splitlines()])
        V, coa = abs(np.linalg.det(cell)) / nat, float(np.linalg.norm(cell[2]) / np.linalg.norm(cell[0]))
    else:                                         # ibrav = 4 (hexagonal): celldm(1) = a in bohr, celldm(3) = c/a
        a = float(re.search(r"celldm\(1\)\s*=\s*([\d.]+)", t).group(1)) * BOHR_A
        coa = float(re.search(r"celldm\(3\)\s*=\s*([\d.]+)", t).group(1))
        V = np.sqrt(3) / 2 * a ** 3 * coa / nat
    return dict(ecutwfc=float(g("ecutwfc")), ecutrho=float(g("ecutrho")), smearing=g("smearing"),
                degauss_Ry=float(g("degauss")), T_smear_K=float(g("degauss")) * RY_K, kmesh="x".join(k),
                nat=nat, V=V, c_over_a=coa)


def src_const(path, name):
    m = re.search(rf"^\s*{name}\s*=\s*([^\n#]+)", open(path).read(), re.M)
    return m.group(1).strip()


def main():
    D = {os.path.basename(f)[:-3]: qe_in(f) for f in glob.glob(os.path.join(paths.DFT_REFS, "*.in"))}
    C = {os.path.basename(f)[:-3]: qe_in(f) for f in glob.glob(os.path.join(paths.EOS_COLD, "*.in"))}
    eos = sorted(v["V"] for k, v in D.items() if k.startswith("eos_"))
    out = dict(dft_refs=dict(
        qe_version=re.search(r"\(QE ([\d.]+)", open(os.path.join(paths.HERE, "val_dft.py")).read()).group(1),
        ecutwfc_Ry=sorted({v["ecutwfc"] for v in D.values()}), ecutrho_Ry=sorted({v["ecutrho"] for v in D.values()}),
        smearing=sorted({v["smearing"] for v in D.values()}),
        T_el_K=sorted({round(v["T_smear_K"], -1) for v in D.values()}),
        kmesh_hcp=sorted({v["kmesh"] for k, v in D.items() if k.startswith(("eos_", "tel_", "coa_", "el_", "ph_"))}),
        kmesh_fcc_bcc=sorted({v["kmesh"] for k, v in D.items() if k.startswith(("fcc_", "bcc_"))}),
        eos_volumes=[round(v, 3) for v in eos], eos_n=len(eos),
        eos_c_over_a=sorted({round(v["c_over_a"], 3) for k, v in D.items() if k.startswith("eos_")}),
        coa_V=sorted({round(v["V"], 3) for k, v in D.items() if k.startswith("coa_")}),
        V_by_kind={kind: sorted({round(v["V"], 3) for k, v in D.items() if k.startswith(kind + "_")})
                   for kind in ("eos", "tel", "fcc", "bcc", "coa", "el", "ph")},
        nat_by_kind={kind: sorted({v["nat"] for k, v in D.items() if k.startswith(kind + "_")})
                     for kind in ("eos", "fcc", "bcc")},
        n_refs=len(D)))
    out["eos_cold"] = dict(ecutwfc_Ry=sorted({v["ecutwfc"] for v in C.values()}),
                           ecutrho_Ry=sorted({v["ecutrho"] for v in C.values()}),
                           smearing=sorted({v["smearing"] for v in C.values()}),
                           degauss_Ry=sorted({v["degauss_Ry"] for v in C.values()}),
                           kmesh=sorted({v["kmesh"] for v in C.values()}),
                           c_over_a=sorted({round(v["c_over_a"], 3) for v in C.values()}), n=len(C))
    out["pseudopotentials"] = {os.path.basename(f): float(re.search(r'z_valence="([^"]+)"', open(f).read()).group(1))
                               for f in sorted(glob.glob(os.path.join(layout.PSEUDO, "*.UPF")))}
    vm = os.path.join(paths.HERE, "val_mace.py")
    vm_txt = open(vm).read()
    out["static"] = dict(strain=float(re.search(r"eps = ([\d.]+)", vm_txt).group(1)),
                         phonon_u_A=float(re.search(r"u = ([\d.]+)", vm_txt).group(1)),
                         density_V_A3=float(re.search(r"rho = 55\.845 \* AMU / ([\d.]+)e-30", vm_txt).group(1)))
    cc = os.path.join(paths.HERE, "cc_check.py")
    nb, rcut, fs, fl = [float(x) for x in src_const(cc, "NB, RCUT, F_SOL, F_LIQ").split(",")]
    nbh = re.search(r"for s in \(([-\d, ]+)\)", open(cc).read()).group(1)
    out["slabs"] = dict(n_slabs=int(nb), q6_rcut_A=rcut, f_solid=fs, f_liquid=fl,
                        pure_neighbours_each_side=max(int(x) for x in nbh.split(",")))
    pv = os.path.join(paths.HERE, "pbe_vs_exp.py")
    pv_txt = open(pv).read()
    out["pbe_vs_exp"] = dict(V_core=float(src_const(pv, "V_CORE")),
                             P_vib_uncertainty_GPa=float(re.search(r"uncertain by about ([\d.]+) GPa", pv_txt).group(1)))
    cells = {os.path.basename(f)[:-5]: json.load(open(f))["N"] for f in glob.glob(os.path.join(layout.CELLS, "*.json"))}
    runs = [re.match(r"coex_hcp_\d+GPa_ti([\d.]+)_s(\d+)_ckpt\.npz$", os.path.basename(f))
            for f in glob.glob(os.path.join(layout.RUNS_TEL, "coex_hcp_*_ckpt.npz"))]
    ti, seeds = sorted({float(m.group(1)) for m in runs}), sorted({int(m.group(2)) for m in runs})
    out["coexistence"] = dict(atoms_per_cell=sorted(set(cells.values())), tinit_frac=ti,
                              n_start_T=len(ti), seeds=seeds, n_seeds=len(seeds))
    lx = [re.match(r"ldaX_([hl])\d+_(PBE|LDA)_Tel(\d+)\.in", os.path.basename(f)).groups()
          for f in glob.glob(os.path.join(layout.LABELS, "ldaX_*.in"))]
    lin = {fn: [open(f).read() for f in glob.glob(os.path.join(layout.LABELS, f"ldaX_*_{fn}_Tel*.in"))]
           for fn in ("PBE", "LDA")}
    upf = lambda txts: sorted({re.search(r"ATOMIC_SPECIES\s*\n\s*\S+\s+\S+\s+(\S+)", x).group(1) for x in txts})
    ls_txt = open(os.path.join(layout.ARC, "ldashift.py")).read()
    P, T = (float(x) for x in re.search(r"^P_GPA, T_K = ([\d.]+), ([\d.]+)", ls_txt, re.M).groups())
    runs, snaps = (int(x) for x in re.search(r"^RUNS, SNAPS = (\d+), (\d+)", ls_txt, re.M).groups())
    out["direct_test_330GPa"] = dict(
        n_hcp=sum(1 for p, fn, t in lx if p == "h" and fn == "PBE"), n_liquid=sum(1 for p, fn, t in lx if p == "l" and fn == "PBE"),
        n_PBE=sum(1 for p, fn, t in lx if fn == "PBE"), n_LDA=sum(1 for p, fn, t in lx if fn == "LDA"),
        T_el_K=sorted({int(t) for *_, t in lx}),
        atoms=sorted({int(re.search(r"\bnat\s*=\s*(\d+)", x).group(1)) for x in lin["PBE"] + lin["LDA"]}),
        P_GPa=P, T_K=T, runs_per_phase=runs, snapshots_per_run=snaps,
        snapshot_gap_ps=float(re.search(r"\((\d+) ps apart\)", ls_txt).group(1)),
        bulk_modulus_strain=float(re.search(r"cell \* \(1 \+ s \* ([\d.]+)\)", ls_txt).group(1)),
        qe_version=sorted({m for f in glob.glob(os.path.join(layout.LABELS, "ldaX*.sbatch"))
                           for m in re.findall(r"QuantumESPRESSO/([\d.]+)", open(f).read())}),
        upf_PBE=upf(lin["PBE"]), upf_LDA=upf(lin["LDA"]),
        input_dft_LDA=sorted({m.group(1) if m else "" for m in (re.search(r"input_dft\s*=\s*'([^']+)'", x) for x in lin["LDA"])}))
    ls = [json.load(open(f)) for f in glob.glob(os.path.join(paths.HERE, "refs", "lda_setup", "*.json"))]
    out["lda_setup"] = dict(V_per_atom=sorted({round(r["V_per_atom"], 3) for r in ls}), n=len(ls))
    print(json.dumps(out, indent=1))
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "settings.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
