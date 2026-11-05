# LDA − PBE shift of T_m at 330 GPa (step X1)

Status: complete. The results are in this folder (`phaseX/ldashift/` on ARC): `frames.xyz`, `md_meta.json`,
`ldashift_result.json`, `ldashift_eval.txt`. The report numbers come from
`melt_curve/validation/ldashift_summary.py`. The scripts are `production_campaign/fe_melt/ldashift.py` and
`production_campaign/arc_scripts/ldashift_md.sbatch`; the inputs are `production_campaign/fe_melt/labels/ldaX_*`.
The commands below ran in `~/fe_melt` on ARC and are kept as they ran.

**Question.** How much of the 104 K gap between Rao & Zhu (LDA+DMFT, 6,225 ± 42 K) and this
work (PBE, 6,329 ± 6 K) comes from LDA versus PBE, and how much remains for the correlation
(DMFT) correction plus their other differences (4–8-atom training cells, potential, protocol)?

**Method.**
1. Take 20 bulk hcp and 20 bulk liquid snapshots (128 atoms) at the PBE melting point
   (330 GPa, 6,329 K) from the Tel potential. Each phase is run in NPT to find its mean volume,
   then in NVT at that volume, with snapshots 1 ps apart.
2. On every snapshot, compute the Mermin free energy at T_el = 6,329 K with PBE and with LDA.
   Both use the same cutoffs, k-points, smearing and PAW file, the official pslibrary 1.0.0 PBE
   dataset; the LDA runs add `input_dft = 'PZ'`. (Until 27 Oct the LDA inputs used the PZ dataset of
   step 0, generated with ld1.x; it gives a divergent SCF in these configurations and is not used,
   see `production_campaign/fe_melt/pseudo/lda_generation/README.md`.)
3. Per phase:
   - first-order term ⟨ΔU⟩/N;
   - configurational second-order term −(β/2)·var(ΔU)/N;
   - volume term −V·ΔP²/(2B), which converts fixed volume to fixed pressure (ΔP = P_LDA − P_PBE
     at the snapshot volume; B from the potential).
4. Shift: δT_m = +(ΔG_liq − ΔG_hcp)/ΔS_m, with ΔS_m = 0.80 k_B. A correction that raises the
   liquid's free energy relative to the solid raises T_m.

Cost: 80 pw.x calculations, about 3,600 core-hours, plus one GPU job of under 3 h.

## 0. Once: the LDA pseudopotential, next to the PBE one (login node)

Superseded on 27 Oct: the LDA runs use the PBE file with `input_dft = 'PZ'`, so this step is not needed.

```bash
cd ~/fe_melt && PD=$(grep -h "pseudo_dir" $(ls labels/hcp-*_Tel4000.in | head -1) | sed "s/.*'\(.*\)'.*/\1/") && PD=${PD/#\~/$HOME} && echo "pseudo_dir: $PD" && \
( cd "$PD" && [ -f Fe.pz-spn-kjpaw_psl.1.0.0.UPF ] || wget -q https://pseudopotentials.quantum-espresso.org/upf_files/Fe.pz-spn-kjpaw_psl.1.0.0.UPF ) && \
ls -l "$PD"/Fe.p*-spn-kjpaw_psl.1.0.0.UPF && grep -m1 -i "functional" "$PD"/Fe.pz-spn-kjpaw_psl.1.0.0.UPF
```

The last line must name an LDA functional (PZ, or SLA PZ). If `wget` is blocked on ARC,
download the file on the Mac and `scp` it into that folder.

## 1. Snapshots and inputs (GPU)

```bash
cd ~/fe_melt && tar xzf ldashift.tgz && mkdir -p logs && sbatch ldashift_md.sbatch
```

## 2. The 80 DFT calculations, when step 1 has finished

```bash
cd ~/fe_melt && tail -6 logs/ldaX_md_*.out && cd labels && sbatch ldaX_array.sbatch && cd ~/fe_melt
```

The log should end with "hcp: V = … B_T = … GPa, 20 snapshots", the same line for the liquid,
and "wrote 80 inputs".

## 3. Result, when the calculations have finished (login node)

```bash
cd ~/fe_melt && echo "finished: $(grep -l 'JOB DONE' labels/ldaX_*.out | wc -l) / 80" && \
source ~/tel_env/bin/activate && python ldashift.py eval
```

The output gives:
- the shift δT_m(LDA − PBE) and its standard error;
- the shift split into a pressure part, −⟨P_LDA − P_PBE⟩(V_liq − V_hcp), and an equal-volume
  structure part. If the structure part is small, a functional change acts on T_m only through
  the pressure it adds at fixed volume, and crystal-only calculations are enough for it;
- the remainder of the Rao & Zhu gap after removing it;
- the potential against PBE on the same 40 snapshots, at the exact state point of the result:
  pressure error and force RMSE for each phase, and the liquid − hcp energy bias as a T_m error,
  split into a pressure part and an equal-volume part. These numbers fill item 1 under Open items
  in the validation doc.

Before computing anything, the script checks that every PBE output used PBE and every LDA
output used LDA. Results are written to `phaseX/ldashift/ldashift_result.json`.
