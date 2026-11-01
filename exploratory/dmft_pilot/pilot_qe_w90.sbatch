#!/bin/bash
# pilot_qe_w90.sbatch — DFT+DMFT pilot, step D0: QE SCF + NSCF and Wannier90 Fe 3d
# orbitals for 2-atom hcp Fe at 7.05 A^3/atom, T = 5,800 K. About 30-60 min.
#   mkdir -p ~/fe_melt/dmft/pilot && cd ~/fe_melt/dmft/pilot  (with the pilot files)
#   sbatch pilot_qe_w90.sbatch
#SBATCH --job-name=dmft_D0
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=D0_%j.out

set -e
module purge
module load QuantumESPRESSO/7.5-foss-2025a
cd $SLURM_SUBMIT_DIR

python3 make_inputs.py scf --labels ~/fe_melt/labels
srun pw.x -nk 8 -i scf.in > scf.out
python3 make_inputs.py w90
srun pw.x -nk 8 -i nscf.in > fe.nscf.out
wannier90.x -pp fe
srun pw2wannier90.x -i pw2wan.in > pw2wan.out
for w in w1 w2; do
  for ext in amn mmn eig nnkp; do ln -sf fe.$ext fe_$w.$ext; done
  ln -sf fe.nscf.out fe_$w.nscf.out
  wannier90.x fe_$w || echo "wannier90 FAILED for window $w: see fe_$w.wout and fe_$w.werr"
  grep -A12 "Final State" fe_$w.wout | tail -12 || true
done
ls -l fe_w?_hr.dat fe_w?_u.mat || true
echo "D0 done. Next: in the TRIQS environment run  python wan_check.py"
