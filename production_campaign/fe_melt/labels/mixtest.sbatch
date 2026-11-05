#!/bin/bash
#SBATCH --job-name=fe_mixtest
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --output=mixtest.out
module purge
module load QuantumESPRESSO/7.5-foss-2025a
ulimit -s unlimited
export OMP_NUM_THREADS=1
cd $SLURM_SUBMIT_DIR
T0=$SECONDS
mpirun -np 48 pw.x -nk 4 -i mixtest.in > mixtest_qe.out
echo "TIMING mixtest $((SECONDS-T0)) s = $(( (SECONDS-T0)*48/3600 )) core-h"
grep -c "JOB DONE" mixtest_qe.out
grep "convergence has been achieved" mixtest_qe.out
grep "^!" mixtest_qe.out | tail -1
