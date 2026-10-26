#!/bin/bash
#SBATCH --job-name=fe_labels
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --time=01:30:00
#SBATCH --array=0-7295%40
#SBATCH --output=lab_%A_%a.out

module purge
module load QuantumESPRESSO/7.5-foss-2025a
cd $SLURM_SUBMIT_DIR
INP=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" q_02)
T0=$SECONDS
mpirun -np $SLURM_NTASKS pw.x -nk 4 -i "$INP" > "${INP%.in}.out"
EL=$((SECONDS-T0))
echo "TIMING $INP ${EL} s on 48 cores = $((EL*48/3600)) core-h"
# G-A4: the first completed task's TIMING line, x 7296 inputs, is the
# projected label cost. Kill-switch: if projection > 1.5 M core-hours,
# STOP the array (scancel on the job id) and invoke the 2-Tel amendment.
