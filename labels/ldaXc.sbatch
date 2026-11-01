#!/bin/bash
#SBATCH --job-name=ldaXc
#SBATCH --account=phys-plasma-dft
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=0
#SBATCH --output=ldaXc_%A_%a.log
module purge
module load QuantumESPRESSO/7.5-foss-2025a && command -v pw.x || { echo "QE module missing on $SLURM_CLUSTER_NAME"; exit 1; }
cd $SLURM_SUBMIT_DIR; mkdir -p claims
[ -e claims/STOP_$SLURM_CLUSTER_NAME ] && { echo "stopped on $SLURM_CLUSTER_NAME"; exit 1; }
if [ "$SLURM_CLUSTER_NAME" = arc ]; then RUN="mpirun -np $SLURM_NTASKS"; else RUN=srun; fi
for INP in $(cat ldaX_list2.txt); do
  S=${INP%.in}
  grep -qs "JOB DONE" "$S.out" && continue
  mkdir "claims/$S" 2>/dev/null || continue
  echo "$SLURM_CLUSTER_NAME job $SLURM_JOB_ID $(date '+%d %b %H:%M')" > "claims/$S/by"
  T0=$SECONDS; $RUN pw.x -nk 4 -i "$INP" > "$S.out"
  if ! grep -q "JOB DONE" "$S.out"; then
    echo "$S failed after $((SECONDS-T0)) s"; rm -rf "claims/$S"
    [ $((SECONDS-T0)) -lt 300 ] && touch claims/STOP_$SLURM_CLUSTER_NAME
  fi
  exit 0
done
echo "no unclaimed inputs on $SLURM_CLUSTER_NAME"
