#!/bin/bash
# fe_phaseE.sbatch — Phase E production array, self-resubmitting (fe_melt_v1)
# 36 runs = 6 cells x 3 initial temperatures x 2 seeds, 200 ps each.
# Each task runs a ~3.5 h segment, checkpoints, and resubmits itself until the
# run reports status "complete". Short segments backfill far better than one
# 11 h block on a saturated 80 GB GPU fleet.
#   sbatch fe_phaseE.sbatch
#SBATCH --job-name=fe_coex
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --exclude=htc-g[003-004],htc-g[032-035],htc-g[037-038],htc-g[041-043],htc-g[048-049]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-35
#SBATCH --output=logs/coex_%A_%a.out

module purge; module load Anaconda3; unset PYTHONPATH
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd ~/fe_melt
M=phaseD/models/fe_ft_Tel6000_f32.model
CELLS=(coex_hcp_100GPa coex_fcc_100GPa coex_hcp_150GPa coex_hcp_200GPa coex_hcp_250GPa coex_hcp_300GPa coex_hcp_330GPa coex_hcp_350GPa)
TI=(0.92 1.00 1.08)
TID=$SLURM_ARRAY_TASK_ID
CELL=${CELLS[$((TID / 6))]}
TIF=${TI[$(( (TID % 6) / 2 ))]}
SEED=$((TID % 2))
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "task $TID -> $CELL tinit=$TIF seed=$SEED"

python3.13 phaseE_coex.py --model $M --cell phaseE/cells/$CELL \
  --tinit-frac $TIF --seed $SEED --ps 200 --wall-hours 3.4 --device cuda
RC=$?
if [ $RC -eq 64 ]; then
  echo "segment done, resubmitting task $TID"
  sbatch --array=$TID ~/fe_melt/fe_phaseE.sbatch
elif [ $RC -ne 0 ]; then
  echo "ERROR rc=$RC for task $TID"
fi
