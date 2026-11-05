#!/bin/bash
#SBATCH --job-name=fe_fill
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-6%7
#SBATCH --output=logs/tether_%A_%a.out

cd ~/scw/fe_melt
PY=""
for c in $(cat env_candidates.txt 2>/dev/null); do
  source "$c" 2>/dev/null
  if python -c "import mace" 2>/dev/null; then PY=python; break; fi
  deactivate 2>/dev/null
done
if [ -z "$PY" ]; then
  for p in python3.13 python3.12 python3.11 python3; do
    if command -v $p >/dev/null && $p -c "import mace" 2>/dev/null; then PY=$p; break; fi
  done
fi
[ -z "$PY" ] && { echo "NO MACE FOUND"; exit 1; }
$PY phaseC_fill.py --model ~/scw/fe_melt/MACE-matpes-pbe-omat-ft.model --job $SLURM_ARRAY_TASK_ID --device cuda
