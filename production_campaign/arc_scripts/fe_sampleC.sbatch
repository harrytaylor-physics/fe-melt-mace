#!/bin/bash
#SBATCH --job-name=fe_sampleC
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-66%12
#SBATCH --output=logs/sampleC_%A_%a.out

cd ~/scw/fe_melt
PY=""
for c in $(cat env_candidates.txt 2>/dev/null); do
  source "$c" 2>/dev/null
  if python -c "import mace" 2>/dev/null; then PY=python; echo "using venv: $c"; break; fi
  deactivate 2>/dev/null
done
if [ -z "$PY" ]; then
  for p in python3.13 python3.12 python3.11 python3; do
    if command -v $p >/dev/null && $p -c "import mace" 2>/dev/null; then PY=$p; echo "using user-site: $p"; break; fi
  done
fi
[ -z "$PY" ] && { echo "NO MACE FOUND in any env"; exit 1; }

MODEL=~/scw/fe_melt/MACE-matpes-pbe-omat-ft.model
TID=$SLURM_ARRAY_TASK_ID
if   [ $TID -lt 28 ]; then PHASE=hcp;    IDX=$TID
elif [ $TID -lt 40 ]; then PHASE=bcc;    IDX=$((TID-28))
elif [ $TID -lt 49 ]; then PHASE=fcc;    IDX=$((TID-40))
else                       PHASE=liquid; IDX=$((TID-49))
fi
echo "task $TID -> $PHASE/$IDX on $(hostname)"
$PY phaseC_sample.py --model $MODEL --phase $PHASE --run-index $IDX --device cuda
