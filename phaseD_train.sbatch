#!/bin/bash
# phaseD_train.sbatch — fine-tune one MACE potential on ARC GPU (fe_melt_v1)
# One submission per electronic temperature:
#   sbatch --export=ALL,TEL=6000 phaseD_train.sbatch     (v0 and primary)
#   sbatch --export=ALL,TEL=4000 phaseD_train.sbatch
#   sbatch --export=ALL,TEL=8000 phaseD_train.sbatch
#
# Design notes (D1):
# - Base: MACE-matpes-pbe-omat-ft.model (per Phase B / amendment A4).
# - multiheads_finetuning is DISABLED: the replay mechanism downloads its
#   replay set at runtime and ARC compute nodes have no internet; and for a
#   single-element production potential, catastrophic forgetting of other
#   chemistry is not a deliverable concern. Logged as the naive-FT decision.
# - Forces weighted 100:1 over energies (forces carry the information; the
#   Mermin energies anchor phase differences).
# - float64 throughout, matching the audit and sampling.

#SBATCH --job-name=fe_ft
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=short
#SBATCH --gres=gpu:h100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/train_Tel%j.out

cd ~/fe_melt
mkdir -p logs phaseD/models
TEL=${TEL:-6000}

# self-locating python env (same block as the sampling jobs)
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

$PY -m mace.cli.run_train \
  --name "fe_ft_Tel${TEL}" \
  --model_dir phaseD/models --log_dir logs --checkpoints_dir phaseD/models \
  --foundation_model ~/fe_melt/MACE-matpes-pbe-omat-ft.model \
  --multiheads_finetuning False \
  --train_file "phaseD/train_Tel${TEL}.xyz" \
  --valid_fraction 0.05 \
  --energy_key REF_energy --forces_key REF_forces \
  --E0s average \
  --energy_weight 1.0 --forces_weight 100.0 \
  --lr 1e-4 --batch_size 4 --valid_batch_size 8 \
  --max_num_epochs 150 --patience 30 \
  --ema --ema_decay 0.999 --amsgrad \
  --scaling rms_forces_scaling \
  --error_table PerAtomRMSE \
  --default_dtype float64 --device cuda --seed 7 --restart_latest

echo "training done: phaseD/models/fe_ft_Tel${TEL}.model"
