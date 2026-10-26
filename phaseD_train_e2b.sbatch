#!/bin/bash
# phaseD_train_e2.sbatch — E2 (A12/A14): Tel-conditioned MACE trained on the
# Tel-4000 + Tel-8000 labels only; judged on the held-out Tel-6000 test set.
#
# Hyperparameters identical to the production Tel-6000 potential
# (fe_ft_Tel6000_f32): same foundation model, loss weights, lr, batch, EMA,
# epochs, patience, float32. Only differences: two Tels of training data, and
# the electronic-temperature embedding (kelvin in; the 1e-4 scale is stored in
# the checkpoint; zero-initialised so the untrained pathway is the base model).
#
# Runs in the isolated venv ~/tel_env (tel-conditioning branch, commit 5168bc0),
# never the ~/.local production install (A12 rule).
#
# Resumable: --restart_latest. If the 12 h limit ends it before epoch 150,
# resubmit the same file and it continues from the last epoch checkpoint.
#   sbatch phaseD_train_e2.sbatch
#SBATCH --job-name=fe_e2
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --exclude=htc-g[003-004],htc-g[032-035],htc-g[037-038],htc-g[041-043],htc-g[048-049]
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=11:55:00
#SBATCH --output=logs/train_e2b_%j.out

module purge; module load Anaconda3; unset PYTHONPATH
source ~/tel_env/bin/activate
export PYTHONNOUSERSITE=1
cd ~/fe_melt
mkdir -p logs phaseD/models
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python -c "import mace; print('mace from', mace.__file__)"

SPEC='{"elec_temp": {"type": "continuous", "per": "graph", "in_dim": 1, "emb_dim": 32, "scale": 1e-4, "zero_init": True}}'

python -m mace.cli.run_train \
  --name fe_tel_e2b \
  --model_dir phaseD/models --log_dir logs --checkpoints_dir phaseD/models \
  --foundation_model ~/fe_melt/MACE-matpes-pbe-omat-ft.model \
  --multiheads_finetuning False \
  --train_file phaseD/train_e2.xyz \
  --valid_fraction 0.05 \
  --energy_key REF_energy --forces_key REF_forces \
  --elec_temp_key elec_temp \
  --embedding_specs "$SPEC" \
  --E0s average \
  --energy_weight 1.0 --forces_weight 100.0 \
  --lr 1e-4 --batch_size 4 --valid_batch_size 8 \
  --max_num_epochs 150 --patience 30 \
  --ema --ema_decay 0.999 --amsgrad \
  --scaling rms_forces_scaling \
  --error_table PerAtomRMSE \
  --default_dtype float32 --device cuda --seed 7 --restart_latest --swa --start_swa 110 --swa_energy_weight 1000 --swa_forces_weight 100

ls -la phaseD/models/ | grep fe_tel_e2b
