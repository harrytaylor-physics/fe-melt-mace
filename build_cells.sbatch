#!/bin/bash
#SBATCH --job-name=fe_cells
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=short
#SBATCH --gres=gpu:h100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-5
#SBATCH --output=logs/cells_%A_%a.out
module purge; module load Anaconda3; unset PYTHONPATH
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd ~/fe_melt
M=phaseD/models/fe_ft_Tel6000_f32.model
P=(300 250 200 150 100 100); B=(hcp hcp hcp hcp hcp fcc)
i=$SLURM_ARRAY_TASK_ID
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
python3.13 phaseE_build.py --model $M --pressure ${P[$i]} --branch ${B[$i]} --device cuda
