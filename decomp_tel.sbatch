#!/bin/bash
# decomp_tel.sbatch — A29 items (1) and (2): first-order T_el / model decomposition
# (decomp_tel.py, three stages) and the 350 GPa energy-conservation check missed on 13 Oct.
#SBATCH --job-name=fe_decomp
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --exclude=htc-g[003-004,041-043]
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=logs/decomp_tel.out
module purge; module load Anaconda3; unset PYTHONPATH
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd ~/fe_melt
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "=== stage tel (tel_env) ==="
( source ~/tel_env/bin/activate; export PYTHONNOUSERSITE=1; python -u decomp_tel.py --stage tel --device cuda 2>&1 | grep -v -i warn ) || exit 1
echo "=== stage ref (production env) ==="
( unset PYTHONNOUSERSITE VIRTUAL_ENV; python3.13 -u decomp_tel.py --stage ref --device cuda 2>&1 | grep -v -i warn ) || exit 1
echo "=== combine ==="
python3.13 -u decomp_tel.py --stage combine
echo "=== energy check 350 GPa (tel_env) ==="
( source ~/tel_env/bin/activate; export PYTHONNOUSERSITE=1; python -u phaseE_analyze_live.py --rundir phaseE/runs_tel --include-drifting --energy-check --echeck-cells coex_hcp_350GPa --model phaseD/models/fe_tel_e3_stagetwo.model --device cuda 2>&1 | sed -n '/ENERGY CONSERVATION/,$p' | grep -v -i -e warn -e reshape -e shape -e elec_temp )
