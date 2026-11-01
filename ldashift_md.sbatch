#!/bin/bash
# ldashift_md.sbatch — step X1a: 40 snapshots (20 hcp, 20 liquid; 128 atoms) at the PBE
# melting point, 330 GPa and 6,329 K, from the Tel potential, plus each phase's bulk
# modulus; then the 80 pw.x inputs (PBE and LDA on every snapshot) in labels/.
#   cd ~/fe_melt && mkdir -p logs && sbatch ldashift_md.sbatch
#SBATCH --job-name=fe_ldaX
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --exclude=htc-g[003-004,041-043]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/ldaX_md_%j.out

module purge; module load Anaconda3; unset PYTHONPATH
source ~/tel_env/bin/activate; export PYTHONNOUSERSITE=1
cd ~/fe_melt
python ldashift.py md
