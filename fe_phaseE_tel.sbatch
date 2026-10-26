#!/bin/bash
# fe_phaseE_tel.sbatch — Phase E with the single Tel-conditioned potential (A26).
# Same 48 runs as the Tel-6000 campaign (8 cells x 3 initial temperatures x
# 2 seeds, 200 ps NVE, self-resubmitting segments), same cells, with the
# electronic temperature of each cell set to that cell's Tel-6000 melting
# temperature (analysis_live_all.json). Runs in ~/tel_env; writes phaseE/runs_tel.
#   sbatch fe_phaseE_tel.sbatch
# 12 Oct: segment wall time is now taken from the job's own remaining time
# (squeue %L minus 0.3 h), so a TimeLimit changed with scontrol while the job
# is pending is respected. 12 Oct (2): each resubmission requests only the time the
# remaining ps need (RATE_PS_PER_H, + 0.5 h start-up, 1-10 h), so one job finishes a run.
#SBATCH --job-name=fe_coxT
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --exclude=htc-g[003-004,041-043]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=0-47
#SBATCH --output=logs/coxT_%A_%a.out

module purge; module load Anaconda3; unset PYTHONPATH
source ~/tel_env/bin/activate; export PYTHONNOUSERSITE=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd ~/fe_melt
M=phaseD/models/fe_tel_e3_stagetwo.model
CELLS=(coex_hcp_100GPa coex_fcc_100GPa coex_hcp_150GPa coex_hcp_200GPa coex_hcp_250GPa coex_hcp_300GPa coex_hcp_330GPa coex_hcp_350GPa)
TI=(0.92 1.00 1.08)
TID=$SLURM_ARRAY_TASK_ID
CELL=${CELLS[$((TID / 6))]}
TIF=${TI[$(( (TID % 6) / 2 ))]}
SEED=$((TID % 2))
TAG=${CELL}_ti${TIF}_s${SEED}
RATE_PS_PER_H=6
TEL=$(python -c "import json; print(round(json.load(open('phaseE/analysis_live_all.json'))['cells']['$CELL']['Tm']))")

# wall time = time left in this job minus 0.3 h (fallback 1.7 h)
LEFT=$(squeue -h -j "$SLURM_JOB_ID" -o %L 2>/dev/null | head -1)
WALL=$(python -c "
import sys
s = sys.argv[1].strip() if len(sys.argv) > 1 else ''
try:
    d = 0
    if '-' in s:
        d, s = s.split('-'); d = int(d)
    p = [int(x) for x in s.split(':')]
    while len(p) < 3: p.insert(0, 0)
    h = d * 24 + p[0] + p[1] / 60 + p[2] / 3600
    print(f'{max(0.3, h - 0.3):.2f}')
except Exception:
    print('1.7')
" "$LEFT")
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "task $TID -> $CELL tinit=$TIF seed=$SEED elec_temp=$TEL  time left $LEFT  wall-hours $WALL"

python phaseE_coex_tel.py --model $M --cell phaseE/cells/$CELL --tinit-frac $TIF --seed $SEED \
  --ps 200 --wall-hours $WALL --elec-temp $TEL --rundir phaseE/runs_tel --device cuda
RC=$?
if [ $RC -eq 64 ]; then
  echo "segment done, resubmitting task $TID"
  NEED=$(python -c "
import math, numpy as np
z = np.load('phaseE/runs_tel/${TAG}_ckpt.npz', allow_pickle=True)
rem = max(0, 200 - int(z['blk']))
print(min(600, max(60, math.ceil((rem / $RATE_PS_PER_H + 0.5) * 60))))
" 2>/dev/null || echo 240)
  echo "remaining time requested: $NEED min"
  ok=0
  for i in 1 2 3; do
    if sbatch --time=$NEED --array=$TID ~/fe_melt/fe_phaseE_tel.sbatch; then ok=1; break; fi
    echo "resubmit attempt $i failed"; sleep 60
  done
  [ $ok -eq 1 ] || echo "RESUBMIT FAILED task $TID"
elif [ $RC -ne 0 ]; then
  echo "ERROR rc=$RC for task $TID"
fi
