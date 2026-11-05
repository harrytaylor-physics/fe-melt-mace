#!/bin/bash
# env_install.sh — TRIQS 3.3 environment for the DFT+DMFT pilot (fe_melt, step D0a).
# Installs TRIQS, TRIQS/DFTTools and TRIQS/CTHYB (conda-forge) and solid_dmft (PyPI,
# matching 3.3 series) into ${DATA:-$HOME}/envs/triqs33, checks MPI and the imports,
# and packs the installed solid_dmft package (~/solid_dmft_params.tgz) so the DMFT
# inputs can be written against this exact version.
# Run on the login node:  bash env_install.sh 2>&1 | tee env_install.log
# If the conda solve is killed for memory, run the same command in an interactive job.
set -e
module purge
module load Anaconda3
source "$(conda info --base)/etc/profile.d/conda.sh"
ENV=${DATA:-$HOME}/envs/triqs33
if [ ! -d "$ENV" ]; then
  conda create -y -p "$ENV" -c conda-forge --override-channels \
    "triqs=3.3" "triqs_dft_tools=3.3" "triqs_cthyb=3.3" h5py scipy toml matplotlib
fi
conda activate "$ENV"
pip install --no-deps "solid_dmft>=3.3,<3.4" || pip install --no-deps "solid_dmft<4"
python - <<'PY'
import importlib
for m in ("triqs", "triqs_dft_tools", "triqs_cthyb", "solid_dmft", "h5py", "numpy"):
    mod = importlib.import_module(m)
    v = getattr(mod, "__version__", None)
    if v is None:
        try:
            v = importlib.import_module(m + ".version").version
        except Exception:
            v = "?"
    print(f"{m:16s} {v}")
PY
mpirun -n 2 python -c "from triqs.utility import mpi; print('MPI rank', mpi.rank, 'of', mpi.size)"
python -c "from triqs_dft_tools.converters.wannier90 import Wannier90Converter; from triqs_cthyb import Solver; print('Wannier90Converter and CTHYB Solver import: ok')"
D=$(python -c "import solid_dmft, os; print(os.path.dirname(solid_dmft.__file__))")
tar czf ~/solid_dmft_params.tgz -C "$(dirname "$D")" "$(basename "$D")"
echo "environment: $ENV"
echo "next: copy ~/solid_dmft_params.tgz and env_install.log to the Mac and attach both in the chat"
