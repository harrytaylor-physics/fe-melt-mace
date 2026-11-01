import sys, torch
S = float(sys.argv[1])
from phaseE_nlcache import enable_nl_cache
enable_nl_cache(skin=S)
from ase.io import read
from mace.calculators import MACECalculator
at = read('phaseE/cells/coex_hcp_300GPa.extxyz')
at.calc = MACECalculator(model_paths='phaseD/models/fe_ft_Tel6000_f32.model', device='cuda', default_dtype='float32', enable_cueq=True)
torch.cuda.reset_peak_memory_stats()
try:
    at.get_forces(); at.get_stress(); at.get_forces()
    print(f'skin {S}: peak {torch.cuda.max_memory_allocated()/2**30:.1f} GiB  {torch.cuda.get_device_name(0)}', flush=True)
except torch.OutOfMemoryError:
    print(f'skin {S}: OOM', flush=True)
