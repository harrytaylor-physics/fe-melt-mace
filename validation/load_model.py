"""load_model.py - the production potential as an ASE calculator, on CPU or GPU.

The model file was saved on a GPU. On a machine without CUDA, torch.jit.load needs
map_location='cpu'; this module patches it before mace is imported.

mace-torch 0.3.16 (upstream) reads the electronic temperature from atoms.info through
info_keys and expects it in units of 10^4 K: set atoms.info["elec_temp"] = T_el / 1e4.
Raw kelvin gives energies that are wrong by orders of magnitude. float32 reproduces the
stored test-set predictions exactly; float64 shifts every energy by a uniform +2.27 meV/atom.
Production scripts on ARC (tel-mace fork) set atoms.info["elec_temp"] = T_el in kelvin; the
fork must therefore apply the 10^4 scaling itself. That is inferred from the exact agreement
of upstream-with-T_el/1e4 against the stored ARC predictions, not yet read from the fork's code.
"""
import os, warnings

warnings.filterwarnings("ignore")
import torch

_jl = torch.jit.load
if not torch.cuda.is_available():
    def _jl_cpu(f, map_location=None, *a, **k):
        return _jl(f, map_location="cpu", *a, **k)
    torch.jit.load = _jl_cpu

from mace.calculators import MACECalculator  # noqa: E402

import paths  # noqa: E402

# atoms.info["elec_temp"] = T_el / TEL_SCALE. 1e4 for upstream mace-torch 0.3.16 (checked against the
# stored ARC predictions). Set FE_TEL_SCALE=1 if the installed mace is a build that scales internally.
TEL_SCALE = float(os.environ.get("FE_TEL_SCALE", "1e4"))


def calc(dtype="float32", device=None, threads=2):
    torch.set_num_threads(threads)
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    return MACECalculator(model_paths=paths.MODEL, device=device, default_dtype=dtype,
                          info_keys={"elec_temp": "elec_temp", "total_spin": "spin", "total_charge": "charge"})
