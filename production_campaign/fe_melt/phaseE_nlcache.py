#!/usr/bin/env python3
"""
phaseE_nlcache.py — Verlet-skin neighbour-list cache for MACE (fe_melt_v1)

MACE's ASE calculator rebuilds the neighbour list from scratch on every force
call. At 12k atoms and ~130 neighbours per atom that is ~1.6M pairs per step,
measured at 257 ms for 6,048 atoms on a 2-core allocation — comparable to the
GPU model evaluation itself, and it does NOT shrink on a faster card.

Standard MD practice is a Verlet skin: build the list with cutoff r_max + skin,
then reuse it until some atom has moved more than skin/2 since the build (at
which point a pair could have entered r_max unseen). Pairs in the skin region
are included but contribute exactly zero, because MACE's radial cutoff function
is identically zero beyond r_max — so this changes no physics, only bookkeeping.

Usage: import and call once, before creating the calculator:
    from phaseE_nlcache import enable_nl_cache
    enable_nl_cache(skin=1.0)

Reports its hit rate on exit if verbose=True.
"""
import numpy as np

_STATE = {"pos": None, "cell": None, "res": None, "cutoff": None,
          "builds": 0, "calls": 0, "skin": 1.0, "verbose": False}


def enable_nl_cache(skin=1.0, verbose=False):
    """Patch MACE's neighbourhood construction with a skin-cached version."""
    import mace.data.atomic_data as ad
    original = ad.get_neighborhood
    _STATE["skin"] = float(skin)
    _STATE["verbose"] = bool(verbose)

    def cached(positions, cutoff, pbc=None, cell=None, true_self_interaction=False):
        s = _STATE
        s["calls"] += 1
        big = float(cutoff) + s["skin"]
        cell_arr = np.asarray(cell) if cell is not None else None
        reuse = (
            s["res"] is not None
            and s["cutoff"] == big
            and s["pos"] is not None
            and positions.shape == s["pos"].shape
            and (cell_arr is None or (s["cell"] is not None
                                      and np.allclose(cell_arr, s["cell"], atol=1e-10)))
        )
        if reuse:
            disp = np.abs(positions - s["pos"])
            # conservative: ignore PBC wrapping, wrapped atoms force a rebuild
            if disp.max() < 0.5 * s["skin"]:
                return s["res"]
        res = original(positions, big, pbc=pbc, cell=cell,
                       true_self_interaction=true_self_interaction)
        s.update(pos=positions.copy(),
                 cell=None if cell_arr is None else cell_arr.copy(),
                 res=res, cutoff=big)
        s["builds"] += 1
        if s["verbose"] and s["builds"] % 50 == 0:
            print(f"[nlcache] {s['builds']} builds / {s['calls']} calls "
                  f"({100*s['builds']/s['calls']:.1f}%)", flush=True)
        return res

    ad.get_neighborhood = cached
    return cached


def nl_cache_stats():
    s = _STATE
    return dict(builds=s["builds"], calls=s["calls"],
                rate=(s["builds"] / s["calls"]) if s["calls"] else float("nan"))
