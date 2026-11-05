#!/usr/bin/env python3
"""
phaseE_q6.py — Lechner–Dellago averaged bond order q̄6 (fe_melt_v1)

Replaces the plain Steinhardt q6, which has no solid/liquid contrast at
melting temperatures (measured on our cells: solid 0.261 vs liquid 0.258 at
300 GPa). The LD average sums each atom's complex q6m vector with its
neighbours' before contracting, cancelling thermal noise while preserving
crystalline coherence. Expected contrast: solid ~0.25-0.35, liquid ~0.05-0.12.

Fully vectorised (bincount over the pair list): ~50x faster than a per-atom
loop, a few seconds for 12k atoms on CPU.

Also a recalibration entry point for cells already built with the old q6:
their structures are valid; only the stored thresholds need replacing. The
sandwich contains both phases, so medians are taken from fractional-z bands
well away from the two interfaces (solid slab occupies z/L in [0, 0.5]):

  python3 phaseE_q6.py recalib phaseE/cells/coex_hcp_300GPa [more prefixes...]
"""
import sys, json
import numpy as np


def q6_per_atom(at, rcut=3.2, device="cpu"):
    from ase.neighborlist import neighbor_list
    try:   # SciPy >= 1.15: sph_harm_y(n, m, polar, azimuth); sph_harm removed in newer releases
        from scipy.special import sph_harm_y
        def sph_harm(m, n, az, pol):          # keep the old argument order used below
            return sph_harm_y(n, m, pol, az)
    except ImportError:
        from scipy.special import sph_harm
    n = len(at)
    i, j, D, r = neighbor_list("ijDd", at, rcut)
    if len(i) == 0:
        return np.zeros(n)
    th = np.arccos(np.clip(D[:, 2] / r, -1.0, 1.0))
    ph = np.arctan2(D[:, 1], D[:, 0])
    cnt = np.bincount(i, minlength=n).astype(float)
    q = np.zeros((n, 13), dtype=complex)
    for k, m in enumerate(range(-6, 7)):
        y = sph_harm(m, 6, ph, th)
        q[:, k] = (np.bincount(i, weights=y.real, minlength=n)
                   + 1j * np.bincount(i, weights=y.imag, minlength=n))
    good = cnt > 0
    q[good] /= cnt[good, None]
    # Lechner-Dellago: average the q6m vectors over {atom + its neighbours}
    Q = q.copy()
    for k in range(13):
        Q[:, k] += (np.bincount(i, weights=q[j, k].real, minlength=n)
                    + 1j * np.bincount(i, weights=q[j, k].imag, minlength=n))
    Q[good] /= (cnt[good] + 1.0)[:, None]
    return np.sqrt(4 * np.pi / 13.0 * (np.abs(Q) ** 2).sum(axis=1))


def recalibrate(prefix):
    from ase.io import read
    at = read(prefix + ".extxyz")
    meta = json.load(open(prefix + ".json"))
    q = q6_per_atom(at)
    zf = at.get_positions()[:, 2] / at.cell.lengths()[2]
    solid_band = (zf > 0.10) & (zf < 0.40)     # solid slab interior
    liquid_band = (zf > 0.60) & (zf < 0.90)    # liquid slab interior
    qs, ql = np.median(q[solid_band]), np.median(q[liquid_band])
    thr = 0.5 * (qs + ql)
    meta.update(q6_solid=float(qs), q6_liquid=float(ql), q6_threshold=float(thr),
                q6_method="LD-qbar6", q6_bands="z/L 0.10-0.40 solid, 0.60-0.90 liquid")
    json.dump(meta, open(prefix + ".json", "w"), indent=1)
    sep = "OK" if qs - ql > 0.10 else "WEAK — inspect before production"
    print(f"{meta['tag']}: qbar6 solid {qs:.3f} / liquid {ql:.3f} "
          f"-> threshold {thr:.3f}   separation {qs-ql:.3f} [{sep}]")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "recalib":
        for p in sys.argv[2:]:
            recalibrate(p)
    else:
        print(__doc__)
