#!/usr/bin/env python3
"""
upf_compare.py - check that two UPF v2 pseudopotential files are the same construction.

Use: python upf_compare.py official.UPF regenerated.UPF

Compares the generation line, functional, valence charge, total pseudo-energy and every
numeric array in the file (radial mesh, local potential, projectors, augmentation
charges, core charge, PAW data). Prints the largest relative difference
(max |a - b| / max |a|) for arrays that differ by more than 1e-6, and ends with MATCH or MISMATCH.

Pass: every array within 1e-4 relative and total pseudo-energy within 1e-4 Ry.
Exit code 0 on MATCH, 1 on MISMATCH, so it can gate the next command with &&.
"""
import re, sys
import numpy as np

TOL_REL, TOL_E = 1e-4, 1e-4


def load(path):
    t = open(path, errors="ignore").read()
    hdr = re.search(r"<PP_HEADER(.*?)/>", t, re.S)
    attrs = dict(re.findall(r'(\w+)\s*=\s*"([^"]*)"', hdr.group(1))) if hdr else {}
    gen = re.search(r"Generated using.*", t)
    arrays = {}
    for tag, body in re.findall(r"<(PP_[A-Z0-9_.]+)[^>]*>([^<]*)</\1>", t):
        try:
            x = np.array(body.replace("D", "E").replace("d", "e").split(), float)
        except ValueError:
            continue
        if x.size > 1:
            arrays[tag] = x
    return attrs, (gen.group(0).strip() if gen else "?"), arrays


def main(a, b):
    ha, ga, A = load(a)
    hb, gb, B = load(b)
    print(f"A: {a}\n   {ga}\n   functional {ha.get('functional', '?').split()}  z_valence {ha.get('z_valence')}")
    print(f"B: {b}\n   {gb}\n   functional {hb.get('functional', '?').split()}  z_valence {hb.get('z_valence')}")
    ok = True
    if ha.get("functional", "").split() != hb.get("functional", "").split():
        print("functional differs"); ok = False
    if float(ha.get("z_valence", 0)) != float(hb.get("z_valence", 0)):
        print("z_valence differs"); ok = False
    try:
        dE = abs(float(ha["total_psenergy"]) - float(hb["total_psenergy"]))
        print(f"total_psenergy: {float(ha['total_psenergy']):.8f} vs {float(hb['total_psenergy']):.8f} Ry  (diff {dE:.2e})")
        ok &= dE < TOL_E
    except KeyError:
        print("total_psenergy missing in one file")
    only = sorted(set(A) ^ set(B))
    if only:
        print(f"arrays in one file only: {only}"); ok = False
    worst, shown = 0.0, 0
    common = sorted(set(A) & set(B))
    print(f"{len(common)} numeric arrays compared; listed below are those differing by more than 1e-6 relative")
    for k in common:
        if A[k].size != B[k].size:
            print(f"{k:28s} size {A[k].size} vs {B[k].size}"); ok = False; continue
        s = np.max(np.abs(A[k]))
        r = np.max(np.abs(A[k] - B[k])) / s if s > 0 else np.max(np.abs(B[k]))
        worst = max(worst, r)
        if r > 1e-6:
            print(f"  {k:28s} size {A[k].size:6d}  max rel diff {r:.2e}"); shown += 1
    ok &= worst < TOL_REL
    print(f"largest relative difference {worst:.2e} (tolerance {TOL_REL:.0e})")
    print("MATCH" if ok else "MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    sys.exit(main(sys.argv[1], sys.argv[2]))
