"""conv_check.py - basis-set and k-point convergence of the static PBE pressure.

hcp Fe at 7.0 A^3/atom, c/a 1.60, Marzari-Vanderbilt 0.01 Ry, QE 6.6, pslibrary 1.0.0 PBE PAW
(refs/convergence/*.json, inputs alongside). Differences are relative to the production setting
used for the static references (ecutwfc 100 Ry, ecutrho 1,000 Ry) at k 16x16x10.
Output: out/convergence.json.
"""
import glob, json, os
import paths

REF = "e100k16"


def main():
    d = os.path.join(os.path.dirname(paths.DFT_REFS), "convergence")
    R = {r["tag"]: r for r in (json.load(open(f)) for f in glob.glob(os.path.join(d, "*.json")))}
    ref = R[REF]
    rows = []
    for t in sorted(R, key=lambda t: (R[t]["ecutwfc"], R[t]["kmesh"])):
        r = R[t]
        rows.append(dict(tag=t, ecutwfc=r["ecutwfc"], ecutrho=r["ecutrho"], kmesh=r["kmesh"], P=r["P"],
                         dP=r["P"] - ref["P"], dE_meV=1e3 * (r["E_per_atom"] - ref["E_per_atom"])))
        print(f"{t:8s} ecut {r['ecutwfc']:5.0f}/{r['ecutrho']:6.0f} Ry  k {r['kmesh']}  P {r['P']:7.2f} GPa  "
              f"P - P(100/1000, k16) {r['P'] - ref['P']:+6.2f} GPa")
    out = dict(rows=rows, reference=REF,
               cutoff_dP_100_to_130=R["e130k16"]["P"] - ref["P"],
               kpoint_dP_12_to_16=ref["P"] - R["e100k12"]["P"])
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "convergence.json"), "w"), indent=1)
    print(f"cutoff 100/1,000 -> 130/1,300 Ry: {out['cutoff_dP_100_to_130']:+.2f} GPa; "
          f"k 12x12x8 -> 16x16x10: {out['kpoint_dP_12_to_16']:+.2f} GPa")


if __name__ == "__main__":
    main()
