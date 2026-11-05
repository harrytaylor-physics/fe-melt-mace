"""train_table.py - architecture, data and final errors of the production potential.

Inputs (paths from ../layout.py):
  production_campaign/arc_scripts/logs/fe_tel_e3_run-7.log         MACE training log: model details, split sizes, final error tables
  production_campaign/arc_scripts/logs/fe_tel_e3_run-7_debug.log   the run's full argument namespace (dtype, embedding spec, weights)
  data/training/dataset_report.txt        frames per T_el (train / held-out test) and per phase
  data/training/train_e3.xyz              the training file: frames per T_el and per phase (LFS; header lines only)
  models/fe_tel_e3_stagetwo.model   model class name, read from the archive (LFS; optional)
Output: out/training.json.
"""
import ast, json, os, re, sys, zipfile
from collections import Counter
import paths

sys.path.insert(0, os.path.dirname(paths.HERE))
import layout  # noqa: E402

TEL = (4000, 6000, 8000)


def log_details(path):
    txt = open(path).read()
    d = {}
    d["r_max_A"] = float(re.findall(r"Radial cutoff: ([\d.]+) A", txt)[-1])
    d["hidden_irreps"] = re.findall(r"Message passing with hidden irreps (\S+?)\)?$", txt, re.M)[-1]
    d["channels"] = int(re.match(r"(\d+)x0e", d["hidden_irreps"]).group(1))
    m = re.findall(r"(\d+) layers, each with correlation order: (\d+) \(body order: (\d+)\) and spherical harmonics up to: l=(\d+)", txt)[-1]
    d["interactions"], d["correlation"], d["body_order"], d["max_ell"] = map(int, m)
    n = re.findall(r"Total number of configurations: train=(\d+), valid=(\d+)", txt)[-1]
    d["n_train_split"], d["n_valid_split"] = int(n[0]), int(n[1])
    d["n_train_file"] = int(re.findall(r"Total Training set \[energy: (\d+)", txt)[-1])
    d["foundation_model"] = os.path.basename(re.findall(r"Using foundation model (\S+) as initial checkpoint", txt)[-1])
    d["stage_two_after_epochs"] = int(re.findall(r"Stage Two \(after (\d+) epochs\)", txt)[-1])
    d["avg_num_neighbors"] = float(re.findall(r"Unusual average number of neighbors: ([\d.]+)", txt)[-1])
    # the last two error tables: stage one (best checkpoint), then stage two (SWA); the production file is stage two
    tabs = re.findall(r"\| train_Default \|\s+([\d.]+)\s+\|\s+([\d.]+)\s+\|\s+([\d.]+)\s+\|\s*\n"
                      r"\| valid_Default \|\s+([\d.]+)\s+\|\s+([\d.]+)\s+\|\s+([\d.]+)\s+\|", txt)
    stage = {}
    for name, t in zip(("stage_one", "stage_two"), tabs[-2:]):
        t = list(map(float, t))
        stage[name] = dict(train_E_meV_atom=t[0], train_F_meV_A=t[1], train_F_rel_pct=t[2],
                           valid_E_meV_atom=t[3], valid_F_meV_A=t[4], valid_F_rel_pct=t[5])
    d["final_errors"] = stage
    d["stage_one_epoch"] = int(re.findall(r"Loaded Stage one model from epoch (\d+)", txt)[-1])
    d["stage_two_epoch"] = int(re.findall(r"Loaded Stage two model from epoch (\d+)", txt)[-1])
    return d


def debug_args(path):
    line = [l for l in open(path) if "Configuration: Namespace(" in l][-1]
    get = lambda k: re.search(rf"\b{k}=('[^']*'|[^,]+)", line).group(1).strip("'")
    spec = ast.literal_eval(re.search(r"embedding_specs=(\{.*?\}\}),", line).group(1))
    return dict(default_dtype=get("default_dtype"), max_num_epochs=int(get("max_num_epochs")),
                batch_size=int(get("batch_size")), seed=int(get("seed")), valid_fraction=float(get("valid_fraction")),
                energy_weight=float(get("energy_weight")), forces_weight=float(get("forces_weight")),
                swa_energy_weight=float(get("swa_energy_weight")), swa_forces_weight=float(get("swa_forces_weight")),
                elec_temp_embedding=spec["elec_temp"])


def dataset_report(path):
    out = {}
    for l in open(path):
        m = re.match(r"Tel (\d+): (\d+) train / (\d+) test\s+by phase: (\{.*\})", l.strip())
        if m:
            out[m.group(1)] = dict(train=int(m.group(2)), test=int(m.group(3)), by_phase=ast.literal_eval(m.group(4)))
    return out


def train_file(path):
    per_tel, per_phase, natoms = Counter(), Counter(), Counter()
    with open(path) as f:
        while True:
            l = f.readline()
            if not l:
                break
            n = int(l)
            h = f.readline()
            per_tel[str(int(float(re.search(r"elec_temp=(\S+)", h).group(1))))] += 1
            per_phase[re.search(r"phase=(\S+)", h).group(1)] += 1
            natoms[n] += 1
            for _ in range(n):
                f.readline()
    tot = sum(per_phase.values())
    return dict(n_frames=tot, per_tel=dict(per_tel), per_phase=dict(per_phase),
                phase_pct={k: 100.0 * v / tot for k, v in per_phase.items()}, atoms_per_frame=sorted(natoms))


def model_class(path):
    with zipfile.ZipFile(path) as z:
        pkl = z.read([n for n in z.namelist() if n.endswith("data.pkl")][0])
    m = re.search(rb"mace\.modules\.models\n(\w+)\n", pkl)
    return m.group(1).decode() if m else None


def lfs_present(path):
    return os.path.exists(path) and open(path, "rb").read(40)[:7] != b"version"


def main():
    out = dict(model=log_details(layout.TRAIN_LOG), args=debug_args(layout.TRAIN_DEBUG_LOG),
               dataset_report=dataset_report(layout.DATASET_REPORT))
    dr = out["dataset_report"]
    out["T_el_K"] = sorted(int(T) for T in dr)
    out["test_frames"] = {T: dr[T]["test"] for T in dr}
    out["train_frames_range"] = [min(dr[T]["train"] for T in dr), max(dr[T]["train"] for T in dr)]
    if lfs_present(layout.TRAIN_SET):
        out["train_file"] = train_file(layout.TRAIN_SET)
        assert out["train_file"]["n_frames"] == out["model"]["n_train_file"]
        assert {int(k): v for k, v in out["train_file"]["per_tel"].items()} == {int(T): dr[T]["train"] for T in dr}
    else:
        sys.exit(f"{layout.TRAIN_SET} is not present (git lfs pull)")
    out["model"]["class"] = model_class(layout.MODEL) if lfs_present(layout.MODEL) else None
    for k, v in sorted(out.items()):
        print(k, json.dumps(v)[:300])
    os.makedirs(paths.OUT, exist_ok=True)
    json.dump(out, open(os.path.join(paths.OUT, "training.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
