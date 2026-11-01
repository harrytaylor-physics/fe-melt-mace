#!/bin/bash
#SBATCH --job-name=fe_stage
#SBATCH --account=phys-plasma-dft
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=fe_stage_%j.log
# stage_fe_melt.sh - read-only snapshot of ~/fe_melt for building the fe-melt-mace repository.
# Works with or without git: if ~/fe_melt is not a git repository, every file is snapshotted
# (subject to the exclusions) and the history files are not written.
#
# Run as a batch job from the home directory:  sbatch stage_fe_melt.sh   (log: ~/fe_stage_<jobid>.log)
# or directly:  bash stage_fe_melt.sh   (prints a timestamped line per step; do not pipe through tail)
# Writes only to $OUT (default $DATA/fe_melt_stage). Nothing in ~/fe_melt is moved, changed or
# committed; git is called with --no-optional-locks so it does not even refresh the index.
# Jobs writing into ~/fe_melt/labels can keep running.
#
# $OUT contents:
#   nested_git_repos.txt      git repositories found inside ~/fe_melt (depth <= 4)
#   fe_melt.bundle            full git history (all refs) of ~/fe_melt, if it is a git repository
#   git_status.txt            untracked (??) and modified ( M) files in the working tree
#   git_log.txt, git_size.txt history listing and object-store size
#   history_large_blobs.txt   blobs over 5 MB anywhere in history
#   worktree.tgz              working-tree files not in git history (all files if no git), with the exclusions below
#   worktree_manifest.tsv     every working-tree file: path, bytes, sha256, included|excluded:<reason>
#   provenance/               model files, per-tensor comparison, training-log lines, model paths in scripts,
#                             tel-mace commit and remote, tel_env packages, QE version, pseudopotential hashes
#   credential_scan.txt       files matching credential-like patterns; excluded from worktree.tgz
#   SUMMARY.txt
# Exclusions from worktree.tgz (hashed in the manifest instead): QE raw outputs (*.out under labels/),
# QE scratch (*.save/, tmp*/ , wfc*), trajectories (*.traj *.dcd *.lammpstrj *.nc), files over 200 MB,
# files listed in credential_scan.txt.
set -u
SRC=${SRC:-$HOME/fe_melt}
DATA=${DATA:?DATA is not set: ARC sets it to the project data directory; or export DATA=<dir>}
OUT=${OUT:-$DATA/fe_melt_stage}
TELSRC=${TELSRC:-$HOME/tel-mace-src}
TELENV=${TELENV:-$HOME/tel_env}
PY=${PY:-python3}
if git --no-optional-locks version >/dev/null 2>&1; then NOL="--no-optional-locks"; else NOL=""; fi   # git >= 2.15
G="git $NOL -C $SRC"

log() { echo "[$(date +%H:%M:%S)] $*"; }
[ -d "$SRC" ] || { echo "no directory $SRC"; exit 1; }
mkdir -p "$OUT/provenance" || exit 1
log "source $SRC -> $OUT ($(du -sh "$SRC" 2>/dev/null | cut -f1))"
find "$SRC" -mindepth 2 -maxdepth 4 -name .git 2>/dev/null > "$OUT/nested_git_repos.txt"
if [ ! -d "$SRC/.git" ]; then
  # no version control: every file goes into the snapshot (subject to the exclusions), no history
  log "no git repository at $SRC: snapshot without history"
  ( cd "$SRC" && find . -type f | sed 's#^\./#?? #' ) > "$OUT/git_status.txt"
  echo "none: $SRC is not a git repository" > "$OUT/git_log.txt"
else

# 1. git history
$G status --porcelain=v1 -uall > "$OUT/git_status.txt"
$G log --all --date=iso --format='%H %ad %an %s' > "$OUT/git_log.txt"
$G count-objects -vH > "$OUT/git_size.txt"
$G for-each-ref --format='%(refname) %(objectname)' >> "$OUT/git_size.txt"
$G rev-list --objects --all | $G cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' \
  | awk '$1=="blob" && $3>5000000 {print $3, $2, $4}' | sort -n -r > "$OUT/history_large_blobs.txt"
$G bundle create "$OUT/fe_melt.bundle" --all 2> "$OUT/bundle.log" && git bundle verify "$OUT/fe_melt.bundle" >> "$OUT/bundle.log" 2>&1
log "history: $(wc -l < "$OUT/git_log.txt") commits, bundle $(du -h "$OUT/fe_melt.bundle" | cut -f1), $(wc -l < "$OUT/history_large_blobs.txt") blobs > 5 MB"
fi

# 2. credential-like files in the working tree (labels/ and data files excluded); kept out of worktree.tgz
log "credential scan"
timeout 900 grep -rIl --exclude=*.xyz --exclude=*.extxyz --exclude=*.npz --exclude=*.npy --exclude=*.model \
  --exclude=*.pt --exclude=*.UPF --exclude=*.out -E "BEGIN [A-Z ]*PRIVATE KEY|password *[:=]|api[_-]?key *[:=]|secret *[:=]|token *[:=]" "$SRC" \
  --exclude-dir=.git --exclude-dir=labels 2>/dev/null > "$OUT/credential_scan.txt"
find "$SRC" -path "$SRC/.git" -prune -o \( -name ".env" -o -name "id_rsa*" -o -name "*.pem" -o -iname "*token*" -o -iname "*secret*" \) -print >> "$OUT/credential_scan.txt"

# 3. working tree: untracked and modified files, filtered, with a full manifest
log "manifest and worktree.tgz (progress every 2,000 files)"
SRC="$SRC" OUT="$OUT" $PY -u - <<'EOF'
import hashlib, os, re, tarfile, time
src, out = os.environ["SRC"], os.environ["OUT"]
st = [l.rstrip("\n") for l in open(os.path.join(out, "git_status.txt"))]
pend = {l[3:].strip('"') for l in st if l[:2] in ("??", " M", "M ", "MM", "A ", "AM")}
cred = {os.path.relpath(l.strip(), src) for l in open(os.path.join(out, "credential_scan.txt")) if l.strip()}
def reason(rel, size):
    p = rel.split("/")
    if rel in cred: return "credential-like"
    if rel.startswith("labels/") and rel.endswith(".out"): return "QE raw output"
    if any(x.endswith(".save") or x.startswith("tmp") for x in p[:-1]) or re.match(r"wfc\d*", p[-1]): return "QE scratch"
    if re.search(r"\.(traj|dcd|lammpstrj|nc)$", rel): return "trajectory"
    if "__pycache__" in p or rel.endswith(".pyc"): return "bytecode"
    if size > 200e6: return "over 200 MB"
    return None
def sha(f):
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
n_in = n_ex = b_in = b_ex = n = 0
t0 = time.time()
with open(os.path.join(out, "worktree_manifest.tsv"), "w") as man, \
     tarfile.open(os.path.join(out, "worktree.tgz"), "w:gz", compresslevel=1) as tar:
    man.write("path\tbytes\tsha256\tgit\tstatus\n")
    for d, dirs, files in os.walk(src):
        dirs[:] = [x for x in dirs if x != ".git"]
        for fn in files:
            f = os.path.join(d, fn); rel = os.path.relpath(f, src)
            if os.path.islink(f) or not os.path.isfile(f): continue
            size = os.path.getsize(f)
            n += 1
            if n % 2000 == 0:
                print(f"  {n} files, {b_in / 1e6:.0f} MB in, {b_ex / 1e6:.0f} MB excluded, {time.time() - t0:.0f} s", flush=True)
            tracked = "snapshot" if rel in pend else "tracked"
            r = reason(rel, size)
            man.write(f"{rel}\t{size}\t{sha(f)}\t{tracked}\t{'excluded:' + r if r else 'ok'}\n")
            if tracked == "tracked": continue
            if r: n_ex += 1; b_ex += size; continue
            tar.add(f, arcname=rel); n_in += 1; b_in += size
print(f"  done in {time.time() - t0:.0f} s: {n_in} snapshot files ({b_in / 1e6:.0f} MB) in worktree.tgz; "
      f"{n_ex} excluded ({b_ex / 1e6:.0f} MB), hashed in worktree_manifest.tsv")
EOF

# 4. model provenance (read-only)
P="$OUT/provenance"
log "model files (search of $DATA limited to depth 4 and 10 min)"
( find "$SRC" -name "*.model" -not -path "$SRC/.git/*"
  timeout 600 find "$DATA" -maxdepth 4 -name "fe_*.model" -not -path "*/fe_melt_stage/*" || echo "# search of $DATA stopped after 10 min" >&2
) 2>"$P/model_search.err" | sort -u > "$P/model_files.txt"
SRC="$SRC" P="$P" $PY -u - <<'EOF'
import hashlib, os, time, zipfile
P = os.environ["P"]; src = os.environ["SRC"]
fs = [l.strip() for l in open(os.path.join(P, "model_files.txt")) if l.strip()]
def tensors(f):
    z = zipfile.ZipFile(f); names = z.namelist()
    return names[0].split("/")[0], {n.split("/", 1)[1]: hashlib.sha256(z.read(n)).hexdigest() for n in names if "/data/" in n}
R = {}
with open(os.path.join(P, "model_files.tsv"), "w") as o:
    o.write("path\tbytes\tmtime\tsha256\tsaved_as\tn_tensors\n")
    for f in fs:
        h = hashlib.sha256(open(f, "rb").read()).hexdigest()
        root, t = tensors(f) if zipfile.is_zipfile(f) else ("not-a-zip", {})
        R[f] = t
        o.write(f"{f}\t{os.path.getsize(f)}\t{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(os.path.getmtime(f)))}\t{h}\t{root}\t{len(t)}\n")
with open(os.path.join(P, "model_tensor_compare.txt"), "w") as o:
    for name in sorted({os.path.basename(f) for f in fs if "fe_" in os.path.basename(f)}):
        same = [f for f in fs if os.path.basename(f) == name]
        o.write(f"{name}: {len(same)} copies, distinct sha256: {len({hashlib.sha256(open(f,'rb').read()).hexdigest() for f in same})}\n")
    # every model in phaseD/models is a production candidate: compare each with every other file
    refs = sorted(f for f in fs if os.path.dirname(f) == os.path.join(src, "phaseD", "models"))
    for r in refs:
        a = R[r]
        o.write(f"\nidentical tensors with {r} ({len(a)} tensors):\n")
        for f in fs:
            b = R[f]
            if b and len(b) == len(a):
                o.write(f"  {sum(1 for k in a if b.get(k) == a[k]):4d}/{len(a)}  {f}\n")
print(f"provenance: {len(fs)} .model files listed in provenance/model_files.tsv")
EOF
log "training logs, model paths, run records"
grep -n -i -E "saving model|stage ?two|swa|KeyError|best model|loading checkpoint" "$SRC"/logs/*.out "$SRC"/logs/*.log 2>/dev/null > "$P/training_log_lines.txt"
grep -rn --exclude-dir=.git --exclude-dir=labels --include=*.py --include=*.sbatch --include=*.sh --include=*.yaml --include=*.yml \
  -E "\.model\b|model_paths|--model" "$SRC" 2>/dev/null > "$P/model_paths_in_scripts.txt"
for d in "$SRC"/phaseE/runs*; do
  [ -d "$d" ] || continue
  echo "== $d"; grep -h -o -E '"model[a-z_]*": *"[^"]*"' "$d"/*.json 2>/dev/null | sort | uniq -c
done > "$P/model_paths_in_run_records.txt"

# 5. tel-mace fork and environment
log "tel-mace and environment"
{
  echo "== $TELSRC"
  git $NOL -C "$TELSRC" log -1 --date=iso --format='HEAD %H %ad %s' 2>&1
  git $NOL -C "$TELSRC" describe --tags --always 2>&1
  git $NOL -C "$TELSRC" remote -v 2>&1
  echo "-- uncommitted changes:"; git $NOL -C "$TELSRC" status --porcelain 2>&1 | head -50
  echo "-- branches:"; git $NOL -C "$TELSRC" branch -a -vv 2>&1 | head -20
} > "$P/tel_mace_git.txt"
grep -rn -E "elec_temp|1e4|10000|10_000" "$TELSRC/mace" --include=*.py 2>/dev/null | head -80 > "$P/tel_mace_elec_temp_lines.txt"
if [ -f "$TELENV/bin/activate" ]; then
  ( source "$TELENV/bin/activate" && timeout 300 python -c "import mace, torch, ase, sys; print('python', sys.version.split()[0]); print('mace', getattr(mace, '__version__', '?'), mace.__file__); print('torch', torch.__version__); print('ase', ase.__version__)" && pip freeze ) > "$P/tel_env.txt" 2>&1
fi
f=$(ls "$SRC"/labels/*.out 2>/dev/null | head -1); [ -n "$f" ] && grep -m3 -E "Program PWSCF|Parallel version|Exchange-correlation" "$f" > "$P/qe_version.txt"
sha256sum "$SRC"/pseudo/*.UPF > "$P/pseudopotentials.sha256" 2>/dev/null
ls "$SRC"/labels/*.sbatch 2>/dev/null | head -3 | xargs -r grep -h -E "module|#SBATCH" > "$P/qe_job_header.txt"

# 6. summary
log "summary"
{
  echo "staged $(date -Iseconds) from $SRC (read-only)"
  du -sh "$SRC"/* 2>/dev/null | sort -h | tail -15
  echo; echo "stage: $(du -sh "$OUT" | cut -f1)"; ls -la "$OUT"
  echo; echo "credential scan: $(wc -l < "$OUT/credential_scan.txt") file(s)"
  echo; cat "$P/model_tensor_compare.txt"
  echo; echo "tel-mace:"; head -3 "$P/tel_mace_git.txt"
} > "$OUT/SUMMARY.txt"
cat "$OUT/SUMMARY.txt"
