"""Repository-wide checks: no absolute machine paths in committed code; the report-number list
matches the ledger; the layout module's paths exist."""
import csv, json, os, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "melt_curve"))
import layout  # noqa: E402

ABS = re.compile(r"/(Users|home|tmp|data|scratch)/[\w.-]+")
CODE = (".py", ".sh", ".sbatch")


def tracked():
    try:
        out = subprocess.run(["git", "ls-files", "-co", "--exclude-standard"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
        return out.split()
    except (OSError, subprocess.CalledProcessError):
        return [os.path.relpath(os.path.join(d, f), ROOT) for d, _, fs in os.walk(ROOT) if ".git" not in d for f in fs]


def test_no_absolute_paths_in_code():
    hits = []
    for f in tracked():
        if f.endswith(CODE):
            for i, l in enumerate(open(os.path.join(ROOT, f), errors="ignore"), 1):
                if ABS.search(l):
                    hits.append(f"{f}:{i}: {l.strip()[:100]}")
    assert not hits, "\n".join(hits)


def lookup(L, key):
    idx = None
    if key.endswith("]"):
        key, i = key[:-1].rsplit("[", 1)
        idx = int(i)
    if key.startswith("keys:"):
        v = L["keys"][key[5:]]
    else:
        v = L["summary"]
        for p in key[8:].split("/"):
            v = v[p]
    return v[idx] if idx is not None else v


def fmt(v):
    return f"{v:.10g}" if isinstance(v, float) else json.dumps(v)


def test_report_numbers_match_ledger():
    L = json.load(open(layout.LEDGER))
    rows = list(csv.DictReader(open(layout.REPORT_NUMBERS), delimiter="\t", quoting=csv.QUOTE_NONE))
    assert rows and list(rows[0]) == ["number_as_printed", "unit", "key", "ledger_value"]
    for r in rows:
        if r["key"]:
            assert fmt(lookup(L, r["key"])) == r["ledger_value"], r
        else:
            assert r["ledger_value"] == "", r


def test_layout_paths_exist():
    for name in ("MODEL", "MODEL_STAGE_ONE", "TRAIN_SET", "PARITY", "TRAIN_LOG", "TEL_CHECK_LOG", "DATASET_REPORT",
                 "ANALYSIS_TEL", "ANALYSIS_T6000", "PN_CORRECT", "DECOMP_TEL", "LITERATURE", "LEDGER", "SPOT_FRAMES"):
        assert os.path.exists(getattr(layout, name)), name
    for p in list(layout.TEST_SETS.values()) + list(layout.MODEL_V1.values()):
        assert os.path.exists(p), p
