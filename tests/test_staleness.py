"""Staleness guard: every committed ledger input that no producer here rebuilds is recorded in
validation/ledger_inputs.tsv with its sha256 and date. Fails when an input's content is not the one
recorded (it changed after the output was built: rebuild the output, then run
validation/record_inputs.py) or when a recorded input is newer than the output built from it."""
import csv, hashlib, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TSV = os.path.join(ROOT, "validation", "ledger_inputs.tsv")
ROWS = list(csv.DictReader(open(TSV), delimiter="\t"))


def sha(p):
    f = os.path.join(ROOT, p)
    data = open(f, "rb").read()
    if data[:7] == b"version":                     # LFS pointer: compare the recorded object id instead
        return data.split(b"sha256:")[1].split()[0].decode()
    return hashlib.sha256(data).hexdigest()


def test_inputs_unchanged_since_outputs_were_built():
    bad = [r["input"] for r in ROWS if sha(r["input"]) != r["input_sha256"]]
    assert not bad, f"inputs changed after their outputs were built (rebuild, then record_inputs.py): {bad[:10]}"


def test_outputs_are_the_recorded_ones():
    bad = sorted({r["output"] for r in ROWS if sha(r["output"]) != r["output_sha256"]})
    assert not bad, f"outputs differ from the record (run validation/record_inputs.py): {bad}"


def test_no_input_newer_than_its_output():
    bad = [(r["input"], r["input_date"], r["output"], r["output_date"]) for r in ROWS if r["input_date"] > r["output_date"]]
    assert not bad, f"input newer than the output built from it: {bad[:10]}"


def test_every_validation_output_is_covered():
    import glob
    outs = {os.path.relpath(f, ROOT) for f in glob.glob(os.path.join(ROOT, "validation", "out", "*.json"))}
    assert outs <= {r["output"] for r in ROWS}
