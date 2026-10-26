"""Root pytest configuration: `pytest` from the repository root collects the test suites
(validation/tests/ and tests/) only.

- validation/stats_test.py is a producer script whose name matches pytest's *_test.py pattern.
  Collecting it imports it, which runs it and rewrites validation/out/test_stats.json.
- The ARC scripts at the root and in labels/, phaseD/, phaseE/ and exploratory/ are kept
  verbatim as provenance; several parse sys.argv or start GPU work at import.
"""
import glob, os

_here = os.path.dirname(os.path.abspath(__file__))
collect_ignore = ["labels", "logs", "phaseD", "phaseE", "phaseX", "exploratory", "traj",
                  "pseudo", "results", "analysis", "tools", "docs", "data", "models"]
collect_ignore += [os.path.basename(f) for f in glob.glob(os.path.join(_here, "*.py"))
                   if os.path.basename(f) != "conftest.py"]
collect_ignore += [os.path.relpath(f, _here) for f in glob.glob(os.path.join(_here, "validation", "*.py"))]
