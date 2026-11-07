"""Limits pytest to the two test suites.

Root pytest configuration: `pytest` from the repository root collects the test suites
(melt_curve/validation/tests/ and tests/) only.

- melt_curve/validation/stats_test.py is a producer script whose name matches pytest's *_test.py
  pattern. Collecting it imports it, which runs it and rewrites melt_curve/validation/out/test_stats.json.
- production_campaign/ and studies/ hold the ARC scripts verbatim, as provenance; several parse
  sys.argv or start GPU work at import.
"""
import glob, os

_here = os.path.dirname(os.path.abspath(__file__))
collect_ignore = ["production_campaign", "studies", "data", "docs", "models",
                  os.path.join("melt_curve", "figures"), os.path.join("melt_curve", "layout.py")]
collect_ignore += [os.path.relpath(f, _here) for f in glob.glob(os.path.join(_here, "melt_curve", "validation", "*.py"))]
