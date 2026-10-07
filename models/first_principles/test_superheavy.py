"""Regression test (referee issue): Track A predict() must not crash for Z > 110 when an s electron
is removed (YS15 QED table stops at Z = 110; F_SE/F_U are now extrapolated quadratically from
Z = 108-110 and flagged as extrapolated), and dirac_point must stay finite for Z alpha >= |kappa|.
Run: py -3.13 models/first_principles/test_superheavy.py
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import relativity as rel  # noqa: E402
import zexp  # noqa: E402

for Z, N in [(120, 2), (119, 1), (111, 1), (118, 2), (118, 118), (120, 118)]:
    v = zexp.predict(Z, N)
    assert math.isfinite(v) and v > 0, (Z, N, v)
    print(f"zexp.predict({Z},{N}) = {v:.3f} eV")
assert math.isfinite(rel.dirac_point(1, -1, 140))
assert abs(zexp.predict(92, 1) - 131819.128) < 0.01   # unchanged inside the table
print("OK")
