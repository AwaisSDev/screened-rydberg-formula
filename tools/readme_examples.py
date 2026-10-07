"""README example table (five ions, predicted vs NIST) from the final model.   py -3.11 tools/readme_examples.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ionization import ie  # noqa: E402

CASES = [("O", 8, "O I"), ("Na", 11, "Na I"), ("Fe", 26, "Fe I"), ("Fe", 10, "Fe XVII"), ("Pb", 82, "Pb I")]
print("| ion | removed | Z_eff | predicted (eV) | NIST (eV) | error |")
print("|---|---|---|---|---|---|")
for sym, N, lab in CASES:
    d = ie(sym, N)
    print(f"| {lab} | {d['removed']} | {d['Zeff']:.4f} | {d['IE_eV']:.3f} | {d['NIST_eV']:.3f} | {d['error_pct']:+.2f} % |")
