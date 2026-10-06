"""results/fp_zexp_coefficients.csv: exact E0, E1, dE0, dE1 and the ab initio screening constant
for every N = 1..110, for two configurations of each isoelectronic sequence:
  'neutral' : NIST ground configuration of the neutral atom Z = N (and of its cation)
  'highZ'   : NIST ground configuration of the most highly charged member in the table (largest Z)
Exact rationals are written when the denominator has < 25 digits; otherwise a float.
"""
import csv
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zexp  # noqa: E402
from zexp import ROOT  # noqa: E402
from build_coeffs import slater_sigma  # noqa: E402
from common.atomdata import load_records, ground_shells, shells_to_str, L_LETTER  # noqa: E402


def exact_str(s):
    if "/" in s or s.lstrip("-").isdigit():
        f = Fraction(s)
        if len(str(f.denominator)) < 25:
            return str(f)
    return ""


def row_for(Z, N, label):
    sh = ground_shells(Z, N)
    shi = ground_shells(Z, N - 1) if N > 1 else []
    a = zexp.config_data(sh)
    c = zexp.coeffs(Z, N, sh)
    return {
        "N": N, "which": label, "Z_ref": Z, "config": shells_to_str(sh), "config_ion": shells_to_str(shi) or "-",
        "removed": f"{c['n']}{L_LETTER[c['l']]}",
        "E0_exact": exact_str(a["E0"]), "E0": float(Fraction(a["E0"])),
        "E1_exact": exact_str(a["E1_sc"]), "E1_singleconfig_Hund": zexp._f(a["E1_sc"]),
        "E1_config_average": zexp._f(a["E_av"]), "E1_Layzer_complex": a["E1_cx"],
        "dE0": c["dE0"], "dE1": c["dE1"], "dE1_singleconfig": c["dE1_sc"],
        "sigma_abinitio": c["sigma"], "sigma_slater": slater_sigma(sh, c["n"], c["l"]),
    }


def main():
    R = load_records()
    maxZ = {}
    for r in R:
        maxZ[r["N"]] = max(maxZ.get(r["N"], 0), r["Z"])
    rows = []
    for N in range(1, 111):
        rows.append(row_for(N, N, "neutral"))
        if maxZ.get(N, N) != N:
            rows.append(row_for(maxZ[N], N, "highZ"))
    out = os.path.join(ROOT, "results", "fp_zexp_coefficients.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    zexp.save_cfg_cache()
    print("wrote", out, len(rows))


if __name__ == "__main__":
    main()
