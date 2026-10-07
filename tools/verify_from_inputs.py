r"""Reference implementation of the Screened Rydberg formula from the paper equation + results/model_inputs.csv.

    py -3.11 tools/verify_from_inputs.py

It imports nothing from models/. Its only inputs are:
  results/model_inputs.csv   per-row inputs (tools/export_model_inputs.py)
  results/pa_params.json     fitted parameters (Table A1 at full precision; pa_hier_rel and pa_bound9)
  CODATA 2018: Ry = 13.605693122994 eV, alpha = 7.2973525693e-3
It evaluates paper eq. (2.2) for every row and compares with the production predictions
(results/pa_hier_rel_predictions.csv, results/pa_bound9_predictions.csv). Exit code 1 if any row differs by more
than 1e-9 relative.

Equation (paper section 2.2):
  Za = max(Z - N + 1, 1)
  T  = sum_g tau_g nu_g + sum_c dtau_c nu_c                      (class deviations: pa_hier_rel only)
  h  = (N - 1) - sigma1 if T >= 0 else sigma1
  D  = T / (Za + kappa_used + |T|/h)     with kappa_used = |kappa| + 0.05; D = 0 when T = 0 (e.g. N = 1)
  Zeff = Z - sigma1 - D
  F    = Dirac/Schroedinger ratio for a point charge Zeff, level (n, j):
         E_D = 1 - [1 + (x/(n - k + sqrt(k^2 - x^2)))^2]^(-1/2),  x = Zeff alpha, k = j + 1/2,  F = 2 n^2 E_D / x^2
  R    = 1 + r_c (Z alpha)^2 / n * (Zeff/Za - 1)                  (pa_hier_rel only; r_c by relativistic class)
  IE   = mu * [Ry Zeff^2/n^2 * F * R + Ry x_l K / n^2] - qedfns * (Zeff/Z)^2
"""
import csv
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RY, ALPHA = 13.605693122994, 7.2973525693e-3
GROUPS = ["same", "in", "core", "df", "out"]


def dirac_ratio(n, j, Ze):
    x = max(Ze * ALPHA, 1e-8)
    k = j + 0.5
    x = min(x, 0.999 * k)
    nr = n - k + math.sqrt(k * k - x * x)
    E = 1.0 - 1.0 / math.sqrt(1.0 + (x / nr) ** 2)
    return E * 2.0 * n * n / (x * x)


def ie_row(r, P, use_dev, use_rel):
    Z, N, n = int(r["Z"]), int(r["N"]), int(r["n"])
    l, j = int(r["l"]), float(r["j"])
    s1 = float(r["sigma1"])
    Za = max(Z - N + 1, 1)
    T = sum(P["tau_" + g] * int(r["nu_" + g]) for g in GROUPS)
    if use_dev:
        T += sum(P[k] * int(r["nuc_" + k[5:]]) for k in P if k.startswith("dtau_"))
    h = max(N - 1 - s1, 1e-9) if T >= 0 else max(s1, 1e-9)
    kap = abs(P["kappa"]) + 0.05
    D = T / (Za + kap + abs(T) / h)
    Ze = max(Z - s1 - D, 1e-3)
    b = RY * Ze ** 2 / n ** 2 * dirac_ratio(n, j, Ze)
    if use_rel:
        b *= 1.0 + P["rel_" + r["rel_class"]] * (Z * ALPHA) ** 2 / n * (Ze / Za - 1.0)
    xl = [0.0, P["x_p"], P["x_d"], P["x_f"]][l]
    b += RY * xl * float(r["K"]) / n ** 2
    return b * float(r["mu"]) - float(r["qedfns_eV"]) * (Ze / Z) ** 2


def main():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "results", "model_inputs.csv"), encoding="utf-8")))
    cands = json.load(open(os.path.join(ROOT, "results", "pa_params.json"), encoding="utf-8"))["candidates"]
    bad = 0
    for name, use_dev, use_rel in (("pa_hier_rel", True, True), ("pa_bound9", False, False)):
        c = cands[name]
        P = dict(zip(c["names"], c["values"]))
        ref = {(int(q["Z"]), int(q["N"])): float(q["IE_pred_eV"]) for q in
               csv.DictReader(open(os.path.join(ROOT, "results", f"{name}_predictions.csv"), encoding="utf-8"))}
        worst, wk = 0.0, None
        for r in rows:
            v = ie_row(r, P, use_dev, use_rel)
            key = (int(r["Z"]), int(r["N"]))
            d = abs(v / ref[key] - 1)
            if d > worst:
                worst, wk = d, key
        ok = worst < 1e-9
        bad += not ok
        print(f"{name:12s} {len(rows)} rows  max |rel diff| vs production = {worst:.2e} at {wk}  "
              f"{'PASS' if ok else 'FAIL'}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
