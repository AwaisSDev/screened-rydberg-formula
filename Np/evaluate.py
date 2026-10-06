"""Common scoring harness for every ionization-energy model.

A model is scored from a predictions CSV with columns  Z,N,IE_pred_eV  (extra columns ignored,
rows may be a subset of the dataset; missing rows are reported as coverage < 100%).

    py -3.13 evaluate.py results/<model>_predictions.csv [--json results/<model>_metrics.json]

Metrics are reported overall and per stratum, because a single global MAE in eV is dominated
by K-shell IEs of super-heavy ions (~10^5 eV). The headline numbers are relative errors.
"""
import argparse
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common.atomdata import load_records, L_LETTER  # noqa: E402


def _stats(pairs):
    """pairs: list of (pred, ref). Returns dict of error statistics."""
    if not pairs:
        return {"n": 0}
    rel = sorted(abs(p - r) / r * 100 for p, r in pairs)
    ab = sorted(abs(p - r) for p, r in pairs)
    n = len(pairs)
    med = lambda v: v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])  # noqa: E731
    return {
        "n": n,
        "MAPE_%": sum(rel) / n,
        "median_APE_%": med(rel),
        "p90_APE_%": rel[min(n - 1, int(0.9 * n))],
        "max_APE_%": rel[-1],
        "MAE_eV": sum(ab) / n,
        "median_AE_eV": med(ab),
        "RMSE_eV": math.sqrt(sum(a * a for a in ab) / n),
        "frac_within_1%": sum(x <= 1 for x in rel) / n,
        "frac_within_5%": sum(x <= 5 for x in rel) / n,
    }


def strata(rec):
    """Yield the stratum names a record belongs to."""
    Z, N = rec["Z"], rec["N"]
    yield "ALL"
    yield f"status={rec['status']}"
    if N == Z:
        yield "first_IE_neutral_atoms"
    if N == 1:
        yield "hydrogen_like"
    if rec["removed"]:
        n, l = rec["removed"]
        yield f"removed_l={L_LETTER[l]}"
        yield f"removed_n={min(n, 7)}"
    if rec["rearranged"]:
        yield "rearranged_config"
    yield "Z<=18" if Z <= 18 else ("19<=Z<=54" if Z <= 54 else "Z>=55")
    yield "N<=10" if N <= 10 else ("11<=N<=36" if N <= 36 else "N>=37")


def evaluate(pred_csv, subset=None):
    preds = {}
    with open(pred_csv, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                v = float(r["IE_pred_eV"])
            except (ValueError, TypeError):
                continue
            if math.isfinite(v):
                preds[(int(r["Z"]), int(r["N"]))] = v
    groups = {}
    recs = load_records()
    if subset is not None:
        recs = [r for r in recs if subset(r)]
    for rec in recs:
        k = (rec["Z"], rec["N"])
        if k not in preds:
            continue
        for s in strata(rec):
            groups.setdefault(s, []).append((preds[k], rec["IE_eV"]))
    out = {s: _stats(v) for s, v in sorted(groups.items())}
    out["_coverage"] = {"predicted": sum(1 for r in recs if (r["Z"], r["N"]) in preds),
                        "dataset": len(recs)}
    return out


def print_table(metrics, title=""):
    if title:
        print(f"\n=== {title} ===")
    cov = metrics.get("_coverage", {})
    print(f"coverage: {cov.get('predicted')}/{cov.get('dataset')}")
    print(f"{'stratum':32s} {'n':>5s} {'MAPE%':>8s} {'medAPE%':>8s} {'p90%':>8s} {'max%':>8s} "
          f"{'MAE eV':>10s} {'<=1%':>6s}")
    for s, m in metrics.items():
        if s.startswith("_") or not m.get("n"):
            continue
        print(f"{s:32s} {m['n']:5d} {m['MAPE_%']:8.3f} {m['median_APE_%']:8.3f} {m['p90_APE_%']:8.3f} "
              f"{m['max_APE_%']:8.2f} {m['MAE_eV']:10.3f} {m['frac_within_1%']:6.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pred_csv")
    ap.add_argument("--json", help="write metrics JSON here")
    a = ap.parse_args()
    m = evaluate(a.pred_csv)
    print_table(m, os.path.basename(a.pred_csv))
    if a.json:
        with open(a.json, "w", encoding="utf-8") as f:
            json.dump(m, f, indent=1)


if __name__ == "__main__":
    main()
