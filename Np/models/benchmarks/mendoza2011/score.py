"""Score the published NRSHM (Mendoza et al. 2011) on the 5847 NIST rows. No refitting.

Writes results/bench_mendoza2011_predictions.csv (primary: total-energy difference, NIST ground configs,
j-split 'low') and variants, the evaluate.py metrics JSONs, and results/bench_mendoza2011_summary.json with
the MAPE on the S1/S2/S3 test rows. Rows whose configurations (N or N-1) need a subshell beyond 5p3/2 have no
published constants and are left out (coverage < 100 %).
"""
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
import model as M  # noqa: E402
from common.atomdata import load_records  # noqa: E402
from evaluate import evaluate  # noqa: E402

VARIANTS = {
    "bench_mendoza2011": dict(method="delta", split="low"),           # primary
    "bench_mendoza2011_stat": dict(method="delta", split="stat"),
    "bench_mendoza2011_xalpha": dict(method="xalpha", split="low"),   # one-electron binding energy
}


def mape(pairs):
    return sum(abs(p - r) / r for p, r in pairs) / len(pairs) * 100 if pairs else float("nan")


def main():
    recs = load_records()
    summary = {}
    for name, kw in VARIANTS.items():
        rows, skipped = [], 0
        for r in recs:
            try:
                v = M.predict(r["Z"], r["N"], **kw)
            except M.NotCovered:
                skipped += 1
                continue
            rows.append((r["Z"], r["N"], v))
        path = os.path.join(ROOT, "results", f"{name}_predictions.csv")
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["Z", "N", "IE_pred_eV"])
            for Z, N, v in rows:
                w.writerow([Z, N, f"{v:.6f}"])
        subprocess.run([sys.executable, os.path.join(ROOT, "evaluate.py"), path, "--json",
                        os.path.join(ROOT, "results", f"{name}_metrics.json")], check=True,
                       stdout=subprocess.DEVNULL)
        pred = {(Z, N): v for Z, N, v in rows}
        m = evaluate(path)

        def sub(sel):
            pr = [(pred[(r["Z"], r["N"])], r["IE_eV"]) for r in recs if sel(r) and (r["Z"], r["N"]) in pred]
            return {"n": len(pr), "n_total": sum(1 for r in recs if sel(r)), "MAPE_%": mape(pr)}
        s = {
            "settings": kw,
            "coverage": m["_coverage"],
            "skipped_beyond_5p": skipped,
            "covered_ALL": {k: m["ALL"][k] for k in ("n", "MAPE_%", "median_APE_%")},
            "neutral": {k: m.get("first_IE_neutral_atoms", {}).get(k) for k in ("n", "MAPE_%", "median_APE_%")},
            "N<=10": {k: m.get("N<=10", {}).get(k) for k in ("n", "MAPE_%")},
            "Z>=55": {k: m.get("Z>=55", {}).get(k) for k in ("n", "MAPE_%")},
            "d_removed": {k: m.get("removed_l=d", {}).get(k) for k in ("n", "MAPE_%")},
            "f_removed": {k: m.get("removed_l=f", {}).get(k) for k in ("n", "MAPE_%")},
            "experimental": {k: m.get("status=experimental", {}).get(k) for k in ("n", "MAPE_%")},
            "H_like": {k: m.get("hydrogen_like", {}).get(k) for k in ("n", "MAPE_%")},
            "ions_not_neutral": sub(lambda r: r["N"] < r["Z"]),
            "S1_test(Z%5==0)": sub(lambda r: r["Z"] % 5 == 0),
            "S2_test(Z>=55)": sub(lambda r: r["Z"] >= 55),
            "S3_test(N%6==0)": sub(lambda r: r["N"] % 6 == 0),
        }
        summary[name] = s
        print(name, json.dumps(s, indent=None, default=float)[:2000])
    with open(os.path.join(ROOT, "results", "bench_mendoza2011_summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1)


if __name__ == "__main__":
    main()
