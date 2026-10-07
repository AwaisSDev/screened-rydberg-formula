"""Write predictions for the zero-parameter analytic 1/Z-expansion variants and score them."""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zexp  # noqa: E402
from zexp import ROOT  # noqa: E402
sys.path.insert(0, ROOT)
import evaluate  # noqa: E402
from common.atomdata import load_records  # noqa: E402

VARIANTS = ["bohr", "zexp2", "zexp2_rel", "zexp_sq", "zexp_sq_rel", "zexp_sq_rel_p4", "zexp"]


def write_preds(name, fn):
    path = os.path.join(ROOT, "results", f"fp_{name}_predictions.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Z", "N", "IE_pred_eV"])
        for r in load_records():
            try:
                v = fn(r["Z"], r["N"])
            except Exception as e:  # noqa: BLE001
                print("fail", r["Z"], r["N"], e)
                v = float("nan")
            w.writerow([r["Z"], r["N"], repr(v)])
    m = evaluate.evaluate(path)
    json.dump(m, open(os.path.join(ROOT, "results", f"fp_{name}_metrics.json"), "w"), indent=1)
    return m


def headline(m):
    keys = ["ALL", "first_IE_neutral_atoms", "hydrogen_like", "N<=10", "11<=N<=36", "N>=37"]
    return {k: round(m[k]["MAPE_%"], 3) for k in keys if k in m}


if __name__ == "__main__":
    out = {}
    for v in (sys.argv[1:] or VARIANTS):
        if v == "zexp":            # headline analytic model = zexp_sq_rel with Lande p = 2
            fn = lambda Z, N: zexp.predict(Z, N, variant="zexp_sq_rel", p=2.0)  # noqa: E731
        elif v == "zexp_sq_rel_p4":
            fn = lambda Z, N: zexp.predict(Z, N, variant="zexp_sq_rel", p=4.0)  # noqa: E731
        else:
            fn = lambda Z, N, v=v: zexp.predict(Z, N, variant=v)  # noqa: E731
        m = write_preds(v, fn)
        out[v] = headline(m)
        print(v, out[v], flush=True)
    zexp.save_cfg_cache()
