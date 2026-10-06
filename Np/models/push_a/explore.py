"""Exploration on the SELECTION splits only (V1, V2, S1, S3). Never computes anything on Z >= 55
as a test set (Z >= 55 rows are only in the S1/S3 training sets, as the protocol prescribes).
Usage: py -3.13 models/push_a/explore.py <tag> [json spec overrides ...]"""
import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import model as M  # noqa: E402
import umodel as U  # noqa: E402

a = U.dataset()
y = a["y"]
Z = a["Z"].astype(int)
N = (a["Z"] - a["Zq"] + 1).astype(int)
SEL = {"V1": (Z <= 36, (Z >= 37) & (Z <= 54)),
       "V2": (Z <= 44, (Z >= 45) & (Z <= 54)),
       "S1": (Z % 5 != 0, Z % 5 == 0),
       "S3": (N % 6 != 0, N % 6 == 0)}


def met(p, m):
    ape = np.abs(p - y) / y * 100
    neu = m & (Z == N)
    return {"MAPE": float(ape[m].mean()), "median": float(np.median(ape[m])),
            "neutral_MAPE": float(ape[neu].mean()) if neu.any() else None}


def run(spec):
    out = {}
    for k, (tr, te) in SEL.items():
        s, th = M.fit(tr, dict(spec))
        out[k] = met(M.evaluate(th, s, a), te)
    out["score"] = float(np.mean([out[k]["MAPE"] for k in SEL]))
    return out


if __name__ == "__main__":
    t = time.time()
    tag = sys.argv[1]
    over = json.loads(sys.argv[2]) if len(sys.argv) > 2 else {}
    spec = M.default_spec(**over)
    r = run(spec)
    line = f"{tag:28s} score {r['score']:7.3f} | " + " ".join(
        f"{k} {r[k]['MAPE']:6.2f}/{r[k]['neutral_MAPE'] if r[k]['neutral_MAPE'] is None else round(r[k]['neutral_MAPE'],1)}" for k in SEL) + f"  np={len(M.param_names(M.auto_tie(spec, np.ones(len(y),bool))))} {time.time()-t:.0f}s"
    print(line)
    with open(os.path.join(_HERE, "explore_log.txt"), "a", encoding="utf-8") as f:
        f.write(line + "  " + json.dumps(over) + "\n")
