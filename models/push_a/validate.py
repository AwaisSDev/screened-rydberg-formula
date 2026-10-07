"""Pre-registered validation protocol for Push A (run: py -3.13 models/push_a/validate.py, ~1 min).

For every candidate (Push A models + the u35 and pocket references, refitted per split with their
own fit functions):
  selection splits  V1 (fit Z<=36 / val 37..54), V2 (fit Z<=44 / val 45..54),
                    S1 (fit Z%5!=0 / test Z%5==0), S3 (fit N%6!=0 / test N%6==0)
  selection score = mean MAPE over V1, V2, S1, S3
  S2 (fit Z<=54 / test Z>=55): computed ONCE, after the candidate list and all hyper-parameters
                    were frozen, for reporting only
  all-data fit      scored with evaluate.py (results/pa_<name>_predictions.csv)
Chosen model = lowest selection score among the Push A candidates.
Writes results/pa_validation.json, results/pa_params.json, results/pa_*_predictions.csv.
"""
import csv
import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)
import model as M  # noqa: E402
import umodel as U  # noqa: E402
import pocket  # noqa: E402
import evaluate as EV  # noqa: E402

RES = os.path.join(_ROOT, "results")
a = U.dataset()
y = a["y"]
Z = a["Z"].astype(int)
N = (a["Z"] - a["Zq"] + 1).astype(int)
SPLITS = {"V1": (Z <= 36, (Z >= 37) & (Z <= 54)),
          "V2": (Z <= 44, (Z >= 45) & (Z <= 54)),
          "S1": (Z % 5 != 0, Z % 5 == 0),
          "S3": (N % 6 != 0, N % 6 == 0),
          "S2": (Z <= 54, Z >= 55)}
SEL = ["V1", "V2", "S1", "S3"]

# ---- frozen candidate list (all choices made on the selection splits only; see NOTES.md)
PA = {
    "pa_hier_rel": M.default_spec(grouping="pocket", rel="fs", hier=True, ridge=1e-4),
    "pa_hier": M.default_spec(grouping="pocket", rel=None, hier=True, ridge=1e-4),
    "pa_bound9": M.default_spec(grouping="pocket", rel=None),
    "pa_bound9_fs0": M.default_spec(grouping="pocket", rel="fs0"),
    "pa_bound14_relfit": M.default_spec(grouping="pocket", rel="fs"),
    "pa_nobound14_relfit": M.default_spec(grouping="pocket", rel="fs", bound=False),
}
U35 = U.default_spec("class", kappa_l=True, S2=True)


def fit_pred(name, train):
    if name in PA:
        s, th = M.fit(train, dict(PA[name]))
        return M.evaluate(th, s, a), (s, th), len(th)
    if name == "ref_u35":
        s = U.auto_tie(U35, train)
        th, _ = U.fit(s, train, ridge=0.02)
        return U.evaluate(th, s, a), (s, th), len(th)
    if name == "ref_pocket":
        th = pocket.fit(train)
        return pocket.evaluate(th, a), (None, th), len(th)
    raise KeyError(name)


def met(p, m):
    ape = np.abs(p - y) / y * 100
    neu = m & (Z == N)
    return {"n": int(m.sum()), "MAPE": float(ape[m].mean()), "median": float(np.median(ape[m])),
            "neutral_MAPE": float(ape[neu].mean()) if neu.any() else None,
            "max": float(ape[m].max())}


def write_pred(name, p):
    path = os.path.join(RES, f"{name}_predictions.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Z", "N", "IE_pred_eV"])
        for zz, nn, v in zip(Z, N, p):
            w.writerow([int(zz), int(nn), repr(float(v))])
    m = EV.evaluate(path)
    with open(path.replace("_predictions.csv", "_metrics.json"), "w", encoding="utf-8") as f:
        json.dump(m, f, indent=1)
    return path, m


def main():
    t0 = time.time()
    out = {"protocol": __doc__, "candidates": {}}
    params = {"candidates": {}}
    names = list(PA) + ["ref_u35", "ref_pocket"]
    # 1) selection splits
    for nm in names:
        r = {}
        for sp in SEL:
            tr, te = SPLITS[sp]
            p, _, _ = fit_pred(nm, tr)
            r[sp] = met(p, te)
        r["selection_score"] = float(np.mean([r[sp]["MAPE"] for sp in SEL]))
        out["candidates"][nm] = r
        print(nm, "score", round(r["selection_score"], 3), f"{time.time()-t0:.0f}s", flush=True)
    chosen = min(PA, key=lambda k: out["candidates"][k]["selection_score"])
    out["chosen"] = chosen
    print("chosen (selection score only):", chosen, flush=True)
    # 2) S2, once, for reporting
    for nm in names:
        tr, te = SPLITS["S2"]
        p, _, _ = fit_pred(nm, tr)
        out["candidates"][nm]["S2"] = met(p, te)
    # 3) all-data fit
    for nm in names:
        p, (s, th), npar = fit_pred(nm, np.ones(len(y), bool))
        path, m = write_pred(nm.replace("ref_", "pa_ref_") if nm.startswith("ref_") else nm, p)
        out["candidates"][nm]["all"] = {"n_params": npar, "MAPE": m["ALL"]["MAPE_%"],
                                        "median": m["ALL"]["median_APE_%"],
                                        "neutral_MAPE": m["first_IE_neutral_atoms"]["MAPE_%"],
                                        "Hlike_MAPE": m["hydrogen_like"]["MAPE_%"], "csv": path}
        if nm in PA:
            params["candidates"][nm] = {"spec": s, "names": M.param_names(s),
                                        "values": list(map(float, th)), "n_params": npar}
        print(nm, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out["candidates"][nm]["all"].items()},
              "S2", round(out["candidates"][nm]["S2"]["MAPE"], 3), flush=True)
    params["chosen"] = chosen
    with open(os.path.join(RES, "pa_params.json"), "w", encoding="utf-8") as f:
        json.dump(params, f, indent=1)
    with open(os.path.join(RES, "pa_validation.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
