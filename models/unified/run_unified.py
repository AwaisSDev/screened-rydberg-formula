"""LEGACY (stage 1 of the unified model). Builds/validates u29, u35, gshm_qed, pocket with the OLD S1/S2/S3 rule.
Since 2026-10-05 the FINAL model is pa_hier_rel (Push A), validated by validate_blind.py; this script no
longer writes results/uni_predictions.csv or uni_validation.json (it writes uni_legacy_validation.json) and
keeps "model" in uni_params.json pointing at the final model. Run: py -3.13 models/unified/run_unified.py (~1 min).

Writes results/uni_* (params, predictions, metrics, validation), results/model_comparison.{csv,md}.
Figures: models/unified/make_figures.py.
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
import umodel as U  # noqa: E402
import pocket  # noqa: E402
import fitlib  # noqa: E402  (Track B)
import gshm  # noqa: E402
sys.path.insert(0, _ROOT)
import evaluate  # noqa: E402

RES = os.path.join(_ROOT, "results")
a = U.dataset()
y = a["y"]
Z, N = fitlib.ZARR, fitlib.NARR
SPL = fitlib.SPLITS
t0 = time.time()


def met(p, mask):
    ape = np.abs(p - y) / y * 100
    neu = mask & (Z == N)
    return {"n": int(mask.sum()), "MAPE": float(ape[mask].mean()), "median": float(np.median(ape[mask])),
            "neutral_MAPE": float(ape[neu].mean()) if neu.any() else None, "n_neutral": int(neu.sum())}


def write_pred(name, p):
    path = os.path.join(RES, f"{name}_predictions.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Z", "N", "IE_pred_eV"])
        for zz, nn, v in zip(Z, N, p):
            w.writerow([int(zz), int(nn), repr(float(v))])
    m = evaluate.evaluate(path)
    with open(path.replace("_predictions.csv", "_metrics.json"), "w", encoding="utf-8") as f:
        json.dump(m, f, indent=1)
    return m


# ------------------------------------------------------------------ unified candidates
CANDS = {"uni_u29": {}, "uni_u35": dict(kappa_l=True, S2=True)}
RIDGES = [0.0, 0.02, 0.05, 0.1, 0.2]
inner_tr, inner_va = Z <= 44, (Z >= 45) & (Z <= 54)
out = {"ridge_selection": {}, "splits": {}, "all": {}, "params": {}}
fits = {}
for name, kw in CANDS.items():
    # ridge chosen WITHOUT looking at any test set: train Z<=44, validate on 45<=Z<=54
    sel = {}
    for lam in RIDGES:
        s = U.auto_tie(U.default_spec("class", **kw), inner_tr)
        th, _ = U.fit(s, inner_tr, ridge=lam)
        sel[lam] = met(U.evaluate(th, s, a), inner_va)["MAPE"]
    lam = min(sel, key=sel.get)
    out["ridge_selection"][name] = {"val_MAPE": {str(k): v for k, v in sel.items()}, "chosen": lam}
    for sp, test in SPL.items():
        s = U.auto_tie(U.default_spec("class", **kw), ~test)
        th, _ = U.fit(s, ~test, ridge=lam)
        out["splits"][f"{name}_{sp}"] = met(U.evaluate(th, s, a), test)
    s = U.default_spec("class", **kw)
    th, _ = U.fit(s, ridge=lam)
    s["ridge"] = lam
    fits[name] = (s, th)
    out["all"][name] = met(U.evaluate(th, s, a), np.ones(len(y), bool))
    out["params"][name] = len(th)
    print(name, "ridge", lam, out["all"][name], {k: round(v["MAPE"], 3) for k, v in out["splits"].items() if k.startswith(name)}, f"{time.time()-t0:.0f}s")

# ------------------------------------------------------------------ Track B (+ Track A QED layer)
with open(os.path.join(RES, "se_params.json"), encoding="utf-8") as f:
    sp_b = json.load(f)["final"]
specB, thB = sp_b["spec"], np.array(sp_b["values"])
aB = gshm.dataset("1e")


def qed_layer(ie, Ze):
    """Track A recoil + QED + FNS layer (parameter-free) added on top of any screened model."""
    return ie * a["mu"] - a["qedfns"] * (np.minimum(Ze, a["Z"]) / a["Z"]) ** 2


ieB, partsB = gshm.evaluate_model(thB, specB, aB, return_parts=True)
pB_qed = qed_layer(ieB, partsB["Ze"])
out["all"]["uni_gshm_qed"] = met(pB_qed, np.ones(len(y), bool))
# Referee fix (major issue 2): Track B's 'final' spec carries sigma_core in [0.8, 1.0], a bound
# added AFTER Track B saw its S2 result. For the held-out numbers used in SELECTION we therefore
# refit with that bound removed ("blind", pre-registered); the bounded (post-hoc) numbers are
# kept under *_posthoc_* for reference only.
specB_blind = fitlib.with_spec(specB, core_bounds=None)
for sp, test in SPL.items():
    for tag, spc in (("", specB_blind), ("_posthoc", specB)):
        s2 = fitlib.auto_tie(spc, ~test)
        th, _ = fitlib.fit(s2, ~test)
        ie2, parts2 = gshm.evaluate_model(th, s2, aB, return_parts=True)
        out["splits"][f"se_gshm{tag}_{sp}"] = met(ie2, test)
        out["splits"][f"uni_gshm_qed{tag}_{sp}"] = met(qed_layer(ie2, parts2["Ze"]), test)
print("B+QED", out["all"]["uni_gshm_qed"], {k: round(v["MAPE"], 3) for k, v in out["splits"].items() if "gshm" in k}, f"{time.time()-t0:.0f}s")

# ------------------------------------------------------------------ pocket formula
thP = pocket.fit()
pP = pocket.evaluate(thP, a)
out["all"]["uni_pocket"] = met(pP, np.ones(len(y), bool))
for sp, test in SPL.items():
    out["splits"][f"uni_pocket_{sp}"] = met(pocket.evaluate(pocket.fit(~test), a), test)
out["pocket_params"] = dict(zip(pocket.NAMES, map(float, thP)))
print("pocket", out["all"]["uni_pocket"], dict(zip(pocket.NAMES, np.round(thP, 4))))

# ------------------------------------------------------------------ choose final model
# pre-registered rule: lowest mean held-out MAPE over S1/S2/S3, using BLIND numbers only
score = {}
for nm in ["uni_u29", "uni_u35", "uni_gshm_qed"]:
    score[nm] = float(np.mean([out["splits"][f"{nm}_{sp}"]["MAPE"] for sp in SPL]))
out["selection_score_mean_heldout_MAPE"] = score
best = min(score, key=score.get)
out["chosen"] = best
print("held-out mean MAPE", score, "->", best)

# predictions
preds = {"uni_u29": U.evaluate(fits["uni_u29"][1], fits["uni_u29"][0], a),
         "uni_u35": U.evaluate(fits["uni_u35"][1], fits["uni_u35"][0], a),
         "uni_gshm_qed": pB_qed, "uni_pocket": pP}
for nm, p in preds.items():
    write_pred(nm, p)
mfinal = write_pred(best, preds[best])  # legacy: no longer overwrites uni_predictions.csv

# parameters
if best.startswith("uni_u"):
    s, th = fits[best]
    names = U.param_names(s)
    params = {"model": "pa_hier_rel", "legacy_stage1_choice": best, "spec": s, "names": names, "values": list(map(float, th)),
              "n_params": len(th)}
else:
    params = {"model": "pa_hier_rel", "legacy_stage1_choice": best, "spec": specB, "values": list(map(float, thB)), "n_params": len(thB)}
params["pocket"] = {"names": pocket.NAMES, "values": list(map(float, thP))}
params["all_candidates"] = {nm: {"spec": fits[nm][0], "names": U.param_names(fits[nm][0]),
                                 "values": list(map(float, fits[nm][1]))} for nm in fits}
with open(os.path.join(RES, "uni_params.json"), "w", encoding="utf-8") as f:
    json.dump(params, f, indent=1)
with open(os.path.join(RES, "uni_legacy_validation.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, indent=1)
print("final ALL", mfinal["ALL"]["MAPE_%"], "neutral", mfinal["first_IE_neutral_atoms"]["MAPE_%"],
      f"{time.time()-t0:.0f}s")
