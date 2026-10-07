"""Pre-registered blind validation of the FINAL model and the references (one re-runnable script).

    py -3.13 models/unified/validate_blind.py            (all models, about 3-5 min)
    py -3.13 models/unified/validate_blind.py --only final ref_u35   (subset; merged into the JSON)

Protocol (identical for every model; refit on each training mask with the model's own fit function):
  V1  fit Z <= 36            validate 37 <= Z <= 54
  V2  fit Z <= 44            validate 45 <= Z <= 54
  S1  fit Z % 5 != 0         test Z % 5 == 0
  S3  fit N % 6 != 0         test N % 6 == 0
  selection score = mean MAPE over V1, V2, S1, S3   (the ONLY quantity used to choose a model)
  S2  fit Z <= 54            test Z >= 55            (BLIND: reported, never used for any choice)
  all-data fit on all 5847 rows, scored with evaluate.py.
For every split: MAPE, median APE and neutral-first-IE MAPE of the test/validation rows.

Models
  final       = pa_hier_rel (Push A, 33 params): exact sigma1 + bounded hierarchical screening remainder
                (models/push_a/model.py, imported read-only; spec frozen by Push A before its S2 run)
  ref_u35     = previous unified final (35 params, ridge 0.02)       models/unified/umodel.py
  ref_u29     = previous unified variant (29 params, ridge 0.0)      models/unified/umodel.py
  ref_pocket  = 8-parameter pocket formula                            models/unified/pocket.py
  ref_gshm    = Track B GSHM final (32 params) WITHOUT the post-hoc sigma_core bound (blind), and
  ref_gshm_core = Track B GSHM core (30 params), same treatment       models/semi_empirical
  ref_slater_total = Slater's rules, total-energy difference (0 params; test MAPE = MAPE on test rows)

Writes results/uni_validation.json, results/uni_validation.md, results/uni_final_params.json,
results/uni_predictions.csv + uni_metrics.json (final, all-data fit), results/uni_ref_*_predictions.csv.
predict() of every model reads only configurations; NIST IEs enter only as fit targets and for scoring.
"""
import argparse
import csv
import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
for _p in (_ROOT, _HERE, os.path.join(_ROOT, "models", "semi_empirical")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import umodel as U  # noqa: E402
import pocket  # noqa: E402
import gshm  # noqa: E402
import fitlib  # noqa: E402
import evaluate as EV  # noqa: E402
from push_a_import import PA  # noqa: E402

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
FINAL = "final"
FINAL_NAME = "pa_hier_rel"
FINAL_SPEC = PA.default_spec(grouping="pocket", rel="fs", hier=True, ridge=1e-4)   # = push_a/validate.py PA['pa_hier_rel']
U35 = U.default_spec("class", kappa_l=True, S2=True)
U29 = U.default_spec("class")
with open(os.path.join(RES, "se_params.json"), encoding="utf-8") as _f:
    _SE = json.load(_f)
GSHM_BLIND = fitlib.with_spec(_SE["final"]["spec"], core_bounds=None)     # referee major issue 2: no post-hoc bound
GSHM_CORE_BLIND = fitlib.with_spec(_SE["core"]["spec"], core_bounds=None)
_aB = None

INFO = {
    FINAL: {"label": "FINAL pa_hier_rel (exact sigma1 + bounded hierarchical remainder)", "n_params": 33, "type": "hybrid"},
    "ref_u35": {"label": "previous unified final u35 (exact sigma1 + linear remainder)", "n_params": 35, "type": "hybrid"},
    "ref_u29": {"label": "previous unified variant u29", "n_params": 29, "type": "hybrid"},
    "ref_pocket": {"label": "pocket formula", "n_params": 8, "type": "fitted"},
    "ref_gshm": {"label": "Track B GSHM final (blind, no post-hoc bound)", "n_params": 32, "type": "fitted"},
    "ref_gshm_core": {"label": "Track B GSHM core (blind, no post-hoc bound)", "n_params": 30, "type": "fitted"},
    "ref_slater_total": {"label": "Slater's rules, total-energy difference", "n_params": 0, "type": "empirical rules (0 params)"},
}
ORDER = list(INFO)


def _slater():
    with open(os.path.join(RES, "se_slater_total_predictions.csv"), encoding="utf-8") as f:
        d = {(int(r["Z"]), int(r["N"])): float(r["IE_pred_eV"]) for r in csv.DictReader(f)}
    return np.array([d[(int(z), int(n))] for z, n in zip(Z, N)])


def fit_pred(name, train):
    """Fit `name` on the rows in `train`; return (predictions for all 5847 rows, fitted params, n_params)."""
    global _aB
    if name == FINAL:
        s, th = PA.fit(train, dict(FINAL_SPEC))
        return PA.evaluate(th, s, a), (s, th), len(th)
    if name in ("ref_u35", "ref_u29"):
        base, lam = (U35, 0.02) if name == "ref_u35" else (U29, 0.0)
        s = U.auto_tie(base, train)
        th, _ = U.fit(s, train, ridge=lam)
        return U.evaluate(th, s, a), (s, th), len(th)
    if name == "ref_pocket":
        th = pocket.fit(train)
        return pocket.evaluate(th, a), (None, th), len(th)
    if name in ("ref_gshm", "ref_gshm_core"):
        if _aB is None:
            _aB = gshm.dataset("1e")
        s = fitlib.auto_tie(GSHM_BLIND if name == "ref_gshm" else GSHM_CORE_BLIND, train)
        th, _ = fitlib.fit(s, train)
        return gshm.evaluate_model(th, s, _aB), (s, th), len(th)
    if name == "ref_slater_total":
        return _slater(), (None, np.zeros(0)), 0
    raise KeyError(name)


def met(p, m):
    ape = np.abs(p - y) / y * 100
    neu = m & (Z == N)
    return {"n": int(m.sum()), "MAPE": float(ape[m].mean()), "median": float(np.median(ape[m])),
            "neutral_MAPE": float(ape[neu].mean()) if neu.any() else None, "n_neutral": int(neu.sum()),
            "max": float(ape[m].max())}


def write_pred(fname, p):
    path = os.path.join(RES, f"{fname}_predictions.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Z", "N", "IE_pred_eV"])
        for zz, nn, v in zip(Z, N, p):
            w.writerow([int(zz), int(nn), repr(float(v))])
    m = EV.evaluate(path)
    with open(path.replace("_predictions.csv", "_metrics.json"), "w", encoding="utf-8") as f:
        json.dump(m, f, indent=1)
    return path, m


def run_model(nm, t0):
    r = {"label": INFO[nm]["label"], "type": INFO[nm]["type"]}
    for sp in SEL:                                   # 1) selection splits
        tr, te = SPLITS[sp]
        p, _, _ = fit_pred(nm, tr)
        r[sp] = met(p, te)
    r["selection_score"] = float(np.mean([r[sp]["MAPE"] for sp in SEL]))
    tr, te = SPLITS["S2"]                            # 2) blind S2, reporting only
    p, _, _ = fit_pred(nm, tr)
    r["S2"] = met(p, te)
    p, (s, th), npar = fit_pred(nm, np.ones(len(y), bool))   # 3) all-data fit
    fname = "uni" if nm == FINAL else "uni_" + nm
    path, m = write_pred(fname, p)
    r["n_params"] = npar
    r["all"] = {"MAPE": m["ALL"]["MAPE_%"], "median": m["ALL"]["median_APE_%"],
                "neutral_MAPE": m["first_IE_neutral_atoms"]["MAPE_%"],
                "Hlike_MAPE": m["hydrogen_like"]["MAPE_%"], "csv": os.path.relpath(path, _ROOT).replace("\\", "/")}
    if nm == FINAL:
        with open(os.path.join(RES, "uni_final_params.json"), "w", encoding="utf-8") as f:
            json.dump({"model": FINAL_NAME, "module": "models/push_a/model.py", "spec": s,
                       "names": PA.param_names(s), "values": list(map(float, th)), "n_params": npar,
                       "note": "all-data fit (5847 NIST rows) written by models/unified/validate_blind.py"},
                      f, indent=1)
    print(f"{nm:18s} score {r['selection_score']:8.3f}  V1 {r['V1']['MAPE']:.3f} V2 {r['V2']['MAPE']:.3f} "
          f"S1 {r['S1']['MAPE']:.3f} S3 {r['S3']['MAPE']:.3f} | S2 {r['S2']['MAPE']:.3f} | all {r['all']['MAPE']:.3f} "
          f"({npar} p)  {time.time()-t0:.0f}s", flush=True)
    return r


def fmt(v, d=2):
    if v is None:
        return "n/a"
    if abs(v) >= 1e4:
        return f"{v:.2e}"
    if abs(v) < 0.1:
        return f"{v:.4f}"
    return f"{v:.{d}f}"


def write_md(out):
    L = ["# Blind validation of the final model (pre-registered protocol)", "",
         "Generated by `py -3.13 models/unified/validate_blind.py`. Every model is refit on each training mask "
         "with its own fit function. Selection score = mean MAPE of V1, V2, S1, S3; S2 (fit Z<=54, test Z>=55) is "
         "BLIND, computed after all choices were frozen, never used for selection. All values in %. "
         "Caveat: S1/S3, as pre-registered, contain Z>=55 rows in train and test (about 76% of their test rows), so heavy-atom "
         "interpolation accuracy did inform the selection score; with S1/S3 restricted to Z<=54 the winner is unchanged "
         "(models/unified/sensitivity_z54.py, results/uni_sensitivity_z54.json).", "",
         "| model | params | type | selection score | V1 | V2 | S1 | S3 | **blind S2** | S2 median | S2 neutral | all-data | all median | all neutral | H-like |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for nm in ORDER:
        if nm not in out["models"]:
            continue
        r = out["models"][nm]
        L.append(f"| {r['label']} | {r['n_params']} | {r['type']} | {fmt(r['selection_score'])} | "
                 + " | ".join(fmt(r[s]["MAPE"]) for s in SEL)
                 + f" | **{fmt(r['S2']['MAPE'])}** | {fmt(r['S2']['median'])} | {fmt(r['S2']['neutral_MAPE'])} | "
                 f"{fmt(r['all']['MAPE'])} | {fmt(r['all']['median'])} | {fmt(r['all']['neutral_MAPE'])} | {fmt(r['all']['Hlike_MAPE'])} |")
    L += ["", "Per-split detail (MAPE / median / neutral-first-IE MAPE, n rows):", ""]
    for nm in ORDER:
        if nm not in out["models"]:
            continue
        r = out["models"][nm]
        L.append(f"- **{nm}**: " + "; ".join(
            f"{s} {fmt(r[s]['MAPE'])} / {fmt(r[s]['median'])} / {fmt(r[s]['neutral_MAPE'])} (n={r[s]['n']})"
            for s in SEL + ["S2"]))
    if out.get("consistency"):
        L += ["", "Consistency with the push agents' own validation runs (|difference| in MAPE points):", ""]
        for k, v in out["consistency"].items():
            L.append(f"- {k}: {v}")
    with open(os.path.join(RES, "uni_validation.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


def consistency(out):
    """Compare with results/pa_validation.json and pb_validation.json for the same models."""
    res = {}
    try:
        pa = json.load(open(os.path.join(RES, "pa_validation.json"), encoding="utf-8"))["candidates"]
        pairs = {FINAL: "pa_hier_rel", "ref_u35": "ref_u35", "ref_pocket": "ref_pocket"}
        for mine, theirs in pairs.items():
            if mine in out["models"] and theirs in pa:
                d = {s: abs(out["models"][mine][s]["MAPE"] - pa[theirs][s]["MAPE"]) for s in SEL + ["S2"]}
                d["all"] = abs(out["models"][mine]["all"]["MAPE"] - pa[theirs]["all"]["MAPE"])
                res[f"{mine} vs push A {theirs}"] = {k: float(f"{v:.2e}") for k, v in d.items()}
    except (OSError, KeyError) as e:
        res["push_a"] = f"not compared: {e}"
    try:
        pb = json.load(open(os.path.join(RES, "pb_validation.json"), encoding="utf-8"))
        pbc = pb.get("candidates", pb)
        for mine, theirs in {"ref_u35": "ref_u35", "ref_pocket": "ref_pocket"}.items():
            if mine in out["models"] and theirs in pbc:
                t = pbc[theirs]
                d = {}
                for s in SEL + ["S2"]:
                    v = t.get(s, {})
                    v = v.get("MAPE") if isinstance(v, dict) else t.get(f"{s}_MAPE")
                    if v is not None:
                        d[s] = abs(out["models"][mine][s]["MAPE"] - v)
                res[f"{mine} vs push B {theirs}"] = {k: float(f"{v:.2e}") for k, v in d.items()}
    except (OSError, KeyError) as e:
        res["push_b"] = f"not compared: {e}"
    # Track B blind S2 from the earlier unified run (run_unified.py, same blind spec)
    if "ref_gshm" in out["models"]:
        res["ref_gshm S2 vs earlier se_gshm_S2 (170.27)"] = float(f"{abs(out['models']['ref_gshm']['S2']['MAPE'] - 170.27):.2e}")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()
    t0 = time.time()
    path = os.path.join(RES, "uni_validation.json")
    out = {}
    if args.only and os.path.exists(path):
        out = json.load(open(path, encoding="utf-8"))
    if "models" not in out:
        out = {"protocol": __doc__, "final_model": FINAL_NAME, "models": {}}
    for nm in (args.only or ORDER):
        out["models"][nm] = run_model(nm, t0)
    out["selection_ranking"] = sorted(((k, v["selection_score"]) for k, v in out["models"].items()),
                                      key=lambda kv: kv[1])
    out["consistency"] = consistency(out)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    write_md(out)
    print(json.dumps(out["consistency"], indent=1))
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
