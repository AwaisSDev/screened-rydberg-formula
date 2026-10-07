"""Pre-registered validation protocol for Push-B candidates and the u35 / pocket references.

    py -3.13 models/push_b/validate.py        (~2 min)

Order (fixed): V1 (fit Z<=36, validate 37..54), V2 (fit Z<=44, validate 45..54), S1 (Z%5), S3 (N%6)
-> selection score = mean MAPE(V1,V2,S1,S3); THEN S2 (fit Z<=54, test Z>=55) once per candidate,
reporting only; THEN all-data fit -> results/pb_<name>_predictions.csv scored with evaluate.py.
The chosen Push-B model (lowest selection score among pb_* candidates) is written to results/pb_params.json.
Writes results/pb_validation.json.
"""
import csv
import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
for _p in (_ROOT, _HERE, os.path.join(_ROOT, "models", "unified")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import model as M  # noqa: E402
import umodel as U  # noqa: E402
import pocket  # noqa: E402
import evaluate  # noqa: E402

RES = os.path.join(_ROOT, "results")
a = U.dataset()
y = a["y"]
Z, N = a["rec_Z"], a["rec_N"]
SEL = {"V1": (Z <= 36, (Z >= 37) & (Z <= 54)), "V2": (Z <= 44, (Z >= 45) & (Z <= 54)),
       "S1": (Z % 5 != 0, Z % 5 == 0), "S3": (N % 6 != 0, N % 6 == 0)}
S2 = (Z <= 54, Z >= 55)
ALL = np.ones(len(y), bool)


def met(p, m):
    ape = np.abs(p - y) / y * 100
    nm = m & (Z == N)
    return {"n": int(m.sum()), "MAPE": float(ape[m].mean()), "median": float(np.median(ape[m])),
            "neutral_MAPE": float(ape[nm].mean()) if nm.any() else None, "n_neutral": int(nm.sum())}


# ---------------------------------------------------------------- candidate definitions
def pb(spec_kw, ridge=0.0):
    def fitpred(mask):
        s, th = M.fit(M.default_spec(**spec_kw), mask, ridge=ridge)
        return M.evaluate(th, s, a), (s, th)
    return fitpred


def u35(mask):
    s = U.auto_tie(U.default_spec("class", kappa_l=True, S2=True), mask)
    th, _ = U.fit(s, mask, ridge=0.02)      # ridge 0.02 = u35's own (inner-split-selected) value
    return U.evaluate(th, s, a), (s, th)


def pocket_ref(mask):
    th = pocket.fit(mask)
    return pocket.evaluate(th, a), ({"model": "pocket"}, th)


# Order and content frozen before any S2 evaluation. Variant history (selection scores only) in NOTES.md.
CANDS = {
    "ref_u35": (u35, 35),
    "ref_pocket": (pocket_ref, 8),
    "pb_clip_pos": (pb(dict(form="clip", rel_form="u35", tpos=True)), None),
    "pb_clip": (pb(dict(form="clip", rel_form="u35", tpos=False)), None),
    "pb_exp_pos": (pb(dict(form="exp", rel_form="sat", tpos=True)), None),
    "pb_logistic_pos": (pb(dict(form="logistic", rel_form="u35", tpos=True)), None),
}


def n_free(spec, th):
    return int(len(th))


def write_pred(name, p):
    path = os.path.join(RES, f"{name}_predictions.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Z", "N", "IE_pred_eV"])
        for zz, nn, v in zip(Z, N, p):
            w.writerow([int(zz), int(nn), repr(float(v))])
    return path, evaluate.evaluate(path)


def main():
    t0 = time.time()
    out = {"protocol": "V1,V2,S1,S3 -> selection score (mean MAPE); S2 computed once after freezing; "
                       "all-data fit scored with evaluate.py", "candidates": {}}
    # stage 1: selection splits for all candidates
    for nm, (fp, _) in CANDS.items():
        r = {}
        for sp, (tr, te) in SEL.items():
            p, _ = fp(tr)
            r[sp] = met(p, te)
        r["selection_score"] = float(np.mean([r[s]["MAPE"] for s in SEL]))
        out["candidates"][nm] = r
        print(nm, {s: round(r[s]["MAPE"], 3) for s in SEL}, "score", round(r["selection_score"], 3),
              f"{time.time()-t0:.0f}s", flush=True)
    pbs = {k: v["selection_score"] for k, v in out["candidates"].items() if k.startswith("pb_")}
    chosen = min(pbs, key=pbs.get)
    out["chosen"] = chosen
    out["chosen_rule"] = "lowest selection score among pb_* candidates (frozen before S2)"
    print("CHOSEN (frozen):", chosen, flush=True)
    # stage 2: S2 once per candidate (reporting only)
    for nm, (fp, _) in CANDS.items():
        p, _ = fp(S2[0])
        out["candidates"][nm]["S2"] = met(p, S2[1])
        print(nm, "S2", {k: round(v, 3) if isinstance(v, float) else v for k, v in out["candidates"][nm]["S2"].items()},
              flush=True)
    # stage 3: all-data fit
    for nm, (fp, npar) in CANDS.items():
        p, (s, th) = fp(ALL)
        path, m = write_pred("pb_" + nm if not nm.startswith("pb_") else nm, p)
        r = out["candidates"][nm]
        r["n_params"] = npar if npar is not None else n_free(s, th)
        r["all"] = {"MAPE": m["ALL"]["MAPE_%"], "median": m["ALL"]["median_APE_%"],
                    "neutral_MAPE": m["first_IE_neutral_atoms"]["MAPE_%"],
                    "Hlike_MAPE": m["hydrogen_like"]["MAPE_%"], "max_APE": m["ALL"]["max_APE_%"],
                    "predictions_csv": os.path.relpath(path, _ROOT).replace("\\", "/")}
        if nm.startswith("pb_"):
            r["spec"] = s
            r["names"] = M.param_names(s)
            r["values"] = list(map(float, th))
            if nm == chosen:
                with open(os.path.join(RES, "pb_params.json"), "w", encoding="utf-8") as f:
                    json.dump({"model": nm, "spec": s, "names": M.param_names(s), "values": list(map(float, th)),
                               "n_params": len(th)}, f, indent=1)
        print(nm, "ALL", round(r["all"]["MAPE"], 3), "neutral", round(r["all"]["neutral_MAPE"], 2),
              "H", r["all"]["Hlike_MAPE"], "npar", r["n_params"], flush=True)
    with open(os.path.join(RES, "pb_validation.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
