"""Sensitivity check requested by the final audit (reporting only, no choice is made here).

The pre-registered selection splits S1 and S3 contain Z >= 55 rows in both train and test
(~76 % of their test rows).  This script recomputes the selection score with every split
restricted to Z <= 54 (train AND test) for the frozen leading candidates, to check whether
the winner depends on heavy-atom interpolation accuracy.  It also counts the Z >= 55 share
of the S1/S3 rows.  Output: results/uni_sensitivity_z54.json.

Run: py -3.13 models/unified/sensitivity_z54.py   (~2-4 min)
"""
import sys, os, json, importlib.util
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "models", "unified"))
import umodel as U
from push_a_import import PA
spec = importlib.util.spec_from_file_location("pbm", os.path.join(ROOT, "models", "push_b", "model.py"))
PB = importlib.util.module_from_spec(spec); spec.loader.exec_module(PB)

a = U.dataset(); y = a["y"]; Z = a["Z"].astype(int); N = (a["Z"] - a["Zq"] + 1).astype(int)
L = Z <= 54
counts = {
    "S1_test_rows": int((Z % 5 == 0).sum()), "S1_test_rows_Zge55": int(((Z % 5 == 0) & ~L).sum()),
    "S1_train_rows_Zge55": int(((Z % 5 != 0) & ~L).sum()),
    "S3_test_rows": int((N % 6 == 0).sum()), "S3_test_rows_Zge55": int(((N % 6 == 0) & ~L).sum()),
    "S3_train_rows_Zge55": int(((N % 6 != 0) & ~L).sum()),
}
print(counts, flush=True)
SPL = {"V1": (Z <= 36, (Z >= 37) & L), "V2": (Z <= 44, (Z >= 45) & L),
       "S1_Zle54": (L & (Z % 5 != 0), L & (Z % 5 == 0)), "S3_Zle54": (L & (N % 6 != 0), L & (N % 6 == 0))}
C = {"pa_hier_rel": ("A", dict(grouping="pocket", rel="fs", hier=True, ridge=1e-4)),
     "pa_hier": ("A", dict(grouping="pocket", rel=None, hier=True, ridge=1e-4)),
     "pa_bound9": ("A", dict(grouping="pocket", rel=None)),
     "pa_bound14_relfit": ("A", dict(grouping="pocket", rel="fs")),
     "pb_clip_pos": ("B", dict(form="clip", rel_form="u35", tpos=True)),
     "pb_exp_pos": ("B", dict(form="exp", rel_form="sat", tpos=True))}
out = {"note": "Selection score recomputed with S1/S3 restricted to Z<=54 (train and test). Reporting only; "
               "the pre-registered choice is unchanged.", "counts": counts, "candidates": {}}
for nm, (w, kw) in C.items():
    r = {}
    for k, (tr, te) in SPL.items():
        if w == "A":
            s, th = PA.fit(tr, PA.default_spec(**kw)); p = PA.evaluate(th, s, a)
        else:
            s, th = PB.fit(PB.default_spec(**kw), tr, ridge=0.0); p = PB.evaluate(th, s, a)
        r[k] = round(float((np.abs(p - y) / y * 100)[te].mean()), 3)
    r["restricted_score"] = round(float(np.mean([r[k] for k in SPL])), 3)
    out["candidates"][nm] = r
    print(nm, r, flush=True)
with open(os.path.join(ROOT, "results", "uni_sensitivity_z54.json"), "w") as f:
    json.dump(out, f, indent=1)
