"""Referee-response computations for the unified stage (run: py -3.13 models/unified/referee_fixes.py, ~1 min).

1. BLIND S2 for the GSHM-based candidate: Track B's 'final' spec carries a sigma_core bound
   [0.8, 1.0] that was added after Track B saw its own S2 result. Here the bound is removed
   (everything else, incl. the ridge prior 0.1 and the auto-tying of unsupported f parameters,
   unchanged) and S1/S2/S3 are refitted and scored, with and without Track A's QED layer.
2. The pre-registered selection rule (lowest mean held-out MAPE over S1/S2/S3) is re-applied
   with the blind S2 numbers.
3. Superheavy sanity check (Z > 110) against published relativistic coupled-cluster estimates.

Writes results/uni_referee.json.
"""
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _HERE)
import umodel as U  # noqa: E402
import fitlib  # noqa: E402
import gshm  # noqa: E402

RES = os.path.join(_ROOT, "results")
a = U.dataset()
y = a["y"]
Z, N = fitlib.ZARR, fitlib.NARR
SPL = fitlib.SPLITS


def met(p, mask):
    ape = np.abs(p - y) / y * 100
    neu = mask & (Z == N)
    return {"n": int(mask.sum()), "MAPE": float(ape[mask].mean()), "median": float(np.median(ape[mask])),
            "neutral_MAPE": float(ape[neu].mean()) if neu.any() else None}


def qed_layer(ie, Ze):
    return ie * a["mu"] - a["qedfns"] * (np.minimum(Ze, a["Z"]) / a["Z"]) ** 2


with open(os.path.join(RES, "se_params.json"), encoding="utf-8") as f:
    specB = json.load(f)["final"]["spec"]
aB = gshm.dataset("1e")
out = {"blind": {}}
for label, spec in (("bounded_posthoc", specB), ("blind_no_core_bound", fitlib.with_spec(specB, core_bounds=None))):
    for sp, test in SPL.items():
        s2 = fitlib.auto_tie(spec, ~test)
        th, _ = fitlib.fit(s2, ~test)
        ie, parts = gshm.evaluate_model(th, s2, aB, return_parts=True)
        out["blind"][f"se_gshm_{label}_{sp}"] = met(ie, test)
        out["blind"][f"uni_gshm_qed_{label}_{sp}"] = met(qed_layer(ie, parts["Ze"]), test)
        print(label, sp, {k: round(v, 3) for k, v in out["blind"][f"uni_gshm_qed_{label}_{sp}"].items() if v})

with open(os.path.join(RES, "uni_validation.json"), encoding="utf-8") as f:
    val = json.load(f)
sc = {}
for nm in ["uni_u29", "uni_u35"]:
    sc[nm] = float(np.mean([val["splits"][f"{nm}_{sp}"]["MAPE"] for sp in SPL]))
sc["uni_gshm_qed (blind S2)"] = float(np.mean([out["blind"][f"uni_gshm_qed_blind_no_core_bound_{sp}"]["MAPE"] for sp in SPL]))
sc["uni_gshm_qed (post-hoc bound)"] = float(np.mean([out["blind"][f"uni_gshm_qed_bounded_posthoc_{sp}"]["MAPE"] for sp in SPL]))
out["selection_score_mean_heldout_MAPE"] = sc
blind_only = {k: v for k, v in sc.items() if "post-hoc" not in k}
out["chosen_blind"] = min(blind_only, key=blind_only.get)
print("selection", sc, "->", out["chosen_blind"])

# superheavy sanity check (neutral first IEs; literature = relativistic coupled-cluster estimates, approx.)
LIT = {(118, 118): 8.9, (120, 120): 5.85, (114, 114): 8.5, (112, 112): 11.97, (119, 119): 4.8}
sys.path.insert(0, _ROOT)
import final  # noqa: E402
sh = {}
for (zz, nn), lit in LIT.items():
    row = {"literature_approx_eV": lit}
    for m in ("uni_gshm_qed", "uni_u35", "uni_pocket"):
        row[m] = final.predict(zz, nn, model=m)
    sh[f"{zz},{nn}"] = row
    print(zz, nn, {k: round(v, 2) for k, v in row.items()})
out["superheavy"] = sh
with open(os.path.join(RES, "uni_referee.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, indent=1)
