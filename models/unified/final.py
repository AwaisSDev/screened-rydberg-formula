"""Final unified predictor ("Screened Rydberg formula"). Default model: pa_hier_rel.

Models (select with model=...):
  pa_hier_rel  [DEFAULT] exact Track A sigma1 + bounded hierarchical screening remainder + Dirac factor,
               relativistic penetration, Hund kink, recoil/QED/FNS (33 fitted params). Chosen by the
               pre-registered score (mean MAPE of V1, V2, S1, S3), never by S2. Code: models/push_a/model.py
               (imported read-only); parameters: results/uni_final_params.json (all-data fit written by
               models/unified/validate_blind.py; identical to results/pa_params.json).
  uni_u35      previous unified final: exact sigma1 + linear remainder (35 params), umodel.py
  uni_u29      previous unified variant (29 params), umodel.py
  uni_pocket   8-parameter pocket formula, pocket.py
  uni_gshm_qed Track B GSHM (32 params) x Track A reduced mass, minus QED + FNS of ns (n<=2) electrons
Aliases: "final" -> pa_hier_rel, "u35" -> uni_u35, "u29", "pocket", "gshm_qed".
No model reads NIST ionization energies at prediction time.
"""
import json
import os
import sys
import warnings
from functools import lru_cache

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _HERE)
import umodel as U  # noqa: E402
import gshm  # noqa: E402
from push_a_import import PA  # noqa: E402

PARAMS = os.path.join(_ROOT, "results", "uni_params.json")
FINAL_PARAMS = os.path.join(_ROOT, "results", "uni_final_params.json")
DEFAULT_MODEL = "pa_hier_rel"
MODELS = ("pa_hier_rel", "uni_u35", "uni_u29", "uni_pocket", "uni_gshm_qed")
ALIASES = {"final": "pa_hier_rel", "u35": "uni_u35", "u29": "uni_u29", "pocket": "uni_pocket",
           "gshm_qed": "uni_gshm_qed", "uni_final": "pa_hier_rel"}


@lru_cache(maxsize=1)
def params():
    """Legacy parameter file (u35, u29, pocket). 'model' names the default model."""
    with open(PARAMS, encoding="utf-8") as f:
        d = json.load(f)
    d["model"] = DEFAULT_MODEL
    return d


@lru_cache(maxsize=1)
def final_params():
    path = FINAL_PARAMS
    if not os.path.exists(path):                        # fallback: Push A's own all-data fit
        with open(os.path.join(_ROOT, "results", "pa_params.json"), encoding="utf-8") as f:
            c = json.load(f)["candidates"]["pa_hier_rel"]
        return c["spec"], np.array(c["values"], float)
    with open(path, encoding="utf-8") as f:
        c = json.load(f)
    return c["spec"], np.array(c["values"], float)


def resolve(model):
    model = ALIASES.get(model, model) if model else DEFAULT_MODEL
    if model not in MODELS:
        raise ValueError(f"unknown model {model!r}; choose from {MODELS} or aliases {sorted(ALIASES)}")
    return model


def predict_many(ZN, shells_list=None, model=None):
    model = resolve(model)
    a = U.build(ZN, shells_list)
    if model == "pa_hier_rel":
        spec, theta = final_params()
        # The bound Za <= Zeff <= Z needs 0 <= sigma1 <= N-1 (checked for all 7021 ground/Madelung
        # configurations with Z <= 118, 0 violations); warn if a user-supplied configuration breaks it.
        nel = a["Z"] - a["Zq"] + 1
        bad = (a["sigma1"] < -1e-9) | (a["sigma1"] > nel - 1 + 1e-9)
        if np.any(bad):
            warnings.warn(f"{int(bad.sum())} configuration(s) have sigma1 outside [0, N-1]; "
                          "the Za <= Zeff <= Z bound is not guaranteed for them", RuntimeWarning)
        return PA.evaluate(theta, spec, a)
    P = params()
    if model == "uni_gshm_qed":
        with open(os.path.join(_ROOT, "results", "se_params.json"), encoding="utf-8") as f:
            b = json.load(f)["final"]
        ie, parts = gshm.evaluate_model(np.array(b["values"]), b["spec"], a, return_parts=True)
        Ze = np.minimum(parts["Ze"], a["Z"])
        return ie * a["mu"] - a["qedfns"] * (Ze / a["Z"]) ** 2
    if model == "uni_pocket":
        import pocket
        return pocket.evaluate(np.array(P["pocket"]["values"]), a)
    c = P["all_candidates"][model]
    return U.evaluate(np.array(c["values"]), c["spec"], a)


def predict(Z, N, shells=None, model=None):
    """Ionization energy (eV) of the N-electron ion of element Z (any Z <= 120, 1 <= N <= Z)."""
    return float(predict_many([(int(Z), int(N))], None if shells is None else [shells], model)[0])
