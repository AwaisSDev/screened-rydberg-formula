"""Fitting utilities for the GSHM (least squares on log IE, held-out splits, metrics)."""
import copy
import os
import sys

import numpy as np
from scipy.optimize import least_squares

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import gshm  # noqa: E402
from common.atomdata import load_records  # noqa: E402

RECS = load_records()
ZARR = np.array([r["Z"] for r in RECS])
NARR = np.array([r["N"] for r in RECS])
YARR = np.array([r["IE_eV"] for r in RECS])

SPLITS = {
    "S1": ZARR % 5 == 0,           # test mask: interpolation in Z
    "S2": ZARR >= 55,              # test mask: extrapolation in Z
    "S3": NARR % 6 == 0,           # test mask: unseen isoelectronic sequences
}


def subset_arrays(a, recmask):
    """Restrict an assembled dataset to the records in recmask (bool over records)."""
    keep_rec = np.where(recmask)[0]
    newidx = -np.ones(a["nrec"], int)
    newidx[keep_rec] = np.arange(len(keep_rec))
    rowmask = recmask[a["idx"]]
    b = {}
    for k, v in a.items():
        if k == "nrec":
            continue
        if k == "_feats":
            b[k] = {kk: vv[recmask] for kk, vv in v.items()}
        elif isinstance(v, np.ndarray) and k.startswith("rec_") or k == "y":
            b[k] = v[recmask]
        elif isinstance(v, np.ndarray) and len(v) == len(a["idx"]):
            b[k] = v[rowmask]
        else:
            b[k] = v
    b["idx"] = newidx[a["idx"][rowmask]]
    b["nrec"] = len(keep_rec)
    return b


def residuals(theta, spec, a, y, lam=0.0, theta_prior=None, scale=None):
    ie = gshm.evaluate_model(theta, spec, a)
    r = np.log(np.maximum(ie, 1e-6 * y) / y)
    if lam > 0:
        r = np.concatenate([r, lam * (theta - theta_prior) / scale])
    return r


def fit(spec, recmask=None, theta0=None, form=None, loss="linear", lam=None, max_nfev=400,
        verbose=0):
    """lam: ridge strength pulling every parameter towards its Slater-like initial value
    (gshm.init_params); default = spec['ridge'] (0 if absent)."""
    if lam is None:
        lam = spec.get("ridge", 0.0)
    form = form or spec["form"]
    a = gshm.dataset(form)
    if recmask is not None:
        a = subset_arrays(a, recmask)
    y = a["y"]
    if theta0 is None:
        theta0 = gshm.init_params(spec)
    prior = gshm.init_params(spec)
    scale = np.ones_like(prior)
    lo, hi = bounds(spec)
    theta0 = np.clip(theta0, lo + 1e-9, hi - 1e-9)
    res = least_squares(residuals, theta0, args=(spec, a, y, lam, prior, scale), loss=loss,
                        f_scale=0.02, x_scale="jac", max_nfev=max_nfev, verbose=verbose,
                        method="trf", bounds=(lo, hi))
    return res.x, res


def bounds(spec):
    """Physical bounds: screening constants in [0, 1.3] (outer-electron classes [0, 1])."""
    names = gshm.param_names(spec)
    lo = np.full(len(names), -np.inf)
    hi = np.full(len(names), np.inf)
    for i, nm in enumerate(names):
        if nm.startswith("sig_"):
            lo[i], hi[i] = 0.0, (1.0 if "out" in nm else 1.3)
            if nm == "sig_core" and spec.get("core_bounds"):
                lo[i], hi[i] = spec["core_bounds"]
        elif nm == "kappa" or nm.startswith("kap_"):
            lo[i], hi[i] = 0.0, 100.0
    return lo, hi


def predict_all(spec, theta, form=None):
    a = gshm.dataset(form or spec["form"])
    return gshm.evaluate_model(theta, spec, a)


def metrics(pred, mask=None):
    y = YARR
    if mask is None:
        mask = np.ones(len(y), bool)
    ape = np.abs(pred - y) / y * 100
    neu = mask & (ZARR == NARR)
    return {
        "n": int(mask.sum()),
        "MAPE": float(ape[mask].mean()),
        "median": float(np.median(ape[mask])),
        "p90": float(np.percentile(ape[mask], 90)),
        "max": float(ape[mask].max()),
        "neutral_MAPE": float(ape[neu].mean()) if neu.any() else float("nan"),
        "n_neutral": int(neu.sum()),
        "H_like_MAPE": float(ape[mask & (NARR == 1)].mean()) if (mask & (NARR == 1)).any() else float("nan"),
    }


def auto_tie(spec, recmask):
    """Tie every parameter that has no support in the training records to its analogue
    (gshm.ANALOG). Returns a new spec (unchanged if everything is supported)."""
    a = gshm.dataset("1e")
    C = a["C"][recmask]
    l = a["l"][recmask]
    tie = dict(spec.get("tie", {}))
    for nm, rep in gshm.ANALOG.items():
        if nm.startswith("sig_"):
            c = nm[4:]
            if c in spec["sigma"] and C[:, gshm.CIDX[c]].sum() == 0:
                tie[nm] = rep
        elif nm == "rel_f" and spec["rel"] and spec.get("rel_pen", True) and not (l == 3).any():
            tie[nm] = rep
        elif nm == "x_f" and spec["exch"] and not (l == 3).any():
            tie[nm] = rep
    return with_spec(spec, tie=tie)


def cv_split(spec, split, theta0=None, loss="linear", lam=None, tie=True):
    test = SPLITS[split]
    if tie:
        spec2 = auto_tie(spec, ~test)
        if spec2.get("tie") != spec.get("tie") and theta0 is not None:
            P = dict(zip(gshm.param_names(spec), theta0))
            theta0 = np.array([P[n] for n in gshm.param_names(spec2)])
        spec = spec2
    th, _ = fit(spec, ~test, theta0=theta0, loss=loss, lam=lam)
    pred = predict_all(spec, th)
    return metrics(pred, test), th, pred


def with_spec(spec, **kw):
    s = copy.deepcopy(spec)
    s.update(kw)
    return s
