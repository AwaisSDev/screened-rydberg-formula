r"""Version-2 candidates of the Screened Rydberg formula (pre-registered in docs/preregistration_v2.md).

Extends the frozen v1 model (models/push_a/model.py, imported read-only) with optional terms:

  rel = "exp"  bounded relativistic penetration factor
               R = exp( r_c (Z alpha)^2 / n * (1 - Za/Zeff) ),   0 <= 1 - Za/Zeff < 1
               -> same first-order Fermi-Segre scaling as v1's 1 + r_c (Z alpha)^2/n (Zeff/Za - 1) near Zeff ~ Za,
               but R > 0 and exp(-|r_c|(Z alpha)^2/n) <= R <= exp(|r_c|(Z alpha)^2/n) for ANY parameters.
  q2 = True    second-order charge dependence inside the bounded remainder:
               T = sum_g tau_g nu_g + sum_c dtau_c nu_c + sum_g tau2_g nu_g / (Za + kappa)
               (vanishes as Za -> infinity, so the exact Z^2 and Z coefficients are kept; the bound is unchanged)
  kappa_l      per-target-l kappa (existing v1 option).
Everything else (sigma1, groups, classes, bound, Dirac factor, Hund term, mu, QED/FNS) is v1.
"""
import os
import sys

import numpy as np
from scipy.optimize import least_squares

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
for _p in (_ROOT, os.path.join(_ROOT, "models", "unified"), os.path.join(_ROOT, "models", "first_principles"),
           os.path.join(_ROOT, "models", "semi_empirical")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from push_a_import import PA  # noqa: E402  (frozen v1 code, read-only)
import umodel as U  # noqa: E402
import gshm  # noqa: E402
from common.atomdata import RYDBERG_EV as RY, ALPHA  # noqa: E402

GROUPS = ["same", "in", "core", "df", "out"]


def param_names(spec):
    names = PA.param_names(dict(spec, tie={}))
    if spec.get("q2"):
        names += ["tau2_" + g for g in GROUPS]
    return [n for n in names if n not in spec.get("tie", {})]


def _resolve(spec, theta):
    P = dict(zip(param_names(spec), theta))
    for _ in range(5):
        for k, v in spec.get("tie", {}).items():
            if v in P:
                P[k] = P[v]
    return P


def zeff(theta, spec, a):
    P = _resolve(spec, theta)
    Z = a["Z"]
    Za = np.maximum(a["Zq"], 1.0)
    N = Z - a["Zq"] + 1
    kap = abs(P["kappa"]) + 0.05
    if spec.get("kappa_l"):
        kap = kap + np.abs(np.array([0.0, P["kappa_p"], P["kappa_d"], P["kappa_f"]]))[a["l"]]
    s1 = a["sigma1"]
    G = PA._group_matrix(spec, a)
    T = G @ np.array([P["tau_" + g] for g in GROUPS])
    if spec.get("dev"):
        Cd = PA._dev_matrix(a)
        T = T + Cd[:, [PA.DEV_IDX[c] for c in spec["dev"]]] @ np.array([P["dtau_" + c] for c in spec["dev"]])
    if spec.get("q2"):
        T = T + (G @ np.array([P["tau2_" + g] for g in GROUPS])) / (Za + kap)
    h = np.where(T >= 0, np.maximum(N - 1 - s1, 1e-9), np.maximum(s1, 1e-9))
    u = np.abs(T) / (h * (Za + kap))
    D = np.sign(T) * h * u / (1.0 + u)
    return np.maximum(Z - s1 - D, 1e-3), P


def evaluate(theta, spec, a, return_Ze=False):
    Ze, P = zeff(theta, spec, a)
    Z = a["Z"]
    Za = np.maximum(a["Zq"], 1.0)
    n = a["n"]
    b = RY * (Ze / n) ** 2 * gshm.dirac_factor(n, a["j"], Ze)
    rel = spec.get("rel")
    if rel:
        r = np.array([P["rel_" + c] for c in PA.REL_CLASSES])[a["rc"]]
        if rel == "exp":
            b = b * np.exp(r * (Z * ALPHA) ** 2 / n * (1.0 - Za / Ze))
        elif rel == "fs":
            b = b * (1.0 + (Z * ALPHA) ** 2 * r * (Ze / Za - 1.0) / n)
        else:
            raise ValueError(rel)
    xl = np.array([0.0, P["x_p"], P["x_d"], P["x_f"]])[a["l"]]
    b = b + RY * xl * a["kink"] / n ** 2
    ie = b * a["mu"] - a["qedfns"] * (Ze / Z) ** 2
    return (ie, Ze) if return_Ze else ie


def auto_tie(spec, mask):
    s = PA.auto_tie(dict(spec, rel=("fs" if spec.get("rel") == "exp" else spec.get("rel"))), mask)
    s["rel"] = spec.get("rel")
    a = U.dataset()
    G = PA._group_matrix(s, {"C": a["C"][mask], "l": a["l"][mask]})
    if s.get("q2"):
        for j, g in enumerate(GROUPS):
            if G[:, j].sum() == 0:
                s["tie"]["tau2_" + g] = "tau2_same"
    return s


def fit(mask=None, spec=None):
    a = U.dataset()
    m = np.ones(a["nrec"], bool) if mask is None else mask
    spec = auto_tie(dict(spec), m)
    if spec.get("hier"):
        cnt = (PA._dev_matrix({"C": a["C"][m], "l": a["l"][m]}) > 0).sum(0)
        spec["dev"] = [c for c in PA.DEV_CLASSES if cnt[PA.DEV_IDX[c]] >= spec.get("min_rows", 1)]
    b = U.subset(a, m)
    y = b["y"]
    names = param_names(spec)
    th0 = np.array([5.0 if nm == "kappa" else (0.3 if nm.startswith("tau_") else 0.0) for nm in names])
    lo = np.array([0.0 if nm.startswith("kappa") else -np.inf for nm in names])
    hi = np.array([100.0 if nm.startswith("kappa") else np.inf for nm in names])
    th0 = np.clip(th0, lo + 1e-9, np.where(np.isfinite(hi), hi - 1e-9, np.inf))
    lam = float(spec.get("ridge", 0.0))
    w = np.sqrt(lam * len(y)) * np.array([1.0 if nm.startswith("dtau_") else 0.0 for nm in names])

    def res(th):
        r = np.log(np.maximum(evaluate(th, spec, b), 1e-6 * y) / y)
        return np.concatenate([r, w * th]) if lam > 0 else r
    r = least_squares(res, th0, x_scale="jac", method="trf", bounds=(lo, hi), max_nfev=2000)
    return spec, r.x
