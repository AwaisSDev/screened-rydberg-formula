r"""Push-B: physics-constrained "Screened Rydberg" model (exact sigma1 + bounded, saturating remainder).

    IE(Z,N) = { Ry Zeff^2/n^2 * D_{n,j}(Zeff) * R_{lj} + Ry x_l K_l(k)/n^2 } * exp(sum_t b_t phi_t) * mu(Z)
              - [dE_QED + dE_FNS](Z,n) (Zeff/Z)^2                     (ns, n <= 2 only; Track A)

Effective charge (form "logistic"; Za = Z - N + 1 is the asymptotic charge):

    D1   = Z - sigma1 - Za  = (N - 1) - sigma1         >= 0   exact first-order "penetration excess"
    X    = sum_c nu_c t_c / (Za + kappa_l)                     charge-dependent relaxation/correlation
    Zeff = Za + D1 * 2 / (1 + exp(X))

  * sigma1: exact first-order screening (Track A, zero parameters).  For X -> 0 (Za -> infinity at
    fixed N) Zeff = Z - sigma1 - D1 X/2 + O(X^3): the exact Z^2 and Z coefficients are preserved and
    the fitted part is O(1/Z), exactly as in u35.
  * Physical bounds hold for ANY parameter values and ANY electron count:  Za <= Zeff <= Za + 2 D1 <= Z
    (D1 <= 0.375 (N-1) for every configuration in the table), i.e. IE >= Ry Za^2/n^2 (non-penetrating
    hydrogenic limit) and Zeff never exceeds the bare charge.  Heavy many-electron atoms (large
    sum nu t) saturate to the screened limit instead of running away.
  * form "exp":   Zeff = Za + D1 exp(-X)   (same bounds when t >= 0).
  * form "clip":  u35's linear Zeff_raw = Z - sigma1 - X, mapped smoothly into [Za, Z].
  * R_{lj}: relativistic penetration factor.  rel_form "sat": exp[(Z alpha)^2 r_{lj} (1 - Za/Zeff)/n]
    (bounded, always positive); rel_form "u35": u35's bracket 1 + (Z alpha)^2 r (Zeff/Za - 1)/n
    (unbounded, can turn negative far outside the training range).
  * f-type parameters are tied to their d analogues when unsupported by the training rows.

Public: evaluate(theta, spec, a); fit(spec, mask, ridge); predict(Z, N, shells=None);
        predict_many(ZN, shells_list=None).  Parameters: results/pb_params.json.  predict() never
        reads NIST IE values (features are built from Z and the configuration only).
"""
import json
import os
import sys
from functools import lru_cache

import numpy as np
from scipy.optimize import least_squares

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
for _p in (_ROOT, os.path.join(_ROOT, "models", "unified"), os.path.join(_ROOT, "models", "first_principles"),
           os.path.join(_ROOT, "models", "semi_empirical")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import umodel as U  # noqa: E402  (read-only import)
import gshm  # noqa: E402
from common.atomdata import ALPHA  # noqa: E402

RY = U.RY
PARAMS_JSON = os.path.join(_ROOT, "results", "pb_params.json")
TCLASSES = list(U.GROUPINGS["class"].keys())          # 18 screening classes (sn_out/out_sp never occur)
CIDX = U.CIDX
# a-priori analogues (deep electrons screen alike regardless of l; inner-shell electrons of a d target)
EXTRA_ANALOG = {"n2_sp_df": "n2_sp_sp", "n1_d_d": "d_near"}
_M = np.zeros((len(U.CLASSES), len(TCLASSES)))
for _j, _c in enumerate(TCLASSES):
    _M[CIDX[_c], _j] = 1.0


def default_spec(**kw):
    s = {"form": "logistic", "kappa_l": True, "rel_form": "sat", "corr": [], "tie": {}, "tpos": False}
    s.update(kw)
    return s


def param_names(spec):
    names = ["t_" + c for c in TCLASSES] + ["kappa"]
    if spec.get("kappa_l"):
        names += ["kappa_p", "kappa_d", "kappa_f"]
    names += ["rel_" + c for c in gshm.REL_CLASSES]
    names += ["x_p", "x_d", "x_f"]
    names += ["b_" + t for t in spec.get("corr", [])]
    return [n for n in names if n not in spec.get("tie", {})]


def init_params(spec):
    return np.array([7.0 if n == "kappa" else (0.2 if n.startswith("t_") else 0.0)
                     for n in param_names(spec)])


def _softclip01(u, beta=12.0):
    """Smooth map R -> (0,1), ~identity on [0.15, 0.85]."""
    return (np.logaddexp(0, beta * u) - np.logaddexp(0, beta * (u - 1.0))) / beta


def zeff(P, spec, a):
    t = np.array([P["t_" + c] for c in TCLASSES])
    Z = a["Z"]
    Za = np.maximum(a["Zq"], 1.0)
    kap = abs(P["kappa"]) + 0.05
    if spec.get("kappa_l"):
        kap = kap + np.abs(np.array([0.0, P["kappa_p"], P["kappa_d"], P["kappa_f"]]))[a["l"]]
    X = ((a["C"] @ _M) @ t) / (Za + kap)
    D1 = np.maximum(Z - a["sigma1"] - Za, 0.0)
    form = spec["form"]
    if form == "logistic":
        return Za + D1 * 2.0 / (1.0 + np.exp(np.clip(X, -50, 50)))
    if form == "exp":
        return Za + D1 * np.exp(-np.clip(X, -50, 50))
    if form == "clip":
        raw = Z - a["sigma1"] - X
        span = Z - Za
        u = np.where(span > 0, (raw - Za) / np.where(span > 0, span, 1.0), 1.0)
        return np.where(span > 0, Za + span * _softclip01(u), Z)
    raise ValueError(form)


def evaluate(theta, spec, a, return_Ze=False):
    P = dict(zip(param_names(spec), theta))
    for k, v in spec.get("tie", {}).items():
        P[k] = v if isinstance(v, float) else P[v]
    Ze = zeff(P, spec, a)
    Z = a["Z"]
    Za = np.maximum(a["Zq"], 1.0)
    n = a["n"]
    r = np.array([P["rel_" + c] for c in gshm.REL_CLASSES])[a["rc"]]
    if spec.get("rel_form", "sat") == "sat":
        R = np.exp(np.clip((Z * ALPHA) ** 2 * r * (1.0 - Za / Ze) / n, -50, 50))
    else:
        R = 1.0 + (Z * ALPHA) ** 2 * r * (Ze / Za - 1.0) / n
    b = RY * (Ze / n) ** 2 * gshm.dirac_factor(n, a["j"], Ze) * R
    xl = np.array([0.0, P["x_p"], P["x_d"], P["x_f"]])[a["l"]]
    b = b + RY * xl * a["kink"] / n ** 2
    ie = np.bincount(a["idx"], weights=a["w"] * b, minlength=a["nrec"])
    if spec.get("corr"):
        ie = ie * np.exp(sum(P["b_" + t] * a["_feats"][t] for t in spec["corr"]))
    ie = ie * a["mu"] - a["qedfns"] * (Ze / Z) ** 2
    return (ie, Ze) if return_Ze else ie


def auto_tie(spec, mask):
    """Tie f-type parameters to their d analogues when the training rows contain no such electrons."""
    a = U.dataset()
    C = a["C"][mask]
    tie = {}
    for c, rep in U.TAU_ANALOG.items():
        if C[:, CIDX[c]].sum() == 0:
            tie["t_" + c] = "t_" + rep
    if spec.get("tie_more", True):
        # a-priori analogues for non-f classes absent from the training rows; anything still
        # unsupported is fixed at 0 (= pure ab-initio sigma1 for that class)
        for c, rep in EXTRA_ANALOG.items():
            if C[:, CIDX[c]].sum() == 0 and C[:, CIDX[rep]].sum() > 0:
                tie["t_" + c] = "t_" + rep
        for c in TCLASSES:
            if C[:, CIDX[c]].sum() == 0 and "t_" + c not in tie:
                tie["t_" + c] = 0.0
            elif "t_" + c in tie and C[:, CIDX[tie["t_" + c][2:]]].sum() == 0:
                tie["t_" + c] = 0.0
    if not (a["l"][mask] == 3).any():
        tie["rel_f"] = "rel_d"
        tie["x_f"] = "x_d"
    s = dict(spec)
    s["tie"] = tie
    return s


def fit(spec, mask=None, ridge=0.0, theta0=None):
    """Least squares on log(IE_pred/IE_NIST) over the rows in mask.  ridge pulls t_c, r, x, b toward 0
    (t = 0 is the pure ab-initio sigma1 model); kappa's are not penalised."""
    a = U.dataset()
    if mask is None:
        mask = np.ones(a["nrec"], bool)
    spec = auto_tie(spec, mask)
    b = U.subset(a, mask)
    y = b["y"]
    names = param_names(spec)
    th0 = init_params(spec) if theta0 is None else np.asarray(theta0, float)
    lo = np.array([0.0 if (n.startswith("kappa") or (spec.get("tpos") and n.startswith("t_"))) else -np.inf
                   for n in names])
    hi = np.array([100.0 if n.startswith("kappa") else np.inf for n in names])
    th0 = np.minimum(np.maximum(th0, lo + 1e-6), hi - 1e-6)
    w = np.array([0.0 if n.startswith("kappa") else ridge for n in names])
    sq = np.sqrt(len(y))

    def res(th):
        r = np.log(np.maximum(evaluate(th, spec, b), 1e-6 * y) / y)
        return np.concatenate([r, sq * w * th]) if ridge > 0 else r
    r = least_squares(res, th0, x_scale="jac", method="trf", bounds=(lo, hi), max_nfev=800)
    return spec, r.x


# ----------------------------------------------------------------------------- public API
@lru_cache(maxsize=1)
def _load_params(path=PARAMS_JSON):
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    return d["spec"], np.array(d["values"], float)


def predict_many(ZN, shells_list=None, params=None):
    spec, theta = params if params is not None else _load_params()
    return evaluate(theta, spec, U.build(ZN, shells_list))


def predict(Z, N, shells=None, params=None):
    """IE (eV) of the N-electron ion of element Z (Z <= 120, 1 <= N <= Z); shells optional."""
    return float(predict_many([(int(Z), int(N))], None if shells is None else [shells], params)[0])
