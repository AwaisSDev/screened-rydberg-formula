r"""Push A: robust / minimal-remainder Screened Rydberg formula.

    IE(Z,N) = { Ry Zeff^2/n^2 * D_{n,j}(Zeff) * R + Ry x_l K_l(k)/n^2 } * mu(Z)
              - [dE_QED + dE_FNS](Z,n) * (Zeff/Z)^2                       (ns, n <= 2 only)

    Zeff = Z - sigma1 - D,      sigma1 = EXACT first-order (1/Z) screening constant (Track A)

    Bounded higher-order screening (the only fitted screening term):
        T = sum_g tau_g nu_g                     (nu_g = electrons in screening group g)
        D = T / (Za + kappa_l + |T| / h),        h = g_max = (N-1) - sigma1   if T >= 0
                                                 h = sigma1                    if T <  0
    so that, for ANY parameter values,  Za <= Zeff <= Z   (Za = Z - N + 1):
        * Zeff >= Za is the non-penetrating hydrogenic limit (IE >= Ry Za^2/n^2 before the
          relativistic factor): penetration can only increase binding;
        * Zeff <= Z: screening cannot be negative overall.
    At large Za (fixed N) D -> T/(Za + kappa_l) = O(1/Z): sigma -> sigma1 exactly, i.e. the
    fitted part only changes the O(Z^0) coefficient and beyond (exact Z^2 and Z terms kept).
    The bound is smooth (rational saturation, no hard clip): D/h = u/(1+u), u = |T|/(h(Za+kappa)).

* D_{n,j}(Zeff): exact Dirac/Schroedinger ratio for a point charge Zeff (no parameters).
* R (relativistic penetration, optional): 1 + r_c (Z alpha)^2 B /n, with r_c per relativistic class
  (s, p1/2, p3/2, d, f) and B a bounded penetration measure (see REL_FORMS).
* x_l K_l(k): Hund exchange kink of the removed electron (Track B), O(Z^0).
* mu(Z), dE_QED, dE_FNS: Track A recoil / QED / finite-size layer (parameter-free).

predict() never reads NIST ionization energies; parameters come from results/pa_params.json.
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

import umodel as U  # noqa: E402  (read-only import: features, sigma1, QED layer)
import gshm  # noqa: E402
from common.atomdata import RYDBERG_EV, ALPHA  # noqa: E402

RY = RYDBERG_EV
CLASSES, CIDX = gshm.CLASSES, gshm.CIDX
PARAMS_JSON = os.path.join(_ROOT, "results", "pa_params.json")

# screening groups (nu_g) -------------------------------------------------------------------
_SAME = [c for c in CLASSES if c.startswith("same_")]
GROUPINGS = {
    # one universal coefficient: T = tau * (N - 1)
    "one": {"all": list(CLASSES)},
    # same subshell / near shells / deep core / outer (n' > n)
    "4group": {"same": _SAME,
               "near": ["sn_in_p", "n1_sp_sp", "n1_sp_d", "n1_sp_f", "d_near", "n1_d_d", "n1_d_f",
                        "f_near", "n1_f_f"],
               "core": ["n2_sp_sp", "n2_sp_df", "core"],
               "out": ["sn_out", "out_sp", "out_d", "out_f"]},
    # pocket-like groups, split by target type (s/p vs d/f)
    "pocket": {"same": _SAME,
               "in": ["sn_in_p", "n1_sp_sp", "n1_sp_d", "n1_sp_f"],
               "core": ["n2_sp_sp", "n2_sp_df"],          # + 'core' class for s,p targets
               "df": ["d_near", "n1_d_d", "n1_d_f", "f_near", "n1_f_f"],   # + 'core' for d,f
               "out": ["sn_out", "out_sp", "out_d", "out_f"]},
    # Track B classes (f analogues tied to d when unsupported)
    "class": {c: [c] for c in CLASSES if c not in ("sn_out", "out_sp")},
}
TAU_ANALOG = U.TAU_ANALOG
REL_CLASSES = gshm.REL_CLASSES          # s, p1, p3, d, f


def default_spec(**kw):
    s = {"grouping": "4group", "kappa_l": False, "rel": "fs", "exch": True, "qed": True,
         "bound": True, "tie": {}}
    s.update(kw)
    return s


def param_names(spec):
    g = GROUPINGS[spec["grouping"]]
    names = ["tau_" + k for k in g] + ["kappa"]
    names += ["dtau_" + c for c in spec.get("dev", [])]
    if spec.get("kappa_l"):
        names += ["kappa_p", "kappa_d", "kappa_f"]
    if spec.get("rel") == "fs1":
        names += ["lam_sp"]
    elif spec.get("rel") == "fsp":
        names += ["lam_df"]
    elif spec.get("rel") and spec.get("rel") != "fs0":
        names += ["rel_" + c for c in REL_CLASSES]
    if spec.get("exch"):
        names += ["x_p", "x_d", "x_f"]
    return [n for n in names if n not in spec.get("tie", {})]


def init_params(spec):
    v = []
    for n in param_names(spec):
        v.append(5.0 if n == "kappa" else (0.3 if n.startswith("tau_") else (1.0 if n == "lam_sp" else 0.0)))
    return np.array(v)


# class-level deviations (hierarchical model): Track B classes with 'core' split by target type
DEV_CLASSES = [c for c in CLASSES if c != "core"] + ["core_sp", "core_df"]
DEV_IDX = {c: i for i, c in enumerate(DEV_CLASSES)}


def _dev_matrix(a):
    C = a["C"]
    sp = a["l"] <= 1
    cols = [C[:, CIDX[c]] for c in CLASSES if c != "core"]
    cols += [np.where(sp, C[:, CIDX["core"]], 0), np.where(sp, 0, C[:, CIDX["core"]])]
    return np.stack(cols, 1)


def _group_matrix(spec, a):
    """Electron counts per screening group, shape (nrec, ngroups)."""
    g = GROUPINGS[spec["grouping"]]
    C = a["C"]
    sp = a["l"] <= 1
    cols = []
    for k, cl in g.items():
        v = C[:, [CIDX[c] for c in cl]].sum(1)
        if spec["grouping"] == "pocket" and k == "core":
            v = v + np.where(sp, C[:, CIDX["core"]], 0)
        if spec["grouping"] == "pocket" and k == "df":
            v = v + np.where(sp, 0, C[:, CIDX["core"]])
        cols.append(v)
    return np.stack(cols, 1)


def _resolve(spec, theta):
    P = dict(zip(param_names(spec), theta))
    tie = spec.get("tie", {})
    for _ in range(5):                      # resolve tie chains (f -> d -> fallback)
        for k, v in tie.items():
            if v in P:
                P[k] = P[v]
    return P


def zeff(theta, spec, a):
    P = _resolve(spec, theta)
    g = GROUPINGS[spec["grouping"]]
    tau = np.array([P["tau_" + k] for k in g])
    Z = a["Z"]
    Za = np.maximum(a["Zq"], 1.0)
    N = Z - a["Zq"] + 1
    kap = abs(P["kappa"]) + 0.05
    if spec.get("kappa_l"):
        kap = kap + np.abs(np.array([0.0, P["kappa_p"], P["kappa_d"], P["kappa_f"]]))[a["l"]]
    s1 = a["sigma1"]
    T = _group_matrix(spec, a) @ tau
    if spec.get("dev"):
        Cd = _dev_matrix(a)
        T = T + Cd[:, [DEV_IDX[c] for c in spec["dev"]]] @ np.array([P["dtau_" + c] for c in spec["dev"]])
    if spec.get("bound", True):
        h = np.where(T >= 0, np.maximum(N - 1 - s1, 1e-9), np.maximum(s1, 1e-9))
        u = np.abs(T) / (h * (Za + kap))
        sat = spec.get("sat", "rat")
        if sat == "rat":
            phi = u / (1.0 + u)
        elif sat == "exp":
            phi = -np.expm1(-u)
        elif sat == "rat2":
            phi = u / np.sqrt(1.0 + u * u)
        elif sat == "tanh":
            phi = np.tanh(u)
        D = np.sign(T) * h * phi
    else:
        D = T / (Za + kap)
    Ze = Z - s1 - D
    return np.maximum(Ze, 1e-3), P


def rel_measure(form, Ze, Za, Z, N):
    """Bounded penetration measure B used in the relativistic factor 1 + r (Z alpha)^2 B / n."""
    if form == "fs":          # Fermi-Segre: dE_rel/E ~ (Z alpha)^2 Zeff/(n Za); hydrogenic part removed
        return Ze / Za - 1.0
    if form == "frac":        # penetration fraction (Zeff - Za)/(Z - Za) in [0, 1]
        return np.where(N > 1, (Ze - Za) / np.maximum(Z - Za, 1e-9), 0.0)
    if form == "sat":         # Fermi-Segre measure saturated: x/(1+x/ (n?)) -> bounded by Z/Za - 1 anyway
        x = Ze / Za - 1.0
        return x / (1.0 + 0.25 * x)
    raise ValueError(form)


def _rel_lambda(rel, P, a):
    sp = (a["l"] <= 1).astype(float)
    if rel == "fs0":                 # parameter-free, s and p only (d, f: no direct penetration)
        return sp
    if rel == "fs1":                 # one fitted scale for s, p
        return P["lam_sp"] * sp
    if rel == "fsp":                 # parameter-free s,p + fitted d/f coefficient
        return sp + P["lam_df"] * (1.0 - sp)
    return np.array([P["rel_" + c] for c in REL_CLASSES])[a["rc"]]     # fsc: per class


def evaluate(theta, spec, a, return_Ze=False):
    Ze, P = zeff(theta, spec, a)
    Z = a["Z"]
    Za = np.maximum(a["Zq"], 1.0)
    N = Z - a["Zq"] + 1
    n = a["n"]
    b = RY * (Ze / n) ** 2 * gshm.dirac_factor(n, a["j"], Ze)
    rel = spec.get("rel")
    if rel in ("fs0", "fs1", "fsc", "fsp"):
        # Fermi-Segre penetration: near-nucleus relativistic energy ~ Ry a^2 Z^2 Za^2/(n*^3 (j+1/2)),
        # n* = n Za/Zeff, minus the part already in the hydrogenic Dirac factor of charge Zeff.
        lam = _rel_lambda(rel, P, a)
        b = b * (1.0 + lam * ALPHA ** 2 * Ze * (Z ** 2 / Za - Ze) / (n * (a["j"] + 0.5)))
    elif rel:
        r = np.array([P["rel_" + c] for c in REL_CLASSES])[a["rc"]]
        b = b * (1.0 + (Z * ALPHA) ** 2 * r * rel_measure(rel, Ze, Za, Z, N) / n)
    if spec.get("exch"):
        xl = np.array([0.0, P["x_p"], P["x_d"], P["x_f"]])[a["l"]]
        b = b + RY * xl * a["kink"] / n ** 2
    ie = b
    if spec.get("qed", True):
        ie = ie * a["mu"] - a["qedfns"] * (Ze / Z) ** 2
    return (ie, Ze) if return_Ze else ie


# ----------------------------------------------------------------------------- fitting
def auto_tie(spec, mask):
    """Tie parameters that the training rows cannot determine to their d analogues."""
    a = U.dataset()
    C = a["C"][mask]
    tie = {}
    if spec["grouping"] == "class":
        for c, rep in TAU_ANALOG.items():
            if C[:, CIDX[c]].sum() == 0:
                tie["tau_" + c] = "tau_" + rep
        fallback = {"n1_d_d": "d_near", "n2_sp_df": "n2_sp_sp", "n1_sp_d": "n1_sp_sp",
                    "out_d": "d_near", "d_near": "n1_sp_sp", "same_d": "same_p", "core": "n2_sp_sp"}
        for c, rep in fallback.items():
            if "tau_" + c not in tie and C[:, CIDX[c]].sum() == 0:
                tie["tau_" + c] = "tau_" + rep
    else:
        G = _group_matrix(spec, {"C": C, "l": a["l"][mask]})
        names = list(GROUPINGS[spec["grouping"]])
        for j, k in enumerate(names):
            if G[:, j].sum() == 0:
                tie["tau_" + k] = "tau_" + names[0]
    if not (a["l"][mask] == 3).any():
        if spec.get("rel") in ("fs", "frac", "sat", "fsc"):
            tie["rel_f"] = "rel_d"
        if spec.get("exch"):
            tie["x_f"] = "x_d"
        if spec.get("kappa_l"):
            tie["kappa_f"] = "kappa_d"
    s = dict(spec)
    s["tie"] = tie
    return s


def fit(mask=None, spec=None, theta0=None):
    """Least squares on log(IE_pred/IE_NIST) over the rows in mask (all if None).
    Returns (spec_with_ties, theta)."""
    a = U.dataset()
    spec = default_spec() if spec is None else spec
    m = np.ones(a["nrec"], bool) if mask is None else mask
    spec = auto_tie(spec, m)
    if spec.get("hier"):
        cnt = (_dev_matrix({"C": a["C"][m], "l": a["l"][m]}) > 0).sum(0)
        spec["dev"] = [c for c in DEV_CLASSES if cnt[DEV_IDX[c]] >= spec.get("min_rows", 1)]
    b = U.subset(a, m)
    y = b["y"]
    names = param_names(spec)
    th0 = init_params(spec) if theta0 is None else theta0
    lo = np.array([0.0 if n.startswith("kappa") else -np.inf for n in names])
    hi = np.array([100.0 if n.startswith("kappa") else np.inf for n in names])
    th0 = np.clip(th0, lo + 1e-9, np.where(np.isfinite(hi), hi - 1e-9, np.inf))
    lam = float(spec.get("ridge", 0.0))
    w = np.sqrt(lam * len(y)) * np.array([1.0 if nm.startswith("dtau_") else 0.0 for nm in names])

    def res(th):
        r = np.log(np.maximum(evaluate(th, spec, b), 1e-6 * y) / y)
        return np.concatenate([r, w * th]) if lam > 0 else r
    r = least_squares(res, th0,
                      x_scale="jac", method="trf", bounds=(lo, hi), max_nfev=2000)
    return spec, r.x


# ----------------------------------------------------------------------------- public API
@lru_cache(maxsize=1)
def _load_params(path=PARAMS_JSON):
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    c = d["candidates"][d["chosen"]]
    return c["spec"], np.array(c["values"], float)


def predict_many(ZN, shells_list=None, params=None):
    spec, theta = params if params is not None else _load_params()
    return evaluate(theta, spec, U.build(ZN, shells_list))


def predict(Z, N, shells=None, params=None):
    """Ionization energy (eV) of the N-electron ion of element Z (any Z <= 120, 1 <= N <= Z).
    shells: optional configuration [(n,l,occ),...]; default NIST ground / Madelung."""
    return float(predict_many([(int(Z), int(N))], None if shells is None else [shells], params)[0])
