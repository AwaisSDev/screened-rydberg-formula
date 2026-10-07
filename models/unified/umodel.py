r"""Unified ionization-energy model (Track A exact leading orders + Track B fitted remainder).

    IE(Z,N) = { Ry Zeff^2/n^2 * D_{n,j}(Zeff) * [1 + r_{lj} (Z alpha)^2 (Zeff/Za - 1)/n]
                + Ry x_l K_l(k)/n^2 } * exp(b_d (Z alpha)^2 [l=d] + b_p3 (Z alpha)^2 [p3/2]) * mu(Z)
              - [dE_QED(Z,n) + dE_FNS(Z,n)] (Zeff/Z)^2          (ns electrons, n <= 2 only)

    Zeff = Z - sigma1(config) - sum_g tau_g nu_g / (Za + kappa),     Za = Z - N + 1

* sigma1      : EXACT first-order screening constant (Track A, zero parameters): Ry (Z-sigma1)^2/n^2
                reproduces the exact Z^2 (Bohr) and Z (first-order e-e repulsion) coefficients.
* tau_g/(Za+kappa): Track B's charge-dependent screening ("penetration/relaxation") term. As
                Z -> infinity at fixed N it is O(1/Z) in Zeff, i.e. it only changes the O(Z^0)
                coefficient dE2 and beyond: the fitted part never touches the exact leading orders.
                nu_g = number of electrons of screening group g (Track B classes, gshm.CLASSES).
* D_{n,j}     : exact Dirac/Schroedinger ratio for a point charge Zeff (Track B, parameter-free).
* r_{lj}      : fitted relativistic penetration coefficient (Track B form).
* x_l K_l(k)  : Hund exchange kink (Track B). Here it is O(Z^0) (no Zeff factor) so the exact O(Z)
                term is preserved; sigma1 already contains the first-order Hund exchange exactly.
* mu(Z)       : reduced-mass factor M/(M+m) (Track A); dE_QED (Yerokhin-Shabaev self-energy table
                + computed Uehling) and dE_FNS (finite nuclear size) from Track A, screened by
                (Zeff/Z)^2 (|psi(0)|^2 scaling).

Public: build(ZN, shells_list) -> feature dict; evaluate(theta, spec, a); fit(spec, mask);
        predict(Z, N, shells=None) -> eV; predict_many(ZN, shells_list=None) -> eV array.
"""
import json
import os
import sys
from functools import lru_cache

import numpy as np
from scipy.optimize import least_squares

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
for _p in (_ROOT, _HERE, os.path.join(_ROOT, "models", "first_principles"),
           os.path.join(_ROOT, "models", "semi_empirical")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from common.atomdata import RYDBERG_EV, ALPHA, HARTREE_EV, load_records  # noqa: E402
import gshm  # noqa: E402  (Track B)
import relativity as rel  # noqa: E402  (Track A)
import abinitio  # noqa: E402

RY = RYDBERG_EV
PARAMS_JSON = os.path.join(_ROOT, "results", "uni_params.json")
CLASSES = gshm.CLASSES
CIDX = gshm.CIDX

# Groupings of screening classes for the fitted remainder tau_g
GROUPINGS = {
    "3group": {"same": gshm.S1_GROUPS["same"], "near": gshm.S1_GROUPS["near"],
               "core": gshm.S1_GROUPS["core"]},
    "class": {c: [c] for c in CLASSES if c not in ("sn_out", "out_sp")},
}
# f-type classes tied to d analogues when unsupported by the training data (as in Track B)
TAU_ANALOG = {"same_f": "same_d", "n1_sp_f": "n1_sp_d", "n1_d_f": "n1_d_d", "f_near": "d_near",
              "n1_f_f": "n1_d_d", "out_f": "out_d"}


# ----------------------------------------------------------------------------- features
def _qed_fns(Z, n, l):
    """QED + finite-size shift (eV, >0 = less bound) of a hydrogenic ns level, n <= 2."""
    if l != 0 or n > 2:
        return 0.0
    # Z > 110: Track A now extrapolates F_SE/F_U quadratically beyond the YS15 table and uses
    # its A ~ 2.5 Z radius formula (flagged as extrapolation; affects only 1s/2s removal).
    fns = rel._FNS.get(f"{n},-1,{Z}")          # Track A's cached FNS shifts (read only)
    if fns is None:
        fns = rel.fns_shift(n, -1, Z)           # computed, never written to Track A's cache
    return (rel.qed_shift(Z, n, 0) + fns) * HARTREE_EV


def _mu(Z):
    _, Mm = rel.nuclear_params(Z)
    return Mm / (1.0 + Mm)


def build(ZN, shells_list=None):
    """Assemble configuration features for a list of (Z, N) (no IE data are read)."""
    ents, sig1, qf, mu = [], [], [], []
    for i, (Z, N) in enumerate(ZN):
        sh = None if shells_list is None else shells_list[i]
        rows, rem, rearr = gshm.record_rows(int(Z), int(N), sh, "1e")
        ents.append((rows, rem, rearr, int(Z), int(N)))
        sig1.append(abinitio.sigma1(Z, N, sh))
        qf.append(_qed_fns(int(Z), rem[0], rem[1]))
        mu.append(_mu(int(Z)))
    a = gshm.assemble(ents)
    a["sigma1"] = np.array(sig1)
    a["qedfns"] = np.array(qf)
    a["mu"] = np.array(mu)
    a["_feats"] = gshm.corr_features(a)
    return a


@lru_cache(maxsize=1)
def dataset():
    """Features for all 5847 NIST rows plus the targets y (used only for fitting/scoring)."""
    recs = load_records()
    a = build([(r["Z"], r["N"]) for r in recs])
    abinitio.save_cache()
    a["y"] = np.array([r["IE_eV"] for r in recs])
    return a


# ----------------------------------------------------------------------------- model
def default_spec(grouping="class", corr=("x2_d", "x2_p3"), **kw):
    s = {"grouping": grouping, "rel": True, "exch": True, "corr": list(corr),
         "qed": True, "tie": {}}
    s.update(kw)
    return s


def param_names(spec):
    g = GROUPINGS[spec["grouping"]]
    names = ["tau_" + k for k in g] + ["kappa"]
    if spec.get("kappa_l"):
        names += ["kappa_p", "kappa_d", "kappa_f"]
    if spec.get("S2"):
        names += ["tau2_" + k for k in GROUPINGS["3group"]]
    if spec["rel"]:
        names += ["rel_" + c for c in gshm.REL_CLASSES]
    if spec["exch"]:
        names += ["x_p", "x_d", "x_f"]
    names += ["b_" + t for t in spec["corr"]]
    return [n for n in names if n not in spec.get("tie", {})]


def init_params(spec):
    return np.array([7.0 if n == "kappa" else 0.0 for n in param_names(spec)])


def evaluate(theta, spec, a, return_Ze=False):
    """Vectorised IE (eV) for an assembled feature dict a (from build()/dataset())."""
    P = dict(zip(param_names(spec), theta))
    for k, v in spec.get("tie", {}).items():
        P[k] = P[v]
    g = GROUPINGS[spec["grouping"]]
    M = np.zeros((len(CLASSES), len(g)))
    for j, (k, cl) in enumerate(g.items()):
        for c in cl:
            M[CIDX[c], j] = 1.0
    tau = np.array([P["tau_" + k] for k in g])
    Za = np.maximum(a["Zq"], 1.0)
    kap = abs(P["kappa"]) + 0.05
    if spec.get("kappa_l"):
        kap = kap + np.array([0.0, P["kappa_p"], P["kappa_d"], P["kappa_f"]])[a["l"]]
    Z = a["Z"]
    Ze = Z - a["sigma1"] - ((a["C"] @ M) @ tau) / (Za + kap)
    if spec.get("S2"):
        g3 = GROUPINGS["3group"]
        M3 = np.zeros((len(CLASSES), 3))
        for j, (k, cl) in enumerate(g3.items()):
            for c in cl:
                M3[CIDX[c], j] = 1.0
        tau2 = np.array([P["tau2_" + k] for k in g3])
        Ze = Ze - ((a["C"] @ M3) @ tau2) / (Za + kap) ** 2
    Ze = np.maximum(Ze, 1e-3)
    n = a["n"]
    b = RY * (Ze / n) ** 2
    if spec["rel"]:
        r = np.array([P["rel_" + c] for c in gshm.REL_CLASSES])[a["rc"]]
        b = b * gshm.dirac_factor(n, a["j"], Ze) * (1.0 + (Z * ALPHA) ** 2 * r * (Ze / Za - 1.0) / n)
    if spec["exch"]:
        xl = np.array([0.0, P["x_p"], P["x_d"], P["x_f"]])[a["l"]]
        b = b + RY * xl * a["kink"] / n ** 2
    ie = np.bincount(a["idx"], weights=a["w"] * b, minlength=a["nrec"])
    if spec["corr"]:
        s = sum(P["b_" + t] * a["_feats"][t] for t in spec["corr"])
        ie = ie * np.exp(s)
    if spec.get("qed", True):
        ie = ie * a["mu"] - a["qedfns"] * (Ze / Z) ** 2
    return (ie, Ze) if return_Ze else ie


# ----------------------------------------------------------------------------- fitting
def subset(a, mask):
    import fitlib  # Track B helper
    b = fitlib.subset_arrays(a, mask)
    for k in ("sigma1", "qedfns", "mu"):
        b[k] = a[k][mask]
    return b


def auto_tie(spec, mask):
    """Tie unsupported f-type parameters to their d analogues (no f electrons in training)."""
    a = dataset()
    C = a["C"][mask]
    tie = {}
    if spec["grouping"] == "class":
        for c, rep in TAU_ANALOG.items():
            if C[:, CIDX[c]].sum() == 0:
                tie["tau_" + c] = "tau_" + rep
    if not (a["l"][mask] == 3).any():
        if spec["rel"]:
            tie["rel_f"] = "rel_d"
        if spec["exch"]:
            tie["x_f"] = "x_d"
    s = dict(spec)
    s["tie"] = tie
    return s


def fit(spec, mask=None, theta0=None, ridge=0.0):
    """Least squares on log(IE_pred/IE_NIST) over the records in mask (all if None).
    ridge > 0 adds ridge*theta (pull towards 0 = pure ab-initio sigma1; kappa excluded)."""
    a = dataset()
    if mask is not None:
        a = subset(a, mask)
    y = a["y"]
    th0 = init_params(spec) if theta0 is None else theta0
    names = param_names(spec)
    lo = np.array([0.0 if n == "kappa" else -np.inf for n in names])
    hi = np.array([100.0 if n == "kappa" else np.inf for n in names])
    w = np.array([0.0 if n == "kappa" else ridge for n in names])

    def res(th):
        r = np.log(np.maximum(evaluate(th, spec, a), 1e-6 * y) / y)
        return np.concatenate([r, w * th]) if ridge > 0 else r
    r = least_squares(res, th0, x_scale="jac", method="trf", bounds=(lo, hi), max_nfev=600)
    return r.x, r


# ----------------------------------------------------------------------------- public API
@lru_cache(maxsize=1)
def _load_params(path=PARAMS_JSON):
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    return d["spec"], np.array(d["values"], float)


def predict_many(ZN, shells_list=None, params=None):
    spec, theta = params if params is not None else _load_params()
    return evaluate(theta, spec, build(ZN, shells_list))


def predict(Z, N, shells=None, params=None):
    """Unified-model ionization energy (eV) of the N-electron ion of element Z.
    shells: optional configuration [(n,l,occ),...]; default NIST ground (Madelung outside table)."""
    return float(predict_many([(int(Z), int(N))], None if shells is None else [shells], params)[0])
