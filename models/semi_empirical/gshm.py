r"""Generalized Screened-Hydrogenic Model (GSHM) for successive ionization energies.

Track B (data-driven / pattern recognition).  Public API:

    predict(Z, N, shells=None, params=None) -> IE in eV

The model (see docs/semi_empirical.md for the full derivation):

    IE = Ry * (Z_eff / n)^2 * F_rel  +  E_x                           (one-electron form)
    Z_eff = Z - S0 - S1 / (Z - S0 + kappa)                           (screening + penetration)
    S0    = sum_c sigma_c * c_c                                     (c_c = # electrons of class c)
    S1    = sum_g tau_g * c_g                                       (q-dependent screening)
    F_rel = D(n, j, Z_eff) * [1 + (Z alpha)^2 * r_{l,j} * (Z_eff/Z_a - 1)/n]
    E_x   = Ry * x_l * K(l, k) * Z_eff / n^2                        (Hund exchange kink)

followed by a small sparse multiplicative correction exp(sum_m beta_m phi_m) whose terms were
selected by cross-validated forward selection (fit_gshm.py).

Every quantity is computed from (Z, N, configuration) only. The NIST IE values are never
read inside predict().
"""
import json
import os
import sys
from functools import lru_cache

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from common.atomdata import RYDBERG_EV, ALPHA, HARTREE_EV, load_records  # noqa: E402
from configs import initial_final, exchange_kink, jj_j, normalize  # noqa: E402

RY = RYDBERG_EV
MC2_EV = 2.0 * RY / ALPHA ** 2           # electron rest energy in eV (511 keV)
LL = "spdf"
PARAMS_JSON = os.path.join(_ROOT, "results", "se_params.json")

# ----------------------------------------------------------------------------- screening classes


def _grp(l):
    return "sp" if l <= 1 else LL[l]


def classify(n, l, n2, l2):
    """Screening class of an electron in (n2,l2) acting on the electron in (n,l).

    Classes (identifiable on ground-state data; full inner s,p shells of d/f targets are
    merged because their populations are always the same constant):
      same_<l>      other electrons of the same subshell
      sn_in_p       p target, s electrons of the same n
      n1_sp_<g>     s/p target, electrons of shell n-1 with l-group g in {sp, d, f}
      n2_sp_<g>     s/p target, electrons of shell n-2, g in {sp, df}
      d_near        d target, s,p electrons of shells n and n-1
      n1_d_d        d target, (n-1)d electrons
      n1_d_f        d target, (n-1)f electrons (4f for 5d)
      f_near        f target, s,p,d electrons of shells n and n-1
      n1_f_f        f target, (n-1)f electrons (4f for 5f)
      core          everything deeper (s/p target: n' <= n-3; d/f target: n' <= n-2)
      sn_out        same n, higher l (only appears in the total-energy form)
      out_<g>       electrons in shells with n' > n (4s acting on 3d, 5s5p6s on 4f, ...)
    """
    if (n2, l2) == (n, l):
        return "same_" + LL[l]
    dn = n - n2
    if dn == 0 and l2 > l:
        return "sn_out"
    if dn < 0:
        return "out_" + _grp(l)
    if l <= 1:
        if dn == 0:
            return "sn_in_p"
        if dn == 1:
            return "n1_sp_" + _grp(l2)
        if dn == 2:
            return "n2_sp_" + ("sp" if l2 <= 1 else "df")
        return "core"
    if l == 2:
        if dn <= 1 and l2 <= 1:
            return "d_near"
        if dn == 1:
            return "n1_d_" + LL[l2]
        return "core"
    # f target
    if dn <= 1 and l2 <= 2:
        return "f_near"
    if dn == 1:
        return "n1_f_f"
    return "core"


CLASSES = ["same_s", "same_p", "same_d", "same_f", "sn_in_p", "n1_sp_sp", "n1_sp_d", "n1_sp_f",
           "n2_sp_sp", "n2_sp_df", "d_near", "n1_d_d", "n1_d_f", "f_near", "n1_f_f", "core",
           "sn_out", "out_sp", "out_d", "out_f"]
CIDX = {c: i for i, c in enumerate(CLASSES)}
OUTER = np.array([c.startswith("out_") or c == "sn_out" for c in CLASSES])

# Groups for the second-order (q-dependent) screening S1
S1_GROUPS = {
    "same": [c for c in CLASSES if c.startswith("same_")],
    "near": ["sn_in_p", "n1_sp_sp", "n1_sp_d", "n1_sp_f", "d_near", "n1_d_d", "n1_d_f",
             "f_near", "n1_f_f"],
    "core": ["n2_sp_sp", "n2_sp_df", "core"],
}
S1_MAT = np.zeros((len(CLASSES), len(S1_GROUPS)))
for _g, (_name, _cl) in enumerate(S1_GROUPS.items()):
    for _c in _cl:
        S1_MAT[CIDX[_c], _g] = 1.0

# relativistic classes of the removed electron
REL_CLASSES = ["s", "p1", "p3", "d", "f"]


def rel_class(l, k):
    if l == 0:
        return 0
    if l == 1:
        return 1 if k <= 2 else 2
    return 3 if l == 2 else 4

# ----------------------------------------------------------------------------- orbital rows


def orbital_rows(Z, shells, targets=None):
    """For each target subshell (n,l) of `shells` return a feature row."""
    out = []
    for n, l, k in shells:
        if targets is not None and (n, l) not in targets:
            continue
        c = np.zeros(len(CLASSES))
        for n2, l2, q2 in shells:
            cnt = q2 - 1 if (n2, l2) == (n, l) else q2
            if cnt > 0:
                c[CIDX[classify(n, l, n2, l2)]] += cnt
        out.append(dict(n=n, l=l, k=k, counts=c, Z=Z, j=jj_j(l, k), rc=rel_class(l, k),
                        kink=exchange_kink(l, k)))
    return out


def record_rows(Z, N, shells=None, form="1e"):
    """Rows (with weights) whose weighted sum of binding energies gives IE.

    form='1e'    : one row, weight +1 (removed electron, initial configuration)
    form='total' : IE = sum_{N config} q b - sum_{N-1 config} q b, identical rows cancelled
    """
    sN, sF, rem, rearr = initial_final(Z, N, shells)
    if form == "1e":
        rows = orbital_rows(Z, sN, targets={rem})
        rows[0]["w"] = 1.0
        return rows, rem, rearr
    rN = orbital_rows(Z, sN)
    rF = orbital_rows(Z, sF)
    for r in rN:
        r["w"] = float(r["k"])
    for r in rF:
        r["w"] = -float(r["k"])
    # cancel identical (target, occupancy, counts) rows
    keyed = {}
    for r in rN + rF:
        key = (r["n"], r["l"], r["k"], tuple(r["counts"]))
        if key in keyed:
            keyed[key]["w"] += r["w"]
        else:
            keyed[key] = dict(r)
    rows = [r for r in keyed.values() if abs(r["w"]) > 1e-12]
    return rows, rem, rearr


def assemble(entries):
    """entries: list of (rows, rem, rearr, Z, N). Returns dict of numpy arrays."""
    idx, n, l, k, w, Zr, j, rc, kink, C = [], [], [], [], [], [], [], [], [], []
    rec_n, rec_l, rec_Z, rec_N, rec_rearr, rec_k, rec_kink = [], [], [], [], [], [], []
    rec_rc, rec_C = [], []
    for i, (rows, rem, rearr, Z, N) in enumerate(entries):
        rr = [r for r in rows if (r["n"], r["l"]) == tuple(rem)]
        if not rr or rr[0]["w"] <= 0:   # total form: removed subshell row of the N config
            sN, _, _, _ = initial_final(Z, N)
            rr = orbital_rows(Z, sN, targets={tuple(rem)}) or rows[:1]
        rec_k.append(rr[0]["k"]); rec_kink.append(rr[0]["kink"]); rec_rc.append(rr[0]["rc"])
        rec_C.append(rr[0]["counts"])
        for r in rows:
            idx.append(i)
            n.append(r["n"]); l.append(r["l"]); k.append(r["k"]); w.append(r["w"])
            Zr.append(Z); j.append(r["j"]); rc.append(r["rc"]); kink.append(r["kink"])
            C.append(r["counts"])
        rec_n.append(rem[0]); rec_l.append(rem[1]); rec_Z.append(Z); rec_N.append(N)
        rec_rearr.append(rearr)
    a = dict(idx=np.array(idx), n=np.array(n, float), l=np.array(l), k=np.array(k, float),
             w=np.array(w), Z=np.array(Zr, float), j=np.array(j), rc=np.array(rc),
             kink=np.array(kink), C=np.array(C), nrec=len(entries),
             rec_n=np.array(rec_n, float), rec_l=np.array(rec_l), rec_Z=np.array(rec_Z, float),
             rec_N=np.array(rec_N, float), rec_rearr=np.array(rec_rearr),
             rec_k=np.array(rec_k, float), rec_kink=np.array(rec_kink), rec_rc=np.array(rec_rc),
             rec_C=np.array(rec_C))
    a["Za"] = a["Z"] - a["C"][:, ~OUTER].sum(1)          # local asymptotic charge seen
    a["Cin"] = a["C"][:, ~OUTER].sum(1)
    a["Zq"] = a["Z"] - a["C"].sum(1)                      # true asymptotic charge q+1
    return a


@lru_cache(maxsize=4)
def dataset(form="1e"):
    recs = load_records()
    ents = []
    for r in recs:
        rows, rem, rearr = record_rows(r["Z"], r["N"], None, form)
        ents.append((rows, rem, rearr, r["Z"], r["N"]))
    a = assemble(ents)
    a["y"] = np.array([r["IE_eV"] for r in recs])
    a["_feats"] = corr_features(a)
    return a

# ----------------------------------------------------------------------------- physics pieces


def dirac_factor(n, j, Zeff):
    """E_Dirac(n, j; Zeff) / E_Schrodinger(n; Zeff) for a point nucleus (>= 1)."""
    x = np.clip(Zeff * ALPHA, 1e-8, None)
    kap = j + 0.5
    x = np.minimum(x, 0.999 * kap)
    gam = np.sqrt(kap ** 2 - x ** 2)
    nr = n - kap + gam
    E = 1.0 - 1.0 / np.sqrt(1.0 + (x / nr) ** 2)      # binding / mc^2
    return E * 2.0 * n ** 2 / x ** 2

# ----------------------------------------------------------------------------- parameter spec


def default_spec():
    """Which blocks are active. Each block lists its parameter names and initial values."""
    return {
        "sigma": [c for c in CLASSES if c not in ("sn_out", "out_sp")],  # screening constants
        "S1": True,                      # q-dependent screening tau_g / (Z - S0 + kappa)
        "rel": True,                     # Dirac(Z_eff) x penetration correction r_{l,j}
        "exch": True,                    # Hund exchange kink x_l
        "qdef": False,                   # quantum defect n* = n - d_l/(1+(Zeff0-1)/kd)
        "corr": [],                      # sparse multiplicative correction terms (names)
        "form": "1e",
    }


SIGMA_INIT = {"same_s": 0.3, "same_p": 0.35, "same_d": 0.35, "same_f": 0.35, "sn_in_p": 0.6,
              "n1_sp_sp": 0.85, "n1_sp_d": 0.85, "n1_sp_f": 0.85, "n2_sp_sp": 1.0,
              "n2_sp_df": 1.0, "d_near": 0.95, "n1_d_d": 1.0, "n1_d_f": 1.0, "f_near": 0.95,
              "n1_f_f": 1.0, "core": 1.0, "sn_out": 0.2, "out_sp": 0.0, "out_d": 0.1,
              "out_f": 0.3}


def param_names(spec):
    names = ["sig_" + c for c in spec["sigma"]]
    if spec["S1"]:
        if spec.get("S1_form", "tau") == "pen":
            names += ["kap_" + g for g in spec.get("S1_names", list(S1_GROUPS))]
        else:
            names += ["tau_" + g for g in spec.get("S1_names", list(S1_GROUPS))] + ["kappa"]
        if spec.get("kappa_l"):
            names += ["kappa_p", "kappa_d", "kappa_f"]
        if spec.get("S2"):
            names += ["tau2_" + g for g in spec.get("S1_names", list(S1_GROUPS))]
    if spec["rel"] and spec.get("rel_pen", True):
        names += ["rel_" + c for c in REL_CLASSES]
    if spec["exch"]:
        names += ["x_p", "x_d", "x_f"]
    if spec["qdef"]:
        names += ["qd_s", "qd_p", "qd_d", "qd_f", "qd_k"]
    if spec.get("relS"):
        names += ["rs_" + c for c in spec["relS"]]
    names += ["b_" + t for t in spec["corr"]]
    tie = spec.get("tie", {})
    return [nm for nm in names if nm not in tie]


# Analogue used when a parameter has no support in a training set (e.g. no f electrons for
# Z <= 54): the f-electron quantity is set equal to the corresponding d-electron quantity.
ANALOG = {"sig_same_f": "sig_same_d", "sig_n1_sp_f": "sig_n1_sp_d", "sig_n1_d_f": "sig_n1_d_d",
          "sig_f_near": "sig_d_near", "sig_n1_f_f": "sig_n1_d_d", "sig_out_f": "sig_out_d",
          "rel_f": "rel_d", "x_f": "x_d"}


def init_params(spec):
    v = []
    for nm in param_names(spec):
        if nm.startswith("sig_"):
            c = nm[4:]
            v.append(SIGMA_INIT[c])
        elif nm == "kappa":
            v.append(7.0)
        elif nm.startswith("kap_"):
            v.append(1.0)
        elif nm.startswith("tau_"):
            v.append(0.0)
        elif nm == "qd_k":
            v.append(1.0)
        else:
            v.append(0.0)
    return np.array(v, float)

# ----------------------------------------------------------------------------- correction library


def corr_features(a):
    """Library of candidate record-level terms phi_m for the sparse multiplicative correction
    IE -> IE * exp(sum_m beta_m phi_m).  Only configuration quantities are used.
    Requires the record-level arrays of the 1e dataset (rows == records)."""
    Z, N, n, l = a["rec_Z"], a["rec_N"], a["rec_n"], a["rec_l"]
    k = a["rec_k"]
    C = a["rec_C"]
    iz = 1.0 / (Z - N + 1.0)                        # 1 / Z_a  (Z_a = q + 1)
    x2 = (Z * ALPHA) ** 2
    rc = a["rec_rc"]
    L = [(l == i).astype(float) for i in range(4)]
    half = (k == 2 * l + 1) & (l > 0)
    full = k == 4 * l + 2
    single = k == 1
    has_open_d_inner = np.zeros(len(Z), bool)
    f = {}
    for i, c in enumerate("spdf"):
        f["invZa_" + c] = L[i] * iz
        f["invZa2_" + c] = L[i] * iz ** 2
    for i, c in enumerate(REL_CLASSES):
        f["x2_" + c] = x2 * (rc == i)
    f["x4"] = x2 ** 2
    f["x4_s"] = x2 ** 2 * L[0]
    for i, c in enumerate("pdf"):
        f["kink_" + c + "_invZa"] = a["rec_kink"] * L[i + 1] * iz
    f["half_invZa"] = half * iz
    f["full_invZa"] = full * iz
    f["single_invZa"] = single * iz
    f["rearr"] = a["rec_rearr"].astype(float)
    f["rearr_invZa"] = a["rec_rearr"] * iz
    f["d_with_outer_s_invZa"] = (L[2] * (C[:, CIDX["out_d"]] > 0)) * iz
    f["d10core_invZa"] = C[:, CIDX["n1_sp_d"]] / 10.0 * iz
    f["f14core_sp_invZa"] = C[:, CIDX["n1_sp_f"]] / 14.0 * iz
    f["f14core_d_invZa"] = C[:, CIDX["n1_d_f"]] / 14.0 * iz
    f["out_f_invZa"] = C[:, CIDX["out_f"]] / 8.0 * iz
    f["NoverZ"] = N / Z
    f["lnZ"] = np.log(Z) / 5.0
    f["qed_s"] = (ALPHA / np.pi) * x2 * np.log(1.0 / np.maximum(x2, 1e-6)) * L[0] / n
    f["n_ge5_invZa"] = (n >= 5) * iz
    f["sameocc_invZa"] = (k - 1) / (4 * l + 2) * iz
    del has_open_d_inner
    many = (N > 1).astype(float)          # one-electron ions are exact (Dirac): no corrections
    return {k_: v * many for k_, v in f.items()}


CORR_LIBRARY = None  # filled lazily: list(corr_features(dataset("1e")))


# ----------------------------------------------------------------------------- model evaluation


def evaluate_model(theta, spec, a, return_parts=False):
    """Vectorized IE prediction for an assembled dataset `a` (from dataset()/assemble())."""
    P = dict(zip(param_names(spec), theta))
    for k_, v_ in spec.get("tie", {}).items():
        P[k_] = P[v_]
    ncls = len(CLASSES)
    sig = np.zeros(ncls)
    for c in spec["sigma"]:
        sig[CIDX[c]] = P["sig_" + c]
    # classes not in the active list keep their fixed defaults (merged/frozen by selection)
    for c in CLASSES:
        if c not in spec["sigma"]:
            sig[CIDX[c]] = spec.get("fixed_sigma", {}).get(c, SIGMA_INIT[c])
    for c, rep in spec.get("sigma_tie", {}).items():   # coarse (Slater-like) class merging
        sig[CIDX[c]] = P["sig_" + rep]
    C = a["C"]
    ZA = "Zq" if spec.get("Za_mode", "local") == "true" else "Za"
    S0 = C @ sig
    Ze0 = a["Z"] - S0
    Ze = Ze0
    if spec["S1"] and spec.get("S1_form", "tau") == "pen":
        # penetration form: the unscreened charge of each group, P_g = sum_c (1-sigma_c) c_c,
        # is only fully felt by a strongly bound electron; near neutrality it is reduced by
        # Z_a / (Z_a + kappa_g):   Z_eff = Z_a + sum_g P_g Z_a/(Z_a + kappa_g)
        S1M = spec.get("S1_mat", S1_MAT)
        groups = spec.get("S1_names", list(S1_GROUPS))
        Pg = (C * (1.0 - sig)) @ S1M                       # records x groups
        kg = np.abs(np.array([P["kap_" + g] for g in groups])) + 0.02
        if spec.get("kappa_l"):
            kg = kg[None, :] * np.exp(np.array([0.0, P["kappa_p"], P["kappa_d"], P["kappa_f"]])[a["l"]])[:, None]
        Za_ = np.maximum(a[ZA], 1.0)
        Ze = Ze0 - (Pg * (kg / (Za_[:, None] + kg))).sum(1)
    elif spec["S1"]:
        S1M = spec.get("S1_mat", S1_MAT)
        groups = spec.get("S1_names", list(S1_GROUPS))
        tau = np.array([P["tau_" + g] for g in groups])
        S1 = (C @ S1M) @ tau
        kap = abs(P["kappa"]) + 0.05
        if spec.get("kappa_l"):
            kap = kap + np.array([0.0, P["kappa_p"], P["kappa_d"], P["kappa_f"]])[a["l"]]
        den = np.maximum(Ze0, 0.05) if spec.get("S1_den", "Za") == "Ze0" else np.maximum(a[ZA], 1.0)
        Ze = Ze0 - S1 / (den + kap)
        if spec.get("S2"):
            tau2 = np.array([P["tau2_" + g] for g in groups])
            Ze = Ze - ((C @ S1M) @ tau2) / (den + kap) ** 2
    if spec.get("relS"):
        # indirect relativistic effect: relativistic contraction of the inner s,p shells
        # increases their screening of the target: dS = rho_l (Z alpha)^2 * (# inner electrons)
        rho = np.zeros(4)
        for c in spec["relS"]:
            rho["spdf".index(c)] = P["rs_" + c]
        Ze = Ze - rho[a["l"]] * (a["Z"] * ALPHA) ** 2 * a["Cin"]
    Ze_pos = np.maximum(Ze, 1e-3)
    n = a["n"]
    neff = n
    if spec["qdef"]:
        d = np.array([P["qd_s"], P["qd_p"], P["qd_d"], P["qd_f"]])[a["l"]]
        neff = n - d / (1.0 + np.maximum(Ze0 - 1.0, 0) / (abs(P["qd_k"]) + 0.05))
    b = RY * (Ze_pos / neff) ** 2
    if spec["rel"]:
        r = (np.array([P["rel_" + c] for c in REL_CLASSES])[a["rc"]]
             if spec.get("rel_pen", True) else 0.0)
        Za = np.maximum(a[ZA], 1.0)
        F = dirac_factor(n, a["j"], Ze_pos) * (1.0 + (a["Z"] * ALPHA) ** 2 * r * (Ze_pos / Za - 1.0) / n)
        b = b * F
    if spec["exch"]:
        xl = np.array([0.0, P["x_p"], P["x_d"], P["x_f"]])[a["l"]]
        xs = {"Ze": Ze_pos, "Za": np.maximum(a[ZA], 1.0), "one": 1.0}[spec.get("exch_scale", "Ze")]
        b = b + RY * xl * a["kink"] * xs / n ** 2
    ie = np.bincount(a["idx"], weights=a["w"] * b, minlength=a["nrec"])
    if spec["corr"]:
        feats = a.get("_feats")
        if feats is None:
            feats = corr_features(a)
        s = np.zeros(a["nrec"])
        for t in spec["corr"]:
            s += P["b_" + t] * feats[t]
        ie = ie * np.exp(s)
    if return_parts:
        return ie, dict(Ze=Ze_pos, S0=S0)
    return ie

# ----------------------------------------------------------------------------- public predict


@lru_cache(maxsize=1)
def _load_params(path=PARAMS_JSON):
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    return d["final"]["spec"], np.array(d["final"]["values"], float)


def predict(Z, N, shells=None, params=None):
    """Ionization energy (eV) of the N-electron ion of element Z.

    shells : optional configuration [(n,l,occ),...] of the N-electron ion; default is the NIST
             ground configuration (Madelung outside the table). With shells=None the final
             (N-1)-electron configuration is the ground configuration too, so rearrangements
             (e.g. Ni 3d8 4s2 -> Ni+ 3d9) are handled by the removed-subshell rule.
    params : optional (spec, theta) tuple; default = final parameters in results/se_params.json
    """
    spec, theta = params if params is not None else _load_params()
    rows, rem, rearr = record_rows(int(Z), int(N), shells, spec["form"])
    a = assemble([(rows, rem, rearr, int(Z), int(N))])
    return float(evaluate_model(theta, spec, a)[0])


def predict_many(ZN, shells_list=None, params=None):
    spec, theta = params if params is not None else _load_params()
    ents = []
    for i, (Z, N) in enumerate(ZN):
        sh = None if shells_list is None else shells_list[i]
        rows, rem, rearr = record_rows(int(Z), int(N), sh, spec["form"])
        ents.append((rows, rem, rearr, int(Z), int(N)))
    return evaluate_model(theta, spec, assemble(ents))


if __name__ == "__main__":
    a = dataset("1e")
    tot = a["C"].sum(0)
    nrec = (a["C"] > 0).sum(0)
    for c, t, m in zip(CLASSES, tot, nrec):
        print(f"{c:12s} electrons={t:9.0f} records={m:5d}")
