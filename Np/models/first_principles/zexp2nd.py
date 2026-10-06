"""Second-order (dE2) completions of the exact two-term 1/Z expansion.

Base (zero fitted parameters, all exact or first-principles):
    B(Z,N) = mu(Z) [Z^2 dE0 + Z dE1] + dR(Z)       (hartree; dR = Dirac+QED+FNS of removed e-, p=2)
Completions:
  sq    : dE2 = dE1^2 / (4 dE0)                     (closed form: IE = dE0 (Z - sigma1)^2)
  lit   : exact literature E2 for N <= 3:  E2(1s^2) = -0.157666429 (Scherr-Knight 1963; Baker et al. 1990),
          E2(1s^2 2s) ~ -0.408165 (Yan, Tambasco & Drake 1998; quoted to ~1e-5, not re-verified here),
          so dE2(He-like) = 0.157666, dE2(Li-like 2s) = 0.250499; others fall back to sq.
  seq   : one fitted dE2 per (N-electron config -> (N-1)-electron config) key  [semi-empirical]
  univ  : universal screening function sigma = sigma1 + sum_{k=1..3} c_{l,k} (N/Z)^k,
          IE = mu dE0 (Z - sigma)^2 + dR; 4 x 3 = 12 fitted parameters  [semi-empirical]
Splits: S1 (Z%5==0 test), S2 (train Z<=54), S3 (N%6==0 test).
"""
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zexp  # noqa: E402
from zexp import ROOT  # noqa: E402
sys.path.insert(0, ROOT)
from common.atomdata import load_records, HARTREE_EV, shells_to_str, ground_shells  # noqa: E402

E2_LIT = {"1s2": -0.157666429, "1s2.2s1": -0.408165}
SPLITS = {"S1": lambda r: r["Z"] % 5 == 0, "S2": lambda r: r["Z"] >= 55, "S3": lambda r: r["N"] % 6 == 0}


def base_terms(Z, N, shells=None, p=2.0):
    c = zexp.coeffs(Z, N, shells)
    mu = zexp.reduced_mass_factor(Z)
    B = mu * (Z * Z * c["dE0"] + Z * c["dE1"]) + zexp.rel_correction(Z, c, p=p)
    return c, mu, B


def key_of(Z, N, shells=None):
    sh = shells if shells is not None else ground_shells(Z, N)
    shi = ground_shells(Z, N - 1) if N > 1 else []
    return shells_to_str(sh) + ">" + (shells_to_str(shi) if shi else "-")


def dE2_lit(Z, N, shells=None):
    sh = shells if shells is not None else ground_shells(Z, N)
    shi = ground_shells(Z, N - 1) if N > 1 else []
    a, b = shells_to_str(sh), (shells_to_str(shi) if shi else "-")
    if a in E2_LIT and (b in E2_LIT or b in ("1s1", "-")):
        return (E2_LIT.get(b, 0.0)) - E2_LIT[a]
    return None


def table():
    """Arrays for all NIST rows."""
    R = load_records()
    out = {k: [] for k in ("Z", "N", "IE", "B", "dE0", "dE1", "sig1", "l", "mu", "key", "lit")}
    for r in R:
        c, mu, B = base_terms(r["Z"], r["N"])
        out["Z"].append(r["Z"]); out["N"].append(r["N"]); out["IE"].append(r["IE_eV"] / HARTREE_EV)
        out["B"].append(B); out["dE0"].append(c["dE0"]); out["dE1"].append(c["dE1"])
        out["sig1"].append(c["sigma"]); out["l"].append(c["l"]); out["mu"].append(mu)
        out["key"].append(key_of(r["Z"], r["N"]))
        out["lit"].append(dE2_lit(r["Z"], r["N"]))
    T = {k: np.array(v) if k not in ("key", "lit") else v for k, v in out.items()}
    T["recs"] = R
    return T


def sq_dE2(T):
    d0, d1 = T["dE0"], T["dE1"]
    return np.where(d0 > 0, d1 ** 2 / (4 * np.where(d0 > 0, d0, 1)), 0.0)


def pred_sq(T):
    return T["B"] + T["mu"] * sq_dE2(T)


def pred_lit(T):
    p = pred_sq(T).copy()
    for i, v in enumerate(T["lit"]):
        if v is not None:
            p[i] = T["B"][i] + T["mu"][i] * v
    return p


def fit_seq(T, train):
    """Per-key dE2 minimising sum of squared relative errors (closed form per key)."""
    params = {}
    y = T["IE"] - T["B"]
    w = 1.0 / T["IE"] ** 2
    keys = np.array(T["key"])
    for k in set(keys[train]):
        m = train & (keys == k)
        params[k] = float(np.sum(w[m] * y[m] * T["mu"][m]) / np.sum(w[m] * T["mu"][m] ** 2))
    return params


def pred_seq(T, params):
    p = pred_sq(T).copy()
    for i, k in enumerate(T["key"]):
        if k in params:
            p[i] = T["B"][i] + T["mu"][i] * params[k]
    return p


def univ_features(T):
    x = T["N"] / T["Z"]
    F = np.zeros((len(x), 12))
    for li in range(4):
        m = np.minimum(T["l"], 3) == li
        for k in range(3):
            F[m, 3 * li + k] = x[m] ** (k + 1)
    return F


def pred_univ(T, theta):
    sig = T["sig1"] + univ_features(T) @ theta
    zo = np.maximum(T["Z"] - sig, 0.05)
    nr = T["mu"] * T["dE0"] * zo ** 2
    rel = T["B"] - T["mu"] * (T["Z"] ** 2 * T["dE0"] + T["Z"] * T["dE1"])
    p = nr + rel
    bad = ~(T["dE0"] > 0) | ~np.isfinite(T["sig1"])
    p[bad] = pred_sq(T)[bad]
    return p


def fit_univ(T, train):
    from scipy.optimize import least_squares
    ok = train & (T["dE0"] > 0) & np.isfinite(T["sig1"])

    def res(th):
        return np.log(np.maximum(pred_univ(T, th)[ok], 1e-6) / T["IE"][ok])
    r = least_squares(res, np.zeros(12), loss="linear", max_nfev=4000)
    return r.x


def metrics(T, pred, mask):
    ape = np.abs(pred[mask] / T["IE"][mask] - 1) * 100
    neu = mask & (T["N"] == T["Z"])
    apn = np.abs(pred[neu] / T["IE"][neu] - 1) * 100
    return {"n": int(mask.sum()), "MAPE_%": float(ape.mean()), "median_APE_%": float(np.median(ape)),
            "neutral_n": int(neu.sum()), "neutral_MAPE_%": float(apn.mean()) if neu.any() else None}


def write_preds(name, T, pred):
    path = os.path.join(ROOT, "results", f"fp_{name}_predictions.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Z", "N", "IE_pred_eV"])
        for Z, N, v in zip(T["Z"], T["N"], pred):
            w.writerow([int(Z), int(N), repr(float(v * HARTREE_EV))])
    return path


def main():
    import evaluate
    T = table()
    R = T["recs"]
    allm = np.ones(len(R), bool)
    val = {}
    # zero-parameter variants
    for name, pr in (("zexp3_lit", pred_lit(T)),):
        path = write_preds(name, T, pr)
        json.dump(evaluate.evaluate(path), open(path.replace("_predictions.csv", "_metrics.json"), "w"), indent=1)
    # held-out validation for fitted variants
    for sp, fn in SPLITS.items():
        test = np.array([fn(r) for r in R])
        train = ~test
        ps = fit_seq(T, train)
        pseq = pred_seq(T, ps)
        th = fit_univ(T, train)
        pun = pred_univ(T, th)
        cover = np.array([k in ps for k in T["key"]])
        val[sp] = {"seq": dict(metrics(T, pseq, test), n_params=len(ps),
                               test_rows_with_fitted_key=int((test & cover).sum())),
                   "univ": dict(metrics(T, pun, test), n_params=12),
                   "sq_zero_param": metrics(T, pred_sq(T), test)}
        print(sp, json.dumps(val[sp]), flush=True)
    # final fits on all data
    ps = fit_seq(T, allm)
    th = fit_univ(T, allm)
    for name, pr in (("zexp3_seq", pred_seq(T, ps)), ("zexp_univ", pred_univ(T, th))):
        path = write_preds(name, T, pr)
        json.dump(evaluate.evaluate(path), open(path.replace("_predictions.csv", "_metrics.json"), "w"), indent=1)
    json.dump({"validation": val, "univ_theta": th.tolist(),
               "univ_theta_names": [f"c_{'spdf'[l]}{k + 1}" for l in range(4) for k in range(3)],
               "seq_n_params_all": len(ps), "seq_dE2": ps},
              open(os.path.join(ROOT, "results", "fp_zexp2nd_validation.json"), "w"), indent=1)
    zexp.save_cfg_cache()


if __name__ == "__main__":
    main()
