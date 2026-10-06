r"""Pocket formula: 8 fitted parameters, hand-calculable (Slater-rule-like counting).

    IE = Ry (Zeff/n)^2 [1 + (Zeff alpha / n)^2 (n/(j+1/2) - 3/4)] + Ry x K_l(k) / n^2
    Zeff = Z - s_same nu_same - s_in nu_in - s_core nu_core - s_df nu_df - s_out nu_out
             - t (N - 1) / (Z - N + 1 + kappa)

nu_same : other electrons in the same subshell (n,l)
nu_in   : s,p target: same-n s electrons (for a p electron) and all (n-1)-shell electrons
nu_core : s,p target: all electrons with n' <= n-2
nu_df   : d,f target: all electrons with n' <= n (other subshells)
nu_out  : electrons with n' > n (e.g. 4s on 3d)
K_l(k)  : Hund kink of the removed electron (p: 0, .6, 1.2, -1.2, -.6, 0 for k = 1..6; d, f
          analogous, = [P(k)-P(k-1)] - 2l(k-1)/(4l+1), P = number of parallel-spin pairs)
The relativistic bracket is the leading (Sommerfeld) Dirac correction for charge Zeff (no
parameters); j = l-1/2 while the subshell holds <= 2l electrons, else l+1/2.
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import umodel as U  # noqa: E402
from umodel import RY, ALPHA, CIDX, CLASSES  # noqa: E402

NAMES = ["s_same", "s_in", "s_core", "s_df", "s_out", "t", "kappa", "x"]


def counts(a):
    """Pocket-group electron counts from the Track B class counts a['C']."""
    C = a["C"]
    sp = a["l"] <= 1

    def g(cl):
        return C[:, [CIDX[c] for c in cl]].sum(1)
    same = g([c for c in CLASSES if c.startswith("same_")])
    inn = g(["sn_in_p", "n1_sp_sp", "n1_sp_d", "n1_sp_f"])
    core_sp = g(["n2_sp_sp", "n2_sp_df"]) + np.where(sp, C[:, CIDX["core"]], 0)
    df = g(["d_near", "n1_d_d", "n1_d_f", "f_near", "n1_f_f"]) + np.where(sp, 0, C[:, CIDX["core"]])
    out = g(["out_sp", "out_d", "out_f", "sn_out"])
    return np.stack([same, inn, core_sp, df, out], 1)


def evaluate(theta, a):
    s, t, kap, x = theta[:5], theta[5], abs(theta[6]) + 0.05, theta[7]
    Za = np.maximum(a["Zq"], 1.0)
    N = a["Z"] - a["Zq"] + 1
    Ze = a["Z"] - counts(a) @ s - t * (N - 1) / (Za + kap)
    Ze = np.maximum(Ze, 1e-3)
    n = a["n"]
    rel = 1.0 + (Ze * ALPHA / n) ** 2 * (n / (a["j"] + 0.5) - 0.75)
    return RY * (Ze / n) ** 2 * rel + RY * x * a["kink"] / n ** 2


def init():
    return np.array([0.35, 0.85, 1.0, 1.0, 0.0, 0.0, 5.0, 0.0])


def fit(mask=None):
    from scipy.optimize import least_squares
    a = U.dataset()
    if mask is not None:
        a = U.subset(a, mask)
    y = a["y"]
    r = least_squares(lambda th: np.log(np.maximum(evaluate(th, a), 1e-6 * y) / y), init(),
                      x_scale="jac", max_nfev=800)
    return r.x


def predict(Z, N, theta, shells=None):
    a = U.build([(int(Z), int(N))], None if shells is None else [shells])
    return float(evaluate(np.asarray(theta), a)[0])
