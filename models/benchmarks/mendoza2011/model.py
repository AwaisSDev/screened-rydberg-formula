"""Relativistic screened hydrogenic model (NRSHM) of Mendoza et al., HEDP 7 (2011) 169.

Published constants only (constants.json, Tables 1-2 of the paper); nothing is refit.

Equations (paper numbering, atomic units, c = 1/alpha):
  (1)  E_T = sum_k P_k eps_k
  (2)  eps_k = c^2 * { [1 + (alpha Q_k / (n - j - 1/2 + sqrt((j+1/2)^2 - (alpha Q_k)^2)))^2]^(-1/2) - 1 }
  (3)  Q_k = Z - sum_k' sigma_kk' (P_k' - delta_kk')
Subshells k = (n, l, j) from 1s1/2 to 5p3/2 only (the paper publishes nothing beyond 5p3/2).

Ionization energy prescriptions
  'delta'  (primary): IE = E_T(N-1) - E_T(N) with NIST ground configurations
           (common.atomdata.ground_shells). The paper fits and tabulates ionization
           energies with the model's total energies (Table 3 is reproduced this way, see
           validate_paper.py) and recommends total-energy differences for isolated atoms
           (Sec. 5.2).
  'xalpha': one-electron (Koopmans-like) binding energy  -dE_T/dP_k  (App. A, eq. A.1-A.2)
           of the subshell that loses the electron.
  'half'  : Slater transition state  -dE_T/dP_k evaluated at P_k - 1/2 ("derivation in 1/2").

j-splitting of a non-relativistic (n,l,q) subshell (paper gives no rule for this):
  'low'  (primary): fill j = l-1/2 first (capacity 2l), then j = l+1/2 -- the jj ground
          configurations the paper prints in Table 4 (e.g. N I 2p1/2^2 2p3/2^1).
  'stat': statistical split q*(2j+1)/(2(2l+1)) (fractional occupations).
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from common.atomdata import ALPHA, HARTREE_EV, ground_shells  # noqa: E402

_C = json.load(open(os.path.join(HERE, "constants.json"), encoding="utf-8"))
LABELS = _C["subshells"]
L_OF = {"s": 0, "p": 1, "d": 2, "f": 3}


def _parse(lab):  # '3d5/2' -> (3, 2, 2.5)
    return int(lab[0]), L_OF[lab[1]], int(lab[2]) / 2.0


KEYS = [_parse(x) for x in LABELS]
IDX = {k: i for i, k in enumerate(KEYS)}
SIG = [[_C["sigma"][a][b] for b in LABELS] for a in LABELS]
CLIGHT = 1.0 / ALPHA


class NotCovered(ValueError):
    pass


def jsplit(shells, mode="low"):
    """[(n,l,q)] -> occupation vector over the 19 nlj subshells."""
    P = [0.0] * len(KEYS)
    for n, l, q in shells:
        if q == 0:
            continue
        if l == 0:
            parts = [(0.5, float(q))]
        elif mode == "low":
            lo = min(q, 2 * l)
            parts = [(l - 0.5, float(lo)), (l + 0.5, float(q - lo))]
        else:
            parts = [(l - 0.5, q * (2 * l) / (2 * (2 * l + 1))), (l + 0.5, q * (2 * l + 2) / (2 * (2 * l + 1)))]
        for j, p in parts:
            if p == 0:
                continue
            key = (n, l, j)
            if key not in IDX:
                raise NotCovered(f"subshell n={n} l={l} j={j} beyond 5p3/2 (no published constants)")
            P[IDX[key]] += p
    return P


def charges(Z, P):
    return [Z - sum(SIG[k][kp] * (P[kp] - (1.0 if kp == k else 0.0)) for kp in range(len(P)))
            for k in range(len(P))]


def eps(k, Q):
    n, l, j = KEYS[k]
    x = ALPHA * Q
    kap = j + 0.5
    if x >= kap:
        raise ValueError("alpha*Q >= j+1/2")
    d = n - kap + math.sqrt(kap * kap - x * x)
    return CLIGHT ** 2 * ((1.0 + (x / d) ** 2) ** -0.5 - 1.0)


def total_energy(Z, P):
    """E_T in hartree (eq. 1); only occupied subshells contribute."""
    Q = charges(Z, P)
    return sum(P[k] * eps(k, Q[k]) for k in range(len(P)) if P[k] != 0)


def binding_xalpha(Z, P, k, h=1e-5):
    """eps_k^Xalpha = dE_T/dP_k (eq. A.1), central difference in the occupation (hartree)."""
    Pp = list(P); Pm = list(P)
    Pp[k] += h; Pm[k] -= h
    return (total_energy(Z, Pp) - total_energy(Z, Pm)) / (2 * h)


def _removed_k(Pi, Pf):
    """Index of the nlj subshell losing the most occupation (outermost on ties)."""
    d = [(Pi[k] - Pf[k], k) for k in range(len(Pi))]
    return max(d)[1]


def predict(Z, N, shells=None, method="delta", split="low", shells_final=None):
    """Ionization energy (eV) of the N-electron ion of element Z."""
    si = shells if shells is not None else ground_shells(Z, N)
    sf = shells_final if shells_final is not None else (ground_shells(Z, N - 1) if N > 1 else [])
    Pi = jsplit(si, split)
    Pf = jsplit(sf, split)
    if method == "delta":
        return (total_energy(Z, Pf) - total_energy(Z, Pi)) * HARTREE_EV
    k = _removed_k(Pi, Pf)
    if method == "xalpha":
        return -binding_xalpha(Z, Pi, k) * HARTREE_EV
    if method == "half":
        P = list(Pi); P[k] -= 0.5
        return -binding_xalpha(Z, P, k) * HARTREE_EV
    raise ValueError(method)


if __name__ == "__main__":
    for Z, N in [(1, 1), (2, 2), (6, 6), (26, 26), (26, 6), (79, 11)]:
        print(Z, N, round(predict(Z, N), 4), "eV")
