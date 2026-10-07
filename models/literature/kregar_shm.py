"""Parameter-free screened hydrogenic model of Kregar / Di Rocco (literature head-to-head).

Published method reproduced here (no fitted constants; never reads NIST IE values):
  * M. Kregar, Phys. Scr. 29, 438 (1984); 31, 246 (1985)          -- splitting of 1/r_ij into one-body terms
  * H. O. Di Rocco, Braz. J. Phys. 22, 227 (1992)                  -- exchange / sub-shell corrections
  * J. Pomarico, D. I. Iriarte, H. O. Di Rocco, Braz. J. Phys. 35, 130 (2005), Eqs. (7), (12)-(15)
  * H. O. Di Rocco, F. Lanzini, Braz. J. Phys. 46, 175 (2016), Eqs. (5)-(9) (definitions of g, f; Dirac energy)

Model.  Each electron of sub-shell i (n_i, l_i, occupation q_i) moves in a hydrogenic orbital of charge
Z_i = Z - sigma_i.  The pair operator is split as 1/r_ij = g_ij/r_i + f_ji/r_j (Kregar), giving

    S_ij = n_i^2/Z_i * Int Int_{r_i > r_j} rho_i(r_i) rho_j(r_j) / r_i     (screening of i by ONE electron of j)

(i inner: S_ij = g_ij, the "external" screening; i outer: S_ij = f_ij, the "internal" screening), and the
exchange / sub-shell correction S_ij -> S_ij (1 - eps_ij), eps_ij = (F0 - {ij}) / F0 with the
average-of-configuration pair energies (Cowan):
    {ij} = F0 - 1/2 sum_k (l k l'; 000)^2 G^k            (i != j)
    {ii} = F0 - (2l+1)/(4l+1) sum_{k>0} (l k l; 000)^2 F^k
Then (Pomarico et al. 2005, Eq. 13)
    sigma_i = sum_{j != i} q_j S_ij (1 - eps_ij) + (q_i - 1) S_ii (1 - eps_ii),
iterated to self-consistency (S depends on all Z_j).  Configuration energy (Eq. 7):
    E = - sum_i q_i Z_i^2 / (2 n_i^2)   [hartree]
which equals <H> in the product of the screened hydrogenic orbitals by construction of the split.

Relativistic variants (the screening itself is always non-relativistic):
  'pauli': + Eq. (15) of Pomarico et al. evaluated with the screened Z_i.  The printed formula shows Z^4;
           we use Z_i^4 (with the bare Z, valence electrons of heavy neutral atoms would get eV-size shifts).
           This is our interpretation and is disclosed.
  'dirac': the orbital term replaced by the (2j+1)-weighted nl-average of the point-nucleus Dirac eigenvalue
           with charge Z_i (Di Rocco & Lanzini 2016 Eq. 9 works with nlj sub-shells; NIST ground configurations
           are given per nl, so we average over j -- a disclosed simplification).
IE(Z, N) = E(Z, N-1) - E(Z, N), both in the NIST ground configuration (configuration-average energies).

Reproduction caveat (see validate_kregar.py): the diagonal parameters k_ii and the total energies of
Pomarico et al. Table 1 are reproduced (<=0.2 %), but their printed Z->inf off-diagonal g/f values (Table 6)
differ from ours by up to 0.05, because the paper evaluates g, f from fitted closed forms (Eq. 14, with
coefficients a_k, b_k given only in Di Rocco 1992, which we could not obtain) and apportions exchange in a way
the 2005 paper does not fully specify.  This implementation follows the published DEFINITIONS exactly.
"""
import math
import os
import sys
from functools import lru_cache

import numpy as np
from scipy.special import eval_genlaguerre, gammaln

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from common.atomdata import load_records, madelung_shells, HARTREE_EV, ALPHA  # noqa: E402

# ---------------------------------------------------------------- radial grid (uniform in x = ln r)
NGRID = int(os.environ.get("KREGAR_NGRID", "700"))
X = np.linspace(math.log(1e-6), math.log(800.0), NGRID)
R = np.exp(X)
W = R * (X[1] - X[0])          # dr = r dx


def P_nl(n, l, Z):
    """Normalised hydrogenic radial function P = r R_nl(r) for charge Z on the grid."""
    rho = 2.0 * Z * R / n
    lognorm = 0.5 * (3 * math.log(2.0 * Z / n) + gammaln(n - l) - math.log(2 * n) - gammaln(n + l + 1))
    with np.errstate(under="ignore", over="ignore", invalid="ignore"):
        P = R * np.exp(lognorm - rho / 2) * rho ** l * eval_genlaguerre(n - l - 1, 2 * l + 1, rho)
    P = np.nan_to_num(P)
    P /= math.sqrt(np.sum(P * P * W))
    return P


@lru_cache(maxsize=None)
def w3j000_sq(l1, k, l2):
    """(l1 k l2; 0 0 0)^2 (zero unless l1+k+l2 even and triangular)."""
    J = l1 + k + l2
    if J % 2 or k < abs(l1 - l2) or k > l1 + l2:
        return 0.0
    g = J // 2
    lf = math.lgamma
    val = (lf(J - 2 * l1 + 1) + lf(J - 2 * k + 1) + lf(J - 2 * l2 + 1) - lf(J + 2)
           + 2 * (lf(g + 1) - lf(g - l1 + 1) - lf(g - k + 1) - lf(g - l2 + 1)))
    return math.exp(val)


def _rk(D, k):
    """Int Int D(r1) D(r2) r<^k / r>^(k+1) for a stack of pair densities D[..., G]."""
    a = D * W
    inner = np.cumsum(a * R ** k, axis=-1) / R ** (k + 1)
    outer = np.cumsum((a / R ** (k + 1))[..., ::-1], axis=-1)[..., ::-1] * R ** k
    Y = inner + outer - a / R            # diagonal point counted once
    return np.sum(a * Y, axis=-1)


def _orbitals(shells, Zeff):
    return np.array([P_nl(n, l, z) for (n, l, _), z in zip(shells, Zeff)])


def _monopole(shells, Zeff, P):
    rho = P * P
    Q = np.cumsum(rho * W, axis=1) - 0.5 * rho * W       # charge of j inside r (midpoint)
    A = (rho / R * W) @ Q.T                              # Int rho_i/r Q_j : region r_i > r_j
    n2z = np.array([s[0] ** 2 for s in shells], float) / np.asarray(Zeff)
    F0 = A + A.T
    S = n2z[:, None] * A
    np.fill_diagonal(S, n2z * np.diag(F0) / 2)
    return S, F0


def exchange_fraction(shells, Zeff):
    P = _orbitals(shells, Zeff)
    _, F0 = _monopole(shells, Zeff, P)
    m = len(shells)
    ls = [s[1] for s in shells]
    D = P[:, None, :] * P[None, :, :]
    Gk = {k: _rk(D, k) for k in range(2 * max(ls) + 1)}   # G^k (i != j) and F^k (i == j)
    eps = np.zeros((m, m))
    for i in range(m):
        for j in range(m):
            li, lj = ls[i], ls[j]
            if i == j:
                ex = (2 * li + 1) / (4 * li + 1) * sum(w3j000_sq(li, k, li) * Gk[k][i, i]
                                                        for k in range(2, 2 * li + 1, 2))
            else:
                ex = 0.5 * sum(w3j000_sq(li, k, lj) * Gk[k][i, j]
                               for k in range(abs(li - lj), li + lj + 1, 2))
            eps[i, j] = ex / F0[i, j]
    return eps


def solve_config(Z, shells, init=None, tol=1e-9, max_iter=400, rounds=3, mix=0.6, hist=6):
    """Self-consistent screened charges Z_i for a configuration [(n,l,q),...].
    init: optional {(n,l): Z_i} warm start (e.g. from the neighbouring ion); the converged fixed point
    does not depend on it.  Exchange fractions are refreshed at the converged charges (rounds=3).
    The fixed-point iteration Z <- Z - sigma(Z) is accelerated by Anderson mixing (numerics only)."""
    shells = [s for s in shells if s[2] > 0]
    if not shells:
        return [], np.zeros(0)
    q = np.array([s[2] for s in shells], float)
    Zeff = np.array([(init or {}).get((n, l), float(Z)) for n, l, _ in shells])

    def g(x, eps):
        S, _ = _monopole(shells, x, _orbitals(shells, x))
        M = S * (1 - eps)
        return np.maximum(Z - (M @ q - np.diag(M)), 0.02)

    for _ in range(rounds):
        eps = exchange_fraction(shells, Zeff)
        X, F = [], []
        x = Zeff.copy()
        for _it in range(max_iter):
            f = g(x, eps) - x
            if np.max(np.abs(f)) < tol * max(1.0, Z):
                break
            X.append(x.copy()); F.append(f.copy())
            X, F = X[-hist:], F[-hist:]
            if len(F) > 1:
                dF = np.array([F[i + 1] - F[i] for i in range(len(F) - 1)]).T
                dX = np.array([X[i + 1] - X[i] for i in range(len(X) - 1)]).T
                gam = np.linalg.lstsq(dF, f, rcond=None)[0]
                xn = x + mix * f - (dX + mix * dF) @ gam
            else:
                xn = x + mix * f
            x = np.maximum(xn, 0.02)
        Zeff = x
    return shells, Zeff


def _dirac_nl_avg(n, l, z):
    """(2j+1)-weighted average over j of the Dirac point-nucleus eigenvalue minus mc^2 [hartree]."""
    a = ALPHA
    out, wsum = 0.0, 0
    for j2 in ([1] if l == 0 else [2 * l - 1, 2 * l + 1]):
        kap = (j2 + 1) / 2
        az = min(a * z, kap - 1e-9)
        gam = math.sqrt(kap * kap - az * az)
        e = (1.0 / math.sqrt(1 + (az / (n - kap + gam)) ** 2) - 1.0) / a ** 2
        out += (j2 + 1) * e
        wsum += j2 + 1
    return out / wsum


def config_energies(Z, shells, init=None):
    """(E_nonrel, E_pauli, E_dirac) in hartree, and {(n,l): Z_i}."""
    shells, Zeff = solve_config(Z, list(shells), init=init)
    En = Ep = Ed = 0.0
    for (n, l, q), z in zip(shells, Zeff):
        en = -q * z * z / (2 * n * n)
        En += en
        Ep += en - q * ALPHA ** 2 / 2 * z ** 4 / n ** 3 * (1 / (l + 0.5) - 0.75 / n - (1 if l == 0 else 0))
        Ed += q * _dirac_nl_avg(n, l, z)
    return (En, Ep, Ed), {(n, l): z for (n, l, _), z in zip(shells, Zeff)}


# ---------------------------------------------------------------- configurations (no IE values used)
@lru_cache(maxsize=1)
def _config_table():
    return {(r["Z"], r["N"]): tuple(r["shells"]) for r in load_records()}


def config(Z, N):
    if N <= 0:
        return ()
    return _config_table().get((Z, N)) or tuple(madelung_shells(N))


def isonuclear_energies(Z):
    """{N: (E_nr, E_pauli, E_dirac)} for N = 0..Z, warm-starting each ion from the previous one."""
    out, init = {0: (0.0, 0.0, 0.0)}, None
    for N in range(1, Z + 1):
        e, zmap = config_energies(Z, config(Z, N), init=init)
        out[N] = e
        init = zmap
    return out


@lru_cache(maxsize=None)
def _iso(Z):
    return isonuclear_energies(Z)


def energies(Z, N):
    return _iso(Z)[N]


VARIANTS = {"nr": 0, "pauli": 1, "dirac": 2}


def predict(Z, N, variant="pauli"):
    k = VARIANTS[variant]
    return (energies(Z, N - 1)[k] - energies(Z, N)[k]) * HARTREE_EV
