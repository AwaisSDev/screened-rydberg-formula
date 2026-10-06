"""Relativistic, nuclear and QED corrections for one-electron (hydrogen-like) binding energies.

All energies in hartree (positive = binding) unless stated.  c = 1/alpha.

Ingredients (each one is first-principles physics; literature input is limited to the
nuclear charge radii / masses and the Mohr-type one-loop self-energy function F_SE(Z alpha)):

1. Dirac point-nucleus energy (Sommerfeld fine-structure formula)
       E_nk = c^2 { [1 + (Z a / (n - |k| + sqrt(k^2 - Z^2 a^2)))^2]^(-1/2) - 1 }
2. Finite nuclear mass (leading relativistic recoil, Barker-Glover / Sapirstein-Yennie):
       B = mu c^2 (1-f) + mu^2 c^2 (1-f)^2 / (2 (M+m)),   mu = m M/(m+M),  f = 1 + E_nk/c^2
3. Finite nuclear size: numerical solution of the radial Dirac equation for a uniformly charged
   sphere, R_sph = sqrt(5/3) R_rms (R_rms from Yerokhin & Shabaev 2015 Table II, which uses
   Angeli & Marinova 2013; fallback R = 0.836 A^(1/3) + 0.570 fm), minus the point result.
4. Vacuum polarisation: Uehling potential (point nucleus)
       V_U(r) = -(2 a /(3 pi)) (Z/r) Int_1^inf dt e^{-2 c r t} (1 + 1/(2t^2)) sqrt(t^2-1)/t^2
   integrated with the exact Dirac 1s density; for ns, n>=2 the tabulated Uehling function is used.
5. One-loop self-energy:  dE_SE = (a/pi) (Z a)^4 / n^3 F_SE(Z a) m c^2,  F_SE(1s), F_SE(2s) from
   Mohr's all-order calculations as tabulated by Yerokhin & Shabaev, J. Phys. Chem. Ref. Data 44,
   033103 (2015) [arXiv:1506.01885], Table II, "SE(pnt)" rows.  For n >= 3 we use F_SE(2s).
   (Self-energy of p and d electrons is an order of magnitude smaller and is neglected.)
"""
import json
import math
import os
import sys

import numpy as np
from numba import njit
from scipy.optimize import brentq
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
from common.atomdata import ALPHA, ELECTRON_MASS_U  # noqa: E402

C = 1.0 / ALPHA
BOHR_FM = 52917.721090  # bohr radius in fm
_QED = {int(k): v for k, v in json.load(open(os.path.join(HERE, "cache", "yerokhin_shabaev_2015_qed.json"))).items()}


def kappa_of(l, j2):
    """kappa from l and 2j."""
    return -(l + 1) if j2 == 2 * l + 1 else l


def dirac_point(n, kappa, Z):
    """Dirac point-nucleus energy (negative, hartree, infinite nuclear mass)."""
    za = Z * ALPHA
    k = abs(kappa)
    if za >= k:      # point-nucleus Dirac equation has no bound solution (Z alpha >= |kappa|)
        za = k * (1.0 - 1e-12)  # clamp to the critical value (only reached for Z >= 138, 1s/2p1/2)
    g = math.sqrt(k * k - za * za)
    f = 1.0 / math.sqrt(1.0 + (za / (n - k + g)) ** 2)
    return C * C * (f - 1.0)


def nuclear_params(Z):
    d = _QED.get(Z)
    if d is None:
        A = round(2.5 * Z)
        return 0.836 * A ** (1 / 3) + 0.570, A * 1822.888 - Z
    return d["R_fm"], d["M_over_m"]


# --------------------------------------------------------------------------- numerical Dirac
@njit(cache=True)
def _V_sphere(r, Z, Rn):
    if r >= Rn:
        return -Z / r
    return -Z / (2 * Rn) * (3.0 - (r / Rn) ** 2)


@njit(cache=True)
def _deriv(r, P, Q, E, kappa, Z, Rn, c):
    V = _V_sphere(r, Z, Rn)
    dP = -kappa / r * P + (E - V + 2 * c * c) / c * Q
    dQ = kappa / r * Q - (E - V) / c * P
    return dP, dQ


@njit(cache=True)
def _shoot(E, kappa, Z, Rn, c, x0, x1, xm, npts):
    """Integrate outward (x0->xm) and inward (x1->xm) in x = ln r with RK4. Return Q/P mismatch."""
    h = (xm - x0) / npts
    r = math.exp(x0)
    l = kappa if kappa > 0 else -kappa - 1
    P = r ** (l + 1)
    V0 = _V_sphere(r, Z, Rn)
    if kappa < 0:
        Q = -(E - V0) * r ** (l + 2) / ((2 * l + 3) * c)
    else:
        Q = (2 * l + 1) * c * r ** l / (E - V0 + 2 * c * c)
    x = x0
    for i in range(npts):
        # RK4 in x: dP/dx = r dP/dr
        r = math.exp(x)
        a1, b1 = _deriv(r, P, Q, E, kappa, Z, Rn, c)
        a1 *= r; b1 *= r
        r2 = math.exp(x + h / 2)
        a2, b2 = _deriv(r2, P + h / 2 * a1, Q + h / 2 * b1, E, kappa, Z, Rn, c)
        a2 *= r2; b2 *= r2
        a3, b3 = _deriv(r2, P + h / 2 * a2, Q + h / 2 * b2, E, kappa, Z, Rn, c)
        a3 *= r2; b3 *= r2
        r4 = math.exp(x + h)
        a4, b4 = _deriv(r4, P + h * a3, Q + h * b3, E, kappa, Z, Rn, c)
        a4 *= r4; b4 *= r4
        P += h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        Q += h / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
        x += h
        if abs(P) > 1e100:
            P *= 1e-100; Q *= 1e-100
    ratio_out = Q / P
    # inward
    lam = math.sqrt(-E * (2 + E / (c * c)))
    r = math.exp(x1)
    P = 1e-30
    Q = -P * math.sqrt(-E / (E + 2 * c * c))
    h = (xm - x1) / npts
    x = x1
    for i in range(npts):
        r = math.exp(x)
        a1, b1 = _deriv(r, P, Q, E, kappa, Z, Rn, c)
        a1 *= r; b1 *= r
        r2 = math.exp(x + h / 2)
        a2, b2 = _deriv(r2, P + h / 2 * a1, Q + h / 2 * b1, E, kappa, Z, Rn, c)
        a2 *= r2; b2 *= r2
        a3, b3 = _deriv(r2, P + h / 2 * a2, Q + h / 2 * b2, E, kappa, Z, Rn, c)
        a3 *= r2; b3 *= r2
        r4 = math.exp(x + h)
        a4, b4 = _deriv(r4, P + h * a3, Q + h * b3, E, kappa, Z, Rn, c)
        a4 *= r4; b4 *= r4
        P += h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        Q += h / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
        x += h
        if abs(P) > 1e100:
            P *= 1e-100; Q *= 1e-100
    return ratio_out - Q / P


def dirac_numeric(n, kappa, Z, R_fm, npts=20000):
    """Energy of the (n, kappa) level for a uniformly charged sphere of rms radius R_fm.
    Only nodeless states (n = |kappa| for kappa<0 ... i.e. 1s, 2p3/2, 3d5/2) and 2s/2p1/2 are
    safely bracketed around the point-nucleus value."""
    Rn = math.sqrt(5.0 / 3.0) * R_fm / BOHR_FM
    E0 = dirac_point(n, kappa, Z)
    x0 = math.log(Rn * 1e-4)
    xm = math.log(1.5 * n * n / Z)
    x1 = math.log((n * n + 25 * n) / Z)
    f = lambda E: _shoot(E, kappa, Z, Rn, C, x0, x1, xm, npts)
    lo, hi = E0 * 1.02, E0 * 0.98
    return brentq(f, lo, hi, xtol=1e-14 * abs(E0), rtol=1e-15)


def fns_shift(n, kappa, Z, R_fm=None):
    """Finite-nuclear-size energy shift (>0, hartree): E(sphere) - E(point), both numerical."""
    if R_fm is None:
        R_fm, _ = nuclear_params(Z)
    e_fin = dirac_numeric(n, kappa, Z, R_fm)
    e_pt = dirac_numeric(n, kappa, Z, 1e-6)
    return e_fin - e_pt


# --------------------------------------------------------------------------- QED
def uehling_1s(Z):
    """<V_Uehling> (hartree, negative) for the Dirac point-nucleus 1s state, computed here."""
    za = Z * ALPHA
    g = math.sqrt(1 - za * za)
    lam = 2 * Z

    # density rho(r) = lam^(2g+1)/Gamma(2g+1) r^(2g) e^(-lam r) (P^2+Q^2 normalised)
    # <V_U> = -(2a/3pi) Z Int_1^inf dt w(t) Int_0^inf rho(r)/r e^{-2 c t r} dr
    #       = -(2a/3pi) Z Int dt w(t) lam^(2g+1) Gamma(2g)/Gamma(2g+1) / (lam + 2 c t)^(2g)
    def w(t):
        return (1 + 1 / (2 * t * t)) * math.sqrt(t * t - 1) / (t * t)

    def integrand(t):
        return w(t) * lam ** (2 * g + 1) / (2 * g) / (lam + 2 * C * t) ** (2 * g)
    val = quad(integrand, 1, np.inf, limit=200, epsabs=0, epsrel=1e-12)[0]
    return -(2 * ALPHA / (3 * math.pi)) * Z * val


def qed_prefactor(Z, n):
    """(alpha/pi) (Z alpha)^4 / n^3 * m c^2  in hartree."""
    return ALPHA / math.pi * (Z * ALPHA) ** 4 / n ** 3 * C * C


_ZMAX_QED = max(_QED)
EXTRAPOLATED_QED_WARNED = set()


def _qed_table(Z, key):
    """Table value of F(Z alpha) for key; for Z beyond the YS15 table (Z > 110) a quadratic
    extrapolation through the last three tabulated points (Z = 108-110) is returned
    (flagged: such values are extrapolated, not literature values)."""
    Z = int(Z)
    if Z in _QED:
        return _QED[Z][key]
    if Z > _ZMAX_QED:
        zs = [_ZMAX_QED - 2, _ZMAX_QED - 1, _ZMAX_QED]
        c = np.polyfit(zs, [_QED[z][key] for z in zs], 2)
        if Z not in EXTRAPOLATED_QED_WARNED:
            EXTRAPOLATED_QED_WARNED.add(Z)
        return float(np.polyval(c, Z))
    raise KeyError(f"no QED table entry for Z={Z}")


def F_SE(Z, n):
    return _qed_table(Z, "SE1s" if n == 1 else "SE2s")


def F_U(Z, n):
    return _qed_table(Z, "Ue1s" if n == 1 else "Ue2s")


def qed_shift(Z, n, l, computed_uehling=True):
    """One-loop QED shift of an ns level (hartree, >0 means less bound). Zero for l > 0."""
    if l > 0:
        return 0.0
    se = qed_prefactor(Z, n) * F_SE(Z, n)
    if n == 1 and computed_uehling:
        vp = uehling_1s(Z)
    else:
        vp = qed_prefactor(Z, n) * F_U(Z, n)
    return se + vp


# --------------------------------------------------------------------------- H-like total
_FNS_CACHE_FILE = os.path.join(HERE, "cache", "fns_shifts.json")
try:
    _FNS = json.load(open(_FNS_CACHE_FILE))
except Exception:
    _FNS = {}


def fns_cached(n, kappa, Z):
    key = f"{n},{kappa},{Z}"
    if key not in _FNS:
        _FNS[key] = fns_shift(n, kappa, Z)
        json.dump(_FNS, open(_FNS_CACHE_FILE, "w"), indent=0)
    return _FNS[key]


def hydrogenic_binding(Z, n=1, l=0, j2=None, recoil=True, fns=True, qed=True):
    """Binding energy (hartree, >0) of a one-electron ion in level (n, l, j)."""
    if j2 is None:
        j2 = 2 * l + 1 if l == 0 else 2 * l - 1
    kappa = kappa_of(l, j2)
    E = dirac_point(n, kappa, Z)          # negative
    if fns and kappa == -1:
        E += fns_cached(n, kappa, Z)
    if recoil:
        _, Mm = nuclear_params(Z)
        mu = Mm / (1 + Mm)
        fm1 = E / (C * C)                  # f - 1
        E = mu * C * C * fm1 - mu * mu * C * C * fm1 ** 2 / (2 * (Mm + 1))
    if qed:
        E += qed_shift(Z, n, l)
    return -E


def dirac_correction(n, l, j2, Z):
    """Relativistic part of the point-nucleus Dirac binding: B_D - Z^2/(2 n^2) (hartree, >0)."""
    return -dirac_point(n, kappa_of(l, j2), Z) - Z * Z / (2.0 * n * n)


if __name__ == "__main__":
    # validation: Uehling vs Yerokhin-Shabaev table
    for Z in (1, 10, 50, 92, 110):
        print(Z, "F_U computed", uehling_1s(Z) / qed_prefactor(Z, 1), "table", F_U(Z, 1))
    # numerical Dirac point vs analytic
    for Z in (1, 50, 92):
        e = dirac_numeric(1, -1, Z, 1e-6)
        print(Z, "numeric pt", e, "analytic", dirac_point(1, -1, Z), "rel", e / dirac_point(1, -1, Z) - 1)
    HARTREE = 27.211386245988
    for Z in (1, 20, 50, 82, 92):
        print(Z, "FNS eV", fns_shift(1, -1, Z) * HARTREE)


# --------------------------------------------------------------------------- closed-form self-energy
LN_K0 = {1: 2.984128556, 2: 2.811769893}   # Bethe logarithms ln k0(ns)


def F_SE_closed(Z, n=1):
    """Low-order Z-alpha expansion of the one-loop self-energy function (no table):
         F = A41 ln(Za)^-2 + A40 + A50 (Za) + (Za)^2 [A62 ln^2(Za)^-2 + A61 ln(Za)^-2]
       A41 = 4/3, A40 = 10/9 - (4/3) ln k0(ns), A50 = 4 pi (139/128 - ln2 / 2), A62 = -1,
       A61(1s) = 28/3 ln2 - 21/20 (the (Za)^2 terms are used for n = 1 only)."""
    za = Z * ALPHA
    L = math.log(za ** -2)
    F = 4.0 / 3.0 * L + (10.0 / 9.0 - 4.0 / 3.0 * LN_K0[min(n, 2)]) \
        + 4 * math.pi * (139.0 / 128.0 - math.log(2) / 2) * za
    if n == 1:
        F += za ** 2 * (-L * L + (28.0 / 3.0 * math.log(2) - 21.0 / 20.0) * L)
    return F


def hydrogenic_binding_closed(Z, n=1, l=0):
    """Fully closed-form/first-principles variant: Dirac + recoil + numerical FNS + computed Uehling
    + Z-alpha-expansion self-energy (no tabulated QED input)."""
    E = -hydrogenic_binding(Z, n, l, qed=False)
    if l == 0:
        E += qed_prefactor(Z, n) * F_SE_closed(Z, n) + (uehling_1s(Z) if n == 1 else qed_prefactor(Z, n) * F_U(Z, n))
    return -E
