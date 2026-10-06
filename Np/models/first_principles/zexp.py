"""Analytic 1/Z-expansion ionization-energy models (Track A, first principles).

    IE(Z, N) = E(Z, N-1) - E(Z, N)
             = Z^2 dE0 + Z dE1 + dE2 + dE3/Z + ...                      [hartree]
    dE0 = E0(ion) - E0(atom)          (= 1/(2 n^2) for removal of an n electron: Bohr)
    dE1 = E1(ion) - E1(atom)          (exact first-order Coulomb repulsion, <0)

The leading-order ab initio screening constant of the removed electron follows from completing
the square, IE ~ dE0 (Z - sigma)^2 with sigma = -dE1 / (2 dE0)  (= -n^2 dE1 when dE0 = 1/(2n^2)).

Relativistic + QED correction for the removed electron (n, l, j):
    dR = [ B_Dirac(n,kappa,Z) - Z^2/(2n^2) ] * s(Z, sigma)  + QED_ns(Z) * s(Z, sigma)
with s = ((Z - sigma)/Z)^p  ('screening' of the relativistic correction).  p = 2 is the
Lande / penetrating-orbital estimate, p = 4 the pure-screening (hydrogenic Z_eff) estimate.

Public API
----------
coeffs(Z, N, shells=None, shells_ion=None) -> dict with dE0, dE1 (best), dE1_sc, dE1_cx, sigma, n, l, j2
predict(Z, N, shells=None, variant="zexp_sq_rel") -> IE in eV
"""
import json
import os
import sys
from fractions import Fraction
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

from common.atomdata import HARTREE_EV, ground_shells, shells_to_str  # noqa: E402
import config_energy as ce  # noqa: E402
import relativity as rel  # noqa: E402

CFG_CACHE = os.path.join(HERE, "cache", "config_E1.json")
_cfg = None


def _load_cfg():
    global _cfg
    if _cfg is None:
        _cfg = json.load(open(CFG_CACHE)) if os.path.exists(CFG_CACHE) else {}
    return _cfg


def save_cfg_cache():
    if _cfg is not None:
        json.dump(_cfg, open(CFG_CACHE, "w"), indent=0)
        ce.save_cache()


def config_data(shells):
    """E0, E_av, E1_sc (exact string if rational), E1_cx for a configuration (cached on disk)."""
    shells = tuple(sorted(tuple(s) for s in shells if s[2] > 0))
    key = shells_to_str(shells) if shells else "-"
    c = _load_cfg()
    if key not in c:
        d = ce.config_coeffs(shells)
        c[key] = {"E0": str(d["E0"]), "E_av": str(d["E_av"]),
                  "E1_sc": str(d["E1_sc"]) if isinstance(d["E1_sc"], Fraction) else repr(float(d["E1_sc"])),
                  "sc_exact": bool(d["sc_exact"]),
                  "E1_cx": d["E1_cx"], "cx_ref_weight": d["cx_ref_weight"]}
    return c[key]


def _f(s):
    return float(Fraction(s)) if "/" in s or s.lstrip("-").isdigit() else float(s)


def removed_subshell(shells, shells_ion):
    a = {(n, l): q for n, l, q in shells}
    b = {(n, l): q for n, l, q in shells_ion}
    lost = {k: a.get(k, 0) - b.get(k, 0) for k in set(a) | set(b)}
    lost = {k: v for k, v in lost.items() if v > 0}
    k = max(lost, key=lambda k: (lost[k], k[0], k[1]))
    return k, a[k]


def ion_shells_default(shells):
    """Remove one electron from the least-bound subshell at first order (highest n, then l)."""
    sh = sorted(shells, key=lambda s: (s[0], s[1]))
    n, l, q = sh[-1]
    out = [list(s) for s in sh]
    out[-1][2] -= 1
    return [tuple(s) for s in out if s[2] > 0]


def coeffs(Z, N, shells=None, shells_ion=None):
    if shells is None:
        shells = ground_shells(Z, N)
    if shells_ion is None:
        shells_ion = ground_shells(Z, N - 1) if N > 1 else []
        if sum(q for *_, q in shells_ion) != N - 1:
            shells_ion = ion_shells_default(shells)
    a = config_data(shells)
    b = config_data(shells_ion) if shells_ion else {"E0": "0", "E_av": "0", "E1_sc": "0", "E1_cx": 0.0}
    dE0 = _f(b["E0"]) - _f(a["E0"])
    dE1_sc = _f(b["E1_sc"]) - _f(a["E1_sc"])
    dE1_av = _f(b["E_av"]) - _f(a["E_av"])
    if a["E1_cx"] is not None and b["E1_cx"] is not None:
        dE1_cx = b["E1_cx"] - a["E1_cx"]
    else:
        dE1_cx = None
    (n, l), q = removed_subshell(shells, shells_ion)
    j2 = 1 if l == 0 else (2 * l - 1 if q <= 2 * l else 2 * l + 1)
    best = dE1_cx if dE1_cx is not None else dE1_sc
    return {"dE0": dE0, "dE1_sc": dE1_sc, "dE1_cx": dE1_cx, "dE1_av": dE1_av, "dE1": best,
            "sigma": -best / (2 * dE0) if dE0 > 0 else float("nan"),
            "n": n, "l": l, "j2": j2, "complex": dE1_cx is not None}


def rel_correction(Z, c, sigma=None, p=2.0, qed=True, fns=True):
    """Relativistic+QED+FNS correction (hartree, added to IE) for the removed electron."""
    n, l, j2 = c["n"], c["l"], c["j2"]
    if sigma is None:
        sigma = c["sigma"]
    if not (sigma == sigma) or c["dE0"] <= 0:   # nan (rearranged, dE0<=0): no screening of dR
        sigma = 0.0
    zo = min(max(Z - sigma, 1.0), Z)
    s = (zo / Z) ** p
    d = rel.dirac_correction(n, l, j2, Z)
    if qed and l == 0:
        d -= rel.qed_shift(Z, n, l)
    if fns and l == 0 and n <= 2:
        d -= rel.fns_cached(n, -1, Z)
    return d * s


def reduced_mass_factor(Z):
    _, Mm = rel.nuclear_params(Z)
    return Mm / (1.0 + Mm)


def predict(Z, N, shells=None, variant="zexp_sq_rel", p=2.0):
    """IE in eV.  Variants (all zero fitted parameters):
       bohr        : Z^2 dE0
       zexp2       : Z^2 dE0 + Z dE1                 (exact two-term series, non-relativistic)
       zexp2_rel   : zexp2 + relativistic/QED correction of removed electron
       zexp_sq     : dE0 (Z - sigma)^2  (closed-form completion of the square, dE2 ~ dE1^2/(4 dE0))
       zexp_sq_rel : zexp_sq + relativistic/QED correction
    """
    c = coeffs(Z, N, shells)
    if variant == "bohr":
        return Z * Z * c["dE0"] * HARTREE_EV
    nr = Z * Z * c["dE0"] + Z * c["dE1"]
    if variant.startswith("zexp_sq"):
        nr += c["dE1"] ** 2 / (4 * c["dE0"]) if c["dE0"] > 0 else 0.0
    if variant.endswith("_rel"):
        nr = nr * reduced_mass_factor(Z) + rel_correction(Z, c, p=p)
    return nr * HARTREE_EV
