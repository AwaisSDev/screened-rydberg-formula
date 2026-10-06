r"""Exact first-order (1/Z) screening of the removed electron, from Track A.

For a configuration C of the N-electron ion and the removed subshell (n, l), define the
frozen-configuration ion C' = C minus one (n,l) electron.  Track A's exact first-order
energies E1(C), E1(C') (rational numbers from hydrogenic Slater integrals, Hund-term angular
algebra, Layzer-complex diagonalisation when one n-shell is open) give

    dE0 = 1/(2 n^2),   dE1 = E1(C') - E1(C),   sigma1 = -dE1 / (2 dE0) = -n^2 dE1

so that  Ry (Z - sigma1)^2 / n^2  reproduces the EXACT Z^2 and Z coefficients of the 1/Z
expansion of IE (in Rydberg units, 1 hartree = 2 Ry).

sigma1 depends only on the configuration (not on Z). Using the frozen configuration C' (instead
of the rearranged NIST ground state of the ion) keeps dE0 = 1/(2 n^2) also for the 63 rows whose
ground configuration reshuffles on ionization.
"""
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
for _p in (_ROOT, os.path.join(_ROOT, "models", "first_principles"),
           os.path.join(_ROOT, "models", "semi_empirical")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import zexp  # noqa: E402  (Track A)
from configs import initial_final, remove_one, normalize  # noqa: E402  (Track B)
from common.atomdata import shells_to_str  # noqa: E402

CACHE = os.path.join(_HERE, "cache", "sigma1.json")
_cache = None


def _load():
    global _cache
    if _cache is None:
        try:
            with open(CACHE, encoding="utf-8") as f:
                _cache = json.load(f)
        except Exception:
            _cache = {}
    return _cache


def save_cache():
    if _cache is not None:
        with open(CACHE, "w", encoding="utf-8") as f:
            json.dump(_cache, f, indent=0)


def sigma1_config(sN, rem):
    """Exact first-order screening constant of the (n,l)=rem electron in configuration sN."""
    sN = normalize(sN)
    key = shells_to_str(sN) + "|" + f"{rem[0]}{rem[1]}"
    c = _load()
    if key not in c:
        ion = remove_one(sN, tuple(rem))
        if ion:
            co = zexp.coeffs(1, sum(q for *_, q in sN), shells=sN, shells_ion=ion)
            c[key] = {"sigma1": float(co["sigma"]), "dE0": float(co["dE0"]),
                      "dE1": float(co["dE1"]), "complex": bool(co["complex"])}
        else:
            c[key] = {"sigma1": 0.0, "dE0": 0.5 / rem[0] ** 2, "dE1": 0.0, "complex": False}
    return c[key]["sigma1"]


def sigma1(Z, N, shells=None):
    """sigma1 for the electron removed from (Z, N) (same removal rule as Track B's GSHM)."""
    sN, _, rem, _ = initial_final(int(Z), int(N), shells)
    return sigma1_config(sN, rem)


if __name__ == "__main__":
    import time
    from common.atomdata import load_records
    t = time.time()
    out = [sigma1(r["Z"], r["N"]) for r in load_records()]
    save_cache()
    print(len(out), "rows", time.time() - t, "s")
    for Z, N in [(2, 2), (3, 3), (8, 8), (11, 11), (26, 26), (55, 55)]:
        print(Z, N, sigma1(Z, N))
