"""Textbook baselines for successive ionization energies IE(Z, N).

All functions have the signature  f(Z, N, shells=None) -> IE in eV.

  bohr                 IE = Ry Z^2 / n^2                         (n of removed electron)
  slater_1e            IE = Ry (Z_eff / n*)^2                    (Slater 1930 rules)
  slater_total         IE = E_S(N-1) - E_S(N),  E_S = -Ry sum_i q_i (Z_eff,i / n*_i)^2
  cr_1e / cr_total     same two variants with Clementi-Raimondi (1963) screening constants
                       (true principal quantum number n, as in their Slater-type orbitals)

Clementi-Raimondi rules (J. Chem. Phys. 38, 2686 (1963)) cover 1s..4p only.  Extension used
here (documented in docs/semi_empirical.md):
  * In the 1s and 2s formulas the 'N3spd + N4sp' outer-electron term is applied to ALL
    electrons with n >= 3.
  * Subshells beyond the CR range are mapped onto the CR formula for the same l with the
    largest available n (s -> 4s, p -> 4p, d -> 3d), shifting every principal quantum number
    by Delta = n - n_analog.  Populations in the explicit CR terms are read from the shifted
    configuration, the CR constant stands for the analog's (full) core, and every electron
    that maps to an invalid orbital (n' - Delta < 1 or l' >= n' - Delta) and is inside the
    target (n' < n, or n' = n and l' < l) screens fully (1.0, Slater's deep-core value).
  * f electrons (no CR formula): same-subshell 0.2693 (the CR d-d value), every electron
    inside (n' < n, or n' = n and l' < 3) 1.0, outer electrons 0.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
for p in (_ROOT, _HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

from common.atomdata import RYDBERG_EV, load_records  # noqa: E402
from configs import initial_final  # noqa: E402

RY = RYDBERG_EV

# ----------------------------------------------------------------------------- Bohr


def bohr(Z, N, shells=None):
    _, _, (n, l), _ = initial_final(Z, N, shells)
    return RY * Z * Z / (n * n)

# ----------------------------------------------------------------------------- Slater


SLATER_NSTAR = {1: 1.0, 2: 2.0, 3: 3.0, 4: 3.7, 5: 4.0, 6: 4.2, 7: 4.3}  # n=7: 4.3 (our choice)


def _slater_group(n, l):
    return (n, 0 if l <= 1 else l - 1)  # (n, sp=0 | d=1 | f=2), ordered as Slater's groups


def slater_sigma(target, shells):
    n, l = target
    gt = _slater_group(n, l)
    s = 0.0
    for n2, l2, q in shells:
        cnt = q - 1 if (n2, l2) == (n, l) else q
        if cnt <= 0:
            continue
        g = _slater_group(n2, l2)
        if g == gt:
            s += cnt * (0.30 if n == 1 else 0.35)
        elif g < gt:
            if l <= 1:
                s += cnt * (0.85 if n2 == n - 1 else 1.0)
            else:
                s += cnt * 1.0
        # groups to the right: 0
    return s


def slater_orbital_energy(target, shells, Z):
    n, _ = target
    zeff = Z - slater_sigma(target, shells)
    return -RY * (zeff / SLATER_NSTAR[n]) ** 2


def slater_total_energy(shells, Z):
    return sum(q * slater_orbital_energy((n, l), shells, Z) for n, l, q in shells)


def slater_1e(Z, N, shells=None):
    sN, _, rem, _ = initial_final(Z, N, shells)
    return -slater_orbital_energy(rem, sN, Z)


def slater_total(Z, N, shells=None):
    sN, sF, _, _ = initial_final(Z, N, shells)
    return slater_total_energy(sF, Z) - slater_total_energy(sN, Z)

# ----------------------------------------------------------------------------- Clementi-Raimondi


def _cr_formula(t, P):
    """Original CR formulas. t = (n,l) in the CR range; P(n,l) -> population."""
    n3up = None
    if t == (1, 0):
        return 0.3 * (P(1, 0) - 1) + 0.0072 * (P(2, 0) + P(2, 1)) + 0.0158 * P("n>=3")
    if t == (2, 0):
        return 1.7208 + 0.3601 * (P(2, 0) - 1 + P(2, 1)) + 0.2062 * P("n>=3")
    if t == (2, 1):
        return (2.5787 + 0.3326 * (P(2, 1) - 1) - 0.0773 * P(3, 0) - 0.0161 * (P(3, 1) + P(4, 0))
                - 0.0048 * P(3, 2) + 0.0085 * P(4, 1))
    if t == (3, 0):
        return (8.4927 + 0.2501 * (P(3, 0) - 1 + P(3, 1)) + 0.0778 * P(4, 0) + 0.3382 * P(3, 2)
                + 0.1978 * P(4, 1))
    if t == (3, 1):
        return 9.3345 + 0.3803 * (P(3, 1) - 1) + 0.0526 * P(4, 0) + 0.3289 * P(3, 2) + 0.1558 * P(4, 1)
    if t == (4, 0):
        return 15.505 + 0.0971 * (P(4, 0) - 1) + 0.8433 * P(3, 2) + 0.0687 * P(4, 1)
    if t == (3, 2):
        return 13.5894 + 0.2693 * (P(3, 2) - 1) - 0.1065 * P(4, 1)
    if t == (4, 1):
        return 24.7782 + 0.2905 * (P(4, 1) - 1)
    raise KeyError(t)
    return n3up


_CR_RANGE = {(1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (3, 2), (4, 1)}
_CR_ANALOG = {0: 4, 1: 4, 2: 3}


def cr_sigma(target, shells):
    n, l = target
    pop = {(a, b): q for a, b, q in shells}
    if l == 3:
        s = 0.2693 * (pop.get(target, 0) - 1)
        s += sum(q for a, b, q in shells if a < n or (a == n and b < 3))
        return s
    if target in _CR_RANGE:
        def P(a, b=None):
            if a == "n>=3":
                return sum(q for x, y, q in shells if x >= 3)
            return pop.get((a, b), 0)
        return _cr_formula(target, P)
    nA = _CR_ANALOG[l]
    D = n - nA
    deep = 0
    for a, b, q in shells:
        am = a - D
        if (am < 1 or b >= am) and (a < n or (a == n and b < l)):
            deep += q

    def P(a, b=None):
        if a == "n>=3":
            return sum(q for x, y, q in shells if x - D >= 3)
        return pop.get((a + D, b), 0)
    return _cr_formula((nA, l), P) + deep


def cr_orbital_energy(target, shells, Z):
    n, _ = target
    zeff = Z - cr_sigma(target, shells)
    return -RY * (zeff / n) ** 2


def cr_total_energy(shells, Z):
    return sum(q * cr_orbital_energy((n, l), shells, Z) for n, l, q in shells)


def cr_1e(Z, N, shells=None):
    sN, _, rem, _ = initial_final(Z, N, shells)
    return -cr_orbital_energy(rem, sN, Z)


def cr_total(Z, N, shells=None):
    sN, sF, _, _ = initial_final(Z, N, shells)
    return cr_total_energy(sF, Z) - cr_total_energy(sN, Z)


BASELINES = {
    "bohr": bohr,
    "slater_1e": slater_1e,
    "slater_total": slater_total,
    "cr_1e": cr_1e,
    "cr_total": cr_total,
}


def main():
    import csv
    import json
    sys.path.insert(0, _ROOT)
    import evaluate
    recs = load_records()
    summary = {}
    for name, f in BASELINES.items():
        path = os.path.join(_ROOT, "results", f"se_{name}_predictions.csv")
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["Z", "N", "IE_pred_eV", "IE_ref_eV"])
            for r in recs:
                w.writerow([r["Z"], r["N"], f"{f(r['Z'], r['N']):.8g}", r["IE_eV"]])
        m = evaluate.evaluate(path)
        with open(os.path.join(_ROOT, "results", f"se_{name}_metrics.json"), "w") as fh:
            json.dump(m, fh, indent=1)
        evaluate.print_table({k: v for k, v in m.items()
                              if k in ("ALL", "first_IE_neutral_atoms", "hydrogen_like", "_coverage")}, name)
        summary[name] = {k: m[k] for k in ("ALL", "first_IE_neutral_atoms", "hydrogen_like")}
    with open(os.path.join(_ROOT, "results", "se_baselines_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)


if __name__ == "__main__":
    main()
