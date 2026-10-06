"""Successive ionization energies of any atom or ion: the Screened Rydberg formula (final model pa_hier_rel;
see docs/unified.md).

    py -3.13 ionization.py Fe              all 26 successive IEs, predicted vs NIST
    py -3.13 ionization.py 26 --N 26       a single value (N = number of electrons before ionization)
    py -3.13 ionization.py 118             beyond the table: Madelung configuration, flagged
    py -3.13 ionization.py Fe --model uni_u35        use another model: pa_hier_rel [default, final],
                                          uni_u35, uni_u29, uni_pocket, uni_gshm_qed (aliases: final, u35, pocket)

    from ionization import ionization_energy, successive_ionization_energies
"""
import argparse
import os
import sys

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_ROOT, "models", "unified"))
sys.path.insert(0, _ROOT)
import final  # noqa: E402
from common.atomdata import get, ground_shells, shells_to_str  # noqa: E402

SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge "
           "As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm "
           "Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U "
           "Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split()


def to_Z(x):
    x = str(x).strip()
    if x.isdigit():
        return int(x)
    low = [s.lower() for s in SYMBOLS]
    if x.lower() not in low:
        raise ValueError(f"unknown element {x!r}")
    return low.index(x.lower()) + 1


def ionization_energy(Z, N=None, shells=None, model=None):
    """IE (eV) to remove one electron from the N-electron ion of element Z (default N = Z,
    i.e. the first IE of the neutral atom). Z may be an int or a symbol."""
    Z = to_Z(Z)
    N = Z if N is None else int(N)
    if not 1 <= N <= Z:
        raise ValueError("need 1 <= N <= Z")
    return final.predict(Z, N, shells, model)


def successive_ionization_energies(Z, model=None):
    """List of the Z successive IEs [IE_1 (N=Z), IE_2 (N=Z-1), ..., IE_Z (N=1)] in eV."""
    Z = to_Z(Z)
    return list(map(float, final.predict_many([(Z, N) for N in range(Z, 0, -1)], model=model)))


def _row(Z, N, pred):
    try:
        r = get(Z, N)
    except KeyError:
        r = None
    stage = Z - N + 1
    cfg = shells_to_str(ground_shells(Z, N))
    if r is None:
        return f"{stage:5d} {N:4d}  {cfg:34s} {pred:14.4f} {'-':>14s} {'-':>8s}  PREDICTION (not in NIST table)"
    err = (pred - r["IE_eV"]) / r["IE_eV"] * 100
    return f"{stage:5d} {N:4d}  {cfg:34s} {pred:14.4f} {r['IE_eV']:14.4f} {err:+8.2f}  {r['status']}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("element", help="symbol or Z")
    ap.add_argument("--N", type=int, help="number of electrons before ionization")
    ap.add_argument("--model", default=None, help=f"one of {final.MODELS} (default {final.DEFAULT_MODEL})")
    a = ap.parse_args()
    Z = to_Z(a.element)
    sym = SYMBOLS[Z - 1] if Z <= len(SYMBOLS) else f"Z{Z}"
    model = final.resolve(a.model)
    print(f"{sym} (Z = {Z}), model {model}; energies in eV")
    print(f"{'stage':>5s} {'N':>4s}  {'configuration (N electrons)':34s} {'predicted':>14s} {'NIST':>14s} {'err %':>8s}  status")
    if a.N is not None:
        print(_row(Z, a.N, ionization_energy(Z, a.N, model=model)))
        return
    ies = successive_ionization_energies(Z, model=model)
    for N, v in zip(range(Z, 0, -1), ies):
        print(_row(Z, N, v))
    if Z > 110:
        print("NOTE: Z > 110 is outside the NIST table: Madelung configurations, QED/FNS coefficients extrapolated beyond the Z<=110 table; "
              "all values are model predictions (extrapolation, larger uncertainty).")


if __name__ == "__main__":
    main()
