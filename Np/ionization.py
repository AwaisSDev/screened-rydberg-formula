"""Successive ionization energies of any atom or ion: the Screened Rydberg formula (final model pa_hier_rel;
see docs/unified.md).

    py -3.13 ionization.py Fe              all 26 successive IEs, predicted vs NIST
    py -3.13 ionization.py 26 --N 26       a single value (N = number of electrons before ionization)
    py -3.13 ionization.py 118             beyond the table: Madelung configuration, flagged
    py -3.13 ionization.py Fe --model uni_u35        use another model: pa_hier_rel [default, final],
                                          uni_u35, uni_u29, uni_pocket, uni_gshm_qed (aliases: final, u35, pocket)

    py -3.13 ionization.py O --N 8 --details          every intermediate quantity + low-confidence warnings
    py -3.13 ionization.py Na --N 11 --details --config "[Ne]3s1"

    from ionization import ie, ionization_energy, successive_ionization_energies
    ie("Fe")            -> dict with IE_eV, NIST_eV, sigma1, nu, T, h, D, Zeff, terms, warnings
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


_CORES = {"[He]": "1s2", "[Ne]": "1s2.2s2.2p6", "[Ar]": "1s2.2s2.2p6.3s2.3p6",
           "[Kr]": "1s2.2s2.2p6.3s2.3p6.3d10.4s2.4p6",
           "[Xe]": "1s2.2s2.2p6.3s2.3p6.3d10.4s2.4p6.4d10.5s2.5p6",
           "[Rn]": "1s2.2s2.2p6.3s2.3p6.3d10.4s2.4p6.4d10.5s2.5p6.4f14.5d10.6s2.6p6"}


def _parse_config(config):
    """'1s2.2s2.2p4', '[Ne].3s1', '[Ne]3s1' or [(n, l, occ), ...] -> [(n, l, occ), ...]."""
    if config is None or isinstance(config, (list, tuple)):
        return config
    from common.atomdata import parse_shells
    s = str(config).replace(" ", "")
    for core, full in _CORES.items():
        if s.startswith(core):
            s = full + "." + s[len(core):].lstrip(".")
    return parse_shells(s.strip("."))


def ie(Z, N=None, config=None):
    """Screened Rydberg formula (final model) with every intermediate quantity and confidence warnings.

    Z: atomic number or symbol; N: electrons before ionization (default Z = neutral atom);
    config: optional configuration of the N-electron ion ('1s2.2s2.2p4', '[Ne]3s1' or [(n, l, occ)]);
    default = NIST ground configuration (Madelung order outside the table).
    Returns a dict: IE_eV, NIST_eV (if tabulated), configuration, removed subshell, sigma1, nu groups, T, h, D,
    Zeff, Rydberg term, Dirac factor F, relativistic bracket R, Hund term, mu, QED/FNS term, warnings.
    """
    import numpy as np
    from configs import initial_final
    import umodel as U
    import gshm
    from push_a_import import PA
    from common.atomdata import RYDBERG_EV as RY, ALPHA
    Z = to_Z(Z)
    N = Z if N is None else int(N)
    if not 1 <= N <= Z:
        raise ValueError("need 1 <= N <= Z")
    shells = _parse_config(config)
    if shells is not None and sum(q for *_, q in shells) != N:
        raise ValueError(f"configuration has {sum(q for *_, q in shells)} electrons, expected N = {N}")
    spec, theta = final.final_params()
    a = U.build([(Z, N)], None if shells is None else [shells])
    P = PA._resolve(spec, theta)
    sN, _, rem, rearr = initial_final(Z, N, shells)
    n, l, j = float(a["n"][0]), int(a["l"][0]), float(a["j"][0])
    Za = max(Z - N + 1, 1)
    s1 = float(a["sigma1"][0])
    groups = ["same", "in", "core", "df", "out"]
    nu = PA._group_matrix(spec, a)[0]
    T = float(nu @ np.array([P["tau_" + g] for g in groups]))
    Cd = PA._dev_matrix(a)[0]
    T += sum(float(Cd[PA.DEV_IDX[c]] * P["dtau_" + c]) for c in spec.get("dev", []))
    h = (N - 1 - s1) if T >= 0 else s1
    kap = abs(P["kappa"]) + 0.05
    D = T / (Za + kap + abs(T) / h) if T else 0.0
    Ze = Z - s1 - D
    F = float(gshm.dirac_factor(n, j, Ze))
    R = 1.0 + P["rel_" + PA.REL_CLASSES[int(a["rc"][0])]] * (Z * ALPHA) ** 2 / n * (Ze / Za - 1.0)
    hund = RY * [0.0, P["x_p"], P["x_d"], P["x_f"]][l] * float(a["kink"][0]) / n ** 2
    mu, qed = float(a["mu"][0]), float(a["qedfns"][0]) * (Ze / Z) ** 2
    val = (RY * Ze ** 2 / n ** 2 * F * R + hund) * mu - qed
    ref = float(PA.evaluate(theta, spec, a)[0])
    assert abs(val / ref - 1) < 1e-9, (val, ref)
    w = []
    if Z >= 55 and Z - N <= 1:
        w.append("heavy neutral or singly charged ion (Z >= 55): weakest regime; neutral first IEs of Z >= 55 had "
                 "21.8 % mean error in the blind Z <= 54 -> Z >= 55 test")
    if Z == N:
        w.append("neutral atom: mean error on neutral first IEs is 7.5 % (all-data fit), up to ~30 %")
    if l == 3:
        w.append("f-electron removal: least accurate class (3.5 % mean error; 2-2.5x too high in extrapolation)")
    if Z > 110:
        w.append("Z > 110: outside the NIST table; Madelung configuration and extrapolated QED/finite-size terms")
    if rearr:
        w.append("rearranged configuration: the ion's ground configuration differs by more than one electron")
    if not 0 <= s1 <= N - 1:
        w.append("sigma1 outside [0, N-1]: the Za <= Zeff <= Z bound is not guaranteed")
    if shells is None and N > 1:
        nxt = final.predict(Z, N - 1)
        if nxt <= val:
            w.append(f"monotonicity violation: IE(N-1) = {nxt:.4f} eV is not larger than this IE")
    try:
        nist = get(Z, N)["IE_eV"] if shells is None else None
    except KeyError:
        nist = None
    L = "spdf"
    return {"Z": Z, "N": N, "IE_eV": val, "NIST_eV": nist,
            "error_pct": (val / nist - 1) * 100 if nist else None,
            "configuration": shells_to_str(sN), "removed": f"{rem[0]}{L[rem[1]]}", "j": j, "Za": Za,
            "sigma1": s1, "nu": dict(zip(groups, map(int, nu))), "T": T, "h": h, "kappa_used": kap, "D": D,
            "Zeff": Ze, "rydberg_eV": RY * Ze ** 2 / n ** 2, "F_dirac": F, "R_rel": R, "hund_eV": hund,
            "mu": mu, "qed_fns_eV": qed, "warnings": w}


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
    ap.add_argument("--details", action="store_true", help="with --N: every intermediate quantity + warnings")
    ap.add_argument("--config", default=None, help="with --details: configuration, e.g. '[Ne]3s1'")
    a = ap.parse_args()
    Z = to_Z(a.element)
    sym = SYMBOLS[Z - 1] if Z <= len(SYMBOLS) else f"Z{Z}"
    model = final.resolve(a.model)
    print(f"{sym} (Z = {Z}), model {model}; energies in eV")
    if not (a.N is not None and a.details):
        print(f"{'stage':>5s} {'N':>4s}  {'configuration (N electrons)':34s} {'predicted':>14s} {'NIST':>14s} {'err %':>8s}  status")
    if a.N is not None and a.details:
        d = ie(Z, a.N, a.config)
        for k, v in d.items():
            if k != "warnings":
                print(f"  {k:14s} {v:.6g}" if isinstance(v, float) else f"  {k:14s} {v}")
        for x in d["warnings"]:
            print("  WARNING:", x)
        return
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
