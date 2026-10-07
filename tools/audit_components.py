r"""Component-by-component audit of the Screened Rydberg formula (paper §2.2, Table A1).

For each requested (Z, N) and each model it prints every intermediate quantity of the paper's equations,
recomputes the final IE from those printed components alone, and asserts that the result equals the
production code (models/push_a/model.py, models/unified/pocket.py) to 1e-9 relative.

    py -3.11 tools/audit_components.py                 # O, Na, N, Ca, Fe, Pb (neutral first IEs)
    py -3.11 tools/audit_components.py 26:25 82:80     # any Z:N pairs

Models:
  pa_hier_rel  final model, 33 parameters (paper "Screened Rydberg formula")
  pa_bound9    9-parameter bounded variant (same equations: 5 group tau_g, kappa, x_p, x_d, x_f;
               no class deviations, no relativistic penetration bracket)
  pocket       8-parameter pocket formula (Appendix A, results/uni_params.json "pocket")

Writes results/audit_components.json and results/audit_components.md.
Notation: D_rem is the bounded screening remainder D of eq. (2.2); F_Dirac is the Dirac/Schroedinger
ratio D_{n,j}(Zeff) (the paper uses the letter D for both).
"""
import json
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _p in (_ROOT, os.path.join(_ROOT, "models", "unified"), os.path.join(_ROOT, "models", "first_principles"),
           os.path.join(_ROOT, "models", "semi_empirical")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import umodel as U  # noqa: E402
import gshm  # noqa: E402
import pocket  # noqa: E402
from push_a_import import PA  # noqa: E402
from common.atomdata import RYDBERG_EV as RY, ALPHA, load_records  # noqa: E402

SYMBOL = {1: "H", 7: "N", 8: "O", 11: "Na", 20: "Ca", 26: "Fe", 82: "Pb"}
LL = "spdf"
GROUPS = ["same", "in", "core", "df", "out"]
DEFAULT = [(8, 8), (11, 11), (7, 7), (20, 20), (26, 26), (82, 82)]


def _nist():
    return {(r["Z"], r["N"]): r for r in load_records()}


def _pa_candidate(name):
    with open(os.path.join(_ROOT, "results", "pa_params.json"), encoding="utf-8") as f:
        c = json.load(f)["candidates"][name]
    return c["spec"], np.array(c["values"], float)


def _cfg_str(shells):
    return ".".join(f"{n}{LL[l]}{q}" for n, l, q in shells)


def audit_pa(name, Z, N, nist):
    spec, theta = _pa_candidate(name)
    a = U.build([(Z, N)])
    P = PA._resolve(spec, theta)
    i = 0
    Zf, n, l, j = float(a["Z"][i]), float(a["n"][i]), int(a["l"][i]), float(a["j"][i])
    Za = max(float(a["Zq"][i]), 1.0)
    s1 = float(a["sigma1"][i])
    nu = PA._group_matrix(spec, a)[i]
    tau = np.array([P["tau_" + g] for g in GROUPS])
    T_groups = float(nu @ tau)
    dev = {}
    if spec.get("dev"):
        Cd = PA._dev_matrix(a)[i]
        for c in spec["dev"]:
            v = Cd[PA.DEV_IDX[c]]
            if v:
                dev[c] = dict(nu=float(v), dtau=float(P["dtau_" + c]), contrib=float(v * P["dtau_" + c]))
    T = T_groups + sum(d["contrib"] for d in dev.values())
    h = (N - 1 - s1) if T >= 0 else s1
    kap = abs(P["kappa"]) + 0.05
    denom = Za + kap + abs(T) / h
    D = T / denom
    Ze = Zf - s1 - D
    ryd = RY * Ze ** 2 / n ** 2
    FD = float(gshm.dirac_factor(n, j, Ze))
    rc_name = PA.REL_CLASSES[int(a["rc"][i])]
    if spec.get("rel") == "fs":
        r_c = float(P["rel_" + rc_name])
        R = 1.0 + r_c * (Zf * ALPHA) ** 2 / n * (Ze / Za - 1.0)
    else:
        r_c, R = None, 1.0
    hyd = ryd * FD * R
    x_l = [0.0, P["x_p"], P["x_d"], P["x_f"]][l]
    K = float(a["kink"][i])
    hund = RY * x_l * K / n ** 2
    mu = float(a["mu"][i])
    qed = float(a["qedfns"][i]) * (Ze / Zf) ** 2
    ie = (hyd + hund) * mu - qed
    code = float(PA.evaluate(theta, spec, a)[0])
    assert abs(ie / code - 1) < 1e-9, (name, Z, N, ie, code)
    ref = nist.get((Z, N))
    return dict(model=name, Z=Z, N=N, n=n, l=LL[l], j=j, k=float(a["k"][i]), Za=Za,
                nu_groups=dict(zip(GROUPS, map(float, nu))),
                tau_groups=dict(zip(GROUPS, map(float, tau))), T_groups=T_groups, class_devs=dev,
                sigma1=s1, T=T, h=h, h_rule="(N-1)-sigma1 (T>=0)" if T >= 0 else "sigma1 (T<0)",
                kappa_fit=float(P["kappa"]), kappa_used=kap, denominator=denom, D_rem=D,
                total_screening=s1 + D, Zeff=Ze, bound_ok=bool(Za <= Ze <= Zf),
                rydberg_term=ryd, F_Dirac=FD, rel_class=rc_name, r_c=r_c, rel_bracket=R,
                hydrogenic_after_rel=hyd, x_l=x_l, K=K, hund_term=hund, mu=mu, qed_fns_term=qed,
                IE=ie, IE_code=code, IE_NIST=ref["IE_eV"] if ref else None,
                err_pct=(ie / ref["IE_eV"] - 1) * 100 if ref else None)


def audit_pocket(Z, N, nist):
    with open(os.path.join(_ROOT, "results", "uni_params.json"), encoding="utf-8") as f:
        th = np.array(json.load(f)["pocket"]["values"], float)
    a = U.build([(Z, N)])
    i = 0
    s, t, kap, x = th[:5], th[5], abs(th[6]) + 0.05, th[7]
    nu = pocket.counts(a)[i]
    Zf, n, j = float(a["Z"][i]), float(a["n"][i]), float(a["j"][i])
    Za = max(float(a["Zq"][i]), 1.0)
    S_groups = float(nu @ s)
    S_charge = t * (N - 1) / (Za + kap)
    Ze = Zf - S_groups - S_charge
    ryd = RY * Ze ** 2 / n ** 2
    R = 1.0 + (Ze * ALPHA / n) ** 2 * (n / (j + 0.5) - 0.75)
    K = float(a["kink"][i])
    hund = RY * x * K / n ** 2
    ie = ryd * R + hund
    code = float(pocket.evaluate(th, a)[0])
    assert abs(ie / code - 1) < 1e-9, ("pocket", Z, N, ie, code)
    ref = nist.get((Z, N))
    return dict(model="pocket", Z=Z, N=N, n=n, l=LL[int(a["l"][i])], j=j, Za=Za,
                nu_groups=dict(zip(GROUPS, map(float, nu))), s_groups=dict(zip(GROUPS, map(float, s))),
                screening_groups=S_groups, t=float(t), kappa_used=float(kap), charge_term=float(S_charge),
                total_screening=S_groups + S_charge, Zeff=Ze, rydberg_term=ryd, sommerfeld_bracket=R,
                x=float(x), K=K, hund_term=hund, IE=ie, IE_code=code,
                IE_NIST=ref["IE_eV"] if ref else None, err_pct=(ie / ref["IE_eV"] - 1) * 100 if ref else None)


def fmt_pa(d, cfg):
    sym = SYMBOL.get(d["Z"], f"Z={d['Z']}")
    lines = [f"### {sym} (Z={d['Z']}, N={d['N']}), model `{d['model']}`",
             f"configuration {cfg}; removed {int(d['n'])}{d['l']} (j={d['j']}, occupancy k={int(d['k'])}); "
             f"Za = Z-N+1 = {d['Za']:g}", "",
             "| quantity | value |", "|---|---|"]
    nu = d["nu_groups"]
    tau = d["tau_groups"]
    lines.append("| nu_g (same, in, core, df, out) | " + ", ".join(f"{nu[g]:g}" for g in GROUPS) + " |")
    lines.append("| tau_g (same, in, core, df, out) | " + ", ".join(f"{tau[g]:.4f}" for g in GROUPS) + " |")
    lines.append(f"| sum tau_g nu_g | {d['T_groups']:.5f} |")
    for c, v in d["class_devs"].items():
        lines.append(f"| + dtau_{c} x nu_{c} | {v['dtau']:+.4f} x {v['nu']:g} = {v['contrib']:+.5f} |")
    lines += [f"| T | {d['T']:.5f} |",
              f"| sigma1 (exact, Track A) | {d['sigma1']:.5f} |",
              f"| h = {d['h_rule']} | {d['h']:.5f} |",
              f"| kappa (fit; used = abs+0.05) | {d['kappa_fit']:.4f}; {d['kappa_used']:.4f} |",
              f"| denominator Za + kappa + abs(T)/h | {d['denominator']:.5f} |",
              f"| D_rem = T / denominator | {d['D_rem']:.5f} |",
              f"| total screening sigma1 + D_rem | {d['total_screening']:.5f} |",
              f"| Zeff = Z - sigma1 - D_rem | {d['Zeff']:.5f} (bound Za<=Zeff<=Z: {d['bound_ok']}) |",
              f"| Rydberg term Ry Zeff^2/n^2 | {d['rydberg_term']:.5f} eV |",
              f"| Dirac factor D_nj(Zeff) | {d['F_Dirac']:.8f} |"]
    if d["r_c"] is not None:
        lines.append(f"| rel. bracket 1 + r_c (Z alpha)^2/n (Zeff/Za - 1), r_{d['rel_class']} = {d['r_c']:.4f} "
                     f"| {d['rel_bracket']:.8f} |")
    else:
        lines.append("| rel. bracket | none in this model (= 1) |")
    lines += [f"| hydrogenic term after relativity | {d['hydrogenic_after_rel']:.5f} eV |",
              f"| Hund term Ry x_l K_l(k)/n^2 (x_l={d['x_l']:.4f}, K={d['K']:g}) | {d['hund_term']:+.5f} eV |",
              f"| reduced-mass factor mu(Z) | {d['mu']:.9f} |",
              f"| QED+FNS term (ns, n<=2 only) | {d['qed_fns_term']:.6f} eV |",
              f"| **IE** | **{d['IE']:.4f} eV** |"]
    if d["IE_NIST"] is not None:
        lines.append(f"| NIST | {d['IE_NIST']:.4f} eV (error {d['err_pct']:+.2f} %) |")
    return "\n".join(lines)


def fmt_pocket(d, cfg):
    sym = SYMBOL.get(d["Z"], f"Z={d['Z']}")
    nu, s = d["nu_groups"], d["s_groups"]
    lines = [f"### {sym} (Z={d['Z']}, N={d['N']}), model `pocket` (8 parameters)",
             f"configuration {cfg}; removed {int(d['n'])}{d['l']} (j={d['j']}); Za = {d['Za']:g}", "",
             "| quantity | value |", "|---|---|",
             "| nu_g (same, in, core, df, out) | " + ", ".join(f"{nu[g]:g}" for g in GROUPS) + " |",
             "| s_g (same, in, core, df, out) | " + ", ".join(f"{s[g]:.4f}" for g in GROUPS) + " |",
             f"| sum s_g nu_g | {d['screening_groups']:.5f} |",
             f"| t (N-1)/(Za + kappa), t={d['t']:.4f}, kappa used={d['kappa_used']:.4f} | {d['charge_term']:.5f} |",
             f"| Zeff | {d['Zeff']:.5f} |",
             f"| Rydberg term Ry Zeff^2/n^2 | {d['rydberg_term']:.5f} eV |",
             f"| Sommerfeld bracket | {d['sommerfeld_bracket']:.8f} |",
             f"| Hund term Ry x K/n^2 (x={d['x']:.4f}, K={d['K']:g}) | {d['hund_term']:+.5f} eV |",
             f"| **IE** | **{d['IE']:.4f} eV** |"]
    if d["IE_NIST"] is not None:
        lines.append(f"| NIST | {d['IE_NIST']:.4f} eV (error {d['err_pct']:+.2f} %) |")
    return "\n".join(lines)


def main(argv):
    pairs = [tuple(int(v) for v in s.split(":")) for s in argv] or DEFAULT
    nist = _nist()
    out, md = [], ["# Component audit of the Screened Rydberg formula", "",
                   "Generated by `tools/audit_components.py`. Every IE below is recomputed from the printed "
                   "components and asserted equal (1e-9 relative) to the production code.", ""]
    for Z, N in pairs:
        cfg = _cfg_str(gshm.initial_final(Z, N)[0]) if hasattr(gshm, "initial_final") else ""
        try:
            from configs import initial_final
            cfg = _cfg_str(initial_final(Z, N)[0])
        except Exception:
            pass
        for name in ("pa_hier_rel", "pa_bound9"):
            d = audit_pa(name, Z, N, nist)
            d["config"] = cfg
            out.append(d)
            md += [fmt_pa(d, cfg), ""]
        d = audit_pocket(Z, N, nist)
        d["config"] = cfg
        out.append(d)
        md += [fmt_pocket(d, cfg), ""]
    summ = ["## Summary (IE in eV)", "", "| ion | NIST | pa_hier_rel (33 p) | pa_bound9 (9 p) | pocket (8 p) |",
            "|---|---|---|---|---|"]
    for Z, N in pairs:
        r = {d["model"]: d for d in out if d["Z"] == Z and d["N"] == N}
        nst = r["pa_hier_rel"]["IE_NIST"]
        cell = lambda m: f"{r[m]['IE']:.3f} ({r[m]['err_pct']:+.2f} %)" if nst else f"{r[m]['IE']:.3f}"
        summ.append(f"| {SYMBOL.get(Z, Z)} Z={Z} N={N} | {nst if nst else '-'} | {cell('pa_hier_rel')} | "
                    f"{cell('pa_bound9')} | {cell('pocket')} |")
    md[4:4] = summ + [""]
    if not argv:
        with open(os.path.join(_ROOT, "results", "audit_components.json"), "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1)
        with open(os.path.join(_ROOT, "results", "audit_components.md"), "w", encoding="utf-8") as f:
            f.write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main(sys.argv[1:])
