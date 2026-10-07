r"""Hand recomputation of the Screened Rydberg formula from the paper text alone (checklist items 1-2).

    py -3.11 tools/hand_recompute.py              # six cases, worked examples, 5847-row precision study
    py -3.11 tools/hand_recompute.py --no-all     # skip the 5847-row study
    py -3.11 tools/hand_recompute.py --json F     # also dump all numbers to F (outside the project)

What it does
------------
For the neutral first IEs of O, Na, N, Ca, Fe, Pb (and Mg2+, the paper's second worked example) it
recomputes every component of three models from the paper's equations and printed parameters only:

  pa_hier_rel  final model, paper sec. 2.2 equation, parameters parsed verbatim from Table A1
  pa_bound9    9-parameter bounded variant: same equations without class deviations and without the
               relativistic bracket. Its parameters are NOT printed in the paper; full precision is read
               from results/pa_params.json (candidates.pa_bound9) and "printed" precision is simulated
               by rounding to 4 decimals, the precision of Table A1.
  pocket       8-parameter pocket formula, paper Appendix A; values parsed verbatim from the Appendix-A
               sentence, full precision from results/uni_params.json "pocket".

Inputs (and nothing else; nothing under models/ is imported or read):
  docs/paper_draft.md                     equations and printed parameter values (parsed, not retyped)
  docs/unified.md                         only to verify the quoted worked-example numbers (sec. 6)
  results/fp_zexp_rows_coefficients.csv   exact sigma1 per (Z, N) (paper sec. 4.6 points to this file)
  results/fp_zexp_coefficients.csv        cross-check of sigma1 and E1 for neutral configurations
  data/nist_ie.csv                        configurations (shells_expanded) and NIST IEs (for errors only)
  results/pa_params.json, results/uni_params.json   full-precision parameter values
  CODATA 2018 constants (stated below)
  results/audit_components.json           ONLY for the final component comparison with the code
  results/{pa_hier_rel,pa_bound9,uni_pocket}_predictions.csv   ONLY for the 5847-row comparison

Where the paper does not define something, the script says so and implements each candidate reading:
  * nu_g groups: G_LIT (literal sec. 2.2 wording), G_OUT (same, but same-n higher-l electrons of an s/p
    target go to "out"), G_DBL (literal, but an s target's same-subshell partner is also counted as a
    "same-n s" electron in "in").
  * nu_c classes (19 class deviations): the paper never defines them. They are implemented from the
    definitions in the code docstring of the class rule (Track B, quoted in CLASS_RULE below), which a
    reader of the paper does not have; the "paper-only" alternative is to omit them (variant NOCLS).
  * mu(Z): the paper gives no formula or masses. MU_A: mu = 1/(1 + m_e/(A m_u)) with A the mass number of
    the isotope (14N, 16O, 23Na, 24Mg, 40Ca, 56Fe, 208Pb; the isotopes listed for these Z in the project's
    Yerokhin-Shabaev table, read once by the reviewer); MU_1: mu = 1.
  * kappa: the equation shows kappa; Table A1 prints the fit value 1.7527 and "1.8027 as used".
  * QED+FNS: not computable from the paper (no table values). It is zero for all seven cases here
    (none removes an ns electron with n <= 2). For the 5847-row study the term is extracted, not computed.
"""
import argparse
import csv
import json
import math
import os
import re
import sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(ROOT, "docs", "paper_draft.md")
UNIFIED = os.path.join(ROOT, "docs", "unified.md")
NIST = os.path.join(ROOT, "data", "nist_ie.csv")
SIG_ROWS = os.path.join(ROOT, "results", "fp_zexp_rows_coefficients.csv")
SIG_TAB = os.path.join(ROOT, "results", "fp_zexp_coefficients.csv")
PA_JSON = os.path.join(ROOT, "results", "pa_params.json")
UNI_JSON = os.path.join(ROOT, "results", "uni_params.json")
AUDIT = os.path.join(ROOT, "results", "audit_components.json")
PRED = {"pa_hier_rel": os.path.join(ROOT, "results", "pa_hier_rel_predictions.csv"),
        "pa_bound9": os.path.join(ROOT, "results", "pa_bound9_predictions.csv"),
        "pocket": os.path.join(ROOT, "results", "uni_pocket_predictions.csv")}

# ----------------------------------------------------------------------------- constants (CODATA 2018)
RY_CODATA = 13.605693122994      # eV, Rydberg energy R_inf h c
ALPHA = 7.2973525693e-3          # fine-structure constant
ME_U = 5.48579909065e-4          # electron mass in u
ALPHA_137 = 1.0 / 137.036        # value used in the worked example of docs/unified.md sec. 6
ISOTOPE_A = {7: 14, 8: 16, 11: 23, 12: 24, 20: 40, 26: 56, 82: 208}

LL = "spdf"
GROUPS = ["same", "in", "core", "df", "out"]
REL = ["s", "p1", "p3", "d", "f"]
CASES = [(8, 8), (11, 11), (7, 7), (20, 20), (26, 26), (82, 82), (12, 10)]
SYM = {7: "N", 8: "O", 11: "Na", 12: "Mg", 20: "Ca", 26: "Fe", 82: "Pb"}
MODELS = ["pa_hier_rel", "pa_bound9", "pocket"]

# Class definitions as written in the docstring of the code's class rule (NOT in the paper):
CLASS_RULE = """
  same_<l>      other electrons of the same subshell
  sn_in_p       p target, s electrons of the same n
  n1_sp_<g>     s/p target, electrons of shell n-1 with l-group g in {sp, d, f}
  n2_sp_<g>     s/p target, electrons of shell n-2, g in {sp, df}
  d_near        d target, s,p electrons of shells n and n-1
  n1_d_d        d target, (n-1)d electrons
  n1_d_f        d target, (n-1)f electrons (4f for 5d)
  f_near        f target, s,p,d electrons of shells n and n-1
  n1_f_f        f target, (n-1)f electrons (4f for 5f)
  core          everything deeper (s/p target: n' <= n-3; d/f target: n' <= n-2), split into
                core_sp / core_df by target type for the deviations
  sn_out        same n, higher l
  out_<g>       electrons in shells with n' > n
"""
DEV19 = ["same_s", "same_p", "same_d", "same_f", "sn_in_p", "n1_sp_sp", "n1_sp_d", "n1_sp_f", "n2_sp_sp",
         "n2_sp_df", "d_near", "n1_d_d", "n1_d_f", "f_near", "n1_f_f", "out_d", "out_f", "core_sp", "core_df"]


def _num(s):
    return float(s.strip().replace("−", "-"))


def _flat(s):
    return " ".join(s.split())


# ----------------------------------------------------------------------------- paper parameters
def read_paper():
    """Parse every printed number the formula needs from docs/paper_draft.md (verbatim)."""
    txt = open(PAPER, encoding="utf-8").read()
    p = {}
    m = re.search(r"converted to eV via Ry = ([0-9.]+) eV", txt)
    p["Ry"] = float(m.group(1))
    m = re.search(r"\| τ_g \(same, in, core, df, out\) \| ([^|]+) \|", txt)
    p["tau"] = dict(zip(GROUPS, [_num(x) for x in m.group(1).split(",")]))
    m = re.search(r"\| κ \| ([0-9.]+) \(([0-9.]+) as used in the code, which applies \\\|κ\\\| \+ 0\.05\) \|", txt)
    p["kappa_fit"], p["kappa_used"] = float(m.group(1)), float(m.group(2))
    m = re.search(r"\| δτ_c \(19 classes\) \| ([^|]+) \|", txt)
    p["dtau"] = {k: _num(v) for k, v in re.findall(r"([a-z][a-z0-9_]*) ([−\-]?[0-9.]+)", m.group(1))}
    m = re.search(r"\| r_c \(s, p½, p3/2, d, f\) \| ([^|]+) \|", txt)
    p["r"] = dict(zip(REL, [_num(x) for x in m.group(1).split(",")]))
    m = re.search(r"\| x_l \(p, d, f\) \| ([^|]+) \|", txt)
    p["x"] = dict(zip("pdf", [_num(x) for x in m.group(1).split(",")]))
    flat = _flat(txt)
    m = re.search(r"The fitted values are s_same ([0-9.]+), s_in ([0-9.]+), s_core ([0-9.]+), s_df ([0-9.]+), "
                  r"s_out ([0-9.]+), t ([0-9.]+), κ ([0-9.]+) \(([0-9.]+) as used\) and x ([0-9.]+)\.", flat)
    v = [float(g) for g in m.groups()]
    p["pocket"] = dict(s=dict(zip(GROUPS, v[:5])), t=v[5], kappa_fit=v[6], kappa_used=v[7], x=v[8])
    assert len(p["dtau"]) == 19 and list(p["dtau"]) == DEV19, p["dtau"]
    return p


def read_full():
    """Full-precision values (results/pa_params.json, results/uni_params.json)."""
    d = json.load(open(PA_JSON, encoding="utf-8"))["candidates"]
    out = {}
    for name in ("pa_hier_rel", "pa_bound9"):
        P = dict(zip(d[name]["names"], d[name]["values"]))
        out[name] = dict(tau={g: P["tau_" + g] for g in GROUPS}, kappa_fit=P["kappa"],
                         kappa_used=abs(P["kappa"]) + 0.05,
                         dtau={c: P["dtau_" + c] for c in DEV19} if "dtau_same_s" in P else {},
                         r={c: P["rel_" + c] for c in REL} if "rel_s" in P else None,
                         x={c: P["x_" + c] for c in "pdf"})
    u = json.load(open(UNI_JSON, encoding="utf-8"))["pocket"]
    P = dict(zip(u["names"], u["values"]))
    out["pocket"] = dict(s={g: P["s_" + g] for g in GROUPS}, t=P["t"], kappa_fit=P["kappa"],
                         kappa_used=abs(P["kappa"]) + 0.05, x=P["x"])
    return out


def paper_params(pp):
    """Parameter sets exactly as printed (pa_bound9: not printed; simulated at Table A1 precision)."""
    return {"pa_hier_rel": dict(tau=pp["tau"], kappa_fit=pp["kappa_fit"], kappa_used=pp["kappa_used"],
                                dtau=pp["dtau"], r=pp["r"], x=pp["x"]),
            "pocket": pp["pocket"]}


def rounded(full, nd):
    """Round a full-precision parameter set to nd decimals (kappa_used = round(|kappa|+0.05))."""
    def r(v):
        return round(v, nd)
    if "s" in full:
        return dict(s={g: r(v) for g, v in full["s"].items()}, t=r(full["t"]), kappa_fit=r(full["kappa_fit"]),
                    kappa_used=r(full["kappa_used"]), x=r(full["x"]))
    return dict(tau={g: r(v) for g, v in full["tau"].items()}, kappa_fit=r(full["kappa_fit"]),
                kappa_used=r(full["kappa_used"]), dtau={c: r(v) for c, v in full["dtau"].items()},
                r=None if full["r"] is None else {c: r(v) for c, v in full["r"].items()},
                x={c: r(v) for c, v in full["x"].items()})


# ----------------------------------------------------------------------------- raw inputs
def parse_shells(s):
    out = {}
    for tok in s.split("."):
        tok = tok.strip()
        if not tok:
            continue
        i = 0
        while tok[i].isdigit():
            i += 1
        n, l, q = int(tok[:i]), LL.index(tok[i]), int(tok[i + 1:] or 1)
        out[(n, l)] = out.get((n, l), 0) + q
    return sorted((n, l, q) for (n, l), q in out.items() if q > 0)


def load_nist():
    recs = {}
    with open(NIST, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            Z, q = int(r["Z"]), int(r["ion_charge"])
            recs[(Z, Z - q)] = dict(shells=parse_shells(r["shells_expanded"]), IE=float(r["IE_eV"]),
                                    status=r["status"])
    return recs


def load_sigma():
    rows = {}
    with open(SIG_ROWS, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            nl = r["removed"]
            rows[(int(r["Z"]), int(r["N"]))] = dict(sigma1=float(r["sigma_abinitio"]), removed=(int(nl[:-1]),
                                                    LL.index(nl[-1])), j=r["j"], rearranged=r["rearranged"] == "1",
                                                    config=r["config"], slater=float(r["sigma_slater"]))
    tab = {}
    with open(SIG_TAB, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["which"] == "neutral":
                tab[int(r["N"])] = dict(sigma1=float(r["sigma_abinitio"]), E1=r["E1_exact"], dE1=float(r["dE1"]))
    return rows, tab


def removed_subshell(shells_N, shells_Nm1):
    """Subshell that loses an electron between the N and N-1 ground configurations. For rearranged ions
    (not covered by the paper) the subshell that lost the most electrons, ties -> outermost (n, then l)."""
    a = {(n, l): q for n, l, q in shells_N}
    b = {(n, l): q for n, l, q in shells_Nm1}
    diff = {k: a.get(k, 0) - b.get(k, 0) for k in set(a) | set(b)}
    lost = {k: d for k, d in diff.items() if d > 0}
    gained = {k: d for k, d in diff.items() if d < 0}
    rearr = bool(gained) or sum(lost.values()) != 1
    if not lost:
        return None, True
    return max(lost, key=lambda k: (lost[k], k[0], k[1])), rearr


# ----------------------------------------------------------------------------- paper building blocks
def jj_j(l, k, after=False):
    """Paper sec. 2.4: j = l - 1/2 while k <= 2l, otherwise l + 1/2 (k = occupancy before removal)."""
    kk = k - 1 if after else k
    if l == 0:
        return 0.5
    return l - 0.5 if kk <= 2 * l else l + 0.5


def rel_class(l, j):
    return ["s", "p1" if j == 0.5 else "p3", "d", "f"][l]


def pairs(l, k):
    """Parallel-spin pairs of l^k under Hund's first rule (maximum spin)."""
    up = min(k, 2 * l + 1)
    dn = k - up
    return up * (up - 1) // 2 + dn * (dn - 1) // 2


def K_l(l, k):
    """Paper sec. 2.2: K_l(k) = [P(k) - P(k-1)] - 2l(k-1)/(4l+1)."""
    if l == 0:
        return 0.0
    return (pairs(l, k) - pairs(l, k - 1)) - 2.0 * l * (k - 1) / (4 * l + 1)


def dirac_ratio(n, j, zeta, alpha=ALPHA):
    """Point-nucleus Dirac binding energy / Schroedinger energy for charge zeta (paper sec. 2.4).
    E_D = mc^2 [1 - (1 + (x/(n - k + sqrt(k^2 - x^2)))^2)^(-1/2)], k = j + 1/2, x = zeta alpha;
    E_S = mc^2 x^2/(2 n^2)."""
    x = zeta * alpha
    k = j + 0.5
    y = (x / (n - k + math.sqrt(k * k - x * x))) ** 2
    s = math.sqrt(1.0 + y)
    return y / (s * (1.0 + s)) * 2.0 * n * n / (x * x)


def nu_groups(shells, n, l, reading="G_LIT"):
    """Paper sec. 2.2 wording: same subshell; inner = same-n s and the (n-1) shell (s/p targets); deeper
    core (s/p); all n' <= n (d/f targets); outer n' > n. Returns (counts, electrons in no group)."""
    g = dict.fromkeys(GROUPS, 0)
    nogroup = 0
    for n2, l2, q in shells:
        same = (n2, l2) == (n, l)
        c = q - 1 if same else q
        if c <= 0:
            continue
        if same:
            g["same"] += c
            if reading == "G_DBL" and l == 0:
                g["in"] += c          # "same-n s" read as including the target's own subshell partner
        elif n2 > n:
            g["out"] += c
        elif l <= 1:
            if n2 == n and l2 == 0:
                g["in"] += c          # same-n s electrons of a p target
            elif n2 == n:             # same n, higher l (e.g. 3d on a 3p target): not in the wording
                if reading == "G_OUT":
                    g["out"] += c
                else:
                    nogroup += c
            elif n2 == n - 1:
                g["in"] += c
            else:
                g["core"] += c
        else:
            if n2 == n and l2 > l and reading == "G_OUT":
                g["out"] += c         # code puts same-n higher-l electrons of d/f targets in "out"
            else:
                g["df"] += c          # literal: all n' <= n
    return g, nogroup


def nu_classes(shells, n, l):
    """19 deviation classes + sn_out + out_sp, from CLASS_RULE (not defined in the paper)."""
    c = defaultdict(float)
    grp = lambda l2: "sp" if l2 <= 1 else LL[l2]  # noqa: E731
    for n2, l2, q in shells:
        same = (n2, l2) == (n, l)
        cnt = q - 1 if same else q
        if cnt <= 0:
            continue
        dn = n - n2
        if same:
            k = "same_" + LL[l]
        elif dn == 0 and l2 > l:
            k = "sn_out"
        elif dn < 0:
            k = "out_" + grp(l)
        elif l <= 1:
            k = ("sn_in_p" if dn == 0 else "n1_sp_" + grp(l2) if dn == 1 else
                 "n2_sp_" + ("sp" if l2 <= 1 else "df") if dn == 2 else "core_sp")
        elif l == 2:
            k = "d_near" if (dn <= 1 and l2 <= 1) else ("n1_d_" + LL[l2] if dn == 1 else "core_df")
        else:
            k = "f_near" if (dn <= 1 and l2 <= 2) else ("n1_f_f" if dn == 1 else "core_df")
        c[k] += cnt
    return dict(c)


CLASS_GROUP = {"same_s": "same", "same_p": "same", "same_d": "same", "same_f": "same", "sn_in_p": "in",
               "n1_sp_sp": "in", "n1_sp_d": "in", "n1_sp_f": "in", "n2_sp_sp": "core", "n2_sp_df": "core",
               "core_sp": "core", "d_near": "df", "n1_d_d": "df", "n1_d_f": "df", "f_near": "df", "n1_f_f": "df",
               "core_df": "df", "out_d": "out", "out_f": "out", "out_sp": "out", "sn_out": "out"}


# ----------------------------------------------------------------------------- one row, all components
def features(Z, N, nist, sig):
    sh = nist[(Z, N)]["shells"]
    if N == 1:
        rem, rearr = (sh[0][0], sh[0][1]), False
    elif (Z, N - 1) not in nist:       # (N-1)-electron ion not tabulated (not covered by the paper):
        rem, rearr = max(((a, b) for a, b, _ in sh), key=lambda t: t), False   # outermost (n, then l)
    else:
        rem, rearr = removed_subshell(sh, nist[(Z, N - 1)]["shells"])
    n, l = rem
    k = dict(((a, b), q) for a, b, q in sh)[rem]
    return dict(Z=Z, N=N, shells=sh, n=n, l=l, k=k, rearr=rearr, Za=Z - N + 1, sigma1=sig[(Z, N)]["sigma1"])


def pa_components(f, P, model, opt):
    """Paper sec. 2.2/2.4 for pa_hier_rel (model='pa_hier_rel') or pa_bound9 (no classes, no bracket)."""
    Z, N, n, l, k, Za = f["Z"], f["N"], f["n"], f["l"], f["k"], f["Za"]
    s1 = round(f["sigma1"], opt.get("sig_nd", 99)) if opt.get("sig_nd") else f["sigma1"]
    Ry, al = opt.get("Ry", RY_CODATA), opt.get("alpha", ALPHA)
    nu, nogroup = nu_groups(f["shells"], n, l, opt.get("groups", "G_OUT"))
    T_groups = sum(P["tau"][g] * nu[g] for g in GROUPS)
    devs = {}
    if model == "pa_hier_rel" and not opt.get("no_classes"):
        cls = nu_classes(f["shells"], n, l)
        for c in DEV19:
            if cls.get(c):
                devs[c] = dict(nu=cls[c], dtau=P["dtau"][c], contrib=cls[c] * P["dtau"][c])
    T = T_groups + sum(d["contrib"] for d in devs.values())
    h = (N - 1) - s1 if T >= 0 else s1
    kap = P["kappa_fit"] if opt.get("kappa_fit") else P["kappa_used"]
    if h == 0.0:                       # N = 1: T = 0 and h = 0 -> 0/0 in the paper; limit D = 0
        denom, D = float("nan"), 0.0
    else:
        denom = Za + kap + abs(T) / h
        D = T / denom
    Ze = Z - s1 - D
    ryd = Ry * Ze ** 2 / n ** 2
    j = jj_j(l, k, after=opt.get("j_after", False))
    F = 1.0 if opt.get("no_dirac") else dirac_ratio(n, j, Ze, al)
    rc = rel_class(l, j)
    if model == "pa_hier_rel":
        R = 1.0 + P["r"][rc] * (Z * al) ** 2 / n * (Ze / Za - 1.0)
        r_c = P["r"][rc]
    else:
        R, r_c = 1.0, None
    hyd = ryd * F * R
    K = K_l(l, k)
    x_l = 0.0 if l == 0 else P["x"]["pdf"[l - 1]]
    hund = Ry * x_l * K / n ** 2
    mu = 1.0 if opt.get("mu1") else 1.0 / (1.0 + ME_U / ISOTOPE_A[Z]) if Z in ISOTOPE_A else 1.0
    qed = 0.0 if not (l == 0 and n <= 2) else float("nan")
    IE = mu * (hyd + hund) - (0.0 if qed != qed else qed)
    return dict(model=model, nu_groups=nu, nogroup=nogroup, T_groups=T_groups, class_devs=devs, T=T, sigma1=s1,
                h=h, kappa_used=kap, denominator=denom, D_rem=D, total_screening=s1 + D, Zeff=Ze,
                rydberg_term=ryd, j=j, F_Dirac=F, rel_class=rc, r_c=r_c, rel_bracket=R, hydrogenic_after_rel=hyd,
                x_l=x_l, K=K, hund_term=hund, mu=mu, qed_fns_term=qed, IE=IE)


def pocket_components(f, P, opt):
    """Paper Appendix A."""
    Z, N, n, l, k, Za = f["Z"], f["N"], f["n"], f["l"], f["k"], f["Za"]
    Ry, al = opt.get("Ry", RY_CODATA), opt.get("alpha", ALPHA)
    nu, nogroup = nu_groups(f["shells"], n, l, opt.get("groups", "G_OUT"))
    S = sum(P["s"][g] * nu[g] for g in GROUPS)
    kap = P["kappa_fit"] if opt.get("kappa_fit") else P["kappa_used"]
    Q = P["t"] * (N - 1) / (Z - N + 1 + kap)
    Ze = Z - S - Q
    j = jj_j(l, k, after=opt.get("j_after", False))
    ryd = Ry * Ze ** 2 / n ** 2
    B = 1.0 + (Ze * al) ** 2 / n ** 2 * (n / (j + 0.5) - 0.75)
    K = K_l(l, k)
    hund = Ry * P["x"] * K / n ** 2
    return dict(model="pocket", nu_groups=nu, nogroup=nogroup, screening_groups=S, kappa_used=kap, charge_term=Q,
                total_screening=S + Q, Zeff=Ze, j=j, rydberg_term=ryd, sommerfeld_bracket=B, K=K,
                hund_term=hund, IE=ryd * B + hund)


def compute(model, f, P, opt):
    return pocket_components(f, P, opt) if model == "pocket" else pa_components(f, P, model, opt)


# ----------------------------------------------------------------------------- comparison with the code
PA_KEYS = ["T_groups", "T", "sigma1", "h", "kappa_used", "denominator", "D_rem", "total_screening", "Zeff",
           "rydberg_term", "F_Dirac", "rel_bracket", "hydrogenic_after_rel", "x_l", "K", "hund_term", "mu",
           "qed_fns_term", "IE"]
PK_KEYS = ["screening_groups", "kappa_used", "charge_term", "total_screening", "Zeff", "rydberg_term",
           "sommerfeld_bracket", "K", "hund_term", "IE"]


def fmt(v, nd=6):
    if v is None:
        return "-"
    if isinstance(v, float) and v != v:
        return "n/a"
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def comp_table(case, hp, hf, code):
    """Markdown rows: quantity | paper-as-printed | code | diff | paper full precision | code | diff."""
    keys = PK_KEYS if hp["model"] == "pocket" else PA_KEYS
    lines = ["| quantity | as printed (P4) | code | P4 - code | full precision (FP) | FP - code |",
             "|---|---|---|---|---|---|"]
    nu_p = ", ".join(f"{hp['nu_groups'][g]:g}" for g in GROUPS)
    nu_c = ", ".join(f"{code['nu_groups'][g]:g}" for g in GROUPS)
    lines.append(f"| nu_g (same, in, core, df, out) | {nu_p} | {nu_c} | {'0' if nu_p == nu_c else 'DIFF'} | "
                 f"{nu_p} | {'0' if nu_p == nu_c else 'DIFF'} |")
    if hp["model"] != "pocket":
        cp = {c: d["nu"] for c, d in hp["class_devs"].items()}
        cc = {c: d["nu"] for c, d in code["class_devs"].items()}
        s_p = ", ".join(f"{c} {v:g}" for c, v in cp.items()) or "none"
        s_c = ", ".join(f"{c} {v:g}" for c, v in cc.items()) or "none"
        lines.append(f"| nu_c (class counts) | {s_p} | {s_c} | {'0' if cp == cc else 'DIFF'} | {s_p} | "
                     f"{'0' if cp == cc else 'DIFF'} |")
        lines.append(f"| j, rel. class | {hp['j']}, {hp['rel_class']} | {code['j']}, {code['rel_class']} | "
                     f"{'0' if (hp['j'], hp['rel_class']) == (code['j'], code['rel_class']) else 'DIFF'} | | |")
    for k in keys:
        a, b, c = hp.get(k), hf.get(k), code.get(k)
        if c is None and k == "rel_bracket":
            c = 1.0
        da = a - c if isinstance(a, float) and isinstance(c, (int, float)) else None
        db = b - c if isinstance(b, float) and isinstance(c, (int, float)) else None
        lines.append(f"| {k} | {fmt(a, 8)} | {fmt(float(c) if c is not None else None, 8)} | "
                     f"{'' if da is None else f'{da:+.2e}'} | {fmt(b, 8)} | {'' if db is None else f'{db:+.2e}'} |")
    return "\n".join(lines)


# ----------------------------------------------------------------------------- worked examples
def _close(x, printed, nd):
    return abs(x - printed) <= 0.5 * 10 ** (-nd) + 1e-12


def worked_examples(res, nist):
    """Every printed number of paper sec. 4.7 and docs/unified.md sec. 6, checked two ways:
    'arith' = the printed arithmetic (re-evaluated from the printed operands) gives the printed result;
    'code' = the code value (audit JSON; Mg2+ = full-precision recomputation, IE checked against the
    production predictions file) rounded to the printed precision equals the printed number."""
    paper = _flat(open(PAPER, encoding="utf-8").read())
    uni = _flat(open(UNIFIED, encoding="utf-8").read())
    O, Op = res[(8, 8)]["code"]["pa_hier_rel"], res[(8, 8)]["code"]["pocket"]
    M, Mp = res[(12, 10)]["FP"]["pa_hier_rel"], res[(12, 10)]["FP"]["pocket"]
    nO, nM = nist[(8, 8)]["IE"], nist[(12, 10)]["IE"]
    a137 = ALPHA_137
    rows = [
        # (doc, quote, label, printed, decimals, arithmetic expression or None, code value)
        ("paper 4.7", "σ₁ = 5.23348 and T = 3.3496", "O sigma1", 5.23348, 5, None, O["sigma1"]),
        ("paper 4.7", "σ₁ = 5.23348 and T = 3.3496", "O T", 3.3496, 4, None, O["T"]),
        ("paper 4.7", "h = 1.76652, so D = 0.71285", "O h", 1.76652, 5, None, O["h"]),
        ("paper 4.7", "h = 1.76652, so D = 0.71285", "O D", 0.71285, 5, None, O["D_rem"]),
        ("paper 4.7", "Z_eff = 2.05367, which lies", "O Zeff", 2.05367, 5, None, O["Zeff"]),
        ("paper 4.7", "The hydrogenic term is 14.3457 eV", "O Ry Zeff^2/n^2", 14.3457, 4, None, O["rydberg_term"]),
        ("paper 4.7", "or 14.3162 eV after relativity", "O after relativity", 14.3162, 4, None,
         O["hydrogenic_after_rel"]),
        ("paper 4.7", "The Hund term is −0.8365 eV", "O Hund", -0.8365, 4, None, O["hund_term"]),
        ("paper 4.7", "IE = 13.479 eV, an error of −1.02 %", "O IE", 13.479, 3, None, O["IE"]),
        ("paper 4.7", "IE = 13.479 eV, an error of −1.02 %", "O error %", -1.02, 2, None, (O["IE"] / nO - 1) * 100),
        ("paper 4.7", "(1s²2s²2p⁴ → 2p³, NIST 13.618 eV)", "O NIST", 13.618, 3, None, nO),
        ("paper 4.7", "pocket calculation gives 9.70 eV (−28.8 %)", "O pocket IE", 9.70, 2, None, Op["IE"]),
        ("paper 4.7", "pocket calculation gives 9.70 eV (−28.8 %)", "O pocket error %", -28.8, 1, None,
         (Op["IE"] / nO - 1) * 100),
        ("paper 4.7", "For Mg²⁺ (NIST 80.144 eV)", "Mg2+ NIST", 80.144, 3, None, nM),
        ("paper 4.7", "the formula gives 78.63 eV (−1.88 %)", "Mg2+ IE", 78.63, 2, None, M["IE"]),
        ("paper 4.7", "the formula gives 78.63 eV (−1.88 %)", "Mg2+ error %", -1.88, 2, None, (M["IE"] / nM - 1) * 100),
        # docs/unified.md sec. 6.1 (oxygen)
        ("unified 6.1", "4 × 1.308370 = **5.23348**", "O sigma1 = 4 x 1.308370", 5.23348, 5, "4*1.308370",
         O["sigma1"]),
        ("unified 6.1", "Slater's rules give 3.45", "O Slater sigma", 3.45, 2, None, 3.45),
        ("unified 6.1", "= 4.3300 − 0.9804 = **3.3496**", "O sum tau_g nu_g", 4.3300, 4, "0.3928*3+0.7879*4",
         O["T_groups"]),
        ("unified 6.1", "= 4.3300 − 0.9804 = **3.3496**", "O class sum", -0.9804, 4, "3*0.0098+2*(-0.1133)+2*(-0.3916)",
         O["T"] - O["T_groups"]),
        ("unified 6.1", "= 4.3300 − 0.9804 = **3.3496**", "O T", 3.3496, 4, "4.3300-0.9804", O["T"]),
        ("unified 6.1", "h = N − 1 − σ₁ = 7 − 5.23348 = 1.76652", "O h", 1.76652, 5, "7-5.23348", O["h"]),
        ("unified 6.1", "u = 3.3496 / [1.76652·(1 + 1.8027)] = 0.67654", "O u", 0.67654, 5,
         "3.3496/(1.76652*(1+1.8027))", O["T"] / (O["h"] * (1 + O["kappa_used"]))),
        ("unified 6.1", "D = h·u/(1+u) = **0.71285**", "O D", 0.71285, 5, "1.76652*0.67654/(1+0.67654)", O["D_rem"]),
        ("unified 6.1", "Z_eff = 8 − 5.23348 − 0.71285 = **2.05367**", "O Zeff", 2.05367, 5, "8-5.23348-0.71285",
         O["Zeff"]),
        ("unified 6.1", "13.6057 × 4.21756/4 = 14.3457 eV", "O Zeff^2", 4.21756, 5, "2.05367**2", O["Zeff"] ** 2),
        ("unified 6.1", "13.6057 × 4.21756/4 = 14.3457 eV", "O Ry Zeff^2/4", 14.3457, 4, "13.6057*4.21756/4",
         O["rydberg_term"]),
        ("unified 6.1", "D_{2,3/2} = 1.000014", "O Dirac factor", 1.000014, 6, None, O["F_Dirac"]),
        ("unified 6.1", "(2.05367 − 1)/2 = 0.997929", "O rel. bracket", 0.997929, 6,
         "1-1.1533*(8/137.036)**2*(2.05367-1)/2", O["rel_bracket"]),
        ("unified 6.1", "Together: 14.3162 eV", "O after relativity", 14.3162, 4, "14.3457*1.000014*0.997929",
         O["hydrogenic_after_rel"]),
        ("unified 6.1", "Ry · 0.2049 · (−1.2)/4 = −0.8365 eV", "O Hund", -0.8365, 4, "13.6057*0.2049*(-1.2)/4",
         O["hund_term"]),
        ("unified 6.1", "(14.3162 − 0.8365) × 0.999966 = **13.479 eV**", "O mu", 0.999966, 6, None, O["mu"]),
        ("unified 6.1", "(14.3162 − 0.8365) × 0.999966 = **13.479 eV**", "O IE", 13.479, 3,
         "(14.3162-0.8365)*0.999966", O["IE"]),
        ("unified 6.1", "an error of **−1.02 %**", "O error %", -1.02, 2, "((14.3162-0.8365)*0.999966/13.618-1)*100",
         (O["IE"] / nO - 1) * 100),
        ("unified 6.1", "= 8 − 5.2553 − 0.8843 = **1.8605**", "O pocket sum s_g nu_g", 5.2553, 4, "0.7499*3+0.7514*4",
         Op["screening_groups"]),
        ("unified 6.1", "= 8 − 5.2553 − 0.8843 = **1.8605**", "O pocket charge term", 0.8843, 4, "1.3725*7/(1+9.865)",
         Op["charge_term"]),
        ("unified 6.1", "= 8 − 5.2553 − 0.8843 = **1.8605**", "O pocket Zeff", 1.8605, 4, "8-5.2553-0.8843",
         Op["Zeff"]),
        ("unified 6.1", "(1 + 1.2·10⁻⁵)", "O pocket bracket - 1 (x1e5)", 1.2, 1, "((1.8605/137.036/2)**2*(2/2-0.75))*1e5",
         (Op["sommerfeld_bracket"] - 1) * 1e5),
        ("unified 6.1", "= 11.775 − 2.074 = **9.70 eV**", "O pocket Rydberg x bracket", 11.775, 3,
         "13.6057*1.8605**2/4*(1+1.2e-5)", Op["rydberg_term"] * Op["sommerfeld_bracket"]),
        ("unified 6.1", "= 11.775 − 2.074 = **9.70 eV**", "O pocket Hund", -2.074, 3, "13.6057*0.508*(-1.2)/4",
         Op["hund_term"]),
        ("unified 6.1", "= 11.775 − 2.074 = **9.70 eV**", "O pocket IE", 9.70, 2, "11.775-2.074", Op["IE"]),
        ("unified 6.1", "**9.70 eV** (−28.8 %)", "O pocket error %", -28.8, 1, "(9.70/13.618-1)*100",
         (Op["IE"] / nO - 1) * 100),
        # docs/unified.md sec. 6.2 (Mg2+)
        ("unified 6.2", "σ₁ = **6.54598**", "Mg2+ sigma1", 6.54598, 5, None, M["sigma1"]),
        ("unified 6.2", "= 5.1156 − 0.9607 = **4.1548**", "Mg2+ sum tau_g nu_g", 5.1156, 4, "0.3928*5+0.7879*4",
         M["T_groups"]),
        ("unified 6.2", "= 5.1156 − 0.9607 = **4.1548**", "Mg2+ class sum", -0.9607, 4, "5*0.0098-2*0.1133-2*0.3916",
         M["T"] - M["T_groups"]),
        ("unified 6.2", "= 5.1156 − 0.9607 = **4.1548**", "Mg2+ T", 4.1548, 4, "5.1156-0.9607", M["T"]),
        ("unified 6.2", "h = 9 − 6.54598 = 2.45402", "Mg2+ h", 2.45402, 5, "9-6.54598", M["h"]),
        ("unified 6.2", "u = 4.1548/[2.45402·(3 + 1.8027)] = 0.35252", "Mg2+ u", 0.35252, 5,
         "4.1548/(2.45402*(3+1.8027))", M["T"] / (M["h"] * (3 + M["kappa_used"]))),
        ("unified 6.2", "D = 2.45402 · 0.35252/1.35252 = **0.63962**", "Mg2+ D", 0.63962, 5,
         "2.45402*0.35252/1.35252", M["D_rem"]),
        ("unified 6.2", "Z_eff = 12 − 6.54598 − 0.63962 = **4.81440**", "Mg2+ Zeff", 4.81440, 5,
         "12-6.54598-0.63962", M["Zeff"]),
        ("unified 6.2", "13.6057 · 23.1784/4 = 78.840 eV", "Mg2+ Zeff^2", 23.1784, 4, "4.81440**2", M["Zeff"] ** 2),
        ("unified 6.2", "13.6057 · 23.1784/4 = 78.840 eV", "Mg2+ Ry Zeff^2/4", 78.840, 3, "13.6057*23.1784/4",
         M["rydberg_term"]),
        ("unified 6.2", "The Dirac factor is 1.000077", "Mg2+ Dirac factor", 1.000077, 6, None, M["F_Dirac"]),
        ("unified 6.2", "(4.8144/3 − 1)/2 = 0.997326", "Mg2+ rel. bracket", 0.997326, 6,
         "1-1.1533*(12/137.036)**2*(4.8144/3-1)/2", M["rel_bracket"]),
        ("unified 6.2", "× 0.999977 = **78.63 eV**", "Mg2+ mu", 0.999977, 6, None, res[(12, 10)]["mu_code"]),
        ("unified 6.2", "× 0.999977 = **78.63 eV**", "Mg2+ IE", 78.63, 2, "78.840*1.000077*0.997326*0.999977",
         res[(12, 10)]["pred"]["pa_hier_rel"]),
        ("unified 6.2", "an error of **−1.88 %**", "Mg2+ error %", -1.88, 2,
         "(78.840*1.000077*0.997326*0.999977/80.144-1)*100", (res[(12, 10)]["pred"]["pa_hier_rel"] / nM - 1) * 100),
        ("unified 6.2", "= 12 − 6.7551 − 0.9602 = **4.2847**", "Mg2+ pocket sum s_g nu_g", 6.7551, 4,
         "0.7499*5+0.7514*4", Mp["screening_groups"]),
        ("unified 6.2", "= 12 − 6.7551 − 0.9602 = **4.2847**", "Mg2+ pocket charge term", 0.9602, 4,
         "1.3725*9/(3+9.865)", Mp["charge_term"]),
        ("unified 6.2", "= 12 − 6.7551 − 0.9602 = **4.2847**", "Mg2+ pocket Zeff", 4.2847, 4, "12-6.7551-0.9602",
         Mp["Zeff"]),
        ("unified 6.2", "4.2847²/4 · 1.00006 = **62.46 eV** (−22.1 %)", "Mg2+ pocket bracket", 1.00006, 5,
         "1+(4.2847/137.036/2)**2*(2/2-0.75)", Mp["sommerfeld_bracket"]),
        ("unified 6.2", "4.2847²/4 · 1.00006 = **62.46 eV** (−22.1 %)", "Mg2+ pocket IE", 62.46, 2,
         "13.6057*4.2847**2/4*1.00006", res[(12, 10)]["pred"]["pocket"]),
        ("unified 6.2", "4.2847²/4 · 1.00006 = **62.46 eV** (−22.1 %)", "Mg2+ pocket error %", -22.1, 1,
         "(62.46/80.144-1)*100", (res[(12, 10)]["pred"]["pocket"] / nM - 1) * 100),
        ("unified 6.2", "`ionization_energy(\"O\")` = 13.4793", "O IE (API)", 13.4793, 4, None, O["IE"]),
        ("unified 6.2", "`ionization_energy(12, 10)` = 78.6333", "Mg2+ IE (API)", 78.6333, 4, None,
         res[(12, 10)]["pred"]["pa_hier_rel"]),
        ("unified 6.2", "gives 9.7013 and 62.4557", "O pocket IE (API)", 9.7013, 4, None, Op["IE"]),
        ("unified 6.2", "gives 9.7013 and 62.4557", "Mg2+ pocket IE (API)", 62.4557, 4, None,
         res[(12, 10)]["pred"]["pocket"]),
        ("unified 6", "μ ≈ 0.99997", "mu (O, Mg)", 0.99997, 5, None, O["mu"]),
    ]
    out = []
    for doc, quote, label, printed, nd, expr, code in rows:
        src = paper if doc.startswith("paper") else uni
        found = _flat(quote) in src
        ar = None if expr is None else eval(expr, {"__builtins__": {}}, {"a137": a137})  # noqa: S307
        out.append(dict(doc=doc, quote=quote, found=found, label=label, printed=printed, nd=nd, expr=expr,
                        arith=ar, arith_ok=None if ar is None else _close(ar, printed, nd), code=code,
                        code_ok=_close(code, printed, nd)))
    return out


# ----------------------------------------------------------------------------- 5847-row study
def load_pred(path):
    d = {}
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            d[(int(r["Z"]), int(r["N"]))] = float(r["IE_pred_eV"])
    return d


def _median(v):
    v = sorted(v)
    return v[len(v) // 2] if len(v) % 2 else 0.5 * (v[len(v) // 2 - 1] + v[len(v) // 2])


def all_rows(nist, sig, full, pp):
    """Every NIST row: (i) does the paper-only recomputation (full precision) reproduce the code's
    prediction? (ii) how much do printed precision and each reading change the IE and the MAPE?
    The paper-only value B omits mu and the QED/FNS term (mu=1, QED=0). mu_Z = code/B is then a
    per-element constant on all rows without 1s/2s removal; a row is 'reproduced' when its ratio equals
    the median of its element to 1e-8. For the variants, the relative change of B is applied to the code
    prediction, so the MAPE shift includes the code's own mu and QED/FNS."""
    keys = sorted(nist)
    feats = {k: features(*k, nist, sig) for k in keys}
    to_j = {"1/2": 0.5, "3/2": 1.5, "5/2": 2.5, "7/2": 3.5}
    gap = {k for k in keys if k[1] > 1 and (k[0], k[1] - 1) not in nist}
    out = dict(n_rows=len(keys),
               n_table_gaps=sum(1 for Z in range(1, 111) for N in range(1, Z + 1) if (Z, N) not in nist),
               n_rows_Nm1_missing=len(gap),
               n_rearranged_diff_rule=sum(1 for k in keys if feats[k]["rearr"]),
               n_rearranged_csv=sum(1 for k in keys if sig[k]["rearranged"]),
               n_removed_differs_from_csv=sum(1 for k in keys if (feats[k]["n"], feats[k]["l"]) != sig[k]["removed"]),
               n_j_differs_from_csv=sum(1 for k in keys if jj_j(feats[k]["l"], feats[k]["k"]) != to_j[sig[k]["j"]]),
               n_sigma_nan=sum(1 for k in keys if sig[k]["sigma1"] != sig[k]["sigma1"]),
               sigma_outside_rows=[k for k in keys if sig[k]["sigma1"] == sig[k]["sigma1"] and
                                   not (0 <= sig[k]["sigma1"] <= k[1] - 1 + 1e-12)])
    out["n_sigma_outside_0_Nm1"] = len(out["sigma_outside_rows"])
    qrows = {k for k in keys if feats[k]["l"] == 0 and feats[k]["n"] <= 2}
    out["n_qed_rows"] = len(qrows)
    out["n_nogroup_rows_literal"] = sum(1 for k in keys if nu_groups(feats[k]["shells"], feats[k]["n"],
                                                                     feats[k]["l"], "G_LIT")[1] > 0)
    out["n_dfout_rows"] = sum(1 for k in keys if feats[k]["l"] >= 2 and any(
        n2 == feats[k]["n"] and l2 > feats[k]["l"] for n2, l2, _ in feats[k]["shells"]))
    pred = {m: load_pred(p) for m, p in PRED.items()}
    nistIE = {k: nist[k]["IE"] for k in keys}
    base = dict(mu1=True, groups="G_LIT")
    res = {}
    for m in MODELS:
        P_full = full[m]
        variants = {"FP": (P_full, {}), "printed": (pp.get(m), {}), "4dp": (rounded(P_full, 4), {}),
                    "5dp": (rounded(P_full, 5), {}), "6dp": (rounded(P_full, 6), {}),
                    "kappa_fit": (P_full, dict(kappa_fit=True)), "G_DBL": (P_full, dict(groups="G_DBL")),
                    "j_after": (P_full, dict(j_after=True))}
        if m == "pa_hier_rel":
            variants["no_classes"] = (P_full, dict(no_classes=True))
        if m != "pocket":
            variants["sig4"] = (P_full, dict(sig_nd=4))
        variants = {v: x for v, x in variants.items() if x[0] is not None}
        comps = {k: compute(m, feats[k], P_full, base) for k in keys}
        B = {"FP": {k: comps[k]["IE"] for k in keys}}
        for v, (P, o) in variants.items():
            if v != "FP":
                B[v] = {k: compute(m, feats[k], P, dict(base, **o))["IE"] for k in keys}
        code = pred[m]
        r = {}
        ok = {k for k in keys if B["FP"][k] == B["FP"][k]}
        mu_z = {}
        if m == "pocket":
            r["FP_vs_code_max_rel"] = max(abs(B["FP"][k] / code[k] - 1) for k in keys)
            notrep = sorted(k for k in keys if abs(B["FP"][k] / code[k] - 1) > 1e-8)
        else:
            r["n_h_zero"] = sum(1 for k in ok if comps[k]["h"] == 0)
            r["n_T_negative"] = sum(1 for k in ok if comps[k]["T"] < 0)
            ratio = defaultdict(list)
            for k in ok - qrows:
                ratio[k[0]].append(code[k] / B["FP"][k])
            mu_z = {z: _median(v) for z, v in ratio.items()}
            notrep = sorted(k for k in ok - qrows if abs(code[k] / B["FP"][k] / mu_z[k[0]] - 1) > 1e-8)
            r["n_Z_with_mu"] = len(mu_z)
            r["mu_Z_range"] = (min(mu_z.values()), max(mu_z.values()))
            good = (ok - qrows) - set(notrep)
            r["reproduced_max_dev"] = max(abs(code[k] / B["FP"][k] / mu_z[k[0]] - 1) for k in good)
            q = sorted(((mu_z[k[0]] * B["FP"][k] - code[k]) / code[k], mu_z[k[0]] * B["FP"][k] - code[k], k)
                       for k in qrows & ok if k[0] in mu_z)
            r["qed_rows_extracted"] = len(q)
            r["qed_rows_not_extractable_Zle4"] = len([k for k in qrows if k[0] not in mu_z])
            qm = max(q, key=lambda x: abs(x[0]))
            qa = max(q, key=lambda x: abs(x[1]))
            r["qed_max_rel_pct"], r["qed_max_rel_row"] = 100 * qm[0], qm[2]
            r["qed_max_abs_eV"], r["qed_max_abs_row"] = qa[1], qa[2]
            r["qed_median_rel_pct"] = 100 * _median([abs(x[0]) for x in q])
            r["qed_n_above_1e-4_rel"] = sum(1 for x in q if abs(x[0]) > 1e-4)
        r["n_not_reproduced"] = len(notrep) + len(set(keys) - ok)
        r["not_reproduced"] = [dict(row=k, csv_rearranged=sig[k]["rearranged"], diff_rule_rearranged=feats[k]["rearr"],
                                    Nm1_missing=k in gap,
                                    ratio_dev=code[k] / B["FP"][k] / mu_z.get(k[0], 1.0) - 1) for k in notrep]
        r["nan_rows"] = sorted(set(keys) - ok)
        stat = sorted(ok - set(notrep))
        r["n_stat_rows"] = len(stat)

        def mape(ie):
            return 100 * sum(abs(ie[k] / nistIE[k] - 1) for k in keys) / len(keys)

        r["MAPE_code"] = mape(code)
        sset = set(stat)
        for v in B:
            if v == "FP":
                continue
            rel = {k: (B[v][k] / B["FP"][k] - 1 if k in sset else 0.0) for k in keys}
            ie_v = {k: code[k] * (1 + rel[k]) for k in keys}
            kmax = max(stat, key=lambda k: abs(rel[k]))
            kabs = max(stat, key=lambda k: abs(code[k] * rel[k]))
            r[v] = dict(max_rel_pct=100 * abs(rel[kmax]), row_max_rel=kmax, max_abs_eV=abs(code[kabs] * rel[kabs]),
                        row_max_abs=kabs, n_rel_gt_1e6=sum(1 for k in stat if abs(rel[k]) > 1e-6),
                        n_rel_gt_1e4=sum(1 for k in stat if abs(rel[k]) > 1e-4),
                        n_rel_gt_1e3=sum(1 for k in stat if abs(rel[k]) > 1e-3),
                        dMAPE=mape(ie_v) - r["MAPE_code"])
        res[m] = r
    out["models"] = res
    return out


def all_rows_md(a):
    L = []
    L.append(f"- rows: {a['n_rows']} (the Z = 1-110 grid has {a['n_table_gaps']} (Z, N) pairs missing from "
             f"data/nist_ie.csv; {a['n_rows_Nm1_missing']} rows have no tabulated (N-1)-electron ion)")
    L.append(f"- rearranged rows: {a['n_rearranged_csv']} by the flag in fp_zexp_rows_coefficients.csv, "
             f"{a['n_rearranged_diff_rule']} when the (N-1) ion is tabulated and the configurations really differ by "
             f"more than one electron")
    L.append(f"- removed subshell (diff of ground configurations; outermost if the N-1 ion is missing) differs from "
             f"the csv 'removed' column on {a['n_removed_differs_from_csv']} rows; jj-rule j differs from the csv 'j' "
             f"on {a['n_j_differs_from_csv']} rows")
    L.append(f"- sigma1 in fp_zexp_rows_coefficients.csv: {a['n_sigma_nan']} NaN; {a['n_sigma_outside_0_Nm1']} "
             f"outside [0, N-1]: {a['sigma_outside_rows']}")
    L.append(f"- rows with 1s/2s removal (QED/FNS term, not computable from the paper): {a['n_qed_rows']}")
    L.append(f"- rows with an electron in no group under the literal wording: {a['n_nogroup_rows_literal']}; "
             f"d/f-target rows with same-n higher-l electrons: {a['n_dfout_rows']}")
    for m, r in a["models"].items():
        L.append(f"\n### {m}\n")
        if m == "pocket":
            L.append(f"- full-precision paper recomputation vs code (all rows): max rel. diff "
                     f"{r['FP_vs_code_max_rel']:.1e}")
        else:
            L.append(f"- h = 0 rows (N = 1, paper formula 0/0): {r['n_h_zero']}; T < 0 rows: {r['n_T_negative']}")
            L.append(f"- mu_Z = code / paper-only value is constant per element on reproduced rows (max dev "
                     f"{r['reproduced_max_dev']:.1e}); range {r['mu_Z_range'][0]:.9f}-{r['mu_Z_range'][1]:.9f} "
                     f"over {r['n_Z_with_mu']} elements")
            L.append(f"- QED/FNS term extracted as mu_Z * B - code on {r['qed_rows_extracted']} rows "
                     f"({r['qed_rows_not_extractable_Zle4']} rows with Z <= 4 not extractable): max "
                     f"{r['qed_max_rel_pct']:.4f} % of IE at {r['qed_max_rel_row']}, max {r['qed_max_abs_eV']:.2f} eV "
                     f"at {r['qed_max_abs_row']}, median {r['qed_median_rel_pct']:.4f} %, "
                     f"{r['qed_n_above_1e-4_rel']} rows above 0.01 %")
        L.append(f"- rows NOT reproduced from the paper + named csv: {r['n_not_reproduced']} "
                 f"(NaN sigma1: {r['nan_rows']})")
        if r["not_reproduced"]:
            c1 = sum(1 for x in r["not_reproduced"] if x["csv_rearranged"])
            c2 = sum(1 for x in r["not_reproduced"] if x["diff_rule_rearranged"])
            c3 = sum(1 for x in r["not_reproduced"] if x["Nm1_missing"])
            mx = max(r["not_reproduced"], key=lambda x: abs(x["ratio_dev"]))
            L.append(f"  - of the non-NaN ones: {c1} flagged rearranged in the csv, {c2} genuinely rearranged, "
                     f"{c3} with the (N-1) ion missing; largest IE deviation {100 * mx['ratio_dev']:+.1f} % at "
                     f"{mx['row']}")
            L.append("  - rows: " + ", ".join(f"{x['row'][0]},{x['row'][1]}" for x in r["not_reproduced"]))
        L.append(f"\nVariant effects on the {r['n_stat_rows']} reproduced rows (relative change applied to the code "
                 f"prediction; MAPE over all 5847 rows, code MAPE {r['MAPE_code']:.4f} %):\n")
        L.append("| variant | max abs change (%) | at (Z,N) | max abs change (eV) | at (Z,N) | rows > 1e-6 | "
                 "rows > 1e-4 | rows > 1e-3 | MAPE shift (% points) |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for v, d in r.items():
            if isinstance(d, dict) and "max_rel_pct" in d:
                L.append(f"| {v} | {d['max_rel_pct']:.2e} | {d['row_max_rel']} | {d['max_abs_eV']:.2e} | "
                         f"{d['row_max_abs']} | {d['n_rel_gt_1e6']} | {d['n_rel_gt_1e4']} | {d['n_rel_gt_1e3']} | "
                         f"{d['dMAPE']:+.5f} |")
    return "\n".join(L)


# ----------------------------------------------------------------------------- main
VARIANTS = [
    ("kappa_fit", "kappa as printed in the fit column (1.7527; pocket 9.815) instead of |kappa|+0.05", dict(kappa_fit=True)),
    ("mu1", "mu(Z) = 1 (paper gives no definition of the mass)", dict(mu1=True)),
    ("Ry_paper", "Ry = 13.6057 eV (paper sec. 2) instead of CODATA 13.605693123", dict(Ry=None)),
    ("alpha137", "alpha = 1/137.036 (unified.md sec. 6) instead of CODATA", dict(alpha=ALPHA_137)),
    ("sig4", "sigma1 rounded to 4 decimals (Table 5 precision)", dict(sig_nd=4)),
    ("no_classes", "class deviations omitted (nu_c never defined in the paper)", dict(no_classes=True)),
    ("G_DBL", "'same-n s' read as also counting the s target's own partner (in += k-1)", dict(groups="G_DBL")),
    ("j_after", "jj rule applied to the occupancy after removal (k-1)", dict(j_after=True)),
    ("no_dirac", "Dirac factor omitted (=1)", dict(no_dirac=True)),
    ("params4dp", "Table A1 / Appendix A precision (4 dp; pocket kappa, x 3 dp) instead of full", dict()),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-all", action="store_true")
    ap.add_argument("--json")
    args = ap.parse_args()
    pp_raw = read_paper()
    full = read_full()
    pp = paper_params(pp_raw)
    # pa_bound9 is not printed: simulate Table A1 precision
    pp["pa_bound9"] = rounded(full["pa_bound9"], 4)
    nist = load_nist()
    sig, sigtab = load_sigma()
    audit = json.load(open(AUDIT, encoding="utf-8"))
    code_by = {(d["Z"], d["N"], d["model"]): d for d in audit}
    preds = {m: load_pred(p) for m, p in PRED.items()}
    P4opt = dict(Ry=pp_raw["Ry"], groups="G_LIT")         # "as printed": paper Ry, literal groups
    FPopt = dict(groups="G_LIT")                          # full precision, CODATA Ry
    res = {}
    for Z, N in CASES:
        f = features(Z, N, nist, sig)
        r = dict(f=f, P4={}, FP={}, code={}, pred={m: preds[m][(Z, N)] for m in MODELS}, var={})
        for m in MODELS:
            r["P4"][m] = compute(m, f, pp[m], P4opt)
            r["FP"][m] = compute(m, f, full[m], FPopt)
            if (Z, N, m) in code_by:
                r["code"][m] = dict(code_by[(Z, N, m)])
            base = r["FP"][m]["IE"]
            for name, _, o in VARIANTS:
                if m != "pa_hier_rel" and name == "no_classes":
                    continue
                if name == "params4dp":
                    val = compute(m, f, pp[m], FPopt)["IE"]
                elif name == "Ry_paper":
                    val = compute(m, f, full[m], dict(FPopt, Ry=pp_raw["Ry"]))["IE"]
                else:
                    val = compute(m, f, full[m], dict(FPopt, **o))["IE"]
                r["var"][(m, name)] = val - base
        r["mu_code"] = preds["pa_hier_rel"][(Z, N)] / (r["FP"]["pa_hier_rel"]["IE"] / r["FP"]["pa_hier_rel"]["mu"])
        res[(Z, N)] = r

    out = []
    pr = out.append
    pr("# Hand recomputation from the paper text (generated by tools/hand_recompute.py)\n")
    pr("## Parameters parsed verbatim from docs/paper_draft.md\n")
    pr(f"- Ry = {pp_raw['Ry']} eV (sec. 2); alpha not printed (CODATA 2018 {ALPHA} used); m_e = {ME_U} u (CODATA 2018)")
    pr(f"- tau_g = {pp_raw['tau']}; kappa fit {pp_raw['kappa_fit']}, as used {pp_raw['kappa_used']}")
    pr(f"- dtau_c = {pp_raw['dtau']}")
    pr(f"- r_c = {pp_raw['r']}; x_l = {pp_raw['x']}")
    pr(f"- pocket = {pp_raw['pocket']}\n")
    pr("## sigma1 inputs (results/fp_zexp_rows_coefficients.csv) and cross-checks\n")
    pr("| ion | config | removed | k | j | Za | sigma1 (rows csv) | sigma1 (neutral table) | rows csv 'removed' agrees |")
    pr("|---|---|---|---|---|---|---|---|---|")
    for (Z, N), r in res.items():
        f = r["f"]
        st = sigtab.get(N, {}).get("sigma1") if Z == N else None
        pr(f"| {SYM[Z]} Z={Z} N={N} | {sig[(Z, N)]['config']} | {f['n']}{LL[f['l']]} | {f['k']} | "
           f"{jj_j(f['l'], f['k'])} (csv {sig[(Z, N)]['j']}) | {f['Za']} | {f['sigma1']:.10f} | "
           f"{'-' if st is None else f'{st:.10f}'} | {sig[(Z, N)]['removed'] == (f['n'], f['l'])} |")
    pr("")
    pr("## Component tables (P4 = as printed: Table A1/App. A values, paper Ry, kappa 'as used', mu from A m_u;"
       " FP = full-precision parameters, CODATA Ry; code = results/audit_components.json)\n")
    maxdiff = defaultdict(float)
    for (Z, N), r in res.items():
        for m in MODELS:
            if m not in r["code"]:
                continue
            pr(f"### {SYM[Z]} (Z={Z}, N={N}), {m}\n")
            pr(comp_table((Z, N), r["P4"][m], r["FP"][m], r["code"][m]))
            pr("")
            maxdiff[(m, "P4")] = max(maxdiff[(m, "P4")], abs(r["P4"][m]["IE"] - r["code"][m]["IE"]))
            maxdiff[(m, "FP")] = max(maxdiff[(m, "FP")], abs(r["FP"][m]["IE"] - r["code"][m]["IE"]))
            maxdiff[(m, "FPrel")] = max(maxdiff[(m, "FPrel")], abs(r["FP"][m]["IE"] / r["code"][m]["IE"] - 1))
    pr("## IE summary (eV)\n")
    pr("| ion | model | NIST | code | as printed (P4) | P4 - code | full precision (FP) | FP - code |")
    pr("|---|---|---|---|---|---|---|---|")
    for (Z, N), r in res.items():
        for m in MODELS:
            c = r["code"][m]["IE"] if m in r["code"] else r["pred"][m]
            pr(f"| {SYM[Z]} {Z},{N} | {m} | {nist[(Z, N)]['IE']:.4f} | {c:.6f} | {r['P4'][m]['IE']:.6f} | "
               f"{r['P4'][m]['IE'] - c:+.6f} | {r['FP'][m]['IE']:.6f} | {r['FP'][m]['IE'] - c:+.2e} |")
    pr("")
    pr("## Sensitivity: IE change (eV) from full precision for each reader choice\n")
    pr("| variant | ion | pa_hier_rel | pa_bound9 | pocket |")
    pr("|---|---|---|---|---|")
    for name, desc, _ in VARIANTS:
        for (Z, N), r in res.items():
            vals = [r["var"].get((m, name)) for m in MODELS]
            pr(f"| {name}: {desc} | {SYM[Z]} {Z},{N} | " + " | ".join("-" if v is None else f"{v:+.4f}" for v in vals)
               + " |")
    pr("")
    wex = worked_examples(res, nist)
    pr("## Worked examples (paper sec. 4.7, docs/unified.md sec. 6)\n")
    pr("| doc | quantity | printed | quote found | printed arithmetic gives | arith OK | code value | code OK |")
    pr("|---|---|---|---|---|---|---|---|")
    for w in wex:
        ar = "" if w["arith"] is None else f"{w['arith']:.7g}"
        ok = "" if w["arith_ok"] is None else str(w["arith_ok"])
        pr(f"| {w['doc']} | {w['label']} | {w['printed']} | {w['found']} | {ar} | {ok} | {w['code']:.7g} | "
           f"{w['code_ok']} |")
    pr("")
    allr = None
    if not args.no_all:
        allr = all_rows(nist, sig, full, paper_params(pp_raw))
        pr("## 5847-row study\n")
        pr(all_rows_md(allr))
    print("\n".join(out))
    if args.json:
        dump = dict(cases={f"{Z},{N}": dict(P4=r["P4"], FP=r["FP"], var={f"{a}|{b}": v for (a, b), v in r["var"].items()},
                                            pred=r["pred"], mu_code=r["mu_code"]) for (Z, N), r in res.items()},
                    worked=wex, all_rows=allr, maxdiff={f"{a}|{b}": v for (a, b), v in maxdiff.items()})
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(dump, f, indent=1, default=str)


if __name__ == "__main__":
    main()
