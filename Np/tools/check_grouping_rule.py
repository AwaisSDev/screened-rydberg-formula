r"""Independent check of the written screening-group / screening-class rule (paper §2.2, Table A2).

The rule (docs/review/grouping_rule.md) is implemented here from its own table, WITHOUT calling
gshm.classify(), PA._group_matrix(), PA._dev_matrix() or pocket.counts() to compute anything. Those
production functions are called only afterwards, to compare.

    export PYTHONPATH=tools/numba_stub    # only on machines without numba
    py -3.11 tools/check_grouping_rule.py

Checks
  0. The class table is a partition: for every target (n,l) and every other subshell (n',l')
     (n, n' = 1..8, l' < n', l, l' <= 3) exactly one table row applies.
  1. All 5847 rows of data/nist_ie.csv (NIST ground configuration; removed subshell from
     configs.initial_final) and all 7021 coverage cases Z = 1..118, N = 1..Z (umodel.build):
     group counts nu_g (5) and class counts nu_c (21 = 19 fitted + sn_out + out_sp) equal the code.
  2. nu_g and nu_c for O, Na, N, Ca, Fe, Pb vs results/audit_components.json.
  3. Paper statements of K_l(k) and of the jj j assignment vs configs.exchange_kink / configs.jj_j.
Nothing is written to disk.
"""
import json
import os
import sys
from fractions import Fraction

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _p in (_ROOT, os.path.join(_ROOT, "models", "unified"), os.path.join(_ROOT, "models", "first_principles"),
           os.path.join(_ROOT, "models", "semi_empirical")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import abinitio  # noqa: E402
abinitio.save_cache = lambda: None        # guard: this check never writes the sigma1 cache
import umodel as U  # noqa: E402
import pocket  # noqa: E402
import configs  # noqa: E402
from push_a_import import PA  # noqa: E402
from common.atomdata import load_records, shells_to_str  # noqa: E402

LL = "spdf"
GROUPS = ["same", "in", "core", "df", "out"]

# ------------------------------------------------------------------ the written rule (Table A2)
# Each row: (class, group, target types, condition on the other electron (n2, l2)).
# dn = n - n2 (dn < 0: other electron has larger n).  Target types: "s", "p", "d", "f".
SP, D_, F_, ANY = ("s", "p"), ("d",), ("f",), ("s", "p", "d", "f")
RULE = [
    # group same: the k-1 other electrons of the removed subshell
    ("same_s",   "same", ("s",), lambda n, l, n2, l2: (n2, l2) == (n, l)),
    ("same_p",   "same", ("p",), lambda n, l, n2, l2: (n2, l2) == (n, l)),
    ("same_d",   "same", ("d",), lambda n, l, n2, l2: (n2, l2) == (n, l)),
    ("same_f",   "same", ("f",), lambda n, l, n2, l2: (n2, l2) == (n, l)),
    # group out: n2 > n, or n2 = n with l2 > l
    ("sn_out",   "out",  ANY,    lambda n, l, n2, l2: n2 == n and l2 > l),
    ("out_sp",   "out",  SP,     lambda n, l, n2, l2: n2 > n),
    ("out_d",    "out",  D_,     lambda n, l, n2, l2: n2 > n),
    ("out_f",    "out",  F_,     lambda n, l, n2, l2: n2 > n),
    # s/p target, inner electrons with n2 >= n-1: group in
    ("sn_in_p",  "in",   SP,     lambda n, l, n2, l2: n2 == n and l2 < l),
    ("n1_sp_sp", "in",   SP,     lambda n, l, n2, l2: n2 == n - 1 and l2 <= 1),
    ("n1_sp_d",  "in",   SP,     lambda n, l, n2, l2: n2 == n - 1 and l2 == 2),
    ("n1_sp_f",  "in",   SP,     lambda n, l, n2, l2: n2 == n - 1 and l2 == 3),
    # s/p target, n2 <= n-2: group core
    ("n2_sp_sp", "core", SP,     lambda n, l, n2, l2: n2 == n - 2 and l2 <= 1),
    ("n2_sp_df", "core", SP,     lambda n, l, n2, l2: n2 == n - 2 and l2 >= 2),
    ("core_sp",  "core", SP,     lambda n, l, n2, l2: n2 <= n - 3),
    # d/f target: every inner electron (n2 < n, or n2 = n with l2 < l): group df
    ("d_near",   "df",   D_,     lambda n, l, n2, l2: n2 in (n, n - 1) and l2 <= 1),
    ("n1_d_d",   "df",   D_,     lambda n, l, n2, l2: n2 == n - 1 and l2 == 2),
    ("n1_d_f",   "df",   D_,     lambda n, l, n2, l2: n2 == n - 1 and l2 == 3),
    ("f_near",   "df",   F_,     lambda n, l, n2, l2: (n2 == n and l2 < l) or (n2 == n - 1 and l2 <= 2)),
    ("n1_f_f",   "df",   F_,     lambda n, l, n2, l2: n2 == n - 1 and l2 == 3),
    ("core_df",  "df",   D_ + F_, lambda n, l, n2, l2: n2 <= n - 2),
]
RULE_CLASSES = [r[0] for r in RULE]
CLASS_GROUP = {r[0]: r[1] for r in RULE}


def group_of(n, l, n2, l2):
    """Five-group rule in prose form (independent of the class table)."""
    if (n2, l2) == (n, l):
        return "same"
    if n2 > n or (n2 == n and l2 > l):
        return "out"                       # outer electron
    if l >= 2:
        return "df"                        # every inner electron of a d or f target
    return "in" if n2 >= n - 1 else "core"  # s or p target


def class_of(n, l, n2, l2):
    """Class from Table A2; asserts that exactly one row applies."""
    t = LL[l]
    hits = [c for c, _, types, cond in RULE if t in types and cond(n, l, n2, l2)]
    if len(hits) != 1:
        raise AssertionError(f"rule is not a partition at target {n}{t}, other {n2}{LL[l2]}: {hits}")
    return hits[0]


def rule_counts(shells, rem):
    """(nu_g dict, nu_c dict) for the electron removed from subshell rem of configuration shells."""
    n, l = rem
    g = dict.fromkeys(GROUPS, 0)
    c = dict.fromkeys(RULE_CLASSES, 0)
    for n2, l2, q2 in shells:
        m = q2 - 1 if (n2, l2) == (n, l) else q2
        if m <= 0:
            continue
        gg, cc = group_of(n, l, n2, l2), class_of(n, l, n2, l2)
        assert CLASS_GROUP[cc] == gg, (rem, (n2, l2), gg, cc)     # classes refine groups
        g[gg] += m
        c[cc] += m
    return g, c

# ------------------------------------------------------------------ comparison with the code


def code_matrices(ZN):
    a = U.build(ZN)
    spec = {"grouping": "pocket"}
    G = PA._group_matrix(spec, a)
    Cd = PA._dev_matrix(a)
    Gp = pocket.counts(a)
    return a, G, Cd, Gp


def compare(label, ZN):
    a, G, Cd, Gp = code_matrices(ZN)
    dev_cols = PA.DEV_CLASSES                       # 21 columns, core split into core_sp / core_df
    assert sorted(dev_cols) == sorted(RULE_CLASSES), (dev_cols, RULE_CLASSES)
    mism_g, mism_c, mism_rem, mism_pocket, sum_bad = [], [], [], [], 0
    Rg = np.zeros((len(ZN), 5))
    Rc = np.zeros((len(ZN), len(dev_cols)))
    rems = []
    for i, (Z, N) in enumerate(ZN):
        sN, _, rem, rearr = configs.initial_final(Z, N)
        rems.append((rem, rearr, sN))
        if (int(a["n"][i]), int(a["l"][i])) != tuple(rem):
            mism_rem.append((Z, N, rem, (a["n"][i], a["l"][i])))
        g, c = rule_counts(sN, rem)
        Rg[i] = [g[k] for k in GROUPS]
        Rc[i] = [c[k] for k in dev_cols]
        if Rg[i].sum() != N - 1:
            sum_bad += 1
        if not np.array_equal(Rg[i], G[i]):
            mism_g.append((Z, N, shells_to_str(sN), rem, Rg[i].tolist(), G[i].tolist()))
        if not np.array_equal(Rc[i], Cd[i]):
            diff = {dev_cols[j]: (Rc[i, j], Cd[i, j]) for j in range(len(dev_cols)) if Rc[i, j] != Cd[i, j]}
            mism_c.append((Z, N, shells_to_str(sN), rem, diff))
        if not np.array_equal(G[i], Gp[i]):
            mism_pocket.append((Z, N))
    print(f"\n== {label}: {len(ZN)} rows ==")
    print(f"removed subshell differs from the built features: {len(mism_rem)}")
    print(f"nu_g (5 groups)  mismatches vs PA._group_matrix(pocket): {len(mism_g)}")
    print(f"nu_c (21 classes) mismatches vs PA._dev_matrix:          {len(mism_c)}")
    print(f"PA._group_matrix vs pocket.counts (code vs code):         {len(mism_pocket)}")
    print(f"rows with sum_g nu_g != N-1:                             {sum_bad}")
    for m in (mism_g[:5] + mism_c[:5]):
        print("   example:", m)
    return Rg, Rc, rems, a


def stats(label, ZN, Rg, Rc, rems):
    dev_cols = PA.DEV_CLASSES
    print(f"\n-- {label}: rows (electrons) per group --")
    for j, g in enumerate(GROUPS):
        print(f"   {g:5s} rows={int((Rg[:, j] > 0).sum()):5d} electrons={int(Rg[:, j].sum()):7d}")
    print(f"-- {label}: rows (electrons) per class, Table A2 order --")
    for c in RULE_CLASSES:
        j = dev_cols.index(c)
        print(f"   {c:9s} {CLASS_GROUP[c]:5s} rows={int((Rc[:, j] > 0).sum()):5d} electrons={int(Rc[:, j].sum()):7d}")
    # edge cases named in the task
    tgt = {"s": 0, "p": 0, "d": 0, "f": 0}
    e = dict(s_with_same_n_p=[], p_with_same_n_d=[], d_with_same_n_f=[], df_same_n_lower_l=0,
             f_with_n1_d=0, n2_gt_n=0, n2_gt_n_sp=[], rearranged=0, sp_deep=0, df_deep=0)
    for i, ((n, l), rearr, sN) in enumerate(rems):
        tgt[LL[l]] += 1
        e["rearranged"] += bool(rearr)
        Z, N = ZN[i]
        sub = {(a, b) for a, b, q in sN}
        if l == 0 and any(a == n and b > 0 for a, b in sub):
            e["s_with_same_n_p"].append((Z, N))
        if l == 1 and any(a == n and b > 1 for a, b in sub):
            e["p_with_same_n_d"].append((Z, N))
        if l == 2 and any(a == n and b > 2 for a, b in sub):
            e["d_with_same_n_f"].append((Z, N))
        if l >= 2 and any(a == n and b < l for a, b in sub):
            e["df_same_n_lower_l"] += 1
        if l == 3 and (n - 1, 2) in sub:
            e["f_with_n1_d"] += 1
        if any(a > n for a, b in sub):
            e["n2_gt_n"] += 1
            if l <= 1:
                e["n2_gt_n_sp"].append((Z, N))
        if any(a <= n - 2 for a, b in sub):
            if l <= 1:
                e["sp_deep"] += 1
            else:
                e["df_deep"] += 1
    print(f"-- {label}: targets by l: {tgt}")
    for k, v in e.items():
        if isinstance(v, list):
            print(f"   {k}: {len(v)} rows {v[:12]}")
        else:
            print(f"   {k}: {v} rows")


# ------------------------------------------------------------------ Hund kink and jj j


def hund_pairs_indep(l, k):
    """Parallel-spin pairs of l^k in the Hund (maximum-S) arrangement: fill 2l+1 spin-up first."""
    up = min(k, 2 * l + 1)
    down = k - up
    return Fraction(up * (up - 1), 2) + Fraction(down * (down - 1), 2)


def check_kink_and_j():
    print("\n== K_l(k) = [P(k) - P(k-1)] - 2l(k-1)/(4l+1)  and  j = l-1/2 (k <= 2l) else l+1/2 ==")
    bad = 0
    for l in range(4):
        row_k, row_j = [], []
        for k in range(1, 4 * l + 3):
            K = (hund_pairs_indep(l, k) - hund_pairs_indep(l, k - 1)) - Fraction(2 * l * (k - 1), 4 * l + 1)
            j = (l - 0.5) if k <= 2 * l else (l + 0.5)
            Kc, jc = configs.exchange_kink(l, k), configs.jj_j(l, k)
            # closed form: (2l+1)(k-1)/(4l+1) up to half filling, antisymmetric K(k) = -K(4l+3-k) above
            kk = k if k <= 2 * l + 1 else 4 * l + 3 - k
            Kcf = Fraction((2 * l + 1) * (kk - 1), 4 * l + 1) * (1 if k <= 2 * l + 1 else -1)
            if Kcf != K:
                bad += 1
                print("   CLOSED-FORM MISMATCH", l, k, K, Kcf)
            if abs(float(K) - Kc) > 1e-12 or abs(j - jc) > 1e-12:
                bad += 1
                print("   MISMATCH", l, k, K, Kc, j, jc)
            row_k.append(str(K) if K.denominator != 1 else str(K.numerator))
            row_j.append(f"{j:g}")
        print(f"   l={LL[l]}: K = {', '.join(row_k)}")
        print(f"        j = {', '.join(row_j)}")
    paper_p = [0, 0.6, 1.2, -1.2, -0.6, 0]
    got_p = [configs.exchange_kink(1, k) for k in range(1, 7)]
    ok_p = all(abs(x - y) < 1e-12 for x, y in zip(paper_p, got_p))
    print(f"   paper p values {paper_p}: code {['%.4g' % v for v in got_p]} -> {'OK' if ok_p else 'MISMATCH'}")
    print(f"   rel class p1/p3 (gshm.rel_class) vs j: "
          f"{'OK' if all((configs.jj_j(1, k) == 0.5) == (__import__('gshm').rel_class(1, k) == 1) for k in range(1, 7)) else 'MISMATCH'}")
    print(f"   K/j mismatches: {bad}")
    return bad + (0 if ok_p else 1)


# ------------------------------------------------------------------ the six case studies


def case_studies():
    print("\n== case studies: rule counts vs results/audit_components.json ==")
    with open(os.path.join(_ROOT, "results", "audit_components.json"), encoding="utf-8") as f:
        aud = json.load(f)
    bad = 0
    for Z, N, sym in [(8, 8, "O"), (11, 11, "Na"), (7, 7, "N"), (20, 20, "Ca"), (26, 26, "Fe"), (82, 82, "Pb")]:
        sN, _, rem, rearr = configs.initial_final(Z, N)
        g, c = rule_counts(sN, rem)
        nz = {k: v for k, v in c.items() if v}
        print(f"   {sym:2s} {shells_to_str(sN)}  removed {rem[0]}{LL[rem[1]]} (k={dict(((a, b), q) for a, b, q in sN)[rem]})"
              f"{' rearranged' if rearr else ''}")
        print(f"      nu_g (same,in,core,df,out) = {[g[k] for k in GROUPS]}")
        print(f"      nu_c = {nz}")
        for d in aud:
            if d["Z"] != Z or d["N"] != N:
                continue
            ag = [int(d["nu_groups"][k]) for k in GROUPS]
            ok = ag == [g[k] for k in GROUPS]
            if d["model"] == "pa_hier_rel":
                ac = {k: int(v["nu"]) for k, v in d["class_devs"].items()}
                rc = {k: v for k, v in nz.items() if k not in ("sn_out", "out_sp")}
                ok = ok and ac == rc
            bad += not ok
            print(f"      audit {d['model']:11s}: nu_g {ag} {'OK' if ok else 'MISMATCH'}")
    return bad


def check_partition():
    n_cells = 0
    for l in range(4):
        for n in range(l + 1, 9):
            for n2 in range(1, 9):
                for l2 in range(0, min(n2, 4)):
                    class_of(n, l, n2, l2)       # raises unless exactly one row applies
                    if (n2, l2) != (n, l):
                        assert CLASS_GROUP[class_of(n, l, n2, l2)] == group_of(n, l, n2, l2)
                    n_cells += 1
    print(f"partition check: {n_cells} (target, other-subshell) cells, each matched by exactly one table row; "
          f"class -> group consistent")


def main():
    check_partition()
    recs = load_records()
    nist = [(r["Z"], r["N"]) for r in recs]
    cov = [(Z, N) for Z in range(1, 119) for N in range(1, Z + 1)]
    Rg, Rc, rems, _ = compare("NIST table (data/nist_ie.csv)", nist)
    stats("NIST table", nist, Rg, Rc, rems)
    with open(os.path.join(_ROOT, "results", "pa_params.json"), encoding="utf-8") as f:
        dev_fit = json.load(f)["candidates"]["pa_hier_rel"]["spec"]["dev"]
    populated = [c for j, c in enumerate(PA.DEV_CLASSES) if (Rc[:, j] > 0).any()]
    print(f"\nclasses populated in the 5847 rows: {len(populated)}; pa_hier_rel fitted dev list: {len(dev_fit)}; "
          f"identical sets: {sorted(populated) == sorted(dev_fit)}; empty: "
          f"{sorted(set(PA.DEV_CLASSES) - set(populated))}")
    Rg2, Rc2, rems2, _ = compare("coverage Z=1..118, N=1..Z", cov)
    stats("coverage", cov, Rg2, Rc2, rems2)
    nb = case_studies()
    kb = check_kink_and_j()
    print(f"\ncase-study mismatches: {nb}; K/j mismatches: {kb}")


if __name__ == "__main__":
    main()
