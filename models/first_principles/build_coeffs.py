"""Compute exact E0/E1 coefficients for every NIST row and write results/fp_zexp_coefficients.csv.

Columns per row (Z, N): configuration of the N- and (N-1)-electron ions, removed subshell, j,
dE0, dE1 (single-configuration Hund term, configuration average, Layzer complex), sigma (ab initio
leading-order screening constant), Slater's empirical screening constant for the same electron.
"""
import csv
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import zexp  # noqa: E402
from zexp import ROOT  # noqa: E402
from common.atomdata import load_records, ground_shells, shells_to_str, L_LETTER  # noqa: E402
import relativity as rel  # noqa: E402

SLATER_GROUPS = [(1, "s"), (2, "sp"), (3, "sp"), (3, "d"), (4, "sp"), (4, "d"), (4, "f"), (5, "sp"),
                 (5, "d"), (5, "f"), (6, "sp"), (6, "d"), (7, "sp")]


def _group(n, l):
    if l <= 1:
        return SLATER_GROUPS.index((n, "sp" if n > 1 else "s"))
    return SLATER_GROUPS.index((n, "d" if l == 2 else "f"))


def slater_sigma(shells, n, l):
    """Slater (1930) screening constant for one electron in subshell (n,l) of configuration."""
    g = _group(n, l)
    s = 0.0
    for n2, l2, q in shells:
        q2 = q - (1 if (n2, l2) == (n, l) else 0)
        g2 = _group(n2, l2)
        if g2 == g:
            s += q2 * (0.30 if n == 1 else 0.35)
        elif g2 < g:
            if l <= 1:
                s += q2 * (0.85 if n2 == n - 1 else (1.0 if n2 < n - 1 else 0.0))
                # same-n d/f groups that precede (none for sp) ; n2>=n with lower group impossible
            else:
                s += q2 * 1.0
    return s


def main():
    t0 = time.time()
    recs = load_records()
    rows = []
    for i, r in enumerate(recs):
        Z, N = r["Z"], r["N"]
        c = zexp.coeffs(Z, N)
        sh = r["shells"]
        shi = ground_shells(Z, N - 1) if N > 1 else []
        rows.append({
            "Z": Z, "N": N, "config": shells_to_str(sh), "config_ion": shells_to_str(shi) if shi else "-",
            "removed": f"{c['n']}{L_LETTER[c['l']]}", "j": f"{c['j2']}/2", "rearranged": int(r["rearranged"]),
            "dE0": repr(c["dE0"]), "dE1_sc": repr(c["dE1_sc"]), "dE1_av": repr(c["dE1_av"]),
            "dE1_cx": "" if c["dE1_cx"] is None else repr(c["dE1_cx"]), "dE1_best": repr(c["dE1"]),
            "sigma_abinitio": repr(c["sigma"]),
            "sigma_slater": repr(slater_sigma(sh, c["n"], c["l"])),
        })
        if i % 500 == 0:
            print(i, Z, N, time.time() - t0, flush=True)
            zexp.save_cfg_cache()
    zexp.save_cfg_cache()
    out = os.path.join(ROOT, "results", "fp_zexp_rows_coefficients.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote", out, time.time() - t0)
    # FNS cache for 1s, 2s
    for Z in range(1, 111):
        rel.fns_cached(1, -1, Z)
        rel.fns_cached(2, -1, Z)
    print("fns done", time.time() - t0)


if __name__ == "__main__":
    main()
