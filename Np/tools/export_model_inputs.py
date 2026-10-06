r"""Export every per-row input of the Screened Rydberg formula to results/model_inputs.csv.

    py -3.11 tools/export_model_inputs.py

The paper prints the equation and the parameters. This table holds the per-row quantities a reader cannot
derive from the paper text alone:
  * sigma1 as used by the model: the exact first-order screening of the removed (n,l) electron, computed for the
    FROZEN configuration (N-electron ground configuration minus one (n,l) electron), not for the NIST ground
    configuration of the ion. The two differ only for rearranged rows.
  * removed subshell, its occupancy k, j (jj filling) and the relativistic class;
  * the five group counts nu_g and the 19 class counts nu_c (the grouping rule of the paper, Table A2);
  * the Hund kink K_l(k);
  * mu(Z) = M/(M + m_e), with the nuclear mass M/m_e of the isotope used in the Yerokhin-Shabaev 2015 table
    (models/first_principles/cache/yerokhin_shabaev_2015_qed.json);
  * the QED + finite-nuclear-size shift (eV) of a 1s or 2s level of charge Z (zero for every other removal).
Inputs come from the production feature builder (models/unified/umodel.build); no NIST ionization energy is read.
tools/verify_from_inputs.py rebuilds the model from this table and the printed equation only.
"""
import csv
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _p in (_ROOT, os.path.join(_ROOT, "models", "unified"), os.path.join(_ROOT, "models", "first_principles"),
           os.path.join(_ROOT, "models", "semi_empirical")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import umodel as U  # noqa: E402
import gshm  # noqa: E402
from push_a_import import PA  # noqa: E402
from configs import initial_final  # noqa: E402
from common.atomdata import load_records  # noqa: E402

LL = "spdf"
GROUPS = ["same", "in", "core", "df", "out"]


def main():
    recs = load_records()
    ZN = [(r["Z"], r["N"]) for r in recs]
    a = U.build(ZN)
    spec = {"grouping": "pocket"}
    G = PA._group_matrix(spec, a)
    C = PA._dev_matrix(a)
    out = os.path.join(_ROOT, "results", "model_inputs.csv")
    head = (["Z", "N", "config_N", "removed", "n", "l", "k", "j", "rel_class", "rearranged", "Za", "sigma1"]
            + ["nu_" + g for g in GROUPS] + ["nuc_" + c for c in PA.DEV_CLASSES] + ["K", "mu", "qedfns_eV"])
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(head)
        for i, (Z, N) in enumerate(ZN):
            sN, _, rem, rearr = initial_final(Z, N)
            cfg = ".".join(f"{n}{LL[l]}{q}" for n, l, q in sN)
            w.writerow([Z, N, cfg, f"{rem[0]}{LL[rem[1]]}", int(a["n"][i]), int(a["l"][i]), int(a["k"][i]),
                        a["j"][i], PA.REL_CLASSES[int(a["rc"][i])], int(bool(rearr)), max(Z - N + 1, 1),
                        repr(float(a["sigma1"][i]))]
                       + [int(v) for v in G[i]] + [int(v) for v in C[i]]
                       + [repr(float(a["kink"][i])), repr(float(a["mu"][i])), repr(float(a["qedfns"][i]))])
    print(f"wrote {out}: {len(ZN)} rows, {len(head)} columns")


if __name__ == "__main__":
    main()
