"""Known failures of the final model, recomputed (paper abstract, sec. 4.3, 4.7, 5.2).

    py -3.11 tools/known_failures.py   -> results/known_failures.json

S2 refit (Z <= 54) of pa_hier_rel, in memory only: first IEs of Tl, Pb, Rn, Lr.
All-data fit (frozen parameters): Og, Rn, Ca, and the 4s2 neutrals' errors. Nothing frozen is written.
"""
import json
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for _p in (_ROOT, os.path.join(_ROOT, "models", "unified"), os.path.join(_ROOT, "models", "first_principles"),
           os.path.join(_ROOT, "models", "semi_empirical")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import numpy as np  # noqa: E402
import umodel as U  # noqa: E402
from push_a_import import PA  # noqa: E402
from common.atomdata import load_records  # noqa: E402

NAMES = {20: "Ca", 38: "Sr", 56: "Ba", 81: "Tl", 82: "Pb", 86: "Rn", 103: "Lr", 118: "Og"}


def main():
    nist = {(r["Z"], r["N"]): r["IE_eV"] for r in load_records()}
    c = json.load(open(os.path.join(_ROOT, "results", "pa_params.json"), encoding="utf-8"))["candidates"]["pa_hier_rel"]
    base = {k: v for k, v in c["spec"].items() if k not in ("tie", "dev")}
    d = U.dataset()
    s2spec, s2th = PA.fit(np.array(d["rec_Z"]) <= 54, dict(base))
    allspec, allth = c["spec"], np.array(c["values"])
    out = {"S2_fit_Zle54": {}, "all_data_fit": {}}
    for Z in (81, 82, 86, 103):
        v = float(PA.evaluate(s2th, s2spec, U.build([(Z, Z)]))[0])
        out["S2_fit_Zle54"][NAMES[Z]] = {"Z": Z, "pred_eV": round(v, 3), "NIST_eV": nist.get((Z, Z))}
    for Z in (20, 38, 56, 86, 118):
        v = float(PA.evaluate(allth, allspec, U.build([(Z, Z)]))[0])
        ref = nist.get((Z, Z))
        out["all_data_fit"][NAMES[Z]] = {"Z": Z, "pred_eV": round(v, 3), "NIST_eV": ref,
                                         "err_pct": round((v / ref - 1) * 100, 2) if ref else None}
    path = os.path.join(_ROOT, "results", "known_failures.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
