r"""Score v2 candidates under docs/preregistration_v2.md.   py -3.11 models/v2/run.py <variant> [<variant> ...]

For each variant: refit on V1, V2, S1, S3 train rows -> selection score (mean test MAPE); all-data fit -> MAPE,
neutral MAPE, physical checks on the 7021 coverage cases (IE finite and > 0, monotonicity violations); then, AFTER
the above, the S2 refit (non-blind history). Appends to results/v2_log.json; never touches v1 files.
"""
import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _HERE)
import model as M  # noqa: E402
import umodel as U  # noqa: E402

BASE = dict(grouping="pocket", kappa_l=False, exch=True, qed=True, bound=True, tie={}, hier=True, ridge=1e-4)
VARIANTS = {
    "v1_baseline": dict(BASE, rel="fs"),                 # must reproduce the pre-registered baseline 2.1473
    "relexp": dict(BASE, rel="exp"),                      # item 1, variant 1
    "relexp_q2": dict(BASE, rel="exp", q2=True),          # item 3 (neutrals), variant 1
    "relexp_kl": dict(BASE, rel="exp", kappa_l=True),     # item 3 (neutrals), variant 2
    "relexp_q2_kl": dict(BASE, rel="exp", q2=True, kappa_l=True),  # item 3, variant 3
}
LOG = os.path.join(_ROOT, "results", "v2_log.json")

a = U.dataset()
y = a["y"]
Z = a["Z"].astype(int)
N = (a["Z"] - a["Zq"] + 1).astype(int)
SPLITS = {"V1": (Z <= 36, (Z >= 37) & (Z <= 54)), "V2": (Z <= 44, (Z >= 45) & (Z <= 54)),
          "S1": (Z % 5 != 0, Z % 5 == 0), "S3": (N % 6 != 0, N % 6 == 0)}
S2 = (Z <= 54, Z >= 55)
_cov = None


def coverage():
    global _cov
    if _cov is None:
        zn = [(z, n) for z in range(1, 119) for n in range(1, z + 1)]
        _cov = (zn, U.build(zn))
    return _cov


def met(p, m):
    ape = np.abs(p - y) / y * 100
    neu = m & (Z == N)
    return {"n": int(m.sum()), "MAPE": round(float(ape[m].mean()), 4), "median": round(float(np.median(ape[m])), 4),
            "neutral_MAPE": round(float(ape[neu].mean()), 3) if neu.any() else None}


def run(name):
    t = time.time()
    spec = VARIANTS[name]
    out = {"variant": name, "spec": {k: v for k, v in spec.items() if k != "tie"}, "splits": {}}
    for k, (tr, te) in SPLITS.items():
        s, th = M.fit(tr, spec)
        out["splits"][k] = met(M.evaluate(th, s, a), te)
    out["selection_score"] = round(float(np.mean([out["splits"][k]["MAPE"] for k in SPLITS])), 4)
    s, th = M.fit(None, spec)
    out["n_params"] = len(th)
    out["all_data"] = met(M.evaluate(th, s, a), np.ones_like(Z, bool))
    zn, ca = coverage()
    ie = M.evaluate(th, s, ca)
    d = {k: v for k, v in zip(zn, ie)}
    out["coverage"] = {"n": len(ie), "finite_and_positive": int(np.sum(np.isfinite(ie) & (ie > 0))),
                       "monotonicity_violations": int(sum(1 for (z, n) in zn if n > 1 and d[(z, n - 1)] <= d[(z, n)])),
                       "Og_I_eV": round(float(d[(118, 118)]), 3), "Lr_I_eV": round(float(d[(103, 103)]), 3)}
    out["params"] = dict(zip(M.param_names(s), [float(v) for v in th]))
    # S2: computed last, non-blind history only (failures known before v2 development)
    s2s, s2th = M.fit(S2[0], spec)
    p2 = M.evaluate(s2th, s2s, a)
    out["S2_nonblind"] = met(p2, S2[1])
    out["S2_nonblind"]["Lr_I_eV"] = round(float(M.evaluate(s2th, s2s, U.build([(103, 103)]))[0]), 3)
    out["seconds"] = round(time.time() - t, 1)
    log = json.load(open(LOG, encoding="utf-8")) if os.path.exists(LOG) else []
    log.append(out)
    json.dump(log, open(LOG, "w", encoding="utf-8"), indent=1)
    sp = out["splits"]
    print(f"{name:14s} p={out['n_params']:2d} score {out['selection_score']:.4f} | V1 {sp['V1']['MAPE']:.3f} V2 "
          f"{sp['V2']['MAPE']:.3f} S1 {sp['S1']['MAPE']:.3f} S3 {sp['S3']['MAPE']:.3f} | all {out['all_data']['MAPE']:.3f} "
          f"neutral {out['all_data']['neutral_MAPE']:.2f} | cov>0 {out['coverage']['finite_and_positive']}/7021 "
          f"mono {out['coverage']['monotonicity_violations']} Og {out['coverage']['Og_I_eV']} | S2(nonblind) "
          f"{out['S2_nonblind']['MAPE']:.2f} neu {out['S2_nonblind']['neutral_MAPE']:.1f} Lr {out['S2_nonblind']['Lr_I_eV']} "
          f"({out['seconds']} s)", flush=True)


if __name__ == "__main__":
    for v in sys.argv[1:]:
        run(v)
