"""Final fitting / validation pipeline for the GSHM (Track B).

Run:  py -3.13 models/semi_empirical/fit_gshm.py
Writes results/se_params.json, results/se_gshm*_predictions.csv, results/se_gshm*_metrics.json,
results/se_validation.json, results/figures/se_gshm_*.png
"""
import csv
import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import gshm  # noqa: E402
import fitlib as F  # noqa: E402
import evaluate  # noqa: E402

RES = os.path.join(_ROOT, "results")
FIG = os.path.join(RES, "figures")

# ridge 0.1 towards the Slater-like initial values (chosen by an inner extrapolation check,
# train Z<=44 / validate 45<=Z<=54) and a physical bound 0.8 <= sigma_core <= 1 on deep-core
# screening (not identifiable from Z<=54 data, where core counts are constant).
CORE = F.with_spec(gshm.default_spec(), S1=True, S1_den="Za", Za_mode="true", ridge=0.1,
                   core_bounds=[0.8, 1.0])
LIBRARY = list(gshm.dataset("1e")["_feats"].keys())


def fsel(train_mask, val_mask, max_terms=6, min_gain=0.01, log=print):
    """Greedy forward selection (orthogonal-matching-pursuit style) of correction terms.

    At each step every candidate phi_t is tried by refitting ONLY the linear correction
    coefficients beta on the log-residuals of the current model (train set, linear least squares);
    the candidate with the lowest validation MAPE is kept and then the whole model (core + all
    chosen beta) is refitted jointly by nonlinear least squares on the train set."""
    feats = gshm.dataset("1e")["_feats"]
    chosen = []
    spec = F.with_spec(CORE, corr=[])
    th, _ = F.fit(spec, train_mask)
    base_pred = F.predict_all(spec, th)
    best = F.metrics(base_pred, val_mask)["MAPE"]
    path = [dict(terms=[], nparam=len(th), val_MAPE=best)]
    log(f"  start val MAPE {best:.4f} ({len(th)} params)")
    y = F.YARR
    for step in range(max_terms):
        r = np.log(y / base_pred)
        trial = []
        for t in LIBRARY:
            if t in chosen:
                continue
            X = feats[t][:, None]
            b, *_ = np.linalg.lstsq(X[train_mask], r[train_mask], rcond=None)
            p = base_pred * np.exp(X @ b)
            trial.append((F.metrics(p, val_mask)["MAPE"], t))
        trial.sort()
        v, t = trial[0]
        if v > best * (1 - min_gain):
            log(f"  stop: best candidate {t} gives {v:.4f} (no >{min_gain:.0%} gain)")
            break
        chosen.append(t)
        spec = F.with_spec(CORE, corr=list(chosen))
        th, _ = F.fit(spec, train_mask, theta0=np.concatenate([th, [0.0]]), max_nfev=300)
        base_pred = F.predict_all(spec, th)
        v2 = F.metrics(base_pred, val_mask)["MAPE"]
        best = min(v2, best)
        path.append(dict(terms=list(chosen), nparam=len(th), val_MAPE=v2))
        log(f"  + {t:22s} val MAPE {v:.4f} (linear screen) -> {v2:.4f} (joint refit)")
    return chosen, path


def param_table(spec, theta, res):
    names = gshm.param_names(spec)
    J = res.jac
    m, p = J.shape
    s2 = float(np.sum(res.fun ** 2) / max(m - p, 1))
    try:
        cov = np.linalg.pinv(J.T @ J) * s2
        unc = np.sqrt(np.clip(np.diag(cov), 0, None))
    except Exception:  # pragma: no cover
        unc = np.full(p, np.nan)
    return [dict(name=n, value=float(v), unc=float(u)) for n, v, u in zip(names, theta, unc)]


def write_preds(name, pred):
    path = os.path.join(RES, f"se_{name}_predictions.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Z", "N", "IE_pred_eV", "IE_ref_eV"])
        for Z, N, p, y in zip(F.ZARR, F.NARR, pred, F.YARR):
            w.writerow([int(Z), int(N), f"{p:.8g}", y])
    m = evaluate.evaluate(path)
    with open(os.path.join(RES, f"se_{name}_metrics.json"), "w") as fh:
        json.dump(m, fh, indent=1)
    return path, m


def main():
    t0 = time.time()
    out = {}
    # ---------------------------------------------------------------- 1. term selection
    # selection pool = S1 train; inner validation = Z % 5 == 2 (disjoint from the S1 test set)
    print("forward selection (train = S1 train minus Z%5==2, validation = Z%5==2)")
    pool = ~F.SPLITS["S1"]
    val = pool & (F.ZARR % 5 == 2)
    chosen, path = fsel(pool & ~val, val)
    out["selection"] = dict(chosen=chosen, path=path)
    FINAL = F.with_spec(CORE, corr=chosen)
    print("time", round(time.time() - t0), "s")

    # ---------------------------------------------------------------- 2. held-out splits
    out["splits"] = {}
    for name, spec in (("core", CORE), ("final", FINAL)):
        for sp in ("S1", "S2", "S3"):
            m, th, _ = F.cv_split(spec, sp)
            out["splits"][f"{name}_{sp}"] = m
            print(f"{name:5s} {sp}: test MAPE {m['MAPE']:.3f} median {m['median']:.3f} "
                  f"neutral {m['neutral_MAPE']:.3f} (n={m['n']}, neutral n={m['n_neutral']})")
    for label, spec, tie in (("core_S2_untied_noprior", F.with_spec(CORE, ridge=0.0, core_bounds=None), False),
                             ("core_S2_tied_noprior", F.with_spec(CORE, ridge=0.0, core_bounds=None), True),
                             ("core_S2_tied_ridge_only", F.with_spec(CORE, core_bounds=None), True)):
        m, _, _ = F.cv_split(spec, "S2", tie=tie)
        out["splits"][label] = m
        print(f"{label:28s}: S2 test MAPE {m['MAPE']:.3f} median {m['median']:.3f} "
              f"neutral {m['neutral_MAPE']:.3f}")
    for sp in ("S1", "S2", "S3"):
        print(sp, "tied parameters:", F.auto_tie(FINAL, ~F.SPLITS[sp]).get("tie", {}))
    print("time", round(time.time() - t0), "s")

    # ---------------------------------------------------------------- 3. ablation (all data)
    abl = {
        "A0 screening only (sigma_c)": F.with_spec(CORE, S1=False, rel=False, exch=False),
        "A1 + q-dependent screening S1": F.with_spec(CORE, rel=False, exch=False),
        "A2 + Dirac factor (no penetration corr.)": F.with_spec(CORE, exch=False, rel_pen=False),
        "A3 + relativistic penetration r_lj": F.with_spec(CORE, exch=False),
        "A4 + Hund exchange x_l (= core GSHM)": CORE,
        "A5 + sparse corrections (= final GSHM)": FINAL,
        "-- core minus S1": F.with_spec(CORE, S1=False),
        "-- core minus relativity": F.with_spec(CORE, rel=False),
        "-- core minus exchange": F.with_spec(CORE, exch=False),
        "-- core with S1 denominator Z_eff0 (not Z_a)": F.with_spec(CORE, S1_den="Ze0"),
        "-- core, total-energy-difference form": None,
    }
    out["ablation"] = {}
    for k, s in abl.items():
        if s is None:
            continue
        th, _ = F.fit(s, max_nfev=1500)
        m = F.metrics(F.predict_all(s, th))
        mS1, _, _ = F.cv_split(s, "S1", theta0=th)
        out["ablation"][k] = dict(nparam=len(th), all=m, S1_test=mS1)
        print(f"{k:48s} p={len(th):3d} all MAPE {m['MAPE']:.3f} med {m['median']:.3f} "
              f"neutral {m['neutral_MAPE']:.3f} | S1 test {mS1['MAPE']:.3f}")
    print("time", round(time.time() - t0), "s")

    # ---------------------------------------------------------------- 4. final fit on all data
    params = {}
    preds = {}
    for name, spec in (("gshm_core", CORE), ("gshm", FINAL)):
        th, res = F.fit(spec, max_nfev=2000)
        pred = F.predict_all(spec, th)
        preds[name] = pred
        _, m = write_preds(name, pred)
        params[name] = dict(spec=spec, values=[float(v) for v in th],
                            table=param_table(spec, th, res), n_params=len(th),
                            metrics_all={k: m[k] for k in ("ALL", "first_IE_neutral_atoms",
                                                            "hydrogen_like")})
        evaluate.print_table(m, name)
    final = dict(params["gshm"])
    doc = {"description": "Generalized Screened-Hydrogenic Model (Track B). 'final' is used by "
                          "models/semi_empirical/gshm.py predict(); 'core' omits the sparse "
                          "correction terms. Uncertainties are formal 1-sigma from the Jacobian "
                          "of the log-residuals (heavy-tailed residuals: indicative only).",
           "final": final, "core": params["gshm_core"]}
    with open(os.path.join(RES, "se_params.json"), "w") as fh:
        json.dump(doc, fh, indent=1)
    with open(os.path.join(RES, "se_validation.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    # sanity: predict() reproduces the vectorized predictions
    gshm._load_params.cache_clear()
    chk = [abs(gshm.predict(int(Z), int(N)) / p - 1) for Z, N, p in
           zip(F.ZARR[::97], F.NARR[::97], preds["gshm"][::97])]
    print("predict() consistency max rel diff:", max(chk))
    print("predict(120, 120) [outside table, Madelung config] =", gshm.predict(120, 120))
    print("predict(26, 26, shells=[(1,0,2),(2,0,2),(2,1,6),(3,0,2),(3,1,6),(3,2,6),(4,0,2)]) =",
          gshm.predict(26, 26, shells=[(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 6),
                                       (3, 2, 6), (4, 0, 2)]))

    # ---------------------------------------------------------------- 5. figures
    figures(preds["gshm"], path)
    print("total time", round(time.time() - t0), "s")


def figures(pred, path):
    y = F.YARR
    lr = np.log(pred / y) * 100
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
    ax = axs[0]
    ax.loglog(y, pred, ".", ms=2, alpha=0.5)
    ax.plot([3, 3e5], [3, 3e5], "k-", lw=0.7)
    ax.set_xlabel("NIST IE (eV)")
    ax.set_ylabel("GSHM IE (eV)")
    ax.set_title("Final GSHM vs NIST (5847 IEs)")
    ax = axs[1]
    sc = ax.scatter(F.ZARR, F.NARR, c=np.clip(lr, -10, 10), cmap="RdBu_r", s=4)
    plt.colorbar(sc, ax=ax, label="100 ln(pred/NIST)  (clipped at +-10)")
    ax.set_xlabel("Z")
    ax.set_ylabel("N")
    ax.set_title("Residual map")
    ax = axs[2]
    neu = F.ZARR == F.NARR
    ax.plot(F.ZARR[neu], y[neu], "k.-", label="NIST")
    ax.plot(F.ZARR[neu], pred[neu], "r.-", label="GSHM")
    ax.set_xlabel("Z")
    ax.set_ylabel("first IE of neutral atom (eV)")
    ax.legend()
    ax.set_title("Neutral-atom first IEs")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "se_gshm_overview.png"), dpi=130)
    plt.close(fig)
    if path and len(path) > 1:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot([p["nparam"] for p in path], [p["val_MAPE"] for p in path], "o-")
        for p in path[1:]:
            ax.annotate(p["terms"][-1], (p["nparam"], p["val_MAPE"]), fontsize=7)
        ax.set_xlabel("number of global parameters")
        ax.set_ylabel("inner-validation MAPE (%)")
        ax.set_title("Forward selection path (learning curve)")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "se_gshm_selection_path.png"), dpi=130)
        plt.close(fig)


if __name__ == "__main__":
    main()
