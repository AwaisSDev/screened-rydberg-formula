r"""One benchmark matrix for every model, built only from files that already exist in results/.

    export PYTHONPATH="D:/Chem-Research/Np/tools/numba_stub"     # only needed for --blind-check
    py -3.11 tools/benchmark_matrix.py                # writes results/benchmark_matrix.{csv,md,json}
    py -3.11 tools/benchmark_matrix.py --blind-check  # also prints the in-memory Z<=54 refit check (writes nothing)

What it does (no model is refit for the matrix):
  * In-sample block: every model's own prediction CSV scored with the shared scorer's statistics
    (evaluate._stats, the function evaluate.py itself uses; non-finite predictions skipped exactly as
    evaluate.evaluate does). Every cell carries its own n, so partial-coverage models (Mendoza 2011:
    5011 rows; LSDA dSCF: 207 rows) are never compared on silently different row sets.
    Each model's ALL cell is asserted equal to evaluate.evaluate(<csv>)["ALL"].
  * Held-out block: V1, V2, S1, S3, selection score and blind S2 (MAPE, median, neutral) COPIED from
    results/uni_validation.json, results/pa_validation.json and results/model_comparison.csv (first file
    that has the model, in that order; the other files are cross-checked). For models in none of those
    files (Kregar/Di Rocco SHM, Mendoza 2011) and for every 0-fitted-parameter model, the same row sets are
    also scored directly ("subset, no fit"); for those models nothing is held out.
  * Common-subset blocks: every model on the 5011 rows covered by Mendoza 2011, and on the 207 LSDA rows.
  * Ratios Slater / {final, pa_bound9, pocket} for the "N times" statements in the paper.

Strata (data/nist_ie.csv via common.atomdata.load_records):
  ALL; neutral = N == Z; H-like = N == 1; status = experimental / semi-empirical / theoretical;
  charge>=3 = Z - N >= 3; N<=10; Z>=55 (in-sample, NOT the blind S2 test).
"""
import argparse
import csv
import json
import math
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import evaluate as EV  # noqa: E402  (shared scorer, imported read-only)
from common.atomdata import load_records  # noqa: E402

RES = os.path.join(_ROOT, "results")
OUT = os.path.join(RES, "benchmark_matrix")

# key, label, prediction file stem, fitted params (as reported), kind, held-out source keys
MODELS = [
    ("bohr", "Bohr Ry Z^2/n^2", "se_bohr", "0", "parameter-free",
     {"mc": "se_bohr"}),
    ("slater", "Slater 1930, total-energy difference", "se_slater_total", "0", "parameter-free",
     {"uv": "ref_slater_total", "mc": "se_slater_total"}),
    ("cr", "Clementi-Raimondi, total-energy difference", "se_cr_total", "0", "parameter-free",
     {"mc": "se_cr_total"}),
    ("zexp", "Exact 1/Z series, completed square + rel/QED (zexp)", "fp_zexp", "0", "parameter-free",
     {"mc": "fp_zexp"}),
    ("gshm", "GSHM final (purely fitted)", "se_gshm", "32", "fitted",
     {"uv": "ref_gshm", "mc": "se_gshm"}),
    ("u35", "u35 (exact sigma1 + linear remainder)", "uni_u35", "35", "fitted",
     {"uv": "ref_u35", "pv": "ref_u35", "mc": "uni_u35"}),
    ("pocket", "pocket formula (8-parameter, Appendix A)", "uni_pocket", "8", "fitted",
     {"uv": "ref_pocket", "pv": "ref_pocket", "mc": "uni_pocket"}),
    ("pa_bound9", "pa_bound9 (9-parameter bounded variant)", "pa_bound9", "9", "fitted",
     {"pv": "pa_bound9", "mc": "pa_bound9"}),
    ("final", "pa_hier_rel = Screened Rydberg formula (final)", "uni", "33", "fitted",
     {"uv": "final", "pv": "pa_hier_rel", "mc": "uni"}),
    ("kregar", "Kregar/Di Rocco SHM + Dirac (our implementation)", "lit_kregar_dirac", "0", "parameter-free",
     {}),
    ("mendoza", "Mendoza 2011 NRSHM, primary (dE, low-j; published constants, not refit)",
     "bench_mendoza2011", "0 refit (published GA-fitted 19x19 sigma matrix)", "published fitted constants",
     {}),
    ("lsda", "LSDA Delta-SCF (own Kohn-Sham solver)", "fp_dft", "0", "parameter-free",
     {"mc": "fp_dft"}),
]
FINAL, SLATER = "final", "slater"

STRATA = [
    ("ALL", lambda r: True),
    ("neutral", lambda r: r["N"] == r["Z"]),
    ("H-like", lambda r: r["N"] == 1),
    ("experimental", lambda r: r["status"] == "experimental"),
    ("semi-empirical", lambda r: r["status"] == "semi-empirical"),
    ("theoretical", lambda r: r["status"] == "theoretical"),
    ("charge>=3", lambda r: r["Z"] - r["N"] >= 3),
    ("N<=10", lambda r: r["N"] <= 10),
    ("Z>=55", lambda r: r["Z"] >= 55),
]
SPLIT_ROWS = {      # test / validation rows of the fixed splits (for "subset, no fit" scoring)
    "V1": lambda r: 37 <= r["Z"] <= 54,
    "V2": lambda r: 45 <= r["Z"] <= 54,
    "S1": lambda r: r["Z"] % 5 == 0,
    "S3": lambda r: r["N"] % 6 == 0,
    "S2": lambda r: r["Z"] >= 55,
}
SEL = ["V1", "V2", "S1", "S3"]


def load_preds(stem):
    """Same reading rule as evaluate.evaluate: rows with a finite IE_pred_eV."""
    path = os.path.join(RES, stem + "_predictions.csv")
    preds = {}
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                v = float(r["IE_pred_eV"])
            except (ValueError, TypeError):
                continue
            if math.isfinite(v):
                preds[(int(r["Z"]), int(r["N"]))] = v
    return path, preds


def stats(preds, recs, pred=None):
    pairs = [(preds[(r["Z"], r["N"])], r["IE_eV"]) for r in recs
             if (r["Z"], r["N"]) in preds and (pred is None or pred(r))]
    s = EV._stats(pairs)
    return {"n": s["n"], "MAPE": s.get("MAPE_%"), "median": s.get("median_APE_%")}


def fmt(v):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "n/a"
    a = abs(v)
    if a >= 1e5:
        return f"{v:.2e}"
    if a >= 100:
        return f"{v:.0f}"
    if a == 0:
        return "0"
    return f"{v:#.3g}"


# ----------------------------------------------------------------------------- held-out copies
def _load_heldout_sources():
    with open(os.path.join(RES, "uni_validation.json"), encoding="utf-8") as f:
        uv = json.load(f)["models"]
    with open(os.path.join(RES, "pa_validation.json"), encoding="utf-8") as f:
        pv = json.load(f)["candidates"]
    with open(os.path.join(RES, "model_comparison.csv"), newline="", encoding="utf-8") as f:
        mc = {r["file"]: r for r in csv.DictReader(f)}
    return uv, pv, mc


def _num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def heldout_from(src, key, uv, pv, mc):
    """Normalised held-out record from one source file."""
    if src in ("uv", "pv"):
        d = (uv if src == "uv" else pv)[key]
        out = {sp: {"n": d[sp]["n"], "MAPE": d[sp]["MAPE"], "median": d[sp].get("median"),
                    "neutral": d[sp].get("neutral_MAPE"), "n_neutral": d[sp].get("n_neutral")}
               for sp in SEL + ["S2"]}
        out["selection"] = d["selection_score"]
        out["source"] = {"uv": "results/uni_validation.json", "pv": "results/pa_validation.json"}[src] \
            + f" [{key}]"
        return out
    r = mc[key]
    out = {sp: {"n": None, "MAPE": _num(r[f"{sp}_MAPE_%"]), "median": None, "neutral": None, "n_neutral": None}
           for sp in SEL + ["S2"]}
    out["S2"]["median"] = _num(r["S2_blind_median_%"])
    out["S2"]["neutral"] = _num(r["S2_blind_neutral_MAPE_%"])
    out["selection"] = _num(r["selection_score_%"])
    out["source"] = f"results/model_comparison.csv [{key}]"
    return out


def subset_noFit(preds, recs):
    """Split columns for a model that is not refit: MAPE on the split's test rows (own coverage n)."""
    out = {}
    for sp, f in SPLIT_ROWS.items():
        s = stats(preds, recs, f)
        neu = stats(preds, recs, lambda r, f=f: f(r) and r["N"] == r["Z"])
        out[sp] = {"n": s["n"], "MAPE": s["MAPE"], "median": s["median"],
                   "neutral": neu["MAPE"], "n_neutral": neu["n"]}
    full = all(out[sp]["n"] == sum(1 for r in recs if SPLIT_ROWS[sp](r)) for sp in SEL)
    out["selection"] = sum(out[sp]["MAPE"] for sp in SEL) / 4 if full else None
    out["source"] = "computed here: MAPE on the split's test rows, no fit (nothing held out)"
    return out


# ----------------------------------------------------------------------------- build
def build():
    recs = load_records()
    nrows = len(recs)
    uv, pv, mc = _load_heldout_sources()
    preds = {}
    files = {}
    for key, label, stem, npar, kind, srcs in MODELS:
        path, p = load_preds(stem)
        preds[key], files[key] = p, os.path.relpath(path, _ROOT).replace("\\", "/")
        ev = EV.evaluate(path)          # consistency with the shared scorer, exactly as evaluate.py runs
        mine = stats(p, recs)
        assert ev["ALL"]["n"] == mine["n"] and abs(ev["ALL"]["MAPE_%"] - mine["MAPE"]) < 1e-12, key

    covered_M = set(preds["mendoza"])
    covered_D = set(preds["lsda"])
    blocks = {
        "in_sample_all_rows": recs,
        "common_mendoza_5011": [r for r in recs if (r["Z"], r["N"]) in covered_M],
        "common_dft_207": [r for r in recs if (r["Z"], r["N"]) in covered_D],
    }
    matrix = {}
    for bname, brecs in blocks.items():
        matrix[bname] = {"rows_in_block": len(brecs), "models": {}}
        for key, *_ in MODELS:
            matrix[bname]["models"][key] = {sname: stats(preds[key], brecs, f) for sname, f in STRATA}

    heldout, checks = {}, []
    for key, label, stem, npar, kind, srcs in MODELS:
        recs_ho = None
        order = [s for s in ("uv", "pv", "mc") if s in srcs]
        if order:
            recs_ho = heldout_from(order[0], srcs[order[0]], uv, pv, mc)
            for other in order[1:]:          # cross-check the copies against each other
                o = heldout_from(other, srcs[other], uv, pv, mc)
                for sp in SEL + ["S2"]:
                    if recs_ho[sp]["MAPE"] is not None and o[sp]["MAPE"] is not None:
                        d = abs(recs_ho[sp]["MAPE"] - o[sp]["MAPE"])
                        checks.append({"model": key, "split": sp, "a": recs_ho["source"], "b": o["source"],
                                       "abs_diff": d})
        noFit = subset_noFit(preds[key], recs)
        if kind != "fitted":
            if recs_ho is not None:          # parameter-free: copied value must equal the subset value
                for sp in SEL + ["S2"]:
                    if recs_ho[sp]["MAPE"] is not None and noFit[sp]["n"]:
                        checks.append({"model": key, "split": sp, "a": recs_ho["source"],
                                       "b": "subset recomputed here", "abs_diff":
                                       abs(recs_ho[sp]["MAPE"] - noFit[sp]["MAPE"])})
                for sp in SEL + ["S2"]:          # fill n / median / neutral the copy lacks
                    for fld in ("n", "n_neutral"):
                        if recs_ho[sp][fld] is None:
                            recs_ho[sp][fld] = noFit[sp][fld]
                    for fld in ("median", "neutral"):
                        if recs_ho[sp][fld] is None and noFit[sp]["n"]:
                            recs_ho[sp][fld] = noFit[sp][fld]
                recs_ho["source"] += "; n, medians and neutral cells filled from the same rows (no fit)"
                if key == "lsda":
                    recs_ho["selection"] = None          # 207 rows only: not comparable (report.py shows n/a)
            else:
                recs_ho = noFit
        else:
            for sp in SEL + ["S2"]:
                if recs_ho[sp]["n"] is None:
                    recs_ho[sp]["n"] = sum(1 for r in recs if SPLIT_ROWS[sp](r))
                if recs_ho[sp]["n_neutral"] is None:
                    recs_ho[sp]["n_neutral"] = sum(1 for r in recs if SPLIT_ROWS[sp](r) and r["N"] == r["Z"])
        heldout[key] = recs_ho

    # Slater ratios (the "N times worse" statements)
    ratios = {}
    for other in (FINAL, "pa_bound9", "pocket"):
        a, b = matrix["in_sample_all_rows"]["models"][SLATER], matrix["in_sample_all_rows"]["models"][other]
        ha, hb = heldout[SLATER], heldout[other]
        ratios[other] = {
            "all-data MAPE": a["ALL"]["MAPE"] / b["ALL"]["MAPE"],
            "all-data median APE": a["ALL"]["median"] / b["ALL"]["median"],
            "neutral first IE MAPE (all-data fit)": a["neutral"]["MAPE"] / b["neutral"]["MAPE"],
            "experimental rows MAPE (all-data fit)": a["experimental"]["MAPE"] / b["experimental"]["MAPE"],
            "charge>=3 MAPE (all-data fit)": a["charge>=3"]["MAPE"] / b["charge>=3"]["MAPE"],
            "blind S2 MAPE": ha["S2"]["MAPE"] / hb["S2"]["MAPE"],
            "blind S2 median APE": ha["S2"]["median"] / hb["S2"]["median"],
            "blind S2 neutral MAPE": ha["S2"]["neutral"] / hb["S2"]["neutral"],
        }

    # row-set pairs that print different numbers for "the same" subset
    allb, mb = blocks["in_sample_all_rows"], blocks["common_mendoza_5011"]
    notM = [r for r in recs if (r["Z"], r["N"]) not in covered_M]
    pairs = {}
    for key in (FINAL, "pa_bound9", "pocket", SLATER):
        pairs[key] = {}
        for sname, f in STRATA:
            pairs[key][sname] = {"all_rows": stats(preds[key], allb, f),
                                 "mendoza_covered": stats(preds[key], mb, f),
                                 "not_mendoza_covered": stats(preds[key], notM, f)}

    meta = {
        "rows_in_dataset": nrows,
        "scorer": "evaluate._stats (shared scorer); finite predictions only, as evaluate.evaluate",
        "models": [{"key": k, "label": lab, "file": files[k], "fitted_params": npar, "kind": kind,
                    "coverage": len(preds[k])} for k, lab, stem, npar, kind, srcs in MODELS],
        "strata": [s for s, _ in STRATA],
        "notes": [
            "In-sample block: fitted models are their all-data fits; Z>=55 there is NOT the blind test.",
            "Held-out block: fitted models copied from the validation files (refit per split). For parameter-free "
            "and not-refit models the split cells are MAPEs on the same test rows; nothing is held out.",
            "Mendoza 2011 constants were GA-fitted by its authors to NIST+FAC data for He- to Eu-like sequences "
            "(with an extension to Z = 68-92), excluding neutral atoms and first ions: its split cells are subsets, "
            "not blind tests.",
            "Mendoza 2011 constants (models/benchmarks/mendoza2011/constants.json): all 361 entries were compared "
            "on 2026-10-06 with the text layer of the deposited article (https://oa.upm.es/11165/2/"
            "INVE_MEM_2011_102088.pdf, Tables 1-2): 0 mismatches (see docs/review/claims_and_matrix.md).",
            "LSDA dSCF covers 207 rows (all ions Z<=18 + neutral Z=19-54); no Z>=55 rows; no selection score.",
            "The pocket formula's V1 (4.4e6, 'fit diverges') is optimizer-path dependent: it converges to 4.35 % "
            "on the Python 3.11 / NumPy 2.4 / SciPy 1.17 stack (docs/review/reproducibility_py311.md).",
        ],
    }
    return {"meta": meta, "matrix": matrix, "heldout": heldout, "heldout_crosschecks": checks,
            "slater_ratios": ratios, "rowset_pairs": pairs}


# ----------------------------------------------------------------------------- writers
def write_csv(res):
    rows = []
    for bname, b in res["matrix"].items():
        for key, cells in b["models"].items():
            for sname, c in cells.items():
                rows.append([bname, key, sname, c["n"], c["MAPE"], c["median"], ""])
    for key, h in res["heldout"].items():
        for sp in SEL + ["S2"]:
            c = h[sp]
            rows.append(["heldout", key, sp, c["n"], c["MAPE"], c["median"], h["source"]])
            if c["neutral"] is not None:
                rows.append(["heldout", key, sp + "_neutral", c["n_neutral"], c["neutral"], "", h["source"]])
        rows.append(["heldout", key, "selection_score", "", h["selection"], "", h["source"]])
    for key, r in res["slater_ratios"].items():
        for q, v in r.items():
            rows.append(["slater_ratio", f"slater/{key}", q, "", v, "", ""])
    with open(OUT + ".csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["block", "model", "subset", "n", "MAPE_%", "median_APE_%", "source"])
        for r in rows:
            w.writerow([("" if x is None else (repr(x) if isinstance(x, float) else x)) for x in r])


def write_md(res):
    labels = {m["key"]: m for m in res["meta"]["models"]}
    L = ["# Benchmark matrix (all models, one scorer)", "",
         "Generated by `py -3.11 tools/benchmark_matrix.py` from the prediction CSVs in `results/` and the "
         "validation files; no model is refit. MAPE in %; every cell shows **MAPE [n]**, where n is the number "
         "of rows that model actually predicts in that cell. Fitted models: all-data fits in the in-sample "
         "blocks. `Z>=55` in the in-sample blocks is NOT the blind test; the blind S2 is in the held-out block.", ""]
    L += ["## Models", "", "| key | model | file | fitted params | coverage |", "|---|---|---|---|---|"]
    for m in res["meta"]["models"]:
        L.append(f"| {m['key']} | {m['label']} | `{m['file']}` | {m['fitted_params']} | {m['coverage']} |")
    titles = {"in_sample_all_rows": "In-sample, all 5847 rows (each model on the rows it covers)",
              "common_mendoza_5011": "Common subset: the 5011 rows covered by Mendoza 2011",
              "common_dft_207": "Common subset: the 207 LSDA Delta-SCF rows (Z<=18 all ions, Z=19-54 neutral)"}
    for bname, b in res["matrix"].items():
        L += ["", f"## {titles[bname]}", "", "MAPE [n]", "",
              "| model | " + " | ".join(res["meta"]["strata"]) + " |",
              "|---" * (len(res["meta"]["strata"]) + 1) + "|"]
        for key, cells in b["models"].items():
            L.append(f"| {key} | " + " | ".join(
                (f"{fmt(c['MAPE'])} [{c['n']}]" if c["n"] else "– [0]") for c in cells.values()) + " |")
        L += ["", "Median APE [n] (same cells)", "",
              "| model | " + " | ".join(res["meta"]["strata"]) + " |",
              "|---" * (len(res["meta"]["strata"]) + 1) + "|"]
        for key, cells in b["models"].items():
            L.append(f"| {key} | " + " | ".join(
                (f"{fmt(c['median'])} [{c['n']}]" if c["n"] else "– [0]") for c in cells.values()) + " |")
    L += ["", "## Held-out (fitted models refit per split; copied from the validation files)", "",
          "V1 fit Z<=36 / validate 37-54; V2 fit Z<=44 / validate 45-54; S1 test Z%5==0; S3 test N%6==0; "
          "selection = mean(V1, V2, S1, S3); S2 fit Z<=54 / test Z>=55 (blind for fitted models). "
          "Cells: MAPE [n]. For parameter-free / not-refit models the cells are MAPEs on the same test rows "
          "(nothing is held out).", "",
          "| model | V1 | V2 | S1 | S3 | selection | S2 | S2 median | S2 neutral | source |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for key, h in res["heldout"].items():
        cells = [f"{fmt(h[sp]['MAPE'])} [{h[sp]['n']}]" for sp in SEL]
        s2 = h["S2"]
        L.append(f"| {key} | " + " | ".join(cells) + f" | {fmt(h['selection'])} | "
                 f"{fmt(s2['MAPE'])} [{s2['n']}] | {fmt(s2['median'])} | "
                 f"{fmt(s2['neutral'])} [{s2['n_neutral']}] | {h['source']} |")
    mx = max((c["abs_diff"] for c in res["heldout_crosschecks"]), default=0.0)
    L += ["", f"Cross-checks between the copied sources (and against the no-fit recomputation for parameter-free "
          f"models): {len(res['heldout_crosschecks'])} comparisons, max |difference| = {mx:.2e} MAPE points.", ""]
    L += ["## Slater / model ratios (for any 'N times' statement)", "",
          "| quantity | Slater / final (33 p) | Slater / pa_bound9 (9 p) | Slater / pocket (8 p) |",
          "|---|---|---|---|"]
    qs = list(next(iter(res["slater_ratios"].values())))
    for q in qs:
        L.append(f"| {q} | " + " | ".join(f"{res['slater_ratios'][k][q]:.2f}" for k in (FINAL, "pa_bound9", "pocket"))
                 + " |")
    L += ["", "## Same subset, different row sets (why two documents print different numbers)", "",
          "| model | subset | all rows | Mendoza-covered rows | rows Mendoza does not cover |", "|---|---|---|---|---|"]
    for key, d in res["rowset_pairs"].items():
        for sname in ("ALL", "neutral", "experimental", "Z>=55", "charge>=3"):
            c = d[sname]
            L.append(f"| {key} | {sname} | " + " | ".join(
                f"{fmt(c[x]['MAPE'])} [{c[x]['n']}]" if c[x]["n"] else "– [0]"
                for x in ("all_rows", "mendoza_covered", "not_mendoza_covered")) + " |")
    L += ["", "## Notes", ""] + [f"- {n}" for n in res["meta"]["notes"]]
    with open(OUT + ".md", "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(L) + "\n")


def write_json(res):
    with open(OUT + ".json", "w", encoding="utf-8", newline="\n") as f:
        json.dump(res, f, indent=1)


# ----------------------------------------------------------------------------- optional refit check
def blind_check():
    """In-memory refit of pa_hier_rel on Z<=54 (the frozen spec), printed only. Nothing is written:
    the sigma1 cache writer is disabled for this process."""
    import numpy as np
    for p in (os.path.join(_ROOT, "models", "unified"), os.path.join(_ROOT, "models", "first_principles"),
              os.path.join(_ROOT, "models", "semi_empirical")):
        if p not in sys.path:
            sys.path.insert(0, p)
    import abinitio
    abinitio.save_cache = lambda: None          # never write models/unified/cache/sigma1.json
    import umodel as U
    from push_a_import import PA
    a = U.dataset()
    y, Z = a["y"], a["Z"].astype(int)
    N = (a["Z"] - a["Zq"] + 1).astype(int)
    recs = load_records()
    spec = PA.default_spec(grouping="pocket", rel="fs", hier=True, ridge=1e-4)   # = push_a/validate.py
    s, th = PA.fit(Z <= 54, dict(spec))
    p = PA.evaluate(th, s, a)
    ape = np.abs(p - y) / y * 100
    te = Z >= 55
    print(f"S2 refit: n_params {len(th)}  MAPE {ape[te].mean():.4f}  median {np.median(ape[te]):.4f}  "
          f"neutral {ape[te & (Z == N)].mean():.4f} (n={int((te & (Z == N)).sum())})  "
          f"ions Z>N {ape[te & (Z > N)].mean():.4f} (n={int((te & (Z > N)).sum())})")
    _, pm = load_preds("bench_mendoza2011")
    inM = np.array([(int(z), int(n)) in pm for z, n in zip(Z, N)])
    print(f"S2 refit on the Mendoza-covered Z>=55 rows: MAPE {ape[te & inM].mean():.4f} (n={int((te & inM).sum())})")
    idx = {(int(z), int(n)): i for i, (z, n) in enumerate(zip(Z, N))}
    for sym, z in (("Lu", 71), ("Tl", 81), ("Pb", 82), ("Rn", 86), ("Lr", 103)):
        i = idx[(z, z)]
        print(f"S2 refit {sym} I: {p[i]:.4f} eV vs NIST {y[i]:.4f}")
    # third IEs with f removal (stage 3: N = Z - 2)
    fr = [i for i, r in enumerate(recs) if te[i] and r["N"] == r["Z"] - 2 and r["removed"] and r["removed"][1] == 3]
    rat = np.array([p[i] / y[i] for i in fr])
    worst = max(fr, key=lambda i: p[i] / y[i])
    print(f"S2 refit third IEs with f removal: n={len(fr)} ratio mean {rat.mean():.3f} median {np.median(rat):.3f} "
          f"range {rat.min():.3f}-{rat.max():.3f}; worst Z={recs[worst]['Z']} {p[worst]:.2f} vs {y[worst]:.2f} eV")
    # all-data (frozen) parameters
    with open(os.path.join(RES, "pa_params.json"), encoding="utf-8") as f:
        c = json.load(f)["candidates"]["pa_hier_rel"]
    fz = (c["spec"], np.array(c["values"], float))
    for sym, z in (("Pb", 82), ("Rn", 86), ("Lr", 103), ("Cn", 112), ("Fl", 114), ("Og", 118), ("E120", 120)):
        print(f"all-data fit {sym} I: {float(PA.predict_many([(z, z)], params=fz)[0]):.4f} eV")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blind-check", action="store_true", help="also print the in-memory Z<=54 refit check")
    a = ap.parse_args()
    res = build()
    write_csv(res)
    write_md(res)
    write_json(res)
    print("wrote", OUT + ".{csv,md,json}")
    if a.blind_check:
        blind_check()


if __name__ == "__main__":
    main()
