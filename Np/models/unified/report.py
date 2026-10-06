"""Model comparison table (results/model_comparison.{csv,md}) and uni_ figures for the FINAL model. ~30 s.

    py -3.13 models/unified/report.py          (run models/unified/validate_blind.py first)

Held-out columns use the pre-registered protocol: V1 (fit Z<=36 / val 37-54), V2 (fit Z<=44 / val 45-54),
S1 (test Z%5==0), S3 (test N%6==0), selection score = mean(V1, V2, S1, S3); S2 (fit Z<=54 / test Z>=55) is BLIND.
Sources: parameter-free models -> MAPE on the test/validation rows directly; refitted models ->
results/uni_validation.json (validate_blind.py), results/pa_validation.json (Push A), results/pb_validation.json
(Push B), results/uni_legacy_validation.json (S1/S2/S3 of uni_gshm_qed), Track A's reported fits.
"""
import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, _ROOT)
import evaluate  # noqa: E402
from common.atomdata import load_records  # noqa: E402

RES = os.path.join(_ROOT, "results")
FIG = os.path.join(RES, "figures")


def _j(name):
    p = os.path.join(RES, name)
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else {}


UV = _j("uni_validation.json").get("models", {})
PAV = _j("pa_validation.json").get("candidates", {})
PBV = _j("pb_validation.json").get("candidates", {})
LEG = _j("uni_legacy_validation.json").get("splits", {})
SEV = _j("se_validation.json").get("splits", {})
SUB = {"V1": lambda r: 37 <= r["Z"] <= 54, "V2": lambda r: 45 <= r["Z"] <= 54,
       "S1": lambda r: r["Z"] % 5 == 0, "S3": lambda r: r["N"] % 6 == 0, "S2": lambda r: r["Z"] >= 55}
SEL = ["V1", "V2", "S1", "S3"]
FP_REPORTED = {"fp_zexp_univ": {"S1": 38.6, "S2": 65.3, "S3": 42.9},
               "fp_zexp3_seq": {"S1": 26.6, "S2": 75.7, "S3": None}}

# (label, predictions-file prefix, n_params, type, held-out source)
MODELS = [
    ("Bohr Ry Z^2/n^2", "se_bohr", 0, "first-principles (baseline)", "zero"),
    ("Slater 1e", "se_slater_1e", 0, "empirical rules (baseline, 0 params)", "zero"),
    ("Slater total-energy diff.", "se_slater_total", 0, "empirical rules (baseline, 0 params)", "zero"),
    ("Clementi-Raimondi 1e", "se_cr_1e", 0, "empirical/HF (baseline, 0 params)", "zero"),
    ("Clementi-Raimondi total", "se_cr_total", 0, "empirical/HF (baseline, 0 params)", "zero"),
    ("A: 1/Z two terms + rel/QED", "fp_zexp2_rel", 0, "first-principles", "zero"),
    ("A: zexp headline (exact O(Z^2,Z) + heuristic O(1) + rel)", "fp_zexp", 0, "first-principles", "zero"),
    ("A: zexp_sq_rel_p4", "fp_zexp_sq_rel_p4", 0, "first-principles", "zero"),
    ("A: universal sigma(N/Z)", "fp_zexp_univ", 12, "hybrid", ("fp", "fp_zexp_univ")),
    ("A: per-sequence dE2 (in-sample)", "fp_zexp3_seq", 409, "hybrid (per configuration)", ("fp", "fp_zexp3_seq")),
    ("A: LSDA Delta-SCF (207 rows)", "fp_dft", 0, "first-principles (DFT)", "zero"),
    ("A: Slater transition state (207 rows)", "fp_dft_ts", 0, "first-principles (DFT)", "zero"),
    ("B: GSHM core (blind)", "se_gshm_core", 30, "fitted", ("uv", "ref_gshm_core")),
    ("B: GSHM final (blind)", "se_gshm", 32, "fitted", ("uv", "ref_gshm")),
    ("U: GSHM + A recoil/QED/FNS", "uni_gshm_qed", 32, "fitted + first-principles layer", ("leg", "uni_gshm_qed")),
    ("U: pocket formula", "uni_pocket", 8, "fitted (hand-calculable)", ("uv", "ref_pocket")),
    ("U: exact sigma1 + remainder (29p)", "uni_u29", 29, "hybrid", ("uv", "ref_u29")),
    ("U: exact sigma1 + remainder (35p) (previous final)", "uni_u35", 35, "hybrid", ("uv", "ref_u35")),
]
for nm in ["pa_hier", "pa_bound9", "pa_bound14_relfit", "pa_bound9_fs0", "pa_nobound14_relfit"]:
    if nm in PAV:
        MODELS.append((f"Push A: {nm}", nm, PAV[nm]["all"]["n_params"], "hybrid", ("pa", nm)))
for nm in ["pb_clip_pos", "pb_exp_pos", "pb_clip", "pb_logistic_pos"]:
    if nm in PBV:
        MODELS.append((f"Push B: {nm}" + (" (Push B choice)" if nm == "pb_clip_pos" else ""), f"pb_{nm[3:]}",
                       PBV[nm].get("n_params", 30), "hybrid", ("pb", nm)))
MODELS.append(("FINAL: pa_hier_rel (Screened Rydberg formula)", "uni", 33, "hybrid", ("uv", "final")))


def heldout(src, s):
    kind, key = src
    if kind == "fp":
        return FP_REPORTED[key].get(s) if s in ("S1", "S2", "S3") else None
    if kind == "leg":
        v = LEG.get(f"{key}_{s}")
        return v["MAPE"] if v else None
    d = {"uv": UV, "pa": PAV, "pb": PBV}[kind].get(key, {})
    v = d.get(s)
    return v["MAPE"] if isinstance(v, dict) else None


rows = []
for label, pre, npar, kind, ho in MODELS:
    path = os.path.join(RES, f"{pre}_predictions.csv")
    mp = os.path.join(RES, f"{pre}_metrics.json")
    m = json.load(open(mp, encoding="utf-8")) if os.path.exists(mp) else evaluate.evaluate(path)

    def g(key, f="MAPE_%"):
        return m.get(key, {}).get(f) if m.get(key, {}).get("n", 0) else None
    row = {"model": label, "file": pre, "type": kind, "n_params": npar, "n_rows": m["ALL"]["n"],
           "MAPE_%": g("ALL"), "median_APE_%": g("ALL", "median_APE_%"),
           "neutral_first_IE_MAPE_%": g("first_IE_neutral_atoms"), "H_like_MAPE_%": g("hydrogen_like")}
    for s in SEL + ["S2"]:
        if ho == "zero":
            v = evaluate.evaluate(path, subset=SUB[s]).get("ALL", {}).get("MAPE_%")
        else:
            v = heldout(ho, s)
        row[f"{s}_MAPE_%"] = v
    vals = [row[f"{s}_MAPE_%"] for s in SEL]
    row["selection_score_%"] = float(np.mean(vals)) if all(v is not None for v in vals) else None
    if pre in ("fp_dft", "fp_dft_ts"):
        # Final-audit fix: DFT rows cover only 207 rows with Z<=54 (V1/V2 = 17/10 neutral rows, no fitting), so a
        # "selection score" for them is not comparable with the fitted models' scores. Not shown.
        row["selection_score_%"] = None
    row["S2_blind_median_%"] = None
    row["S2_blind_neutral_MAPE_%"] = None
    if ho != "zero" and ho[0] in ("uv", "pa", "pb"):
        d = {"uv": UV, "pa": PAV, "pb": PBV}[ho[0]].get(ho[1], {}).get("S2", {})
        row["S2_blind_median_%"], row["S2_blind_neutral_MAPE_%"] = d.get("median"), d.get("neutral_MAPE")
    elif ho == "zero":
        e = evaluate.evaluate(path, subset=SUB["S2"])
        row["S2_blind_median_%"] = e.get("ALL", {}).get("median_APE_%")
        row["S2_blind_neutral_MAPE_%"] = e.get("first_IE_neutral_atoms", {}).get("MAPE_%") if e.get("first_IE_neutral_atoms", {}).get("n") else None
    row["S2_posthoc_bound_MAPE_%"] = {"se_gshm_core": SEV.get("core_S2", {}).get("MAPE"),
                                      "se_gshm": SEV.get("final_S2", {}).get("MAPE"),
                                      "uni_gshm_qed": LEG.get("uni_gshm_qed_posthoc_S2", {}).get("MAPE")}.get(pre)
    rows.append(row)

cols = list(rows[0].keys())
with open(os.path.join(RES, "model_comparison.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)


def fmt(v):
    if v is None:
        return "n/a"
    if isinstance(v, str):
        return v
    if isinstance(v, int):
        return str(v)
    if abs(v) >= 1e4:
        return f"{v:.1e}"
    return f"{v:.4f}" if v < 0.1 else (f"{v:.2f}" if v < 100 else f"{v:.0f}")


hdr = ["model", "type", "n_params", "n_rows", "MAPE_%", "median_APE_%", "neutral_first_IE_MAPE_%", "H_like_MAPE_%",
       "V1_MAPE_%", "V2_MAPE_%", "S1_MAPE_%", "S3_MAPE_%", "selection_score_%", "S2_MAPE_%",
       "S2_blind_neutral_MAPE_%", "S2_posthoc_bound_MAPE_%"]
md = ["# Model comparison (all tracks, pre-registered protocol)", "",
      "Generated by `py -3.13 models/unified/report.py`. All-data columns: evaluate.py on the NIST table (5847 rows unless "
      "n_rows says otherwise; fitted models fitted on all rows). Held-out columns (test/validation MAPE after refitting on the "
      "training part): **V1** fit Z<=36 / validate 37-54; **V2** fit Z<=44 / validate 45-54; **S1** test Z%5==0; **S3** test "
      "N%6==0; **selection score** = mean(V1, V2, S1, S3), the only quantity used to choose the final model; **S2** fit Z<=54 / "
      "test Z>=55 is the BLIND extrapolation test, computed once per candidate after all choices were frozen. For parameter-free "
      "models the held-out columns are the MAPE on those rows (DFT rows: only the 207 computed rows, Z<=54). Track A's fitted "
      "variants: S1/S2/S3 as reported by Track A (V1/V2 not run). uni_gshm_qed: S1/S2/S3 from the earlier unified run "
      "(V1/V2 not run; its fit is identical to GSHM final, so expect GSHM final's values). GSHM rows are BLIND: the sigma_core in "
      "[0.8,1] bound that Track B chose after seeing its S2 result is removed; the last column shows S2 WITH that post-hoc bound, "
      "for reference only. The pocket formula diverges on V1 (unbounded Zeff), so its selection score is meaningless. "
      "DFT rows (LSDA Delta-SCF, Slater transition state): their V1/V2/S1/S3 cells are MAPEs on the few computed rows only "
      "(207 rows in all, Z<=54; V1 and V2 contain just 17 and 10 neutral rows) and nothing is fitted, so no selection score "
      "is shown for them (n/a): it would not be comparable with the fitted models' scores. "
      "Caveat: S1 and S3, as pre-registered, contain Z>=55 rows in train and test (913/1188 S1 and 703/928 S3 test rows), "
      "so heavy-atom interpolation accuracy did inform the selection score; only S2 was never used for any choice. With S1/S3 "
      "restricted to Z<=54 the winner is unchanged (results/uni_sensitivity_z54.json: pa_hier_rel 1.844 < pa_hier 1.989 < "
      "pa_bound9 2.094 < pb_clip_pos 2.118). "
      "Type: first-principles = no fitted parameters; fitted = all screening fitted to data; hybrid = exact Track A first-order "
      "screening sigma1 (and relativistic/QED layers) plus a fitted remainder.", "",
      "| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
for r in rows:
    cells = [fmt(r[h]) for h in hdr]
    if r["file"] == "uni":
        cells = [f"**{c}**" for c in cells]
    md.append("| " + " | ".join(cells) + " |")
fin = UV.get("final", {})
md += ["", "Final model = **pa_hier_rel** (Push A), the lowest pre-registered selection score among all push candidates "
       f"({fmt(fin.get('selection_score'))}; Push B's choice pb_clip_pos {fmt(PBV.get('pb_clip_pos', {}).get('selection_score'))}; "
       f"previous final u35 {fmt(UV.get('ref_u35', {}).get('selection_score'))}). Its blind S2 is {fmt(fin.get('S2', {}).get('MAPE'))}% "
       f"(median {fmt(fin.get('S2', {}).get('median'))}%, neutral first IE {fmt(fin.get('S2', {}).get('neutral_MAPE'))}%).",
       "", "Note: pb_exp_pos has a lower blind S2 (6.03%) than pa_hier_rel but a worse selection score (3.73); promoting it "
       "would be selection on S2, so it is not the final model."]
with open(os.path.join(RES, "model_comparison.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(md) + "\n")
print("\n".join(md))

# ------------------------------------------------------------------ figures (final model)
recs = load_records()
Z = np.array([r["Z"] for r in recs]); N = np.array([r["N"] for r in recs])
y = np.array([r["IE_eV"] for r in recs]); L = np.array([r["removed"][1] for r in recs])
pred = {}
for nm in ["uni", "uni_u35", "uni_pocket"]:
    with open(os.path.join(RES, f"{nm}_predictions.csv"), encoding="utf-8") as f:
        d = {(int(r["Z"]), int(r["N"])): float(r["IE_pred_eV"]) for r in csv.DictReader(f)}
    pred[nm] = np.array([d[(z, n)] for z, n in zip(Z, N)])
p = pred["uni"]
cl = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
TITLE = "Screened Rydberg formula (pa_hier_rel, 33 params)"

fig, ax = plt.subplots(figsize=(6, 6))
for l in range(4):
    k = L == l
    ax.loglog(y[k], p[k], ".", ms=2, color=cl[l], label="removed " + "spdf"[l])
ax.plot([3, 3e5], [3, 3e5], "k-", lw=0.7)
ax.set_xlabel("NIST IE (eV)"); ax.set_ylabel("predicted IE (eV)")
ax.set_title(TITLE + "\nall 5847 IEs, all-data fit", fontsize=10); ax.legend(markerscale=5)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "uni_parity.png"), dpi=150); plt.close(fig)

k = Z == N
o = np.argsort(Z[k])
fig, ax = plt.subplots(figsize=(10, 4.5))
ax.plot(Z[k][o], y[k][o], "k.-", lw=1, label="NIST")
ax.plot(Z[k][o], p[k][o], "o-", ms=3, lw=1, color="#d62728", label="final pa_hier_rel (33p)")
ax.plot(Z[k][o], pred["uni_u35"][k][o], "s-", ms=2, lw=0.7, color="#1f77b4", alpha=0.6, label="previous u35 (35p)")
ax.plot(Z[k][o], pred["uni_pocket"][k][o], "-", lw=0.7, color="#2ca02c", alpha=0.7, label="pocket (8p)")
ax.set_xlabel("Z"); ax.set_ylabel("first ionization energy (eV)"); ax.set_ylim(0, 30); ax.legend()
ax.set_title("Neutral-atom first IE (all-data fits)")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "uni_first_IE.png"), dpi=150); plt.close(fig)

fig, axs = plt.subplots(1, 4, figsize=(16, 4))
for ax, (sym, zz) in zip(axs, [("C", 6), ("Fe", 26), ("Xe", 54), ("U", 92)]):
    k = Z == zz
    st = zz - N[k] + 1
    o = np.argsort(st)
    ax.semilogy(st[o], y[k][o], "k.-", label="NIST")
    ax.semilogy(st[o], p[k][o], "o", mfc="none", color="#d62728", ms=4, label="final")
    ax.set_title(sym); ax.set_xlabel("ionization stage (Z - N + 1)")
axs[0].set_ylabel("IE (eV)"); axs[0].legend()
fig.suptitle(TITLE + ": successive IEs", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "uni_successive.png"), dpi=150); plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 4.5))
for l in range(4):
    k = L == l
    ax.plot(Z[k] + 0.25 * (l - 1.5), (p[k] - y[k]) / y[k] * 100, ".", ms=2.5, color=cl[l], label="removed " + "spdf"[l])
ax.axhline(0, color="k", lw=0.6); ax.set_ylim(-25, 25)
ax.set_xlabel("Z"); ax.set_ylabel("(pred - NIST)/NIST  (%)"); ax.legend(markerscale=5)
ax.set_title(TITLE + ": relative residuals, all-data fit (clipped at +-25%)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "uni_residuals.png"), dpi=150); plt.close(fig)
print("figures written")
