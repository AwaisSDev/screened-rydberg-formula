"""H-like (N = 1) accuracy in layers; writes results/fp_hlike_*_predictions.csv + metrics."""
import csv, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from zexp import ROOT
sys.path.insert(0, ROOT)
import evaluate
import relativity as rel
from common.atomdata import load_records, HARTREE_EV

LAYERS = {
    "hlike_a_dirac": lambda Z: rel.hydrogenic_binding(Z, recoil=False, fns=False, qed=False),
    "hlike_b_dirac_recoil_fns": lambda Z: rel.hydrogenic_binding(Z, qed=False),
    "hlike_c_plus_qed_table": lambda Z: rel.hydrogenic_binding(Z),
    "hlike_d_plus_qed_closedform": lambda Z: rel.hydrogenic_binding_closed(Z),
}
out = {}
for name, fn in LAYERS.items():
    path = os.path.join(ROOT, "results", f"fp_{name}_predictions.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["Z", "N", "IE_pred_eV"])
        for r in load_records():
            if r["N"] == 1:
                w.writerow([r["Z"], 1, repr(fn(r["Z"]) * HARTREE_EV)])
    m = evaluate.evaluate(path)
    json.dump(m, open(path.replace("_predictions.csv", "_metrics.json"), "w"), indent=1)
    h = m["hydrogen_like"]
    out[name] = {"MAPE_%": h["MAPE_%"], "median_%": h["median_APE_%"], "max_%": h["max_APE_%"]}
    print(name, out[name])
for Z in (1, 10, 30, 50, 70, 92, 110):
    print(Z, "F_SE table", rel.F_SE(Z, 1), "closed", round(rel.F_SE_closed(Z, 1), 4))
json.dump(out, open(os.path.join(ROOT, "results", "fp_hlike_layers_summary.json"), "w"), indent=1)
