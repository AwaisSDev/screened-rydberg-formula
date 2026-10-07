"""Paper case-study table (O, Na, N, Ca, Fe, Pb x three models) from results/audit_components.json.

    py -3.11 tools/audit_components.py && py -3.11 tools/case_study_table.py   -> results/case_studies.md
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SYM = {8: "O", 11: "Na", 7: "N", 20: "Ca", 26: "Fe", 82: "Pb"}
LAB = {"pa_hier_rel": "final (33 p)", "pa_bound9": "pa_bound9 (9 p)", "pocket": "pocket (8 p)"}
G = ["same", "in", "core", "df", "out"]


def main():
    d = json.load(open(os.path.join(ROOT, "results", "audit_components.json"), encoding="utf-8"))
    L = ["| atom (removed) | model | σ₁ | ν_g (same, in, core, df, out) | T | h | D | Z_eff | Ry Z_eff²/n² (eV) "
         "| relativistic factor | Hund term (eV) | IE (eV) | NIST (eV) | error |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in d:
        nu = ", ".join(str(int(r["nu_groups"][g])) for g in G)
        head = f"| {SYM[r['Z']]} ({int(r['n'])}{r['l']}) | {LAB[r['model']]} |"
        tail = (f" {r['Zeff']:.4f} | {r['rydberg_term']:.4f} | {{rel}} | {r['hund_term']:+.4f} | **{r['IE']:.3f}** "
                f"| {r['IE_NIST']:.3f} | {r['err_pct']:+.2f} % |")
        if r["model"] == "pocket":
            mid = f" – | {nu} | – | – | {r['total_screening']:.4f}* |"
            L.append(head + mid + tail.format(rel=f"{r['sommerfeld_bracket']:.5f}"))
        else:
            mid = f" {r['sigma1']:.4f} | {nu} | {r['T']:.4f} | {r['h']:.4f} | {r['D_rem']:.4f} |"
            L.append(head + mid + tail.format(rel=f"{r['F_Dirac'] * r['rel_bracket']:.5f}"))
    L += ["", "*Pocket formula: total screening Σ s_g ν_g + t(N−1)/(Z_a+κ); it has no σ₁ and no bound. "
          "Relativistic factor: F_{n,j}(Z_eff)·R for the σ₁ models (R = 1 for pa_bound9), Sommerfeld bracket for the "
          "pocket formula. μ(Z) and QED/FNS (zero here) are applied as in eq. (2.2)."]
    out = os.path.join(ROOT, "results", "case_studies.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
