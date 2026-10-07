"""Insert the validation and head-to-head numbers into docs/literature.md (placeholders {{VALIDATION}}, {{HEAD2HEAD}})."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RES = os.path.join(ROOT, "results")
DOC = os.path.join(ROOT, "docs", "literature.md")


def main():
    v = json.load(open(os.path.join(RES, "lit_kregar_validation.json")))
    t1 = v["table1_total_energy_Ht"]
    t1rows = " ".join(f"{el} {d['ref_Ht']:.2f}/{d['pauli']:.2f};" for el, d in t1.items())
    t1dev = max(abs(d["pauli"] - d["ref_Ht"]) / d["ref_Ht"] * 100 for d in t1.values())
    t2 = v["table2_Ar_isonuclear"]
    t3 = v["table3_Ar_isoelectronic"]
    val = (
        f"- Diagonal Z → ∞ parameters k_ii (Table 6 of [Pomarico2005]): reproduced to all printed digits, e.g. "
        f"1s 0.3125, 2s 0.3008, 2p 0.3492, 3d 0.3765. Off-diagonal g/f: mean |Δ| = {v['table6_mean_abs_dev']:.3f}, "
        f"max {v['table6_max_abs_dev']:.3f}. Our values follow the definitions exactly; the paper uses fitted closed forms "
        f"and an exchange apportioning that the 2005 paper does not fully specify.\n"
        f"- Total configuration energies (Table 1, hartree, published/ours with the Pauli correction): {t1rows} "
        f"max deviation {t1dev:.2f} %.\n"
        f"- Ar isonuclear IEs (Table 2): mean deviation from the published model values "
        f"{v['table2_mean_rel_dev_%']:.1f} %. Inner-shell IEs agree to <1 %; valence IEs are 3–5 eV higher in our "
        f"implementation (Ar I 18.96 vs printed 14.72 eV; experiment 15.76 eV). Ar-like sequence (Table 3): mean "
        f"deviation {v['table3_mean_rel_dev_%']:.1f} %. So our reproduction matches the published model where screening "
        f"is dominated by k_ii and the Z → ∞ limit. It is less faithful for valence shells of near-neutral ions, where the "
        f"off-diagonal g/f (fitted closed forms in the original) matter. Treat the near-neutral numbers in §2.2 as "
        f"'the Kregar/Di Rocco model as defined', not as the authors' code.")
    rows = json.load(open(os.path.join(RES, "lit_comparison.json")))
    md = open(os.path.join(RES, "lit_comparison.md"), encoding="utf-8").read().split("\n", 3)[3]
    best_free = min((r for r in rows if r["kind"] == "free" and "Kregar" in r["model"]), key=lambda r: r["ALL"])
    sl = next(r for r in rows if r["model"].startswith("Slater"))
    pa = next((r for r in rows if "pa_hier_rel" in r["file"]), None)
    pk = next((r for r in rows if "pocket" in r["file"]), None)

    def f(x):
        return "n/a" if x is None else (f"{x:.3g}" if x < 1000 else f"{x:.0f}")
    h2h = (
        "MAPE in %, from evaluate.py on the identical 5847 NIST rows (`models/literature/compare.py`, "
        "`results/lit_comparison.{md,json}`). "
        "Parameter-free models: S1/S2/S3 = MAPE on the test rows. Fitted project models: refit held-out values from "
        "`docs/unified.md`. The final project model may still change (another agent is integrating it); "
        "`pa_hier_rel` is the current final candidate and will be re-scored when frozen.\n\n" + md + "\n\n"
        "**Reading the table.**\n"
        f"- Best literature variant: {best_free['model']}. It reaches {f(best_free['ALL'])} % all-data MAPE "
        f"(median {f(best_free['median_APE_%'])} %), {f(best_free['first_IE_neutral_atoms'])} % on neutral first IEs "
        f"and {f(best_free['hydrogen_like'])} % on H-like ions. Its blind-equivalent S2 (Z ≥ 55) score is "
        f"{f(best_free['S2_test'])} %.\n"
        f"- Slater's rules on the same rows: {f(sl['ALL'])} % all-data, {f(sl['first_IE_neutral_atoms'])} % neutral.\n"
        + (f"- Project pocket formula (8 parameters): {f(pk['ALL'])} % all-data, blind S2 {f(pk['S2_test'])} %.\n" if pk else "")
        + (f"- Project final candidate `pa_hier_rel` (33 parameters): {f(pa['ALL'])} % all-data, "
           f"{f(pa['first_IE_neutral_atoms'])} % neutral, blind S2 {f(pa['S2_test'])} %.\n" if pa else "")
        + "- A parameter-free published model cannot be compared with a fitted one on parameter count alone. The fair "
          "statements are:\n"
          "  1. on highly charged ions (N/Z small), the parameter-free SHM and the project's exact-σ₁ models both work, "
          "because the Z → ∞ limit is built in;\n"
          "  2. on near-neutral atoms, every parameter-free hydrogenic model degrades. The fitted remainder is what "
          "makes the project's formula useful there, at the cost of 8–35 parameters;\n"
          "  3. the project's blind S2 numbers come from extrapolating a fit, while the literature model has nothing "
          "to extrapolate. The comparison tests whether fitting buys more than it risks. Against the literature "
          f"model's {f(best_free['S2_test'])} % on the same Z ≥ 55 rows: "
        + "; ".join(f"{name} ({f(val)} %) is {'better' if val < best_free['S2_test'] else 'worse'}"
                    for name, val in [("`pa_hier_rel`", 6.60), ("pocket", 11.3), ("`u35`", 57.2), ("GSHM", 170.0)])
        + ".\n")
    s = open(DOC, encoding="utf-8").read()
    s = s.replace("{{VALIDATION}}", val).replace("{{HEAD2HEAD}}", h2h)
    open(DOC, "w", encoding="utf-8", newline="\n").write(s)
    print(val)
    print(h2h)


if __name__ == "__main__":
    main()
