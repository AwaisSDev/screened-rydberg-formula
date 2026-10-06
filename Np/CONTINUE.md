# CONTINUE HERE: project handoff

## ▶ 2026-10-06 review pass (D:\ machine). READ THIS FIRST

**Environment.** This machine has no Python 3.13. Use `py -3.11` with
`export PYTHONPATH="D:/Chem-Research/Np/tools/numba_stub"`, a no-op numba stand-in with identical numerics. Frozen
hashes are intact; models/push_a/model.py matches once CRLF is normalised to LF. The refits reproduce the final model
to ≤0.005 percentage points. The unbounded references do not: the pocket formula's V1 converges to 4.35 % here
instead of diverging. Details: `docs/review/reproducibility_py311.md`.

**Task.** Work through a reviewer's 25-item checklist. Items 1–19 must be done before anyone else sees the paper;
20–21 are optional; 22–25 are the user's own actions. At the user's request, all agents were stopped at 16:08 to
save session limit. **Continue without spawning agents unless the user asks for them.**

| item | status | where |
|---|---|---|
| 1 component printout (O, Na, N, Ca, Fe, Pb) | ✅ | `tools/audit_components.py` → `results/audit_components.{md,json}`; asserts equality with the code |
| 2 hand recompute, paper only | ⚠ partial | `tools/hand_recompute.py`. The agent was stopped mid-diagnosis: the pocket formula matches the code to 6e-15 on all rows, but the paper-only pa_* versions disagreed for some Z (NaNs). Check whether the paper text is insufficient. |
| 3 Ca/Fe gaps, Pb σ₁ 76.2 | ✅ | `docs/review/cheat_sheet_reconciliation.md`. No model gives +0.24/+0.39 eV; the cheat sheet's Ca/Fe lines are not reproducible (true pa_bound9: Ca 4.32 eV, Fe 7.30 eV). Pb σ₁ = 63.72; "76.2" is total screening or Slater's σ. |
| 4 cheat sheet = pa_bound9; add its parameters to Table A1 | ✅ identified, ⏳ paper edit | same file; parameters in results/pa_params.json |
| 5 one grouping rule | ✅ | `docs/review/grouping_rule.md` (paper-ready text + Table A2); `tools/check_grouping_rule.py`: 0 mismatches on 5847 rows and 7021 Z ≤ 118 configurations |
| 6 rename models (pa_bound9 ≠ pocket) | ⏳ | canonical names: "Screened Rydberg formula (final, pa_hier_rel, 33 p)", "bounded 9-parameter variant (pa_bound9)", "pocket formula (8 p, uni_pocket)" |
| 7 remove "destroyed"/"best" | ⏳ | only "best Push B candidate" (paper l.265) is in the repo; "destroyed" is in the user's own material |
| 8 Mendoza numbers | ✅ sourced | constants transcribed from the user-supplied open-access PDF (oa.upm.es/11165); 6 printed tables reproduced (`models/benchmarks/mendoza2011/NOTES.md`). Keep the numbers, but with row sets (5011 covered rows). Paper §4.5 still wrongly says "not obtained". |
| 9 reconcile matrix vs paper | ✅ matrix, ⏳ paper | `results/benchmark_matrix.{md,csv,json}` (`tools/benchmark_matrix.py`, n in every cell); row-set issues in cheat_sheet_reconciliation.md §5 |
| 10 Slater ratio | ⏳ | mean: 11.8/1.87 = 6.3× (final), 11.8/2.90 = 4.1× (pa_bound9); blind S2: 11.6/6.60 = 1.8× |
| 11 authors, affiliation, acknowledgements | ⏳ user | |
| 12 NIST ASD version + access date | ⏳ | refs agent stopped mid-check; access date not recorded anywhere, user must supply |
| 13 verify references | ⏳ | refs-part1 (A–K) was ~77 % done when stopped, no report written; its transcript is in the session's `subagents/workflows/wf_28a7f90b-4d3/agent-a2b967d2350b4c32d.jsonl` |
| 14–16 σ₁-table novelty, ML search, superheavy values | ⏳ not started | |
| 17, 18 Zenodo DOI, AI statement | ⏳ user (after journal choice) | |
| 19 disclose failures | ⏳ | paper §4.3/§4.7 already list Lr, Og, Pb, Tl; add them to the abstract and conclusions |
| 20 Mendoza on the 5847 rows | ✅ (2026-10-05) | 5011/5847 covered: 2.82 % vs final 1.62 % on the same rows; Mendoza is the most accurate of the three on N ≤ 10 (0.40 vs 0.44) |
| 21 Rodrigues 2004 Dirac–Fock | ⏳ optional | |

**UPDATE 2026-10-06 16:30: the user's PHASE 1 (items 1–4 of their new plan) is DONE in `docs/paper_draft.md`, and
the PDF is rebuilt (`py -3.14 tools/build_paper_pdf.py`; the build script needs Python ≥ 3.12).**
- **Case studies.** New Table 6 (O, N, Na, Ca, Fe, Pb × final / pa_bound9 / pocket), from
  `tools/case_study_table.py`.
- **Hand-recompute gap closed.** The paper never stated four things: the ν_c classes, that σ₁ uses the frozen
  configuration, μ(Z), and the QED/FNS values. The fix:
  - `results/model_inputs.csv` (`tools/export_model_inputs.py`);
  - `tools/verify_from_inputs.py`, which uses no model code and reproduces both models on all 5847 rows to 7e-16;
  - §2.2 rewritten with the grouping rule, the removed-subshell rule, K_l, F_{n,j} (renamed from D_{n,j}), μ and QED.
- **Naming.** A naming table in §2.6; the pa_bound9 parameters added to Appendix A; Table A2 holds the class
  definitions.
- **Paper edits.**
  - Mendoza §4.5 with Table 4a (5011 rows, with caveats).
  - §3.4 reproducibility, including the pocket-formula V1 caveat (†).
  - The Slater ratios: 6.3× / 4.1× / 1.8×.
  - "best" changed to "lowest-scoring".
  - Failures added to the abstract, §5.2 and the conclusions; numbers from `tools/known_failures.py`.
- `docs/literature.md`: the Mendoza lines are fixed.

**2026-10-06 ~17:30: "do all the still open things" pass.**
- **Paper is in final form** (IEEE-style numbered citations and reference list, built by `tools/ieee_refs.py`;
  "Abstract—" / "Index Terms—" front matter). It has no draft banner and 0 TODO/VERIFY markers.
- **References.** The 14 flagged references were checked against Crossref (DOIs added). Layzer1967 is not in
  Crossref and was dropped; DiRocco1992 is labelled "as cited in Pomarico2005".
- **Other paper fixes.**
  - NIST ASD v5.12, accessed on or before 2026-10-05.
  - The NumPy/SciPy versions are stated as not recorded.
  - Code availability points to the GitHub URL.
  - Acknowledgements: None.
- **Phase 3 done** (except the web calculator):
  - `ionization.ie(Z, N, config)` returns all intermediates and warnings; CLI `--details --config`.
  - `README.md`; its examples come from `tools/readme_examples.py`.
  - `tests/test_paper_examples.py`: 7 tests pass.
- **Zenodo.** `CITATION.cff` at the git root. The user must link GitHub to Zenodo and make a release to get a DOI.
- **The journal** is still the user's choice, after human review.

**2026-10-06 ~17:00: reviewer fix list 1–9 applied.**
- Author: Muhammad Awais, independent researcher, mawais9171@gmail.com.
- §4.5 split into Kregar (4.5.1) and Mendoza (4.5.2).
- Table 6 verdict fixed.
- Pb cancellation added to §5.2.
- N worked example added.
- Abstract: 57.2 % / 170 %.
- Unverifiable claims softened (the 88 % claim, the σ₁-table novelty, ML priority).
- Cheat sheet: `docs/cheat_sheet_pa_bound9.md` (`tools/cheat_sheet_bound9.py`).
- PDF rebuilt.
- Remaining paper markers: the NIST version and date, NumPy/SciPy versions, the Zenodo DOI, acknowledgements, and
  the reference [VERIFY] tags.
- **The pre-registration commit is waiting for the user's OK. §3.3 already says "committed afterwards".**

**PHASE 2, round 1 done (16:35).** Four pre-registered variants were scored and all rejected (`models/v2/NOTES.md`);
v1 stays final. The paper discloses this in §3.3 (one sentence). §4.5 Table 4a is now a same-row comparison of the
33 p model, pa_bound9 and Mendoza (`tools/compare_mendoza.py`). Claim: "competitive", not "better". Mendoza has 331
non-zero published constants.
(Original note: PHASE 2 had not started at that point.) Its pre-registration is written in `docs/preregistration_v2.md`; the user must commit it
(baseline score in this environment 2.1473). Still open from the reviewer list: 12–16 (NIST version, references,
σ₁ novelty, ML search, superheavy values) and the user-only items.

**Earlier next steps (kept for reference):**
1. **Replace the Ca, Fe and Na case studies with real code output, labelled by model.** The user's earlier Ca/Fe case
   studies (+2.9 %, +4.8 %, the "+0.24/+0.39 eV" gaps, and "Hund/Spin correction = 8.28") were AI-written numbers,
   not model output. Never reuse them. Real values (`results/audit_components.md`), IE in eV with error vs NIST:

   | ion (NIST) | final pa_hier_rel (33 p) | pa_bound9 (9 p) | pocket (8 p) |
   |---|---|---|---|
   | Na I (5.139) | 6.203 (+20.70 %) | 5.414 (+5.36 %) | 6.287 (+22.33 %) |
   | Ca I (6.113) | 5.027 (−17.77 %) | 4.323 (−29.28 %) | 4.925 (−19.43 %) |
   | Fe I (7.902) | 7.611 (−3.69 %) | 7.299 (−7.64 %) | 8.386 (+6.12 %) |

   O, N and Pb are in the same file. Print the full component chains from `tools/audit_components.py`. Fe's 3d⁶
   electrons belong in the **in** group, not df.
2. Item 2 hand-recompute gap: find which pa_* term the paper text doesn't specify (start from `tools/hand_recompute.py`).
   Then either describe it fully in §2.2 or point the reader to `tools/audit_components.py` as the reference
   implementation.
3. Edit `docs/paper_draft.md` with the ✅ findings (items 4–10, 19, Mendoza §4.5, grouping rule,
Table A1/A2, the "D" vs Dirac-factor notation clash, the reproducibility note), then rebuild the PDF with
`py -3.11 tools/build_paper_pdf.py`.

> **03:05 orchestrator check.** I verified the main workflow results independently:
> - frozen fingerprints intact (`handoff/FROZEN.txt`); `uni_predictions.csv` byte-identical to the frozen winner;
> - evaluate.py: 1.874% overall, 7.515% neutral atoms, 4.63% on the 311 experimental rows;
> - coverage test passes.
>
> **STILL RUNNING:** the competitor benchmark `wf_1a66b31c-6b0` (Faussurier 1997 / More 1982 / Mendoza 2011), started
> 02:56. See `handoff/BENCHMARK_RUNNING.md`. When it finishes, merge `docs/benchmarks.md` into the paper draft (Table 4)
> and `docs/literature.md`.

Last updated: **2026-10-05 02:59 (PKT)** by Claude (Opus 5.5), final fix pass (end of the final workflow `wf_2cf67b7c-ef5`).

> **The final workflow has finished** (script: `handoff/final_push_workflow.js`). All four stages are done:
> 1. **Push** (`models/push_a`, `models/push_b`): winner picked by the **pre-registered** selection score, mean MAPE of
>    V1 (fit Z≤36, validate Z 37–54), V2 (fit Z≤44, validate 45–54), S1 and S3. **S2 (fit Z≤54 → Z≥55) was never used
>    for any choice.** Caveat: S1/S3, as pre-registered, contain Z≥55 rows in train and test (~76 % of their test
>    rows), so heavy-atom *interpolation* accuracy did inform the choices; with S1/S3 restricted to Z≤54 the winner is
>    unchanged (`results/uni_sensitivity_z54.json`).
> 2. **Integrate** (final model, `validate_blind.py`, `tests/test_coverage.py`, `docs/unified.md`) and **Literature**
>    (`docs/literature.md`, `docs/references.bib`, Kregar/Di Rocco head-to-head).
> 3. **Final audit** (no critical, 1 major, 6 minor issues) and **paper draft** (`docs/paper_draft.md`).
> 4. **Final fix pass**: all audit issues fixed; see `docs/unified.md` §9 "Audit responses".

If you are a new Claude session and the user says **"continue"**: read this whole file, then start at
**§6 "What to do next"**.

---

## 1. The goal

Find a general analytical equation for the **successive ionization energies** IE(Z, N) of any atom or ion. Inputs are
the atomic number Z and the electron configuration. Use first principles where possible and pattern recognition on
experimental data otherwise, with shielding, e⁻–e⁻ repulsion, quantum, relativistic and QED corrections, and DFT.
Bohr's 13.6·Z²/n² eV is exact only for one-electron ions; Slater's rules are empirical.

The user wants:
- the work split across **two agents running in parallel**;
- speed;
- a paper at the end.

Ultracode (multi-agent workflows) was on. Paper name decided: **"Screened Rydberg formula"**. The user rejected
eponyms ("Awais formula") and acronyms ("ASH").

## 2. Environment (important)

- Project root: `E:\Dr. Navin\Np` on Windows 10. It is not a git repository.
- **Run everything with `py -3.13`.** The default `python` is 3.14 with no packages. Python 3.13 already has numpy,
  scipy, pandas, matplotlib, sympy, scikit-learn and numba. Don't create venvs.
- The machine has 4 CPU cores, so a workflow runs at most 2 agents at once.
- **Lesson:** never `SendMessage` to an agent running inside a Workflow. It spawns a duplicate that overwrites files.
  This happened once and was cleaned up.
- **Lesson:** workflow script files must have LF line endings (CRLF is rejected as "control characters"). Write them
  in binary mode from Python.

## 3. Layout

```
data/nist_ie.csv              NIST ASD successive IEs: 5847 rows, Z = 1..110, all ion stages (clean_nist.py rebuilds it from nist_ie_raw.csv)
common/atomdata.py            SHARED loader: load_records(), get(Z,N), ground_shells(Z,N), madelung_shells(N), constants
evaluate.py                   SHARED scorer: py -3.13 evaluate.py results/<name>_predictions.csv  (CSV columns Z,N,IE_pred_eV)
models/first_principles/      Track A: exact 1/Z expansion, Dirac/recoil/finite-size/QED, own Kohn-Sham LSDA solver → docs/first_principles.md
models/semi_empirical/        Track B: baselines + fitted GSHM (Generalized Screened-Hydrogenic Model)            → docs/semi_empirical.md
models/unified/               Unified model: umodel.py, final.py (dispatcher), pocket.py, run_unified.py, report.py, referee_fixes.py → docs/unified.md
ionization.py                 Top-level API/CLI: py -3.13 ionization.py Fe | py -3.13 ionization.py 118 | ionization_energy(Z,N), successive_ionization_energies(Z)
results/                      *_predictions.csv, *_metrics.json (prefixes fp_, se_, uni_), model_comparison.{csv,md}, uni_validation.json, figures/*.png
handoff/builds.json           Final reports of Track A and Track B (structured)
handoff/referee_report.json   Adversarial referee report on both tracks (structured)
handoff/unify_result.json     Unifier's final report (structured)
handoff/unify_workflow.js     Workflow script: Unify || Referee → Final fix. Its final-fix prompt has the 7 REQUIRED items (see §6)
```

Fixed held-out splits, used by every fitted model:
- **S1:** test = rows with Z % 5 == 0.
- **S2 (blind extrapolation):** train on Z ≤ 54, test on Z ≥ 55.
- **S3:** test = rows with N % 6 == 0.

## 4. Results so far

All numbers are MAPE over the 5847 NIST rows unless noted.

| Model | Fitted params | All ions | Neutral 1st IE | H-like | Blind S2 (Z≥55) | S1 / S3 |
|---|---|---|---|---|---|---|
| Bohr | 0 | 1363% | 22878% | 6.0% | 1517% | — |
| Slater's rules (total-energy difference) | 0 | 11.8% | 51.8% | 6.0% | 11.6% | — |
| Track A `fp_zexp` (exact 1/Z + rel/QED, closed form) | 0 | 69.8% | ~1100% | **0.0012%** | 77.8% | — |
| Track A LSDA ΔSCF DFT (207 rows only) | 0 | 1.32% | 3.30% | — | n/a | — |
| Track B GSHM (final) | 32 | 1.57% | 7.6% | 0.19% | **170–178% blind** (19% only with a post-hoc bound) | 1.59 / 1.70 |
| Unified pocket formula | 8 | 4.68% | 16.8% | 1.16% | **11.3%** | 4.62 / 5.24 |
| Unified previous final `uni_u35` (exact σ₁ + linear remainder) | 35 | **1.48%** | **6.43%** | **0.0012%** | 57.2% (median ~2%) | 1.52 / 1.65 |
| **FINAL `pa_hier_rel`** (exact σ₁ + bounded hierarchical remainder; pre-registered V1/V2/S1/S3 winner, score 2.148) | 33 | 1.87% | 7.52% | **0.0012%** | **6.60% blind** (median 1.62%, neutral 21.8%) | 1.91 / 2.00 |

`results/uni_predictions.csv` is now the FINAL `pa_hier_rel` model (all-data fit), verified 2026-10-05 02:23 with evaluate.py: 1.874% overall,
0.940% median, 7.515% neutral, 0.001% H-like. Full protocol table: `results/uni_validation.md` (validate_blind.py).

The final model was chosen by the **pre-registered selection score** (mean of V1, V2, S1, S3): pa_hier_rel 2.148 <
pb_clip_pos 2.370 < pa_hier 2.374 < pa_bound9 2.514 < … < u35 7.76. The earlier stage-1 rule (mean of S1, S2, S3;
u35 20.1, u29 47.3, gshm_qed 57.8) is **superseded**: it used S2 for selection, so it was not blind.

Key findings:
- **Track A.** The exact first-order screening constants σ₁ = −n²·ΔE₁ (rationals built from hydrogenic Slater
  integrals) are tabulated for N = 1..110 in `results/fp_zexp_coefficients.csv`. They act as ab-initio replacements
  for Slater's σ. Validated: He 5/8; Li-like E₁ = 5965/5832.
- **Track A, H-like.** Error by layer: Dirac 0.19%, then recoil + finite nuclear size 0.113%, then one-loop QED
  0.0012%. **Caveat:** the QED values come from the Yerokhin–Shabaev 2015 tables, and NIST's H-like references come
  from the same theory, so that layer shows consistency, not independent validation.
- **Track A, DFT.** Agreement with NIST to 1e-6 Ha is against NIST's LDA tables, with the same functional. That
  verifies the implementation only. The physical LDA error for Ne is about 0.7 Ha.
- **Track B.** The penetration charge follows `p_inf − tau/(Za + kappa)`, with one universal kappa = 8.2 and
  Za = Z − N + 1. **This is a rediscovery of Edlén-type isoelectronic behaviour, not a new law.** The referee
  flagged this.
- **Thesis for the paper.** Derived screening (Track A) covers what data fitting cannot extrapolate (Track B), and the
  reverse. Evidence: blind S2 went from 178% to 57% when exact σ₁ was added. **That is not yet under the ~10% bar**
  needed for a strong paper. The median heavy-element error is only about 2%; the mean is driven by a few heavy
  **neutral and near-neutral** atoms (7s, 6d and 5f elements) that blow up (u35 S2 neutral MAPE ~3800%). The
  8-parameter pocket formula gets 11.3% blind, which suggests the extra fitted terms extrapolate badly.

Referee verdict (handoff/referee_report.json): both tracks honest and reproducible, no leakage, no critical issues.
- **Major (fixed):** Track A crashed for Z > 110. relativity.py now extrapolates the QED and finite-size tables.
  `ionization.py 118` and `ionization_energy(120, 2)` work.
- **Major (addressed):** Track B's 19% S2 used a post-hoc bound. Blind S2 is now reported, and the final model was
  re-selected on blind scores.
- **Minor:** Track A's "zero-parameter" formula contains heuristic modelling choices. GSHM core extrapolates better
  than GSHM final. Structural leakage: Track B's model form was discovered on all rows.

## 5. Status checklist

- [x] Data, shared loader and scorer
- [x] Track A build: verified
- [x] Track B build: verified
- [x] Referee: done (handoff/referee_report.json)
- [x] Unify: done (handoff/unify_result.json, docs/unified.md, ionization.py, model_comparison.*)
- [x] (superseded by the integration stage below) **Final fix pass: INTERRUPTED at 01:37**, about 4 minutes into a 30-minute time box. What it finished:
  - the Z > 110 crash fix (relativity.py, test_superheavy.py);
  - a re-run of the blind S2 for all candidates;
  - re-selection of `uni_u35` as final;
  - regenerated uni_* results, model_comparison and figures.

  What it probably did **not** finish:
  - updating docs/unified.md, which still describes the earlier gshm_qed choice (check it);
  - a "Referee responses" section;
  - a "Literature positioning / novelty" section;
  - the Edlén framing;
  - `models/unified/validate_blind.py`;
  - `tests/test_coverage.py`;
  - the structured summary.
- [x] Headline blind-test push (S2 < ~10%): DONE, pre-registered winner pa_hier_rel, blind S2 6.60%. **Push A done (2026-10-05)**. `models/push_a/` chose `pa_hier_rel`
  (33 params, selection score 2.148 = mean of V1/V2/S1/S3). Blind S2 is **6.60%** (median 1.62%, neutral 21.8%).
  All-data: 1.87%, neutral 7.5%, H-like 0.0012%. Details: models/push_a/NOTES.md, results/pa_validation.json.
  Push B was compared in the integration stage (pb_clip_pos 2.370, second).
  - [~] Approach B (push_b, 2026-10-05): models/push_b/{model,validate}.py, NOTES.md, results/pb_*. Chosen by pre-registered V1/V2/S1/S3 score: pb_clip_pos (30 params; Za<=Zeff<=Z smooth clip, t>=0): score 2.37 (u35 7.76), all-data 1.745%, neutral 7.49%, blind S2 13.9% (median 1.56%, neutral 248%). Misses 10%. Non-selected pb_exp_pos had S2 6.0% (cannot be switched to: S2-selection).
- [x] **Integration stage (2026-10-05 02:23)**: final model = `pa_hier_rel` (models/unified/final.py default; ionization.py CLI unchanged,
  `--model` option). `models/unified/validate_blind.py` reproduces Push A/B numbers exactly (diff 0.0) for final, u35,
  pocket; also GSHM final/core blind (170 / 178.5), u29, Slater total. `results/model_comparison.{csv,md}` now has V1/V2/S1/S3/
  selection/blind S2 for all models + push candidates; uni_* figures regenerated (`models/unified/report.py`).
  docs/unified.md rewritten (headline table, equation + parameter table + lineage, protocol history, referee responses,
  Edlén framing, worked O / Mg III examples, limitations; 'Literature positioning and novelty' -> docs/literature.md).
  run_unified.py is now legacy (writes uni_legacy_validation.json, no longer overwrites uni_predictions.csv).
- [x] Coverage test Z = 1..118 (`py -3.13 tests/test_coverage.py`): 7021/7021 finite and > 0, 0 crashes, PASS;
  22 monotonicity violations (Pt–Bi N≈60–62; Rf–Ds at NIST configuration jumps). Lr negative (−3.3 eV) in the Z≤54
  S2 fit (negative r_c); Og 3.66 eV vs ≈8.9 lit. Documented in docs/unified.md §7.
- [x] **Literature + head-to-head stage (2026-10-05 ~03:05)**: docs/literature.md, docs/references.bib (each entry
  tagged VERIFIED/CITED-IN/STANDARD/UNVERIFIED). More 1982 / Faussurier 1997 constants NOT obtainable (not fabricated).
  Implemented the parameter-free Kregar/Di Rocco SHM (models/literature/kregar_shm.py; validate_kregar.py vs
  Pomarico 2005 tables; run_lit.py; compare.py; fill_doc.py). On all 5847 rows: Dirac variant 12.5% MAPE (median
  1.7%), neutral 224%, Z>=55 13.6%, charge>=3 4.7%; results/lit_*. Project pa_hier_rel beats it (1.87%, blind S2 6.6%).
- [x] **Final audit (2026-10-05 ~02:47)**: independent re-implementation of the protocol reproduced every headline
  number (selection 2.148, blind S2 6.603, all-data 1.874; 33 params; no S2 leakage; coverage PASS). No critical
  issues; 1 major (false claim that no Z≥55 row influenced any choice) and 6 minor.
- [x] **Paper draft (2026-10-05)**: `docs/paper_draft.md` + `results/figures/paper_blind_s2_by_Z.png`,
  `paper_selection_vs_blind.png`. Open [TODO]/[VERIFY] marks remain (see §6).
- [x] **Final fix pass (2026-10-05 02:59)**: all audit issues fixed, numbers unchanged; responses in `docs/unified.md` §9.
  - Z≥55 wording corrected in unified.md §1/§3/§7, paper abstract/§3.2/Fig. 2 caption, push_b/NOTES.md, this file,
    model_comparison.md and uni_validation.md; new robustness check `models/unified/sensitivity_z54.py` →
    `results/uni_sensitivity_z54.json` (S1/S3 restricted to Z≤54: pa_hier_rel 1.844 < pa_hier 1.989 < pa_bound9
    2.094 < pb_clip_pos 2.118; reporting only).
  - Exploration size now "about 57 variants (43 Push A + ~14 Push B)" everywhere; f-removal S2 error "2.0–2.5×
    (mean 2.2×, 23 rows)"; literature.md no longer "provisional", SHM described as "our implementation", 224 %
    neutral caveat; DFT rows' selection score shown as n/a in model_comparison (report.py); bound condition
    0 ≤ σ₁ ≤ N−1 stated (docs + paper) and checked with a RuntimeWarning in `final.predict_many`.
  - Paper: experimental-row metrics filled in (311 rows: 4.63 %, median 2.86 %).
  - Re-runs after the fixes: validate_blind.py reproduces uni_validation.json exactly and uni_predictions.csv byte for
    byte; evaluate.py 1.874 / 0.940 / 7.515 / 0.001; test_coverage PASS (7021/7021 > 0, 22 violations);
    `ionization.py Fe` and `ionization.py 118` OK.
- [ ] Remaining for the user (not agent work): final verdict read-through, venue choice, author verification of
  every term and reference, resolving the paper's [TODO]/[VERIFY] marks.

## 6. Final headline (verified 2026-10-05 02:59) and what remains

| model | params | all-data MAPE (median) | neutral 1st IE | H-like | V1 | V2 | S1 | S3 | selection | **blind S2** (median / neutral) |
|---|---|---|---|---|---|---|---|---|---|---|
| **FINAL pa_hier_rel** ("Screened Rydberg formula") | 33 | 1.87 (0.94) | 7.52 | 0.0012 | 3.09 | 1.60 | 1.91 | 2.00 | **2.15** | **6.60** (1.62 / 21.8) |
| Push B pb_clip_pos | 30 | 1.75 | 7.49 | 0.0012 | 4.04 | 1.84 | 1.76 | 1.85 | 2.37 | 13.9 (1.56 / 248) |
| previous final u35 | 35 | 1.48 | 6.43 | 0.0012 | 26.4 | 1.48 | 1.52 | 1.65 | 7.76 | 57.2 (2.00 / 3831) |
| pocket formula | 8 | 4.68 | 16.8 | 1.16 | 4.4e6 (diverges) | 2.55 | 4.62 | 5.24 | – | 11.3 (2.04 / 12.3) |
| Track B GSHM final (blind) | 32 | 1.57 | 7.61 | 0.19 | 34.1 | 1.65 | 1.59 | 1.70 | 9.77 | 170 (3.51 / 7393) |
| Slater total-energy | 0 | 11.8 | 51.8 | 6.00 | 13.4 | 14.4 | 11.6 | 12.4 | 12.9 | 11.6 (8.07 / 50.7) |
| Kregar/Di Rocco SHM (our implementation, Dirac) | 0 | 12.5 | 224 | 0.19 | – | – | – | – | – | 13.6 on Z≥55 rows |

Sources: `results/uni_validation.md`, `results/model_comparison.md`, `docs/literature.md`. Full write-up:
`docs/unified.md` (technical) and `docs/paper_draft.md` (paper).

If the user says **"continue"**, the project work is done; what remains is the user's:
1. **Read the honest verdict**: a solid, honest paper, not a breakthrough. Most original piece: the exact σ₁ table
   (Track A). The 1/(Z_a+κ) remainder is an Edlén-type rediscovery. Blind S2 6.6 % beats Slater (11.6 %) and our
   implementation of the Kregar/Di Rocco SHM (13.6 %), but heavy neutral atoms remain weak (21.8 % blind).
2. **Venue choice**: J. Phys. B, Atoms (MDPI), Eur. J. Phys. / Am. J. Phys. (pedagogical angle, pocket formula),
   HEDP, JQSRT.
3. **Author verification**: every [VERIFY]/[TODO] in `docs/paper_draft.md` (NIST ASD version and access date,
   Rodrigues 2004 Dirac-Fock comparison, ML prior-work search, novelty of the σ₁ table vs Safronova / Layzer 1964,
   DOI for the code archive, acknowledgements, author details, every STANDARD/CITED-IN/UNVERIFIED reference in
   `docs/references.bib`, the approximate superheavy literature values in `docs/unified.md` §7). Disclose AI
   assistance; be able to defend every term.
4. **Optional further work** (each must be chosen on V1/V2/S1/S3 only, and disclosed as post-audit):
   - a physically bounded relativistic term (fitted r_c are negative; Lr −3.33 eV in the Z≤54 fit; Og 3.66 eV is
     implausibly low);
   - better f-electron removal (S2 third IEs 2.0–2.5× too high);
   - extending LSDA ΔSCF to all neutral atoms (only 207 rows now);
   - the lighter 9-parameter pa_bound9 (selection 2.51, S2 7.80 %) as a paper alternative;
   - publishing the paper as an artifact page if the user wants a shareable view.
