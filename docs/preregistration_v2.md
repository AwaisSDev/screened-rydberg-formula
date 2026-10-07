# Pre-registration for model version 2 (Phase 2 improvements)

**Written 2026-10-06, before any v2 model was built, fitted or scored.** It must be committed to git before the
first v2 experiment, so that the commit timestamp proves the order. Any later change to it must be a new, dated section; the original text stays.

## 1. What stays frozen

- `pa_hier_rel` is version 1, the paper's final model. Its files are not modified:
  - `models/push_a/model.py`;
  - `results/pa_params.json`;
  - `results/pa_hier_rel_predictions.csv`;
  - `results/uni_*`;
  - the hashes in `docs/FROZEN.txt`.
- Every v2 model lives in a new folder, `models/v2/`, with new result-file prefixes `v2_*`.

## 2. Selection score (the only quantity used to choose)

**score = mean MAPE (%) over the four held-out splits V1, V2, S1, S3**, each refit with the candidate's own fitting
routine on the training rows:

| split | train | test |
|---|---|---|
| V1 | Z ≤ 36 | 37 ≤ Z ≤ 54 |
| V2 | Z ≤ 44 | 45 ≤ Z ≤ 54 |
| S1 | Z mod 5 ≠ 0 | Z mod 5 = 0 |
| S3 | N mod 6 ≠ 0 | N mod 6 = 0 |

All splits use the 5847 rows of `data/nist_ie.csv`, scored with the definitions of `evaluate.py`.

**Baseline.** The baseline is `pa_hier_rel` refit by the same code in the same software environment. On this machine
(Python 3.11.9, NumPy 2.4.4, SciPy 1.17.1) that is **2.1473**. On the frozen 3.13 run it was 2.1481. Every
comparison uses scores computed in the same environment.

## 3. Acceptance rule (fixed now)

A candidate change is **accepted** only if all of the following hold:

1. Its selection score is lower than the baseline by **at least 0.02 percentage points**. Smaller differences count
   as ties, and a tie is resolved in favour of the simpler model (fewer parameters, then the frozen one).
2. **Hard physical checks** on the all-data fit, evaluated on all 7021 coverage cases (Z = 1–118, N = 1–Z):
   - every IE is finite and > 0;
   - the number of monotonicity violations IE(Z, N−1) > IE(Z, N) is **not larger** than v1's (22).
3. Its parameter count is reported. Any parameter added beyond v1's 33 must be justified physically in the notes.

Changes are evaluated **one at a time against the current accepted model**, in the order of §4, and accepted
cumulatively. A change that fails is dropped and reported in `models/v2/NOTES.md` with its score.

## 4. Planned changes (in this order)

1. **A physically bounded relativistic term.** v1's fitted r_c are all negative and Lr comes out negative. The form
   must be bounded for any parameter values so that IE stays > 0, and it should reduce to the Fermi–Segrè scaling at
   small Zα.
2. **f-electron removal.**
3. **Heavy p-block neutrals.**
4. **Monotonicity.** Either fix the 22 violations in the model or flag them in the API output. Flagging is not a
   model change and needs no score.

At most **5 variants per item** will be scored, including failed ones. Every variant scored is logged with its
score, so the number of models tried is public.

## 5. S2 is no longer blind

The S2 split (fit Z ≤ 54, predict Z ≥ 55) and its failures (Pb, Tl, Rn, Lr, f removal) are already known. **S2 is
not used for any v2 decision.** It is reported only afterwards, labelled "non-blind (failures known before
development)".

## 6. Honesty rules

- No number in the paper, README or notes is typed from memory. Each one comes from a script output that is named
  next to it.
- Every failed or rejected variant is reported.
- No tuning against the all-data fit, S2 or individual elements.
- If no change passes §3, v1 stays the final model and the paper says so.
