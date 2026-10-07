# v2 experiments (pre-registered in docs/preregistration_v2.md; hash in PREREG_HASH.txt)

**Harness.** `py -3.11 models/v2/run.py <variant>`. Results are appended to `results/v2_log.json`.

**Baseline.** v1 refit in this harness scores **2.1476**; the pre-registered value is 2.1473 from validate_blind.
The acceptance threshold is 2.1273 (the stricter baseline minus 0.02).

## Round 1 (2026-10-06), all scored variants

The four held-out columns are each split's neutral-atom MAPE (%).

| variant | item | params | selection score | V1 neutral | V2 neutral | S1 neutral | S3 neutral | all-data neutral | all-data | mono | Og (eV) | S2 non-blind (Lr) | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| v1_baseline | – | 33 | 2.1476 | 12.5 | 14.2 | 6.71 | 7.95 | 7.52 | 1.874 | 22 | 3.66 | 6.60 (−3.32) | reference |
| relexp | 1 (bounded rel.) | 33 | 2.4919 | 21.4 | 6.41 | 7.06 | 7.93 | 8.67 | 1.987 | 19 | 9.34 | 9.10 (+5.96) | **rejected** (score) |
| relexp_q2 | 3 (neutrals) | 38 | 2.1778 | 36.4 | 9.48 | 8.37 | 8.95 | 9.35 | 1.779 | 22 | 10.24 | 7.20 (+6.76) | **rejected** (score) |
| relexp_kl | 3 (neutrals) | 36 | 2.4919 | 21.4 | 6.41 | 7.06 | 7.93 | 8.67 | 1.986 | 19 | 9.34 | 9.10 (+5.96) | **rejected** (score; κ_l stays at 0) |
| relexp_q2_kl | 3 (neutrals) | 41 | 2.3568 | 36.3 | 8.48 | 7.65 | 7.41 | 8.25 | 1.767 | 22 | 9.88 | 102.9 | **rejected** (score) |

## Findings

- **The bounded relativistic factor (relexp) does what physics asks, but the score says no.** It is
  R = exp(r_c(Zα)²(1 − Z_a/Z_eff)/n), which stays positive for any parameter values. Compared with v1:
  - Lr is positive in the Z ≤ 54 fit;
  - Og comes out at 9.3 eV, below its congener Rn (NIST 10.75), instead of 3.7;
  - monotonicity violations fall from 22 to 19.

  But the selection score is worse (2.49 vs 2.15), driven by V1 (3.98 vs 3.08), and so is the all-data neutral
  MAPE. Rejected under the pre-registered rule.
- **The second-order charge term (q2) does not help neutral atoms in extrapolation.** It improves the in-sample and
  interpolation errors: all data 1.78 %, S1 1.83 %. But held-out neutral error gets much worse in V1 (36 %), so the
  term overfits the Z-dependence of near-neutral screening. Rejected.
- **Variants used so far:** item 1, 1 of 5; item 3, 3 of 5.
- **Why neutral gains cannot pass this score.** Neutral atoms are about 2 % of the test rows, so the pre-registered
  score is almost blind to them. Any neutral-specific gain must not cost the ions anything, or it cannot pass. A
  neutral-specific criterion was not pre-registered. Adding one now must be done as a dated amendment *before* any
  further variant is scored.

**Verdict after round 1: v1 (pa_hier_rel) remains the final model.**
