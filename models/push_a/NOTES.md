# Push A: robust minimal-remainder Screened Rydberg formula

Files: `model.py` (predict / predict_many / fit), `validate.py` (full protocol, writes `results/pa_validation.json`,
`results/pa_params.json`, `results/pa_*_predictions.csv`), `explore.py` + `explore_log.txt` (exploration log, run on the
selection splits V1/V2/S1/S3 only).

Run: `py -3.13 models/push_a/validate.py` (~25 s).

## Protocol (pre-registered, identical for all candidates)
V1 fit Z<=36 / val 37..54; V2 fit Z<=44 / val 45..54; S1 fit Z%5!=0 / test Z%5==0; S3 fit N%6!=0 / test N%6==0.
Selection score = mean MAPE(V1, V2, S1, S3). S2 (fit Z<=54 / test Z>=55) was computed once, in the single
`validate.py` run made after the candidate list and every hyper-parameter were frozen. All-data fit scored with
evaluate.py. The references u35 (ridge 0.02, kappa_l, S2 terms) and pocket (8 params) were refitted per split with
their own fit functions.

## Final formula (chosen: `pa_hier_rel`, lowest selection score)

    IE = { Ry Zeff^2/n^2 * D_{n,j}(Zeff) * [1 + r_c (Z alpha)^2 (Zeff/Za - 1)/n] + Ry x_l K_l(k)/n^2 } * mu(Z)
         - [dE_QED + dE_FNS](Z, n) (Zeff/Z)^2                       (ns electrons with n <= 2 only)
    Zeff = Z - sigma1 - D
    T    = sum_{g in pocket groups} tau_g nu_g + sum_{c in Track-B classes} dtau_c nu_c
    D    = T / (Za + kappa + |T|/h),  h = (N-1) - sigma1 if T >= 0, else h = sigma1

- Z: nuclear charge. N: number of electrons. Za = Z - N + 1 (asymptotic charge). (n, l, j, k): removed subshell,
  j from jj filling (l-1/2 first), k = its occupancy.
- sigma1: exact first-order (1/Z) screening constant of the removed electron (Track A, `models/unified/abinitio.py`,
  0 parameters). It keeps the exact Z^2 and Z coefficients of the 1/Z expansion.
- D: the only fitted screening term. Algebraically D/h = u/(1+u) with u = |T|/(h (Za+kappa)). So **Za <= Zeff <= Z
  holds for every parameter value**. Zeff >= Za is the non-penetrating hydrogenic limit, and Zeff <= Z means no net
  anti-screening. The saturation is smooth (no clip). As Za -> infinity at fixed N, D -> T/(Za+kappa) = O(1/Z), so
  sigma -> sigma1.
- nu_g: electron counts in the pocket groups same / in / core / df / out (pocket.py definitions). nu_c: counts in
  the Track B screening classes, with 'core' split by target type (s,p vs d,f). The class deviations dtau_c carry a
  ridge penalty of 1e-4 * n_rows * dtau^2, which shrinks them towards the pocket-group value. A class with no
  training electrons gets no parameter, so it inherits its pocket-group value exactly. The f classes therefore fall
  back to physically grouped values, not arbitrary ones.
- D_{n,j}(Zeff): exact point-Dirac / Schroedinger ratio (parameter-free). The bracket with r_c (s, p1/2, p3/2, d, f)
  is the Fermi-Segre-type penetration correction. r_f is tied to r_d when the training set has no f removals.
- x_l K_l(k): Hund exchange kink (Track B), x_f tied to x_d when there are no f rows.
- mu, dE_QED, dE_FNS: Track A reduced-mass / QED / finite-size layer (parameter-free). H-like MAPE is 0.0012%.

Parameter count (all-data fit): 5 tau_g + kappa + 19 dtau_c + 5 r_c + 3 x_l = **33**. The 19 dtau are
ridge-shrunk, so they are less than 19 free degrees of freedom, but all 33 are counted.

Parameter values (all-data): see `results/pa_params.json` -> candidates.pa_hier_rel. The main values are
tau_same 0.393, tau_in 0.788, tau_core 2.224, tau_df 2.549, tau_out 10.57, kappa 1.753, r_s -0.454, r_p1 -0.813,
r_p3 -1.153, r_d -0.843, r_f -1.467, x_p 0.205, x_d 0.440, x_f 0.461.

## Results (MAPE %; from results/pa_validation.json)

| candidate | params | V1 | V2 | S1 | S3 | **selection** | S2 blind (median / neutral) | all-data (median / neutral / H-like) |
|---|---|---|---|---|---|---|---|---|
| **pa_hier_rel** (chosen) | 33 | 3.09 | 1.60 | 1.91 | 2.00 | **2.148** | **6.60** (1.62 / 21.8) | 1.87 (0.94 / 7.52 / 0.0012) |
| pa_hier (no fitted rel) | 28 | 3.47 | 1.53 | 2.16 | 2.34 | 2.374 | 8.47 (1.78 / 35.2) | 2.16 (0.99 / 9.65 / 0.0012) |
| pa_bound9 (pocket groups + bound, no rel) | 9 | 2.16 | 1.95 | 2.90 | 3.05 | 2.514 | 7.80 (1.69 / 28.0) | 2.90 (1.39 / 11.99 / 0.0012) |
| pa_bound14_relfit | 14 | 3.85 | 2.83 | 2.46 | 2.48 | 2.905 | 5.53 (1.92 / 47.2) | 2.41 (1.31 / 8.42 / 0.0012) |
| pa_bound9_fs0 (parameter-free Fermi-Segre rel) | 9 | 2.73 | 2.58 | 3.21 | 3.39 | 2.980 | 8.49 (2.50 / 41.2) | 3.24 (1.78 / 15.5 / 0.0012) |
| pa_nobound14_relfit (ablation: no bound) | 14 | 6.61 | 3.16 | 2.91 | 3.15 | 3.957 | 9.43 (3.22 / 105.7) | 2.88 (1.81 / 13.1 / 0.0012) |
| ref_u35 | 35 | 26.39 | 1.48 | 1.52 | 1.65 | 7.760 | 57.2 (2.00 / 3831) | 1.48 (0.73 / 6.43 / 0.0012) |
| ref_pocket | 8 | 4.4e6 (fit diverges) | 2.55 | 4.62 | 5.24 | 1.1e6 | 11.33 (2.04 / 12.3) | 4.68 (3.04 / 16.8 / 1.16) |

## What each change bought on the selection score (explore_log.txt)
- 4-group remainder + bound + fitted rel (13 p): 6.57. Without the bound: 6.28. The 4-group form is a poor base.
- Pocket grouping (same / in / core / df / out) instead: 2.91 (14 p). Removing the bound from it: 3.96. **The bound is
  worth about 1.05 points**, mostly on V1 (6.6 -> 3.9) and the neutral validation rows.
- Dropping the fitted relativistic term: 2.51 (9 p). The parameter-free Fermi-Segre term (fs0): 2.98, worse. Per-class
  or fraction-based relativistic forms were also worse (3.1-3.6).
- Saturation shape: rational u/(1+u) 2.51, exp 2.92, u/sqrt(1+u^2) 2.95, tanh 3.15. Rational kept.
- Hierarchical class deviations with ridge (shrink to pocket groups): ridge 0 diverges on V1 (60.9). Ridge 3e-5 gave
  2.57, 1e-4 gave 2.37, 3e-4 gave 2.42, 1e-3 gave 2.51, 1e-2 gave 2.54. **1e-4 chosen.**
- Hierarchical + fitted rel (Fermi-Segre form, as in u35): **2.148** (chosen). With ridge 3e-4: 2.157. With per-class
  "fsc" or "frac" rel forms: 2.68 / 2.86. Adding kappa_l gave no gain (2.374).

## Coverage (chosen model, after freezing)
Z = 1..120, N = 1..Z gives 7260 predictions: all finite and > 0. There are 22 monotonicity violations
(IE(N+1) > IE(N)).
Superheavy neutral 7p atoms come out low: Og (118) 3.66 eV against about 8.9 eV in the literature. This is a known
weakness, not fixed, because fixing it now would be post-hoc.

## Disclosures
- Before designing anything I knew from earlier project notes that u35's S2 failures are heavy near-neutral
  atoms (Po, Hs, No, Ac, Tl, Ta, Sg, Hf, Md). The idea of bounding Zeff between Za and Z (and so protecting the neutral
  rows) was suggested by the task brief with that knowledge in the background. Its *selection* was made on V1/V2/S1/S3
  only (ablation above).
- I also knew that the pocket formula extrapolates relatively well to S2 (11.3%). That is why the pocket groups were
  tried first.
- The S2 numbers above come from the single validate.py run, made after the candidate list was frozen. No change was
  made after seeing them.
- The ref_pocket reference fit diverges on V1 with pocket.fit's own settings (unbounded Zeff, init values). It is
  reported as is and was not repaired.
- The fitted r_c are negative (Dirac(Zeff) over-binds and r_c absorbs other O((Z alpha)^2) effects), so the relativistic
  bracket is not physically bounded. It stayed positive over Z <= 120 but makes the superheavy 7p neutrals too low.
- The S2 neutral first-IE MAPE is still 21.8%, so heavy neutral atoms remain the weak point. The all-data neutral
  MAPE is 7.5%, against 6.4% for u35.
- Track B's class definitions and the Hund kink were discovered earlier on all rows (structural leakage already
  flagged by the referee). That applies here too.
