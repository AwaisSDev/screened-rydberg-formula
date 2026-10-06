# Push-B: a physics-constrained full model (exact σ₁ plus a bounded remainder)

Reproduce with `py -3.13 models/push_b/validate.py` (about 40 s). It writes `results/pb_validation.json`,
`results/pb_params.json` and `results/pb_<candidate>_predictions.csv`.

## Protocol (pre-registered, the same for every candidate)

The selection score is the mean MAPE over four splits:
- V1: fit Z ≤ 36, validate on 37–54.
- V2: fit Z ≤ 44, validate on 45–54.
- S1: hold out Z % 5 = 0.
- S3: hold out N % 6 = 0.

All choices were made on this score only. Before S2 was run, `validate.py` fixed the candidate list and chose the
final model in code (lowest score among `pb_*`). S2 (fit Z ≤ 54, test Z ≥ 55) was computed **once** per candidate,
for reporting only. The u35 and pocket references were refit for every split with their own fit functions. u35
used its own ridge of 0.02, which an earlier session chose on the inner split.

## Final formula (chosen candidate `pb_clip_pos`, 30 fitted parameters)

    IE = { Ry Zeff²/n² · D_{n,j}(Zeff) · [1 + (Zα)² r_{lj} (Zeff/Za − 1)/n] + Ry x_l K_l(k)/n² } · μ(Z)
         − [ΔE_QED + ΔE_FNS](Z, n) · (Zeff/Z)²          (ns electrons with n ≤ 2 only; Track A, no parameters)

    Za     = Z − N + 1                                   asymptotic charge
    Zraw   = Z − σ₁(config) − Σ_c ν_c t_c / (Za + κ + κ_l),   with t_c ≥ 0
    Zeff   = Za + (Z − Za) · h((Zraw − Za)/(Z − Za)),     h(u) = [ln(1+e^{βu}) − ln(1+e^{β(u−1)})]/β,  β = 12
             (Zeff = Z for N = 1)

Symbols:
- σ₁ is the exact first-order screening constant (Track A, rational, no parameters).
- ν_c counts the other electrons in screening class c (the 18 Track B classes).
- D_{n,j} is the exact Dirac/Schrödinger ratio for a point charge Zeff.
- K_l(k) is the Hund exchange kink.
- μ is the reduced-mass factor.
- h is a smooth clip of the real line onto (0, 1); it is close to the identity on about [0.15, 0.85].

Constraints, decided from physics and confirmed only on V1, V2, S1 and S3:
1. **Za ≤ Zeff ≤ Z, smoothly, for any parameters and any electron count.** This means IE ≥ Ry Za²/n², the
   non-penetrating hydrogenic limit. On Z ≤ 54 data, the hydrogenic bound is broken by only 2 of the 1,485 rows,
   and only at the non-relativistic level.
2. **t_c ≥ 0.** Beyond first order, relaxation and correlation can only add screening relative to the hydrogenic σ₁
   value. For neutral atoms σ₁ badly under-screens: for Rb, D1 = N − 1 − σ₁ = 7.6 while about 1.2 is needed. This
   bound stops the fit from cancelling large positive and negative class terms, which is what blew up extrapolation
   into a new shell.
3. **No data, no parameter.** f-class t values are tied to their d analogues (as in u35). Two more a-priori analogues
   are used: n2_sp_df → n2_sp_sp and n1_d_d → d_near. Any class still unsupported by the training rows is fixed at
   t = 0, which is the pure ab-initio σ₁. With all data every class is supported, so these ties only act inside the
   split fits.
4. The O(1/Z) structure of u35 is kept: as Za → ∞ at fixed N, Zeff → Z − σ₁ + O(1/Z). The exact Z² and Z
   coefficients are preserved, and H-like ions are unchanged (0.0012%).

Dropped from u35:
- the second-order τ₂/(Za+κ)² terms;
- the exp(b (Zα)²) correction factors.

Together these remove 5 parameters.

### Parameter table (all-data fit, `results/pb_params.json`)

| param | value | param | value | param | value |
|---|---|---|---|---|---|
| t_same_s | 0.0600 | t_n2_sp_df | 2.0069 | kappa | 7.4057 |
| t_same_p | 0.8659 | t_d_near | 1.8894 | kappa_p | 0.0000 (at bound) |
| t_same_d | 1.7159 | t_n1_d_d | 2.1835 | kappa_d | 0.8867 |
| t_same_f | 1.7749 | t_n1_d_f | 3.1934 | kappa_f | 0.0001 |
| t_sn_in_p | 0.6947 | t_f_near | 3.3469 | rel_s | 0.1861 |
| t_n1_sp_sp | 1.1587 | t_n1_f_f | 1.5773 | rel_p1 | −0.4683 |
| t_n1_sp_d | 1.9888 | t_core | 2.4114 | rel_p3 | −0.9166 |
| t_n1_sp_f | 1.9947 | t_out_d | 1.7013 | rel_d | −0.9006 |
| t_n2_sp_sp | 1.6252 | t_out_f | 0.0000 (at bound) | rel_f | −1.8430 |
| x_p | 0.2051 | x_d | 0.2366 | x_f | 0.3143 |

There are 30 free parameters in total. Two of them sit at their bounds (kappa_p and t_out_f) but are still counted.

## What each change bought on the selection score (V1/V2/S1/S3 only)

| variant | V1 | V2 | S1 | S3 | score |
|---|---|---|---|---|---|
| ref u35 (35 params, unbounded) | 26.39 | 1.48 | 1.52 | 1.65 | 7.76 |
| ref pocket (8 params) | 4.4e6 (unbounded Zeff blow-up) | 2.55 | 4.63 | 5.24 | 1.1e6 |
| logistic Zeff = Za + 2D1/(1+e^X), t free | 90.0 | 1.73 | 1.73 | 1.86 | 23.8 |
| logistic, t ≥ 0 (before the analogue ties) | 15.4 | 2.32 | 1.99 | 2.04 | 5.44 |
| exp Zeff = Za + D1 e^{−X}, t ≥ 0, sat-rel (before the ties) | 13.9 | 2.20 | 1.86 | 1.90 | 4.97 |
| clip, t free, u35 rel (before the ties) | 9.13 | 1.75 | 1.76 | 1.84 | 3.62 |
| clip, t free, sat-rel | 447.6 | 1.81 | 1.74 | 1.80 | 113 |
| ridge toward σ₁ (0.0003–0.01 · √n), clip | worse at every strength tried (4.1–6.2) | | | | |
| clip + u35 exp(b(Zα)²) corr terms | 10.33 | 1.94 | 1.79 | 1.87 | 3.98 |
| **clip, t ≥ 0, u35 rel, with the ties (final)** | **4.04** | **1.84** | **1.76** | **1.85** | **2.37** |
| clip, t free, with the ties | 10.42 | 1.75 | 1.76 | 1.84 | 3.94 |
| exp, t ≥ 0, with the ties | 8.94 | 2.20 | 1.86 | 1.90 | 3.73 |
| logistic, t ≥ 0, u35 rel, with the ties | 9.73 | 2.35 | 2.05 | 2.14 | 4.07 |

The ridge was dropped: every strength tried made the score worse, so the final ridge is 0.

## Results (validate.py; MAPE / median APE / neutral first-IE MAPE, in %)

| candidate | params | V1 | V2 | S1 | S3 | score | **S2 blind** | all-data MAPE | all-data median | all-data neutral | H-like |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ref_u35 | 35 | 26.39 | 1.48 | 1.52 | 1.65 | 7.76 | 57.2 / 2.00 / 3831 | 1.484 | 0.728 | 6.43 | 0.0012 |
| ref_pocket | 8 | 4.4e6 | 2.55 | 4.63 | 5.24 | 1.1e6 | 11.3 / 2.04 / 12.3 | 4.684 | 3.04 | 16.8 | 1.16 |
| **pb_clip_pos (chosen)** | 30 | 4.04 | 1.84 | 1.76 | 1.85 | **2.37** | **13.9 / 1.56 / 248** | 1.745 | 0.788 | 7.49 | 0.0012 |
| pb_clip | 30 | 10.42 | 1.75 | 1.76 | 1.84 | 3.94 | 65.5 / 1.64 / 4132 | 1.744 | 0.787 | 7.46 | 0.0012 |
| pb_exp_pos | 30 | 8.94 | 2.20 | 1.86 | 1.90 | 3.73 | 6.0 / 2.27 / 40.4 | 1.864 | 0.774 | 10.36 | 0.0012 |
| pb_logistic_pos | 30 | 9.73 | 2.35 | 2.05 | 2.14 | 4.07 | 10.1 / 2.51 / 44.9 | 2.076 | 0.829 | 10.67 | 0.0012 |

**Honest verdict:** the selected model, pb_clip_pos, cuts blind S2 from 57.2% to 13.9%. Its median is the best of
all candidates at 1.56%. It does **not** reach the ~10% target. It still has neutral-atom overshoots (S2 neutral MAPE
248%), because its upper bound is Z, which is far too loose for heavy neutral atoms. The bounded-saturation forms
pb_exp_pos (6.0%) and pb_logistic_pos (10.1%) do much better on S2, but they scored worse on the pre-registered
selection splits, so they were **not** chosen. Switching to them now would be selecting on S2. The protocol forbids
that, and it would repeat the inflation the referee flagged earlier.

## Disclosures

- Before any work started, the task description told me that u35 fails on heavy near-neutral atoms (Po, Hs, No, Ac,
  Tl, Ta, Sg, Hf, Md; 6p/6d/7s/5d removal). That framing motivated looking at bounds and saturation at all. The
  specific forms (bounds, t ≥ 0, ties) were chosen from physical arguments plus V1, V2, S1 and S3 diagnostics. I
  looked at no Z ≥ 55 rows before S2.
- The t ≥ 0 constraint and the extra analogue ties were adopted after diagnosing the worst **V1** rows (Rb–Xe
  neutrals; t_core, t_n2_sp_df and t_n1_d_d took large values of opposite sign). That is allowed by the protocol
  (V1 is a selection split), but these are data-informed choices, not pure a-priori ones.
- About 14 variants were tried on the selection score (table above). Choosing the minimum among them carries some
  optimism into the selection score itself.
- S1 and S3 include Z ≥ 55 rows in both training and test sets (about 76 % of their test rows); the protocol defines
  them that way. So heavy-atom *interpolation* accuracy, through the S1/S3 part of the selection score, did inform the
  choice of model form. Only S2 (fit Z ≤ 54 → Z ≥ 55) was never used for any choice. (Corrected after the final audit;
  an earlier version of this note said no Z ≥ 55 row influenced the choice. With S1/S3 restricted to Z ≤ 54,
  pb_clip_pos scores 2.118 and stays behind pa_hier_rel's 1.844: `results/uni_sensitivity_z54.json`.)
- The pocket reference's V1 is a genuine blow-up (its Zeff is unbounded below and above). Its selection score is
  therefore huge. This does not reflect its S2 robustness.
- Track A's LSDA ΔSCF code was **not** used (time box).
- The `clip` upper bound Z is loose. A tighter physical upper bound (for example Za + D1, which pb_exp_pos uses) would
  presumably help heavy neutral atoms. pb_exp_pos's S2 result has now been seen, so any future model that adopts this
  idea must disclose that it was inspired by S2.
