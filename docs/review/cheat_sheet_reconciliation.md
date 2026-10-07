# Cheat sheet, case studies and benchmark matrix vs the code (checklist items 3, 4, 9)

The user pasted the cheat sheet's contents on 2026-10-06. The code values below come from `tools/audit_components.py`
(`results/audit_components.md`) and from a sweep of every model and split refit in the repo (scratch script, nothing
written to `results/`).

## 1. Which model is the cheat sheet? `pa_bound9`

The cheat sheet's nine parameters are the `pa_bound9` all-data fit in `results/pa_params.json`, rounded:

| cheat sheet | `pa_bound9` |
|---|---|
| τ_same 0.235 | 0.2347 |
| τ_in 0.597 | 0.5971 |
| τ_core 2.956 | 2.9563 |
| τ_df 2.951 | 2.9515 |
| τ_out 23.71 | 23.7087 |
| κ 2.26 (2.21 + 0.05) | 2.2091 (+0.05) |
| x_p 0.360 | 0.3602 |
| x_d 0.424 | 0.4237 |
| x_f 0.821 | 0.8215 |

`pa_bound9` is **not** the paper's 8-parameter "pocket formula" (Appendix A, `uni_pocket`).
- `pa_bound9` uses the pocket formula's five electron groups, but it is the final model's equation (exact σ₁ +
  bounded remainder D + Dirac factor + Hund term + μ) with 9 parameters. It has no class deviations and no
  relativistic bracket.
- The pocket formula has no σ₁ and no bound. It uses Z_eff = Z − Σ s_g ν_g − t(N−1)/(Z_a+κ) and a Sommerfeld bracket.

The cheat sheet must not call `pa_bound9` a "pocket" model. Its parameters belong in the paper next to Table A1.

## 2. Line by line

| case | cheat sheet | `pa_bound9` code | verdict |
|---|---|---|---|
| Na 3s | σ₁ 7.786, ν_in 8, ν_core 2, T 10.69, h 2.214, Z_eff 1.893, IE 5.41 (NIST 5.14) | σ₁ 7.78564, ν_in 8, ν_core 2, T 10.68923, h 2.21436, Z_eff 1.89247, IE 5.4144 (+5.36 %) | **matches** |
| Ca 4s | Z_eff 2.668, IE 6.29 (NIST 6.113) | ν_same 1, ν_in 8, ν_core 10; T 34.5745; Z_eff 2.25478; IE **4.3234 (−29.3 %)** | **does not match** |
| Fe 4s | ν_df 6, Z_eff 3.045, IE 8.28 (NIST 7.90) | ν_same 1, ν_in 14 (3s²3p⁶ **and 3d⁶**), ν_core 10, ν_df 0; T 38.157; Z_eff 2.92962; IE **7.2989 (−7.6 %)** | **does not match** |

**Fe grouping error.** The cheat sheet puts the 3d⁶ electrons in ν_df. In the model the df group exists only for d
and f **targets**. For an s or p target, every electron of shell n − 1, whatever its l, is in the **in** group
(class `n1_sp_d`). With the cheat sheet's own grouping (ν_in 8, ν_df 6) and `pa_bound9` parameters you get
T = 52.28, D = 4.268, Z_eff = 2.547 and IE ≈ 5.5 eV. So the printed Z_eff 3.045 doesn't follow from its own counts
either.

## 3. The "+0.24 eV (Ca) and +0.39 eV (Fe) above the bare Rydberg term"

Using the cheat sheet's own Z_eff:
- Ca: Ry·Z_eff²/n² = 13.6057 · 2.668²/16 = 6.053 eV, so 6.29 − 6.053 = **+0.237 eV**.
- Fe: 13.6057 · 3.045²/16 = 7.884 eV, so 8.28 − 7.884 = **+0.395 eV**.

In `pa_bound9`, and in the paper's equation, a 4s electron gets nothing on top of Ry·Z_eff²/n² except:
- the Dirac factor, +0.006 % to +0.009 % (≈ +0.0004 to +0.0007 eV);
- μ(Z), about −10⁻⁵.

The Hund term is zero for s electrons (K_0 = 0, and there is no x_s). There is no QED or FNS term for n = 4. So
**no term in `pa_bound9` or in the paper produces +0.24 or +0.39 eV.** The reviewer is right that the paper has no
s-type term that could explain them.

### Sweep of every model in the repo: IE (Z_eff) [IE − Ry·Z_eff²/n²]

| model | Ca I | Fe I |
|---|---|---|
| **cheat sheet** | 6.290 (2.668) [+0.237] | 8.280 (3.045) [+0.395] |
| `pa_bound9` (all-data; V1/V2/S1/S3/S2 refits all within ±0.4 eV and [+0.000]) | 4.323 (2.255) [+0.000] | 7.299 (2.930) [+0.001] |
| `pa_bound9_fs0` (parameter-free Fermi–Segrè term on s, p) | 4.035 (2.166) [+0.046] | 6.991 (2.831) [+0.174] |
| final `pa_hier_rel` (r_s < 0 lowers s IEs) | 5.027 (2.436) [−0.017] | 7.611 (3.004) [−0.062] |
| `pa_hier_rel` refits (V1, V2, S1, S3, S2) | 4.97–5.71, [−0.017 to −0.184] | 7.49–7.61, [−0.06 to −0.55] |
| `pa_hier` | 5.019 (2.429) [+0.000] | 7.502 (2.970) [+0.001] |
| `pa_bound14_relfit` / `pa_nobound14_relfit` | 4.596 / 5.402, [−0.018 / −0.017] | 7.662 / 9.062, [−0.076 / −0.071] |
| `pb_clip_pos` | 6.458 (2.753) [+0.012] | 8.331 (3.124) [+0.030] |
| u35 / u29 | 6.711 / 6.865, [+0.019 / −0.004] | 8.231 / 8.664, [+0.045 / −0.010] |
| pocket formula | 4.925 (2.407) [+0.000] | 8.386 (3.140) [+0.001] |
| **Track B GSHM final / core** | **6.137 (2.666) [+0.092]** / 6.117 (2.663) [+0.089] | **8.028 (3.026) [+0.243]** / 8.033 (3.028) [+0.236] |

**Answer.** No model reproduces both lines. The closest is the **Track B GSHM** (32 parameters, `docs/semi_empirical.md`):
- its Ca Z_eff 2.666 is within 0.002 of the cheat sheet's 2.668;
- it is the only model with a sizeable *positive* term for an s electron above the Rydberg term, from its fitted
  relativistic-penetration bracket 1 + r_s(Zα)²(Z_eff/Z_a − 1)/n with r_s > 0 (the final model's fitted r_s is
  negative, −0.454);
- that term is worth +0.09 eV (Ca) and +0.24 eV (Fe), not +0.24 and +0.39;
- its IEs, 6.137 and 8.028, are not the cheat sheet's 6.29 and 8.28.

So the Ca and Fe lines look like a hand-made mix: a GSHM-like Z_eff plus an extra term that matches no code path, under
`pa_bound9`'s parameter header. They must not be published. The true `pa_bound9` values are Ca 4.32 eV (−29.3 %) and
Fe 7.30 eV (−7.6 %). Alkaline-earth and other 4s² neutrals are a known weak spot of the bounded models; the paper's
§5.2 should say so. The unbounded models (GSHM, u35, `pb_clip_pos`) do better on Ca but extrapolate catastrophically
(blind S2 57–170 %).

## 4. Pb: "σ₁ = 76.2" vs about 63 in Figure 3

- The code's σ₁ for Pb I (…6s²6p², 6p removed) is **63.724**, the same as Figure 3.
- The quantities near 76 are:
  - the total screening σ₁ + D: 76.87 (final), 77.04 (`pa_bound9`);
  - the pocket formula's total screening, 76.60;
  - Slater's σ, 76.35;
  - Figure 3's "experimental" σ = Z − n√(IE/Ry) = 77.57.
- The case study almost certainly printed the **total** screening, or Slater's σ, under the label σ₁.
- `pa_bound9` gives Pb I 9.37 eV (+26.3 %); the final model gives 8.01 eV (+8.0 %).

## 5. Benchmark matrix (33 p / 9 p / Mendoza) vs the paper: row sets

| matrix row | as pasted | same-row values (n) | problem |
|---|---|---|---|
| all ions | 1.87 / 2.90 / – | 1.874 / 2.897 (5847 rows; this is all rows, neutrals included) | label it "all rows" |
| Mendoza-covered | 1.62 / 2.63 / 2.82 | 1.623 / 2.634 / 2.822 (5011) | consistent |
| neutral 1st IE | 7.5 / 12.0 / 27.1 | **8.37 / 11.45 / 27.09 (54 covered neutrals)**; 7.52 / 11.99 are on all 108 neutrals | **mixed row sets**: 108 rows for 33 p and 9 p, 54 rows for Mendoza |
| experimental subset | 4.4 / 7.2 / 11.3 | 4.374 / 7.169 / 11.276 (236 covered experimental rows) | consistent in itself, but the paper prints 4.63 % (33 p) on all **311** experimental rows; 9 p on 311 rows is 7.54 % |
| N ≤ 10 | 0.44 / 0.51 / 0.40 | 0.439 / 0.514 / 0.396 (1048, all covered) | consistent; **Mendoza is the most accurate of the three here**, and the paper must say so |
| one-electron | 0.001 / 0.001 / 0.19 | 0.0012 / 0.0012 / 0.190 (110) | consistent |

**Rule for the paper.**
- Every cell states its n.
- Cross-model comparisons use the same rows.
- 4.63 % (311 rows) and 4.37 % (236 Mendoza-covered rows) are both correct, but on different rows; print n beside each.
- The Mendoza 11.3 % on experimental rows coincides with the pocket formula's blind S2 of 11.3 %. Never print either
  without its label.
