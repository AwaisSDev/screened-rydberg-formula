# Screening groups ν_g and classes ν_c: one written rule (checklist item 5)

Date: 2026-10-06. Checker: `tools/check_grouping_rule.py`. It implements the rule below from its own table and calls
`gshm.classify` / `PA._group_matrix` / `PA._dev_matrix` / `pocket.counts` only to compare. Environment:
`py -3.11` with `PYTHONPATH=tools/numba_stub`, run time about 95 s. Nothing is written to disk.

**Result.** The rule reproduces the code exactly, with **0 mismatches**:
- all 5847 NIST rows (n = 5847);
- all coverage cases Z = 1–118, N = 1–Z (n = 7021);
- both the 5 group counts ν_g and the 21 class counts ν_c.

It gives O, Na, N, Ca, Fe and Pb exactly as in `results/audit_components.json`. The same five groups serve the final
model (pa_hier_rel), the 9-parameter bounded variant (pa_bound9) and the 8-parameter pocket formula.

What the current §2.2 text gets wrong or leaves out (details in §6):
1. It does not say where same-n, higher-l electrons go. The code puts them in *out*.
2. It writes "all n′ ≤ n" for d/f targets. That wrongly includes same-n, higher-l electrons.
3. It never defines the 19 class names used in Table A1.

None of the three changes a number, because the affected cases occur in 0 of the 5847 rows.

---

## 1. Paper-ready text: replaces the ν_g and ν_c bullets of §2.2

The two bullets "**ν_g** counts the other electrons …" and "**ν_c and δτ_c** are finer screening classes …" are
replaced by these three:

```markdown
- **ν_g (five screening groups).** Let (n, l) be the removed subshell. Each of the other N − 1 electrons, in a
  subshell (n′, l′), is counted in exactly one group. *same*: the k − 1 other electrons of (n, l). *out*: outer
  electrons, n′ > n, or n′ = n with l′ > l. All other electrons (n′ < n, or n′ = n with l′ < l) are inner. Inner
  electrons of an s or p target form *in* if n′ ≥ n − 1 (the ns electrons of an np target and the whole (n − 1)
  shell, including its d and f electrons) and *core* if n′ ≤ n − 2. Inner electrons of a d or f target, from the
  same-n lower-l subshells down to 1s, all form *df*. Hence Σ_g ν_g = N − 1. Examples, (ν_same, ν_in, ν_core, ν_df,
  ν_out): O 2p⁴: (3, 4, 0, 0, 0); Na 3s: (0, 8, 2, 0, 0); Fe 3d⁶4s², 4s removed: (1, 14, 10, 0, 0), since 3d⁶ is in
  the (n − 1) shell; Pb 6p²: (1, 20, 60, 0, 0), with 6s²5s²5p⁶5d¹⁰ in *in* and all n′ ≤ 4, including 4f¹⁴, in
  *core*. Each group has a coefficient τ_g.
- **ν_c and δτ_c (screening classes, Table A2).** The class of an electron depends only on the target type (s/p, d
  or f), on n − n′ and on l′. Each class lies inside one group, and an electron of class c in group g contributes
  τ_g + δτ_c to T. The deviations are shrunk toward their group value by a ridge penalty (10⁻⁴ per row). Two of the
  21 classes, sn_out (n′ = n, l′ > l) and out_sp (n′ > n, s/p target), are empty in all 5847 rows and carry no
  δτ_c, which leaves the 19 of Table A1. The 9-parameter variant pa_bound9 and the pocket formula (Appendix A) use
  the same five groups without classes.
- **Removed subshell.** (n, l) is the subshell whose occupancy drops from the N-electron configuration (NIST
  ground configuration; Madelung order if the ion is not tabulated) to the NIST ground configuration of the
  (N − 1)-electron ion. For a rearranging ion (e.g. V 3d³4s² → V⁺ 3d⁴) it is the subshell that loses most electrons
  (ties: larger n, then larger l). If the (N − 1)-electron ion is not tabulated, or the configuration is supplied
  by the user, the outermost subshell (largest n, then largest l) is removed.
```

Optional sentence for the first bullet, to be checked against [Slater1930] before use (not checked against the
1930 paper in this session): "The partition has the structure of Slater's rules (the (n − 1) shell separate from
deeper shells for s/p electrons; all inner electrons together for d/f electrons), except that the ns electrons of an
np target are counted with the (n − 1) shell."

## 2. Paper-ready Table A2 (goes directly under Table A1)

```markdown
**Table A2.** Screening classes c and groups g of §2.2. (n, l) is the removed subshell and (n′, l′) the subshell of
another electron. Each class lies inside one group. "rows" is the number of the 5847 NIST rows with ν_c > 0.
δτ_c is from Table A1, and τ_g + δτ_c is the coefficient of one electron of that class in T.

| group g | class c | target l | other electron (n′, l′) | rows | δτ_c | τ_g + δτ_c |
|---|---|---|---|---|---|---|
| same | same_s | s | (n′, l′) = (n, l) | 534 | −0.1552 | 0.2376 |
| same | same_p | p | (n′, l′) = (n, l) | 1713 | 0.0098 | 0.4026 |
| same | same_d | d | (n′, l′) = (n, l) | 1720 | 0.0004 | 0.3932 |
| same | same_f | f | (n′, l′) = (n, l) | 717 | 0.1457 | 0.5385 |
| in | sn_in_p | p | n′ = n, l′ = s | 2073 | −0.1133 | 0.6745 |
| in | n1_sp_sp | s, p | n′ = n − 1, l′ = s, p | 2912 | −0.3916 | 0.3963 |
| in | n1_sp_d | s, p | n′ = n − 1, l′ = d | 1248 | 0.3009 | 1.0888 |
| in | n1_sp_f | s, p | n′ = n − 1, l′ = f | 386 | 0.2042 | 0.9921 |
| core | n2_sp_sp | s, p | n′ = n − 2, l′ = s, p | 2083 | −0.2195 | 2.0048 |
| core | n2_sp_df | s, p | n′ = n − 2, l′ = d, f | 668 | −0.0622 | 2.1621 |
| core | core_sp | s, p | n′ ≤ n − 3 | 1311 | 0.2828 | 2.5070 |
| df | d_near | d | n′ = n or n − 1, l′ = s, p | 1938 | −0.6962 | 1.8526 |
| df | n1_d_d | d | n′ = n − 1, l′ = d | 1078 | −0.2689 | 2.2799 |
| df | n1_d_f | d | n′ = n − 1, l′ = f | 391 | 0.3286 | 2.8773 |
| df | f_near | f | n′ = n or n − 1, l′ = s, p, d | 778 | 0.5476 | 3.0964 |
| df | n1_f_f | f | n′ = n − 1, l′ = f | 201 | −0.0918 | 2.4569 |
| df | core_df | d, f | n′ ≤ n − 2 | 2716 | 0.1812 | 2.7300 |
| out | out_d | d | n′ > n | 15 | −0.0176 | 10.5512 |
| out | out_f | f | n′ > n | 294 | 0.0173 | 10.5861 |
| out | sn_out | any | n′ = n, l′ > l | 0 | none | 10.5688 |
| out | out_sp | s, p | n′ > n | 0 | none | 10.5688 |
```

The τ_g and δτ_c values were read from `results/pa_params.json` (candidate `pa_hier_rel`) through `PA._resolve`. All
19 δτ_c and all 5 τ_g agree with the values printed in Table A1 to 4 decimals. The τ_g + δτ_c column is their sum,
rounded to 4 decimals.

### 2a. Smaller edits that go with the rule (every old text checked to occur exactly once in docs/paper_draft.md)

- **§2.2, K_l(k) bullet** (values verified in §5 below). The whole bullet becomes:

  ```markdown
  - **K_l(k) = [P(k) − P(k−1)] − 2l(k−1)/(4l+1)** is the change in the number of parallel-spin pairs relative to a
    statistical average (Hund kink); P(k) is the number of parallel-spin pairs of l^k under Hund's first rule.
    Equivalently K_l(k) = (2l+1)(k−1)/(4l+1) for k ≤ 2l+1 and K_l(k) = −K_l(4l+3−k) above half filling, so K_s ≡ 0,
    K_p = 0, 0.6, 1.2, −1.2, −0.6, 0 for k = 1…6, K_d = (5/9)(0, 1, 2, 3, 4, −4, −3, −2, −1, 0) and
    K_f = (7/13)(0, 1, …, 6, −6, …, −1, 0). Its amplitude x_l is fitted for l = p, d, f; s electrons have no Hund
    term.
  ```
- **Table A1.** Change the row label "| δτ_c (19 classes) |" to "| δτ_c (19 classes, defined in Table A2) |".
- **Appendix A, pocket formula.** Change "The fitted values are s_same 0.7499" to "The ν_g are the five group
  counts of §2.2 (Table A2). The fitted values are s_same 0.7499".
- **§3.1.** Change "beyond Z = 110 the Madelung order is used." to "for ions not in the table the Madelung order is
  used. These are all ions with Z > 110 and 258 ions with Z = 104–110. They enter no fit or score, only the
  7021-configuration checks (§2.2, §4.7) and predictions outside the table." (see §7.2).

## 3. Verification (`tools/check_grouping_rule.py`)

**Inputs.** The configuration and removed subshell come from `configs.initial_final(Z, N)`, which is the code's own
rule. Every count is then made by the independent rule: Table A2 as a list of conditions, plus a separate prose-form
`group_of()`. For each electron the checker asserts:
- exactly one table row applies;
- the class lies in the group that `group_of()` gives.

**Code side of the comparison.** `umodel.build(ZN)` → `PA._group_matrix({"grouping": "pocket"}, a)` (ν_g),
`PA._dev_matrix(a)` (ν_c, 21 columns, core split into core_sp / core_df) and `pocket.counts(a)`.
`abinitio.save_cache` is disabled inside the checker, so `models/unified/cache/sigma1.json` is never rewritten
(git status confirmed it unchanged).

| check | n | mismatches |
|---|---|---|
| Partition: every (target n, l; other n′, l′) cell, with n, n′ = 1–8, l′ < n′ and l, l′ ≤ 3, matched by exactly one row of Table A2; class → group consistent | 676 cells | 0 |
| NIST rows: removed subshell of rule input vs built features | 5847 | 0 |
| NIST rows: ν_g (5 groups) vs `PA._group_matrix` | 5847 | 0 |
| NIST rows: ν_c (21 classes) vs `PA._dev_matrix` | 5847 | 0 |
| NIST rows: `PA._group_matrix` vs `pocket.counts` (code vs code) | 5847 | 0 |
| NIST rows: Σ_g ν_g ≠ N − 1 | 5847 | 0 |
| Coverage Z = 1–118, N = 1–Z: removed subshell, ν_g, ν_c, pocket.counts, Σ_g ν_g = N − 1 | 7021 | 0 each |
| Classes populated in the 5847 rows vs the pa_hier_rel `spec["dev"]` list | 19 vs 19 | identical sets |

Empty classes are sn_out and out_sp, both with 0 rows in the NIST table and 0 rows in the coverage set.

**Targets by l.**
- NIST: s 1058, p 2073, d 1938, f 778 (sum 5847).
- Coverage: s 1887, p 2418, d 1938, f 778 (sum 7021).

**Edge cases (rows in which each occurs).**

| edge case | where the rule puts it (class) | NIST rows (n = 5847) | coverage rows (n = 7021) |
|---|---|---|---|
| same n, higher l (any target) | out (sn_out) | 0 | 0 |
| s target with same-n p electrons | out (sn_out) | 0 | 0 |
| p target with same-n d electrons | out (sn_out) | 0 | 0 |
| d target with same-n f electrons | out (sn_out) | 0 | 0 |
| d/f target with same-n lower-l electrons | df (d_near / f_near) | 2716 | 2716 |
| f target with (n − 1)d electrons | df (f_near; there is no separate n1_f_d class) | 778 (every f target) | 778 |
| deep core n′ ≤ n − 2, s/p target | core (n2_sp_sp / n2_sp_df at n − 2; core_sp below) | 2083 | 3170 |
| deep core n′ ≤ n − 2, d/f target | df (core_df) | 2716 | 2716 |
| n′ > n, d/f target | out (out_d 15 rows, out_f 294) | 309 | 309 |
| n′ > n, s/p target | out (out_sp) | 0 | 0 |
| rearranged per `configs.initial_final` | (removed-subshell rule) | 40 | 58 |

The 15 out_d rows are d removal with outer s electrons:
- Y I (4d with 5s²);
- Lu I, Hf I (5d with 6s²);
- Ac I, Th I, Pa I, U I, Cm I (6d with 7s²);
- 7 rows of Rf–Hs, one of which, Hs N = 106, has 7s¹.

Electrons per group in the NIST rows: same 18 792, in 38 890, core 35 241, df 110 073, out 1719.

## 4. The six case studies (neutral first IEs)

Configurations are the NIST ground configurations; k is the occupancy of the removed subshell. ν_g agrees with
`results/audit_components.json` for all three models (pa_hier_rel, pa_bound9, pocket). Nonzero ν_c agrees with the
`class_devs` of pa_hier_rel.

| ion | configuration (outer part) | removed (k) | ν_g (same, in, core, df, out) | nonzero ν_c |
|---|---|---|---|---|
| O | 1s² 2s² 2p⁴ | 2p (4) | 3, 4, 0, 0, 0 | same_p 3, sn_in_p 2, n1_sp_sp 2 |
| Na | [Ne] 3s¹ | 3s (1) | 0, 8, 2, 0, 0 | n1_sp_sp 8, n2_sp_sp 2 |
| N | 1s² 2s² 2p³ | 2p (3) | 2, 4, 0, 0, 0 | same_p 2, sn_in_p 2, n1_sp_sp 2 |
| Ca | [Ar] 4s² | 4s (2) | 1, 8, 10, 0, 0 | same_s 1, n1_sp_sp 8, n2_sp_sp 8, core_sp 2 |
| Fe | [Ar] 3d⁶ 4s² | 4s (2) | 1, 14, 10, 0, 0 | same_s 1, n1_sp_sp 8, n1_sp_d 6, n2_sp_sp 8, core_sp 2 |
| Pb | [Xe] 4f¹⁴ 5d¹⁰ 6s² 6p² | 6p (2) | 1, 20, 60, 0, 0 | same_p 1, sn_in_p 2, n1_sp_sp 8, n1_sp_d 10, n2_sp_sp 8, n2_sp_df 24, core_sp 28 |

How the Pb counts arise:
- **in** (20): 6s² (sn_in_p) + 5s²5p⁶ (n1_sp_sp) + 5d¹⁰ (n1_sp_d).
- **core** (60): 4s²4p⁶ (n2_sp_sp) + 4d¹⁰4f¹⁴ (n2_sp_df) + all n′ ≤ 3 (core_sp, 28).

How the Fe counts arise:
- The 4s target puts 3d⁶ into **in**, as n1_sp_d. It does not go into df: df applies only to d/f *targets*.

None of these six ions has a nonzero ν_df or ν_out. All six remove an s or p electron.

## 5. Hund kink K_l(k) and jj assignment of j (paper §2.2 and §2.4)

The checker recomputes P(k) independently in exact rational arithmetic. P(k) is the number of parallel-spin pairs of
l^k when the 2l + 1 spin-up orbitals are filled first. It then compares with `configs.exchange_kink` and
`configs.jj_j` for l = 0–3 and every k = 1 … 4l + 2.

**Result: 0 mismatches (32 (l, k) pairs: 2 + 6 + 10 + 14).**
- The paper's p values 0, 0.6, 1.2, −1.2, −0.6, 0 are correct.
- The statement "j = l − ½ while k ≤ 2l, otherwise l + ½" is correct. It gives j = ½ for s.
- The relativistic class p½/p3/2 (`gshm.rel_class`) is consistent with j.

| l | K_l(k), k = 1 … 4l+2 | j, k = 1 … 4l+2 |
|---|---|---|
| s | 0, 0 | ½, ½ |
| p | 0, 3/5, 6/5, −6/5, −3/5, 0 | ½, ½, 3/2, 3/2, 3/2, 3/2 |
| d | 0, 5/9, 10/9, 5/3, 20/9, −20/9, −5/3, −10/9, −5/9, 0 | 3/2 (k ≤ 4), 5/2 (k ≥ 5) |
| f | 0, 7/13, 14/13, 21/13, 28/13, 35/13, 42/13, −42/13, …, −7/13, 0 | 5/2 (k ≤ 6), 7/2 (k ≥ 7) |

The checker also verified this closed form for all 32 (l, k) pairs:

K_l(k) = (2l+1)(k−1)/(4l+1) for k ≤ 2l+1, and K_l(k) = −K_l(4l+3−k) above half filling.

The 2l(k−1)/(4l+1) term is the change in the statistical (configuration-average) number of parallel pairs. That
number is C(k,2)·2l/(4l+1), because 2l/(4l+1) is the fraction of spin-orbital pairs with equal spin.

K_s ≡ 0, so no model has an s-type Hund term. This bears on checklist item 3: Ca and Fe both remove 4s, so their
Hund term is exactly zero in all three models.

## 6. Current paper text vs code

1. **Same n, higher l.** §2.2 now reads "all n′ ≤ n for d/f targets; outer n′ > n". In the code (`gshm.classify`
   → `sn_out`) every electron with n′ = n and l′ > l is *outer*, whatever the target type. "n′ ≤ n" therefore
   wrongly includes, for example, 4f electrons on a 4d target. The case also covers np electrons on an ns target,
   which the old description left unassigned. It occurs in 0 of 5847 rows and 0 of 7021 coverage cases, so no
   number changes.
2. **Class names.** Table A1 lists 19 δτ_c names that the paper never defines. The code has 21 classes: the 20 of
   `gshm.CLASSES`, with `core` split by target type into core_sp / core_df. Two of them (sn_out, out_sp) are empty
   in the data. `PA.fit` keeps only classes with at least one training row (`min_rows` = 1), which is why 19 are
   fitted on all data.

   In held-out fits the class list is whatever is populated in the training rows. With the rule, the S2 training
   rows (Z ≤ 54) populate only 13 classes. Six f-type classes are empty there: same_f, n1_sp_f, n1_d_f, f_near,
   n1_f_f and out_f. A Z ≤ 54 fit therefore has 13 δτ_c rather than 19. This follows from `PA.fit`; the S2 fit
   itself was not re-run here.
3. **Code comments with the same imprecision.** Both are in shared or frozen code and were not changed:
   - `models/unified/pocket.py` docstring: "nu_df : d,f target: all electrons with n' <= n (other subshells)".
     The counts are right; the comment is not.
   - The `configs.jj_j` docstring says it "Returns (j, kappa_sign)", but the function returns j only.

## 7. Side findings outside item 5 (for the owners of §3.1)

1. **"63 rearranged ions" (§3.1).** 63 is the `rearranged` flag of `common/atomdata.py`, which is also the
   `rearranged_config` stratum of `evaluate.py`. 23 of those 63 are Z = 104–110 rows whose (N − 1)-electron ion is
   missing from the NIST table. For them, `atomdata._removed_subshell` compares with an empty configuration, which
   marks them "rearranged" and returns a meaningless `removed` subshell. The breakdown by Z is Z = 104–108: 3 rows
   each, Z = 109 and 110: 4 rows each. One example is Rf N = 4, 1s²2s², flagged because Rf N = 3 is not tabulated.

   The model's own rule (`configs.initial_final`) flags 40 rows. Every one of them is also among the 63: 40 are
   common to both rules, 23 are atomdata-only and 0 are configs-only. V I → V II (3d³4s² → 3d⁴) and Ni I → Ni II
   are among the 40, both with 4s removed. All 23 atomdata-only rows have the (N − 1)-electron ion missing. In 16
   rows (all Z ≥ 104) the `removed` field of atomdata differs from the subshell the model removes. The 63-row
   stratum therefore contains 23 rows that are not rearrangements. `evaluate.py` and `common/` are shared
   infrastructure and were not touched. The paper should either say this or quote 40.
2. **Madelung configurations (§3.1).** "beyond Z = 110 the Madelung order is used" is incomplete. 258 (Z, N) with
   Z = 104–110 are also absent from the table: Z = 104–108 have 35 each, Z = 109 has 41 and Z = 110 has 42. The
   coverage test therefore uses Madelung configurations for 258 + 916 = 1174 of its 7021 cases (7021 = 5847 + 258 +
   916). In 25 of these the N-electron ion is Madelung but the (N − 1)-electron ion is tabulated (N = 3, 47, 79,
   103 or 106; Z = 104–110). The removed subshell then comes from diffing a Madelung against a NIST configuration,
   and the code flags 18 of the 25 as "rearranged", for example Rf N = 47 and N = 79. That is why the coverage set
   has 58 rearranged cases (40 + 18). The 5847 fitted rows all use NIST configurations, so no fitted or scored
   number is affected.

## 8. Reproduce

```bash
export PYTHONPATH="D:/Chem-Research/Np/tools/numba_stub"   # only where numba is missing
py -3.11 tools/check_grouping_rule.py                      # or py -3.13 where available
```
