# Track B: a generalized screened-hydrogenic formula for successive ionization energies

**Track B: data-driven / pattern recognition.** Code is in `models/semi_empirical/`, results in `results/se_*`, and figures in `results/figures/se_*.png`.

All numbers below were produced by the scripts named in each section and scored with the shared `evaluate.py` on the NIST ASD table `data/nist_ie.csv` (5847 successive IEs, $Z=1\dots110$, all charge states).

## 0. Summary

* **What it is:** a 32-parameter closed-form, screened-Dirac-hydrogenic formula. It has class-wise screening plus an ion-charge-dependent penetration term $\tau/(Z_a+\kappa)$, an exact Dirac factor with $l,j$-dependent penetration corrections, and a Hund-exchange kink term.
* **What it covers:** every successive IE of every element from configuration input alone, with no per-element parameters.
* **Accuracy on all 5847 NIST IEs:** **1.57 % MAPE, 0.74 % median**, 59 % of rows within 1 %.
  * Neutral-atom first IEs: 7.6 % MAPE (0.57 eV MAE).
  * H-like ions: 0.19 %.
* **Comparison with the best textbook baseline** (Slater total-energy difference): 11.8 % overall and 51.8 % on neutral atoms.
* **Held-out tests:**
  * S1 (unseen Z): 1.59 % test MAPE.
  * S3 (unseen isoelectronic sequences): 1.70 %.
  * S2 (train $Z\le54$, test $Z\ge55$): median 1.5 %, but MAPE 19 % because heavy near-neutral atoms fail. The f-electron and deep-core parameters cannot be learned from light atoms (Section 6).
* **Key discovery:** the expansion variable for penetration and relaxation is the *ionic* charge $Z_a=q+1$, not $Z$ or $Z_{\rm eff}$. Writing the second-order screening as $\tau_g/(Z_a+\kappa)$ with one universal $\kappa\approx8.2$ lowers the error from 9.3 % to 1.6 %.

## 1. Notation

| symbol | meaning |
|---|---|
| $Z$ | nuclear charge |
| $N$ | number of electrons **before** ionization; ion charge $q = Z-N$ |
| $Z_a = Z-N+1 = q+1$ | asymptotic charge felt by the departing electron (charge of the remaining ion) |
| $(n,l)$ | principal and orbital quantum numbers of the removed electron (subshell that loses an electron between the NIST ground configurations of the $N$- and $(N-1)$-electron ions) |
| $k$ | occupancy of the removed subshell in the $N$-electron configuration ($1\le k\le 4l+2$) |
| $\nu_c$ | number of **other** electrons (all electrons except the departing one) that belong to screening class $c$ (Table 2) |
| $\mathrm{Ry}$ | Rydberg energy, 13.605693 eV |
| $\alpha$ | fine-structure constant |
| $j$ | total angular momentum of the removed electron, assigned by $jj$ filling: the first $2l$ electrons of an $l$ subshell are $j=l-\tfrac12$, the rest $j=l+\tfrac12$ |
| $\kappa_j = j+\tfrac12$ | Dirac quantum number magnitude |
| MAPE | mean absolute percentage error $\frac{100}{M}\sum \lvert \mathrm{IE}_{\rm pred}-\mathrm{IE}_{\rm NIST}\rvert/\mathrm{IE}_{\rm NIST}$ |

## 2. Baselines: textbook formulas applied to successive IEs

Script: `models/semi_empirical/baselines.py`. Predictions are in `results/se_<name>_predictions.csv` and metrics in `results/se_<name>_metrics.json`. All rows use the same removed-subshell logic as the new model.

* **Bohr:** $\mathrm{IE}=\mathrm{Ry}\,Z^2/n^2$.
* **Slater (1930), one-electron:** $\mathrm{IE}=\mathrm{Ry}\,(Z-s)^2/n^{*2}$. This uses Slater's grouping (1s)(2s,2p)(3s,3p)(3d)(4s,4p)(4d)(4f)… and his constants 0.35 (0.30 for 1s), 0.85, 1.00. The effective quantum numbers are $n^*=1,2,3,3.7,4.0,4.2$, and I chose $4.3$ for $n=7$.
* **Slater, total-energy difference:** $\mathrm{IE}=E_S(N-1)-E_S(N)$ with $E_S=-\mathrm{Ry}\sum_i q_i (Z_{{\rm eff},i}/n^*_i)^2$.
* **Clementi–Raimondi (1963), both variants:** these use their screening formulas $\sigma(1s)\dots\sigma(4p)$ with the true $n$, as in their Slater-type orbitals. Their rules stop at 4p, so I extended them as follows (also documented in the module docstring):
  1. In the 1s/2s formulas, the "$n\ge3$" outer term counts *all* electrons with $n\ge 3$.
  2. A subshell beyond the CR range is mapped to the CR formula for the same $l$ with the largest available $n$ (s, p → 4s, 4p; d → 3d), shifting every principal quantum number by $\Delta=n-n_{\rm analog}$. Inner electrons that map below the CR range screen with 1.0.
  3. For f electrons, which have no CR formula, the same-subshell constant is 0.2693 (the CR d–d value), everything inside screens with 1.0, and outer electrons screen with 0.

**Table 1. Baselines on all 5847 IEs.**

| model | params | MAPE % | median APE % | neutral-atom first IE MAPE % (108) | H-like MAPE % (110) |
|---|---|---|---|---|---|
| Bohr $Z^2/n^2$ | 0 | 1362.9 | 171.7 | 22878 | 6.00 |
| Slater, one-electron | 0 (textbook) | 27.56 | 8.60 | 118.7 | 6.00 |
| Slater, $\Delta E$ total | 0 (textbook) | 11.77 | 7.49 | 51.8 | 6.00 |
| Clementi–Raimondi, one-electron | 0 (literature) | 31.21 | 11.27 | 180.0 | 6.00 |
| Clementi–Raimondi, $\Delta E$ total | 0 (literature) | 119.6 | 18.66 | 2281 | 6.00 |
| **GSHM core (this work)** | 30 | 1.618 | 0.795 | 7.71 | 0.19 |
| **GSHM final (this work)** | 32 | **1.573** | **0.737** | **7.61** | **0.19** |

Observations:

* Slater's total-energy difference is clearly better than the one-electron version (11.8 % vs 27.6 % overall). The difference automatically includes the relaxation of the remaining electrons, which are less screened after ionization.
* Clementi–Raimondi is worse than Slater for successive IEs, and catastrophically so in the total-energy form. The CR constants were fitted to *orbital exponents* of neutral atoms, and include small positive screening of inner orbitals by outer electrons (for example $0.0158$ per $n\ge3$ electron for 1s). Removing one valence electron changes every inner-shell energy by $\approx 2\,\mathrm{Ry}\,Z_{\rm eff}\,\Delta\sigma$. For a heavy atom that is tens of eV per inner electron, which swamps a 5–10 eV first IE.
* All of these baselines are exact only for hydrogen-like ions, and even there they miss relativity: H-like MAPE is 6 % and reaches 19 % at $Z=110$.

## 3. Pattern recognition in the NIST data

Script: `models/semi_empirical/explore.py`. Figures are listed in Section 9.

1. **Shell structure** (`se_explore_ie_vs_N.png`). At fixed $Z$, $\log \mathrm{IE}$ vs $N$ is a staircase that jumps at closed shells: $N=2,10,18,28,36,46,54,\dots$
2. **Moseley linearity along isoelectronic sequences** (`se_explore_moseley.png`, `results/se_moseley_fits.csv`). For every $N$, $y=n\sqrt{\mathrm{IE}/\mathrm{Ry}}$ (the "experimental $Z_{\rm eff}$") is almost exactly linear in $Z$ once the ion is a few times ionized. Fitting $y=aZ+b$ on $2\le q$, $Z\le 60$ gives:
   * a slope $a=1.006$ for 1s, 1.01–1.015 for $n=2$, 1.03 for 3s/3p, and 1.05–1.07 for 3d. The 3d slopes above 1 are a first sign of relativistic growth and of the $Z$-dependence of penetration.
   * an intercept $\sigma_N=-b/a$, which is the screening of the other $N-1$ electrons. $\sigma_N$ is a nearly piecewise-linear function of $N$ with slope $\approx 0.6$–$0.75$ inside a subshell (same-subshell screening) and $\approx 0.85$–$1$ when a new shell opens. This is the empirical basis of the class-wise constant screening $S_0=\sum_c\sigma_c\nu_c$.
3. **Penetration depends on the ion charge** (`se_explore_penetration.png`). The unscreened charge $p=Z_{\rm eff}-Z_a$ is not constant along a sequence: it grows with $q$, quickly at first and then slowly. For the Na sequence (N = 11, 3s), $p=0.84,1.15,1.34,1.46,1.56$ for $q=0\dots4$ and 2.11 at $q=20$. The pattern is the same for every $N$. The functional form $p(q)=p_\infty-\tau/(Z_a+\kappa)$ describes it with a common $\kappa\approx 8$. This was the single most important discovery of the analysis: it lowers the overall MAPE from ≈ 9 % to ≈ 2 % (Table 4).
4. **Half-filled-subshell kinks** (`se_explore_kinks.png`). $\mathrm{IE}/(q+1)^2$ along 2p, 3d and 4f filling shows the p³/p⁴ (N/O), d⁵/d⁶ and f⁷/f⁸ drops at every charge state. The drop has the sign and shape of the change in the number of parallel-spin pairs (Hund's rule) relative to a statistical average.
5. **4s/3d and 6s/5d/4f competition** (`se_explore_removal_map.png`). The removed subshell switches between $ns$ and $(n-1)d$ / $(n-2)f$ near neutrality. 63 ions rearrange on ionization. The model uses the NIST ground configurations, so the removed subshell is known.
6. **Relativity** (`se_explore_relativity.png`). Dividing NIST IEs by the non-relativistic Moseley line shows that the excess for s and p$_{1/2}$ electrons grows like the Dirac factor in $(Z\alpha)^2$, reaching +20–25 % at $Z\approx110$. d and f electrons show the opposite trend (−7 % to −20 %). That is the *indirect* relativistic effect: contracted inner s, p shells screen d and f electrons more. A single Dirac factor cannot describe both, so the model multiplies the Dirac factor by an $l,j$-dependent penetration correction.

## 4. The new model: Generalized Screened-Hydrogenic Model (GSHM)

Code: `models/semi_empirical/gshm.py`. The entry point is `predict(Z, N, shells=None)`, which returns eV. Parameters come from `results/se_params.json`.

### 4.1 Final formula

$$
\boxed{\;\mathrm{IE}(Z,N)=\Big[\underbrace{\mathrm{Ry}\,\frac{Z_{\rm eff}^2}{n^2}\,D_{n\kappa_j}(Z_{\rm eff})\,\Big(1+\frac{r_{lj}\,(Z\alpha)^2}{n}\Big(\frac{Z_{\rm eff}}{Z_a}-1\Big)\Big)}_{\text{screened Dirac hydrogen}}\;+\;\underbrace{\mathrm{Ry}\,x_l\,K_l(k)\,\frac{Z_{\rm eff}}{n^2}}_{\text{exchange / Hund}}\Big]\;\exp\Big(\sum_{m}\beta_m\phi_m\Big)\;}
$$

**Effective charge (two-level screening):**

$$
Z_{\rm eff}^{(0)} = Z-\sum_{c}\sigma_c\,\nu_c ,\qquad
Z_{\rm eff} = Z_{\rm eff}^{(0)} - \frac{\tau_{\rm same}\,\nu_{\rm same}+\tau_{\rm near}\,\nu_{\rm near}+\tau_{\rm core}\,\nu_{\rm core}}{Z_a+\kappa}.
$$

* The first-order screening constants $\sigma_c$ are charge-independent. They give the Moseley intercepts.
* The second term is a **charge-dependent screening (penetration) correction**. The other electrons are grouped into "same" (same subshell), "near" (same shell or the shell just below) and "core" (deeper) electrons. Their combined screening is enhanced by $\tau_g/(Z_a+\kappa)$, which is large near neutrality ($Z_a=1$) and vanishes for highly charged ions ($Z_a\to\infty$, where $Z_{\rm eff}\to Z_{\rm eff}^{(0)}$, the exact $1/Z$ screening limit).
* Physically, the departing electron of a weakly charged ion has a diffuse orbital that sees the core more completely shielded (less penetration). This is the $1/Z$-expansion structure $E=Z^2E_0+ZE_1+E_2+E_3/Z+\dots$ resummed into a Padé-like form in the *ionic* charge.

**Relativistic factor.** $D$ is the exact ratio of the Dirac point-nucleus binding energy to the Schrödinger one for a hydrogenic level of charge $\zeta=Z_{\rm eff}$:

$$
D_{n\kappa}(\zeta)=\frac{2n^2}{(\zeta\alpha)^2}\left[1-\left(1+\frac{(\zeta\alpha)^2}{\big(n-\kappa+\sqrt{\kappa^2-(\zeta\alpha)^2}\big)^2}\right)^{-1/2}\right],\qquad D\to 1+\frac{(\zeta\alpha)^2}{n^2}\Big(\frac{n}{\kappa}-\frac34\Big).
$$

It is exact for H-like ions, which have no free parameter. For many-electron ions, the factor $1+\frac{r_{lj}(Z\alpha)^2}{n}(Z_{\rm eff}/Z_a-1)$ accounts for the penetrating part of the orbital feeling the full nuclear relativistic potential (direct effect, $r>0$) or the relativistically contracted core (indirect effect, $r<0$ for d, f). It vanishes for H-like ions ($Z_{\rm eff}=Z_a$) and for $Z\alpha\to0$. There are five classes of the removed electron: s, p$_{1/2}$, p$_{3/2}$, d, f.

**Exchange (Hund) term.** Let $P_l(k)$ be the number of parallel-spin pairs in $l^k$ under Hund's first rule. With $m=2l+1$, $u=\min(k,m)$ and $d=k-u$, $P_l(k)=\binom{u}{2}+\binom{d}{2}$. The kink function is the exchange stabilization lost on removing one electron, relative to the statistical (configuration-average) value:

$$K_l(k)=\big[P_l(k)-P_l(k-1)\big]-\frac{2l\,(k-1)}{4l+1}.$$

For p this gives $k=1\dots6$: $0,\,0.6,\,1.2,\,-1.2,\,-0.6,\,0$. It is positive up to half filling (IE raised) and negative just after it (p⁴, d⁶, f⁸: IE lowered), which is exactly the N/O, Cr/Mn/Fe and Eu/Gd kink pattern. $x_l\,\mathrm{Ry}\,Z_{\rm eff}/n^2$ is the corresponding exchange integral per pair: it scales like $\langle 1/r\rangle\propto Z_{\rm eff}/n^2$. For s, $K=0$.

**Sparse correction.** $\exp(\sum_m\beta_m\phi_m)$ is a small multiplicative correction. The terms $\phi_m$ were chosen by cross-validated forward selection from a library of 33 physically motivated candidates (Section 5). Selected terms: $\phi_1=(Z\alpha)^2\,[l=2]$ and $\phi_2=(Z\alpha)^2\,[\mathrm{p}_{3/2}]$, with $\beta_1=-0.0247$ and $\beta_2=-0.0163$. These are small extra relativistic corrections for d and p$_{3/2}$ electrons (the spin-orbit-raised level is slightly less bound than the jj-assigned Dirac factor implies).

### 4.2 Screening classes

**Table 2.** The class of each other electron relative to the removed electron $(n,l)$. Ground-state data cannot distinguish the full inner s,p shells seen by a d or f electron (their populations are always the same constants), so those shells are merged.

| class $c$ | removed electron | screening electrons | S1 group |
|---|---|---|---|
| same_s / same_p / same_d / same_f | any | other electrons of the same subshell | same |
| sn_in_p | p | s electrons of the same $n$ | near |
| n1_sp_sp, n1_sp_d, n1_sp_f | s or p | shell $n-1$: s,p / d / f electrons | near |
| n2_sp_sp, n2_sp_df | s or p | shell $n-2$: s,p / d,f electrons | core |
| d_near | d | s,p electrons of shells $n$ and $n-1$ | near |
| n1_d_d, n1_d_f | d | $(n-1)$d / $(n-1)$f electrons | near |
| f_near | f | s,p,d electrons of shells $n$ and $n-1$ | near |
| n1_f_f | f | $(n-1)$f electrons (4f seen by 5f) | near |
| core | any | all deeper electrons ($n'\le n-3$ for s,p; $n'\le n-2$ for d,f) | core |
| out_d, out_f | d, f | electrons with $n'>n$ (4s seen by 3d; 5s,5p,6s seen by 4f…) | – (enter $Z_a$) |

### 4.3 Parameters

**Table 3.** Final parameters, refitted on all 5847 IEs (`results/se_params.json`, block `final`). There are **32 fitted global parameters**: 18 screening constants, 3 tau plus kappa, 5 relativistic, 3 exchange and 2 sparse-correction coefficients. There are no per-element or per-ion parameters. The uncertainty is the formal 1-sigma from the Jacobian.

| parameter | value | 1 sigma | meaning |
|---|---|---|---|
| $\sigma_{\rm same\_s}$ | 0.5908 | 0.0084 | other electron in the same s subshell |
| $\sigma_{\rm same\_p}$ | 0.6956 | 0.0053 | same p subshell |
| $\sigma_{\rm same\_d}$ | 0.7914 | 0.0043 | same d subshell |
| $\sigma_{\rm same\_f}$ | 0.8800 | 0.0044 | same f subshell |
| $\sigma_{\rm sn\_in\_p}$ | 0.7361 | 0.0051 | p target: s electrons of the same shell |
| $\sigma_{\rm n1\_sp\_sp}$ | 0.7110 | 0.0022 | s/p target: s,p electrons of shell n-1 |
| $\sigma_{\rm n1\_sp\_d}$ | 0.7595 | 0.0021 | s/p target: d electrons of shell n-1 |
| $\sigma_{\rm n1\_sp\_f}$ | 0.7394 | 0.0023 | s/p target: f electrons of shell n-1 |
| $\sigma_{\rm n2\_sp\_sp}$ | 0.8182 | 0.0029 | s/p target: s,p electrons of shell n-2 |
| $\sigma_{\rm n2\_sp\_df}$ | 0.8335 | 0.0020 | s/p target: d,f electrons of shell n-2 |
| $\sigma_{\rm d\_near}$ | 0.7945 | 0.0014 | d target: s,p electrons of shells n, n-1 |
| $\sigma_{\rm n1\_d\_d}$ | 0.7684 | 0.0019 | d target: (n-1)d electrons |
| $\sigma_{\rm n1\_d\_f}$ | 0.8224 | 0.0022 | d target: (n-1)f electrons (4f14 under 5d) |
| $\sigma_{\rm f\_near}$ | 0.8644 | 0.0011 | f target: s,p,d electrons of shells n, n-1 |
| $\sigma_{\rm n1\_f\_f}$ | 0.8282 | 0.0032 | f target: (n-1)f electrons |
| $\sigma_{\rm core}$ | 0.8616 | 0.0019 | deep core (s,p: n' <= n-3; d,f: n' <= n-2); bounded to [0.8, 1] |
| $\sigma_{\rm out\_d}$ | 0.4497 | 0.0135 | d target: electrons in higher shells (e.g. 4s for 3d) |
| $\sigma_{\rm out\_f}$ | 0.7437 | 0.0041 | f target: electrons in higher shells (5s5p6s for 4f) |
| $\tau_{\rm same}$ | 0.6721 | 0.0514 | charge-dependent extra screening, same-subshell electrons |
| $\tau_{\rm near}$ | 1.6737 | 0.0266 | charge-dependent extra screening, "near" electrons |
| $\tau_{\rm core}$ | 1.1751 | 0.0283 | charge-dependent extra screening, core electrons |
| $\kappa$ | 8.2155 | 0.1134 | offset in the $1/(Z_a+\kappa)$ penetration law |
| $r_{\rm s}$ | 1.7068 | 0.0969 | relativistic penetration coefficient, s |
| $r_{\rm p1}$ | -0.0265 | 0.0711 | same, p$_{1/2}$ |
| $r_{\rm p3}$ | -0.5888 | 0.0490 | same, p$_{3/2}$ |
| $r_{\rm d}$ | 0.3592 | 0.0874 | same, d (indirect effect) |
| $r_{\rm f}$ | 0.1228 | 0.0901 | same, f |
| $x_p$ | 0.1294 | 0.0127 | Hund exchange per parallel pair, p |
| $x_d$ | 0.1021 | 0.0097 | same, d |
| $x_f$ | 0.0629 | 0.0103 | same, f |
| $\beta_{\rm x2\_d}$ | -0.0247 | 0.0034 | sparse correction: $\phi=(Z\alpha)^2[l=2]$ |
| $\beta_{\rm x2\_p3}$ | -0.0163 | 0.0030 | sparse correction: $\phi=(Z\alpha)^2[\mathrm{p}_{3/2}]$ |

Fixed, not fitted: $\mathrm{Ry}$ and $\alpha$ (CODATA 2018), the Slater-like prior values for the ridge, the bound $0.8\le\sigma_{\rm core}\le1$, and the jj and Hund rules.

**How the parameters were fitted:** ridge-regularized nonlinear least squares on $\ln(\mathrm{IE}_{\rm pred}/\mathrm{IE}_{\rm NIST})$, $\chi^2=\sum_i\ln^2(\cdot)+\lambda^2\lVert\theta-\theta_0\rVert^2$ with $\lambda=0.1$ and $\theta_0$ the Slater-like start values ($\sigma$: 0.3 to 1.0, $\tau=0$, $\kappa=7$, $r=x=\beta=0$). $\lambda$ was chosen by an inner extrapolation check (train $Z\le44$, validate $45\le Z\le54$): $\lambda=0$ gives 2.19 % and $\lambda=0.1$ gives 1.31 % on that check. The ridge leaves the all-data fit unchanged (1.617 to 1.618 %).

## 5. Discovery process and ablation

All variants were fitted with `scipy.optimize.least_squares` (TRF, physical bounds $0\le\sigma_c\le1.3$, $0\le\sigma_{\rm out}\le1$, $0\le\kappa\le100$; the final model also uses ridge $\lambda=0.1$ and $0.8\le\sigma_{\rm core}\le1$) on the residuals $\ln(\mathrm{IE}_{\rm pred}/\mathrm{IE}_{\rm NIST})$, so relative errors are weighted equally from 4 eV to 200 keV. The sequence of discoveries was:

1. **Class-wise constant screening only** (Slater-like, 18 classes). MAPE 10.9 %. The residuals were strongly structured in the ion charge $q$: too low near neutrality, too high for intermediate $q$.
2. **Relativity:** the Dirac factor plus five penetration coefficients. This brings H-like ions from 6 % to 0.19 % (no parameter acts on them) and the overall MAPE to ≈ 9 %.
3. **Charge-dependent screening $\tau_g/(Z_a+\kappa)$:**
   * With $Z_{\rm eff}^{(0)}$ in the denominator, a "Padé in $Z$", the MAPE is 3.5 %.
   * With the *ionic charge* $Z_a=q+1$ in the denominator, it is 1.66 %.
   * The asymptotic charge is the correct expansion variable for the outermost electron. This choice alone halves the error.
4. **Hund exchange kinks:** three parameters. These fix the p³/p⁴ and d⁵/d⁶ steps.
5. **Sparse residual corrections:** forward selection (below).

Rejected or neutral variants. These were measured during exploration on all data, without the final ridge and bound, and were not re-run:

| variant | params | MAPE % | neutral MAPE % |
|---|---|---|---|
| quantum defect $n^*=n-\delta_l/(1+(Z^{(0)}_{\rm eff}-1)/k)$ on top of step 3 (old denominator) | 35 | 3.02 | 16.9 |
| per-class $\tau_c$ (20 instead of 3) | 45 | 1.58 | 7.5 |
| indirect relativistic screening $\rho_l(Z\alpha)^2\nu_{\rm in}$ | 34 | 1.62 | 7.0 |
| second order $\tau^{(2)}/(Z_a+\kappa)^2$ | 33 | 1.54 | 8.0 |
| $l$-dependent $\kappa_l$ | 33 | 1.65 | 7.2 |
| "penetration-reduction" form $Z_{\rm eff}=Z^{(0)}_{\rm eff}-\sum_g P_g\kappa_g/(Z_a+\kappa_g)$ | 29 | 2.79 | 15.0 |
| total-energy-difference form $\sum_{N-1}q\,b-\sum_N q\,b$ with the same ingredients (1st-order $Z_{\rm eff}^{(0)}$ in S1) | 32 | 2.92 | 15.3 |

The extra parameters in these variants did not buy enough to be kept within the 40-parameter budget. Notably, the **total-energy-difference** form is *worse* than the one-electron form once the charge-dependent screening is present (2.9 % vs 1.6 %). The relaxation it adds is already represented by the $Z_a$-dependent screening, and summing over all inner orbitals amplifies errors in inner-shell energies. Its fit is also ~100× slower.

**Table 4. Ablation** (`results/se_validation.json`, key `ablation`). Each row is refitted on all data. The "S1 test" column is the held-out MAPE on $Z\equiv0 \pmod 5$ after refitting on the S1 training set.

| model | params | MAPE % (all) | median % | neutral first IE MAPE % | S1 test MAPE % |
|---|---|---|---|---|---|
| A0 screening only (sigma_c) | 18 | 11.904 | 8.620 | 66.67 | 11.722 |
| A1 + q-dependent screening S1 | 22 | 3.196 | 1.695 | 8.93 | 3.224 |
| A2 + Dirac factor (no penetration corr.) | 22 | 1.942 | 0.908 | 8.48 | 1.934 |
| A3 + relativistic penetration r_lj | 27 | 1.630 | 0.773 | 7.94 | 1.634 |
| A4 + Hund exchange x_l (= core GSHM) | 30 | 1.618 | 0.795 | 7.71 | 1.633 |
| A5 + sparse corrections (= final GSHM) | 32 | 1.573 | 0.737 | 7.61 | 1.587 |
| -- core minus S1 | 26 | 9.341 | 6.777 | 48.10 | 9.211 |
| -- core minus relativity | 25 | 3.156 | 1.676 | 8.75 | 3.182 |
| -- core minus exchange | 27 | 1.630 | 0.773 | 7.94 | 1.634 |
| -- core with S1 denominator Z_eff0 (not Z_a) | 30 | 3.207 | 1.617 | 18.60 | 3.192 |

What each term buys:

* The **charge-dependent screening** is the dominant ingredient: removing it raises the MAPE from 1.62 to 9.34 %.
* The ionic charge $Z_a$ is the right variable in its denominator: using $Z^{(0)}_{\rm eff}$ instead gives 3.21 %.
* **Relativity** is worth a factor of 2 (3.16 to 1.62 %).
* The **exchange term** barely changes the global MAPE (1.630 to 1.618 %), because kinks affect few rows. It reduces the neutral-atom error from 7.94 to 7.71 % and is needed for the N/O and Cr/Mn/Fe shapes.
* The **sparse corrections** gain 0.045 points overall and 0.05 on S1.

**Forward selection of correction terms.**

* Library: 33 terms (`gshm.corr_features`): $l$-resolved $1/Z_a$ and $1/Z_a^2$; class-resolved $(Z\alpha)^2$; $(Z\alpha)^4$; kink$\times1/Z_a$; half-filled, full and single-electron indicators $\times1/Z_a$; rearrangement indicators; d$^{10}$/f$^{14}$ core indicators $\times1/Z_a$; outer-electron count $\times1/Z_a$; $N/Z$; $\ln Z$; a QED-like $\frac{\alpha}{\pi}(Z\alpha)^2\ln(Z\alpha)^{-2}/n$ for s; $(n\ge5)/Z_a$; same-subshell occupancy $/Z_a$.
* Selection set: the S1 training records. Inner validation: the $Z\equiv2 \pmod 5$ subset of them, so the S1 test set never influences selection.
* Each step screens every candidate by linear least squares of its coefficient on the current log-residuals, keeps the candidate with the lowest inner-validation MAPE, and then jointly refits the whole model. It stops when the gain is < 1 %, or after 6 terms.

Path (`results/figures/se_gshm_selection_path.png`, a parameters-vs-held-out-error learning curve):

| step | terms | params | inner-validation MAPE % |
|---|---|---|---|
| 0 | (core) | 30 | 1.6447 |
| 1 | x2_d | 31 | 1.6228 |
| 2 | x2_d, x2_p3 | 32 | 1.6060 |

The next-best candidate, $(Z\alpha)^4$, gave < 1 % relative gain, so selection stopped at 32 parameters. Held-out error is flat from 30 to 32 parameters (S1 test 1.633 to 1.587 %), and train and test errors agree within 0.02 points on S1 and S3. There is **no sign of over-fitting** for interpolation.

## 6. Held-out validation

The standard splits were used with fits on train only:

* **S1**, interpolation in $Z$: test $Z\equiv0 \pmod 5$.
* **S2**, extrapolation: train $Z\le54$, test $Z\ge55$.
* **S3**, unseen isoelectronic sequences: test $N\equiv0 \pmod 6$.

**Table 5.** Test-set metrics. "core" is the GSHM without the sparse corrections (fully clean). "final" includes them.

| split | model | n test | test MAPE % | median APE % | p90 % | neutral first-IE MAPE % (n) |
|---|---|---|---|---|---|---|
| S1 | core | 1188 | 1.633 | 0.784 | 3.97 | 6.89 (21) |
| S1 | final | 1188 | 1.587 | 0.728 | 3.79 | 6.70 (21) |
| S2 | core | 4362 | 19.234 | 1.534 | 39.39 | 259.80 (54) |
| S2 | final | 4362 | 18.712 | 3.600 | 37.78 | 234.81 (54) |
| S3 | core | 928 | 1.707 | 0.794 | 4.56 | 8.36 (18) |
| S3 | final | 928 | 1.700 | 0.801 | 4.43 | 8.12 (18) |

**S2 (extrapolation from $Z\le54$ to $Z\ge55$) is the honest failure mode.** S1 and S3 test errors equal the training error, so the model interpolates in Z and generalizes to unseen isoelectronic sequences. S2 is different. Its *median* error is small (1.5 %), but its mean is dominated by heavy neutral and near-neutral atoms (7s, 6d and 5f electrons), where the error reaches hundreds of %. These are the diagnosed causes and the measured history (core model, S2 test):

| S2 variant | test MAPE % | median % | neutral MAPE % |
|---|---|---|---|
| no constraints | 850.97 | 2.39 | 35499.1 |
| + f parameters tied to their d analogues | 390.73 | 2.16 | 19380.0 |
| + ridge lambda = 0.1 | 178.50 | 1.41 | 7953.1 |
| + bound 0.8 <= sigma_core <= 1 (final) | 19.23 | 1.53 | 259.8 |

1. $Z\le54$ contains no f electrons, so eight parameters ($\sigma$ of the f classes, $r_f$, $x_f$) have zero support. In every split they are set equal to their d analogue (`gshm.ANALOG`, applied automatically by `fitlib.auto_tie`; only S2 triggers it).
2. In $Z\le54$ the deep-core class always has the same electron count for a given target subshell (2 or 10). $\sigma_{\rm core}$ is therefore not identifiable and drifts to 0.3, and every heavy atom then gets a huge $Z_{\rm eff}$.
3. The bound $0.8\le\sigma_{\rm core}\le1$ encodes the physics that deep shells screen almost fully (Slater: 1.00). **Disclosure:** I added this bound *after* seeing variant 3 fail on S2, so the final S2 row is not fully blind. The bound changes neither the all-data fit nor S1 or S3.
4. The two sparse $(Z\alpha)^2$ corrections, fitted on $Z\le54$ where $(Z\alpha)^2<0.16$, extrapolate poorly (the S2 median rises from 1.5 to 3.6 %). The core model is the better extrapolator.

**Caveat on "final":** the *choice* of correction terms was made on the S1 training pool, which contains part of the S2 and S3 test sets. The coefficients are always refitted on each split's training set only. The "core" rows have no such issue.

## 7. Final model on all data (refit on all 5847 IEs)

Predictions are in `results/se_gshm_predictions.csv` (final) and `results/se_gshm_core_predictions.csv`, with metrics in `results/se_gshm_metrics.json` and `results/se_gshm_core_metrics.json`. Output of `evaluate.py`:

```
=== gshm ===
coverage: 5847/5847
stratum                              n    MAPE%  medAPE%     p90%     max%     MAE eV   <=1%
11<=N<=36                         2275    1.164    0.625    2.880    35.72     27.648   0.64
19<=Z<=54                         1314    1.700    0.906    3.813    35.72     11.844   0.53
ALL                               5847    1.573    0.737    3.929    35.72     36.436   0.59
N<=10                             1048    0.698    0.332    1.029    21.57    112.559   0.90
N>=37                             2524    2.305    1.307    5.386    35.27     12.749   0.41
Z<=18                              171    3.092    1.903    7.709    21.57      3.641   0.37
Z>=55                             4362    1.476    0.683    3.841    25.60     45.129   0.61
first_IE_neutral_atoms             108    7.612    6.443   16.358    35.72      0.572   0.09
hydrogen_like                      110    0.190    0.118    0.501     0.91    233.372   1.00
rearranged_config                   63    4.712    2.756   12.467    17.44     42.632   0.27
removed_l=d                       1929    0.982    0.509    2.151    35.27     10.386   0.74
removed_l=f                        794    2.675    1.849    6.081    17.44     20.572   0.28
removed_l=p                       2066    1.544    0.758    3.791    35.72     27.650   0.58
removed_l=s                       1058    1.881    0.789    4.717    20.31    112.991   0.56
removed_n=1                        219    0.392    0.308    0.885     1.49    386.226   0.93
removed_n=2                        829    0.779    0.337    1.324    21.57     40.263   0.89
removed_n=3                       1632    0.828    0.422    1.838    23.44     21.868   0.76
removed_n=4                       1917    1.868    1.106    4.105    35.72     26.403   0.46
removed_n=5                        986    2.412    1.502    4.793    25.60      7.994   0.34
removed_n=6                        224    4.243    3.642    8.472    16.36      3.862   0.16
removed_n=7                         40    5.117    4.854    9.426    20.17      0.680   0.15
status=experimental                311    5.415    3.950   12.130    35.72      1.921   0.14
status=semi-empirical              919    2.146    1.181    4.877    25.60     11.796   0.42
status=theoretical                4617    1.200    0.607    3.220    20.17     43.665   0.65
```

The core model (30 parameters) scores 1.618 % MAPE, 0.795 % median and 7.71 % on neutral first IEs (`results/se_gshm_core_metrics.json`).

### Hard cases

* **H-like ions:** 0.19 % MAPE, with no fitted parameter acting on them, because the Dirac factor is exact for a point nucleus. The remaining error, up to 0.9 % at $Z=110$, is finite nuclear size plus QED, which are not modelled. For example U$^{91+}$: 131816 eV (NIST) vs 132280 eV (+0.35 %).
* **Highly charged heavy ions** are excellent:
  * U with N = 3, 10, 11 and 47: +0.48 %, -0.14 %, -0.01 % and +0.04 %.
  * Overall, $Z\ge55$ scores 1.48 % MAPE and 0.68 % median.
* **Neutral-atom first IEs** (108): 7.6 % MAPE, 6.4 % median, with an absolute MAE of only 0.57 eV. That compares with 51.8 % for Slater total-energy difference and 22878 % for Bohr. Examples (NIST vs GSHM):

  | | He | Li | N | O | Ne | Na | Ar | K | Cr | Fe |
  |---|---|---|---|---|---|---|---|---|---|---|
  | NIST (eV) | 24.59 | 5.39 | 14.53 | 13.62 | 21.56 | 5.14 | 15.76 | 4.34 | 6.77 | 7.90 |
  | GSHM (eV) | 24.31 | 5.04 | 12.57 | 13.59 | 21.94 | 5.93 | 16.36 | 4.67 | 6.00 | 8.03 |
  | error | -1.1 % | -6.6 % | -13.5 % | -0.2 % | +1.7 % | +15.4 % | +3.8 % | +7.5 % | -11.3 % | +1.6 % |

  | | Cu | Kr | Rb | Pd | Xe | Cs | La | Gd | Au | Rn | Fr | U |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | NIST (eV) | 7.73 | 14.00 | 4.18 | 8.34 | 12.13 | 3.89 | 5.58 | 6.15 | 9.23 | 10.75 | 4.07 | 6.19 |
  | GSHM (eV) | 7.56 | 14.87 | 4.52 | 9.15 | 11.38 | 3.78 | 4.96 | 6.04 | 8.60 | 9.41 | 3.96 | 5.46 |
  | error | -2.1 % | +6.2 % | +8.2 % | +9.8 % | -6.2 % | -2.8 % | -11.1 % | -1.7 % | -6.8 % | -12.5 % | -2.9 % | -11.8 % |

  The alkali minima and noble-gas maxima are reproduced qualitatively everywhere (`se_gshm_overview.png`, right panel). The N/O ordering is not: the exchange kink is too weak for 2p$^3$ at $q=0$.
* **Rearranged configurations** (63): 4.7 % MAPE.
* **f electrons:** 2.7 % MAPE. **n = 6 and n = 7 electrons:** 4.2 % and 5.1 %. The weakest strata are the near-neutral lanthanides and actinides.

## 8. Limitations

* **Neutral atoms and low charge states remain the weak spot.** Their error is several times the overall MAPE. Some of the largest errors are structural:
  * Ga, In and Tl, the first p electron after a filled d¹⁰ shell.
  * B, C, Al: penetration of the first p electron.
  * Y, Nb, Mo and other 4d atoms, where the configuration-average screening ignores the s/d energy competition.
  * The q = 0–2 lanthanide and actinide ions.

  The model does not predict *which* configuration is the ground state. It takes the NIST ground configurations as input (Madelung outside the table). With user-supplied `shells`, the removed electron is the outermost one (largest $n$, then $l$). That rule reproduces the NIST removed subshell for 94.4 % of the table, so for nf / (n+1)s competition it may pick the wrong electron.
* **Term-averaged energies only.** No fine-structure or multiplet splitting beyond the single Hund kink term. The NIST IEs refer to the lowest *levels*, so a residual error of a few tenths of an eV near neutrality is irreducible in this form.
* **Several heavy, highly charged values in NIST are themselves theoretical** (Dirac–Fock + QED) with stated uncertainties of up to ~1–2 % for some actinide and superheavy ions. Some of the largest "errors" (for example the $N=69$ sequence for $Z\ge104$, where NIST's ground configuration 4f¹⁴5d¹ loses a 4f electron) are probably configuration-assignment issues rather than model failures.
* **Parameter uncertainties are formal.** They come from the Jacobian with heavy-tailed log-residuals, and some parameters are strongly correlated (for example the $\tau_g$ with $\kappa$, and $\sigma_{\rm core}$ with $\tau_{\rm core}$). Individual values should not be over-interpreted.
* **Skipped for time** (orchestrator request to wrap up):
  * per-split term selection, so the "final" S2/S3 rows have mild selection leakage, and "core" is the clean reference;
  * a finer learning curve (only the selection path was produced);
  * a dedicated per-case analysis of lanthanides and actinides;
  * combining the total-energy form with the $Z_a$-denominator screening (only tested with the older variant);
  * bootstrapped parameter uncertainties;
  * extra figures (relativistic residual vs $(Z\alpha)^2$ for the final model).
* **QED** is not modelled explicitly beyond whatever the $(Z\alpha)^2$ terms absorb. Neither is the finite nuclear size (relevant for 1s at $Z>90$).

## 9. Files

* Code (`models/semi_empirical/`):
  * `configs.py`: removed subshell, Hund pairs, $jj$ filling.
  * `baselines.py`: Bohr, Slater and CR.
  * `explore.py`: figures and Moseley fits.
  * `gshm.py`: the model and `predict()`.
  * `fitlib.py`: fitting and splits.
  * `fit_gshm.py`: selection, validation, ablation and final fit.
* Parameters: `results/se_params.json`, which has `final` and `core` blocks with names, values, formal 1σ and the spec.
* Validation and ablation: `results/se_validation.json`. The fit log is `results/se_fit_log.txt`.
* Figures (`results/figures/`):
  * `se_explore_ie_vs_N.png`
  * `se_explore_moseley.png`
  * `se_explore_penetration.png`
  * `se_explore_kinks.png`
  * `se_explore_neutral.png`
  * `se_explore_relativity.png`
  * `se_explore_removal_map.png`
  * `se_gshm_overview.png`
  * `se_gshm_selection_path.png`

### For the unification stage

* `gshm.predict(Z, N, shells=None, params=None)` takes `params=(spec, theta)`, so Track A's exact leading coefficients can be imposed. For example, fix $\sigma_{\rm same\_s}$ to the exact $1/Z$ two-electron value.
* `gshm.dirac_factor(n, j, Zeff)` is the parameter-free relativistic factor.
* `gshm.CLASSES` and `gshm.orbital_rows()` give the screening-class occupation vectors for any configuration.
* `fitlib.fit(spec, recmask)` fits any spec on a training mask.
