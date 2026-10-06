# A screened Rydberg formula for the successive ionization energies of all atoms and ions

**[Author name(s), affiliation(s), corresponding e-mail: to be supplied by the author]**

*Manuscript draft, 2026-10-05. Built from the project results files at the time of writing. Numbers may be
updated by the final audit/fix pass; every number below is traceable to a file in `results/` or `docs/`
(see the list of sources returned with this draft). Items marked **[TODO]** or **[VERIFY]** must be resolved
before submission.*

---

## Abstract

We present a single closed-form expression, the *screened Rydberg formula*, for the successive ionization energy
IE(Z, N) of any atom or ion. It needs only the nuclear charge Z and the ground electron configuration. The formula
keeps the hydrogenic form Ry Z_eff²/n². Its effective charge is split into two parts. The first is an exact,
parameter-free first-order screening constant σ₁ = −n²ΔE₁ from the 1/Z perturbation expansion; it is a rational
number for every configuration, and we tabulate it for N = 1–110. The second is a bounded, fitted higher-order
remainder that depends on the ion charge in an Edlén-type form. A saturation construction guarantees
Z − N + 1 ≤ Z_eff ≤ Z for any parameter values (given 0 ≤ σ₁ ≤ N − 1, which holds for every configuration we
use) and keeps the exact Z² and Z coefficients of the 1/Z series. One-electron
Dirac, recoil, finite-nuclear-size and QED corrections are added without fitted parameters. The final model has 33
fitted global parameters and no per-element or per-ion parameters. On all 5847 successive ionization energies in the
NIST Atomic Spectra Database (Z = 1–110), its mean absolute percentage error (MAPE) is 1.87 % (median 0.94 %). That is
7.52 % for neutral-atom first ionization energies and 0.0012 % for hydrogen-like ions; the H-like figure is consistent
with, but not independent of, the QED theory behind the NIST values. The model was chosen by a pre-registered score
that never used the blind extrapolation test (its interpolation splits do contain heavy elements). In a single blind
extrapolation (fit Z ≤ 54, predict Z = 55–110), it gives 6.60 % MAPE
(median 1.62 %). The same test gives 11.6 % for Slater's rules, 13.6 % for a published parameter-free
screened-hydrogenic model scored on the same rows, and 57–178 % for earlier fitted variants without the bound.
Extrapolation to heavy neutral atoms remains the weak point (21.8 % MAPE on their first ionization energies).

**Keywords:** ionization energy; screening constants; 1/Z expansion; screened hydrogenic model; isoelectronic
sequences; Slater's rules; blind validation.

---

## 1. Introduction

The Bohr energy Ry Z²/n² is exact only for a non-relativistic one-electron ion. For every other atom and ion the
ionization energy (IE) is controlled by electron–electron repulsion, which has no closed-form solution for N ≥ 2. Two
broad routes exist.

**Screening rules.** Slater [Slater1930] replaced Z by an effective charge Z − σ, with σ given by simple counting
rules fitted by hand to atomic data. Clementi and Raimondi [ClementiRaimondi1963; ClementiRaimondiReinhardt1967]
obtained σ from optimised self-consistent-field orbital exponents of neutral atoms. Both are designed for orbitals of
neutral atoms, not for the full set of successive IEs. Scored as IE formulas on all NIST ions with a total-energy
difference, Slater's rules give 11.8 % MAPE and Clementi–Raimondi 120 % (§4).

**Perturbation theory in 1/Z.** Layzer [Layzer1959; LayzerEtAl1964; Layzer1967] showed that the non-relativistic
energy of a fixed configuration is an asymptotic series E = Z²E₀ + ZE₁ + E₂ + …. E₀ is hydrogenic and E₁ is a
rational combination of hydrogenic Slater integrals, which fixes the exact Z → ∞ limit of the screening constant.
Higher orders need continuum sums [ScherrKnight1963; DalgarnoStewart1958]. Along isoelectronic sequences, Edlén
[Edlen1964] systematised the smooth variation of screening with ion charge, using expansions in 1/(ζ + s). Z-expansion
codes with relativistic corrections were developed for selected sequences [Safronova1993].

**Screened hydrogenic models (SHMs).** Starting with Mayer [Mayer1947] and More [More1982], plasma-physics codes use
screened hydrogenic energies with tabulated or fitted screening constants: l-splitting [Faussurier1997;
Faussurier2008], analytical potentials [Martel1998; Rubiano2002], relativistic constants fitted by genetic algorithm
[Mendoza2011], and a parameter-free, self-consistent (Z, N)-dependent screening [Kregar1984; Kregar1985; DiRocco1992;
Pomarico2005; LanziniDiRocco2015; DiRoccoLanzini2016]. Average-atom and kinetics codes value these for speed and
coverage of all ions [Rozsnyai1972; Chung2005; Crilly2023].

**Ab initio methods.** Koopmans' theorem [Koopmans1934], ΔSCF Hartree–Fock and Kohn–Sham DFT [KohnSham1965;
Kotochigova1997], correlated non-relativistic energies [Chakravorty1993] and Dirac–Fock total energies for all
ground configurations up to Z = 118 [Rodrigues2004] give more accurate numbers. They are numerical procedures, not
formulas.

**The gap.** We found no single closed-form expression that (i) covers every ion of every element from the
configuration alone, (ii) reproduces the exact large-Z behaviour of the 1/Z expansion, (iii) reports every fitted
parameter, and (iv) is validated on all NIST successive IEs with held-out tests fixed in advance, including a blind
extrapolation to heavier elements. The SHM papers we read report accuracy on their fitting data. For example, about
88 % of the database of [Mendoza2011] lies within ±10 %, and [Pomarico2005] reports selected sequences. This paper
tries to fill that gap. It does not claim a new law of atomic physics; §5.1 states what is rediscovered.

**Contributions.**
1. A table of the exact first-order screening constants σ₁ for every NIST ground configuration, N = 1–110, as a
   parameter-free, Z → ∞-exact counterpart to Slater's σ (§2.1, §4.6).
2. The screened Rydberg formula: σ₁ plus a bounded Edlén-type remainder with 33 global parameters (§2.2–2.4).
3. A validation protocol with a pre-registered selection score and a single blind extrapolation test, applied
   identically to the final model, its ancestors and published baselines (§3, §4).
4. A head-to-head on the same 5847 rows against a published parameter-free SHM, which we re-implemented from its
   definitions (§4.5).

---

## 2. Theory

Atomic units are used in §2.1, with energies converted to eV via Ry = 13.6057 eV. Z is the nuclear charge, N the
number of electrons before ionization, and Z_a = Z − N + 1 the charge seen by the departing electron. (n, l) is the
removed subshell and k its occupancy.

### 2.1 Exact first-order screening from the 1/Z expansion

Scaling r → ρ/Z turns the non-relativistic Hamiltonian into H = Z²[H₀ + Z⁻¹V], with H₀ hydrogenic and
V = Σ 1/ρ_ij. Rayleigh–Schrödinger perturbation theory in 1/Z then gives, for a fixed configuration and term,
E(Z, N) = Z²E₀ + ZE₁ + E₂ + …, and the ionization energy is

$$\mathrm{IE}(Z,N)=Z^2\Delta E_0+Z\,\Delta E_1+\Delta E_2+\dots,\qquad \Delta E_k=E_k(N-1)-E_k(N).$$

- **Zeroth order.** ΔE₀ = 1/(2n²), which is exactly Bohr's formula.
- **First order.** E₁ = ⟨V⟩ over the Z = 1 hydrogenic state. Every hydrogenic radial integral F^k and G^k is an exact
  rational number, because each product P_a P_c is a polynomial times e^{−βr} with rational β. Examples:
  F⁰(1s,1s) = 5/8, F⁰(1s,2s) = 17/81, G⁰(1s,2s) = 16/729. We evaluate them in exact rational arithmetic and combine
  them with exact angular algebra.
- **Open shells.** For open subshells, E₁ is evaluated for the Hund ground term by projecting onto the
  highest-weight (S, L) subspace. Where hydrogenic degeneracy couples configurations (the Layzer complex, e.g.
  1s²2s² with 1s²2p²), V is diagonalised within the complex. For two or more open subshells, a high-spin coupled term
  is used; this is the one approximation in σ₁.

Completing the square in the first two terms defines the first-order screening constant

$$\mathrm{IE}\simeq \Delta E_0\,(Z-\sigma_1)^2,\qquad \sigma_1=-\frac{\Delta E_1}{2\Delta E_0}=-n^2\Delta E_1 .$$

σ₁ is Layzer's Z → ∞ screening constant. It has no adjustable parameter and depends only on the configuration, not on
Z. Two checks against known values:

- **He-like:** E₁(1s²) = 5/8, so σ₁(He) = 0.6250 (Slater: 0.30).
- **Li-like:** E₁(1s²2s) = 5965/5832 = 1.0228052 = 5/8 + 2·17/81 − 16/729, so ΔE₁ = −0.397805 and σ₁ = 1.5912.
- **Be-like:** the complex value E₁ = 1.5592742 agrees with Layzer's.

The exact series fails for near-neutral atoms (§4.6). The expansion parameter is effectively N/Z. For a neutral atom,
Z − σ is only 1–3, so a 10 % error in σ becomes a 100–1000 % error in IE. A closed form for the higher orders does
not exist: E₂ already requires sums over the hydrogenic continuum. For example, the completed-square guess
ΔE₂ = ΔE₁²/(4ΔE₀) is 24–26 % too large for He- and Li-like ions. This is why the remainder in §2.2 has to be fitted.

### 2.2 The screened Rydberg form with a bounded remainder

The final model ("pa_hier_rel") is

$$
\mathrm{IE}(Z,N)=\mu(Z)\Big\{\mathrm{Ry}\,\frac{Z_\mathrm{eff}^2}{n^2}\,D_{n,j}(Z_\mathrm{eff})\Big[1+r_c\,\frac{(Z\alpha)^2}{n}\Big(\frac{Z_\mathrm{eff}}{Z_a}-1\Big)\Big]+\mathrm{Ry}\,\frac{x_l\,K_l(k)}{n^2}\Big\}-\big[\Delta E_\mathrm{QED}+\Delta E_\mathrm{FNS}\big]_{Z,n}\Big(\frac{Z_\mathrm{eff}}{Z}\Big)^2 ,
$$

where the last term applies only to the removal of an ns electron with n ≤ 2, and

$$
Z_\mathrm{eff}=Z-\sigma_1(\mathcal C)-D,\qquad
D=\frac{T}{Z_a+\kappa+|T|/h},\qquad
h=\begin{cases}(N-1)-\sigma_1 & T\ge 0\\ \sigma_1 & T<0\end{cases},\qquad
T=\sum_{g}\tau_g\,\nu_g+\sum_{c}\delta\tau_c\,\nu_c .
$$

The terms are:

- **ν_g** counts the other electrons in five screening groups (same subshell; inner, i.e. same-n s and the (n−1)
  shell, for s/p targets; deeper core for s/p targets; all n′ ≤ n for d/f targets; outer n′ > n). Each group has a
  coefficient τ_g.
- **ν_c and δτ_c** are finer screening classes, with deviations shrunk toward their group value by a ridge penalty
  (10⁻⁴ per row).
- **K_l(k) = [P(k) − P(k−1)] − 2l(k−1)/(4l+1)** is the change in the number of parallel-spin pairs relative to a
  statistical average (Hund kink). For p electrons it is 0, 0.6, 1.2, −1.2, −0.6, 0 for k = 1…6. Its amplitude x_l
  is fitted per l.
- **μ(Z)** is the reduced-mass factor.

The saturation form of D has two properties that hold for any parameter values.

1. **The bound Z_a ≤ Z_eff ≤ Z.** Write D = sign(T)·h·u/(1+u) with u = |T|/[h(Z_a+κ)]. Then |D| < h, so Z_eff can
   neither fall below the fully screened, non-penetrating limit Z_a nor exceed the bare charge Z. Since
   h = (N−1) − σ₁ or σ₁, this requires 0 ≤ σ₁ ≤ N − 1. That is a checked property of the configurations, not a
   proven one: it holds for all 7021 ground/Madelung configurations with Z ≤ 118, with no violations, and the code
   warns if a user-supplied configuration breaks it.
2. **The exact asymptotics.** At fixed N and Z_a → ∞, D → T/(Z_a+κ) = O(1/Z). The expansion of IE therefore begins
   Ry[Z² − 2Zσ₁ + O(1)]/n², and its Z² and Z coefficients are those of the exact 1/Z series. The fitted part only
   represents ΔE₂ and higher orders: relaxation, correlation and penetration at low ion charge.

### 2.3 Charge-dependent penetration as an Edlén-type term

The remainder's dependence on 1/(Z_a + κ) was found empirically (docs/semi_empirical.md). Along every isoelectronic
sequence, the excess charge p = Z_eff − Z_a grows with ion charge q as p∞ − τ/(Z_a + κ), with one κ common to all
sequences. For the Na sequence (3s), for example, p = 0.84, 1.15, 1.34, 1.46, 1.56 for q = 0–4 and 2.11 at q = 20.
Introducing this term reduced the error of the purely fitted model from about 9 % to about 2 % MAPE. This is an
independent rediscovery of the isoelectronic regularity that Edlén formalised [Edlen1964]. It is also consistent with
Layzer's theory, in which the screening constant is σ₀ + σ₁′/Z + …. We do not present it as a new law. What is
specific here is narrower: one universal denominator across all sequences, combined with exact σ₁ and with the
saturation bound of §2.2.

### 2.4 Relativistic, recoil, finite-size and QED corrections

- **D_{n,j}(ζ)** is the exact ratio of the point-nucleus Dirac binding energy to the Schrödinger energy for charge ζ,
  with j assigned by jj filling: l − ½ while k ≤ 2l, otherwise l + ½.
- **The bracket 1 + r_c(Zα)²(Z_eff/Z_a − 1)/n** is a Fermi–Segrè-type correction for the penetration of the
  valence electron into the region of full nuclear charge. Its five coefficients r_c (classes s, p½, p3/2, d, f) are
  fitted.
- **One-electron terms.** For one-electron ions we use the Dirac energy with a numerically solved finite-nuclear-size
  (FNS) shift, Barker–Glover recoil, the Uehling vacuum polarisation computed over the Dirac 1s density, and the
  one-loop self-energy function F_SE(Zα) from the all-order tabulation of Yerokhin and Shabaev [YerokhinShabaev2015].
  The 2 × 110 tabulated values are external theory inputs, not fitted parameters. In the many-electron formula, the
  QED and FNS shifts are applied only to 1s and 2s removal, scaled by (Z_eff/Z)².
- **No closed-form substitute for the table.** The closed-form low-order Zα expansion of F_SE matches the table at
  Z = 1 but diverges for Z ≳ 15. It gives 0.354 % MAPE on H-like ions, worse than no QED.

### 2.5 Density functional theory as a physics check

We also wrote a radial Kohn–Sham LSDA solver (Slater exchange + VWN5 correlation). It reproduces the NIST LDA
reference total energies [Kotochigova1997] to about 10⁻⁶ hartree (e.g. Ne −128.233481). Because that reference uses
the same functional, the agreement verifies the implementation only, not physical accuracy. ΔSCF ionization energies
from this solver have no fitted parameters and serve as an independent, non-formula check (§4.2). DFT does not enter
the formula.

### 2.6 Parameter count

The final model has 33 fitted global parameters:

| block | count |
|---|---|
| group screening τ_g | 5 |
| saturation denominator κ | 1 |
| class deviations δτ_c (ridge-shrunk; counted in full) | 19 |
| relativistic penetration r_c | 5 |
| Hund amplitudes x_l (p, d, f) | 3 |

There are no per-element or per-ion parameters, and the prediction code never reads a NIST ionization energy.
Fitted values are in Table A1 (Appendix) and `results/uni_final_params.json`.

---

## 3. Data and validation protocol

### 3.1 Data

The reference data are the 5847 successive ionization energies of the NIST Atomic Spectra Database [NISTASD] for
Z = 1–110, all charge states, with NIST ground configurations (`data/nist_ie.csv`). The database flags each value as
experimental (311 rows), semi-empirical (919) or theoretical (4617). Every metric below is computed over all 5847
rows unless stated otherwise. By status, the final model's all-data fit gives (evaluate.py, `status=` strata):
experimental 4.63 % MAPE (median 2.86 %, n = 311), semi-empirical 2.20 % (median 0.83 %, n = 919), theoretical
1.62 % (median 0.92 %, n = 4617). The experimental rows are mostly neutral atoms and low-charge ions (98 neutral, 193 of
311 with charge ≤ 2, median charge 2), the hardest regime
for the formula, so their higher error mirrors the neutral-atom weakness rather than a disagreement with experiment
specifically. **[TODO: insert the ASD version and access date.]**

Subsets used below:
- 108 neutral-atom first IEs;
- 110 hydrogen-like ions;
- 63 "rearranged" ions, whose ground configuration changes by more than one electron on ionization (e.g.
  V 3d³4s² → V⁺ 3d⁴).

The configuration of each ion is an input. Inside the table it is the NIST ground configuration; beyond Z = 110 the
Madelung order is used.

The error measure is MAPE = (100/M) Σ |IE_pred − IE_NIST| / IE_NIST, together with the median absolute percentage
error. All scores come from one shared scorer (`evaluate.py`). Fits minimise squared log-ratios ln(IE_pred/IE_NIST),
so a 4 eV and a 100 keV ionization energy carry equal relative weight.

### 3.2 Held-out splits

Every fitted model is refit on the training rows of each split with its own fitting routine:

| split | training rows | test / validation rows | n test rows | role |
|---|---|---|---|---|
| V1 | Z ≤ 36 | 37 ≤ Z ≤ 54 | 819 | selection |
| V2 | Z ≤ 44 | 45 ≤ Z ≤ 54 | 495 | selection |
| S1 | Z mod 5 ≠ 0 | Z mod 5 = 0 (unseen elements) | 1188 | selection |
| S3 | N mod 6 ≠ 0 | N mod 6 = 0 (unseen isoelectronic sequences) | 928 | selection |
| **S2** | Z ≤ 54 | **Z ≥ 55** | 4362 | **blind**, reported once |

**Selection score.** The score is the mean MAPE over V1, V2, S1 and S3. It was fixed before the final round of
model development and was the only quantity used to choose the final model. S2 (fit Z ≤ 54 → Z ≥ 55) was never
used for any choice: it was computed once per candidate, after its specification was frozen. The selection splits
themselves are not free of heavy elements, however. As pre-registered, S1 and S3 contain Z ≥ 55 rows in both training
and test sets (913 of the 1188 S1 test rows and 703 of the 928 S3 test rows, about 76 %), so heavy-atom
*interpolation* accuracy did inform the design choices; only heavy-atom *extrapolation* was kept blind. As a
robustness check (reporting only, after the audit), we recomputed the selection score with S1 and S3 restricted to
Z ≤ 54 in both training and test. The winner is unchanged: pa_hier_rel 1.844, then the 28-parameter variant without
the relativistic term 1.989, the 9-parameter bounded variant 2.094 and the best Push B candidate 2.118. This checks
the ranking of the frozen candidates, not the exploration path that produced them. The final model's numbers
were reproduced by an independent validation script, which matched the developing agent's own run with a difference
of 0.0 in every split.

### 3.3 History of the protocol (disclosure)

The protocol above was the third. In earlier stages of this project:

1. A purely fitted model reported S2 = 19 % using a bound on core screening chosen after seeing S2. Without that
   bound, the blind value is 170 %.
2. Later, an intermediate unified model was chosen on a score that included S2.

Both are reported here only as non-blind history (§4.3). Other disclosures:

- The developers of the final round knew which heavy near-neutral atoms the previous model failed on, and that an
  8-parameter model extrapolated well. This motivated the direction of the search, though not the choice, which was
  made on the selection score.
- About 57 variants were explored in that round (43 logged by one agent, about 14 by the other). The minimum
  selection score over them is therefore optimistic.
- The screening classes and the Hund term were designed earlier on all rows, including Z ≥ 55. This is structural
  leakage that we cannot remove retrospectively.

---

## 4. Results

### 4.1 Baselines versus the final model

**Table 1.** All values are MAPE in %. "All", "neutral" and "H-like" come from all-data fits. V1, V2, S1, S3 and S2
are held-out values after refitting. For parameter-free models nothing is fitted, so the split columns are simply the
error on those rows. Sources: `results/model_comparison.md`, `results/uni_validation.md`.

| model | fitted params | all (median) | neutral 1st IE | H-like | V1 | V2 | S1 | S3 | selection score | **blind S2** (median / neutral) |
|---|---|---|---|---|---|---|---|---|---|---|
| Bohr Ry Z²/n² | 0 | 1363 (172) | 2.3·10⁴ | 6.00 | 1050 | 1045 | 1317 | 1399 | 1203 | 1517 |
| Slater 1930, total-energy difference | 0 | 11.8 (7.48) | 51.8 | 6.00 | 13.4 | 14.4 | 11.6 | 12.4 | 12.9 | 11.6 (8.07 / 50.7) |
| Clementi–Raimondi, total-energy difference | 0 | 120 (18.7) | 2281 | 6.00 | 87.1 | 77.3 | 116 | 119 | 99.9 | 134 |
| Exact 1/Z series, completed square + rel/QED ("zexp") | 0 | 69.8 (7.20) | 1105 | 0.0012 | 53.4 | 54.0 | 68.3 | 74.3 | 62.5 | 77.8 (– / 1772) |
| GSHM final (purely fitted, blind) | 32 | 1.57 (0.74) | 7.61 | 0.19 | 34.1 | 1.65 | 1.59 | 1.70 | 9.77 | 170 (3.51 / 7393) |
| exact σ₁ + linear remainder ("u35") | 35 | 1.48 (0.73) | 6.43 | 0.0012 | 26.4 | 1.48 | 1.52 | 1.65 | 7.76 | 57.2 (2.00 / 3831) |
| pocket formula (hand-calculable) | 8 | 4.68 (3.04) | 16.8 | 1.16 | 4.4·10⁶ (fit diverges) | 2.55 | 4.62 | 5.24 | – | 11.3 (2.04 / 12.3) |
| lighter bounded variant ("pa_bound9") | 9 | 2.90 (1.39) | 12.0 | 0.0012 | 2.16 | 1.95 | 2.90 | 3.05 | 2.51 | 7.80 (– / 28.0) |
| **Screened Rydberg formula (final, "pa_hier_rel")** | **33** | **1.87 (0.94)** | **7.52** | **0.0012** | **3.09** | **1.60** | **1.91** | **2.00** | **2.15** | **6.60 (1.62 / 21.8)** |

The final model has the lowest selection score of all candidates (2.15). The next were the other developer's choice,
"pb_clip_pos" (30 parameters), at 2.37 and the 28-parameter hierarchical variant "pa_hier" at 2.37 (blind S2 13.9 %
and 8.47 %, respectively). In-distribution, the bound costs a little
accuracy: 1.87 % against 1.48 % for u35, and 7.52 % against 6.43 % on neutral atoms. In return, it removes the
catastrophic extrapolation failures.

### 4.2 Errors by stratum

**Table 2.** All-data fit, MAPE in % (source: `results/lit_comparison.md` / `docs/literature.md` §2.2).

| model | N ≤ 10 | 11 ≤ N ≤ 36 | N ≥ 37 | Z ≥ 55 | removed d | removed f |
|---|---|---|---|---|---|---|
| Slater 1930 | 4.89 | 9.38 | 16.8 | 11.6 | 16.2 | 12.9 |
| pocket formula (8 p) | 2.37 | 2.76 | 7.38 | 4.7 | 3.67 | 11.5 |
| GSHM final (32 p) | 0.698 | 1.16 | 2.31 | 1.48 | 0.982 | 2.67 |
| u35 (35 p) | 0.463 | 1.08 | 2.28 | 1.58 | 1.17 | 2.58 |
| **final (33 p)** | **0.439** | **1.43** | **2.87** | **1.92** | **1.52** | **3.51** |

- **Charge.** The error grows with the number of electrons and falls with ion charge. On ions of charge ≥ 3 the
  final model gives 1.61 %.
- **Where theory alone suffices.** For few-electron, highly charged ions, the parameter-free exact series is already
  good: median 0.7 % for N ≤ 10, and 0.48 % MAPE for N/Z ≤ 0.2 with two terms.
- **Neutral atoms, DFT check.** The parameter-free LSDA ΔSCF solver reaches 3.30 % on the 54 neutral first IEs with
  Z ≤ 54, and 0.80 % on all 171 ions with Z ≤ 18. It is better than every closed form on neutral atoms, but it is a
  numerical procedure and was run on only 207 rows.

Figures (existing): parity plot `results/figures/uni_parity.png`; first IEs against Z
`results/figures/uni_first_IE.png`; successive IEs of selected elements `results/figures/uni_successive.png`;
residuals `results/figures/uni_residuals.png`.

### 4.3 Blind extrapolation (S2)

With a fit on Z ≤ 54 only, the final model predicts the 4362 ionization energies of Z = 55–110 with 6.60 % MAPE and a
1.62 % median. The ancestors of the formula show what each ingredient contributes:

| model | what it adds | blind S2 MAPE (%) | S2 neutral (%) |
|---|---|---|---|
| GSHM final (fitted screening, no exact σ₁) | – | 170 | 7393 |
| u35 | exact σ₁ as the large-Z anchor, linear remainder | 57.2 | 3831 |
| final | + saturation bound Z_a ≤ Z_eff ≤ Z, hierarchical shrinkage | **6.60** | **21.8** |
| Slater 1930 (0 p, nothing extrapolated) | – | 11.6 | 50.7 |

- **What fixed extrapolation.** Exact σ₁ alone was not enough. The decisive step was confining the fitted part to a
  bounded remainder. In a 14-parameter ablation, the bound improved the selection score from 3.96 to 2.91; its
  blind S2 changed from 9.43 % to 5.53 %, and that S2 change was seen only after the choice.
- **Ions versus neutral atoms.** For ions (Z > N) the blind S2 MAPE is 6.4 %. For neutral atoms it is 21.8 %.
- **Per element.** Figure 1 shows the per-element median error on the S2 test rows.
- **Does the selection score predict extrapolation?** Figure 2 plots selection score against blind S2 for all fitted
  candidates. Candidates with a low selection score generally extrapolate well, but the relation is not monotone. A
  candidate with a worse selection score (pb_exp_pos, 3.73) has a lower blind S2 (6.03 %) than the final model. It
  was not promoted, because doing so would be selection on S2.

![Figure 1](../results/figures/paper_blind_s2_by_Z.png)

**Figure 1.** Blind extrapolation. Each point is the median absolute percentage error over all ion stages of one
element Z = 55–110. The three fitted models are fitted on Z ≤ 54 only; Slater's rules and the Kregar/Di Rocco SHM have
no fitted parameters. Generated from refits with the project's validation code; the S2 MAPEs reproduce Table 1
(final 6.603 %, u35 57.233 %, pocket 11.334 %, Slater 11.638 %).

![Figure 2](../results/figures/paper_selection_vs_blind.png)

**Figure 2.** Selection score (mean of V1, V2, S1, S3; no S2 data, although S1 and S3 contain Z ≥ 55 rows) against blind S2 MAPE for every fitted
candidate in `results/model_comparison.csv`. The pocket formula is omitted because its V1 fit diverges. The dashed
line is Slater's rules (0 parameters).

**Failures of the blind fit**, all reported here and none corrected, since correcting them now would be post hoc:
- **Heavy p-block neutrals.** Their first IEs are badly underestimated. Pb is predicted at 1.19 eV against 7.42 eV,
  Tl at 1.73 against 6.11, and Rn at 6.10 against 10.75.
- **Lr is negative.** Lr is predicted at −3.33 eV. The fitted r_c are negative, and the relativistic bracket is not
  bounded, so the Z_eff bound alone does not guarantee a positive IE once Zα is large.
- **f removal.** Third IEs of lanthanides and actinides (4f/5f removal) are about 2.0–2.5 times too high (mean
  2.2, median 2.1, 23 rows; worst Er²⁺ at 55.65 eV against 22.7 eV), because the Z ≤ 54 training set contains no
  f electrons.

### 4.4 Hydrogen-like ions: correction layers

**Table 3.** 110 H-like ions, Z = 1–110 (`docs/first_principles.md` §2.1).

| layer | MAPE (%) | median (%) | max (%) |
|---|---|---|---|
| Bohr Ry Z² | 6.0 | 4.2 | 19.5 |
| (a) Dirac, point nucleus | 0.190 | 0.118 | 0.91 |
| (b) + recoil + finite nuclear size | 0.113 | 0.109 | 0.248 |
| (c) + one-loop QED (Uehling computed; F_SE from [YerokhinShabaev2015]) | 0.00117 | 0.000154 | 0.0096 |
| (d) as (c) but with the closed-form Zα expansion of F_SE | 0.354 | 0.191 | 1.18 |

**Caveat.** The NIST H-like reference values are themselves computed from the same QED theory. Layer (c) therefore
shows consistency with that source, not independent validation. The remaining 10⁻⁵–10⁻⁴ residual at high Z is the
omitted two-loop QED, nuclear-polarisation and recoil-QED terms.

### 4.5 Head-to-head with a published screened hydrogenic model

Of the published SHMs, only the parameter-free Kregar/Di Rocco model could be fully specified from articles we were
able to read [Pomarico2005; DiRoccoLanzini2016]. We implemented it from its definitions:
- screening from hydrogenic densities with an exchange correction, iterated to self-consistency;
- E = −Σ q_i Z_i²/2n_i²;
- non-relativistic, Pauli and Dirac variants.

Two other cases:
- The constants of More [More1982] and Faussurier et al. [Faussurier1997] are in papers we could not access, and no
  verifiable reprint of the tables was found. We did not reconstruct them from memory.
- The constants of Mendoza et al. [Mendoza2011] were also not obtained.

**Fidelity of our implementation.**
- **Matches the published numbers:** same-shell Z → ∞ screening constants reproduce every printed digit (e.g. 1s
  0.3125, 2p 0.3492), and total energies agree with the published table within 0.71 %.
- **Differs:** cross-shell constants differ by 0.011 on average (at most 0.053). The original uses fitted closed
  forms whose coefficients we could not obtain [DiRocco1992]. Valence IEs of near-neutral ions come out 3–5 eV higher
  than the printed model values; Ar I, for example, is 18.96 eV here against 14.72 eV printed.

Our near-neutral numbers are therefore "the model as defined", not the authors' code.

**Table 4.** Same 5847 rows, same scorer, MAPE in % (`results/lit_comparison.md`). For the 0-parameter models the S2
column is the error on the Z ≥ 55 rows. For the fitted models it is the blind S2 refit.

| model | params | all (median) | neutral 1st IE | H-like | charge ≥ 3 | S2 (Z ≥ 55) |
|---|---|---|---|---|---|---|
| Kregar/Di Rocco SHM, non-relativistic | 0 | 13.6 (2.34) | 230 | 6 | – | 15.1 |
| Kregar/Di Rocco SHM + Pauli | 0 | 12.6 (1.79) | 225 | 1.16 | – | 13.9 |
| Kregar/Di Rocco SHM + Dirac | 0 | 12.5 (1.71) | 224 | 0.19 | 4.68 | 13.6 |
| Slater 1930 | 0 | 11.8 (7.48) | 51.8 | 6 | 10.2 | 11.6 |
| pocket formula | 8 | 4.68 (3.04) | 16.8 | 1.16 | – | 11.3 |
| **final** | **33** | **1.87 (0.94)** | **7.52** | **0.0012** | **1.61** | **6.6** |

The 224 % neutral figure is for our implementation, built from the published definitions. It is not fully faithful
for valence shells: it gives 18.96 eV for Ar I where the authors print 14.72 eV, and the authors' printed valence IEs
are 3–5 eV lower than ours, so this figure probably overstates the error of the authors' own model on neutral atoms.

The parameter-free SHM works well for ionised species, with a median of 1.54 % for charge ≥ 3. It fails near
neutrality: 37 of its 5847 predictions are zero or negative, all for ions of charge 0–4 with an open d or f shell.
The fair conclusion has three parts:

1. Where the Z → ∞ limit dominates, parameter-free and exact-σ₁ models both work.
2. Near neutrality, every parameter-free hydrogenic model degrades. The fitted remainder is what makes the formula
   useful there, at the price of 8–33 parameters.
3. The bounded fitted formula extrapolates better (6.6 %) than the parameter-free SHM scores on the same Z ≥ 55 rows
   (13.6 %). The unbounded fitted models extrapolate worse (57 % and 170 %).

A comparison with Dirac–Fock ionization energies from [Rodrigues2004] on the same rows has not yet been done.
**[TODO]**

### 4.6 The exact first-order screening constants compared with Slater's

**Table 5.** σ₁ for neutral-atom ground configurations, compared with Slater's σ. Excerpt; the full table for
N = 1–110 is in `results/fp_zexp_coefficients.csv`, with values for all 5847 rows in
`results/fp_zexp_rows_coefficients.csv`.

| N | configuration | removed | E₁ (exact) | ΔE₁ | σ₁ (ab initio) | σ (Slater) |
|---|---|---|---|---|---|---|
| 2 | 1s² | 1s | 5/8 | −0.625000 | 0.6250 | 0.30 |
| 3 | 1s²2s | 2s | 5965/5832 | −0.397805 | 1.5912 | 1.70 |
| 4 | 1s²2s² | 2s | 586373/373248 | −0.536469 | 2.1459 | 2.05 |
| 6 | 1s²2s²2p² (³P) | 2p | 10957679/3359232 | −0.931338 | 3.7254 | 2.75 |
| 8 | 1s²2s²2p⁴ | 2p | 4754911/839808 | −1.308370 | 5.2335 | 3.45 |
| 10 | 1s²2s²2p⁶ | 2p | 2455271/279936 | −1.636495 | 6.5460 | 4.15 |
| 11 | [Ne]3s | 3s | 9.635901 | −0.865071 | 7.7856 | 8.80 |
| 18 | [Ne]3s²3p⁶ | 3p | 17.980333 | −1.407535 | 12.6678 | 11.25 |
| 19 | [Ar]4s | 4s | 18.859971 | −0.879639 | 14.0742 | 16.80 |
| 29 | [Ar]3d¹⁰4s | 4s | 39.317547 | −1.347129 | 21.5541 | 25.30 |
| 36 | [Ar]3d¹⁰4s²4p⁶ | 4p | 50.145538 | −1.667284 | 26.6765 | 27.75 |

Notes:
- E₁ is the Hund-term single-configuration value. For N = 4 (Be) the Layzer-complex value 1.5592742 is used in σ₁.
  Decimal E₁ entries are exact rationals stored in the project cache.
- σ₁ = −n²ΔE₁.

Slater's rules describe the neutral atom; σ₁ is the exact screening of the bare-nucleus limit. The two disagree in a
systematic way (Figure 3):
- **Inner electrons:** the hydrogenic same-shell screening is stronger than Slater's 0.35 per electron (Ne 2p: 6.55
  against 4.15).
- **Valence electrons of heavy atoms:** σ₁ is much smaller than Slater's σ and than the screening implied by the
  measured IE (Na 3s: 7.79 against Slater 8.80 and an empirical 9.16; Cs 6s: 42.8 against Slater 52.8).

This difference is exactly what the fitted remainder D supplies near neutrality.

![Figure 3](../results/figures/fp_sigma_abinitio_vs_slater.png)

**Figure 3.** Exact first-order screening σ₁ against Slater's σ for neutral-atom configurations (existing project
figure).

### 4.7 Worked example and coverage

For the first IE of oxygen (1s²2s²2p⁴ → 2p³, NIST 13.618 eV):
1. σ₁ = 5.23348 and T = 3.3496.
2. h = 1.76652, so D = 0.71285 and Z_eff = 2.05367, which lies between Z_a = 1 and Z = 8.
3. The hydrogenic term is 14.3457 eV, or 14.3162 eV after relativity.
4. The Hund term is −0.8365 eV.
5. IE = 13.479 eV, an error of −1.02 %.

The analogous 8-parameter pocket calculation gives 9.70 eV (−28.8 %). For Mg²⁺ (NIST 80.144 eV), the formula gives
78.63 eV (−1.88 %). Full arithmetic is in `docs/unified.md` §6.

**Coverage test.** The public API was evaluated for every Z = 1–118 and every N = 1–Z: 7021 values. All are finite
and positive, with no crashes. There are 22 violations of the monotonicity IE(Z, N−1) > IE(Z, N):
- 6 at Pt–Bi, caused by 4f/5s ordering;
- 16 at Rf–Ds, where the input ground configuration jumps between neighbouring ions (e.g. Rf N = 68 is 4f¹²6s², N = 69
  is 4f¹⁴5d¹).

Above Z = 110 the predictions are qualitative. The 7p neutrals are too low: Og is predicted at 3.66 eV, implausibly
below its lighter congener Rn (NIST 10.75 eV), and the cause is again the negative r_c. **[VERIFY: if published
relativistic coupled-cluster values for Cn, Fl and Og are to be quoted, cite them from a checked source; the
project's approximate values are from memory and must not be printed as literature.]**

---

## 5. Discussion

### 5.1 What is new and what is rediscovered

**Rediscovered, and credited as such:**
- the 1/Z expansion and its exact first-order term (σ₁ is Layzer's Z → ∞ screening constant [Layzer1959], computed
  with textbook Slater-integral algebra [Cowan1981; FroeseFischer1997]);
- the screened hydrogenic form [Slater1930; Mayer1947; More1982], and self-consistent (Z, N)-dependent screening
  [Kregar1984; Pomarico2005];
- the 1/(Z_a + κ) remainder, which is Edlén-type isoelectronic behaviour [Edlen1964];
- the one-electron Dirac, recoil, finite-size and QED corrections [YerokhinShabaev2015; NISTASD];
- LSDA ΔSCF [KohnSham1965; Kotochigova1997].

**Plausibly new, and claimed with care:**
1. **A complete σ₁ table.** We found no published table of exact first-order screening constants covering every
   NIST ground configuration for N = 1–110. **[VERIFY against the Z-expansion literature, e.g. Safronova et al. and
   Layzer et al. 1964, before claiming this.]**
2. **A blind-tested hybrid formula.** A single closed form joins exact σ₁ to a bounded fitted remainder. It is
   validated on all 5847 NIST successive IEs, with a pre-registered selection score and a blind extrapolation from
   Z ≤ 54 to Z ≥ 55. In the SHM literature we read, we found no blind out-of-range extrapolation test.
3. **A published SHM scored on every NIST row**, re-implemented and scored on all 5847 rows (§4.5).
4. **A teaching rule.** The 8-parameter pocket formula is a hand-calculable successor to Slater's rules: 4.68 %
   against 11.8 % on all rows, and 11.3 % against 11.6 % on blind S2.

**Not claimed:**
- a "first formula for all ionization energies", since SHMs and Dirac–Fock tables already cover all ions;
- a "first-principles formula", since only σ₁ and the one-electron corrections are first principles and the final
  model has 33 fitted parameters;
- priority over machine-learning work. We found no ML study of successive IEs of all atomic ions, but a systematic
  search is still required. **[TODO]**

**The main lesson is methodological.** Fitting alone interpolates to about 1.6 % and extrapolates catastrophically
(170 %). Exact theory alone extrapolates smoothly but is about 70 % off overall. Fixing the Z² and Z terms by theory
and confining the fit to a bounded O(Z⁰) remainder gives about 2 % inside the data and 6.6 % in blind extrapolation.

### 5.2 Where the formula fails

1. **Heavy near-neutral atoms.** The blind neutral first-IE MAPE is 21.8 %, and the all-data neutral MAPE is 7.5 %
   (6.4 % for u35). Valence screening in a neutral atom is an all-order, non-perturbative effect. The parameter-free
   LSDA ΔSCF does better on neutral atoms (3.3 % for Z ≤ 54).
2. **An unbounded relativistic term.** The fitted r_c are all negative (−0.45 to −1.47), opposite to the direct
   relativistic contraction of s and p½ electrons. The bracket is not bounded, which produces the negative Lr value
   of the Z ≤ 54 fit and the low superheavy 7p IEs. A physically bounded relativistic term is the obvious next step.
   It must be chosen on V1/V2/S1/S3 only.
3. **f-electron removal.** This is the least accurate class: 3.5 % on all data, and about 2.0–2.5 times too high
   (mean 2.2) for third IEs in blind S2.
4. **Multiplet structure and rearranged configurations.** The formula has one Hund term per removed electron and no
   term-dependent multiplet energies. Ions whose ground configuration rearranges on ionization are described by
   single-configuration inputs. This is also the origin of the 22 monotonicity violations.
5. **Parameter economy.** 19 of the 33 parameters are shrunk class deviations. The 9-parameter bounded variant
   (selection 2.51, blind S2 7.80 %, all data 2.90 %) is a reasonable lighter alternative. It was not the
   pre-registered winner.
6. **Reference data.** Most NIST values are labelled theoretical or semi-empirical (4617 + 919 of 5847). Agreement
   with them is partly agreement with other theory.

### 5.3 Why no exact closed form exists for N ≥ 2

With V = Σ 1/r_ij, the many-electron Schrödinger equation does not separate. E(Z) is analytic in 1/Z only up to a
critical charge. For He, 1/Z_c ≈ 1.0975. Beyond first order, every coefficient E_k (k ≥ 2) is an infinite sum over the
hydrogenic continuum, with no known closed form even for two electrons; it is known only numerically
[ScherrKnight1963]. A closed formula for all ions must therefore combine exact low orders with an approximation for
the rest. The question is only which approximation, and how honestly it is validated. Our answer keeps the exact
orders exact and confines the fit to the remainder, under a bound that keeps it from extrapolating wildly.

---

## 6. Conclusions

- **The formula.** The screened Rydberg formula gives every successive ionization energy of every element from the
  configuration alone, with 33 global parameters.
- **Accuracy.** 1.87 % MAPE (median 0.94 %) on 5847 NIST values; 7.5 % for neutral atoms; 0.0012 % for H-like ions.
- **Blind extrapolation.** A fit on Z ≤ 54 predicts Z = 55–110 with 6.60 % MAPE, under the 11.6 % of Slater's rules and
  the 13.6 % of a published parameter-free screened hydrogenic model. This holds because the large-Z behaviour is
  fixed exactly by first-order perturbation theory and the fitted remainder is bounded.
- **Remaining defects.** Heavy neutral atoms, f-electron removal and the sign of the fitted relativistic coefficients.
- **Most reusable result.** The table of exact rational first-order screening constants σ₁ for N = 1–110.

---

## Data and code availability

All code, data and results are in the project repository **[TODO: archive with a DOI, e.g. Zenodo, and insert the
link]**:
- the NIST table (`data/nist_ie.csv`) and the shared scorer (`evaluate.py`);
- the exact first-order coefficients (`results/fp_zexp_coefficients.csv`, `results/fp_zexp_rows_coefficients.csv`);
- the final model and its parameters (`models/push_a/model.py`, `models/unified/final.py`,
  `results/uni_final_params.json`);
- predictions for all rows (`results/uni_predictions.csv`);
- the validation script (`models/unified/validate_blind.py`);
- the literature re-implementation (`models/literature/kregar_shm.py`).

A command-line and Python interface (`ionization.py`) evaluates the formula for any Z ≤ 118 and N ≤ Z.

## AI-assistance disclosure

The computations, code, analysis and a draft of this text were produced with AI agents (Anthropic Claude) working
under the author's direction. The author designed and directed the study and is responsible for its content.
**Before submission, the author must:**
- independently verify every number against the results files;
- verify every reference against the original publication (several bibliographic fields are marked unverified in
  `docs/references.bib`);
- word this statement according to the target journal's policy.

An AI system is not listed as an author.

## Acknowledgements

**[TODO]**

---

## References

Only entries from `docs/references.bib` are used. Fields marked there as STANDARD, CITED-IN or UNVERIFIED must be
checked against the publisher before submission. Most DOIs are omitted deliberately for that reason.

- [ClementiRaimondi1963] E. Clementi and D. L. Raimondi, Atomic screening constants from SCF functions, J. Chem. Phys. 38, 2686–2689 (1963).
- [ClementiRaimondiReinhardt1967] E. Clementi, D. L. Raimondi and W. P. Reinhardt, Atomic screening constants from SCF functions. II. Atoms with 37 to 86 electrons, J. Chem. Phys. 47, 1300–1307 (1967).
- [Chakravorty1993] S. J. Chakravorty, S. R. Gwaltney, E. R. Davidson, F. A. Parpia and C. Froese Fischer, Ground-state correlation energies for atomic ions with 3 to 18 electrons, Phys. Rev. A 47, 3649–3670 (1993).
- [Chung2005] H.-K. Chung, M. H. Chen, W. L. Morgan, Yu. Ralchenko and R. W. Lee, FLYCHK: generalized population kinetics and spectral model for rapid spectroscopic analysis for all elements, High Energy Density Phys. 1, 3–12 (2005).
- [Cowan1981] R. D. Cowan, *The Theory of Atomic Structure and Spectra* (University of California Press, Berkeley, 1981).
- [Crilly2023] A. J. Crilly et al., SpK: a fast atomic and microphysics code for the high-energy-density regime, High Energy Density Phys. (2023), doi:10.1016/j.hedp.2023.101053, arXiv:2211.16464.
- [DalgarnoStewart1958] A. Dalgarno and A. L. Stewart, A perturbation calculation of properties of the helium iso-electronic sequence, Proc. R. Soc. Lond. A 247, 245–259 (1958). [VERIFY]
- [DiRocco1992] H. O. Di Rocco, Braz. J. Phys. 22, 227 (1992). [VERIFY title]
- [DiRoccoLanzini2016] H. O. Di Rocco and F. Lanzini, Breit and quantum electrodynamics energy contributions in multielectron atoms from the relativistic screened hydrogenic model, Braz. J. Phys. 46, 175–183 (2016), doi:10.1007/s13538-015-0397-9.
- [Edlen1964] B. Edlén, Atomic spectra, in *Handbuch der Physik* vol. 27, ed. S. Flügge (Springer, Berlin, 1964), pp. 80–220. [VERIFY pages]
- [Faussurier1997] G. Faussurier, C. Blancard and A. Decoster, New screening coefficients for the hydrogenic ion model including l-splitting for fast calculations of atomic structure in plasmas, J. Quant. Spectrosc. Radiat. Transfer 58, 233 (1997). [VERIFY]
- [Faussurier2008] G. Faussurier, C. Blancard and P. Renaudin, Equation of state of dense plasmas using a screened-hydrogenic model with l-splitting, High Energy Density Phys. 4, 114–123 (2008).
- [FroeseFischer1997] C. Froese Fischer, T. Brage and P. Jönsson, *Computational Atomic Structure: An MCHF Approach* (Institute of Physics Publishing, Bristol, 1997).
- [KohnSham1965] W. Kohn and L. J. Sham, Self-consistent equations including exchange and correlation effects, Phys. Rev. 140, A1133–A1138 (1965).
- [Koopmans1934] T. Koopmans, Über die Zuordnung von Wellenfunktionen und Eigenwerten zu den einzelnen Elektronen eines Atoms, Physica 1, 104–113 (1934).
- [Kotochigova1997] S. Kotochigova, Z. H. Levine, E. L. Shirley, M. D. Stiles and C. W. Clark, Local-density-functional calculations of the energy of atoms, Phys. Rev. A 55, 191–199 (1997).
- [Kregar1984] M. Kregar, Phys. Scr. 29, 438 (1984). [VERIFY title]
- [Kregar1985] M. Kregar, Phys. Scr. 31, 246 (1985). [VERIFY title]
- [LanziniDiRocco2015] F. Lanzini and H. O. Di Rocco, Screening parameters for the relativistic hydrogenic model, High Energy Density Phys. 17, 240–247 (2015).
- [Layzer1959] D. Layzer, On a screening theory of atomic spectra, Ann. Phys. (N.Y.) 8, 271 (1959). [VERIFY pages]
- [LayzerEtAl1964] D. Layzer et al., Ann. Phys. (N.Y.) 29, 101 (1964). [VERIFY title and authors]
- [Layzer1967] D. Layzer, Int. J. Quantum Chem. 1S, 45 (1967). [VERIFY title]
- [Martel1998] P. Martel, J. G. Rubiano, J. M. Gil, L. Doreste and E. Mínguez, J. Quant. Spectrosc. Radiat. Transfer 60, 623 (1998). [VERIFY title]
- [Mayer1947] H. Mayer, Methods of opacity calculations, Los Alamos Scientific Laboratory report LA-647 (1947).
- [Mendoza2011] M. A. Mendoza, J. G. Rubiano, J. M. Gil, R. Rodríguez, R. Florido, P. Martel and E. Mínguez, A new set of relativistic screening constants for the screened hydrogenic model, High Energy Density Phys. 7, 169–179 (2011). [VERIFY volume/pages]
- [More1982] R. M. More, Electronic energy levels in dense plasmas, J. Quant. Spectrosc. Radiat. Transfer 27, 345–357 (1982).
- [NISTASD] A. Kramida, Yu. Ralchenko, J. Reader and NIST ASD Team, NIST Atomic Spectra Database (ver. 5.x), National Institute of Standards and Technology, Gaithersburg, MD, https://physics.nist.gov/asd. [Insert exact version and access date]
- [Pomarico2005] J. Pomarico, D. I. Iriarte and H. O. Di Rocco, An efficient screening approach to be used in plasma modeling and ion-surface collision experiments, Braz. J. Phys. 35(1), 130–135 (2005).
- [Rodrigues2004] G. C. Rodrigues, P. Indelicato, J. P. Santos, P. Patté and F. Parente, Systematic calculation of total atomic energies of ground state configurations, At. Data Nucl. Data Tables 86, 117–233 (2004), doi:10.1016/j.adt.2003.11.005.
- [Rozsnyai1972] B. F. Rozsnyai, Relativistic Hartree-Fock-Slater calculations for arbitrary temperature and matter density, Phys. Rev. A 5, 1137–1149 (1972). [VERIFY pages]
- [Rubiano2002] J. G. Rubiano, R. Rodríguez, J. M. Gil, F. H. Ruano, P. Martel and E. Mínguez, A screened hydrogenic model using analytical potentials, J. Quant. Spectrosc. Radiat. Transfer 72, 575 (2002).
- [Safronova1993] U. I. Safronova et al., Phys. Scr. 47, 364 (1993). [VERIFY title and authors]
- [ScherrKnight1963] C. W. Scherr and R. E. Knight, Two-electron atoms III. A sixth-order perturbation study of the 1¹S ground state, Rev. Mod. Phys. 35, 436–442 (1963). [VERIFY pages]
- [Slater1930] J. C. Slater, Atomic shielding constants, Phys. Rev. 36, 57–64 (1930).
- [YerokhinShabaev2015] V. A. Yerokhin and V. M. Shabaev, Lamb shift of n = 1 and n = 2 states of hydrogen-like atoms, 1 ≤ Z ≤ 110, J. Phys. Chem. Ref. Data 44, 033103 (2015).

---

## Appendix A. Fitted parameters (all-data fit)

**Table A1.** Source: `results/uni_final_params.json`, which is identical to `results/pa_params.json`.

| block | values |
|---|---|
| τ_g (same, in, core, df, out) | 0.3928, 0.7879, 2.2243, 2.5488, 10.5688 |
| κ | 1.7527 (1.8027 as used in the code, which applies \|κ\| + 0.05) |
| δτ_c (19 classes) | same_s −0.1552, same_p 0.0098, same_d 0.0004, same_f 0.1457, sn_in_p −0.1133, n1_sp_sp −0.3916, n1_sp_d 0.3009, n1_sp_f 0.2042, n2_sp_sp −0.2195, n2_sp_df −0.0622, d_near −0.6962, n1_d_d −0.2689, n1_d_f 0.3286, f_near 0.5476, n1_f_f −0.0918, out_d −0.0176, out_f 0.0173, core_sp 0.2828, core_df 0.1812 |
| r_c (s, p½, p3/2, d, f) | −0.4543, −0.8129, −1.1533, −0.8426, −1.4672 |
| x_l (p, d, f) | 0.2049, 0.4402, 0.4612 |

**Pocket formula (8 parameters, all-data fit, `results/uni_params.json`).**

$$Z_\mathrm{eff} = Z - \sum_g s_g\nu_g - t\,\frac{N-1}{Z-N+1+\kappa},\qquad
\mathrm{IE}=\mathrm{Ry}\,\frac{Z_\mathrm{eff}^2}{n^2}\Big[1+\frac{(Z_\mathrm{eff}\alpha)^2}{n^2}\Big(\frac{n}{j+\tfrac12}-\frac34\Big)\Big]+\mathrm{Ry}\,\frac{x\,K_l(k)}{n^2}$$

The fitted values are s_same 0.7499, s_in 0.7514, s_core 0.8432, s_df 0.8408, s_out 1.0497, t 1.3725, κ 9.815
(9.865 as used) and x 0.508.
