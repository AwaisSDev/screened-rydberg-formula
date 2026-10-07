# A screened Rydberg formula for the successive ionization energies of all atoms and ions

**Muhammad Awais**

Independent researcher. Correspondence: mawais9171@gmail.com



---

**Abstract—**We present a closed-form expression, the *screened Rydberg formula*, for the successive ionization
energy IE(Z, N) of any atom or ion. Its inputs are the nuclear charge Z and the ground electron configuration. The
formula keeps the hydrogenic form Ry Z_eff²/n² and splits the effective charge into two parts. The first is the
first-order screening constant σ₁ = −n²ΔE₁ of the 1/Z perturbation expansion. It has no adjustable parameter, it is a
rational number for every configuration, and we tabulate it for N = 1–110. The second is a fitted higher-order
remainder that depends on the ion charge in an Edlén-type form. A saturation construction keeps
Z − N + 1 ≤ Z_eff ≤ Z for any parameter values (given 0 ≤ σ₁ ≤ N − 1, which holds for every configuration we use)
and preserves the exact Z² and Z coefficients of the 1/Z series. One-electron Dirac, recoil, finite-nuclear-size and
QED corrections enter without fitted parameters. The final model has 33 global parameters and none per element or
ion. On the 5847 successive ionization energies in the NIST Atomic Spectra Database (Z = 1–110) its mean absolute
percentage error (MAPE) is 1.87 % (median 0.94 %). It is 7.52 % for the first ionization energies of neutral atoms and
0.0012 % for hydrogen-like ions; the H-like figure agrees with, but is not independent of, the QED theory behind the
NIST values. We chose the model with a pre-registered score that never used the extrapolation test, although its
interpolation splits contain heavy elements. Fitted on Z ≤ 54 and used to predict Z = 55–110, the model gives 6.60 %
MAPE (median 1.62 %). On the same rows Slater's rules give 11.6 %, a published parameter-free screened-hydrogenic
model 13.6 %, and two earlier fitted variants without the bound 57.2 % and 170 % (Table 1). Heavy neutral atoms
remain the weak point, with 21.8 % MAPE on their first ionization energies in that test. The Z ≤ 54 fit puts Pb and
Tl at about 1.2 and 1.7 eV (NIST 7.42 and 6.11 eV) and makes Lr negative, and the all-data fit puts oganesson
(3.66 eV) far below radon (NIST 10.75 eV). We report these failures and have not corrected them.

**Index Terms—** ionization energy; screening constants; 1/Z expansion; screened hydrogenic model; isoelectronic
sequences; Slater's rules; blind validation.

---

## 1. Introduction

The Bohr energy Ry Z²/n² is exact only for a non-relativistic one-electron ion. For every other atom and ion the
ionization energy (IE) depends on electron–electron repulsion, which has no closed-form solution for N ≥ 2. Four
lines of work bear on the problem.

Slater [1] replaced Z by an effective charge Z − σ, with σ given by counting rules fitted by hand to atomic data.
Clementi and Raimondi [2], [3] obtained σ from optimised self-consistent-field orbital exponents of neutral atoms.
Both were designed for orbitals of neutral atoms, not for the full set of successive IEs. Scored as IE formulas on
all NIST ions with a total-energy difference, Slater's rules give 11.8 % MAPE and Clementi–Raimondi 120 % (§4).

Layzer [4], [5] showed that the non-relativistic energy of a fixed configuration is an asymptotic series
E = Z²E₀ + ZE₁ + E₂ + …. E₀ is hydrogenic, and E₁ is a rational combination of hydrogenic Slater integrals, which
fixes the Z → ∞ limit of the screening constant exactly. Higher orders need continuum sums [6], [7]. Along
isoelectronic sequences, Edlén [8] described the smooth variation of screening with ion charge by expansions in
1/(ζ + s). Z-expansion codes with relativistic corrections were developed for selected sequences [9].

Screened hydrogenic models (SHMs), starting with Mayer [10] and More [11], are used in plasma-physics codes. Their
screening constants have been tabulated or fitted with l-splitting [12], [13], derived from analytical potentials
[14], [15], fitted by a genetic algorithm with relativistic subshells [16], or computed self-consistently as functions
of Z and N without parameters [17]–[22]. Average-atom and kinetics codes use SHMs because they are fast and cover all
ions [23]–[25].

Ab initio methods give more accurate numbers: Koopmans' theorem [26], ΔSCF Hartree–Fock and Kohn–Sham DFT [27], [28],
correlated non-relativistic energies [29], and Dirac–Fock total energies for all ground configurations up to Z = 118
[30]. They are numerical procedures, not formulas.

We found no closed-form expression that covers every ion of every element from the configuration alone, reproduces
the large-Z behaviour of the 1/Z expansion exactly, reports every fitted parameter, and is validated on all NIST
successive IEs with held-out tests fixed in advance, including an extrapolation to heavier elements. The SHM papers
we read report accuracy on their fitting data or, like [20], on selected sequences. Scored on our rows (§4.5), the
constants of Mendoza et al. give 2.82 % MAPE on the 5011 rows they cover. This paper addresses that gap. It does not
propose a new law of atomic physics, and §5.1 lists what is rediscovered.

The paper makes four contributions. It tabulates σ₁ for every NIST ground configuration with N = 1–110, as a
parameter-free counterpart to Slater's σ that is exact as Z → ∞ (§2.1, §4.6). It defines the screened Rydberg
formula, σ₁ plus a bounded Edlén-type remainder with 33 global parameters (§2.2–2.4). It applies one validation
protocol, with a pre-registered selection score and a single extrapolation test, to the final model, its ancestors
and the published baselines (§3, §4). And it compares the formula with two published SHMs, re-implemented and scored
with the same scorer (§4.5): the parameter-free Kregar/Di Rocco model on all 5847 rows, and the fitted constants of
Mendoza et al. (2011) on the 5011 rows they cover.

---

## 2. Theory

Section 2.1 uses atomic units; energies are converted to eV with Ry = 13.6057 eV. Z is the nuclear charge, N the
number of electrons before ionization, and Z_a = Z − N + 1 the charge seen by the departing electron. (n, l) is the
removed subshell and k its occupancy.

### 2.1 Exact first-order screening from the 1/Z expansion

Scaling r → ρ/Z turns the non-relativistic Hamiltonian into H = Z²[H₀ + Z⁻¹V], with H₀ hydrogenic and
V = Σ 1/ρ_ij. Rayleigh–Schrödinger perturbation theory in 1/Z then gives, for a fixed configuration and term,
E(Z, N) = Z²E₀ + ZE₁ + E₂ + …, and the ionization energy is

{{L98}}

At zeroth order ΔE₀ = 1/(2n²), which is Bohr's formula. At first order E₁ = ⟨V⟩ over the Z = 1 hydrogenic state.
Every hydrogenic radial integral F^k and G^k is a rational number, because each product P_a P_c is a polynomial times
e^{−βr} with rational β; for example F⁰(1s,1s) = 5/8, F⁰(1s,2s) = 17/81 and G⁰(1s,2s) = 16/729. We evaluate the
integrals in exact rational arithmetic and combine them with exact angular algebra. For open subshells E₁ is
evaluated for the Hund ground term by projecting onto the highest-weight (S, L) subspace. Where hydrogenic degeneracy
couples configurations (the Layzer complex, e.g. 1s²2s² with 1s²2p²), V is diagonalised within the complex. For two
or more open subshells we use a high-spin coupled term, which is the one approximation in σ₁.

Completing the square in the first two terms defines the first-order screening constant

{{L112}}

σ₁ is Layzer's Z → ∞ screening constant. It has no adjustable parameter and depends on the configuration, not on Z.
For He-like ions E₁(1s²) = 5/8, so σ₁(He) = 0.6250 (Slater: 0.30). For Li-like ions
E₁(1s²2s) = 5965/5832 = 1.0228052 = 5/8 + 2·17/81 − 16/729, so ΔE₁ = −0.397805 and σ₁ = 1.5912. For Be-like ions the
complex value E₁ = 1.5592742 agrees with Layzer's.

The series fails for near-neutral atoms (§4.6), because its expansion parameter is effectively N/Z. For a neutral
atom Z − σ is only 1–3, so a 10 % error in σ becomes a 100–1000 % error in IE. The higher orders have no closed form:
E₂ already requires sums over the hydrogenic continuum. The completed-square guess ΔE₂ = ΔE₁²/(4ΔE₀), for example, is
24–26 % too large for He- and Li-like ions. For this reason the remainder in §2.2 is fitted.

### 2.2 The screened Rydberg form with a bounded remainder

The final model (code name pa_hier_rel) is

{{L130-132}}

where the last term applies only to the removal of an ns electron with n ≤ 2, and

{{L136-141}}

Each of the other N − 1 electrons, in a subshell (n′, l′), is counted in exactly one of five screening groups ν_g.
The group *same* holds the k − 1 other electrons of (n, l). The group *out* holds outer electrons, with n′ > n or with
n′ = n and l′ > l. All remaining electrons are inner. For an s or p target, inner electrons with n′ ≥ n − 1 form *in*
(the ns electrons of an np target plus the whole (n − 1) shell, including its d and f electrons) and those with
n′ ≤ n − 2 form *core*. For a d or f target, every inner electron, from the same-n lower-l subshells down to 1s, is
in *df*. Hence Σ_g ν_g = N − 1, and each group has a coefficient τ_g. As (ν_same, ν_in, ν_core, ν_df, ν_out), O 2p⁴
gives (3, 4, 0, 0, 0) and Na 3s gives (0, 8, 2, 0, 0). Fe 3d⁶4s² with 4s removed gives (1, 14, 10, 0, 0), because 3d⁶
belongs to the (n − 1) shell. Pb 6p² gives (1, 20, 60, 0, 0), with 6s²5s²5p⁶5d¹⁰ in *in* and every electron with
n′ ≤ 4, including 4f¹⁴, in *core*.

The screening classes ν_c (Table A2) refine the groups. An electron's class depends on the target type (s/p, d or
f), on n − n′ and on l′. Each class lies inside one group, and an electron of class c in group g contributes
τ_g + δτ_c to T. A ridge penalty (10⁻⁴ per row) shrinks the deviations toward their group value. Two of the 21
classes are empty in all 5847 rows and carry no δτ_c: sn_out (n′ = n, l′ > l) and out_sp (n′ > n, s/p target). That
leaves the 19 of Table A1. The bounded 9-parameter variant (pa_bound9) and the pocket formula (Appendix A) use the
five groups without classes. An independent implementation of the rule reproduces the code's group and class counts
with 0 mismatches on all 5847 rows and on all 7021 configurations with Z ≤ 118 (`tools/check_grouping_rule.py`).

The removed subshell (n, l) is the one whose occupancy drops from the N-electron configuration to the NIST ground
configuration of the (N − 1)-electron ion. The N-electron configuration is the NIST ground configuration, or the
Madelung order if the ion is not tabulated. For a rearranging ion (e.g. V 3d³4s² → V⁺ 3d⁴) it is the subshell that
loses most electrons, with ties going to larger n and then larger l. If the (N − 1)-electron ion is not tabulated,
or the user supplies the configuration, the outermost subshell (largest n, then largest l) is removed.

The Hund term uses K_l(k) = [P(k) − P(k−1)] − 2l(k−1)/(4l+1), the change in the number of parallel-spin pairs
relative to a statistical average. P(k) is the number of parallel-spin pairs of l^k under Hund's first rule.
Equivalently K_l(k) = (2l+1)(k−1)/(4l+1) for k ≤ 2l+1 and K_l(k) = −K_l(4l+3−k) above half filling, so K_s ≡ 0,
K_p = 0, 0.6, 1.2, −1.2, −0.6, 0 for k = 1…6, K_d = (5/9)(0, 1, 2, 3, 4, −4, −3, −2, −1, 0) and
K_f = (7/13)(0, 1, …, 6, −6, …, −1, 0). The amplitude x_l is fitted for l = p, d, f; s electrons have no Hund term.

Three inputs need a precise definition. σ₁ is evaluated for the *frozen* configuration, the N-electron configuration
minus one (n, l) electron, and not for the ground configuration of the ion. The two differ only for the 63
rearranged rows, and the choice keeps ΔE₀ = 1/(2n²) on those rows. The reduced-mass factor is μ(Z) = M/(M + m_e),
where M is the nuclear mass of the isotope used in the QED tabulation [31]. The QED and finite-nuclear-size shift is
applied only to 1s and 2s removal (427 rows); it is the one-electron shift of §2.4 for charge Z, scaled by
(Z_eff/Z)², and reaches at most 0.9 % of the IE (at Z = 110, N = 2). For N = 1, T = 0 and so D = 0.

A reader cannot reconstruct the frozen-configuration σ₁, μ(Z) or the QED/FNS shift from the text alone, so we publish
every per-row input in `results/model_inputs.csv`: removed subshell, k, j, σ₁, ν_g, ν_c, K, μ and QED/FNS. A short
script that uses only this table, the equations above and the parameters of Table A1 (`tools/verify_from_inputs.py`,
which imports no model code) reproduces the production code on all 5847 rows to 7·10⁻¹⁶ relative, for the final
model and for pa_bound9.

The saturation form of D has two properties for any parameter values. First, it bounds the effective charge,
Z_a ≤ Z_eff ≤ Z. Writing D = sign(T)·h·u/(1+u) with u = |T|/[h(Z_a+κ)] gives |D| < h, so Z_eff can fall neither
below the fully screened, non-penetrating limit Z_a nor above the bare charge Z. Because h = (N−1) − σ₁ or σ₁, this
requires 0 ≤ σ₁ ≤ N − 1. We have checked rather than proved that condition: it holds for all 7021 ground and Madelung
configurations with Z ≤ 118, and the code warns if a user-supplied configuration breaks it. Second, it keeps the
exact asymptotics. At fixed N and Z_a → ∞, D → T/(Z_a+κ) = O(1/Z), so the expansion of IE begins
Ry[Z² − 2Zσ₁ + O(1)]/n² and its Z² and Z coefficients are those of the exact 1/Z series. The fitted part represents
only ΔE₂ and higher orders, that is relaxation, correlation and penetration at low ion charge.

### 2.3 Charge-dependent penetration as an Edlén-type term

We found the remainder's dependence on 1/(Z_a + κ) empirically (docs/semi_empirical.md). Along every isoelectronic
sequence the excess charge p = Z_eff − Z_a grows with ion charge q as p∞ − τ/(Z_a + κ), with one κ for all
sequences. For the Na sequence (3s), p = 0.84, 1.15, 1.34, 1.46, 1.56 for q = 0–4 and 2.11 at q = 20. Adding the term
reduced the error of the purely fitted model from about 9 % to about 2 % MAPE. This is an independent rediscovery of
the isoelectronic regularity that Edlén formalised [8], and it is consistent with Layzer's theory, in which the
screening constant is σ₀ + σ₁′/Z + …. What is specific here is narrower: one denominator for all sequences, combined
with exact σ₁ and with the saturation bound of §2.2.

### 2.4 Relativistic, recoil, finite-size and QED corrections

F_{n,j}(ζ) is the exact ratio of the point-nucleus Dirac binding energy to the Schrödinger energy for charge ζ. With
x = ζα and κ = j + ½, F = (2n²/x²)·{1 − [1 + (x/(n − κ + √(κ² − x²)))²]^(−1/2)}. j follows jj filling, l − ½ while
k ≤ 2l and l + ½ otherwise. We write F to avoid confusion with the screening remainder D.

The bracket 1 + r_c(Zα)²(Z_eff/Z_a − 1)/n is a Fermi–Segrè-type correction for the penetration of the valence
electron into the region of full nuclear charge. Its five coefficients r_c, for the classes s, p½, p3/2, d and f,
are fitted.

For one-electron ions we use the Dirac energy with a numerically solved finite-nuclear-size (FNS) shift,
Barker–Glover recoil, the Uehling vacuum polarisation computed over the Dirac 1s density, and the one-loop self-energy
function F_SE(Zα) from the all-order tabulation of Yerokhin and Shabaev [31]. The 2 × 110 tabulated values are
external theory inputs, not fitted parameters. In the many-electron formula the QED and FNS shifts are applied only to
1s and 2s removal, scaled by (Z_eff/Z)². The closed-form low-order Zα expansion of F_SE cannot replace the table: it
matches the table at Z = 1 but diverges for Z ≳ 15, and it gives 0.354 % MAPE on H-like ions, worse than no QED.

### 2.5 Density functional theory as a physics check

We also wrote a radial Kohn–Sham LSDA solver (Slater exchange plus VWN5 correlation). It reproduces the NIST LDA
reference total energies [28] to about 10⁻⁶ hartree (e.g. Ne −128.233481). That reference uses the same functional,
so the agreement verifies the implementation, not the physics. ΔSCF ionization energies from this solver have no
fitted parameters and serve as an independent check that is not a formula (§4.2). DFT does not enter the formula.

### 2.6 Parameter count

The final model has 33 fitted global parameters:

{{L249-255}}

There are no per-element or per-ion parameters, and the prediction code never reads a NIST ionization energy. The
fitted values are in Table A1 (Appendix) and `results/uni_final_params.json`; the class definitions are in Table A2.

Each model name in this paper refers to one model:

{{L263-267}}

pa_bound9 uses the pocket formula's five electron groups but is a different model.

---

## 3. Data and validation protocol

### 3.1 Data

The reference data are the 5847 successive ionization energies of the NIST Atomic Spectra Database [32] for
Z = 1–110, all charge states, with NIST ground configurations (`data/nist_ie.csv`). We took them from ASD version
5.12 [32], the current version at the time, on or before 5 October 2026; the project log records results computed
from these data on that date, and the original download timestamp was not kept. The database flags 311 values as
experimental, 919 as semi-empirical and 4617 as theoretical. Every metric below uses all 5847 rows unless stated
otherwise. By status, the final model's all-data fit gives 4.63 % MAPE on experimental rows (median 2.86 %,
n = 311), 2.20 % on semi-empirical rows (median 0.83 %, n = 919) and 1.62 % on theoretical rows (median 0.92 %,
n = 4617) (evaluate.py, `status=` strata). The experimental rows are mostly neutral atoms and low-charge ions (98
neutral, 193 of 311 with charge ≤ 2, median charge 2), the hardest regime for the formula, so their higher error
reflects the neutral-atom weakness rather than a disagreement with experiment.

Three subsets recur below: the 108 neutral-atom first IEs, the 110 hydrogen-like ions, and 63 rearranged ions whose
ground configuration changes by more than one electron on ionization (e.g. V 3d³4s² → V⁺ 3d⁴).

The configuration of each ion is an input. Inside the table it is the NIST ground configuration. For ions outside
it we use the Madelung order: all ions with Z > 110, and 258 ions with Z = 104–110 that the table lacks. These ions
enter no fit and no score, only the 7021-configuration checks (§4.7) and predictions outside the table.

The error measure is MAPE = (100/M) Σ |IE_pred − IE_NIST| / IE_NIST, together with the median absolute percentage
error. All scores come from one shared scorer (`evaluate.py`). Fits minimise squared log-ratios ln(IE_pred/IE_NIST),
so a 4 eV and a 100 keV ionization energy carry equal relative weight.

### 3.2 Held-out splits

Every fitted model is refit on the training rows of each split with its own fitting routine:

{{L305-311}}

The selection score is the mean MAPE over V1, V2, S1 and S3. It was fixed before the final round of model
development and was the only quantity used to choose the final model. S2 (fit Z ≤ 54, predict Z ≥ 55) was never used
for a choice: we computed it once per candidate, after the candidate's specification was frozen. The selection splits
do contain heavy elements. As pre-registered, S1 and S3 have Z ≥ 55 rows in both training and test sets (913 of the
1188 S1 test rows and 703 of the 928 S3 test rows, about 76 %), so heavy-atom *interpolation* accuracy informed the
design choices and only heavy-atom *extrapolation* stayed blind. As a robustness check, made after the audit and for
reporting only, we recomputed the selection score with S1 and S3 restricted to Z ≤ 54 in training and test. The
ranking is unchanged: pa_hier_rel 1.844, the 28-parameter variant without the relativistic term 1.989, the bounded
9-parameter variant 2.094 and the lowest-scoring Push B candidate 2.118. This checks the ranking of the frozen
candidates, not the exploration that produced them. An independent validation script reproduced the final model's
numbers, matching the developing agent's own run with a difference of 0.0 in every split.

### 3.3 History of the protocol

The protocol above is the project's third. In an earlier stage a purely fitted model reported S2 = 19 % using a bound
on core screening chosen after seeing S2; without that bound its S2 value is 170 %. Later, an intermediate unified
model was chosen on a score that included S2. We report both only as non-blind history (§4.3).

Several other facts bear on how blind the final choice was. The developers of the final round knew which heavy
near-neutral atoms the previous model failed on, and that an 8-parameter model extrapolated well; this shaped the
direction of the search, though the choice itself was made on the selection score. About 57 variants were explored in
that round (43 logged by one agent and about 14 by the other), so the minimum selection score over them is
optimistic. The screening classes and the Hund term were designed earlier on all rows, including Z ≥ 55, which is
structural leakage that cannot be removed now.

After the first draft of this paper we tried four more pre-registered variants (`docs/preregistration_v2.md`,
`models/v2/NOTES.md`): a relativistic factor that cannot change sign, a second-order charge term, and a κ per orbital
type. The sign-preserving relativistic factor fixes the signs of Lr and Og, but no variant improved the selection
score (2.15 → 2.18–2.49), so the model was left unchanged. That brings the number of explored variants to about 61.
The pre-registration file was fingerprinted (sha256, `models/v2/PREREG_HASH.txt`) before that round and committed to
git afterwards.

### 3.4 Computational reproducibility

The frozen run used Python 3.13 with NumPy and SciPy; their exact versions were not recorded. We repeated every refit
in a second software environment (Python 3.11.9, NumPy 2.4.4, SciPy 1.17.1; `docs/review/reproducibility_py311.md`).
The final model reproduces to within 0.005 percentage points in every split: selection score 2.1473 vs 2.1481, blind
S2 6.6030 vs 6.6031 %, and 1.8742 % on all data in both. The unbounded reference models depend on the least-squares
path. The pocket formula's V1 fit, which diverges in the frozen run (MAPE 4.4·10⁶ %), converges in the second
environment to 4.35 % (selection score 4.19). The u29 V1 fit converges in the frozen run and diverges in the second,
and u35's blind S2 moves from 57.2 % to 56.5 %. Their V1 entries in Table 1 therefore describe the optimizer as much
as the model, and we do not use the pocket formula's V1 divergence as evidence against it.

---

## 4. Results

### 4.1 Baselines versus the final model

Every cross-model number in this section comes from one generated matrix with the row count in each cell,
`results/benchmark_matrix.md` (`tools/benchmark_matrix.py`), and comparisons between models always use the same rows.

**Table 1.** MAPE in %. "All", "neutral" and "H-like" are from all-data fits; V1, V2, S1, S3 and S2 are held-out values
after refitting. Parameter-free models are not fitted, so their split columns are the error on those rows. Sources:
`results/model_comparison.md`, `results/uni_validation.md`.

{{L379-389}}

† Optimizer-path dependent (§3.4): the V1 fit diverges in the frozen run and converges to 4.35 % in a second
software environment.

Slater's rules have an all-row MAPE of 11.8 %, 6.3 times that of the final model (1.87 %) and 4.1 times that of
pa_bound9 (2.90 %). Their median is 7.48 % against 0.94 %, and their blind S2 error is 11.6 % against 6.60 %, a factor
of 1.8 (ratios from `results/benchmark_matrix.md`).

The final model has the lowest selection score of all candidates (2.15). The next were the other developer's choice,
pb_clip_pos (30 parameters), at 2.37, and the 28-parameter hierarchical variant pa_hier at 2.37, with blind S2 errors
of 13.9 % and 8.47 %. The bound costs some accuracy inside the data, 1.87 % against 1.48 % for u35 and 7.52 % against
6.43 % on neutral atoms, and in return it removes the catastrophic extrapolation failures.

### 4.2 Errors by stratum

**Table 2.** All-data fit, MAPE in % (source: `results/lit_comparison.md` / `docs/literature.md` §2.2).

{{L410-416}}

The error grows with the number of electrons and falls with ion charge; on ions of charge ≥ 3 the final model gives
1.61 %. For few-electron, highly charged ions the parameter-free exact series is already good, with a median of 0.7 %
for N ≤ 10 and 0.48 % MAPE for N/Z ≤ 0.2 with two terms. On neutral atoms the parameter-free LSDA ΔSCF solver reaches
3.30 % on the 54 first IEs with Z ≤ 54, and 0.80 % on all 171 ions with Z ≤ 18. It beats every closed form on neutral
atoms, but it is a numerical procedure and was run on only 207 rows.

Further figures are in the repository: a parity plot (`results/figures/uni_parity.png`), first IEs against Z
(`results/figures/uni_first_IE.png`), successive IEs of selected elements (`results/figures/uni_successive.png`) and
residuals (`results/figures/uni_residuals.png`).

### 4.3 Blind extrapolation (S2)

Fitted on Z ≤ 54 only, the final model predicts the 4362 ionization energies of Z = 55–110 with 6.60 % MAPE and a
1.62 % median. Its ancestors show what each ingredient contributes:

{{L435-440}}

Exact σ₁ alone did not fix extrapolation; confining the fitted part to a bounded remainder did. In a 14-parameter
ablation the bound improved the selection score from 3.96 to 2.91, and its blind S2 changed from 9.43 % to 5.53 %, a
change we saw only after the choice. For ions (Z > N) the blind S2 MAPE is 6.4 %; for neutral atoms it is 21.8 %.
Figure 1 shows the per-element median error on the S2 test rows.

Figure 2 plots the selection score against blind S2 for all fitted candidates. Candidates with a low selection score
generally extrapolate well, but the relation is not monotone. One candidate with a worse selection score (pb_exp_pos,
3.73) has a lower blind S2 (6.03 %) than the final model. We did not promote it, because that would be selection on
S2.

{{L452-452}}

**Figure 1.** Blind extrapolation. Each point is the median absolute percentage error over all ion stages of one
element Z = 55–110. The three fitted models are fitted on Z ≤ 54 only; Slater's rules and the Kregar/Di Rocco SHM have
no fitted parameters. The S2 MAPEs, computed from refits with the project's validation code, reproduce Table 1 (final
6.603 %, u35 57.233 %, pocket 11.334 %, Slater 11.638 %).

{{L459-459}}

**Figure 2.** Selection score (mean of V1, V2, S1, S3; no S2 data, although S1 and S3 contain Z ≥ 55 rows) against
blind S2 MAPE for every fitted candidate in `results/model_comparison.csv`. The pocket formula is omitted because its
V1 fit diverged in the frozen run; that divergence depends on the optimizer path (§3.4), and with the converged V1 its
selection score would be 4.19. The dashed line is Slater's rules (0 parameters).

The blind fit has three clear failures, which we report and have not corrected, since a correction now would be post
hoc. Heavy p-block neutrals come out far too low (`results/known_failures.json`, refit on Z ≤ 54): Pb at 1.20 eV
against 7.42 eV, Tl at 1.74 against 6.11 and Rn at 6.10 against 10.75. The frozen run gave 1.19 and 1.73; refits in
the two software environments differ here by about 0.01 eV. Lr comes out at −3.33 eV: the fitted r_c are negative and
the relativistic bracket is not bounded, so the Z_eff bound alone does not keep the IE positive once Zα is large.
Third IEs of lanthanides and actinides (4f/5f removal) are about 2.0–2.5 times too high (mean 2.2, median 2.1, 23
rows; worst Er²⁺ at 55.65 eV against 22.7 eV), because the Z ≤ 54 training set contains no f electrons.

### 4.4 Hydrogen-like ions: correction layers

**Table 3.** 110 H-like ions, Z = 1–110 (`docs/first_principles.md` §2.1).

{{L485-491}}

The NIST H-like reference values are themselves computed from the same QED theory, so layer (c) shows consistency
with that source rather than independent validation. The remaining 10⁻⁵–10⁻⁴ residual at high Z is the omitted
two-loop QED, nuclear-polarisation and recoil-QED terms.

### 4.5 Head-to-head with published screened hydrogenic models

We could specify two published SHMs fully from articles we were able to read: the parameter-free Kregar/Di Rocco
model [20], [22], implemented from its definitions, and the relativistic model of Mendoza et al. [16], implemented
from its published 19 × 19 table of fitted screening constants. The constants of More [11] and of Faussurier et al.
[12] are in papers we could not access, and we found no verifiable reprint of their tables. We did not reconstruct
them from memory.

#### 4.5.1 Kregar/Di Rocco SHM (parameter-free; all 5847 rows)

Our implementation follows the published definitions: screening from hydrogenic densities with an exchange
correction, iterated to self-consistency; energies E = −Σ q_i Z_i²/2n_i²; and non-relativistic, Pauli and Dirac
variants. Its same-shell Z → ∞ screening constants reproduce every printed digit (e.g. 1s 0.3125, 2p 0.3492), and its
total energies agree with the published table within 0.71 %. Its cross-shell constants differ by 0.011 on average (at
most 0.053), because the original uses fitted closed forms whose coefficients we could not obtain [19]. Valence IEs
of near-neutral ions come out 3–5 eV higher than the printed model values; for Ar I we get 18.96 eV against 14.72 eV
printed. Our near-neutral numbers therefore describe the model as defined, not the authors' code.

**Table 4.** Same 5847 rows, same scorer, MAPE in % (`results/lit_comparison.md`). For the 0-parameter models the S2
column is the error on the Z ≥ 55 rows; for the fitted models it is the blind S2 refit.

{{L525-532}}

Because of the valence discrepancy above, the 224 % neutral-atom figure probably overstates the error of the authors'
own model on neutral atoms. The parameter-free SHM works well for ionised species, with a median of 1.54 % for
charge ≥ 3. It fails near neutrality: 37 of its 5847 predictions are zero or negative, all for ions of charge 0–4 with
an open d or f shell. Where the Z → ∞ limit dominates, parameter-free and exact-σ₁ models both work. Near neutrality
every parameter-free hydrogenic model degrades, and the fitted remainder is what makes the formula useful there, at
the price of 8–33 parameters. The bounded fitted formula extrapolates better (6.6 %) than the parameter-free SHM scores
on the same Z ≥ 55 rows (13.6 %), while the unbounded fitted models extrapolate worse (57 % and 170 %).

#### 4.5.2 Mendoza et al. 2011 (published constants; 5011 covered rows)

We transcribed the relativistic nlj screening constants of [16] (the 19 × 19 matrix σ_kk′, subshells 1s½ to 5p3/2)
from Tables 1 and 2 of the authors' open-access deposit of the article (https://oa.upm.es/11165/), and implemented the
model as published: Dirac energies of screened charges, Q_k = Z − Σ_k′ σ_kk′(P_k′ − δ_kk′), and
IE = E_T(N−1) − E_T(N) with NIST ground configurations. The implementation reproduces six of the paper's printed
tables to their rounding, the 84 IEs of its Table 3 to ≤ 0.023 % and its Tables 4–8 to ≤ 0.16 %
(`models/benchmarks/mendoza2011/`). We fitted no parameter.

Four caveats apply. The authors fitted the constants with a genetic algorithm to NIST and FAC energies of
isoelectronic sequences from He to Eu, with Z up to 92, which overlaps our test rows, so no split is held out for
this model. Neutral atoms and singly charged ions were excluded from their fit. The constants were optimised for
ionization and excitation energies together, so an IE-only score is not what they were optimised for. And the tables
stop at 5p3/2, so 836 rows (Z ≥ 55, ground configurations with 5d, 5f, 6s, 6p, 6d or 7s electrons) have no
prediction. The comparison is therefore restricted to the 5011 covered rows.

**Table 4a.** The 5011 rows covered by Mendoza et al.; same scorer, all three models fitted to all data (in-sample).
Cells: mean / median absolute percentage error [n]. Source: `tools/compare_mendoza.py` → `results/compare_mendoza.md`.

{{L573-586}}

On the ions it covers, the 33-parameter formula is competitive with the published constants of Mendoza et al., not
decisively better. Its mean error is lower (1.62 % vs 2.82 %), its median is similar (0.78 % vs 0.85 %), and it has
far fewer large misses (99th percentile 11 % vs 28 %). Mendoza et al. are more accurate for few-electron ions (N ≤ 10:
0.40 % vs 0.44 %) and have the lower median for 11 ≤ N ≤ 36. Row by row, the formula is closer to NIST on 51.3 % of
the covered rows. The 9-parameter variant matches their mean on ions (2.54 % vs 2.56 %) with a larger median.

The two models differ in size and scope. Mendoza et al. publish a 19 × 19 matrix with 331 non-zero constants, fitted
by the authors; we count these as published numbers rather than independent degrees of freedom, and they are the only
model-specific numbers our implementation of their IE prescription uses. The formula has 33 global parameters and was
also validated on held-out splits. Their model also yields excitation energies and orbital properties and is used in
plasma codes, while ours gives only ionization energies. The neutral-atom row lies outside the range they fitted and
is shown for completeness. The experimental-row value here (236 rows) differs from the all-row value of §3.1 (4.63 %
on 311 rows) only because the row sets differ.

A comparison with Dirac–Fock ionization energies [30] on the same rows is left for future work.

### 4.6 The exact first-order screening constants compared with Slater's

**Table 5.** σ₁ for neutral-atom ground configurations, compared with Slater's σ. Excerpt; the full table for
N = 1–110 is in `results/fp_zexp_coefficients.csv`, with values for all 5847 rows in
`results/fp_zexp_rows_coefficients.csv`.

{{L619-631}}

E₁ is the Hund-term single-configuration value, except for N = 4 (Be), where σ₁ uses the Layzer-complex value
1.5592742. Decimal E₁ entries are exact rationals stored in the project cache, and σ₁ = −n²ΔE₁.

Slater's rules describe the neutral atom, while σ₁ is the exact screening in the bare-nucleus limit, and the two
disagree in a systematic way (Figure 3). For inner electrons the hydrogenic same-shell screening is stronger than
Slater's 0.35 per electron (Ne 2p: 6.55 against 4.15). For valence electrons of heavy atoms σ₁ is much smaller than
Slater's σ and than the screening implied by the measured IE (Na 3s: 7.79 against Slater 8.80 and an empirical 9.16;
Cs 6s: 42.8 against Slater 52.8). The fitted remainder D supplies this difference near neutrality.

{{L647-647}}

**Figure 3.** Exact first-order screening σ₁ against Slater's σ for neutral-atom configurations.

### 4.7 Worked examples and coverage

For the first IE of oxygen (1s²2s²2p⁴ → 2p³, NIST 13.618 eV), σ₁ = 5.23348 and T = 3.3496. With h = 1.76652 this
gives D = 0.71285 and Z_eff = 2.05367, which lies between Z_a = 1 and Z = 8. The hydrogenic term is 14.3457 eV, or
14.3162 eV after relativity, and the Hund term is −0.8365 eV, so IE = 13.479 eV, an error of −1.02 %.

Nitrogen (1s²2s²2p³ → 2p², NIST 14.534 eV) has a half-filled shell. Here σ₁ = 4.37867 and T = 2.9470; with
h = 1.62133, D = 0.63782 and Z_eff = 1.98351. The hydrogenic term is 13.3822 eV, or 13.3626 eV after relativity. The
Hund term is +0.8365 eV, since K_p(3) = +1.2 has the opposite sign to oxygen's, and IE = 14.199 eV, an error of
−2.31 %.

For Mg²⁺ (NIST 80.144 eV) the formula gives 78.63 eV (−1.88 %).

**Table 6.** Case studies: first IEs of six neutral atoms, with every intermediate quantity, for the three named
models. Printed by `tools/audit_components.py` and `tools/case_study_table.py` (`results/case_studies.md`); the script
asserts that the components reproduce the production code to 10⁻⁹. D is the screening remainder. The relativistic
factor is F_{n,j}·R for the σ₁ models (R = 1 for pa_bound9) and the Sommerfeld bracket for the pocket formula.

{{L675-696}}

For s removal (Na, Ca, Fe) the IE is the screened Rydberg term times relativistic factors within about 1 % of unity;
no model in this paper adds an s-type correction. The final model is the most accurate of the three on O, N, Ca, Fe
and Pb, pa_bound9 only on Na, and the pocket formula on none. Both bounded models fail on Ca (−17.8 % and −29.3 %).
For Pb, σ₁ = 63.72 (Figure 3), while the total screening σ₁ + D is 76.87.

We evaluated the public API for every Z = 1–118 and every N = 1–Z, 7021 values. All are finite and positive, with no
crashes. There are 22 violations of the monotonicity IE(Z, N−1) > IE(Z, N): 6 at Pt–Bi, caused by 4f/5s ordering,
and 16 at Rf–Ds, where the input ground configuration jumps between neighbouring ions (e.g. Rf N = 68 is 4f¹²6s² and
N = 69 is 4f¹⁴5d¹).

Above Z = 110 the predictions are qualitative. The 7p neutrals come out too low: Og is predicted at 3.66 eV, far below
its lighter congener Rn (NIST 10.75 eV), again because of the negative r_c.

---

## 5. Discussion

### 5.1 What is new and what is rediscovered

Much of the formula is known. The 1/Z expansion and its exact first-order term are Layzer's (σ₁ is his Z → ∞
screening constant [4]), computed with textbook Slater-integral algebra [33], [34]. The screened hydrogenic form goes
back to [1], [10], [11], and self-consistent (Z, N)-dependent screening to [17], [20]. The 1/(Z_a + κ) remainder is
Edlén-type isoelectronic behaviour [8]. The one-electron Dirac, recoil, finite-size and QED corrections come from
[31], [32], and LSDA ΔSCF from [27], [28].

Four things may be new. Z-expansion work [4], [5], [9] computed first-order energies for selected configurations and
isoelectronic sequences, and we are not aware of a table of exact first-order screening constants for every NIST
ground configuration with N = 1–110; we offer ours as a complete tabulation, not as a new quantity. The formula joins
exact σ₁ to a bounded fitted remainder in one closed form and is validated on all 5847 NIST successive IEs, with a
pre-registered selection score and an extrapolation from Z ≤ 54 to Z ≥ 55; we found no out-of-range extrapolation
test in the SHM literature we read. Two published SHMs are scored on the same rows: the Kregar/Di Rocco model on all
5847, and the constants of Mendoza et al. on the 5011 they cover, where the formula is competitive rather than
decisively better (§4.5). Finally, the 8-parameter pocket formula is a hand-calculable successor to Slater's rules,
with 4.68 % against 11.8 % on all rows and 11.3 % against 11.6 % on blind S2.

We do not claim a first formula for all ionization energies, since SHMs and Dirac–Fock tables already cover all ions,
or a first-principles formula, since only σ₁ and the one-electron corrections are first principles and the final model
has 33 fitted parameters. We did not search the machine-learning literature on ionization energies systematically and
make no claim relative to it.

Fitting alone interpolates to about 1.6 % and extrapolates catastrophically (170 %). Exact theory alone extrapolates
smoothly but is about 70 % off overall. Fixing the Z² and Z terms by theory and confining the fit to a bounded O(Z⁰)
remainder gives about 2 % inside the data and 6.6 % in blind extrapolation.

### 5.2 Where the formula fails

Heavy near-neutral atoms, alkaline earths and noble gases are the main failure. The blind neutral first-IE MAPE is
21.8 %, and the all-data neutral MAPE is 7.5 % (6.4 % for u35). Even with all data fitted
(`results/known_failures.json`) the ns² alkaline earths are too low (Ca −17.8 %, Sr −14.8 %, Ba −8.4 %), Rn is too low
by 24.6 %, and Na is too high by 20.7 % (Table 6). Pb's small all-data error (+8.0 %) is partly a cancellation. Its
Rydberg term alone is 9.96 eV (+34 % against NIST 7.42 eV); the relativistic bracket, 0.80 because r_p½ is negative,
lowers it to 7.97 eV (+7.4 %), and the final IE is 8.01 eV (Table 6; `results/audit_components.md`). A term of the
wrong physical sign is cancelling an overestimate. Valence screening in a neutral atom is an all-order,
non-perturbative effect, and the parameter-free LSDA ΔSCF does better on neutral atoms (3.3 % for Z ≤ 54).

The relativistic term is unbounded. The fitted r_c are all negative (−0.45 to −1.47), opposite to the direct
relativistic contraction of s and p½ electrons. Because the bracket is not bounded, the Z ≤ 54 fit gives a negative
Lr and low superheavy 7p IEs. A bounded relativistic term chosen on V1/V2/S1/S3 alone is the obvious next step; the
first attempt (§3.3) fixed the signs but lowered the selection score.

f-electron removal is the least accurate class, with 3.5 % on all data and third IEs about 2.0–2.5 times too high
(mean 2.2) in blind S2.

The formula has one Hund term per removed electron and no term-dependent multiplet energies, and ions whose ground
configuration rearranges on ionization are described by single-configuration inputs. These limits also cause the 22
monotonicity violations.

19 of the 33 parameters are shrunk class deviations. The bounded 9-parameter variant pa_bound9 (selection 2.51,
blind S2 7.80 %, all data 2.90 %; parameters in Appendix A) is a reasonable alternative with fewer parameters, though
it was not the pre-registered winner.

Most NIST values are labelled theoretical or semi-empirical (4617 + 919 of 5847), so agreement with them is partly
agreement with other theory.

### 5.3 Why no exact closed form exists for N ≥ 2

With V = Σ 1/r_ij the many-electron Schrödinger equation does not separate. E(Z) is analytic in 1/Z only up to a
critical charge; for He, 1/Z_c ≈ 1.0975. Beyond first order every coefficient E_k (k ≥ 2) is an infinite sum over the
hydrogenic continuum, with no known closed form even for two electrons, and is known only numerically [6]. A closed
formula for all ions must therefore combine exact low orders with an approximation for the rest. The choice is which
approximation, and how it is validated. We keep the exact orders exact and confine the fit to the remainder, under a
bound that prevents it from extrapolating wildly.

---

## 6. Conclusions

The screened Rydberg formula gives every successive ionization energy of every element from the configuration
alone, with 33 global parameters. On 5847 NIST values its MAPE is 1.87 % (median 0.94 %), 7.5 % for neutral atoms
and 0.0012 % for H-like ions. Fitted on Z ≤ 54, it predicts Z = 55–110 with 6.60 % MAPE, below the 11.6 % of Slater's
rules and the 13.6 % of a published parameter-free screened hydrogenic model, because first-order perturbation theory
fixes the large-Z behaviour exactly and the fitted remainder is bounded.

The formula still fails in known places. In blind S2, Pb comes out at about 1.2 eV against 7.42 and Tl at about
1.7 eV against 6.11. Lr is negative in the Z ≤ 54 fit, and Og (3.66 eV) lies far below Rn. Alkaline-earth neutrals
are too low, f-electron removal is too high in extrapolation, and all fitted relativistic coefficients are negative.

The most reusable result is the table of exact rational first-order screening constants σ₁ for N = 1–110.

---

## Data and code availability

All code, data and results are in the project repository, https://github.com/AwaisSDev/Chem-Research (folder `Np/`).
It contains the NIST table (`data/nist_ie.csv`) and the shared scorer (`evaluate.py`); the exact first-order
coefficients (`results/fp_zexp_coefficients.csv`, `results/fp_zexp_rows_coefficients.csv`); the final model and its
parameters (`models/push_a/model.py`, `models/unified/final.py`, `results/uni_final_params.json`); predictions for all
rows (`results/uni_predictions.csv`); the validation script (`models/unified/validate_blind.py`); the literature
re-implementations (`models/literature/kregar_shm.py`, `models/benchmarks/mendoza2011/`); every per-row model input
(`results/model_inputs.csv`) with the stand-alone reference implementation (`tools/verify_from_inputs.py`), which
reproduces the production predictions to 7·10⁻¹⁶; and the component audit (`tools/audit_components.py`, Table 6) and
benchmark matrix with row counts in every cell (`tools/benchmark_matrix.py` → `results/benchmark_matrix.md`). A
command-line and Python interface (`ionization.py`) evaluates the formula for any Z ≤ 118 and N ≤ Z.

## AI-assistance disclosure

The computations, code, analysis and a draft of this text were produced with AI agents (Anthropic Claude) working
under the author's direction. The author designed and directed the study and is responsible for the content. Every
number in the paper is generated by a script in the repository and can be regenerated from it. Bibliographic data
were checked against Crossref. No AI system is listed as an author.

## Acknowledgements

None.

{{TAIL}}
