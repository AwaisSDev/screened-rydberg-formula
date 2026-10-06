# Prior art, head-to-head comparison and novelty assessment

Stage: literature + head-to-head. Written 2026-10-05.
Code: `models/literature/`. Results: `results/lit_*`. Bibliography: `docs/references.bib`.

Every bibliographic entry in `references.bib` has a verification tag:
- **VERIFIED-PDF**: we read the article itself.
- **VERIFIED-WEB**: confirmed from a publisher or repository record.
- **CITED-IN**: taken from the reference list of an article we read.
- **STANDARD**: a classic citation given from memory and not re-checked online in this stage.
- **UNVERIFIED**: some field could not be confirmed.

Nothing has been invented. Check every STANDARD and UNVERIFIED field against the publisher before submission. Most
DOIs are deliberately left out for that reason.

---

## 1. Prior art

### 1.1 Empirical screening rules

**Slater (1930)** [Slater1930] gave rules for the screening constant σ and an effective quantum number n\*. The
orbital energy is −Ry(Z−σ)²/n\*². The rules were fitted by hand to atomic data and are not derived. In this project
they serve as the 0-parameter baseline. Using the total-energy difference, they give 11.8 % MAPE on all 5847 NIST
rows and 51.8 % on neutral first IEs (`results/model_comparison.md`).

**Clementi and Raimondi (1963), and Clementi, Raimondi and Reinhardt (1967)**
[ClementiRaimondi1963; ClementiRaimondiReinhardt1967] replaced Slater's hand-made σ. Their σ are orbital exponents
optimised in SCF (minimal-basis Hartree-Fock) calculations of neutral atoms, with a fitted additive formula for σ.
The constants are designed to describe orbitals, not energies. Used as an IE formula, they do worse than Slater:
31 % (one-electron form) and 120 % (total-energy form).

### 1.2 The 1/Z (Z-dependent) expansion and isoelectronic sequences

**Layzer's Z-dependent theory** [Layzer1959; LayzerEtAl1964; Layzer1967] writes the non-relativistic energy of a
fixed configuration as

E = Z²E₀ + ZE₁ + E₂ + E₃/Z + …

where:
- E₀ is a sum of hydrogenic terms;
- E₁ is the first-order Coulomb interaction, a rational combination of hydrogenic Slater integrals F^k and G^k;
- E₂ and higher terms need sums over the hydrogenic continuum.

The leading relativistic terms scale as α²Z⁴ (Layzer-Bahcall-type extensions; see the Breit-Pauli discussion in
[Pomarico2005]). High-order perturbation series for two-electron ions are classical [ScherrKnight1963;
DalgarnoStewart1958]. For many-electron isoelectronic sequences, Z-expansion codes with relativistic corrections
were developed by Safronova and co-workers (e.g. [Safronova1993]). Layzer's framework shows that a screening
constant σ = σ₀ + σ₁/Z + … has an exact Z → ∞ limit σ₀, fixed by E₁.

**This project's Track A σ₁ is exactly this σ₀.** The project calls it σ₁ = −n²ΔE₁, built from hydrogenic Slater
integrals and Racah algebra. The physics is Layzer's (1959). What may be new is the *engineering*:
- the closed rational values are tabulated for every NIST ground configuration with N = 1..110
  (`results/fp_zexp_coefficients.csv`);
- the table is used as the exact O(Z²) and O(Z) part of a formula fitted to all ions.

**Edlén** [Edlen1964] systematised isoelectronic-sequence regularities in the *Handbuch der Physik*. These include
the smooth variation of screening parameters and quantum defects with Z along a sequence, expansions in 1/(Z − s),
and their use to interpolate and extrapolate ionization energies. A referee has already pointed out that the
remainder form τ/(Z_a + κ) used in the unified models is Edlén-type behaviour. It should be presented as a
rediscovery, not as a new law (`docs/unified.md`, referee item 5).

### 1.3 Screened hydrogenic models (SHM) and average-atom codes

| model | screening constants obtained from | relativistic? | fitted? | status in this stage |
|---|---|---|---|---|
| Mayer 1947 [Mayer1947] | first SHM for opacity, n-shell screening | no | partly | historical; constants not obtained |
| More 1982 [More1982] | σ(n, m) table for n-shells; a secondary source states it was fitted to ~800 Hartree-Fock-Slater IPs of 30 elements (not checked in the original) | no | yes | **constants not obtained**: paper is paywalled, and no open source reproducing the table could be verified |
| Perrot 1989 [Perrot1989] | SHM refinements | – | – | cited only |
| Faussurier, Blancard, Decoster 1997 [Faussurier1997] | σ(nl, n′l′) with l-splitting, fitted | no | yes | **coefficients not obtained** (paywalled; SpK [Crilly2023] uses them but does not print them) |
| Faussurier, Blancard, Renaudin 2008 [Faussurier2008] | SHM with l-splitting for the equation of state of dense plasmas | – | yes | cited only |
| Rubiano et al. 2002 [Rubiano2002]; Martel et al. 1998 [Martel1998] | SHM from analytical potentials | – | partly | cited only |
| Mendoza et al. 2011 [Mendoza2011] | nlj-dependent relativistic screening constants, genetic-algorithm fit to NIST energies plus FAC-calculated IPs and excitation energies | yes | **yes** | constants transcribed from the authors' open-access deposit (oa.upm.es/11165) and validated against six printed tables; scored on the 5011 rows it covers: 2.82 % MAPE (`models/benchmarks/mendoza2011/`, paper §4.5) |
| **Kregar 1984/85; Di Rocco 1992; Pomarico, Iriarte & Di Rocco 2005; Lanzini & Di Rocco 2015; Di Rocco & Lanzini 2016** [Kregar1984; Kregar1985; DiRocco1992; Pomarico2005; LanziniDiRocco2015; DiRoccoLanzini2016] | **parameter-free**: σ computed from hydrogenic densities by splitting 1/r_ij into one-body terms, made self-consistent (Z, N)-dependent, with exchange corrections | Pauli (2005), Dirac + Breit + QED (2015/16) | **no** | **implemented and scored on all 5847 rows (§2)** |

Average-atom and HEDP kinetics codes use SHM-type or Hartree-Fock-Slater atomic data. Examples are the
relativistic Hartree-Fock-Slater average atom [Rozsnyai1972], FLYCHK [Chung2005] and SpK [Crilly2023]. In these
codes the SHM is valued for **speed and coverage of all ions**, not for spectroscopic accuracy. This is the
community whose problem the "Screened Rydberg formula" most directly addresses.

### 1.4 Ab initio benchmarks

- **Koopmans' theorem** [Koopmans1934]. A Hartree-Fock orbital energy approximates the IE but neglects orbital
  relaxation and correlation. Delta-SCF (ΔSCF: the difference of two self-consistent total energies) includes
  relaxation and is generally more accurate for atoms. We did not find, in this stage, a single verified study that
  tabulates Koopmans-HF or ΔSCF-HF IEs against NIST for all ions, so **we quote no numbers**.
- **Kohn-Sham DFT** [KohnSham1965]. NIST's LDA reference [Kotochigova1997] gives total energies of neutral atoms
  (Z = 1–92) and singly charged cations to 1 µhartree *numerical* precision. That precision is the convergence of the
  calculation, not agreement with experiment. The difference of the two energies is the first-IE ΔSCF in LDA. This
  project's own LSDA ΔSCF solver gives 1.32 % MAPE on 207 rows and 3.30 % on neutral first IEs
  (`docs/first_principles.md`), consistent with the usual LDA accuracy level.
- **Correlated non-relativistic energies.** Exact or near-exact energies for 3–18 electron ions are given by
  [Chakravorty1993]. They fix E₂ and higher terms of the 1/Z series empirically in the non-relativistic limit.
- **Dirac-Fock for all ions.** Rodrigues et al. (2004) [Rodrigues2004] tabulate Dirac-Fock total energies of ground
  configurations for the Li to Db (3–105 electron) isoelectronic series, for all Z up to 118. IEs follow from
  differences of these energies. This is the natural *ab initio* competitor for "all ions" coverage. We did not
  download the tables in this stage. **Recommended follow-up:** score Rodrigues et al. IEs on the same 5847 rows; it
  is a clean, citable comparison.
- **H-like ions.** NIST ASD [NISTASD] and the QED tabulation of Yerokhin & Shabaev (2015) [YerokhinShabaev2015].
  The project's Dirac + recoil + finite-size + one-loop QED layer reaches 0.0012 % MAPE on H-like ions. This is a
  re-implementation of standard theory and is not new.

### 1.5 Machine learning and symbolic regression

Several targeted searches in this stage found ML work on **molecular** ionization potentials and atomization
energies (e.g. [Rupp2012]). They did **not** find a verified paper that fits or discovers a closed formula for
**successive ionization energies of all atomic ions**. This is not proof that none exists. Before submission, a
Google Scholar / Web of Science search for "ionization energies" together with machine learning, symbolic regression
or neural network, restricted to atoms and ions, is required. The paper should not claim priority over ML work.

### 1.6 Accuracies as published, or measured by us on the same rows

| work | what was measured | accuracy | source |
|---|---|---|---|
| Slater 1930 | IE as a total-energy difference, all 5847 NIST rows | 11.8 % MAPE (neutral 51.8 %) | measured here (`se_slater_total`) |
| Clementi-Raimondi 1963/67 | same | 120 % MAPE (one-electron form 31 %) | measured here (`se_cr_*`) |
| Mendoza et al. 2011 | their fitting database (NIST levels + FAC IPs and excitation energies) | ~88 % of data within ±10 % | abstract of [Mendoza2011] |
| Pomarico et al. 2005 (Kregar/Di Rocco SHM) | Ar isonuclear IEs (18 values) | ratio to experiment 0.93–1.03 | Table 2 of [Pomarico2005] |
| Pomarico et al. 2005 | first IPs of closed-subshell atoms, Z ≤ 54 | stated to be within about a factor of two (their Fig. 1) | text of [Pomarico2005] |
| Kregar/Di Rocco SHM, **our implementation** | all 5847 rows | 12.5 % MAPE, median 1.71 %; charge ≥ 3: 4.68 % | measured here (§2) |
| Kotochigova et al. 1997 (LDA) | total energies, Z ≤ 92, neutral atoms and +1 ions | 1 µhartree numerical precision (not accuracy versus experiment) | [Kotochigova1997] |
| This project, LSDA ΔSCF | 207 rows | 1.32 % MAPE, neutral 3.30 % | `docs/first_principles.md` |
| Rodrigues et al. 2004 (Dirac-Fock) | total energies, 3–105 electrons, Z ≤ 118 | not scored in this stage | [Rodrigues2004] |

---

## 2. Head-to-head: published parameter-free SHM vs this project (same 5847 rows, same scorer)

### 2.1 What could and could not be implemented

- **More (1982)** and **Faussurier et al. (1997)** constants could **not** be obtained reliably. Both papers are
  paywalled, and no open document reproducing the tables could be verified. Following the honesty rules, **we did
  not reconstruct them from memory**, so there is no head-to-head against them.
- **Mendoza et al. (2011)**: constants obtained from the open-access deposit and implemented; see paper §4.5 and `models/benchmarks/mendoza2011/NOTES.md`.
- **The Kregar / Di Rocco parameter-free SHM** is fully specified by its published definitions:
  - the 1/r_ij split, Eqs. (5)–(7) of [DiRoccoLanzini2016] and Eqs. (12)–(13) of [Pomarico2005];
  - Cowan's average-of-configuration exchange [Cowan1981];
  - E = −Σ q_i Z_i²/2n_i², Eq. (7) of [Pomarico2005];
  - the Pauli correction, Eq. (15) of [Pomarico2005].

  It has **no fitted constants**, so it is the cleanest published "screened hydrogenic formula from first
  principles" against which to test this project. It is implemented in `models/literature/kregar_shm.py`.

  Implementation details and deviations (all disclosed in the code docstring):
  1. g_ij and f_ji are evaluated by exact numerical radial integration of hydrogenic densities. The paper uses fitted
     closed forms (Eq. 14) whose coefficients are only in [DiRocco1992], which we could not obtain.
  2. Eq. (15) is printed with Z⁴; we use the screened Z_i⁴.
  3. The Dirac variant averages the nlj Dirac eigenvalues over j, because NIST configurations are per nl.
  4. The Anderson-accelerated self-consistency and the warm starts are purely numerical choices; they do not change
     the fixed point.

**Validation against the published numbers** (`models/literature/validate_kregar.py`, written to
`results/lit_kregar_validation.json`):
- Diagonal Z → ∞ parameters k_ii (Table 6 of [Pomarico2005]): reproduced to all printed digits, e.g. 1s 0.3125, 2s 0.3008, 2p 0.3492, 3d 0.3765. Off-diagonal g/f: mean |Δ| = 0.011, max 0.053. The largest deviations are the 2s/2p–3d and 2s/2p–4s/4p entries, where
our values are lower (for example f(2s←1s) is 0.686 vs 0.692 printed, and 3d←2s is 0.723 vs 0.776). Our values
follow the definitions exactly. The paper instead uses fitted closed forms and an exchange apportioning that the
2005 paper does not fully specify. Some printed entries even exceed the exchange-free monopole value, which the
definition with ε ≥ 0 cannot give.
- Total configuration energies (Table 1, hartree, published/ours with the Pauli correction): He 2.85/2.85; Be 14.59/14.69; Ne 128.44/128.36; Mg 200.06/200.04; Ar 528.78/529.82; Ca 680.24/681.33; Zn 1793.04/1792.43; Kr 2787.80/2787.20; Sr 3177.06/3177.08; Cd 5580.65/5582.60; Xe 7423.77/7426.31; max deviation 0.71 %.
- Ar isonuclear IEs (Table 2): mean deviation from the published model values 4.6 %. Inner-shell IEs agree to <1 %; valence IEs are 3–5 eV higher in our implementation (Ar I 18.96 vs printed 14.72 eV; experiment 15.76 eV). Ar-like sequence (Table 3): mean deviation 11.9 %. So our reproduction matches the published model where screening is dominated by k_ii and the Z → ∞ limit. It is less faithful for valence shells of near-neutral ions, where the off-diagonal g/f (fitted closed forms in the original) matter. Treat the near-neutral numbers in §2.2 as 'the Kregar/Di Rocco model as defined', not as the authors' code.

### 2.2 Results

MAPE in %, from evaluate.py on the identical 5847 NIST rows (`models/literature/compare.py`, `results/lit_comparison.{md,json}`). Parameter-free models: S1/S2/S3 = MAPE on the test rows. Fitted project models: refit held-out values from `docs/unified.md`. `pa_hier_rel` is the frozen final model; its numbers below were re-checked after freezing by the final audit and are unchanged.


| model | n_params | ALL | median_APE_% | first_IE_neutral_atoms | hydrogen_like | N<=10 | 11<=N<=36 | N>=37 | Z>=55 | removed_l=d | removed_l=f | S1_test | S2_test | S3_test |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Kregar/Di Rocco SHM, non-rel. (literature, this work's implementation) | 0 | 13.6 | 2.34 | 230 | 6 | 4.75 | 5.23 | 24.9 | 15.1 | 4.03 | 14.1 | 13.4 | 15.1 | 13.6 |
| Kregar/Di Rocco SHM + Pauli corr. (literature) | 0 | 12.6 | 1.79 | 225 | 1.16 | 2.49 | 4.38 | 24.2 | 13.9 | 3.78 | 14.5 | 12.4 | 13.9 | 12.4 |
| Kregar/Di Rocco SHM + nl-avg Dirac (literature) | 0 | 12.5 | 1.71 | 224 | 0.19 | 2.13 | 4.29 | 24.1 | 13.6 | 3.82 | 14.3 | 12.2 | 13.6 | 12.2 |
| Slater 1930 rules, total-energy difference | 0 | 11.8 | 7.48 | 51.8 | 6 | 4.89 | 9.38 | 16.8 | 11.6 | 16.2 | 12.9 | 11.6 | 11.6 | 12.4 |
| Clementi-Raimondi 1963/67, total-energy difference | 0 | 120 | 18.7 | 2281 | 6 | 5.95 | 32.4 | 245 | 134 | 85.1 | 136 | 116 | 134 | 119 |
| Track A zexp (exact O(Z^2,Z) + heuristic O(1) + rel/QED) | 0 | 69.8 | 7.2 | 1105 | 0.0012 | 1.71 | 15.4 | 147 | 77.8 | 57.6 | 126 | 68.3 | 77.8 | 74.3 |
| Project: pocket formula | 8 | 4.68 | 3.04 | 16.8 | 1.16 | 2.37 | 2.76 | 7.38 | 4.7 | 3.67 | 11.5 | 4.62 | 11.3 | 5.24 |
| Project: Track B GSHM final | 32 | 1.57 | 0.737 | 7.61 | 0.19 | 0.698 | 1.16 | 2.31 | 1.48 | 0.982 | 2.67 | 1.59 | 170 | 1.7 |
| Project: u35 (exact sigma1 + fitted remainder) | 35 | 1.48 | 0.728 | 6.43 | 0.0012 | 0.463 | 1.08 | 2.28 | 1.58 | 1.17 | 2.58 | 1.52 | 57.2 | 1.65 |
| Project: pa_hier_rel (final model) | 33 | 1.87 | 0.94 | 7.52 | 0.0012 | 0.439 | 1.43 | 2.87 | 1.92 | 1.52 | 3.51 | 1.91 | 6.6 | 2 |


**Reading the table.**
- Best literature variant: Kregar/Di Rocco SHM + nl-avg Dirac (literature). It reaches 12.5 % all-data MAPE (median 1.71 %), 224 % on neutral first IEs and 0.19 % on H-like ions. Caveat on the neutral figure: our implementation is not fully faithful for valence shells (Ar I 18.96 eV against the 14.72 eV printed by the authors; the authors' printed valence IEs are 3–5 eV lower than ours, §2.1), so 224 % probably overstates how badly the authors' own model does on neutral atoms. Its blind-equivalent S2 (Z ≥ 55) score is 13.6 %.
- Where the literature SHM fails:
  - neutral first IEs: median error 178 %;
  - "rearranged" rows: 164 % MAPE;
  - 37 of the 5847 predictions are non-positive. All are ions of charge 0–4 with an open 3d, 4d, 4f or 5f shell,
    for example Ni II, Pd I, Er III, Cf III and Lu I. Our interpretation, not tested further: the open d/f shell is
    over-bound relative to the outer s electrons by nodeful hydrogenic orbitals with a large Z_i, so the
    configuration-average energy difference changes sign.
- Where it works: on ions with charge ≥ 3 its MAPE is 4.68 % (median 1.54 %), against 10.2 % for Slater and 1.61 %
  for `pa_hier_rel`. On N ≤ 10 it gives 2.1 %.
- Slater's rules on the same rows: 11.8 % all-data, 51.8 % neutral, but a 7.5 % median. Slater's MAPE is slightly
  lower than the SHM's only because Slater never blows up near neutrality; the SHM is far better on typical
  (ionized) rows.
- Project pocket formula (8 parameters): 4.68 % all-data, blind S2 11.3 %.
- Project final model `pa_hier_rel` (33 parameters): 1.87 % all-data, 7.52 % neutral, blind S2 6.6 %.
- A parameter-free published model cannot be compared with a fitted one on parameter count alone. The fair statements are:
  1. on highly charged ions (N/Z small), the parameter-free SHM and the project's exact-σ₁ models both work, because the Z → ∞ limit is built in;
  2. on near-neutral atoms, every parameter-free hydrogenic model degrades. The fitted remainder is what makes the project's formula useful there, at the cost of 8–35 parameters;
  3. the project's blind S2 numbers come from extrapolating a fit, while the literature model has nothing to extrapolate. The comparison tests whether fitting buys more than it risks. Against the literature model's 13.6 % on the same Z ≥ 55 rows: `pa_hier_rel` (6.6 %) is better; pocket (11.3 %) is better; `u35` (57.2 %) is worse; GSHM (170 %) is worse.


---

## 3. Novelty assessment

### Re-derivations of known physics (must be credited, not claimed)

1. **The 1/Z expansion and its exact first-order term.** σ₁ = −n²ΔE₁ is Layzer's (1959) Z → ∞ screening constant.
   Hydrogenic F^k/G^k algebra for E₁ is textbook material [Cowan1981; FroeseFischer1997].
2. **The screened hydrogenic form** IE ≈ Ry(Z − σ)²/n². This goes back to Slater (1930), Mayer (1947) and More
   (1982); self-consistent (Z, N)-dependent screening is Kregar/Di Rocco.
3. **Edlén-type isoelectronic behaviour.** The remainder τ/(Z_a + κ) and the smooth screening along sequences are
   Edlén-type regularities.
4. **Dirac + recoil + finite-size + QED for H-like ions.** This is standard theory (NIST ASD; Yerokhin & Shabaev
   2015).
5. **LSDA ΔSCF.** This is standard DFT (Kotochigova et al. 1997 is the NIST reference implementation).

### Plausibly new (claims to make carefully)

1. **An exact first-order screening table for every NIST ground configuration, N = 1..110**, presented as an
   *ab initio replacement for Slater's rules*. The values are rational and parameter-free, and they carry the exact
   Z → ∞ limit. We found no published table with this coverage. Isoelectronic Z-expansion work, such as Safronova's,
   covers selected sequences and levels; the SHM literature uses fitted or iterated σ. **Caveat:** this must be
   checked against the Z-expansion literature (Safronova et al., Layzer et al. 1964) before claiming it.
2. **One hybrid closed-form formula** (exact σ₁ plus a small fitted remainder, 33 parameters in the final
   model `pa_hier_rel`) **validated on all 5847 NIST successive IEs with pre-registered held-out tests**: S1, S3,
   V1, V2 and a *blind* extrapolation from Z ≤ 54 to Z ≥ 55. In the SHM papers we read, accuracy is reported as
   "fraction within ±10 %" on a fitting database (Mendoza 2011: ~88 %) or on selected sequences (Pomarico 2005).
   We saw no blind out-of-range extrapolation test. The S2 protocol and the reporting of every fitted parameter are
   a methodological contribution.
3. **A quantitative head-to-head of a published parameter-free SHM on all 5847 rows** (§2). We are not aware of a
   published all-ions NIST score for the Kregar/Di Rocco SHM. Results:
   - our implementation of it gives 12.5 % MAPE (median 1.7 %) and 13.6 % on the Z ≥ 55 rows;
   - the project's final model gives 1.87 % all-data and 6.6 % *blind* S2;
   - the 8-parameter pocket formula gives 4.68 % and 11.3 %.

   So the project's fitted formulas beat our implementation (from the published definitions) of the best parameter-free
   published SHM we found, including on blind
   extrapolation. The exception is the earlier `u35` and GSHM models, whose blind S2 (57 %, 170 %) is far worse than
   the SHM's 13.6 %. That contrast is itself a useful, honest result.
4. **The 8-parameter "pocket" formula**: a hand-calculable rule with 4.7 % all-data MAPE and 11 % blind S2, against
   Slater's 11.8 %. It is a candidate *pedagogical* "Slater's rules upgrade" for chemistry education.

### Not new and must not be claimed

- "First formula for all ionization energies": SHMs (More 1982; Faussurier 1997; Mendoza 2011; Kregar/Di Rocco) and
  Dirac-Fock tables (Rodrigues 2004) already cover all ions.
- "First-principles formula": the final model has 33 fitted parameters. Only σ₁ and the H-like layers are first
  principles.

---

## 4. Recommended venues

| venue | why | risk |
|---|---|---|
| **Journal of Physics B: Atomic, Molecular and Optical Physics** | Natural home for atomic-structure methodology; accepts semi-empirical formulas with careful validation. | Referees will demand comparison with Dirac-Fock (Rodrigues 2004) and with the SHM literature: the §2 head-to-head plus the Rodrigues follow-up. |
| **High Energy Density Physics** or **JQSRT** | The SHM community (More, Faussurier, Mendoza, Di Rocco) publishes here. Fast, all-ion IEs are exactly what plasma kinetics and opacity codes need. | Must show value for plasma use (speed, smoothness, coverage up to Z = 118) and compare with the community's SHMs. |
| **Atomic Data and Nuclear Data Tables** | If the centrepiece is the σ₁ table for all configurations plus formula predictions for all ions with Z ≤ 118. | ADNDT expects comprehensive tables and an error analysis; the formula alone is a weaker fit. |
| **Atoms (MDPI, open access)** | Fast review, open access, receptive to atomic-data methods papers. | Lower prestige. |
| **Journal of Chemical Education** (pocket formula only) | "Beyond Slater's rules" teaching paper: an 8-parameter hand-calculable rule, 2.5× better than Slater on all ions. | Must be written for teachers, not as a research claim. |

Physical Review A is **not** recommended unless the σ₁ table turns out to be new *and* is presented with new physics
insight. The fitted-remainder formula alone is unlikely to meet PRA's novelty bar.

---

## 5. Required disclosures for the paper

1. **AI assistance.** Code, analysis and draft text were produced with an AI assistant (Anthropic Claude) under the
   author's direction. Most publishers (IOP, APS, Elsevier, MDPI) require a statement in Methods or
   Acknowledgements, and an AI cannot be listed as an author. Check the target journal's current policy and word the
   statement accordingly. The author is responsible for verifying every reference: see the tags in `references.bib`.
2. **Data provenance.** `data/nist_ie.csv` is from NIST ASD and includes values flagged experimental, semi-empirical
   and theoretical (`status` column). Report the ASD version and access date, and give metrics on the experimental
   subset too.
3. **Fitted-parameter count.** State it for every model. "Parameter-free" applies only to σ₁, the H-like layers,
   LSDA ΔSCF and the literature Kregar SHM.
4. **Design decisions influenced by seen results.**
   - The post-hoc σ_core bound in Track B was chosen after seeing S2. It is reported only as non-blind.
   - Push A's starting point was informed by the pocket formula's known S2 result (`docs/unified.md`).
   - In this literature stage, the choice of the Kregar/Di Rocco SHM as the head-to-head model was made *before*
     scoring it, because it was the only published SHM whose full specification could be obtained. No variant was
     selected on the basis of its score. All three variants (nr, Pauli, Dirac) are reported.
5. **Reproduction caveat** for the literature model (§2.1). Our Kregar SHM follows the published definitions. It is
   not the authors' own code, and its off-diagonal screening parameters differ from their printed Z → ∞ table.
