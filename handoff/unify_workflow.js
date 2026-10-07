export const meta = {
  name: 'ionization-energy-unify',
  description: 'Unify the two finished ionization-energy tracks into one formula while a referee audits both, then one final fix pass',
  phases: [
    { title: 'Unify + Referee', detail: 'unified hybrid model || adversarial referee of both tracks' },
    { title: 'Final fix', detail: 'apply referee findings, smoke-test, literature positioning' },
  ],
}

const ROOT = 'E:/Dr. Navin/Np'
const SCRATCH = 'C:/Users/HP/AppData/Local/Temp/claude/E--Dr--Navin-Np/43fc5e11-e216-4a38-838f-413efbd43896/scratchpad'

const SHARED = `PROJECT: Derive/approximate a generalizable analytical equation for the successive ionization energies IE(Z, N) of any atom or ion, taking atomic number Z and the electron configuration as input, as accurately as possible. Bohr's IE = 13.6 Z^2/n^2 eV is exact only for one-electron ions (and only non-relativistically); Slater's rules are empirical. We want a model derived from first principles and/or by pattern recognition on experimental data, incorporating shielding, electron-electron repulsion, quantum corrections (exchange, correlation, relativistic, QED), and density-functional theory where useful.

The work is split between two parallel tracks:
 - Track A (first principles): exact 1/Z perturbation expansion + relativistic/QED corrections + a home-made Kohn-Sham DFT Delta-SCF atomic solver.
 - Track B (data-driven): baselines (Bohr, Slater, Clementi-Raimondi) + a new generalized screened-hydrogenic formula discovered by pattern recognition / sparse regression on the NIST data, with held-out validation.
A later unification stage will merge them (likely: exact leading coefficients from Track A + fitted higher-order screening from Track B), so make your key functions importable and well documented.

ENVIRONMENT
 - Project root: ${ROOT} (Windows 10; Git Bash and PowerShell tools available). Not a git repo.
 - Python: ALWAYS run with "py -3.13" (has numpy 2.2, scipy 1.15, pandas 2.3, matplotlib 3.10, sympy 1.14, scikit-learn 1.7, numba 0.61). Do not create venvs or pip install. Use matplotlib's Agg backend.
 - The machine has only 4 CPU cores, shared with one other agent running concurrently: use at most 2 processes in parallel; prefer vectorized numpy / numba; keep any single long job under ~40 minutes (run long jobs in the background and poll, checkpoint partial results to disk).
 - Data: data/nist_ie.csv = NIST ASD successive ionization energies, 5847 rows, Z = 1..110, all ion stages. Columns include Z, N, IE_eV, status (experimental | semi-empirical | theoretical; NIST lists most highly-charged-ion values as Dirac-Fock+QED theoretical values, which are accurate references), shells_expanded (e.g. 1s2.2s2.2p6.3s1).
 - Loader: common/atomdata.py. load_records() returns a tuple of dicts with keys Z, N, ion_charge, symbol, shells [(n,l,occ),...], IE_eV, status, removed = (n,l) of the subshell losing the electron (from diffing the N and N-1 ground configurations), rearranged (True for the 63 ions whose ground configuration reshuffles on ionization, e.g. Ni 3d8 4s2 -> Ni+ 3d9). Also get(Z,N), ground_shells(Z,N) (NIST configuration if in table, else Madelung), madelung_shells(N), constants HARTREE_EV, RYDBERG_EV, ALPHA, ELECTRON_MASS_U.
 - Scoring: "py -3.13 evaluate.py results/<name>_predictions.csv --json results/<name>_metrics.json". The predictions CSV needs columns Z,N,IE_pred_eV. It prints MAPE / median / p90 / max relative error and MAE overall and by strata (neutral-atom first IEs, H-like, removed l, removed n, Z range, N range, status, rearranged). For held-out evaluation, import evaluate and call evaluate.evaluate(pred_csv, subset=lambda rec: ...).
 - Reference point: naive Bohr IE = Ry Z^2/n^2 (n of removed electron) gives overall MAPE 1365% and neutral-atom first-IE MAPE 22878%; even the 110 H-like ions have 6% MAPE because of relativity at high Z.

RULES
 - Do NOT modify data/, common/, or evaluate.py (shared; the other agent depends on them). Put helpers in your own directory.
 - Write only inside your own track directory (given below), results/ (files prefixed with your track prefix), results/figures/ (prefix too), and docs/ (your own file).
 - Units: eV for outputs; Hartree internally is fine. Relative error is the headline metric (IEs span ~4 eV to ~2e5 eV), so fit in log space or on relative residuals.
 - Honesty over hype: report real numbers, count every fitted parameter, no per-element lookup tables disguised as formulas, no reading NIST IE values inside predict(), no test-set leakage. If something does not work, say so and show the numbers.
 - Standard held-out splits for anything fitted (use exactly these so the tracks are comparable):
     S1 interpolation in Z: test = records with Z % 5 == 0; train = rest.
     S2 extrapolation in Z: train = Z <= 54; test = Z >= 55.
     S3 unseen isoelectronic sequences: test = records with N % 6 == 0; train = rest.
   Report test-set MAPE, median APE, and neutral-atom first-IE MAPE for each split.
 - Every model must be callable as predict(Z, N, shells=None) -> IE in eV (shells defaults to common.atomdata.ground_shells(Z, N)), so it works for ANY (Z, N), not only rows in the table.
 - The derivation document is the main intellectual product: write it in Markdown with LaTeX math ($...$ and $$...$$), define every symbol, explain every coefficient, include result tables and figure paths. Figures go to results/figures/<prefix>*.png.
 - Work autonomously to completion; do not ask questions. Your final message is returned as structured data to the orchestrator.`

const TRACK_A = `YOU ARE TRACK A - FIRST PRINCIPLES.
Your directory: models/first_principles/   Results prefix: fp_   Doc: docs/first_principles.md

Goal: derive, as far as rigorously possible, a closed-form IE(Z, N, configuration), and quantify exactly how far pure theory gets before empirical input is needed.

1. Exact 1/Z (Hylleraas-Layzer) perturbation expansion. Scaling r -> r/Z gives H = Z^2 [H0 + Z^-1 V_ee], so for a fixed configuration/term E(Z,N) = Z^2 E0 + Z E1 + E2 + E3/Z + ... (Hartree), and
     IE(Z,N) = E(Z,N-1) - E(Z,N) = Z^2 dE0 + Z dE1 + dE2 + ...
   - E0 = -sum_i q_i/(2 n_i^2); so dE0 = 1/(2 n^2) for the removed electron: this term IS Bohr's formula.
   - E1 = <V_ee> over hydrogenic Z=1 orbitals. Build it EXACTLY from hydrogenic radial Slater integrals F^k(nl,n'l') and G^k(nl,n'l') (rational numbers; compute exactly with sympy and cache to disk) and the angular coefficients: Slater configuration-average energy (Cowan, Theory of Atomic Structure and Spectra, ch. 6) plus the term correction for the Hund's-rule ground LS term of open subshells (at least for one open s/p/d/f subshell; for two open subshells use configuration average + a documented approximation, e.g., high-spin coupling). Discuss Layzer complex degeneracy (same-n states degenerate in the hydrogenic limit, e.g. 1s2 2s2 mixes with 1s2 2p2 at first order for Be-like ions) and handle it where it matters (Be-like, Mg-like ...) by diagonalising within the complex if feasible.
   - Validate E1 against known exact values: He-like 1s2: E1 = 5/8 = 0.625; Li-like 1s2 2s: E1 = 5965/972 = 6.136831; He-like E2 = -0.157666 (literature); check others you can find.
   - Result: an analytic formula IE = Ry [ Z^2/n^2 - 4 dE1' Z + ... ] (be careful with Hartree vs Rydberg factors) with exact Z^2 and Z coefficients for every configuration: a first-principles derivation of screening. The implied leading-order screening constant sigma (IE ~ Ry (Z - sigma)^2/n^2 to first order in 1/Z) follows from dE1; tabulate it and compare with Slater's empirical sigma.
   - Second order dE2 is not closed-form in general (needs sums over the continuum). Explore: (i) exact/literature values for small N; (ii) a physically motivated closed-form estimate; (iii) one fitted parameter per isoelectronic sequence N (clearly labelled semi-empirical, validated on the S1 split, and S3 is then impossible: say so) and/or a smooth universal function of the configuration (validated on S1/S2/S3).
   - Quantify where the truncated series (2 terms, 3 terms) works (by N, by Z/N ratio). The expansion is asymptotic in N/Z; near-neutral heavy atoms are hard. Show this with numbers and plots (IE along isoelectronic sequences vs NIST).
2. Relativistic + QED corrections. Implement the Dirac-Coulomb hydrogenic energy E(n,kappa) for the removed electron with an appropriate (effective) charge and/or the alpha^2 Z^4 Breit-Pauli expansion (mass-velocity, Darwin, spin-orbit), using the j of the ground level (NIST ground_level designation is available in the CSV) where useful, finite-nuclear-mass (reduced-mass) correction, and an estimate of the leading QED self-energy (Lamb shift) for s electrons of high-Z ions (e.g., Mohr's F(Z alpha) tabulation or a fitted Mohr-type formula, cited). Show the H-like sequence (110 rows, Z up to 110) is reproduced to <= 1e-4 relative error (state what you achieve) vs 6% MAPE for naive Bohr. Then combine with the 1/Z expansion for He-like, Li-like ... ('relativistic Z-expansion'): quantify how good it is for few-electron highly-charged ions.
3. Density-functional theory. Write your own radial (spherically averaged) Kohn-Sham solver for atoms and ions: logarithmic radial grid, robust eigen-solver (Numerov shooting with node counting, or finite-difference/Chebyshev generalized eigenproblem), Hartree potential by Poisson integration, LDA/LSDA exchange-correlation (Slater exchange + VWN5 correlation, and PZ81 if you like; spin-polarised for open shells with Hund's-rule occupations; spherically averaged fractional occupations allowed), robust SCF mixing (Anderson/Broyden, or linear with adaptive damping). Compute Delta-SCF IE = E_KS(N-1) - E_KS(N) for each ion; also the Slater transition-state / Janak estimate IE ~ -eps_HOMO at half removed occupation, a near-closed-form link to orbital energies.
   - Validate total energies against NIST 'Atomic Reference Data for Electronic Structure Calculations' (Kotochigova, Levine, Shirley, Stiles, Clark), non-spin-polarised LDA (VWN): He -2.834836, Be -14.447209, Ne -128.233481, Ar -525.946195 Hartree; aim for ~1e-5 Ha agreement. WebFetch/WebSearch of NIST pages is allowed for more reference numbers.
   - Add a scalar-relativistic option (Koelling-Harmon or ZORA-type radial equation) if feasible so heavy-element inner-shell IEs are sensible; validate against the NIST ScRLDA reference energies if you implement it.
   - Run Delta-SCF for all ions with Z <= 36 (666 rows) at minimum, and all neutral first IEs Z <= 86 (or as far as runtime allows); extend further if fast. If time permits, implement one improvement that matters for IEs (Perdew-Zunger self-interaction correction for the removed orbital, or PBE GGA) and quantify the gain.
   - Predictions: results/fp_dft_predictions.csv (rows you computed; say which), results/fp_zexp_predictions.csv (all 5847 rows, the analytic expansion model including relativistic corrections), plus any variants (e.g. fp_zexp2_, fp_zexp_rel_). Score each with evaluate.py and save metrics JSONs.
4. docs/first_principles.md: complete derivation (Hamiltonian scaling, perturbation series, Slater integral formulas, angular coefficients, term corrections, complex degeneracy, relativistic and QED formulas, KS equations and functionals); table of exact E0, E1, dE1 and the implied ab-initio screening constants for ground configurations N = 1..36 in the doc, with the full N = 1..110 table in results/fp_zexp_coefficients.csv; accuracy tables; plots; and an honest discussion: what is truly first-principles, what needed empirical input, why no exact closed form exists for N >= 2 (non-separable Coulomb many-body problem), and what the best closed-form approximation is.`

const TRACK_B = `YOU ARE TRACK B - DATA-DRIVEN / PATTERN RECOGNITION.
Your directory: models/semi_empirical/   Results prefix: se_   Doc: docs/semi_empirical.md

Goal: the most accurate, compact, generalizable closed-form formula for IE(Z, N, configuration), discovered from patterns in the NIST data, physically interpretable (screening, penetration / quantum defect, exchange / pairing, relativity), and validated out of sample.

1. Baselines on all 5847 rows (score each with evaluate.py, save predictions): (a) Bohr Ry Z^2/n^2; (b) Slater's rules Z_eff with Slater's n* (1, 2, 3, 3.7, 4.0, 4.2), in two variants: one-electron IE = Ry (Z_eff/n*)^2, and the proper total-energy difference IE = E_S(N-1) - E_S(N) with E_S = -Ry sum_i q_i (Z_eff,i/n*_i)^2; (c) Clementi-Raimondi (1963) screening-constant rules (their empirical formula for sigma in terms of subshell populations; extend sensibly beyond its original range and say how), same two variants. These quantify how bad textbook approaches are on SUCCESSIVE IEs.
2. Data exploration / pattern recognition, with figures: IE vs N at fixed Z (shell jumps); IE vs Z along isoelectronic sequences (Moseley-type plots: n*sqrt(IE/Ry) is nearly linear in Z; extract slope and intercept, i.e. 1/n_eff and screening sigma per sequence, and study how they depend on the configuration); half-filled/full-subshell kinks (p3/p4, d5/d6, f7/f8: exchange / spin-pairing energy); 4s/3d and 6s/5d/4f competition; relativistic growth at high Z (residual vs (Z alpha)^2).
3. Build a new model. Suggested backbone (improve freely):
     IE(Z,N) = Ry (Z - S)^2/(n - delta)^2 * F_rel + E_x
   with S = sum_j sigma(i<-j) q_j over all other electrons (including the others in the same subshell), sigma depending only on the relative (n,l) of the removed electron i and the screening electron j (inner shell, same n lower l, same subshell, same n higher l, outer shells...), delta a quantum defect depending on l (and possibly on the degree of ionization q/Z), F_rel a Dirac-type relativistic factor in Z alpha, E_x an exchange / spin-pairing term (e.g., proportional to (Z - S)/n^2 times the change in the number of parallel-spin pairs under Hund's rule; this should capture the p3/p4, d5/d6 kinks). Also try a total-energy-difference form IE = E(N-1) - E(N) with screened-hydrogenic total energies (this includes orbital relaxation automatically). Compare forms.
   Then do residual analysis and sparse model discovery: build a library of physically motivated candidate terms (powers of Z and of (Z - S), 1/n, 1/n^3, occupancy features, q/Z, l indicators, half/full-shell indicators, (Z alpha)^2 Z^2/n^3, ...) and select terms with cross-validated sparse regression (OMP / LASSO / forward selection on the train split) to minimise held-out relative error. Keep the final model compact (target <= 40 global parameters; no per-element or per-ion parameters; report the count) and fit with scipy least_squares on log(IE) or relative residuals.
4. Validation on S1/S2/S3 (fit on train only; report test metrics) plus a learning curve (number of parameters vs held-out error) to show it is not over-fitted. Also refit on all data for the final parameters and report full-dataset metrics from evaluate.py. Look specifically at the hard cases: neutral-atom first IEs (alkali minima, noble-gas maxima, transition metals, lanthanides), rearranged configurations, highly-charged heavy ions.
5. Deliverables: models/semi_empirical/*.py with predict(Z, N, shells=None); parameters in results/se_params.json (names, values, uncertainties if possible); predictions results/se_<variant>_predictions.csv (all 5847 rows each); metrics JSONs; figures results/figures/se_*.png; docs/semi_empirical.md with the final formula written out in full (every symbol, every parameter value in a table), physical interpretation of each term, the discovery process (which patterns were found and what each term buys: an ablation table), baseline comparison table, held-out results, honest limitations.`

const BUILD_SCHEMA = {
  type: 'object',
  properties: {
    track: { type: 'string' },
    headline_formula: { type: 'string', description: 'The main closed-form result, plain text / LaTeX, with symbols defined' },
    summary: { type: 'string', description: 'What was built and found (8-20 sentences)' },
    files: { type: 'array', items: { type: 'string' } },
    metrics_markdown: { type: 'string', description: 'Markdown table: each model variant x {overall MAPE, median APE, neutral first-IE MAPE, H-like MAPE, S1/S2/S3 test MAPE where applicable, n_params}' },
    n_fitted_parameters: { type: 'integer' },
    limitations: { type: 'array', items: { type: 'string' } },
    notes_for_other_track: { type: 'array', items: { type: 'string' } },
  },
  required: ['track', 'headline_formula', 'summary', 'files', 'metrics_markdown', 'n_fitted_parameters', 'limitations', 'notes_for_other_track'],
}

const REVIEW_SCHEMA = {
  type: 'object',
  properties: {
    verified_claims: { type: 'array', items: { type: 'string' } },
    issues: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          severity: { type: 'string', enum: ['critical', 'major', 'minor'] },
          title: { type: 'string' },
          evidence: { type: 'string' },
          suggested_fix: { type: 'string' },
        },
        required: ['severity', 'title', 'evidence', 'suggested_fix'],
      },
    },
    reproduced_metrics: { type: 'string' },
    overall_verdict: { type: 'string' },
    notes_for_unification: { type: 'array', items: { type: 'string' } },
  },
  required: ['verified_claims', 'issues', 'reproduced_metrics', 'overall_verdict', 'notes_for_unification'],
}

const FIX_SCHEMA = {
  type: 'object',
  properties: {
    ...BUILD_SCHEMA.properties,
    issues_addressed: { type: 'array', items: { type: 'string' }, description: 'issue title -> what was done (or rebuttal)' },
  },
  required: [...BUILD_SCHEMA.required, 'issues_addressed'],
}

const TRACKS = [
  { key: 'A', name: 'first-principles', brief: TRACK_A, dir: 'models/first_principles/', prefix: 'fp_', doc: 'docs/first_principles.md' },
  { key: 'B', name: 'data-driven', brief: TRACK_B, dir: 'models/semi_empirical/', prefix: 'se_', doc: 'docs/semi_empirical.md' },
]

// Build results are passed in from the main session (the two track agents ran outside this workflow)
const builds = args.builds
if (!Array.isArray(builds) || builds.length !== 2) throw new Error('args.builds must be [trackA, trackB]')

// ---------------- FAST MODE: Unify || Referee, then one Final fix ----------------
const FAST = `WORKING DIRECTORY: first run cd to the project root E:/Dr. Navin/Np; all relative paths in this brief are relative to it. SPEED MODE: the user wants results fast. Hard time-box: finish within about 35 minutes. No long compute jobs (nothing over ~5 minutes); reuse the tracks' existing predictions CSVs, caches and fitted parameters wherever possible instead of recomputing. Prefer a working, honest, well-documented result over extra variants.`
const LI_NOTE = `Note: the correct exact first-order energy of Li-like 1s2 2s is E1 = 5965/5832 = 1.022805 Hartree (= 5/8 + 2*17/81 - 16/729); an earlier brief mistakenly said 5965/972, which is a typo, not a defect in the code.`

phase('Unify + Referee')
const [unified, referee] = await parallel([
  () => agent(`${SHARED}

${FAST}

UNIFICATION STAGE. Both tracks have finished their builds. Read docs/first_principles.md, docs/semi_empirical.md, the code in models/first_principles/ and models/semi_empirical/, and results/.
Track A summary (JSON):
${JSON.stringify(builds[0], null, 1)}
Track B summary (JSON):
${JSON.stringify(builds[1], null, 1)}
${LI_NOTE}
An independent referee is reviewing both tracks concurrently; a short fix pass follows you, so leave the code clean and the doc accurate.

Your directory: models/unified/   Results prefix: uni_   Doc: docs/unified.md. You may import (not modify) both tracks' modules.

Build the final unified model and a top-level API ionization.py at the project root:
  py -3.13 ionization.py Fe            -> table of all 26 successive IEs, predicted vs NIST (with % error)
  py -3.13 ionization.py 26 --N 26     -> a single value
  py -3.13 ionization.py 118           -> works beyond the table (Madelung configuration), clearly flagged as a prediction
  from ionization import ionization_energy, successive_ionization_energies
Aim: the most accurate single closed-form equation that is maximally first-principles: e.g. keep the exact Z^2 and Z coefficients (Bohr term + exact first-order electron-electron repulsion dE1, i.e. ab-initio screening) and the Dirac/QED corrections from Track A, and model only the remainder (second-order correlation / relaxation) with Track B's compact fitted terms, so the fitted part is small and smooth. Compare against simply using Track B's best model (+ Track A's relativistic factor if B lacks one); pick the better on held-out S1/S2/S3 and say so honestly. Report the parameter count. Also give a 'pocket formula' (<= 10 parameters, hand-calculable) with its accuracy.
Deliver: results/uni_predictions.csv (all 5847 rows) + uni metrics JSON; results/model_comparison.csv and results/model_comparison.md comparing all models from both tracks + baselines + unified with the same columns (overall MAPE, median APE, neutral first-IE MAPE, H-like MAPE, S1/S2/S3 test MAPE where applicable, n_params, first-principles vs fitted); figures results/figures/uni_parity.png (log-log), uni_first_IE.png (first IE vs Z, predicted vs NIST, Z = 1..110), uni_successive.png (successive IEs for C, Fe, Xe, U), uni_residuals.png (relative residual vs Z coloured by removed l); docs/unified.md with the final equation(s) in full, parameter table, derivation lineage (which term came from where), validation, limitations, and a worked hand example (1st IE of O).`,
    { label: 'unify', phase: 'Unify + Referee', schema: BUILD_SCHEMA, effort: 'high' }),
  () => agent(`${SHARED}

${FAST} Your own time-box is about 25 minutes.

You are a single ADVERSARIAL REFEREE for BOTH tracks (a unification agent is working concurrently in models/unified/; ignore that directory). Find what is WRONG, overstated, irreproducible or leaky, and confirm what is right. A claim is unverified until you reproduce it.
Track A (first principles: models/first_principles/, docs/first_principles.md, results/fp_*) summary (JSON):
${JSON.stringify(builds[0], null, 1)}
Track B (data-driven: models/semi_empirical/, docs/semi_empirical.md, results/se_*) summary (JSON):
${JSON.stringify(builds[1], null, 1)}
${LI_NOTE}
Check, quickly but for real:
 - Re-score the main predictions CSVs of each track with evaluate.py and compare with the claimed numbers.
 - Leakage: fitted parameters seeing test-split rows (S1: Z%5==0 test; S2: train Z<=54; S3: N%6==0 test); per-element/per-ion parameters; NIST IE lookups inside predict(); predict(Z,N) for out-of-table inputs (Z=118 N=118, Z=120 N=2).
 - Physics spot checks with your own short calculations: E1(He-like) = 5/8, E1(Li-like) = 5965/5832, a Dirac 1s energy for Z=92, one or two NIST LDA total energies vs the DFT solver (He -2.834836, Ne -128.233481 Ha) if the solver is fast, Slater's rules on a textbook example, Hartree/Rydberg factor-of-2 consistency.
 - Doc claims vs actual numbers; anything labelled first-principles that contains fitted numbers.
 - VALIDATION-CLAIM AUDIT (mandatory, state answers explicitly in verified_claims or issues):
   (a) DFT: which reference energies is Track A's solver compared against (expected: NIST Atomic Reference Data, Kotochigova et al., LDA/LSD with Slater exchange + VWN5, i.e. the SAME functional) and with which functional? Confirm that 1e-6 Ha agreement is framed only as implementation verification, not physical accuracy; check that LDA's physical error is reported (total energy vs exact non-relativistic, e.g. Ne exact -128.9376 Ha vs LDA -128.2335 Ha; Delta-SCF IE errors vs NIST experiment). Look for unit / reference-energy bugs (Hartree vs Rydberg, spin-polarised vs unpolarised reference).
   (b) H-like circularity: Track A's QED uses per-Z self-energy F(Z alpha) values and nuclear radii tabulated by Yerokhin & Shabaev 2015 (models/first_principles/cache/yerokhin_shabaev_2015_qed.json). Determine whether NIST's H-like reference IEs derive from the same theory; if so the ~12 ppm H-like agreement is consistency, not independent validation. Require accuracy reported in layers: Dirac only; + recoil + finite nuclear size (analytic); + tabulated QED. Count the per-Z tabulated inputs honestly as external theory inputs.
 - Quick literature positioning (one or two WebSearch calls at most): what is known (Slater 1930, Clementi-Raimondi 1963, Layzer 1959 1/Z expansion, screened-hydrogenic models of Mayer/More/Faussurier, Delta-SCF DFT accuracy) and whether anything here is genuinely new. Put it in notes_for_unification.
Do NOT edit project files; scratch only under ${SCRATCH}/referee/ . Severity: critical (wrong result, leakage, broken code), major (significant error or unsupported claim), minor.`,
    { label: 'referee:both-tracks', phase: 'Unify + Referee', schema: REVIEW_SCHEMA, effort: 'high' }),
])
if (referee) log(`Referee: ${referee.issues.filter(x => x.severity === 'critical').length} critical, ${referee.issues.filter(x => x.severity === 'major').length} major, ${referee.issues.filter(x => x.severity === 'minor').length} minor`)

phase('Final fix')
const finalFix = await agent(`${SHARED}

${FAST} Your own time-box is about 40 minutes.

FINAL FIX PASS over the whole project (models/unified/, ionization.py, docs/unified.md, results/model_comparison.*, and the track code/docs where an issue lives).
Unified-stage summary (JSON): ${JSON.stringify(unified, null, 1)}
Referee report on both tracks (JSON): ${JSON.stringify(referee, null, 1)}
${LI_NOTE}
1. Fix every critical and major issue (code, quick re-run, re-score with evaluate.py, update docs and model_comparison) or rebut it with evidence in a 'Referee responses' section of docs/unified.md. If a track's issue changes its numbers, propagate it to model_comparison and the unified model.
2. Check that the unified model is not affected by any leakage the referee found; re-verify its held-out S1/S2/S3 numbers if the referee questioned them.
3. Smoke-test: py -3.13 ionization.py Fe ; py -3.13 ionization.py 118 ; py -3.13 evaluate.py results/uni_predictions.csv .
4. Add a 'Literature positioning and novelty' section to docs/unified.md using the referee's notes (honest: what is a re-derivation of known physics, what is new here).
5. HEADLINE BLIND TEST (mandatory; this is the paper's key number). Run a truly BLIND S2: fit the unified model's fitted part ONLY on Z <= 54 rows, with NO constraint, bound, prior or term choice that was selected by looking at Z >= 55 results. A constraint is allowed only if justified a priori by physics, e.g. fixing the deep-core and f-class screening to Track A's exact first-order sigma1 values (zexp.coeffs) or using them as ridge priors. Score on Z >= 55 and report MAPE, median APE, p90 and neutral first-IE MAPE, side by side with Track B's blind S2 (178.5% MAPE) and its post-hoc 19.2%. Also report S1 and S3 test numbers for the unified model. Put it in a re-runnable script (models/unified/validate_blind.py) writing results/uni_validation.json, and make this table the first result in docs/unified.md. State plainly whether the derived screening fixes the heavy-element extrapolation. If it does not, say so.
6. COVERAGE CLAIM Z = 1..118 (mandatory). Fix the Z > 110 crash: Track A's QED/nuclear tables end at Z = 110, so extrapolate F_SE and nuclear radii sensibly or fall back gracefully, and document how. Write tests/test_coverage.py, runnable with 'py -3.13 tests/test_coverage.py', that calls the public ionization.py API for EVERY Z = 1..118 and EVERY N = 1..Z, and asserts that each IE is finite and > 0. It also counts monotonicity violations, i.e. cases where a successive IE fails to increase as electrons are removed (IE(Z, N-1) <= IE(Z, N)), lists them, and exits non-zero only on crashes or non-finite/negative values. Run it and report the results in docs/unified.md, then state the claim 'ionization.py covers Z = 1-118, all ion stages' only if the test passes.
7. FRAMING: describe Track B's tau/(Za+kappa) penetration law as an independent data-driven REDISCOVERY of Edlen-type isoelectronic behaviour (cite Edlen), not as a new law. The paper's thesis is 'derived first-order screening (Track A) covers where data-driven fitting cannot extrapolate (Track B), and vice versa'. Support it with the blind S2 numbers.
Return the final unified summary (all fields reflect the final state) plus issues_addressed.`,
  { label: 'final-fix', phase: 'Final fix', schema: FIX_SCHEMA, effort: 'high' })

return { builds, referee, unified, finalFix }
