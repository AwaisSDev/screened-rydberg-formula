export const meta = {
  name: 'ionization-final-push',
  description: 'Blind-extrapolation push (2 approaches) -> integrate || literature -> audit || paper draft -> final fix',
  phases: [
    { title: 'Push', detail: 'two independent approaches to get blind heavy-element (S2) error under ~10% with a pre-registered selection rule' },
    { title: 'Integrate + Literature', detail: 'winner becomes final model, tests, docs || prior-art research and head-to-head vs published screened-hydrogenic models' },
    { title: 'Audit + Paper', detail: 'adversarial audit of final numbers and protocol || paper-style manuscript draft' },
    { title: 'Final fix', detail: 'apply audit findings everywhere, run tests, update CONTINUE.md' },
  ],
}

const SCRATCH = 'C:/Users/HP/AppData/Local/Temp/claude/E--Dr--Navin-Np/43fc5e11-e216-4a38-838f-413efbd43896/scratchpad'

const CONTEXT = `PROJECT CONTEXT (read carefully)
Project root: E:/Dr. Navin/Np (Windows 10; Git Bash + PowerShell). FIRST cd to the project root; all relative paths are relative to it. Read CONTINUE.md and CLAUDE.md there first: they describe the goal, layout, results so far and rules.
Goal: a general closed-form ("Screened Rydberg formula") for successive ionization energies IE(Z,N) of every atom/ion from Z and the electron configuration, as accurate and as first-principles as possible, validated honestly. The user will write a paper from this.
Environment: ALWAYS run Python as "py -3.13" (numpy, scipy, pandas, matplotlib(Agg), sympy, scikit-learn, numba available). No venvs, no pip installs. 4 CPU cores shared with one other agent: at most 2 processes in parallel; no single job longer than ~10 minutes.
Shared infrastructure (do NOT modify): data/nist_ie.csv (5847 NIST rows, Z=1..110), common/atomdata.py (load_records, get, ground_shells, madelung_shells, constants), evaluate.py ("py -3.13 evaluate.py results/<name>_predictions.csv [--json ...]"; in Python evaluate.evaluate(csv, subset=lambda rec: ...)).
Existing work: models/first_principles/ (Track A: exact first-order screening sigma1 via zexp.coeffs / models/unified/abinitio.py, Dirac/recoil/finite-size/QED in relativity.py, own Kohn-Sham LSDA solver ks_atom.py), models/semi_empirical/ (Track B: GSHM gshm.py, fitlib.py), models/unified/ (umodel.py = exact sigma1 + fitted remainder "u35"/"u29"; pocket.py = 8-parameter pocket formula; final.py dispatcher; run_unified.py; report.py), ionization.py (top-level API/CLI), results/ (predictions, metrics, model_comparison.md, uni_validation.json), docs/ (first_principles.md, semi_empirical.md, unified.md [currently STALE, has a banner]), handoff/ (builds.json, referee_report.json, unify_result.json).
Current numbers (MAPE): u35 (35 params) all-data 1.48%, neutral first IE 6.43%, H-like 0.0012%, S1 1.52%, S3 1.65%, blind S2 57.2% (median 2.0%, neutral 3831% = catastrophic blow-ups for heavy near-neutral atoms). pocket (8 params) all-data 4.68%, S1 4.62%, S3 5.24%, blind S2 11.3% (neutral 12.3%: robust). Track B GSHM blind S2 170%. Worst u35 near-neutral heavy rows: Po, Hs, No, Ac, Tl, Ta, Sg, Hf, Md (6p/6d/7s/5d removal).
Honesty rules: report real numbers; count every fitted parameter; no per-element/per-ion parameters; predict() must never read NIST IE values; no test-set leakage; disclose any design decision influenced by results you have seen. The referee found earlier that a constraint chosen after looking at S2 inflated a headline number; never do that again.`

const PROTOCOL = `PRE-REGISTERED VALIDATION PROTOCOL (mandatory, identical for every candidate model):
 - V1 (inner extrapolation into new shells): fit on Z <= 36, validate on 37 <= Z <= 54.
 - V2 (inner extrapolation): fit on Z <= 44, validate on 45 <= Z <= 54.
 - S1 (interpolation): fit on rows with Z % 5 != 0, test on Z % 5 == 0.
 - S3 (unseen isoelectronic sequences): fit on N % 6 != 0, test on N % 6 == 0.
 - SELECTION SCORE = mean of the four MAPEs (V1, V2, S1, S3). All design choices, hyper-parameters, term selection and the choice of your final candidate must use ONLY these (and data with Z <= 54 for anything extrapolation-related). Rows with Z >= 55 must never influence any choice.
 - S2 (blind extrapolation): fit on Z <= 54, test on Z >= 55. Compute it ONCE per candidate, AFTER all choices are frozen, purely for reporting. Never iterate on it.
 - Also report the all-data fit (fit on all 5847 rows) scored with evaluate.py: MAPE, median, neutral first-IE MAPE, H-like MAPE.
 - For every split report MAPE, median APE and neutral first-IE MAPE of the test/validation rows.
Include the existing u35 and pocket models as reference candidates evaluated with exactly this protocol (refit them per split using their existing fit functions).`

const PUSH_SCHEMA = {
  type: 'object',
  properties: {
    approach: { type: 'string' },
    summary: { type: 'string' },
    formula: { type: 'string', description: 'Final chosen formula in full, every symbol defined' },
    candidates: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          name: { type: 'string' },
          n_params: { type: 'integer' },
          V1_MAPE: { type: 'number' }, V2_MAPE: { type: 'number' }, S1_MAPE: { type: 'number' }, S3_MAPE: { type: 'number' },
          selection_score: { type: 'number' },
          all_MAPE: { type: 'number' }, all_median: { type: 'number' }, all_neutral_MAPE: { type: 'number' }, all_Hlike_MAPE: { type: 'number' },
          S2_MAPE: { type: 'number' }, S2_median: { type: 'number' }, S2_neutral_MAPE: { type: 'number' },
          predictions_csv: { type: 'string', description: 'all-data-fit predictions, all 5847 rows' },
          module: { type: 'string', description: 'python module + function to call: predict(Z, N, shells=None)' },
        },
        required: ['name', 'n_params', 'V1_MAPE', 'V2_MAPE', 'S1_MAPE', 'S3_MAPE', 'selection_score', 'all_MAPE', 'all_median', 'all_neutral_MAPE', 'all_Hlike_MAPE', 'S2_MAPE', 'S2_median', 'S2_neutral_MAPE', 'predictions_csv', 'module'],
      },
    },
    chosen: { type: 'string', description: 'name of your candidate with the lowest selection_score (never chosen by S2)' },
    disclosures: { type: 'array', items: { type: 'string' } },
    files: { type: 'array', items: { type: 'string' } },
  },
  required: ['approach', 'summary', 'formula', 'candidates', 'chosen', 'disclosures', 'files'],
}

const PUSHES = [
  { key: 'A', dir: 'models/push_a/', prefix: 'pa_', title: 'ROBUST-MINIMAL', brief: `APPROACH A: ROBUST / MINIMAL REMAINDER. Start from the exact first-order screening sigma1 (Track A; models/unified/abinitio.py) plus the 8-parameter pocket formula's robust structure, and add only terms that improve the SELECTION SCORE. Ideas (use physics, not test peeking):
 - A-priori physical bounds on the effective charge: Z_a <= Zeff <= Z (Z_a = Z - N + 1). The lower bound is the non-penetrating hydrogenic limit: penetration (quantum defect >= 0) can only increase binding, so IE >= Ry Z_a^2/n^2. Implement as smooth saturating transforms, not hard clips that kill gradients.
 - Bounded, saturating forms for the charge-dependent penetration/relaxation term (e.g. tau*nu/(Za+kappa) replaced by a form that cannot exceed a physically sensible fraction of the screening electrons), monotonic in the number of screening electrons.
 - Exact sigma1 everywhere it is available as the large-Z limit (sigma -> sigma1 as Z -> infinity at fixed N), so the fitted part only describes O(N/Z) corrections.
 - Few parameters; tie f-class parameters to d analogues when unsupported; avoid (Z alpha)^2 exponential correction terms that extrapolate badly unless they win on V1/V2.
 - Parameter-free Dirac factor + Track A recoil/QED/finite-size layer (so H-like stays ~0.001%).
Goal: blind S2 MAPE under ~10% while keeping all-data and S1/S3 MAPE as low as possible (ideally <= 2.5%).` },
  { key: 'B', dir: 'models/push_b/', prefix: 'pb_', title: 'PHYSICS-CONSTRAINED FULL', brief: `APPROACH B: PHYSICS-CONSTRAINED FULL MODEL. Start from the u35 structure (models/unified/umodel.py: exact sigma1 + fitted remainder) which is the most accurate in-distribution, and remove its extrapolation failure modes using physics constraints decided a priori:
 - Diagnose (on Z <= 54 data and V1/V2 only, plus physical reasoning about what happens as N/Z, n, l grow) which terms can run away: e.g. tau_g*nu_g/(Za+kappa) with large electron counts for heavy neutral atoms (Za = 1) driving Zeff to the clip floor or far above Z; exp((Z alpha)^2) correction terms; unsupported f-class parameters.
 - Constrain: physical bounds Z_a <= Zeff <= Z via smooth transforms (non-penetrating hydrogenic limit: IE >= Ry Za^2/n^2); sign/magnitude bounds on tau from screening physics (a screening electron screens at most 1 charge unit and at least 0); per-electron (not per-group-total) penetration so the term saturates with electron count; f-class tied to d analogues or fixed by Track A exact sigma1 values; ridge toward the pure ab-initio sigma1 model with strength chosen on V1/V2.
 - Optionally use Track A's LSDA Delta-SCF code (models/first_principles/ks_atom.py, run_dft.py) only as a physics prior for Z <= 54 (never fit to or look at Z >= 55 DFT or NIST values for choices).
Goal: blind S2 MAPE under ~10% while keeping all-data MAPE near u35's 1.48% and neutral first-IE MAPE near 6.4%.` },
]

phase('Push')
log('Two independent approaches in parallel; selection by pre-registered inner-validation score (no Z>=55 information)')
const pushes = await parallel(PUSHES.map(p => () => agent(`${CONTEXT}

${PROTOCOL}

${p.brief}

YOUR AREA: write code only in ${p.dir} and results files prefixed results/${p.prefix} (figures results/figures/${p.prefix}*.png). Import (never modify) models/unified, models/first_principles, models/semi_empirical, common, evaluate.py. Another agent is concurrently working on a different approach in a different directory; do not touch its files.
DELIVERABLES: (1) ${p.dir}model.py with predict(Z, N, shells=None) -> eV (works for any Z up to 120, any N <= Z; params loaded from results/${p.prefix}params.json), predict_many, fit(mask) ; (2) ${p.dir}validate.py that re-runs the full protocol (V1, V2, S1, S3, then S2 once, plus all-data fit) for every candidate including u35 and pocket references and writes results/${p.prefix}validation.json; (3) results/${p.prefix}<candidate>_predictions.csv (all-data fit, all 5847 rows) for each candidate you report; (4) ${p.dir}NOTES.md: what you tried, what each change bought on the selection score, final formula, parameter table, disclosures.
Time box: about 35 minutes. Return the structured result (all candidates with complete metrics; chosen = lowest selection_score among YOUR candidates; disclosures must be candid, including any idea that was inspired by having already seen which heavy atoms failed).`,
  { label: `push:${p.title}`, phase: 'Push', schema: PUSH_SCHEMA, effort: 'high' })))

const cands = []
pushes.forEach((r, i) => {
  if (!r) { log(`Push ${PUSHES[i].key} failed`); return }
  r.candidates.forEach(c => cands.push({ ...c, push: PUSHES[i].key }))
})
const complete = cands.filter(c => [c.V1_MAPE, c.V2_MAPE, c.S1_MAPE, c.S3_MAPE].every(x => typeof x === 'number' && isFinite(x)))
complete.forEach(c => { c.score_recomputed = (c.V1_MAPE + c.V2_MAPE + c.S1_MAPE + c.S3_MAPE) / 4 })
complete.sort((a, b) => a.score_recomputed - b.score_recomputed)
const winner = complete[0] || null
log(winner ? `Pre-registered winner: ${winner.name} (push ${winner.push}) score ${winner.score_recomputed.toFixed(3)}; blind S2 ${winner.S2_MAPE.toFixed(2)}%` : 'No complete candidate; keeping u35')
const ranking = complete.map(c => ({ name: c.name, push: c.push, n_params: c.n_params, score: +c.score_recomputed.toFixed(4), V1: c.V1_MAPE, V2: c.V2_MAPE, S1: c.S1_MAPE, S3: c.S3_MAPE, S2_blind: c.S2_MAPE, S2_neutral: c.S2_neutral_MAPE, all: c.all_MAPE, neutral: c.all_neutral_MAPE, module: c.module, predictions_csv: c.predictions_csv }))

const INTEGRATE_SCHEMA = {
  type: 'object',
  properties: {
    summary: { type: 'string' },
    final_model: { type: 'string' },
    final_formula: { type: 'string' },
    headline_table_markdown: { type: 'string', description: 'final + references: params, all-data MAPE, neutral, H-like, V1, V2, S1, S3, blind S2' },
    coverage_test_result: { type: 'string' },
    files: { type: 'array', items: { type: 'string' } },
    open_issues: { type: 'array', items: { type: 'string' } },
  },
  required: ['summary', 'final_model', 'final_formula', 'headline_table_markdown', 'coverage_test_result', 'files', 'open_issues'],
}
const LIT_SCHEMA = {
  type: 'object',
  properties: {
    summary: { type: 'string' },
    novelty_assessment: { type: 'string' },
    head_to_head_markdown: { type: 'string', description: 'published models implemented and scored on the same 5847 rows, or why not possible' },
    venues: { type: 'array', items: { type: 'string' } },
    files: { type: 'array', items: { type: 'string' } },
  },
  required: ['summary', 'novelty_assessment', 'head_to_head_markdown', 'venues', 'files'],
}

phase('Integrate + Literature')
const [integ, lit] = await parallel([
  () => agent(`${CONTEXT}

${PROTOCOL}

INTEGRATION STAGE. Two push agents produced candidates. The PRE-REGISTERED winner (lowest selection score = mean of V1, V2, S1, S3; chosen by the orchestrator's script, never by S2) is:
${JSON.stringify(winner, null, 1)}
Full ranking of all candidates: ${JSON.stringify(ranking, null, 1)}
Push A report: ${JSON.stringify(pushes[0], null, 1)}
Push B report: ${JSON.stringify(pushes[1], null, 1)}
(If the winner is null, keep u35 as final.)

Make the winner the project's FINAL model. You own models/unified/, ionization.py, tests/, results/uni_* and results/model_comparison.*, docs/unified.md. Do not modify models/push_a, models/push_b, models/first_principles, models/semi_empirical (import only). Do NOT write docs/literature.md or models/literature/ (another agent is writing them concurrently).
 1. models/unified/final.py must dispatch to the winner by default (keep u35, u29, pocket, gshm_qed available as options); ionization.py uses it by default and keeps the CLI: "py -3.13 ionization.py Fe", "py -3.13 ionization.py 26 --N 26", "py -3.13 ionization.py 118", --model option. Regenerate results/uni_predictions.csv (winner, all-data fit, all 5847 rows) and score it with evaluate.py.
 2. models/unified/validate_blind.py: one re-runnable script that applies the full protocol (V1, V2, S1, S3, blind S2, all-data) to the final model and the references (u35, pocket, Track B GSHM, Slater total-energy baseline) and writes results/uni_validation.json + results/uni_validation.md. Its numbers must match the push agents' numbers for the same models; if they differ, investigate and report.
 3. tests/test_coverage.py ("py -3.13 tests/test_coverage.py"): calls the public ionization.py API for every Z = 1..118 and every N = 1..Z; asserts finite and > 0; counts and lists monotonicity violations (IE(Z, N-1) <= IE(Z, N), i.e. a successive IE not increasing as electrons are removed); exits non-zero only for crashes or non-finite/non-positive values. Run it; report the result.
 4. results/model_comparison.csv/.md: all models (baselines, Track A, Track B, previous unified variants, all push candidates, final) with the same columns incl. V1, V2, S1, S3, blind S2, n_params, type (first-principles / fitted / hybrid). Regenerate the uni_* figures for the final model (parity, first IE vs Z, successive IEs for C, Fe, Xe, U, residuals vs Z).
 5. Rewrite docs/unified.md completely (remove the STALE banner): (a) headline validation table FIRST (blind S2 next to Track B's 170-178% and the earlier post-hoc 19%), with a plain statement whether derived screening fixed heavy-element extrapolation; (b) the final equation in full with a parameter table and derivation lineage (which term came from Track A exact theory, which from Track B, which from the push); (c) the selection-protocol history, candidly: the earlier unified model was selected using S2 (a mild leak), now replaced by the pre-registered V1/V2/S1/S3 rule; disclosures from the push agents; (d) a 'Referee responses' section answering every issue in handoff/referee_report.json; (e) framing: Track B's tau/(Za+kappa) law is an independent rediscovery of Edlen-type isoelectronic behaviour, not a new law; the thesis is 'derived first-order screening covers where data fitting cannot extrapolate, and vice versa'; (f) a worked hand example (first IE of O, and third IE of Mg) with the pocket formula and the final formula; (g) limitations; (h) leave a section heading 'Literature positioning and novelty' containing only: 'See docs/literature.md' (the literature agent fills that file).
Time box: about 30 minutes.`, { label: 'integrate', phase: 'Integrate + Literature', schema: INTEGRATE_SCHEMA, effort: 'high' }),
  () => agent(`${CONTEXT}

LITERATURE + HEAD-TO-HEAD STAGE. You own ONLY docs/literature.md, docs/references.bib, models/literature/ and results files prefixed results/lit_ . Another agent is concurrently editing models/unified, ionization.py, docs/unified.md and results/uni_*: do not touch those.
 1. Research prior art with WebSearch/WebFetch (cite precisely; never invent references; mark anything you could not verify): Slater 1930 screening rules; Clementi & Raimondi 1963 / Clementi, Raimondi & Reinhardt 1967; Layzer 1959 and later 1/Z expansion work (e.g. Dalgarno, Stewart; Safronova et al. Z-expansion for isoelectronic sequences); Edlen's isoelectronic-sequence formulas (Handbuch der Physik 1964 and related) and the screening-parameter/quantum-defect regularities along sequences; screened hydrogenic models (Mayer 1947, More 1982 JQSRT, Faussurier, Blancard, Renaudin, Lanzini et al., relativistic SHM e.g. Faussurier 1997, Mendoza/Rubiano; average-atom codes); Hartree-Fock Koopmans and Delta-SCF accuracy for IEs; LDA/GGA Delta-SCF IE benchmarks; Rodrigues et al. 2004 (ADNDT) Dirac-Fock IEs for all ions; NIST ASD (Kramida et al.) and Yerokhin & Shabaev 2015 for H-like; Kotochigova et al. 1997 NIST LDA reference; Chakravorty & Davidson 1993 exact non-relativistic energies; recent machine-learning / symbolic-regression approaches to ionization energies.
 2. HEAD-TO-HEAD (the key 'is it publishable' test): if you can obtain a published screened-hydrogenic model's screening constants in full (e.g. More 1982 screening-constant formula, or Faussurier et al. 1997 relativistic SHM tables), implement it in models/literature/ with predict(Z, N), write results/lit_<model>_predictions.csv for all 5847 rows, score with evaluate.py, and compare with the project's models (results/model_comparison.md has current numbers; the final model may still be changing, so compare against u35, pocket and Track B GSHM, and note the final model will be added later). If the constants cannot be obtained reliably, say so explicitly and do not fabricate them.
 3. docs/literature.md: prior-art survey; table of published accuracies where available; an honest novelty assessment (what in this project is a re-derivation of known physics: 1/Z expansion, Edlen-type behaviour, screened hydrogenic form; what is plausibly new: e.g. an exact first-order screening table for all ground configurations Z=1..110 as an ab-initio replacement of Slater's rules, a single validated hybrid formula over all 5847 ions with blind extrapolation tests, specific functional forms); recommended venues with reasons; required disclosures (AI assistance). docs/references.bib with every cited work.
Time box: about 30 minutes.`, { label: 'literature', phase: 'Integrate + Literature', schema: LIT_SCHEMA, effort: 'high' }),
])

const AUDIT_SCHEMA = {
  type: 'object',
  properties: {
    verified_claims: { type: 'array', items: { type: 'string' } },
    issues: { type: 'array', items: { type: 'object', properties: { severity: { type: 'string', enum: ['critical', 'major', 'minor'] }, title: { type: 'string' }, evidence: { type: 'string' }, suggested_fix: { type: 'string' } }, required: ['severity', 'title', 'evidence', 'suggested_fix'] } },
    reproduced_headline: { type: 'string' },
    overall_verdict: { type: 'string' },
  },
  required: ['verified_claims', 'issues', 'reproduced_headline', 'overall_verdict'],
}
const PAPER_SCHEMA = {
  type: 'object',
  properties: { summary: { type: 'string' }, files: { type: 'array', items: { type: 'string' } }, numbers_used: { type: 'string', description: 'list of every number in the draft and the results file it came from' } },
  required: ['summary', 'files', 'numbers_used'],
}

phase('Audit + Paper')
const [audit, paper] = await parallel([
  () => agent(`${CONTEXT}

${PROTOCOL}

FINAL ADVERSARIAL AUDIT. Default to skepticism; a claim is unverified until you reproduce it. Do not edit project files; scratch only under ${SCRATCH}/final_audit/ .
Integration summary: ${JSON.stringify(integ, null, 1)}
Literature summary: ${JSON.stringify(lit, null, 1)}
Pre-registered ranking: ${JSON.stringify(ranking, null, 1)}
Check:
 1. Re-run models/unified/validate_blind.py (or your own independent implementation of the protocol) and confirm the final model's V1, V2, S1, S3, blind S2 and all-data numbers; re-score results/uni_predictions.csv with evaluate.py.
 2. Protocol integrity: read the push and unified code; confirm that no fitted value, hyper-parameter, bound, term choice or candidate selection used rows with Z >= 55 (other than the single final S2 report). Confirm the winner really has the lowest pre-registered selection score. Confirm predict() never reads NIST IE values; no per-element/per-ion parameters; parameter counts are right.
 3. Run tests/test_coverage.py; run ionization.py for H, He, C, O, Na, Fe, Ni, Cu, Xe, Gd, W, Au, U, 118, and check sanity (signs, magnitudes, monotonic successive IEs, H-like values).
 4. Doc consistency: every number in docs/unified.md headline table, results/model_comparison.md and docs/literature.md head-to-head must match re-computed values (spot-check >= 20 numbers). Check framing honesty (no 'new law', no overclaiming 'first principles' for fitted parts, H-like QED consistency caveat, DFT validation framed as implementation check).
 5. Literature head-to-head: if a published model was implemented, spot-check its implementation against the cited source description.
Severity: critical (wrong numbers, leakage, broken code), major (unsupported claim, protocol violation), minor.`, { label: 'final-audit', phase: 'Audit + Paper', schema: AUDIT_SCHEMA, effort: 'high' }),
  () => agent(`${CONTEXT}

PAPER DRAFT STAGE. Write docs/paper_draft.md: a journal-style manuscript draft titled along the lines of 'A screened Rydberg formula for the successive ionization energies of all atoms and ions'. You own ONLY docs/paper_draft.md and results/figures/paper_*.png. Read docs/unified.md, docs/first_principles.md, docs/semi_empirical.md, docs/literature.md, docs/references.bib, results/model_comparison.md, results/uni_validation.json/.md.
Integration summary: ${JSON.stringify(integ, null, 1)}
Literature summary: ${JSON.stringify(lit, null, 1)}
Structure: Abstract (with the honest headline numbers incl. blind extrapolation); 1 Introduction (Bohr, Slater, Clementi-Raimondi, Layzer, Edlen, screened-hydrogenic models, DFT; the gap: a single validated closed form over all ions with blind tests); 2 Theory (exact first-order screening from the 1/Z expansion with the He 5/8 and Li 5965/5832 examples; the screened Rydberg form; charge-dependent penetration as an Edlen-type term; relativistic Dirac factor, recoil, finite size, QED; DFT used as a physics check); 3 Data and validation protocol (NIST ASD 5847 rows, status categories, pre-registered V1/V2/S1/S3 selection, single blind S2); 4 Results (tables: baselines vs final; per-stratum errors; blind extrapolation; H-like layers with the Yerokhin-Shabaev consistency caveat; head-to-head with published models if available; the exact sigma1 table excerpt vs Slater); 5 Discussion (what is new vs rediscovered; where it fails: heavy near-neutral atoms, multiplet structure, rearranged configurations; why no exact closed form exists for N >= 2); 6 Conclusions; Data and code availability; AI-assistance disclosure statement (the work was carried out with AI agents under the author's direction; the author must verify); References (from docs/references.bib only: never invent references). Use figures that exist in results/figures/ (reference by path) and make 1-2 paper_*.png figures if useful.
Every number must come from a results file or doc: list each number and its source in numbers_used. Where the final numbers may still change after the audit, keep them exactly as in the current results files; the final fix pass will update them. Time box: about 25 minutes.`, { label: 'paper-draft', phase: 'Audit + Paper', schema: PAPER_SCHEMA, effort: 'high' }),
])
if (audit) log(`Audit: ${audit.issues.filter(x => x.severity === 'critical').length} critical, ${audit.issues.filter(x => x.severity === 'major').length} major, ${audit.issues.filter(x => x.severity === 'minor').length} minor`)

phase('Final fix')
const FINAL_SCHEMA = {
  type: 'object',
  properties: {
    summary: { type: 'string' },
    final_headline_markdown: { type: 'string', description: 'final verified headline table' },
    issues_addressed: { type: 'array', items: { type: 'string' } },
    remaining_limitations: { type: 'array', items: { type: 'string' } },
    files_changed: { type: 'array', items: { type: 'string' } },
  },
  required: ['summary', 'final_headline_markdown', 'issues_addressed', 'remaining_limitations', 'files_changed'],
}
const finalFix = await agent(`${CONTEXT}

${PROTOCOL}

FINAL FIX PASS. Apply the audit findings across the whole project (code, results, docs/unified.md, results/model_comparison.*, docs/literature.md, docs/paper_draft.md). Fix every critical and major issue (or rebut with evidence in an 'Audit responses' section of docs/unified.md); fix minor ones where cheap. Never change the pre-registered selection or use Z >= 55 rows for any choice.
Audit report: ${JSON.stringify(audit, null, 1)}
Integration summary: ${JSON.stringify(integ, null, 1)}
Paper draft summary: ${JSON.stringify(paper, null, 1)}
Then: re-run models/unified/validate_blind.py, py -3.13 evaluate.py results/uni_predictions.csv, py -3.13 tests/test_coverage.py, py -3.13 ionization.py Fe, py -3.13 ionization.py 118; make sure docs/unified.md, results/model_comparison.md and docs/paper_draft.md all carry the same final numbers. Finally update CONTINUE.md: set 'Last updated' to the time this run finished (read it with 'date'), mark the completed checklist items, record the final headline numbers and the remaining open items (e.g. venue choice, author verification, optional further work). Time box: about 25 minutes.`, { label: 'final-fix', phase: 'Final fix', schema: FINAL_SCHEMA, effort: 'high' })

return { ranking, winner, pushes, integ, lit, audit, paper, finalFix }
