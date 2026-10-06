export const meta = {
  name: 'shm-competitor-benchmark',
  description: 'Find published screened-hydrogenic model constants (Faussurier 1997, More 1982, Mendoza 2011) from open sources, implement, verify, benchmark vs frozen formula',
  phases: [
    { title: 'Hunt + implement', detail: 'one agent per competitor: locate verbatim constants in open-access sources, implement, validate vs published examples, score on 5847 ions' },
    { title: 'Verify', detail: 'independent adversarial re-transcription and re-implementation check per competitor found' },
    { title: 'Report', detail: 'combined benchmark table with fidelity caveats' },
  ],
}

const SCRATCH = 'C:/Users/HP/AppData/Local/Temp/claude/E--Dr--Navin-Np/43fc5e11-e216-4a38-838f-413efbd43896/scratchpad'

const CONTEXT = `PROJECT CONTEXT
Project root: E:/Dr. Navin/Np (Windows 10; Git Bash + PowerShell). FIRST cd there; relative paths are relative to it. Skim CONTINUE.md for background.
Goal of the project: a closed-form 'Screened Rydberg formula' for successive ionization energies IE(Z,N) of all atoms/ions; the final model is FROZEN: results/pa_hier_rel_predictions.csv (33 fitted params; all-data MAPE 1.87%, blind heavy-element S2 6.6%) and the 9-parameter version results/pa_bound9_predictions.csv (2.90% all-data, blind S2 7.8%). NEVER modify these files, models/push_a/, models/unified/, ionization.py, docs/unified.md, docs/paper_draft.md, docs/literature.md (another workflow may still be finishing them).
Environment: run Python as "py -3.13" (numpy, scipy, pandas, matplotlib Agg, sympy available). No venvs, no pip installs.
Data/scoring (do not modify): data/nist_ie.csv (5847 NIST rows, Z=1..110), common/atomdata.py (load_records(): dicts with Z, N, shells [(n,l,occ)], IE_eV, removed (n,l), rearranged; ground_shells(Z,N); constants HARTREE_EV, RYDBERG_EV, ALPHA), evaluate.py ("py -3.13 evaluate.py results/<name>_predictions.csv --json results/<name>_metrics.json", CSV columns Z,N,IE_pred_eV).
Already benchmarked: Kregar/Di Rocco parameter-free SHM (models/literature/, results/lit_kregar_*_predictions.csv, results/lit_comparison.md: best variant 12.5% MAPE all rows, 2.13% for N<=10, 224% neutral first IE). Its fidelity check against the paper was only partial; learn from that: validate against the paper's own printed numbers.
RULES (strict):
 - Never fabricate or reconstruct constants from memory. Use only values read from a source you actually opened (WebFetch/WebSearch; read PDFs/pages online). Do NOT download or save PDFs or other files to disk; transcribe only the numeric constants you need into JSON with provenance (URL, document, table/equation number, page).
 - A table quoted in a secondary source counts only if it is stated to be the original constants verbatim; cross-check against a second independent source or against numbers printed in the original paper (example energies, total energies, IEs) that your implementation must reproduce.
 - Copyright: do not copy paper text; numeric constants, equations (re-typeset) and short citations only.
 - Report honestly: if constants cannot be obtained or verified, say 'not obtainable' with the list of places you searched. A missing benchmark is acceptable; a fake one is not.
 - IE definition: use the prescription the original paper recommends for ionization energies; if it gives none, use the total-energy difference IE = E(N-1) - E(N) with NIST ground configurations (common.atomdata ground_shells), and also report the one-electron (orbital-energy) variant if cheap. State which is primary BEFORE scoring.`

const COMPETITORS = [
  { key: 'faussurier1997', title: 'Faussurier 1997', brief: `TARGET: G. Faussurier, C. Blancard, A. Decoster, "New screening coefficients for the hydrogenic ion model including l-splitting for fast calculations of atomic structure in plasmas", J. Quant. Spectrosc. Radiat. Transfer 58 (1997) 233 (verify the citation). Screening coefficients sigma(nl, n'l') with l-splitting, fitted to atomic-structure calculations; the strongest fitted competitor (fitted vs fitted). Places to search: OSTI.gov, IAEA INIS (inis.iaea.org: CEA reports often have open full text), HAL / theses.fr (CEA-affiliated theses often reproduce the table), arXiv, Semantic Scholar / CORE / Unpaywall open copies, later open-access papers that tabulate or reuse these coefficients (e.g. relativistic SHM papers, average-atom / opacity code documentation), ResearchGate public previews.` },
  { key: 'more1982', title: 'More 1982', brief: `TARGET: R. M. More, "Electronic energy levels in dense plasmas", J. Quant. Spectrosc. Radiat. Transfer 27 (1982) 345 (verify the citation). Classic screened hydrogenic model with screening constants sigma(n, m) depending on principal quantum numbers only, plus its energy expression. Places to search: OSTI.gov (LLNL preprint/report version, UCRL number), IAEA INIS, arXiv, open-access later papers / theses / code manuals that reproduce More's screening-constant table (e.g. FLYCHK or average-atom documentation, Lanzini/Di Rocco, Faussurier, Mendoza papers comparing to More), Semantic Scholar / CORE / Unpaywall.` },
  { key: 'mendoza2011', title: 'Mendoza 2011', brief: `TARGET: the relativistic screened hydrogenic model of M. A. Mendoza, J. G. Rubiano, J. M. Gil, R. Rodriguez, R. Florido, P. Martel, E. Minguez (likely "A new set of relativistic screening constants for the screened hydrogenic model", High Energy Density Physics 7 (2011) 169; VERIFY the exact citation first). Screening constants with nlj splitting. Places to search: the ULPGC institutional repository (accedacris.ulpgc.es), Universidad de Las Palmas theses, arXiv, OSTI, INIS, Semantic Scholar / CORE / Unpaywall, ResearchGate previews.` },
]

const HUNT_SCHEMA = {
  type: 'object',
  properties: {
    model: { type: 'string' },
    citation_verified: { type: 'string' },
    constants_found: { type: 'boolean' },
    sources: { type: 'array', items: { type: 'object', properties: { url: { type: 'string' }, what: { type: 'string' }, verbatim_original: { type: 'boolean' } }, required: ['url', 'what', 'verbatim_original'] } },
    searched_but_not_found: { type: 'array', items: { type: 'string' } },
    implementation_validation: { type: 'string', description: 'which published numbers were reproduced and the deviations' },
    ie_definition: { type: 'string' },
    deviations_from_paper: { type: 'array', items: { type: 'string' } },
    metrics_markdown: { type: 'string', description: 'MAPE/median/neutral/N<=10/Z>=55/d/f/H-like on 5847 rows, plus S1/S2/S3 test-row MAPE' },
    files: { type: 'array', items: { type: 'string' } },
    confidence: { type: 'string', enum: ['high', 'medium', 'low', 'not-applicable'] },
    notes: { type: 'string' },
  },
  required: ['model', 'citation_verified', 'constants_found', 'sources', 'searched_but_not_found', 'implementation_validation', 'ie_definition', 'deviations_from_paper', 'metrics_markdown', 'files', 'confidence', 'notes'],
}
const VERIFY_SCHEMA = {
  type: 'object',
  properties: {
    model: { type: 'string' },
    transcription_check: { type: 'string' },
    implementation_check: { type: 'string' },
    reproduced_metrics: { type: 'string' },
    issues: { type: 'array', items: { type: 'object', properties: { severity: { type: 'string', enum: ['critical', 'major', 'minor'] }, title: { type: 'string' }, evidence: { type: 'string' } }, required: ['severity', 'title', 'evidence'] } },
    verdict: { type: 'string', enum: ['trustworthy', 'usable-with-caveats', 'not-trustworthy'] },
  },
  required: ['model', 'transcription_check', 'implementation_check', 'reproduced_metrics', 'issues', 'verdict'],
}

phase('Hunt + implement')
const results = await pipeline(COMPETITORS,
  c => agent(`${CONTEXT}

${c.brief}

YOUR AREA: models/benchmarks/${c.key}/ and results files prefixed results/bench_${c.key}_ . Steps:
 1. Verify the exact citation. Search the open web thoroughly for the original constants (tables) and the model's energy/IE equations. Record every place searched.
 2. If found: transcribe the constants into models/benchmarks/${c.key}/constants.json with full provenance per table; implement the model in models/benchmarks/${c.key}/model.py with predict(Z, N, shells=None) -> eV; write models/benchmarks/${c.key}/validate_paper.py reproducing numbers PRINTED in the original paper or the verbatim-quoting source (e.g. example level energies, total energies, IEs, screened charges) and report the deviations.
 3. Score on all 5847 rows: results/bench_${c.key}_predictions.csv (+ variants if the paper allows), evaluate.py metrics JSON; also MAPE on the S1 (Z%5==0), S2 (Z>=55), S3 (N%6==0) test rows (parameter-free / published-constant models are not refit).
 4. models/benchmarks/${c.key}/NOTES.md: sources, equations, choices, deviations, validation, results, confidence.
If the constants cannot be obtained/verified within the time box, stop and return constants_found=false with the search log. Time box: about 30 minutes.`,
    { label: `hunt:${c.key}`, phase: 'Hunt + implement', schema: HUNT_SCHEMA, effort: 'high' }),
  (hunt, c) => {
    if (!hunt || !hunt.constants_found) { log(`${c.title}: constants not obtained; skipping verification`); return { hunt, verify: null } }
    return agent(`${CONTEXT}

INDEPENDENT ADVERSARIAL VERIFICATION of the ${c.title} benchmark implementation. Another agent claims to have found and implemented it:
${JSON.stringify(hunt, null, 1)}
Default to skepticism. Do not edit models/benchmarks/${c.key}/ ; scratch only under ${SCRATCH}/verify_${c.key}/ .
 1. Open the cited sources yourself and RE-TRANSCRIBE at least 25 constants (or all, if the table is small) independently; compare with models/benchmarks/${c.key}/constants.json.
 2. Check the implemented equations against the source (energy expression, screening sum including same-shell factor, any relativistic terms, units, IE definition).
 3. Re-run validate_paper.py and confirm the reproduction of printed values; independently reproduce at least 3 printed numbers with your own short code.
 4. Re-score results/bench_${c.key}_predictions.csv with evaluate.py; spot-check predict() on 10 rows against the CSV.
 5. Check that the comparison is fair: intended domain of the model, configuration/IE definition consistent with the project, no accidental use of NIST IE values.
Verdict: trustworthy / usable-with-caveats / not-trustworthy. Time box: about 20 minutes.`,
      { label: `verify:${c.key}`, phase: 'Verify', schema: VERIFY_SCHEMA, effort: 'high' }).then(v => ({ hunt, verify: v }))
  })

phase('Report')
const report = await agent(`${CONTEXT}

BENCHMARK REPORT. Results per competitor (hunt + independent verification):
${JSON.stringify(results, null, 1)}

Write docs/benchmarks.md and results/bench_comparison.md (+ .csv). You own only these files. Include:
 1. A comparison table on the same 5847 NIST rows, same scorer: Bohr, Slater (total-energy difference), Kregar/Di Rocco variants (results/lit_kregar_*_predictions.csv), every competitor found here (with its verification verdict), and the project's frozen 9-parameter (results/pa_bound9_predictions.csv) and 33-parameter (results/pa_hier_rel_predictions.csv) formulas. Columns: n_params, fitted?, ALL MAPE, median, within-5% fraction, N<=10, 11<=N<=36, N>=37, Z>=55, removed d, removed f, neutral first IE, H-like. Compute every number yourself with evaluate.py.
 2. Fair held-out comparison: for published-constant models the S1/S2/S3 test-row MAPE; for the project's formulas the refit held-out numbers (9p: S1 2.90, S3 3.05, blind S2 7.80; 33p: S1 1.91, S3 2.00, blind S2 6.60; source models/push_a/NOTES.md and results/pa_validation.json; verify them from the JSON).
 3. Fidelity caveats per competitor (sources, verbatim or not, reproduction of printed values, deviations, verdict), intended domain (plasma / highly charged ions), parameter-count asymmetry (fitted vs parameter-free).
 4. A plain verdict paragraph: on which strata the project's formulas beat or lose to each competitor, and whether the 'fitted vs fitted' comparison (Faussurier) could be done.
 5. The list of competitors that could not be obtained, with where we searched (so the user can request them later).
Time box: about 15 minutes. Return a short summary with the headline table (markdown) in your final message.`, { label: 'report', phase: 'Report', effort: 'high' })

return { results, report }
