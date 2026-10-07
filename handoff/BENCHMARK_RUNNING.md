# Competitor benchmark: a second workflow running in parallel

Started **2026-10-05 02:56 PKT**. Workflow `wf_1a66b31c-6b0`; script saved as `handoff/benchmark_workflow.js`.

The user asked to benchmark the frozen Screened Rydberg formula against the published screened-hydrogenic models
(SHMs). Priority: Faussurier et al. 1997 (JQSRT 58, 233), then More 1982 (JQSRT 27, 345), then Mendoza et al. 2011
(HEDP). The user could not get the paywalled PDFs, so the agents search open-access sources: OSTI, IAEA INIS, HAL,
theses, and papers that quote the tables verbatim.

Rules the agents follow:
- Never reconstruct constants from memory.
- Never download files to disk.
- Prove each implementation by reproducing numbers printed in the paper.
- Report "not obtainable" if a source can't be found.

Stages:
1. **Hunt and implement**, one agent per model. Output goes to `models/benchmarks/<model>/` and `results/bench_<model>_*`.
2. **Independent verification** of each model that was found.
3. **Report:** `docs/benchmarks.md` and `results/bench_comparison.{md,csv}`.

The workflow does not touch the frozen model or `docs/unified.md`, `docs/paper_draft.md` or `docs/literature.md`.
After it finishes, merge the benchmark table into the paper draft and the literature doc, and record it in
`CONTINUE.md`.
