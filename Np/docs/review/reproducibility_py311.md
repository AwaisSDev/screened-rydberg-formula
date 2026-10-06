# Cross-environment reproducibility check (2026-10-06)

**Why.** The project moved from `E:\Dr. Navin\Np` (Windows 10, Python 3.13 with numba) to `D:\Chem-Research\Np` (Windows 11).
The new machine has no Python 3.13. Python 3.11.9 has NumPy 2.4.4 and SciPy 1.17.1 but no numba, so the runs here use
`py -3.11` with the no-op stand-in `tools/numba_stub` on `PYTHONPATH`. The numba-decorated functions are only the FNS and
DFT solvers; their outputs are cached in `models/first_principles/cache/`.

**What was run.** `models/unified/validate_blind.py` ran on a scratch copy of the project, so no project file was
overwritten. It refits every model on every split (V1, V2, S1, S3, blind S2, all data). Run time was 20 s.

## Frozen artefacts are intact

- `results/pa_hier_rel_predictions.csv` and `results/pa_params.json` match the sha256 in `handoff/FROZEN.txt`.
- `models/push_a/model.py` matches as well once git's CRLF checkout is normalised back to LF.
- `evaluate.py results/uni_predictions.csv` gives 1.874 % overall, 0.940 % median, 7.515 % neutral, 0.001 % H-like
  and 4.631 % on the 311 experimental rows, the same as before.

## Refit results: frozen 3.13 run vs this 3.11 run

| model / split | frozen (3.13) | this machine (3.11) | comment |
|---|---|---|---|
| final pa_hier_rel: selection score | 2.1481 | 2.1473 | robust |
| final: V1 / V2 / S1 / S3 | 3.085 / 1.603 / 1.909 / 1.995 | 3.081 / 1.603 / 1.909 / 1.995 | robust |
| final: **blind S2** (median / neutral) | **6.6031** (1.616 / 21.76) | **6.6030** (1.617 / 21.74) | robust |
| final: all-data | 1.8742 | 1.8742 | robust |
| u35: blind S2 | 57.23 | 56.50 | optimizer-path dependent |
| u35: V1 | 26.39 | 26.34 | |
| u29: V1 | 21.97 | **55 914 (diverges)** | optimizer-path dependent |
| **pocket: V1** | **4.4·10⁶ (diverges)** | **4.35** | the paper's "fit diverges" is not a property of the model |
| pocket: selection score | 1.1·10⁶ | 4.19 | |
| pocket: blind S2 | 11.334 | 11.334 | robust |
| GSHM final: blind S2 | 170.27 | 170.24 | robust |
| Slater total: S2 | 11.638 | 11.638 | no fit |

## Consequences for the paper

1. Every headline number of the final model reproduces to within 0.005 percentage points on a different NumPy/SciPy
   stack. This supports the claim that the bounded model is numerically stable.
2. Two printed reference entries are artefacts of the original optimizer path, not properties of the models:
   - the pocket formula's V1 entry "4.4·10⁶ (fit diverges)", which converges to 4.35 % here;
   - u29's V1.

   The paper should either report both outcomes ("diverges or 4.35 % depending on the least-squares path") or drop the
   pocket formula's V1 and selection-score entries. It must not use the V1 divergence as evidence against the pocket
   formula. Figure 2 omits the pocket formula "because its V1 fit diverges", and that caption needs the same caveat.
3. Add a computational-details sentence naming the software versions of the frozen run. They are not recorded
   anywhere in the repo; the author must take them from the original machine (Python 3.13.x, NumPy, SciPy versions).
   Also note that the 3.11 re-run reproduces the final model's numbers.
