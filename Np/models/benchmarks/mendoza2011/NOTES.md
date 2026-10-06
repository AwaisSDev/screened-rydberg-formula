# Benchmark: Mendoza et al. 2011 relativistic screened hydrogenic model (NRSHM)

## Citation (verified)
M.A. Mendoza, J.G. Rubiano, J.M. Gil, R. Rodríguez, R. Florido, P. Martel, E. Mínguez,
"A new set of relativistic screening constants for the screened hydrogenic model",
*High Energy Density Physics* **7**(3) (2011) 169–179, doi:10.1016/j.hedp.2011.04.006.
Checked against three sources: the Crossref record, the UPM repository record (https://oa.upm.es/11165/), and
the title page of the deposited article.

## Source of the constants
- **The original article.** It is an author-deposited scan with an OCR text layer, held as green open access in the
  Universidad Politécnica de Madrid repository: https://oa.upm.es/11165/2/INVE_MEM_2011_102088.pdf. The link was
  found through the Semantic Scholar `openAccessPdf` field.
- **Constants.** Table 1 (journal p. 171) and Table 2 (journal p. 172) give the full 19×19 matrix σ_kk' for k,
  k' = 1s1/2 … 5p3/2. It was transcribed from the OCR layer, and every row was checked by eye against the page
  image. Each row has the expected 9 and 10 columns. The constants are in `constants.json`, with provenance.
- **Not used as a source.** The ULPGC repository (accedacris.ulpgc.es/handle/10553/45566 and 10553/45552)
  returned HTTP 403 to every automated fetch.
- **File handling.** The PDF was read online through WebFetch. WebFetch keeps its own temporary cache of the binary.
  Nothing was saved into the project, and the cache copy was removed after transcription.

## Equations (paper numbering, atomic units, c = 1/α)
- (1) E_T = Σ_k P_k ε_k
- (2) ε_k = c²{[1 + (αQ_k / (n − j − ½ + √((j+½)² − (αQ_k)²)))²]^(−½) − 1}, which is the Dirac eigenvalue.
- (3) Q_k = Z − Σ_k' σ_kk' (P_k' − δ_kk'). The row is the screened subshell k and the column is the screening
  subshell k'.
- App. A: the Xα binding energy is ε_k^Xα = ∂E_T/∂P_k. The model computes it by numerical differentiation of
  eq. (1).
- App. B: ε^HF = ε^Xα − σ_kk/r_k, with r_k = n̄√(n̄² + α²Q²)/Q and n̄ = n − |κ| + √(κ² − α²Q²).

## IE prescription (chosen before scoring)
**Primary: IE = E_T(N−1) − E_T(N)**, using NIST ground configurations from `common.atomdata.ground_shells`.
- The paper fits and tabulates ionization energies this way, and Table 3 is reproduced exactly with it.
- Sec. 5.2 of the paper recommends total-energy differences for isolated atoms.

**Variants:**
- `_xalpha`: the one-electron binding energy, −ε^Xα, of the subshell that loses the electron.
- `_stat`: statistical j-split.

**j-splitting.** The primary rule fills j = l−½ before j = l+½. This gives the jj configurations the paper prints in
Table 4, and it is the rule that reproduces Tables 3–8. The `stat` variant splits the occupations as
q(2j+1)/(2(2l+1)).

**Fe-like row of Table 3.** It is reproduced only with the NIST ground configurations (for example Co II 3d⁸).
Madelung filling gives deviations of up to 270 %. So the paper also used NIST ground configurations.

## Validation against numbers printed in the paper (`validate_paper.py` → `results/bench_mendoza2011_validation.json`)
| Printed quantity | Rows | Max deviation |
|---|---|---|
| Table 3: IE (a.u.) of the C-, Ne-, Al-, Ar- and Fe-like sequences, Z = 6–32 | 84 | 0.023 % (rounding) |
| Table 4: N I–VI and O I–VII configuration-average energies | 13 | 0.007 % |
| Table 5: Na-like total energies, Z = 18–79 | 13 | 0.006 % |
| Table 6: Ne-like binding energies, Xα and HF columns, Z = 16, 26, 66, 91 | 16 | 0.05 % and 0.02 % |
| Table 7: Au²⁵⁺ (4f⁸) level energies, 1s–4f7/2 (NRSHM column = Xα) | 16 | 0.16 % (rounding to 0.1 a.u.) |
| Table 8: Fe²⁺, Fe⁶⁺, Fe⁹⁺ and Fe¹²⁺ binding energies (eV), Xα and HF columns | 24 | 0.006 % and 0.005 % |

**Not reproduced.** Table 6's "der 1/2" column (Slater transition state) does not match. Our transition-state
reading deviates by up to 13 %, and the paper does not spell out its definition. That column is not used for any
IE.

**Not covered by any printed check.** No printed value exercises the 5s/5p rows of the matrix, or the 5s/5p columns
that only matter when those subshells are occupied. They were checked only by eye against the scan.

## Results on the 5847 NIST rows (no refit)
**Coverage: 5011 of 5847 rows.** The paper publishes nothing beyond 5p3/2. The 836 rows whose ground configuration
needs 5d, 5f, 6s, 6p, 6d or 7s are not predicted. These are all at Z ≥ 55 and N ≥ 55.

| Variant | MAPE, covered rows | Median APE | Neutral atoms (54) | N≤10 | Z≥55 | S1 test | S2 test | S3 test |
|---|---|---|---|---|---|---|---|---|
| **delta, low j (primary)** | **2.82 %** | 0.85 % | 27.1 % | 0.40 % | 2.74 % | 2.83 % (1016/1188) | 2.74 % (3526/4362) | 3.02 % (793/928) |
| delta, stat j | 3.40 % | 1.41 % | 25.6 % | 2.06 % | 3.50 % | 3.43 % | 3.50 % | 3.64 % |
| Xα one-electron | 5.22 % | 2.28 % | 86.6 % | 3.32 % | 3.82 % | 5.12 % | 3.82 % | 5.19 % |

**Primary variant, by group:**
- By removed subshell: s 4.14 %, p 2.23 %, d 1.69 %, f 5.75 %.
- Hydrogen-like ions: 0.19 %.
- Experimental-status rows: 11.3 %.
- N ≥ 37: 5.68 %.

**Same 5011 rows for the frozen models** (`results/bench_mendoza2011_common_subset.json`, in-sample):

| Model | All covered rows | Neutral atoms | N≤10 | Experimental rows |
|---|---|---|---|---|
| pa_hier_rel | 1.62 % | 8.37 % | 0.44 % | 4.37 % |
| pa_bound9 | 2.63 % | 11.45 % | 0.51 % | 7.17 % |

## Caveats
- **S1/S2/S3 are not held out for this model.** The NRSHM constants were fitted by a genetic algorithm to 61,350
  NIST + FAC energies of isoelectronic sequences from He to Eu, with Z up to 92. Those data include heavy ions, so
  for this model the S1/S2/S3 numbers are only subsets of the rows, not blind tests. The paper also left neutral atoms
  and singly charged ions out of the fit, which explains the large error on neutral atoms.
- **Incomplete coverage.** About 14 % of the rows (heavy, near-neutral ions with n ≥ 5 d/f or n ≥ 6 shells) are not
  predicted. No constants for those shells were invented.

**Confidence: high.** The constants come from the original paper, and six of its printed tables are reproduced to
rounding.
