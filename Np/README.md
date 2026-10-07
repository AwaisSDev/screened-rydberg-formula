# Screened Rydberg formula for successive ionization energies

A closed-form formula for the ionization energy IE(Z, N) of any atom or ion. It needs only Z and the electron
configuration, which defaults to the NIST ground configuration.

**Author:** Muhammad Awais (independent researcher, Multan, Pakistan; mawais9171@gmail.com).

**Paper:** `docs/paper_draft.pdf` (source `docs/paper_draft.md`).

## The formula (final model, 33 global parameters)

    IE   = mu(Z) * { Ry Zeff^2/n^2 * F_nj(Zeff) * [1 + r_c (Z alpha)^2/n (Zeff/Za - 1)] + Ry x_l K_l(k)/n^2 }
           - [dE_QED + dE_FNS](Z, n) * (Zeff/Z)^2                 (QED/FNS term: 1s and 2s removal only)
    Zeff = Z - sigma1 - D,      D = T / (Za + kappa + |T|/h),      Za = Z - N + 1
    T    = sum_g tau_g nu_g + sum_c dtau_c nu_c,      h = (N-1) - sigma1 if T >= 0 else sigma1

The terms:
- **σ₁** is the *exact* first-order (1/Z) screening constant, a rational number for each configuration.
- **D** is a bounded fitted remainder. It guarantees Za ≤ Zeff ≤ Z.
- **ν_g / ν_c** are electron counts in 5 screening groups and 19 classes (paper §2.2, Table A2).
- **F** is the exact Dirac/Schrödinger ratio.
- **K_l** is the Hund-kink term.

Parameters: paper Table A1, or `results/pa_params.json`. Every per-row input, including σ₁, μ and QED/FNS, is in
`results/model_inputs.csv`. `tools/verify_from_inputs.py` rebuilds the model from that table and the printed equation
alone.

## Install

You need Python 3.11+ with numpy and scipy (`pip install numpy scipy`). numba is optional: without it, set
`PYTHONPATH=tools/numba_stub`. The numerics are identical; the relevant results are cached.

## Use

```bash
python ionization.py Fe                       # all 26 successive IEs, predicted vs NIST
python ionization.py O --N 8 --details        # every intermediate quantity + low-confidence warnings
python ionization.py Na --N 11 --details --config "[Ne]3s1"
python ionization.py 118                      # beyond the NIST table (Madelung configuration, flagged)
```

```python
from ionization import ie, ionization_energy, successive_ionization_energies
d = ie("Fe")      # dict: IE_eV, NIST_eV, sigma1, nu, T, h, D, Zeff, rydberg_eV, F_dirac, R_rel, hund_eV, warnings
```

## Five examples checked against NIST

These come from `python tools/readme_examples.py`.

| ion | removed | Z_eff | predicted (eV) | NIST (eV) | error |
|---|---|---|---|---|---|
| O I | 2p | 2.0537 | 13.479 | 13.618 | -1.02 % |
| Na I | 3s | 2.0266 | 6.203 | 5.139 | +20.70 % |
| Fe I | 4s | 3.0039 | 7.611 | 7.902 | -3.69 % |
| Fe XVII | 2p | 19.2513 | 1258.687 | 1262.700 | -0.32 % |
| Pb I | 6p | 5.1339 | 8.012 | 7.417 | +8.02 % |

Accuracy on all 5847 NIST successive IEs (Z = 1–110):
- 1.87 % mean error, 0.94 % median;
- 7.5 % on neutral first IEs;
- 0.0012 % on H-like ions;
- 6.60 % in a blind extrapolation (fit on Z ≤ 54, predict Z ≥ 55).

## Limitations (read before use)

- **Neutral and near-neutral atoms are the weak spot.** Na is +21 %, Ca −18 % and Rn −25 %. Heavy neutrals had a
  21.8 % mean error in the blind test.
- **f-electron removal** is the least accurate class: 3.5 % mean error, and 2–2.5× too high in extrapolation.
- **Z > 110** is outside the NIST table; values there are qualitative. Og comes out at 3.66 eV, far below Rn.
- **All fitted relativistic coefficients are negative.** For Pb, part of the accuracy comes from this term
  cancelling an overestimate.
- **22 monotonicity violations** (Pt–Bi, Rf–Ds); `ie()` flags them.
- **On the ions Mendoza et al. (2011) cover, the formula is competitive with their model, not decisively better**
  (paper §4.5).
- **`ie()` warns** for heavy neutrals, f removal, Z > 110, rearranged configurations and monotonicity violations.

## Tests

```bash
python tests/test_paper_examples.py     # worked examples of the paper + reference implementation
python tests/test_coverage.py           # all 7021 ions with Z <= 118: finite, positive
```

## Citation

```
M. Awais, "A screened Rydberg formula for the successive ionization energies of all atoms and ions,"
manuscript, 2026. Code and data: https://github.com/AwaisSDev/Chem-Research (folder Np/).
```

Data: NIST Atomic Spectra Database, version 5.12 (Kramida, Ralchenko, Reader and NIST ASD Team),
doi:10.18434/T4W30F.

AI assistance: code, analysis and draft text were produced with Anthropic Claude under the author's direction. The
author is responsible for the content.
