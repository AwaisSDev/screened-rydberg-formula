# Track A: ionization energies from first principles

**Question.** How far does pure theory get toward a closed-form $IE(Z,N,\text{configuration})$, and where does empirical input become unavoidable?

**Answer in one paragraph.**
1. The non-relativistic Hamiltonian, scaled by $r\to r/Z$, gives the exact asymptotic series $IE = Z^2\Delta E_0 + Z\Delta E_1 + \Delta E_2 + \dots$ (hartree).
2. $\Delta E_0 = 1/(2n^2)$ is exactly Bohr's formula. $\Delta E_1$ is an exact rational number for every configuration, built here from hydrogenic Slater integrals and Racah algebra. It yields a first-principles screening constant, $\sigma_1 = -\Delta E_1/(2\Delta E_0)$, tabulated for $N = 1..110$.
3. With Dirac + recoil + finite-size + QED corrections for the removed electron, the two-term (or completed-square) series gives a closed-form formula with zero fitted parameters, $IE \approx \mu\,\Delta E_0 (Z-\sigma_1)^2 + \Delta_{\rm rel}$.
   - It is accurate to about 0.2-0.9% for few-electron highly charged ions ($N/Z \lesssim 0.3$).
   - It fails for near-neutral atoms (median error about 1000% at $N = Z$). This is not a coding defect: the series is asymptotic in $N/Z$.
4. No exact closed form exists for $N\ge 2$. $\Delta E_2$ requires continuum sums, and fitting it per sequence only moves the problem into data. Even then the fitted term cannot extrapolate to neutral atoms (held-out neutral MAPE about 750%).
5. The first-principles route that works across the whole range is my own Kohn-Sham LSDA $\Delta$SCF solver. It gives 0.8% MAPE for all ions with $Z\le 18$ and 3.3% for neutral first IEs with $Z\le 54$, with no fitted parameters, but it is a numerical procedure, not a formula.

Code: `models/first_principles/`. Results: `results/fp_*`. Figures: `results/figures/fp_*.png`. Every model is callable as `predict(Z, N, shells=None)` (see §8).

**Note on a brief typo.** The Li-like first-order coefficient is $E_1(1s^2 2s) = 5965/5832 = 1.0228052$, which equals $5/8 + 2\cdot 17/81 - 16/729$. The task brief's "5965/972" was a typo, confirmed by the orchestrator. The code reproduces $5965/5832$ exactly, so $\Delta E_1(\text{Li-like}) = 0.625 - 1.022805 = -0.397805$.

---

## 1. The 1/Z (Hylleraas-Layzer) expansion

### 1.1 Hamiltonian scaling

The non-relativistic, infinite-nuclear-mass Hamiltonian of $N$ electrons and nuclear charge $Z$ in atomic units is

$$H = \sum_{i=1}^{N}\left(-\tfrac12\nabla_i^2 - \frac{Z}{r_i}\right) + \sum_{i<j}\frac{1}{r_{ij}}.$$

Substitute $\mathbf r_i = \boldsymbol\rho_i/Z$. Then $\nabla_r^2 = Z^2\nabla_\rho^2$ and $1/r = Z/\rho$, so

$$H = Z^2\Big[\underbrace{\sum_i\left(-\tfrac12\nabla_{\rho_i}^2 - \tfrac1{\rho_i}\right)}_{H_0} + \frac1Z\underbrace{\sum_{i<j}\frac1{\rho_{ij}}}_{V}\Big].$$

Rayleigh-Schrödinger perturbation theory in $\lambda = 1/Z$ gives, for a fixed configuration and term,

$$E(Z,N) = Z^2E_0 + Z E_1 + E_2 + \frac{E_3}{Z} + \dots \quad\text{(hartree)},$$

where:
- $E_0$ is the eigenvalue of $H_0$, a sum of hydrogenic $Z=1$ energies;
- $E_1 = \langle V\rangle$ is evaluated over the $Z=1$ hydrogenic zeroth-order state, using degenerate perturbation theory when needed (§1.4);
- $E_2 = -\sum_{k\neq 0}|\langle k|V|0\rangle|^2/(E^{(0)}_k-E^{(0)}_0)$ includes the hydrogenic continuum.

The ionization energy is

$$IE(Z,N) = E(Z,N-1) - E(Z,N) = Z^2\Delta E_0 + Z\Delta E_1 + \Delta E_2 + \dots,\qquad \Delta E_k \equiv E_k(N-1)-E_k(N).$$

### 1.2 Zeroth order: Bohr

$$E_0 = -\sum_{\text{subshells}} \frac{q_{nl}}{2n^2},\qquad \Delta E_0 = \frac{1}{2n^2}\ \ \text{(removal of one electron from shell } n).$$

So $Z^2\Delta E_0$ in hartree equals $Ry\,Z^2/n^2$ in Rydberg units. This is Bohr's formula, and it is exactly the leading term.

For the rearranged ions (e.g. V $3d^34s^2 \to$ V$^+$ $3d^4$), the ground configurations of the N- and (N-1)-electron ions differ by more than one electron. Then $\Delta E_0 = E_0(\text{ion}) - E_0(\text{atom})$ is computed literally, e.g. $1/32 - 1/18 + 1/32 = 0.00694$ for V. These rows are where the expansion is least meaningful.

### 1.3 First order: exact hydrogenic Slater integrals

With $P_{nl}(r) = rR_{nl}(r)$ the $Z=1$ hydrogen radial functions,

$$P_{nl}(r) = \mathcal N_{nl}\Big(\frac{2}{n}\Big)^{l} r^{l+1}e^{-r/n}L^{(2l+1)}_{n-l-1}(2r/n),\qquad \mathcal N^2_{nl} = \Big(\frac2n\Big)^3\frac{(n-l-1)!}{2n\,(n+l)!}.$$

The radial two-electron integrals are

$$R^k(ab;cd) = \int_0^\infty\!\!\int_0^\infty P_a(r_1)P_c(r_1)P_b(r_2)P_d(r_2)\,\frac{r_<^k}{r_>^{k+1}}\,dr_1dr_2,$$

with direct $F^k(a,b) = R^k(ab;ab)$ and exchange $G^k(a,b) = R^k(ab;ba)$. Every $P_aP_c$ is a polynomial times $e^{-\beta r}$ with rational $\beta = 1/n_a + 1/n_c$. Using

$$\int_0^x t^s e^{-bt}dt = \frac{s!}{b^{s+1}}\Big[1-e^{-bx}\sum_{j\le s}\frac{(bx)^j}{j!}\Big],$$

every $F^k$ and $G^k$ is an exact rational number. The code (`slater_integrals.py`) evaluates them in `fractions.Fraction` arithmetic and caches them in `cache/slater_integrals.json`. For charge $Z$ every integral scales as $Z\cdot R^k(Z{=}1)$.

Checks:
- $F^0(1s,1s) = 5/8$, $F^0(1s,2s) = 17/81$, $G^0(1s,2s) = 16/729$;
- $F^0(1s,2p) = 59/243$, $G^1(1s,2p) = 112/2187$, $F^0(2s,2s) = 77/512$, $F^2(2p,2p) = 45/512$;
- a general $G^1(3d,4p)$ agrees with brute-force quadrature to $3\times10^{-12}$.

**Angular algebra** (`angular.py`, Condon-Shortley phases). The Coulomb coefficients are

$$c^k(lm,l'm') = (-1)^m\sqrt{(2l+1)(2l'+1)}\begin{pmatrix}l&k&l'\\0&0&0\end{pmatrix}\begin{pmatrix}l&k&l'\\-m&m-m'&m'\end{pmatrix}.$$

The spin-orbital matrix elements are

$$\langle ab|r_{12}^{-1}|cd\rangle = \delta_{\sigma_a\sigma_c}\delta_{\sigma_b\sigma_d}\delta_{m_a+m_b,m_c+m_d}\sum_k c^k(l_am_a,l_cm_c)\,c^k(l_dm_d,l_bm_b)\,R^k(ab;cd).$$

These are evaluated exactly (sympy Wigner 3j) for diagonal energies.

**Configuration-average energy** (Slater; Cowan, *Theory of Atomic Structure and Spectra*, eq. 6.38):

$$E_{\rm av} = \sum_{a}\binom{q_a}{2}\Big[F^0(aa) - \frac{2l_a+1}{4l_a+1}\sum_{k>0}\begin{pmatrix}l_a&k&l_a\\0&0&0\end{pmatrix}^2F^k(aa)\Big] + \sum_{a<b}q_aq_b\Big[F^0(ab)-\frac12\sum_k\begin{pmatrix}l_a&k&l_b\\0&0&0\end{pmatrix}^2G^k(ab)\Big].$$

**Term correction (Hund's-rule ground term).** For closed subshells, $E_1$ equals the configuration average. For open subshells, the code builds all determinants with fixed subshell occupations and with $M_S = S$, $M_L = L$ of the Hund term: maximum $S$, then maximum $L$, with all open subshells coupled high-spin. It then projects onto the highest-weight subspace, the null space of $S_+$ and $L_+$, which is the pure $|S L M_S{=}S M_L{=}L\rangle$ space. Finally it diagonalises $V$ there (`config_energy.E1_single_config`).
- When that space is a single determinant, the energy is an exact rational: $E_1 = E_{\rm core} + \sum_{o}u_{\rm core}(o) + \sum_{o<o'}(J_{oo'} - \delta_{\sigma\sigma'}K_{oo'})$. This covers $p^1$, $p^2\,{}^3P$, $p^3\,{}^4S$, $d^5\,{}^6S$, and all single-open-subshell Hund terms of $p$ and $d$.
- Example: C $2p^2\ {}^3P$ gives $E_1 = 10957679/3359232 = 3.261960$, versus the configuration average 3.272506.
- For two or more open subshells (e.g. $3d^64s^2$ or $4f^n5d6s^2$), the high-spin coupled term is used. That is the documented approximation for multiple open subshells.

### 1.4 Layzer complex degeneracy

In the hydrogenic limit, all configurations with the same number of electrons in each principal shell $n$ are degenerate: $E_0$ depends only on $n$. Examples are $1s^22s^2$ and $1s^22p^2$ (${}^1S$, Be-like), and $3s^2$ and $3p^2$ (Mg-like). First-order perturbation theory must therefore diagonalise $V$ within the whole complex. The code (`config_energy.E1_complex`) does this whenever exactly one $n$-shell is open and the complex has at most 60,000 determinants. The space is restricted to the $S$, $L$ and parity of the reference Hund term.
- Be-like: $E_1$ falls from the single-configuration $586373/373248 = 1.571001$ to $1.559274$ (complex). This matches Layzer's known $E_1 = 1.5592742$ for Be ${}^1S$. The ground eigenvector has 94.9% weight on $2s^2$ and 5.1% on $2p^2$.
- B-like: $2.334449 \to 2.327527$. C-like ${}^3P$: $3.261960 \to 3.258865$. Mg-like: $10.567378 \to 10.563805$.

This first-order mixing is the origin of the well-known Be-like $2s^2$-$2p^2$ near-degeneracy correlation. When the complex is not available (two open $n$ shells, or too large), the single-configuration value is used. The CSV columns `dE1_cx`/`E1_Layzer_complex` show which applies.

**Validation of first order:**

| quantity | this work (exact) | literature |
|---|---|---|
| He-like $E_1(1s^2)$ | $5/8$ | $5/8$ |
| Li-like $E_1(1s^22s)$ | $5965/5832 = 1.0228052$ | $1.0228052$ (brief's 5965/972 was a typo) |
| Be-like $E_1(1s^22s^2)$, complex | 1.5592742 | 1.5592742 (Layzer 1959) |
| He-like $E_2$ | not computed (continuum sum) | $-0.157666429$ (Scherr-Knight 1963; Baker et al. 1990) |

### 1.5 The first-principles screening constant

Completing the square in the first two terms gives

$$IE \simeq \Delta E_0\,(Z-\sigma_1)^2,\qquad \sigma_1 = -\frac{\Delta E_1}{2\Delta E_0} = -n^2\Delta E_1\ \ (\text{when }\Delta E_0 = 1/2n^2).$$

This reproduces the series exactly through order $Z$ and implies the specific closed-form guess $\Delta E_2^{\rm sq} = \Delta E_1^2/(4\Delta E_0)$. $\sigma_1$ is the exact leading-order (bare-nucleus-limit) screening of the removed electron. Slater's 1930 rules are an empirical fit to the opposite limit, the neutral atom.

Comparison for neutral-atom configurations (full table in §6; figure `results/figures/fp_sigma_abinitio_vs_slater.png`):
- Inner electrons: the ab initio value is larger than Slater's. For He 1s, 0.625 vs 0.30. For Ne 2p, 6.55 vs 4.15: hydrogenic same-shell screening is stronger than Slater's 0.35 per electron.
- Valence electrons of heavy neutral atoms: it is far smaller. Na 3s gives 7.79 vs Slater 8.80 vs the "experimental" $Z - n\sqrt{IE/Ry} = 9.16$. Cs 6s gives 42.8 vs 52.8. The reason is that bare-nucleus hydrogenic orbitals put the outer electron at the same length scale as the core, so it is barely screened.

Real screening of the valence electron in neutral atoms is a strongly non-perturbative, all-order effect. This is the quantitative reason the expansion fails there.

**Exact checks of $\Delta E_2$.** Known values are $E_2(\text{He-like}) = -0.157666429$ and $E_2(\text{Li-like }1s^22s) \approx -0.408165$ (Yan, Tambasco & Drake 1998; quoted to about $10^{-5}$, not re-verified here). These give $\Delta E_2 = 0.157666$ (He-like) and $0.250499$ (Li-like). The completed-square guess gives $0.1953$ and $0.3165$: the right sign and magnitude, but 24-26% too large.

### 1.6 Second order: what is and is not possible in closed form

$E_2$ is a sum over the full hydrogenic spectrum, continuum included, of squared first-order couplings. It equals the second-order pair (and, for open shells, single-excitation) energies. It has no finite closed form for any $N\ge 2$. Even for He, $E_2$ is known only numerically, to many digits. Four options were explored:

1. **Literature values** (zero fitted parameters): exact $E_2$ for He-like and Li-like, in model `fp_zexp3_lit`. He-like MAPE falls from 0.29% to 0.22% (max error 4.6% to 0.49%). Li-like falls from 0.69% to 0.29%.
2. **Closed-form estimate** $\Delta E_2 = \Delta E_1^2/(4\Delta E_0)$, i.e. $(Z-\sigma_1)^2$ (model `fp_zexp`, zero parameters). It reduces the overall MAPE from 653% (two terms) to 70%.
3. **One fitted $\Delta E_2$ per sequence/configuration key** (semi-empirical). The key is the pair (N-electron configuration, (N-1)-electron configuration). Weighted least squares on relative error gives a closed-form value per key.
   - S1 is meaningful (1145 of 1188 test rows have a fitted key).
   - S2 is partially meaningful.
   - **S3 is impossible by construction**: unseen sequences have no key. The model falls back to option 2 and the S3 numbers equal option 2's.
4. **A smooth universal function** with 12 fitted parameters: $\sigma = \sigma_1 + \sum_{k=1}^3 c_{l,k}(N/Z)^k$ with separate coefficients for removed $l = s,p,d,f$. It is fitted in log space and validated on all splits.

---

## 2. Relativistic, nuclear and QED corrections

### 2.1 One-electron ions (`relativity.py`)

- **Dirac point nucleus** (Sommerfeld): $E_{n\kappa} = c^2\Big\{\big[1+\big(\tfrac{Z\alpha}{n-|\kappa|+\sqrt{\kappa^2-Z^2\alpha^2}}\big)^2\big]^{-1/2}-1\Big\}$, with $c = 1/\alpha$.
- **Finite nuclear size (FNS):** the radial Dirac equation is solved numerically (RK4 shooting on a log grid) for a uniformly charged sphere of radius $R_{\rm sph} = \sqrt{5/3}\,R_{\rm rms}$, minus the same calculation for a point nucleus. The $R_{\rm rms}$ values (Angeli-Marinova, via Yerokhin & Shabaev 2015 Table II) are external nuclear data. The numerical point-nucleus solution reproduces the analytic Dirac energy to $10^{-13}$. The 1s shift is 0.0146 eV at Z=20, 1.96 eV at Z=50 and 199.0 eV at Z=92.
- **Recoil** (finite nuclear mass, Barker-Glover form): $B = \mu c^2(1-f) + \mu^2c^2(1-f)^2/[2(M+m)]$ with $f = 1 + E_{n\kappa}/c^2$ and $\mu = mM/(m+M)$.
- **QED, vacuum polarisation:** the Uehling potential
  $V_U(r) = -\frac{2\alpha Z}{3\pi r}\int_1^\infty dt\,e^{-2crt}\big(1+\frac{1}{2t^2}\big)\frac{\sqrt{t^2-1}}{t^2}$
  is integrated analytically over the Dirac 1s density, with one numerical $t$ quadrature. It reproduces the tabulated Uehling function $F_U$ to all 6 printed digits for Z=1, 10, 50, 92, 110.
- **QED, one-loop self-energy:** $\Delta E_{SE} = \frac{\alpha}{\pi}\frac{(Z\alpha)^4}{n^3}F_{SE}(Z\alpha)\,mc^2$. $F_{SE}(1s)$ and $F_{SE}(2s)$ are taken per $Z$ from Mohr's all-order results as tabulated by Yerokhin & Shabaev, *J. Phys. Chem. Ref. Data* 44, 033103 (2015) [arXiv:1506.01885], Table II, "SE(pnt)" rows. **These 2×110 tabulated values are external theory inputs, not fitted parameters.** $n \ge 3$ uses $F_{SE}(2s)$. Self-energy of $p$ and $d$ electrons is neglected.
- **Closed-form alternative for the self-energy (no table):** the low-order $Z\alpha$ expansion
  $F = A_{41}\ln(Z\alpha)^{-2} + A_{40} + A_{50}Z\alpha + (Z\alpha)^2[A_{62}\ln^2(Z\alpha)^{-2} + A_{61}\ln(Z\alpha)^{-2}]$,
  with $A_{41} = 4/3$, $A_{40} = 10/9 - \frac43\ln k_0$ ($\ln k_0(1s) = 2.984128556$, $\ln k_0(2s) = 2.811769893$), $A_{50} = 4\pi(139/128 - \frac12\ln2)$, $A_{62} = -1$ and $A_{61}(1s) = \frac{28}{3}\ln2 - \frac{21}{20}$ (`relativity.F_SE_closed`).

**H-like accuracy in layers** (110 ions, Z = 1..110; `results/fp_hlike_*`; figure `results/figures/fp_hlike_error.png`):

| layer | MAPE % | median APE % | max APE % |
|---|---|---|---|
| naive Bohr $Ry\,Z^2$ | 6.0 | 4.2 | 19.5 (Z=110) |
| (a) Dirac point nucleus only | 0.190 | 0.118 | 0.91 |
| (b) + recoil + finite nuclear size (numerical Dirac) | 0.113 | 0.109 | 0.248 |
| (c) + one-loop QED (Uehling computed; $F_{SE}$ from the Yerokhin-Shabaev table) | **0.00117** | **0.000154** | **0.0096** (Z=110) |
| (d) + QED with closed-form $Z\alpha$-expansion self-energy | 0.354 | 0.191 | 1.18 |

- Layer (c) reaches the brief's $\le 10^{-4}$ target for every Z: 12 ppm mean, 96 ppm worst.
- **Caveat:** NIST's H-like reference IEs come from the same Yerokhin-Shabaev QED theory. Layer (c)'s agreement is therefore *consistency with that source*, not independent validation. The residual about $10^{-5}$-$10^{-4}$ at high Z is the omitted two-loop QED, nuclear polarisation and recoil-QED terms, which are tabulated there but not used here.
- Layer (d) is honest about the closed form. The $Z\alpha$ expansion matches the table at Z=1 (10.318 vs 10.317) but diverges for $Z\gtrsim 15$ (Z=30: 3.56 vs 2.55; Z=92: 6.09 vs 1.49). It is worse than no QED at all for heavy ions. No useful closed-form all-order self-energy exists: the all-order function must be computed numerically or tabulated.

### 2.2 Many-electron ions: relativistic Z-expansion

For the removed electron (quantum numbers $n,l$; $j$ from Hund's rule), the correction added to the non-relativistic series is

$$\Delta_{\rm rel} = \Big[B_{\rm Dirac}(n,\kappa,Z) - \frac{Z^2}{2n^2} - \Delta E_{\rm QED} - \Delta E_{\rm FNS}\Big]\Big(\frac{Z-\sigma_1}{Z}\Big)^{p}.$$

The non-relativistic part is multiplied by the reduced-mass factor $\mu = M/(M+m)$.
- $p = 2$ is the Landé penetrating-orbital scaling, $\Delta E_{\rm rel}\propto Z_i^2Z_o^2/n^3$. It is the default, with no fitted parameter.
- $p = 4$ is pure hydrogenic $Z_{\rm eff}^4$ scaling.
- Screening uses the first-order $\sigma_1$, so no fitted parameter enters.

Not included: the relativistic correction to the electron-electron interaction (Breit, and the $\alpha^2Z^3$ relativistic part of $E_1$). This is what limits the few-electron sequences at high Z. The plateau at 0.2-0.5% for $Z\gtrsim 60$ is visible in `results/figures/fp_zexp_isoelectronic.png`.

---

## 3. Kohn-Sham DFT (own radial solver, `ks_atom.py`)

**Equations** (hartree units). The spherically averaged, spin-polarised KS equations are

$$\Big[-\tfrac12\frac{d^2}{dr^2} + \frac{l(l+1)}{2r^2} + v_{s\sigma}(r)\Big]P_{nl\sigma} = \varepsilon_{nl\sigma}P_{nl\sigma},\qquad v_{s\sigma} = -\frac Zr + v_H[\rho] + v_{xc,\sigma}[\rho_\uparrow,\rho_\downarrow],$$

with $\rho_\sigma = \sum f_{nl\sigma}P^2_{nl\sigma}/(4\pi r^2)$. Fractional occupations give spherical averaging; open subshells are filled high-spin (Hund).

The total energy is $E = T_s + E_{\rm ext} + E_H + E_{xc}$, with $T_s = \sum f\varepsilon - \sum_\sigma\int\rho_\sigma v_{s\sigma}$.

**Numerics:**
- Grid: logarithmic, $r_i = r_{\min}e^{ih}$, $r_{\min} = 10^{-6}/Z$, $r_{\max} = 60$.
- Radial solver: with $P = r^{1/2}y$ the equation becomes $y'' = [(l+\tfrac12)^2 + 2r^2(v-\varepsilon)]y$, solved by Numerov.
- Eigenvalues: bisection on the Sturm node count. Outward integration is stopped once $h^2g/12 > 0.2$; this fixed a spurious-node instability for deep states.
- Eigenvectors: outward/inward Numerov matched at the classical turning point.
- Hartree potential: radial Poisson equation by inhomogeneous Numerov, with the boundary condition $r v_H\to N_{el}$.
- Functionals: LDA/LSDA = Slater exchange + VWN5 correlation, with $f(\zeta)$ spin interpolation and spin stiffness.
- Mixing: Pulay/DIIS on the potential.
- Scalar-relativistic Koelling-Harmon option (`scalar_rel=True`): implemented but not used for the production numbers, per the time box. The two-component KH equations are integrated by RK4 in $\ln r$. They reproduce the exact Dirac 1s energy for H-like Z=92 to $10^{-13}$. Without the MacDonald-Vosko relativistic exchange correction used by NIST's ScRLDA, total energies differ from NIST ScRLDA by 0.02 Ha (O) to 51 Ha (Rn), so this mode is not validated against ScRLDA.

**What the validation proves.** Agreement with the NIST "Atomic Reference Data for Electronic Structure Calculations" (Kotochigova, Levine, Shirley, Stiles & Clark), which uses the same Slater + VWN functional, **only verifies the implementation, not physical accuracy**.

| atom | this solver | NIST (same functional) | difference |
|---|---|---|---|
| He LDA | -2.834836 | -2.834836 | $3.8\times10^{-7}$ |
| Be LDA | -14.447209 | -14.447209 | $-4.7\times10^{-7}$ |
| Ne LDA | -128.233481 | -128.233481 | $-2.7\times10^{-7}$ |
| Ar LDA | -525.946195 | -525.946195 | $7\times10^{-8}$ |
| H LSD / LDA | -0.478671 / -0.445671 | -0.478671 / -0.445671 | $<10^{-6}$ |
| O, O$^+$ LSD | -74.527410, -74.016721 | -74.527410, -74.016721 | $<10^{-6}$ |
| Fe, Fe$^+$ LSD | -1261.223291, -1260.927746 | identical | $<10^{-6}$ |
| Kr, Rn LDA | -2750.147941, -21861.346870 | -2750.147940, -21861.346869 | $\sim10^{-6}$ |

Orbital eigenvalues also agree, e.g. Ar 1s gives $-113.800134$ in both. The production grid ($h = 0.009$) changes total energies by less than $10^{-6}$ Ha (Ar, Kr).

**Physical error of LDA.**
- Ne: LDA total energy $-128.2335$ Ha vs the exact non-relativistic $-128.9376$ Ha (Chakravorty & Davidson 1993). The error is $+0.70$ Ha (0.55%), mostly from LDA exchange.
- H: LSD gives $-0.4787$ vs exact $-0.5$, a self-interaction error of 4.3%. As a result, $\Delta$SCF IEs of one-electron ions are 1-4% too low.
- Most of the absolute total-energy error cancels in $\Delta$SCF differences. The IE errors against NIST data are below.

**Delta-SCF and Slater transition-state results.** All ions with $Z\le18$ (171 rows) plus neutral first IEs for $Z = 19..54$: 207 rows, all converged, run time about 150 s (`run_dft.py 18 54`).
- $\Delta$SCF: $IE = E_{KS}(N-1) - E_{KS}(N)$, using the NIST ground configurations of both ions.
- Transition state: $IE \approx -\varepsilon_{\rm HOMO}$ at half occupation of the removed spin-orbital (Janak's theorem plus the trapezoid rule). This is a near-closed-form link between IE and one orbital energy.

| model (`results/fp_*`) | all 207 | neutral first IE (54) | neutral, experimental status only (53) | ions $Z\le18$ (171) | H-like (18) |
|---|---|---|---|---|---|
| LSDA $\Delta$SCF (`dft`) | 1.32 / 0.73 | **3.30 / 2.09** | 3.29 / 2.04 (MAE 0.28 eV) | **0.80 / 0.55** | 1.36 / 0.93 |
| Slater transition state (`dft_ts`) | 1.16 / 0.55 | 3.00 / 1.77 | 3.03 / 1.75 (MAE 0.25 eV) | 0.66 / 0.49 | 0.82 / 0.66 |
| $-\varepsilon_{\rm HOMO}$ (`dft_homo`) | 17.7 / 11.3 | 39.1 / 38.4 | | 13.4 / 9.1 | 9.8 / 5.6 |

Entries are MAPE % / median APE %. Figures: `results/figures/fp_dft_neutral_first_IE.png`, `results/figures/fp_dft_ions_Zle18.png`.

- $-\varepsilon_{\rm HOMO}$ underestimates IEs by about 40% because LDA's potential decays exponentially, not as $-1/r$ (self-interaction).
- The transition-state value is marginally *better* than full $\Delta$SCF, a known LSDA error cancellation.
- Not done (per the time box): PZ-SIC, GGA, scalar-relativistic $\Delta$SCF for heavy atoms, and ions with $Z > 18$.

---

## 4. Accuracy of the analytic models on all 5847 rows

Entries are MAPE % / median APE % from `evaluate.py`. All variants without "fit" have **zero fitted parameters**.

| model | ALL | neutral 1st IE | H-like | N≤10 | 11≤N≤36 | N≥37 | removed n=1 | rearranged |
|---|---|---|---|---|---|---|---|---|
| `bohr`: $Z^2\Delta E_0$ | 1340 / 172 | 22500 / 21100 | 6.0 / 4.2 | 37.6 / 12.8 | 302 / 85 | 2820 / 669 | 7.6 / 5.2 | 2340 / 1100 |
| `zexp2`: $Z^2\Delta E_0 + Z\Delta E_1$ | 658 / 33 | 12500 / 11500 | 6.0 / 4.2 | 11.6 / 3.5 | 113 / 15 | 1420 / 231 | 6.1 / 4.3 | 1140 / 374 |
| `zexp2_rel`: + relativistic/QED | 653 / 31 | 12500 / 11500 | 0.0012 / 0.0002 | 7.2 / 0.46 | 109 / 9.5 | 1410 / 229 | 0.24 / 0.008 | 1130 / 395 |
| `zexp_sq`: $\Delta E_0(Z-\sigma_1)^2$ | 67.4 / 9.3 | 1060 / 1010 | 6.0 / 4.2 | 5.5 / 2.3 | 15.5 / 3.4 | 140 / 32 | 6.0 / 4.2 | 194 / 112 |
| **`zexp`**: $\mu\Delta E_0(Z-\sigma_1)^2 + \Delta_{\rm rel}$, $p=2$ (headline) | **69.8 / 7.2** | 1100 / 1030 | **0.0012 / 0.0002** | **1.71 / 0.71** | 15.4 / 3.4 | 147 / 35 | 0.15 / 0.010 | 207 / 128 |
| `zexp_sq_rel_p4`: same with $p = 4$ | 66.7 / 5.9 | 1060 / 1010 | 0.0012 / 0.0002 | 1.44 / 0.21 | 14.5 / 2.5 | 141 / 32 | 0.10 / 0.010 | 200 / 113 |
| `zexp3_lit`: + literature $E_2$ for N≤3 | 69.7 / 7.2 | 1100 / 1030 | 0.0012 / 0.0002 | 1.66 / 0.69 | 15.4 / 3.4 | 147 / 35 | 0.11 / 0.010 | 207 / 128 |
| `zexp_univ`: 12-parameter universal $\sigma(N/Z)$ (fit, in-sample) | 38.7 / 4.1 | 819 / 698 | 3.4 / 0.41 | 4.6 / 2.2 | 8.3 / 3.7 | 80 / 13 | 3.7 / 0.74 | 163 / 103 |
| `zexp3_seq`: 409 fitted $\Delta E_2$ (in-sample only, see note) | 4.7 / 2.3 | 8.0 / ~0 | 0.0012 / 0.0002 | 0.58 / 0.48 | 2.9 / 2.2 | 8.1 / 4.7 | 0.11 / 0.008 | 6.3 / 0.71 |

The `zexp3_seq` in-sample numbers are not a model quality: many configuration keys occur only once (neutral atoms), so the fit becomes a per-row lookup there. Use the held-out numbers below.

**Few-electron highly charged ions** (the regime where theory alone works):

| subset | `zexp2_rel` (2 terms) | `zexp` (completed square) | `zexp3_lit` (exact $E_2$) |
|---|---|---|---|
| He-like, N=2 (109) | 0.48 / 0.23 / max 17 | 0.29 / 0.22 / max 4.6 | **0.22 / 0.20 / max 0.49** |
| Li-like, N=3 (101) | 2.5 / 0.25 / max 135 | 0.69 / 0.25 / max 25 | **0.29 / 0.20 / max 8.1** |
| N≤10 and Z≥3N (938) | 0.75 / 0.40 / max 6.9 | 0.65 / 0.56 / max 2.1 | 0.64 / 0.55 / max 2.1 |
| N/Z ≤ 0.2 (1170) | **0.48 / 0.38 / max 2.1** | 0.90 / 0.91 / max 3.1 | 0.90 / 0.91 |
| N/Z > 0.8 (1219) | 2900 / 1060 | 297 / 133 | 297 / 133 |

Entries are MAPE / median / max APE %. Figures: `results/figures/fp_zexp_error_vs_N_over_Z.png` (median APE vs N/Z for 1, 2 and "3" terms) and `results/figures/fp_zexp_isoelectronic.png` (APE along the N = 2, 3, 4, 10, 11, 18, 28, 36, 46, 54 sequences).

**Where the truncated series works:**
- Two terms: median error below 1% only for $N/Z\lesssim 0.2$.
- Completed square: median error below 3% for $N/Z\lesssim 0.3$ and below 10% for $N/Z\lesssim 0.5$.
- At $N/Z\to1$ every truncation fails by orders of magnitude.
- For pure high-Z ions (N/Z ≤ 0.2) the plain two-term series beats the completed square. There the true $\Delta E_2$ is smaller than $\Delta E_1^2/4\Delta E_0$ (He-like 0.158 vs 0.195) and the missing Breit/relativistic $E_1$ terms dominate.

## 5. Held-out validation of the fitted second-order variants

Standard splits: S1 = test $Z\%5{=}0$; S2 = train $Z\le54$, test $Z\ge55$; S3 = test $N\%6{=}0$. Entries are test MAPE / median APE / neutral first-IE MAPE, in %. Source: `results/fp_zexp2nd_validation.json`.

| model | params (fit on train) | S1 | S2 | S3 |
|---|---|---|---|---|
| zero-parameter `zexp` (reference) | 0 | 68.3 / 7.0 / 1092 | 77.8 / 7.4 / 1772 | 74.3 / 8.3 / 1171 |
| per-key $\Delta E_2$ (`zexp3_seq`) | 366 (S1), 78 (S2), 343 (S3) | 26.6 / 2.5 / 757 | 75.7 / 3.9 / 1772 (2609 of 4362 test rows had a key) | impossible: 0 test rows have a key, so it equals the reference |
| universal $\sigma = \sigma_1 + \sum c_{l,k}(N/Z)^k$ (`zexp_univ`) | 12 | 38.6 / 4.0 / 851 | 65.3 / 4.3 / 1415 | 42.9 / 4.6 / 789 |

Fitted universal coefficients (all data): $c_{s} = (-6.90, 21.77, -14.11)$, $c_p = (-6.14, 20.39, -12.80)$, $c_d = (17.40, -46.33, 35.45)$, $c_f = (31.23, -78.32, 61.21)$, for powers $k = 1, 2, 3$ of $N/Z$.

Conclusion: second-order information helps the bulk (median APE about 7% down to about 4%), but **no smooth correction to the 1/Z series rescues neutral atoms**. A screened-hydrogenic data-driven form (Track B) is needed there, or DFT.

## 6. Table of exact coefficients, ground configurations N = 1..36

Neutral-atom (Z = N) NIST ground configurations. Columns:
- $E_0$, $E_1$: total-energy coefficients (hartree, Z=1 scaling). $E_1$ is the Hund-term single-configuration value: an exact rational where printed, otherwise a float whose exact rational is in `cache/config_E1.json`.
- $\Delta E_0$, $\Delta E_1$: IE coefficients. "Best" uses the Layzer complex when available.
- $\sigma_1 = -\Delta E_1/(2\Delta E_0)$.

The full N = 1..110 table, with both neutral and highest-charge configurations, is in `results/fp_zexp_coefficients.csv`. Per-row coefficients for all 5847 table rows are in `results/fp_zexp_rows_coefficients.csv`.

| N | configuration (Z = N) | removed | $E_0$ | $E_1$ | $\Delta E_0$ | $\Delta E_1$ | $\sigma_1$ ab initio | $\sigma$ Slater |
|---|---|---|---|---|---|---|---|---|
| 1 | 1s1 | 1s | -1/2 | 0 | 0.500000 | 0.000000 | 0.0000 | 0.00 |
| 2 | 1s2 | 1s | -1 | 5/8 | 0.500000 | -0.625000 | 0.6250 | 0.30 |
| 3 | 1s2.2s1 | 2s | -9/8 | 5965/5832 | 0.125000 | -0.397805 | 1.5912 | 1.70 |
| 4 | 1s2.2s2 | 2s | -5/4 | 586373/373248 | 0.125000 | -0.536469 | 2.1459 | 2.05 |
| 5 | 1s2.2s2.2p1 | 2p | -11/8 | 1960489/839808 | 0.125000 | -0.768252 | 3.0730 | 2.40 |
| 6 | 1s2.2s2.2p2 | 2p | -3/2 | 10957679/3359232 | 0.125000 | -0.931338 | 3.7254 | 2.75 |
| 7 | 1s2.2s2.2p3 | 2p | -13/8 | 2437421/559872 | 0.125000 | -1.094668 | 4.3787 | 3.10 |
| 8 | 1s2.2s2.2p4 | 2p | -7/4 | 4754911/839808 | 0.125000 | -1.308370 | 5.2335 | 3.45 |
| 9 | 1s2.2s2.2p5 | 2p | -15/8 | 11982943/1679616 | 0.125000 | -1.472432 | 5.8897 | 3.80 |
| 10 | 1s2.2s2.2p6 | 2p | -2 | 2455271/279936 | 0.125000 | -1.636495 | 6.5460 | 4.15 |
| 11 | [Ne]3s1 | 3s | -37/18 | 9.635901 | 0.055556 | -0.865071 | 7.7856 | 8.80 |
| 12 | [Ne]3s2 | 3s | -19/9 | 10.567378 | 0.055556 | -0.927904 | 8.3511 | 9.15 |
| 13 | [Ne]3s2.3p1 | 3p | -13/6 | 11.630399 | 0.055556 | -1.062197 | 9.5598 | 9.50 |
| 14 | [Ne]3s2.3p2 | 3p | -20/9 | 12.758090 | 0.055556 | -1.127250 | 10.1452 | 9.85 |
| 15 | [Ne]3s2.3p3 | 3p | -41/18 | 13.950451 | 0.055556 | -1.192313 | 10.7308 | 10.20 |
| 16 | [Ne]3s2.3p4 | 3p | -7/3 | 15.229075 | 0.055556 | -1.277698 | 11.4993 | 10.55 |
| 17 | [Ne]3s2.3p5 | 3p | -43/18 | 16.572369 | 0.055556 | -1.342608 | 12.0835 | 10.90 |
| 18 | [Ne]3s2.3p6 | 3p | -22/9 | 17.980333 | 0.055556 | -1.407535 | 12.6678 | 11.25 |
| 19 | [Ar]4s1 | 4s | -713/288 | 18.859971 | 0.031250 | -0.879639 | 14.0742 | 16.80 |
| 20 | [Ar]4s2 | 4s | -361/144 | 19.776882 | 0.031250 | -0.916910 | 14.6706 | 17.15 |
| 21 | [Ar]3d1.4s2 | 4s | -41/16 | 21.492902 | 0.031250 | -0.964077 | 15.4252 | 18.00 |
| 22 | [Ar]3d2.4s2 | 4s | -377/144 | 23.286948 | 0.031250 | -1.011243 | 16.1799 | 18.85 |
| 23 | [Ar]3d3.4s2 | (rearr.) | -385/144 | 25.164339 | 0.006944 | -0.209804 | 15.1059 | 19.70 |
| 24 | [Ar]3d5.4s1 | 4s | -793/288 | 28.011097 | 0.031250 | -1.111297 | 17.7808 | 21.05 |
| 25 | [Ar]3d5.4s2 | 4s | -401/144 | 29.163840 | 0.031250 | -1.152742 | 18.4439 | 21.40 |
| 26 | [Ar]3d6.4s2 | 4s | -409/144 | 31.310090 | 0.031250 | -1.199074 | 19.1852 | 22.25 |
| 27 | [Ar]3d7.4s2 | (rearr.) | -139/48 | 33.534366 | 0.006944 | -0.236912 | 17.0577 | 23.10 |
| 28 | [Ar]3d8.4s2 | (rearr.) | -425/144 | 35.841988 | 0.006944 | -0.247065 | 17.7887 | 23.95 |
| 29 | [Ar]3d10.4s1 | 4s | -97/32 | 39.317547 | 0.031250 | -1.347129 | 21.5541 | 25.30 |
| 30 | [Ar]3d10.4s2 | 4s | -49/16 | 40.701948 | 0.031250 | -1.382912 | 22.1266 | 25.65 |
| 31 | [Ar]3d10.4s2.4p1 | 4p | -99/32 | 42.182528 | 0.031250 | -1.480091 | 23.6815 | 26.00 |
| 32 | [Ar]3d10.4s2.4p2 | 4p | -25/8 | 43.698057 | 0.031250 | -1.515169 | 24.2427 | 26.35 |
| 33 | [Ar]3d10.4s2.4p3 | 4p | -101/32 | 45.248537 | 0.031250 | -1.550480 | 24.8077 | 26.70 |
| 34 | [Ar]3d10.4s2.4p4 | 4p | -51/16 | 46.845921 | 0.031250 | -1.597384 | 25.5581 | 27.05 |
| 35 | [Ar]3d10.4s2.4p5 | 4p | -103/32 | 48.478254 | 0.031250 | -1.632334 | 26.1173 | 27.40 |
| 36 | [Ar]3d10.4s2.4p6 | 4p | -13/4 | 50.145538 | 0.031250 | -1.667284 | 26.6765 | 27.75 |

"(rearr.)" marks V, Co and Ni, whose cation ground configuration is $3d^{n+1}$. For them $\Delta E_0 = E_0(3d^{n+1}) - E_0(3d^n4s^2) = 1/32 - 1/18 + 1/32$ is not a single-electron Bohr term, and $\sigma_1$ is not a screening constant in the usual sense.

## 7. Honest discussion

**What is truly first-principles.**
- $\Delta E_0$, $\Delta E_1$ and $\sigma_1$: exact rational numbers from the Coulomb Hamiltonian, for every configuration. The only approximation is high-spin coupling for several open subshells.
- The Layzer-complex diagonalisation.
- Dirac energies, recoil and the finite-size shift (given nuclear radii), and the Uehling vacuum polarisation.
- The KS-LSDA solver. The functional itself (VWN5) is a fit to quantum Monte Carlo electron-gas data, but nothing is fitted to atoms.

**What needed external or empirical input.**
- Nuclear rms radii and masses (nuclear data).
- The all-order self-energy function $F_{SE}(Z\alpha)$: 2×110 tabulated values from Mohr / Yerokhin-Shabaev. These are external theory inputs, not fitted parameters. The closed-form $Z\alpha$ expansion that replaces them fails above $Z\approx15$.
- The literature $E_2$ for He- and Li-like ions.
- The choice $p = 2$ in the relativistic screening factor. This is physics-motivated (Landé); $p = 4$ scores slightly better and is reported as a variant, not adopted.
- The fitted second-order variants: 409 or 12 parameters, as labelled.

**Why no exact closed form exists for N ≥ 2.** With $V = \sum 1/r_{ij}$ the Schrödinger equation does not separate. $E(Z)$ is analytic in $1/Z$ only up to a critical charge. Already for He ($\lambda_c = 1/Z_c \approx 1.0975$), the coefficients $E_k$ for $k\ge2$ are infinite sums over the continuum with no closed form. Each order adds a new transcendental number (Bethe-logarithm-like sums for QED, pair sums for correlation).

**Why the expansion fails for neutral atoms.** The expansion parameter is effectively the ratio of electron-electron to electron-nucleus interaction, $\sim N/Z$. For a neutral atom it is of order 1, and the valence electron sees a nearly completely screened nucleus that perturbation theory around the bare-nucleus orbitals cannot reach at any low order. Quantitatively, $\sigma_1$ for the valence electron is 10-20% below the effective screening (Na 7.8 vs 9.2; Cs 42.8 vs 51.8). Because $IE\propto(Z-\sigma)^2$ and $Z-\sigma\approx 1$-$3$, a 10% error in $\sigma$ becomes a 100-1000% error in IE.

**Best closed-form approximation (Track A's answer).**

$$\boxed{IE(Z,N) \simeq \mu(Z)\,\frac{(Z-\sigma_1)^2}{2n^2} + \Big[B_{\rm Dirac}(n\kappa;Z) - \frac{Z^2}{2n^2} - \Delta E_{\rm QED} - \Delta E_{\rm FNS}\Big]\Big(\frac{Z-\sigma_1}{Z}\Big)^2,\qquad \sigma_1 = -n^2\Delta E_1\ \text{(exact rational)}}$$

- Hartree units; multiply by 27.211386 for eV.
- Zero fitted parameters.
- Exact through $O(Z)$, with all-order one-electron relativity.
- Valid for $N/Z\lesssim0.3$: median 0.7% for N≤10, 0.2% He-like, 12 ppm H-like.
- The exact $\sigma_1$ supplies the leading coefficients for Track B's screened-hydrogenic fit. $\sigma = \sigma_1 + O(N/Z)$ must reduce to $\sigma_1$ as $Z\to\infty$.

For near-neutral atoms the best first-principles route is numerical: LSDA $\Delta$SCF gives 3.3% MAPE on neutral first IEs ($Z\le54$) and 0.8% on all ions with $Z\le18$, with no parameters.

## 8. Code map and API

| file | content |
|---|---|
| `slater_integrals.py` | exact hydrogenic $F^k$, $G^k$, $R^k$ (Fraction), disk cache |
| `angular.py` | $c^k$ coefficients and 3j symbols (exact) |
| `config_energy.py` | $E_0$, $E_{\rm av}$, Hund-term $E_1$ (determinants + highest-weight projection), Layzer complex $E_1$ |
| `relativity.py` | Dirac energies, numerical Dirac FNS, recoil, Uehling, tabulated $F_{SE}$, closed-form $F_{SE}$, `hydrogenic_binding(Z, n, l, j2)` |
| `zexp.py` | `coeffs(Z, N, shells=None)` returns $\Delta E_0$, $\Delta E_1$ (sc, avg, complex), $\sigma_1$, $n$, $l$, $j$; `predict(Z, N, shells=None, variant='zexp_sq_rel', p=2)` returns IE in eV, with variants `bohr`, `zexp2`, `zexp2_rel`, `zexp_sq`, `zexp_sq_rel` |
| `zexp2nd.py` | literature $E_2$; per-key and universal second-order fits; split validation |
| `ks_atom.py` | radial KS LDA/LSDA solver: `total_energy(Z, shells, spin=True, ...)`; `Atom(...).scf()`; KH scalar-relativistic option |
| `run_dft.py` | $\Delta$SCF, transition state and HOMO for chosen rows: `py -3.13 run_dft.py ZMAX_ALL ZMAX_NEUTRAL` |
| `run_zexp.py`, `hlike_layers.py`, `coeff_table.py`, `build_coeffs.py`, `make_figures.py` | drivers |

All `results/fp_*` prediction files come from the current code version. They were regenerated from a clean state after the duplicate-agent episode, and the DFT run was re-run from scratch.
