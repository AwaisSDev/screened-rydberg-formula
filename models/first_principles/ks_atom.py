"""Spherical (radial) Kohn-Sham DFT solver for atoms and ions, written from scratch.

Equations (hartree atomic units).  For each spin sigma and orbital (n, l) with occupation f:

    [-1/2 d^2/dr^2 + l(l+1)/(2 r^2) + v_s,sigma(r)] P_nl,sigma(r) = eps_nl,sigma P_nl,sigma(r)
    v_s,sigma = -Z/r + v_H[rho] + v_xc,sigma[rho_up, rho_dn]
    rho_sigma(r) = sum f P^2 / (4 pi r^2)        (spherically averaged, fractional occupations)

Numerics
  * logarithmic grid r_i = r_min exp(i h); with P = r^{1/2} y the radial equation becomes
    y'' = g(x) y,  g = (l+1/2)^2 + 2 r^2 (v - eps), solved by Numerov (O(h^4)).
  * eigenvalues by bisection on the Sturm node count of the outward Numerov solution
    (# sign changes = # eigenvalues below eps), eigenvector by outward/inward integration
    matched at the outermost classical turning point.
  * Hartree potential from the radial Poisson equation U'' = -4 pi r rho (U = r v_H) by
    inhomogeneous Numerov, with the homogeneous solution a r added to impose U(inf) = N_el.
  * LDA/LSDA: Slater exchange + VWN5 correlation (Vosko, Wilk, Nusair 1980, Ceperley-Alder fit),
    spin interpolation with f(zeta) and the spin stiffness alpha_c.
  * SCF: Pulay/Anderson (DIIS) mixing of the potential.
  * Optional scalar-relativistic (Koelling-Harmon) mode, scalar_rel=True: the two-component
    radial equations dP/dr = 2 M c Q + P/r, dQ/dr = -Q/r + [l(l+1)/(2 M c r^2) + (v - eps)/c] P,
    M = 1 + (eps - v)/(2 c^2) (mass-velocity + Darwin, spin-orbit dropped), integrated by RK4 in
    x = ln r; eigenvalue by bisection on the node count of P; density from P^2 + Q^2.
    No relativistic correction to the XC functional (NIST ScRLDA includes MacDonald-Vosko).

Total energy  E = T_s + E_ext + E_H + E_xc,  T_s = sum f eps - sum_sigma Int rho_sigma v_s,sigma.
"""
import math

import numpy as np
from numba import njit

ALPHA = 7.2973525693e-3
C_LIGHT = 1.0 / ALPHA


# ============================================================================ grid
class Grid:
    def __init__(self, Z, rmin=None, rmax=60.0, h=0.0045):
        if rmin is None:
            rmin = 1e-6 / Z
        M = int(math.log(rmax / rmin) / h) + 1
        if M % 2 == 0:
            M += 1
        self.h = math.log(rmax / rmin) / (M - 1)
        self.x = np.arange(M) * self.h + math.log(rmin)
        self.r = np.exp(self.x)
        self.M = M
        # Simpson weights in x for int f dx
        w = np.ones(M)
        w[1:-1:2] = 4
        w[2:-1:2] = 2
        self.wx = w * self.h / 3.0
        self.wr = self.wx * self.r          # int f dr = sum wr f

    def integrate(self, f):
        return float(np.dot(self.wr, f))


# ============================================================================ Numerov kernels
@njit(cache=True)
def _count_nodes(V, r, h, l, E):
    M = r.shape[0]
    h12 = h * h / 12.0
    lp = (l + 0.5) ** 2
    y0 = r[0] ** (l + 0.5)
    y1 = r[1] ** (l + 0.5)
    g0 = lp + 2 * r[0] ** 2 * (V[0] - E)
    g1 = lp + 2 * r[1] ** 2 * (V[1] - E)
    nodes = 0
    for i in range(1, M - 1):
        g2 = lp + 2 * r[i + 1] ** 2 * (V[i + 1] - E)
        if h12 * g2 > 0.3 and V[i + 1] > E:
            # deep in the classically forbidden region: Numerov step unstable, no more nodes
            break
        y2 = (2 * y1 * (1 + 5 * h12 * g1) - y0 * (1 - h12 * g0)) / (1 - h12 * g2)
        if (y2 < 0) != (y1 < 0) and y2 != 0:
            nodes += 1
        if abs(y2) > 1e150:
            y2 *= 1e-150
            y1 *= 1e-150
        y0, y1 = y1, y2
        g0, g1 = g1, g2
    return nodes


@njit(cache=True)
def _eigen(V, r, h, l, k, Elo, Ehi, tol):
    """Find eigenvalue with k nodes by bisection between Elo (count<=k) and Ehi (count>k)."""
    for it in range(200):
        Em = 0.5 * (Elo + Ehi)
        if _count_nodes(V, r, h, l, Em) > k:
            Ehi = Em
        else:
            Elo = Em
        if Ehi - Elo < tol * max(1.0, abs(Em)):
            break
    return 0.5 * (Elo + Ehi)


@njit(cache=True)
def _wavefunction(V, r, h, l, E):
    """Numerov outward to turning point, inward from decay region; return normalised P(r)."""
    M = r.shape[0]
    h12 = h * h / 12.0
    lp = (l + 0.5) ** 2
    g = lp + 2 * r * r * (V - E)
    # outermost classical turning point
    it = M - 2
    while it > 2 and g[it] > 0:
        it -= 1
    if it < 10:
        it = 10
    # inward start: where integral of sqrt(g) dx from turning point exceeds 60
    acc = 0.0
    ist = M - 1
    for i in range(it, M):
        if g[i] > 0:
            acc += math.sqrt(g[i]) * h
        if acc > 60.0 or h12 * g[i] > 0.2:
            ist = i
            break
    y = np.zeros(M)
    y[0] = r[0] ** (l + 0.5)
    y[1] = r[1] ** (l + 0.5)
    for i in range(1, it + 1):
        y[i + 1] = (2 * y[i] * (1 + 5 * h12 * g[i]) - y[i - 1] * (1 - h12 * g[i - 1])) / (1 - h12 * g[i + 1])
        if abs(y[i + 1]) > 1e150:
            for j in range(i + 2):
                y[j] *= 1e-150
    yout = y[it]
    yin = np.zeros(M)
    if ist >= M - 1:
        ist = M - 1
        yin[ist] = 0.0
        yin[ist - 1] = 1e-200
    else:
        yin[ist] = 1e-200
        yin[ist - 1] = 1e-200 * math.exp(math.sqrt(g[ist]) * h)
    for i in range(ist - 1, it, -1):
        yin[i - 1] = (2 * yin[i] * (1 + 5 * h12 * g[i]) - yin[i + 1] * (1 - h12 * g[i + 1])) / (1 - h12 * g[i - 1])
        if abs(yin[i - 1]) > 1e150:
            for j in range(i - 1, ist + 1):
                yin[j] *= 1e-150
    sc = yout / yin[it]
    for i in range(it + 1, ist + 1):
        y[i] = yin[i] * sc
    for i in range(ist + 1, M):
        y[i] = 0.0
    P = y * np.sqrt(r)
    # normalise: int P^2 dr = int P^2 r dx  (Simpson)
    s = 0.0
    for i in range(M):
        w = 1.0 if (i == 0 or i == M - 1) else (4.0 if i % 2 == 1 else 2.0)
        s += w * P[i] * P[i] * r[i]
    s *= h / 3.0
    return P / math.sqrt(s)


@njit(cache=True)
def _poisson(rho, r, h, Q):
    """U(r) = r v_H(r) from U'' = -4 pi r rho with U(0)=0, U(inf)=Q. Uses U = r^{1/2} w."""
    M = r.shape[0]
    h12 = h * h / 12.0
    # w'' = w/4 - 4 pi r^{5/2} rho  ->  w'' = g w + s
    s = -4 * math.pi * r ** 2.5 * rho
    w = np.zeros(M)
    for i in range(1, M - 1):
        w[i + 1] = (2 * w[i] * (1 + 5 * h12 * 0.25) - w[i - 1] * (1 - h12 * 0.25)
                    + h12 * (s[i + 1] + 10 * s[i] + s[i - 1])) / (1 - h12 * 0.25)
    U = w * np.sqrt(r)
    a = (Q - U[M - 1]) / r[M - 1]
    return U + a * r



# ============================================================================ scalar-relativistic (Koelling-Harmon)
@njit(cache=True)
def _half_values(V, r):
    """V at x_{i+1/2} by 4-point interpolation of the smooth function r V(r)."""
    M = r.shape[0]
    u = V * r
    Vh = np.empty(M - 1)
    for i in range(M - 1):
        if i == 0 or i >= M - 2:
            uh = 0.5 * (u[i] + u[i + 1])
        else:
            uh = (-u[i - 1] + 9 * u[i] + 9 * u[i + 1] - u[i + 2]) / 16.0
        Vh[i] = uh / math.sqrt(r[i] * r[i + 1])
    return Vh


@njit(cache=True)
def _kh_rhs(r, P, Q, E, V, ll, c):
    """Koelling-Harmon radial equations in x = ln r (P = r g, Q = r f):
       dP/dx = 2 M c r Q + P,   dQ/dx = -Q + r [l(l+1)/(2 M c r^2) + (V - E)/c] P,
       M = 1 + (E - V)/(2 c^2)."""
    Mm = 1.0 + (E - V) / (2 * c * c)
    dP = 2 * Mm * c * r * Q + P
    dQ = -Q + r * (ll / (2 * Mm * c * r * r) + (V - E) / c) * P
    return dP, dQ


@njit(cache=True)
def _kh_start(r0, V0, E, l, Z, c):
    ll = l * (l + 1)
    za = Z / c
    gam = math.sqrt(max(ll + 1 - za * za, 0.01))
    Mm = 1.0 + (E - V0) / (2 * c * c)
    P = r0 ** gam
    Q = (gam - 1) * r0 ** (gam - 1) / (2 * Mm * c)
    return P, Q


@njit(cache=True)
def _kh_count(V, Vh, r, h, l, E, Z, c):
    M = r.shape[0]
    ll = l * (l + 1)
    h12 = h * h / 12.0
    P, Q = _kh_start(r[0], V[0], E, l, Z, c)
    nodes = 0
    for i in range(M - 1):
        if h12 * 2 * r[i + 1] ** 2 * (V[i + 1] - E) > 0.2:
            break
        rm = math.sqrt(r[i] * r[i + 1])
        a1, b1 = _kh_rhs(r[i], P, Q, E, V[i], ll, c)
        a2, b2 = _kh_rhs(rm, P + 0.5 * h * a1, Q + 0.5 * h * b1, E, Vh[i], ll, c)
        a3, b3 = _kh_rhs(rm, P + 0.5 * h * a2, Q + 0.5 * h * b2, E, Vh[i], ll, c)
        a4, b4 = _kh_rhs(r[i + 1], P + h * a3, Q + h * b3, E, V[i + 1], ll, c)
        Pn = P + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        Q = Q + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
        if (Pn < 0) != (P < 0) and Pn != 0:
            nodes += 1
        P = Pn
        if abs(P) > 1e150:
            P *= 1e-150
            Q *= 1e-150
    return nodes


@njit(cache=True)
def _kh_eigen(V, Vh, r, h, l, k, Elo, Ehi, tol, Z, c):
    for it in range(200):
        Em = 0.5 * (Elo + Ehi)
        if _kh_count(V, Vh, r, h, l, Em, Z, c) > k:
            Ehi = Em
        else:
            Elo = Em
        if Ehi - Elo < tol * max(1.0, abs(Em)):
            break
    return 0.5 * (Elo + Ehi)


@njit(cache=True)
def _kh_wavefunction(V, Vh, r, h, l, E, Z, c):
    M = r.shape[0]
    ll = l * (l + 1)
    g = (l + 0.5) ** 2 + 2 * r * r * (V - E)
    it = M - 2
    while it > 2 and g[it] > 0:
        it -= 1
    if it < 10:
        it = 10
    h12 = h * h / 12.0
    acc = 0.0
    ist = M - 1
    for i in range(it, M):
        if g[i] > 0:
            acc += math.sqrt(g[i]) * h
        if acc > 60.0 or h12 * g[i] > 0.2:
            ist = i
            break
    P = np.zeros(M)
    Q = np.zeros(M)
    P[0], Q[0] = _kh_start(r[0], V[0], E, l, Z, c)
    for i in range(it):
        rm = math.sqrt(r[i] * r[i + 1])
        p, q = P[i], Q[i]
        a1, b1 = _kh_rhs(r[i], p, q, E, V[i], ll, c)
        a2, b2 = _kh_rhs(rm, p + 0.5 * h * a1, q + 0.5 * h * b1, E, Vh[i], ll, c)
        a3, b3 = _kh_rhs(rm, p + 0.5 * h * a2, q + 0.5 * h * b2, E, Vh[i], ll, c)
        a4, b4 = _kh_rhs(r[i + 1], p + h * a3, q + h * b3, E, V[i + 1], ll, c)
        P[i + 1] = p + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        Q[i + 1] = q + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
        if abs(P[i + 1]) > 1e150:
            for j in range(i + 2):
                P[j] *= 1e-150
                Q[j] *= 1e-150
    Pin = np.zeros(M)
    Qin = np.zeros(M)
    Mm = 1.0 + (E - V[ist]) / (2 * c * c)
    kap = math.sqrt(max(-2 * Mm * E, 1e-12))
    Pin[ist] = 1e-200
    Qin[ist] = -(kap + 1.0 / r[ist]) * Pin[ist] / (2 * Mm * c)
    for i in range(ist, it, -1):
        rm = math.sqrt(r[i] * r[i - 1])
        p, q = Pin[i], Qin[i]
        hh = -h
        a1, b1 = _kh_rhs(r[i], p, q, E, V[i], ll, c)
        a2, b2 = _kh_rhs(rm, p + 0.5 * hh * a1, q + 0.5 * hh * b1, E, Vh[i - 1], ll, c)
        a3, b3 = _kh_rhs(rm, p + 0.5 * hh * a2, q + 0.5 * hh * b2, E, Vh[i - 1], ll, c)
        a4, b4 = _kh_rhs(r[i - 1], p + hh * a3, q + hh * b3, E, V[i - 1], ll, c)
        Pin[i - 1] = p + hh / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        Qin[i - 1] = q + hh / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
        if abs(Pin[i - 1]) > 1e150:
            for j in range(i - 1, ist + 1):
                Pin[j] *= 1e-150
                Qin[j] *= 1e-150
    sc = P[it] / Pin[it]
    for i in range(it + 1, ist + 1):
        P[i] = Pin[i] * sc
        Q[i] = Qin[i] * sc
    s = 0.0
    for i in range(M):
        w = 1.0 if (i == 0 or i == M - 1) else (4.0 if i % 2 == 1 else 2.0)
        s += w * (P[i] * P[i] + Q[i] * Q[i]) * r[i]
    s = math.sqrt(s * h / 3.0)
    return P / s, Q / s


# ============================================================================ XC: Slater + VWN5
_VWN = {
    "P": (0.0310907, -0.10498, 3.72744, 12.9352),
    "F": (0.01554535, -0.32500, 7.06042, 18.0578),
    "A": (-1.0 / (6 * math.pi ** 2), -0.0047584, 1.13107, 13.0045),
}
FPP0 = 4.0 / (9.0 * (2 ** (1 / 3) - 1))  # f''(0) = 1.709921


def _vwn_G(x, A, x0, b, c):
    X = x * x + b * x + c
    X0 = x0 * x0 + b * x0 + c
    Q = math.sqrt(4 * c - b * b)
    at = np.arctan(Q / (2 * x + b))
    G = A * (np.log(x * x / X) + 2 * b / Q * at
             - b * x0 / X0 * (np.log((x - x0) ** 2 / X) + 2 * (b + 2 * x0) / Q * at))
    dG = A * (2 / x - (2 * x + b) / X - 4 * b / ((2 * x + b) ** 2 + Q * Q)
              - b * x0 / X0 * (2 / (x - x0) - (2 * x + b) / X - 4 * (b + 2 * x0) / ((2 * x + b) ** 2 + Q * Q)))
    return G, dG


def lsda_xc(rho_up, rho_dn, exchange_only=False):
    """Return (eps_xc per electron, v_up, v_dn) for LSDA (Slater + VWN5)."""
    rho = rho_up + rho_dn
    tiny = 1e-30
    rho_s = np.maximum(rho, tiny)
    ru = np.maximum(rho_up, 0.0)
    rd = np.maximum(rho_dn, 0.0)
    # exchange
    cx = -(3.0 / 4.0) * (6.0 / math.pi) ** (1 / 3)
    ex_dens = cx * (ru ** (4 / 3) + rd ** (4 / 3))           # energy density
    vx_up = -(6.0 / math.pi * ru) ** (1 / 3)
    vx_dn = -(6.0 / math.pi * rd) ** (1 / 3)
    if exchange_only:
        return ex_dens / rho_s, vx_up, vx_dn
    rs = (3.0 / (4 * math.pi * rho_s)) ** (1 / 3)
    x = np.sqrt(rs)
    zeta = np.clip((ru - rd) / rho_s, -1.0, 1.0)
    eP, dP = _vwn_G(x, *_VWN["P"])
    eF, dF = _vwn_G(x, *_VWN["F"])
    eA, dA = _vwn_G(x, *_VWN["A"])
    f = ((1 + zeta) ** (4 / 3) + (1 - zeta) ** (4 / 3) - 2) / (2 ** (4 / 3) - 2)
    df = (4 / 3) * ((1 + zeta) ** (1 / 3) - (1 - zeta) ** (1 / 3)) / (2 ** (4 / 3) - 2)
    z4 = zeta ** 4
    ec = eP + eA * f / FPP0 * (1 - z4) + (eF - eP) * f * z4
    # d/drs = (1/(2x)) d/dx
    dec_dx = dP + dA * f / FPP0 * (1 - z4) + (dF - dP) * f * z4
    dec_drs = dec_dx / (2 * x)
    dec_dz = (eA / FPP0 * (df * (1 - z4) - 4 * zeta ** 3 * f)
              + (eF - eP) * (df * z4 + 4 * zeta ** 3 * f))
    vc_common = ec - rs / 3 * dec_drs
    vc_up = vc_common + (1 - zeta) * dec_dz
    vc_dn = vc_common - (1 + zeta) * dec_dz
    eps = ex_dens / rho_s + ec
    return eps, vx_up + vc_up, vx_dn + vc_dn


# ============================================================================ SCF
def occupations(shells, spin=True, remove=None, remove_amount=1.0):
    """List of (n, l, f_up, f_dn). Hund high-spin filling inside each subshell if spin=True.
    remove=(n,l): take remove_amount electrons out of that subshell (from the spin channel that
    loses an electron in the high-spin filling)."""
    occ = []
    for n, l, q in shells:
        q = float(q)
        deg = 2 * l + 1
        if remove is not None and (n, l) == tuple(remove):
            if spin:
                up = min(q, deg)
                dn = q - up
                if dn > 0:
                    dn -= remove_amount
                else:
                    up -= remove_amount
                occ.append((n, l, up, dn))
                continue
            q -= remove_amount
        if spin:
            up = min(q, deg)
            occ.append((n, l, up, q - up))
        else:
            occ.append((n, l, q / 2, q / 2))
    return [o for o in occ if o[2] + o[3] > 1e-12]


class Atom:
    def __init__(self, Z, occ, spin=True, grid=None, xc="lsda", scalar_rel=False, sic_orbital=None):
        self.Z = Z
        self.occ = occ
        self.spin = spin
        self.g = grid or Grid(Z)
        self.Nel = sum(o[2] + o[3] for o in occ)
        self.xc = xc
        self.scalar_rel = scalar_rel

    def _solve_orbitals(self, Vup, Vdn):
        g = self.g
        res = {}
        for (n, l, fu, fd) in self.occ:
            for s, V, f in ((0, Vup, fu), (1, Vdn, fd)):
                if f <= 0:
                    continue
                lp = (l + 0.5) ** 2
                Elo = 1.5 * float(np.min(V + lp / (2 * g.r ** 2))) - 1e-6
                k = n - l - 1
                if self.scalar_rel:
                    Vh = _half_values(V, g.r)
                    E = _kh_eigen(V, Vh, g.r, g.h, l, k, Elo, 50.0, 1e-13, float(self.Z), C_LIGHT)
                    P, Q = _kh_wavefunction(V, Vh, g.r, g.h, l, E, float(self.Z), C_LIGHT)
                    res[(n, l, s)] = (E, P, f, Q)
                else:
                    E = _eigen(V, g.r, g.h, l, k, Elo, 50.0, 1e-13)
                    P = _wavefunction(V, g.r, g.h, l, E)
                    res[(n, l, s)] = (E, P, f, None)
        return res

    def _density(self, orbs):
        r = self.g.r
        ru = np.zeros_like(r)
        rd = np.zeros_like(r)
        for (n, l, s), (E, P, f, Q) in orbs.items():
            d = f * P * P / (4 * math.pi * r * r)
            if Q is not None:
                d = d + f * Q * Q / (4 * math.pi * r * r)
            if s == 0:
                ru += d
            else:
                rd += d
        return ru, rd

    def _potential(self, ru, rd):
        r = self.g.r
        rho = ru + rd
        U = _poisson(rho, r, self.g.h, float(self.Nel))
        vH = U / r
        if not self.spin:
            eps, vu, vd = lsda_xc(rho / 2, rho / 2, exchange_only=(self.xc == "x"))
        else:
            eps, vu, vd = lsda_xc(ru, rd, exchange_only=(self.xc == "x"))
        vext = -self.Z / r
        return vext + vH + vu, vext + vH + vd, vH, eps

    def scf(self, tol=1e-9, maxit=300, beta=0.35, hist=7, verbose=False):
        g = self.g
        r = g.r
        # initial guess: Thomas-Fermi-like screened Coulomb
        Zeff = self.Z - (self.Nel - 1) * (1 - np.exp(-r * 1.5 * self.Z ** (1 / 3)))
        V0 = -np.maximum(Zeff, 1.0 if self.Nel < self.Z + 0.5 else 0.5) / r
        Vin = np.concatenate([V0, V0])
        Xs, Fs = [], []
        Eold = 0.0
        M = g.M
        w = g.wr * r * r  # weight for residual norm
        for it in range(maxit):
            orbs = self._solve_orbitals(Vin[:M], Vin[M:])
            ru, rd = self._density(orbs)
            vu, vd, vH, eps = self._potential(ru, rd)
            Vout = np.concatenate([vu, vd])
            F = Vout - Vin
            # energy
            E = self._energy(orbs, Vin[:M], Vin[M:], ru, rd, vH, eps)
            res = math.sqrt(np.dot(np.concatenate([w, w]), F * F))
            if verbose:
                print(it, E, res)
            if res < tol and abs(E - Eold) < tol * 10:
                break
            Eold = E
            Xs.append(Vin.copy())
            Fs.append(F.copy())
            if len(Xs) > hist:
                Xs.pop(0)
                Fs.pop(0)
            if len(Xs) >= 2:
                # Pulay DIIS
                m = len(Fs)
                W = np.concatenate([w, w])
                B = np.empty((m + 1, m + 1))
                for i in range(m):
                    for j in range(m):
                        B[i, j] = np.dot(W, Fs[i] * Fs[j])
                B[m, :] = -1
                B[:, m] = -1
                B[m, m] = 0
                rhs = np.zeros(m + 1)
                rhs[m] = -1
                try:
                    cvec = np.linalg.solve(B, rhs)[:m]
                    Xb = sum(ci * x for ci, x in zip(cvec, Xs))
                    Fb = sum(ci * f for ci, f in zip(cvec, Fs))
                    Vin = Xb + beta * Fb
                except np.linalg.LinAlgError:
                    Vin = Vin + beta * F
            else:
                Vin = Vin + beta * F
        self.converged = res < tol * 100
        self.iterations = it
        self.orbs = orbs
        self.E = E
        self.res = res
        return E

    def _energy(self, orbs, vin_u, vin_d, ru, rd, vH, eps):
        g = self.g
        r = g.r
        four_pi_r2 = 4 * math.pi * r * r
        band = sum(v[2] * v[0] for v in orbs.values())
        Ts = band - g.integrate(four_pi_r2 * (ru * vin_u + rd * vin_d))
        rho = ru + rd
        Eext = g.integrate(four_pi_r2 * rho * (-self.Z / r))
        EH = 0.5 * g.integrate(four_pi_r2 * rho * vH)
        Exc = g.integrate(four_pi_r2 * rho * eps)
        self.parts = dict(Ts=Ts, Eext=Eext, EH=EH, Exc=Exc, band=band)
        return Ts + Eext + EH + Exc

    def eigenvalue(self, n, l, s=None):
        if s is None:
            cands = [(k, v) for k, v in self.orbs.items() if k[:2] == (n, l)]
            # minority channel if occupied (the one losing an electron in high-spin removal)
            k, v = max(cands, key=lambda kv: kv[0][2]) if len(cands) > 1 and cands[-1][1][2] > 0 else cands[0]
            return v[0]
        return self.orbs[(n, l, s)][0]


def total_energy(Z, shells, spin=True, xc="lsda", scalar_rel=False, remove=None, remove_amount=1.0,
                 return_atom=False, **kw):
    if not shells or sum(q for *_, q in shells) - (remove_amount if remove else 0) <= 1e-12:
        return (0.0, None) if return_atom else 0.0
    occ = occupations(shells, spin=spin, remove=remove, remove_amount=remove_amount)
    a = Atom(Z, occ, spin=spin, xc=xc, scalar_rel=scalar_rel)
    E = a.scf(**kw)
    return (E, a) if return_atom else E


if __name__ == "__main__":
    import time
    refs = {"He": (2, [(1, 0, 2)], -2.834836), "Be": (4, [(1, 0, 2), (2, 0, 2)], -14.447209),
            "Ne": (10, [(1, 0, 2), (2, 0, 2), (2, 1, 6)], -128.233481),
            "Ar": (18, [(1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 6)], -525.946195)}
    for name, (Z, sh, ref) in refs.items():
        t = time.time()
        E, a = total_energy(Z, sh, spin=False, return_atom=True)
        print(f"{name}: E_LDA = {E:.6f}  NIST = {ref:.6f}  diff = {E - ref:.2e}  it={a.iterations} "
              f"t={time.time() - t:.2f}s", {k: round(v[0], 6) for k, v in a.orbs.items()})
    # LSD hydrogen: NIST -0.478671 ; LDA -0.445671
    print("H LSD", total_energy(1, [(1, 0, 1)], spin=True), "(NIST -0.478671)")
    print("H LDA", total_energy(1, [(1, 0, 1)], spin=False), "(NIST -0.445671)")
