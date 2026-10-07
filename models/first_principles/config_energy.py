"""Exact zeroth- and first-order energies of the hydrogenic 1/Z (Hylleraas-Layzer) expansion.

After the scaling r -> r/Z the non-relativistic N-electron Hamiltonian is
H = Z^2 [ H0 + Z^{-1} V ],  H0 = sum_i (-1/2 nabla_i^2 - 1/r_i),  V = sum_{i<j} 1/r_ij,
so for a fixed configuration and term  E(Z) = Z^2 E0 + Z E1 + E2 + E3/Z + ...  (hartree), with

    E0 = -sum_i q_i / (2 n_i^2)
    E1 = < V >   over hydrogenic Z=1 orbitals (first-order, degenerate perturbation theory).

This module computes E1 three ways:

  * E_av  : Slater/Cowan configuration-average energy (exact rational)
  * E1_sc : single-configuration energy of the Hund's-rule LS term (max S, then max L);
            exact rational when the highest-weight (M_S=S, M_L=L) space is one determinant,
            otherwise lowest eigenvalue of V in the highest-weight subspace of the configuration.
            For several open subshells this is high-spin (all open-subshell spins parallel) coupling.
  * E1_cx : Layzer complex: lowest eigenvalue of V within the full hydrogenic complex (all
            configurations with the same number of electrons in each n-shell, e.g. 1s2 2s2 + 1s2 2p2
            for Be-like) restricted to the S, L, parity of the reference Hund term.  Only computed
            when exactly one n-shell is incompletely filled and the determinant space is small.

All quantities are in hartree for Z = 1 (multiply E1 by Z for charge Z).
"""
import itertools
import math
from fractions import Fraction
from functools import lru_cache

import numpy as np

from slater_integrals import Fk, Gk, Rk_float, save_cache
from angular import ck_diag, ck_sq, ck_float, threej0_sq

MAX_COMPLEX_COMBOS = 60000


# ----------------------------------------------------------------------------- basics
def E0(shells):
    return -sum(Fraction(q, 2 * n * n) for n, l, q in shells)


def is_closed(l, q):
    return q == 2 * (2 * l + 1)


def E_average(shells):
    """Cowan (Theory of Atomic Structure and Spectra, eq. 6.38) configuration-average <V>, exact."""
    E = Fraction(0)
    sh = [(n, l, q) for n, l, q in shells if q > 0]
    for i, (n, l, q) in enumerate(sh):
        a = (n, l)
        if q > 1:
            t = Fk(a, a, 0)
            for k in range(2, 2 * l + 1, 2):
                t -= Fraction(2 * l + 1, 4 * l + 1) * threej0_sq(l, k, l) * Fk(a, a, k)
            E += Fraction(q * (q - 1), 2) * t
        for n2, l2, q2 in sh[i + 1:]:
            b = (n2, l2)
            t = Fk(a, b, 0)
            for k in range(abs(l - l2), l + l2 + 1):
                if (l + l2 + k) % 2 == 0:
                    t -= Fraction(1, 2) * threej0_sq(l, k, l2) * Gk(a, b, k)
            E += q * q2 * t
    return E


def core_energy(closed):
    """<V> among closed subshells (single determinant == configuration average)."""
    return E_average(closed)


def u_core(a, closed):
    """Interaction of one electron in subshell a=(n,l) with all closed subshells (exact)."""
    n, l = a
    u = Fraction(0)
    for n2, l2, q2 in closed:
        b = (n2, l2)
        t = Fk(a, b, 0)
        for k in range(abs(l - l2), l + l2 + 1):
            if (l + l2 + k) % 2 == 0:
                t -= Fraction(1, 2) * threej0_sq(l, k, l2) * Gk(a, b, k)
        u += q2 * t
    return u


# ----------------------------------------------------------------------------- spin orbitals
class SpinOrbitals:
    """Spin-orbital index <-> (n, l, m, ms2) with ms2 = +1 / -1."""

    def __init__(self, subshells):
        self.orbs = []
        for n, l in subshells:
            for m in range(l, -l - 1, -1):
                for s in (1, -1):
                    self.orbs.append((n, l, m, s))
        self.index = {o: i for i, o in enumerate(self.orbs)}
        self.M = len(self.orbs)


def _sign_remove(det, j):
    """det: sorted tuple. Return (sign, newdet) for a_j det, or (0,None)."""
    if j not in det:
        return 0, None
    pos = det.index(j)
    return (-1) ** pos, det[:pos] + det[pos + 1:]


def _sign_add(det, k):
    if k in det:
        return 0, None
    pos = 0
    while pos < len(det) and det[pos] < k:
        pos += 1
    return (-1) ** pos, det[:pos] + (k,) + det[pos:]


def _apply_one(det, k, j):
    """a+_k a_j det -> (sign, newdet)."""
    s1, d1 = _sign_remove(det, j)
    if not s1:
        return 0, None
    s2, d2 = _sign_add(d1, k)
    if not s2:
        return 0, None
    return s1 * s2, d2


def _two_e(so, p, q, r, s):
    """<pq|rs> (physicist; p,r electron 1) over spin orbitals, float."""
    n1, l1, m1, s1 = so.orbs[p]
    n2, l2, m2, s2 = so.orbs[q]
    n3, l3, m3, s3 = so.orbs[r]
    n4, l4, m4, s4 = so.orbs[s]
    if s1 != s3 or s2 != s4 or m1 + m2 != m3 + m4:
        return 0.0
    v = 0.0
    for k in range(max(abs(l1 - l3), abs(l2 - l4)), min(l1 + l3, l2 + l4) + 1):
        if (l1 + l3 + k) % 2 or (l2 + l4 + k) % 2:
            continue
        a = ck_float(l1, m1, l3, m3, k)
        b = ck_float(l4, m4, l2, m2, k)
        if a == 0 or b == 0:
            continue
        v += a * b * Rk_float((n1, l1), (n2, l2), (n3, l3), (n4, l4), k)
    return v


def _pair_exact(so, p, q):
    """Exact J_pq - delta_spin K_pq for two occupied spin orbitals."""
    n1, l1, m1, s1 = so.orbs[p]
    n2, l2, m2, s2 = so.orbs[q]
    a, b = (n1, l1), (n2, l2)
    J = Fraction(0)
    for k in range(0, 2 * min(l1, l2) + 1, 2):
        J += ck_diag(l1, m1, k) * ck_diag(l2, m2, k) * Fk(a, b, k)
    K = Fraction(0)
    if s1 == s2:
        for k in range(abs(l1 - l2), l1 + l2 + 1):
            if (l1 + l2 + k) % 2 == 0:
                K += ck_sq(l1, m1, l2, m2, k) * Gk(a, b, k)
    return J - K


def det_energy_exact(so, det):
    E = Fraction(0)
    for i, p in enumerate(det):
        for q in det[i + 1:]:
            E += _pair_exact(so, p, q)
    return E


class TwoBody:
    """Antisymmetrised <pq||rs> table over a spin-orbital set, grouped by (r,s)."""

    def __init__(self, so):
        self.so = so
        M = so.M
        pairs = [(p, q) for p in range(M) for q in range(p + 1, M)]
        key = lambda pq: (so.orbs[pq[0]][2] + so.orbs[pq[1]][2], so.orbs[pq[0]][3] + so.orbs[pq[1]][3])
        groups = {}
        for pq in pairs:
            groups.setdefault(key(pq), []).append(pq)
        self.table = {}
        for rs in pairs:
            lst = []
            r, s = rs
            for p, q in groups[key(rs)]:
                v = _two_e(so, p, q, r, s) - _two_e(so, p, q, s, r)
                if abs(v) > 1e-15:
                    lst.append((p, q, v))
            self.table[rs] = lst


def build_H(so, tb, dets, onebody):
    """V matrix (+ diagonal one-body onebody[orb]) in the given determinant list."""
    idx = {d: i for i, d in enumerate(dets)}
    H = np.zeros((len(dets), len(dets)))
    for j, D in enumerate(dets):
        H[j, j] += sum(onebody[o] for o in D)
        Dset = set(D)
        for a in range(len(D)):
            r = D[a]
            for b in range(a + 1, len(D)):
                s = D[b]
                # a_s a_r |D>  (apply a_r first)
                s1, d1 = _sign_remove(D, r)
                s2, d2 = _sign_remove(d1, s)
                sg = s1 * s2
                for p, q, v in tb.table[(r, s)]:
                    if (p in Dset and p not in (r, s)) or (q in Dset and q not in (r, s)):
                        continue
                    t1, e1 = _sign_add(d2, q)
                    t2, e2 = _sign_add(e1, p)
                    i = idx.get(e2)
                    if i is not None:
                        H[i, j] += sg * t1 * t2 * v
    return H


def _raising(so, dets, kind):
    """Matrix of S+ (kind='S') or L+ (kind='L') from dets to image dets."""
    img = {}
    entries = []
    for j, D in enumerate(dets):
        for o in D:
            n, l, m, s = so.orbs[o]
            if kind == "S":
                if s != -1:
                    continue
                t = so.index[(n, l, m, 1)]
                coef = 1.0
            else:
                if m == l:
                    continue
                t = so.index[(n, l, m + 1, s)]
                coef = math.sqrt(l * (l + 1) - m * (m + 1))
            sg, nd = _apply_one(D, t, o)
            if sg:
                i = img.setdefault(nd, len(img))
                entries.append((i, j, sg * coef))
    A = np.zeros((len(img), len(dets)))
    for i, j, v in entries:
        A[i, j] += v
    return A


def highest_weight_basis(so, dets):
    A = np.vstack([_raising(so, dets, "S"), _raising(so, dets, "L")])
    if A.shape[0] == 0:
        return np.eye(len(dets))
    u, sv, vt = np.linalg.svd(A, full_matrices=True)
    rank = int((sv > 1e-9).sum())
    return vt[rank:].T


def _qnums(so, D):
    ms2 = sum(so.orbs[o][3] for o in D)
    ml = sum(so.orbs[o][2] for o in D)
    par = sum(so.orbs[o][1] for o in D) % 2
    return ms2, ml, par


def _subshell_dets(so, n, l, q):
    orbs = [so.index[(n, l, m, s)] for m in range(l, -l - 1, -1) for s in (1, -1)]
    return [tuple(sorted(c)) for c in itertools.combinations(orbs, q)]


def hund_target(open_shells):
    """(2*S, L) of the Hund (max S, then max L) term for parallel-coupled open subshells."""
    S2 = 0
    L = 0
    for n, l, q in open_shells:
        nup = min(q, 2 * l + 1)
        ndn = q - nup
        S2 += nup - ndn
        ms = list(range(l, -l - 1, -1))
        L += sum(ms[:nup]) + sum(ms[:ndn])
    return S2, L


# ----------------------------------------------------------------------------- E1 drivers
@lru_cache(maxsize=None)
def E1_single_config(shells):
    """Hund-term energy of a single configuration. Returns dict with exact pieces."""
    shells = tuple(s for s in shells if s[2] > 0)
    closed = [s for s in shells if is_closed(s[1], s[2])]
    opn = [s for s in shells if not is_closed(s[1], s[2])]
    Ecore = core_energy(closed)
    Eav = E_average(shells)
    out = {"E_av": Eav, "E_core": Ecore}
    if not opn:
        out.update(E1=Ecore, exact=True, dim=1)
        return out
    so = SpinOrbitals([(n, l) for n, l, q in opn])
    S2, L = hund_target(opn)
    # enumerate dets with fixed subshell occupations and M_S = S, M_L = L
    per = []
    for n, l, q in opn:
        per.append([(d, _qnums(so, d)) for d in _subshell_dets(so, n, l, q)])
    dets = []

    def rec(i, cur, ms2, ml):
        if i == len(per):
            if ms2 == S2 and ml == L:
                dets.append(tuple(sorted(cur)))
            return
        for d, (a, b, _) in per[i]:
            rec(i + 1, cur + list(d), ms2 + a, ml + b)
    rec(0, [], 0, 0)
    ucore = {o: u_core(so.orbs[o][:2], closed) for o in range(so.M)}
    if len(dets) == 1:
        D = dets[0]
        E = Ecore + sum(ucore[o] for o in D) + det_energy_exact(so, D)
        out.update(E1=E, exact=True, dim=1, S2=S2, L=L)
        return out
    tb = TwoBody(so)
    H = build_H(so, tb, dets, {o: float(v) for o, v in ucore.items()})
    Q = highest_weight_basis(so, dets)
    w = np.linalg.eigvalsh(Q.T @ H @ Q)
    out.update(E1=float(Ecore) + float(w[0]), exact=False, dim=len(dets), hw_dim=Q.shape[1], S2=S2, L=L)
    return out


@lru_cache(maxsize=None)
def E1_complex(shells):
    """Layzer-complex first-order energy (None if not applicable / too large)."""
    shells = tuple(s for s in shells if s[2] > 0)
    nq = {}
    for n, l, q in shells:
        nq[n] = nq.get(n, 0) + q
    open_n = [n for n, q in nq.items() if q < 2 * n * n]
    if len(open_n) == 0:
        sc = E1_single_config(shells)
        return {"E1": float(sc["E1"]), "dim": 1, "hw_dim": 1, "ref_weight": 1.0}
    if len(open_n) != 1:
        return None
    n0 = open_n[0]
    q0 = nq[n0]
    if math.comb(2 * n0 * n0, q0) > MAX_COMPLEX_COMBOS:
        return None
    closed = [s for s in shells if s[0] != n0]
    opn = [s for s in shells if not is_closed(s[1], s[2])]
    S2, L = hund_target(opn)
    par = sum(l * q for n, l, q in shells if n == n0) % 2
    so = SpinOrbitals([(n0, l) for l in range(n0)])
    dets = []
    for c in itertools.combinations(range(so.M), q0):
        if _qnums(so, c) == (S2, L, par):
            dets.append(c)
    Ecore = core_energy(closed)
    ucore = {o: float(u_core(so.orbs[o][:2], closed)) for o in range(so.M)}
    tb = _tb_cache(n0)
    H = build_H(so, tb, dets, ucore)
    Q = highest_weight_basis(so, dets)
    Hq = Q.T @ H @ Q
    w, v = np.linalg.eigh(Hq)
    # weight of the reference configuration in the ground eigenvector
    ref_occ = {(n, l): q for n, l, q in shells if n == n0}
    c = Q @ v[:, 0]
    wref = 0.0
    for i, D in enumerate(dets):
        occ = {}
        for o in D:
            k = so.orbs[o][:2]
            occ[k] = occ.get(k, 0) + 1
        if all(occ.get((n0, l), 0) == ref_occ.get((n0, l), 0) for l in range(n0)):
            wref += c[i] ** 2
    return {"E1": float(Ecore) + float(w[0]), "dim": len(dets), "hw_dim": Q.shape[1], "ref_weight": wref}


_TB = {}


def _tb_cache(n0):
    if n0 not in _TB:
        _TB[n0] = TwoBody(SpinOrbitals([(n0, l) for l in range(n0)]))
    return _TB[n0]


def config_coeffs(shells):
    """Summary dict for a configuration: E0, E_av, E1_sc (+exact flag), E1_cx (or None)."""
    shells = tuple(sorted(tuple(s) for s in shells if s[2] > 0))
    if not shells:
        return {"E0": Fraction(0), "E_av": Fraction(0), "E1_sc": Fraction(0), "sc_exact": True,
                "E1_cx": 0.0, "cx_ref_weight": 1.0}
    sc = E1_single_config(shells)
    cx = E1_complex(shells)
    return {"E0": E0(shells), "E_av": sc["E_av"], "E1_sc": sc["E1"], "sc_exact": sc["exact"],
            "E1_cx": None if cx is None else cx["E1"],
            "cx_ref_weight": None if cx is None else cx["ref_weight"]}


if __name__ == "__main__":
    import time
    t = time.time()
    tests = {
        "He 1s2": ((1, 0, 2),),
        "Li 1s2 2s": ((1, 0, 2), (2, 0, 1)),
        "Li 1s2 2p": ((1, 0, 2), (2, 1, 1)),
        "Be 1s2 2s2": ((1, 0, 2), (2, 0, 2)),
        "B 2s2 2p": ((1, 0, 2), (2, 0, 2), (2, 1, 1)),
        "C 2s2 2p2": ((1, 0, 2), (2, 0, 2), (2, 1, 2)),
        "N 2s2 2p3": ((1, 0, 2), (2, 0, 2), (2, 1, 3)),
        "Ne 2s2 2p6": ((1, 0, 2), (2, 0, 2), (2, 1, 6)),
        "Mg-like 3s2": ((1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2)),
        "Fe 3d6 4s2": ((1, 0, 2), (2, 0, 2), (2, 1, 6), (3, 0, 2), (3, 1, 6), (3, 2, 6), (4, 0, 2)),
    }
    for name, sh in tests.items():
        c = config_coeffs(sh)
        print(f"{name:14s} E0={float(c['E0']):10.5f} E_av={float(c['E_av']):.6f} "
              f"E1_sc={c['E1_sc']} ({float(c['E1_sc']):.6f}) E1_cx={c['E1_cx']} w={c['cx_ref_weight']}")
    print("time", time.time() - t)
    save_cache()
