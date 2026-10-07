"""Exact hydrogenic (Z = 1) radial Slater integrals, in rational arithmetic.

For a one-electron orbital nl of hydrogen (Z = 1, atomic units) the reduced radial function is

    P_nl(r) = r R_nl(r) = N_nl (2/n)^l  r^(l+1) e^(-r/n) L_{n-l-1}^{(2l+1)}(2r/n),
    N_nl^2  = (2/n)^3 (n-l-1)! / (2n (n+l)!)

(associated Laguerre polynomial in the modern convention).  The general two-electron radial
integral is

    R^k(ab; cd) = Int Int P_a(r1) P_c(r1) P_b(r2) P_d(r2) r_<^k / r_>^(k+1) dr1 dr2

with F^k(a,b) = R^k(ab; ab) (direct) and G^k(a,b) = R^k(ab; ba) (exchange).  Because each
P_a P_c is (polynomial) x exp(-beta r) with rational beta = 1/n_a + 1/n_c, every R^k is an exact
rational number times sqrt(N_a^2 N_b^2 N_c^2 N_d^2); for F^k and G^k the square root disappears
and the result is rational.  We evaluate with fractions.Fraction (exact), using

    int_0^x t^s e^{-bt} dt = s!/b^{s+1} [1 - e^{-bx} sum_{j<=s} (bx)^j/j!]
    int_x^inf t^s e^{-bt} dt = s!/b^{s+1} e^{-bx} sum_{j<=s} (bx)^j/j!

For a hydrogenic ion of charge Z every integral scales as R^k(Z) = Z R^k(1).

Results are cached in cache/slater_integrals.json (keys "F|n l|n' l'|k", "G|...").
"""
import json
import math
import os
from fractions import Fraction
from functools import lru_cache

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "cache", "slater_integrals.json")


def _fact(n):
    return math.factorial(n)


@lru_cache(maxsize=None)
def poly_P(n, l):
    """Return (coeffs, beta, Nsq): P_nl(r) = sqrt(Nsq) * sum_i coeffs[i] r^i * exp(-r/n)."""
    m = n - l - 1
    alpha = 2 * l + 1
    two_n = Fraction(2, n)
    coeffs = {}
    for i in range(m + 1):
        c = Fraction((-1) ** i * math.comb(m + alpha, m - i), _fact(i)) * two_n ** i
        coeffs[i + l + 1] = c * two_n ** l
    Nsq = two_n ** 3 * Fraction(_fact(n - l - 1), 2 * n * _fact(n + l))
    maxp = max(coeffs)
    arr = [coeffs.get(i, Fraction(0)) for i in range(maxp + 1)]
    return tuple(arr), Fraction(1, n), Nsq


def _pmul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            if y:
                out[i + j] += x * y
    return out


def _gint(q, g):
    """int_0^inf r^q e^{-g r} dr, q >= 0 integer."""
    if q < 0:
        raise ValueError("divergent integral")
    return Fraction(_fact(q)) / g ** (q + 1)


def _Rk_raw(f1, b1, f2, b2, k):
    """Unnormalised double integral for densities f1(r1) e^{-b1 r1}, f2(r2) e^{-b2 r2}."""
    tot = Fraction(0)
    bb = b1 + b2
    for m, A in enumerate(f1):
        if A == 0:
            continue
        for p, B in enumerate(f2):
            if B == 0:
                continue
            # term 1: r1 < r2 :  r2^{-(k+1)} int_0^{r2} r1^{m+k}
            s = m + k
            pref = A * B * Fraction(_fact(s)) / b1 ** (s + 1)
            t = _gint(p - k - 1, b2)
            acc = Fraction(0)
            bj = Fraction(1)
            for j in range(s + 1):
                acc += bj / _fact(j) * _gint(p - k - 1 + j, bb)
                bj *= b1
            tot += pref * (t - acc)
            # term 2: r1 > r2 : r2^k int_{r2}^inf r1^{m-k-1}
            s = m - k - 1
            pref = A * B * Fraction(_fact(s)) / b1 ** (s + 1)
            acc = Fraction(0)
            bj = Fraction(1)
            for j in range(s + 1):
                acc += bj / _fact(j) * _gint(p + k + j, bb)
                bj *= b1
            tot += pref * acc
    return tot


@lru_cache(maxsize=None)
def Rk_general(a, b, c, d, k):
    """R^k(ab; cd) for hydrogenic Z=1 orbitals a=(n,l) ... Returns (rational_part, Nsq_product)
    such that R^k = rational_part * sqrt(Nsq_product)."""
    pa, ba, Na = poly_P(*a)
    pb, bb_, Nb = poly_P(*b)
    pc, bc, Nc = poly_P(*c)
    pd, bd, Nd = poly_P(*d)
    f1 = _pmul(pa, pc)
    f2 = _pmul(pb, pd)
    raw = _Rk_raw(f1, ba + bc, f2, bb_ + bd, k)
    return raw, Na * Nb * Nc * Nd


def Rk_float(a, b, c, d, k):
    raw, nsq = Rk_general(tuple(a), tuple(b), tuple(c), tuple(d), k)
    return float(raw) * math.sqrt(float(nsq))


_cache = None


def _load():
    global _cache
    if _cache is None:
        _cache = {}
        if os.path.exists(CACHE):
            with open(CACHE, encoding="utf-8") as f:
                _cache = {k: Fraction(v) for k, v in json.load(f).items()}
    return _cache


def save_cache():
    if _cache is None:
        return
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    with open(CACHE, "w", encoding="utf-8") as f:
        json.dump({k: str(v) for k, v in sorted(_cache.items())}, f, indent=0)


def _key(kind, a, b, k):
    # F and G are symmetric in a <-> b
    a, b = sorted([tuple(a), tuple(b)])
    return f"{kind}|{a[0]} {a[1]}|{b[0]} {b[1]}|{k}", a, b


def Fk(a, b, k):
    """Exact direct Slater integral F^k(a,b) (Fraction, hartree, Z=1)."""
    key, a, b = _key("F", a, b, k)
    c = _load()
    if key not in c:
        raw, nsq = Rk_general(a, b, a, b, k)
        c[key] = raw * nsq_sqrt_rational(nsq)
    return c[key]


def Gk(a, b, k):
    """Exact exchange Slater integral G^k(a,b) (Fraction, hartree, Z=1)."""
    key, a, b = _key("G", a, b, k)
    c = _load()
    if key not in c:
        raw, nsq = Rk_general(a, b, b, a, k)
        c[key] = raw * nsq_sqrt_rational(nsq)
    return c[key]


def nsq_sqrt_rational(nsq):
    """nsq = N_a^2 N_b^2 N_a^2 N_b^2 is a perfect square of a rational: return its sqrt."""
    num, den = nsq.numerator, nsq.denominator
    rn, rd = math.isqrt(num), math.isqrt(den)
    assert rn * rn == num and rd * rd == den, nsq
    return Fraction(rn, rd)


if __name__ == "__main__":
    print("F0(1s,1s) =", Fk((1, 0), (1, 0), 0), "(expect 5/8)")
    print("F0(1s,2s) =", Fk((1, 0), (2, 0), 0), "(expect 17/81)")
    print("G0(1s,2s) =", Gk((1, 0), (2, 0), 0), "(expect 16/729)")
    print("F0(1s,2p) =", Fk((1, 0), (2, 1), 0), "(expect 59/243)")
    print("G1(1s,2p) =", Gk((1, 0), (2, 1), 1), "(expect 112/6561)")
    print("F0(2s,2s) =", Fk((2, 0), (2, 0), 0), "(expect 77/512)")
    print("F2(2p,2p) =", Fk((2, 1), (2, 1), 2), "(expect 45/512)")
    # numerical cross-check of one general integral
    import numpy as np
    r = np.linspace(1e-8, 200, 400001)

    def P(n, l):
        from scipy.special import genlaguerre
        Nn = math.sqrt((2 / n) ** 3 * math.factorial(n - l - 1) / (2 * n * math.factorial(n + l)))
        return Nn * r * (2 * r / n) ** l * np.exp(-r / n) * genlaguerre(n - l - 1, 2 * l + 1)(2 * r / n)
    a, b = (3, 2), (4, 1)
    k = 1
    rho1 = P(*a) * P(*b)
    from scipy.integrate import cumulative_trapezoid as ct
    inner = ct(rho1 * r ** k, r, initial=0)
    outer = ct((rho1 / r ** (k + 1))[::-1], r[::-1], initial=0)[::-1] * -1
    Y = inner / r ** (k + 1) + outer * r ** k
    num = np.trapezoid(rho1 * Y, r)
    print("G1(3d,4p) exact", float(Gk(a, b, 1)), "numeric", num)
    save_cache()
